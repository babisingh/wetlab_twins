from __future__ import annotations

import asyncio
import logging
import os

import httpx

from aixbio.models.structure import StructureResult

logger = logging.getLogger(__name__)

ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

_TIMEOUT = 120
_MAX_RETRIES = 3
_BACKOFF_BASE = 2.0

# ESMFold works best for sequences between 1 and 400 residues.
_MAX_SEQUENCE_LENGTH = 400


def _parse_plddt(pdb_text: str) -> float:
    """Average B-factor from ESMFold PDB output; B-factor column holds pLDDT (0-100)."""
    scores: list[float] = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                scores.append(float(line[60:66].strip()))
            except (ValueError, IndexError):
                pass
    return round(sum(scores) / len(scores), 2) if scores else 0.0


async def fold_sequence(
    chain_id: str,
    aa_sequence: str,
    output_dir: str,
) -> tuple[StructureResult, str | None]:
    """Fold aa_sequence using the ESMFold public API.

    Returns (StructureResult, pdb_file_path).  pdb_file_path is None on failure.
    pLDDT values come from the B-factor column of the returned PDB (0-100 scale).

    ESMFold folds the *actual pipeline-produced sequence*, so this is a genuine
    validation of what the pipeline output represents structurally — unlike
    looking up a pre-computed AFDB entry keyed on the input UniProt ID.
    """
    if len(aa_sequence) > _MAX_SEQUENCE_LENGTH:
        logger.warning(
            f"Chain {chain_id} is {len(aa_sequence)} residues — exceeds ESMFold "
            f"recommended limit of {_MAX_SEQUENCE_LENGTH}. Skipping ESMFold."
        )
        return _failed_result(chain_id), None

    last_error: Exception | None = None
    pdb_text: str | None = None

    for attempt in range(_MAX_RETRIES):
        try:
            async with httpx.AsyncClient(timeout=_TIMEOUT) as client:
                resp = await client.post(
                    ESMFOLD_API,
                    content=aa_sequence,
                    headers={"Content-Type": "text/plain"},
                )
                resp.raise_for_status()
                pdb_text = resp.text
            break
        except httpx.HTTPStatusError as exc:
            last_error = exc
            if exc.response.status_code < 500:
                logger.warning(f"ESMFold API 4xx for chain {chain_id}: {exc}")
                return _failed_result(chain_id), None
            _log_retry(chain_id, attempt, exc)
            await asyncio.sleep(_BACKOFF_BASE * (2**attempt))
        except (httpx.TimeoutException, httpx.ConnectError) as exc:
            last_error = exc
            _log_retry(chain_id, attempt, exc)
            await asyncio.sleep(_BACKOFF_BASE * (2**attempt))

    if pdb_text is None:
        logger.warning(
            f"ESMFold failed for chain {chain_id} after {_MAX_RETRIES} attempts: {last_error}"
        )
        return _failed_result(chain_id), None

    os.makedirs(output_dir, exist_ok=True)
    pdb_path = os.path.join(output_dir, f"esm_{chain_id}.pdb")
    with open(pdb_path, "w") as f:
        f.write(pdb_text)

    plddt = _parse_plddt(pdb_text)
    logger.info(f"ESMFold: chain {chain_id} ({len(aa_sequence)} aa), pLDDT={plddt:.1f}")

    return StructureResult(
        id=chain_id,
        plddt_mean=plddt,
        rmsd_to_ref=None,
        perplexity=None,
        structure_file=pdb_path,
        method="esmfold",
    ), pdb_path


def _failed_result(chain_id: str) -> StructureResult:
    return StructureResult(
        id=chain_id,
        plddt_mean=0.0,
        rmsd_to_ref=None,
        perplexity=None,
        structure_file="",
        method="esmfold_failed",
    )


def _log_retry(chain_id: str, attempt: int, exc: Exception) -> None:
    wait = _BACKOFF_BASE * (2**attempt)
    logger.warning(
        f"ESMFold attempt {attempt + 1}/{_MAX_RETRIES} failed for {chain_id}: {exc}. "
        f"Retrying in {wait:.1f}s..."
    )