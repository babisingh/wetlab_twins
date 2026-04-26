from __future__ import annotations

import logging
import os
import tempfile

import httpx
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.Polypeptide import protein_letters_3to1

from aixbio.models.structure import StructureResult

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.com/api"

_TIMEOUT = 60
_MAX_RETRIES = 3
_BACKOFF_BASE = 1.0

# Sequence aligner shared across calls — stateless, thread-safe after setup.
_ALIGNER = PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.match_score = 2
_ALIGNER.mismatch_score = -1
_ALIGNER.open_gap_score = -4
_ALIGNER.extend_gap_score = -0.5


def _load_structure(path: str, entry_id: str):
    """Parse CIF or PDB based on file extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True).get_structure(entry_id, path)
    return PDBParser(QUIET=True).get_structure(entry_id, path)


def _get_chain_residues(structure) -> dict[str, list]:
    """Return {chain_id: [standard residues]} from the first model."""
    result: dict[str, list] = {}
    for model in structure:
        for chain in model:
            residues = [r for r in chain if r.id[0] == " " and "CA" in r]
            if residues:
                result[chain.id] = residues
        break  # first model only
    return result


def _residues_to_seq(residues: list) -> str:
    return "".join(protein_letters_3to1.get(r.resname, "X") for r in residues)


def _best_ref_chain(
    query_aa: str,
    ref_chains: dict[str, list],
) -> tuple[str, list] | tuple[None, None]:
    """Select the reference chain whose sequence best aligns to query_aa."""
    best_id: str | None = None
    best_score = -float("inf")
    for chain_id, residues in ref_chains.items():
        ref_seq = _residues_to_seq(residues)
        try:
            score = _ALIGNER.score(query_aa, ref_seq)
        except Exception:
            continue
        if score > best_score:
            best_score = score
            best_id = chain_id
    if best_id is None:
        return None, None
    return best_id, ref_chains[best_id]


def _aligned_ca_pairs(
    query_seq: str,
    ref_seq: str,
    query_res: list,
    ref_res: list,
) -> tuple[list, list]:
    """Return paired (query_CA, ref_CA) atom lists for sequence-aligned positions."""
    try:
        alignments = list(_ALIGNER.align(query_seq, ref_seq))
    except Exception:
        return [], []
    if not alignments:
        return [], []

    aln = alignments[0]
    q_ca: list = []
    r_ca: list = []
    for (qs, qe), (rs, re) in zip(aln.aligned[0], aln.aligned[1]):
        for qi, ri in zip(range(qs, qe), range(rs, re)):
            if qi < len(query_res) and ri < len(ref_res):
                q_ca.append(query_res[qi]["CA"])
                r_ca.append(ref_res[ri]["CA"])
    return q_ca, r_ca


def _download_reference_cif(pdb_id: str, dest_dir: str) -> str | None:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    dest_path = os.path.join(dest_dir, f"{pdb_id.upper()}.cif")
    try:
        with httpx.Client(timeout=_TIMEOUT) as client:
            resp = client.get(url)
            resp.raise_for_status()
            with open(dest_path, "wb") as f:
                f.write(resp.content)
        return dest_path
    except (httpx.HTTPStatusError, httpx.ConnectError, httpx.TimeoutException) as exc:
        logger.warning(f"Could not download reference PDB {pdb_id}: {exc}")
        return None


def _compute_rmsd(
    query_struct_path: str,
    reference_pdb_id: str,
    query_entry_id: str,
    query_aa: str,
) -> float | None:
    """Compute Cα RMSD using sequence-guided alignment.

    Fixes three bugs in the original implementation:
    1. Multi-chain reference PDBs: selects the chain whose sequence best
       matches query_aa rather than concatenating all chains.
    2. Length mismatch / signal peptide offsets: uses pairwise alignment
       to match only corresponding residues before superimposing.
    3. Wrong atom selection: only paired Cα atoms from the alignment are
       passed to Superimposer, never a raw positional slice.
    """
    try:
        tmpdir = tempfile.mkdtemp(prefix="pdb_ref_")
        ref_path = _download_reference_cif(reference_pdb_id, tmpdir)
        if ref_path is None:
            return None

        query_struct = _load_structure(query_struct_path, query_entry_id)
        ref_struct = _load_structure(ref_path, reference_pdb_id)

        query_chains = _get_chain_residues(query_struct)
        ref_chains = _get_chain_residues(ref_struct)

        if not query_chains or not ref_chains:
            return None

        query_res = next(iter(query_chains.values()))
        query_seq = _residues_to_seq(query_res)

        _, ref_res = _best_ref_chain(query_aa, ref_chains)
        if ref_res is None:
            return None
        ref_seq = _residues_to_seq(ref_res)

        q_ca, r_ca = _aligned_ca_pairs(query_seq, ref_seq, query_res, ref_res)

        if len(q_ca) < 3:
            logger.warning(
                f"Too few aligned Cα atoms for RMSD ({len(q_ca)}) — "
                f"query={len(query_res)} res, ref={len(ref_res)} res"
            )
            return None

        sup = Superimposer()
        sup.set_atoms(r_ca, q_ca)
        return round(float(sup.rms), 3)

    except Exception:
        logger.warning(
            f"RMSD computation failed for {query_entry_id} vs {reference_pdb_id}",
            exc_info=True,
        )
        return None


async def _fetch_predictions(uniprot_id: str) -> list[dict]:
    import asyncio

    url = f"{ALPHAFOLD_API}/prediction/{uniprot_id}"
    last_error: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            async with httpx.AsyncClient(timeout=_TIMEOUT) as client:
                resp = await client.get(url)
                resp.raise_for_status()
                return resp.json()
        except (httpx.TimeoutException, httpx.ConnectError, httpx.HTTPStatusError) as exc:
            last_error = exc
            if isinstance(exc, httpx.HTTPStatusError) and exc.response.status_code < 500:
                raise
            wait = _BACKOFF_BASE * (2**attempt)
            logger.warning(
                f"AlphaFold API attempt {attempt + 1}/{_MAX_RETRIES} failed for "
                f"{uniprot_id}: {exc}. Retrying in {wait:.1f}s..."
            )
            await asyncio.sleep(wait)

    raise RuntimeError(
        f"Failed to fetch AlphaFold predictions for {uniprot_id} "
        f"after {_MAX_RETRIES} attempts"
    ) from last_error


async def _download_cif(cif_url: str, dest_path: str) -> None:
    async with httpx.AsyncClient(timeout=_TIMEOUT) as client:
        resp = await client.get(cif_url)
        resp.raise_for_status()
        with open(dest_path, "wb") as f:
            f.write(resp.content)


async def predict_structure(
    chain_id: str,
    uniprot_id: str,
    aa_sequence: str,
    reference_pdb: str | None = None,
) -> StructureResult:
    """Fetch a pre-computed AlphaFold DB structure for uniprot_id.

    aa_sequence is the actual mature chain sequence produced by the pipeline.
    It is used to select the correct chain in a multimer reference PDB and to
    align residues before Cα superimposition.
    """
    predictions = await _fetch_predictions(uniprot_id)

    if not predictions:
        logger.warning(f"No AlphaFold predictions found for {uniprot_id}")
        return StructureResult(
            id=chain_id,
            plddt_mean=0.0,
            rmsd_to_ref=None,
            perplexity=None,
            structure_file="",
            method="afdb_fallback",
        )

    pred = predictions[0]
    plddt_mean = pred.get("globalMetricValue", 0.0)
    cif_url = pred.get("cifUrl", "")
    entry_id = pred.get("entryId", f"AF-{uniprot_id}-F1")

    cif_path = ""
    rmsd = None

    if cif_url:
        output_dir = os.path.join("output", "structures")
        os.makedirs(output_dir, exist_ok=True)
        cif_path = os.path.join(output_dir, f"{entry_id}.cif")

        if not os.path.exists(cif_path):
            await _download_cif(cif_url, cif_path)
            logger.info(f"Downloaded AlphaFold structure to {cif_path}")

        if reference_pdb:
            rmsd = _compute_rmsd(cif_path, reference_pdb, entry_id, aa_sequence)

    logger.info(
        f"AFDB prediction for {uniprot_id} (chain '{chain_id}'): "
        f"pLDDT={plddt_mean:.1f}, RMSD={rmsd}"
    )

    return StructureResult(
        id=chain_id,
        plddt_mean=plddt_mean,
        rmsd_to_ref=rmsd,
        perplexity=None,
        structure_file=cif_path,
        method="afdb",
    )