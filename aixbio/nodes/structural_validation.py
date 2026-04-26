from __future__ import annotations

import asyncio
import logging

from aixbio.models.structure import StructureReport, StructureResult
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.alphafold import predict_structure as afdb_predict
from aixbio.tools.esmfold import fold_sequence as esmfold_fold

logger = logging.getLogger(__name__)

REFERENCE_PDBS = {
    "Insulin": "4INS",
    "hGH": "1HGU",
    "EPO": "1BUY",
    "Interferon": "1AU1",
    "Chymosin": "4CMS",
    "tPA": "1A5H",
}

_STRUCTURES_DIR = "output/structures"


async def _validate_chain(
    chain,
    ref_pdb: str | None,
    uniprot_id: str,
) -> StructureResult:
    """Fold chain.aa_sequence and optionally compute RMSD vs reference.

    Priority:
    1. ESMFold — folds the actual mature chain sequence the pipeline produced.
    2. AFDB fallback — used only when ESMFold is unavailable; pLDDT reflects
       the precursor, not the mature chain (logged as a warning).
    """
    result, struct_path = await esmfold_fold(
        chain_id=chain.id,
        aa_sequence=chain.aa_sequence,
        output_dir=_STRUCTURES_DIR,
    )

    if result.method == "esmfold_failed":
        logger.warning(
            f"ESMFold unavailable for chain {chain.id}; falling back to AFDB "
            f"(pLDDT will reflect precursor {uniprot_id}, not the mature chain)"
        )
        return await afdb_predict(
            chain_id=chain.id,
            uniprot_id=uniprot_id,
            aa_sequence=chain.aa_sequence,
            reference_pdb=ref_pdb,
        )

    rmsd: float | None = None
    if ref_pdb and struct_path:
        from aixbio.tools.alphafold import _compute_rmsd
        rmsd = _compute_rmsd(
            query_struct_path=struct_path,
            reference_pdb_id=ref_pdb,
            query_entry_id=f"ESM-{chain.id}",
            query_aa=chain.aa_sequence,
        )

    return StructureResult(
        id=result.id,
        plddt_mean=result.plddt_mean,
        rmsd_to_ref=rmsd,
        perplexity=result.perplexity,
        structure_file=result.structure_file,
        method=result.method,
    )


def structural_validation(state: PipelineState) -> dict:
    protein = state["protein_record"]
    if protein is None:
        return {"structure_report": None}

    ref_pdb: str | None = None
    for name, pdb_id in REFERENCE_PDBS.items():
        if name.lower() in protein.name.lower():
            ref_pdb = pdb_id
            break

    results = []
    for chain in protein.chains:
        result = asyncio.run(
            _validate_chain(
                chain=chain,
                ref_pdb=ref_pdb,
                uniprot_id=protein.uniprot_id,
            )
        )
        results.append(result)

    return {"structure_report": StructureReport(chains=tuple(results))}