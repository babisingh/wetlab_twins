from __future__ import annotations

from aixbio.tools.biosafety import screen_protein
from aixbio.state.pipeline_state import PipelineState


def biosafety_screen(state: PipelineState) -> dict:
    """Layer 2 + 3 biosafety check, runs after sequence retrieval.

    If the protein name or any chain sequence matches a select-agent blocklist,
    the pipeline status is set to 'biosafety_rejected' so the router can halt
    immediately without writing any output artifacts.
    """
    protein = state.get("protein_record")
    if protein is None:
        return {"pipeline_status": "biosafety_rejected",
                "warnings": ["Biosafety screen: no protein record returned by retrieval step."]}

    chains_aa = [c.aa_sequence for c in protein.chains]
    result = screen_protein(protein.name, chains_aa)

    if not result.safe:
        return {
            "pipeline_status": "biosafety_rejected",
            "biosafety_result": result,
            "warnings": [
                f"BIOSAFETY BLOCK — {result.reason}"
            ],
        }

    return {"biosafety_result": result}
