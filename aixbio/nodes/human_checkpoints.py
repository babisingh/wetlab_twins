from __future__ import annotations

from langgraph.types import interrupt

from aixbio.state.pipeline_state import PipelineState


def human_checkpoint_chains(state: PipelineState) -> dict:
    protein = state["protein_record"]
    if protein is None:
        return {"pipeline_status": "halted"}

    reasoning = state.get("chain_extraction_reasoning", "")
    rec = state.get("host_recommendation")

    interrupt_data: dict = {
        "stage": "chain_extraction_review",
        "protein": protein.name,
        "uniprot_id": protein.uniprot_id,
        "chains_extracted": [
            {
                "id": c.id,
                "length": c.length,
                "sequence_preview": c.aa_sequence[:30] + ("..." if len(c.aa_sequence) > 30 else ""),
            }
            for c in protein.chains
        ],
        "agent_reasoning": reasoning,
        "question": "Approve chain extraction? Reply 'approve', 'reject', or provide corrections as JSON.",
    }

    if rec is not None:
        interrupt_data["host_recommendation"] = {
            "primary_host": rec.primary_host,
            "confidence": rec.confidence,
            "reasoning": rec.reasoning,
            "alternative_hosts": list(rec.alternative_hosts),
            "caveats": list(rec.caveats),
        }

    decision = interrupt(interrupt_data)

    if decision == "reject":
        return {"pipeline_status": "halted"}

    return {"human_approval_chain_extraction": True}


def human_checkpoint_plasmid(state: PipelineState) -> dict:
    chain_results = state.get("chain_results", [])

    decision = interrupt({
        "stage": "plasmid_review",
        "chain_results": [
            {
                "chain_id": cr["chain_id"],
                "validation_passed": cr["validation_passed"],
                "remediation_rounds": cr["remediation_rounds"],
                "insert_size": cr["insert_size"],
            }
            for cr in chain_results
        ],
        "question": "Approve plasmids? Reply 'approve' or 'reject'.",
    })

    if decision == "reject":
        return {"pipeline_status": "halted"}

    return {"human_approval_plasmid": True}
