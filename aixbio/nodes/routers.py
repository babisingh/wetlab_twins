from __future__ import annotations

from typing import Literal

from langgraph.types import Send

from aixbio.state.chain_state import ChainSubgraphState
from aixbio.state.pipeline_state import PipelineState


def _validation_core(state: ChainSubgraphState) -> str:
    validation = state["chain_validation"]
    if validation is None:
        return "halt_pipeline"
    if validation.passed:
        return "package_result"
    for check in validation.checks:
        if check.name == "back_translation" and not check.passed:
            return "halt_pipeline"
    if state["remediation_attempt"] >= state["max_remediation_attempts"]:
        if state.get("enable_escalation") and not state.get("escalation_used"):
            return "escalation_agent"
        return "package_result_failed"
    return "remediation_agent"


def validation_router(state: ChainSubgraphState) -> str:
    return _validation_core(state)


def revalidation_router(state: ChainSubgraphState) -> str:
    return _validation_core(state)


def escalation_router(state: ChainSubgraphState) -> Literal[
    "apply_fixes", "codon_optimization",
    "package_result_escalated", "halt_pipeline",
]:
    d = state.get("escalation_decision")
    if d is None or d.kind == "give_up":
        return "halt_pipeline"
    if d.kind == "incompatible":
        return "package_result_escalated"
    if d.kind == "apply_plan":
        return "apply_fixes" if d.actions else "package_result_escalated"
    if d.kind == "change_strategy":
        return "codon_optimization"
    return "halt_pipeline"


def fan_out_to_chains(state: PipelineState) -> list[Send]:
    protein = state["protein_record"]
    if protein is None:
        return []

    sends = []
    for chain in protein.chains:
        chain_state: ChainSubgraphState = {
            "chain": chain,
            "host_organism": state["host_organism"],
            "avoid_sites": state["avoid_sites"],
            "tag_type": state["tag_type"],
            "protease_site": state["protease_site"],
            "vector": state["vector"],
            "cloning_sites": state["cloning_sites"],
            "protein_record": protein,
            "solubility_result": None,
            "optimized_dna": None,
            "cassette": None,
            "plasmid": None,
            "chain_validation": None,
            "remediation_attempt": 0,
            "max_remediation_attempts": state.get("max_remediation_attempts", 3),
            "failed_checks": (),
            "remediation_plan": None,
            "remediation_history": [],
            "enable_escalation": state.get("enable_escalation", False),
            "escalation_used": False,
            "escalation_decision": None,
            "decision_log": [],
            "warnings": [],
        }
        sends.append(Send("chain_processing", chain_state))
    return sends


def biosafety_router(state: PipelineState) -> str:
    if state.get("pipeline_status") == "biosafety_rejected":
        return "__end__"
    return "host_selection"


def structural_router(
    state: PipelineState,
) -> Literal["structural_validation", "protocol_generation", "__end__"]:
    if state.get("run_structural_validation"):
        return "structural_validation"
    if state.get("run_protocol_generation"):
        return "protocol_generation"
    return "__end__"


def post_structural_router(
    state: PipelineState,
) -> Literal["protocol_generation", "__end__"]:
    if state.get("run_protocol_generation"):
        return "protocol_generation"
    return "__end__"
