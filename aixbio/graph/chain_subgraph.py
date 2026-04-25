from __future__ import annotations

from langgraph.graph import END, START, StateGraph

from aixbio.nodes.cassette_assembly import cassette_assembly
from aixbio.nodes.codon_optimization import codon_optimization
from aixbio.nodes.escalation_agent import escalation_agent
from aixbio.nodes.merge_results import (
    halt_pipeline,
    package_result,
    package_result_escalated,
    package_result_failed,
)
from aixbio.nodes.plasmid_assembly import plasmid_assembly
from aixbio.nodes.remediation_agent import apply_fixes, remediation_agent
from aixbio.nodes.routers import escalation_router, revalidation_router, validation_router
from aixbio.nodes.sequence_validation import sequence_validation
from aixbio.state.chain_state import ChainSubgraphState

_VALIDATION_MAP = {
    "package_result": "package_result",
    "halt_pipeline": "halt_pipeline",
    "package_result_failed": "package_result_failed",
    "remediation_agent": "remediation_agent",
    "escalation_agent": "escalation_agent",
}

_ESCALATION_MAP = {
    "apply_fixes": "apply_fixes",
    "codon_optimization": "codon_optimization",
    "package_result_escalated": "package_result_escalated",
    "halt_pipeline": "halt_pipeline",
}


def build_chain_subgraph() -> StateGraph:
    g = StateGraph(ChainSubgraphState)

    # Steps 2-4: linear pipeline
    g.add_node("codon_optimization", codon_optimization)
    g.add_node("cassette_assembly", cassette_assembly)
    g.add_node("plasmid_assembly", plasmid_assembly)

    # Step 5: validation
    g.add_node("sequence_validation", sequence_validation)

    # Remediation loop
    g.add_node("remediation_agent", remediation_agent)
    g.add_node("apply_fixes", apply_fixes)
    g.add_node("reassemble_cassette", cassette_assembly)
    g.add_node("reassemble_plasmid", plasmid_assembly)
    g.add_node("revalidate", sequence_validation)

    # Escalation
    g.add_node("escalation_agent", escalation_agent)
    g.add_node("package_result_escalated", package_result_escalated)

    # Terminal nodes
    g.add_node("package_result", package_result)
    g.add_node("package_result_failed", package_result_failed)
    g.add_node("halt_pipeline", halt_pipeline)

    # Linear pipeline edges
    g.add_edge(START, "codon_optimization")
    g.add_edge("codon_optimization", "cassette_assembly")
    g.add_edge("cassette_assembly", "plasmid_assembly")
    g.add_edge("plasmid_assembly", "sequence_validation")

    # First validation routing
    g.add_conditional_edges(
        "sequence_validation",
        validation_router,
        _VALIDATION_MAP,
    )

    # Remediation cycle
    g.add_edge("remediation_agent", "apply_fixes")
    g.add_edge("apply_fixes", "reassemble_cassette")
    g.add_edge("reassemble_cassette", "reassemble_plasmid")
    g.add_edge("reassemble_plasmid", "revalidate")

    g.add_conditional_edges(
        "revalidate",
        revalidation_router,
        _VALIDATION_MAP,
    )

    # Escalation routing
    g.add_conditional_edges(
        "escalation_agent",
        escalation_router,
        _ESCALATION_MAP,
    )

    # Terminal edges
    g.add_edge("package_result", END)
    g.add_edge("package_result_failed", END)
    g.add_edge("package_result_escalated", END)
    g.add_edge("halt_pipeline", END)

    return g


def compile_chain_subgraph():
    return build_chain_subgraph().compile()
