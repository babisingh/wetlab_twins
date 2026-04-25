from __future__ import annotations

from langgraph.checkpoint.memory import MemorySaver
from langgraph.checkpoint.serde.jsonplus import JsonPlusSerializer
from langgraph.graph import END, START, StateGraph

from aixbio.graph.chain_subgraph import compile_chain_subgraph
from aixbio.nodes.human_checkpoints import (
    human_checkpoint_chains,
    human_checkpoint_plasmid,
)
from aixbio.nodes.merge_results import merge_all_chain_results
from aixbio.nodes.routers import fan_out_to_chains, structural_router
from aixbio.nodes.sequence_retrieval import sequence_retrieval_agent
from aixbio.nodes.structural_validation import structural_validation
from aixbio.state.pipeline_state import PipelineState


def build_main_graph() -> StateGraph:
    g = StateGraph(PipelineState)

    # Step 1: Sequence Retrieval (LLM agent)
    g.add_node("sequence_retrieval", sequence_retrieval_agent)

    # Human checkpoint: review chain extraction
    g.add_node("human_checkpoint_chains", human_checkpoint_chains)

    # Per-chain subgraph (fan-out target)
    g.add_node("chain_processing", compile_chain_subgraph())

    # Merge fan-out results
    g.add_node("merge_results", merge_all_chain_results)

    # Human checkpoint: review plasmids
    g.add_node("human_checkpoint_plasmid", human_checkpoint_plasmid)

    # Step 6: Structural Validation (optional)
    g.add_node("structural_validation", structural_validation)

    # Edges
    g.add_edge(START, "sequence_retrieval")
    g.add_edge("sequence_retrieval", "human_checkpoint_chains")

    # Fan-out: human_checkpoint_chains -> N x chain_processing
    g.add_conditional_edges(
        "human_checkpoint_chains",
        fan_out_to_chains,
        ["chain_processing"],
    )

    # All chain_processing branches merge into merge_results
    g.add_edge("chain_processing", "merge_results")
    g.add_edge("merge_results", "human_checkpoint_plasmid")

    # Conditional: run structural validation or finish
    g.add_conditional_edges(
        "human_checkpoint_plasmid",
        structural_router,
        {
            "structural_validation": "structural_validation",
            "__end__": END,
        },
    )

    g.add_edge("structural_validation", END)

    return g


_ALLOWED_MSGPACK_MODULES = [
    ("aixbio.models.protein", "Chain"),
    ("aixbio.models.protein", "ProteinRecord"),
    ("aixbio.models.validation", "CheckResult"),
    ("aixbio.models.validation", "ChainValidation"),
    ("aixbio.models.validation", "ValidationReport"),
    ("aixbio.models.dna", "DNAChain"),
    ("aixbio.models.dna", "OptimizedDNA"),
    ("aixbio.models.dna", "CassetteElement"),
    ("aixbio.models.dna", "CassetteChain"),
    ("aixbio.models.dna", "CassetteDNA"),
    ("aixbio.models.plasmid", "PlasmidChain"),
    ("aixbio.models.plasmid", "PlasmidRecord"),
    ("aixbio.models.audit", "AgentDecision"),
    ("aixbio.models.remediation", "RemediationAction"),
    ("aixbio.models.remediation", "PlannedFix"),
    ("aixbio.models.remediation", "RemediationPlan"),
    ("aixbio.models.structure", "StructureResult"),
    ("aixbio.models.structure", "StructureReport"),
    ("aixbio.models.escalation", "EscalationApplyPlan"),
    ("aixbio.models.escalation", "EscalationIncompatible"),
    ("aixbio.models.escalation", "EscalationChangeStrategy"),
    ("aixbio.models.escalation", "EscalationGiveUp"),
]


def compile_pipeline(checkpointer=None):
    if checkpointer is None:
        serde = JsonPlusSerializer(allowed_msgpack_modules=_ALLOWED_MSGPACK_MODULES)
        checkpointer = MemorySaver(serde=serde)
    return build_main_graph().compile(checkpointer=checkpointer)
