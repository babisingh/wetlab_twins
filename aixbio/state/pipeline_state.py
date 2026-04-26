from __future__ import annotations

from typing import Annotated, Literal, TypedDict

from aixbio.models.audit import AgentDecision
from aixbio.models.biosafety import BiosafetyResult
from aixbio.models.host import HostRecommendation
from aixbio.models.protein import ProteinRecord
from aixbio.models.remediation import RemediationAction
from aixbio.models.structure import StructureReport
from aixbio.models.validation import ChainValidation, CheckResult, ValidationReport



def merge_chain_validations(
    existing: tuple[ChainValidation, ...], new: tuple[ChainValidation, ...]
) -> tuple[ChainValidation, ...]:
    return existing + new


def append_log(existing: list, new: list) -> list:
    """Reducer for LangGraph Annotated fields.

    Each node returns only *its own* new items (e.g. remediation actions from
    this round). The reducer automatically accumulates them by concatenation.
    This means nodes should NOT return the full history — just the delta.
    """
    return existing + new


class ChainProcessingResult(TypedDict):
    chain_id: str
    optimized_dna_sequence: str
    cassette_dna: str
    genbank_file: str
    insert_size: int
    validation_passed: bool
    checks: tuple[CheckResult, ...]
    remediation_rounds: int
    remediation_history: tuple[RemediationAction, ...]
    # Solubility / inclusion-body prediction (from chain subgraph step 1)
    solubility_score: float | None
    solubility_reasoning: str
    disulfide_risk: bool
    # Status distinguishes pass, fail, and max-retries-exceeded outcomes
    status: Literal[
        "passed", "failed", "max_retries_exceeded",
        "host_incompatible", "escalation_failed",
    ]


class PipelineState(TypedDict):
    # Inputs
    compound_id: str
    host_organism: str
    avoid_sites: tuple[str, ...]
    tag_type: str
    protease_site: str
    vector: str
    cloning_sites: tuple[str, ...]
    run_structural_validation: bool
    max_remediation_attempts: int
    enable_escalation: bool

    # Step 1 output
    protein_record: ProteinRecord | None
    chain_extraction_reasoning: str

    # Biosafety screen (runs after sequence_retrieval, before host selection)
    biosafety_result: BiosafetyResult | None

    # Host recommendation (runs after sequence_retrieval, before fan-out)
    host_recommendation: HostRecommendation | None

    # Per-chain results (fan-out merge target)
    chain_results: Annotated[list[ChainProcessingResult], append_log]

    # Step 5 aggregate
    validation_report: ValidationReport | None

    # Step 6 output (optional)
    structure_report: StructureReport | None

    # Protocol generation output (optional, requires --protocol flag)
    run_protocol_generation: bool
    protocol: str | None

    # Control flow
    pipeline_status: Literal[
        "running", "completed", "failed", "halted", "awaiting_human",
        "biosafety_rejected",
    ]

    # Auditability
    decision_log: Annotated[list[AgentDecision], append_log]
    warnings: Annotated[list[str], append_log]

    # Human-in-the-loop
    human_approval_chain_extraction: bool | None
    human_approval_plasmid: bool | None
