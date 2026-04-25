from __future__ import annotations

from typing import Annotated, TypedDict

from aixbio.models.audit import AgentDecision
from aixbio.models.dna import CassetteChain, DNAChain
from aixbio.models.escalation import EscalationDecision
from aixbio.models.plasmid import PlasmidChain
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.models.remediation import RemediationAction, RemediationPlan
from aixbio.models.validation import ChainValidation, CheckResult
from aixbio.state.pipeline_state import ChainProcessingResult, append_log


class ChainSubgraphState(TypedDict):
    # Context from outer graph
    chain: Chain
    host_organism: str
    avoid_sites: tuple[str, ...]
    tag_type: str
    protease_site: str
    vector: str
    cloning_sites: tuple[str, ...]
    protein_record: ProteinRecord

    # Per-chain outputs (built up through the subgraph)
    optimized_dna: DNAChain | None
    cassette: CassetteChain | None
    plasmid: PlasmidChain | None
    chain_validation: ChainValidation | None

    # Remediation control
    remediation_attempt: int
    max_remediation_attempts: int
    failed_checks: tuple[CheckResult, ...]
    remediation_plan: RemediationPlan | None
    # Remediation history uses append_log reducer: each `apply_fixes` call
    # returns only *this round's* actions, and the reducer accumulates them.
    # Node code should NOT return the full history — just the new delta.
    remediation_history: Annotated[list[RemediationAction], append_log]

    # Escalation (opt-in, single-fire)
    enable_escalation: bool
    escalation_used: bool
    escalation_decision: EscalationDecision | None

    # Output for outer graph
    chain_results: Annotated[list[ChainProcessingResult], append_log]

    # Auditability
    decision_log: Annotated[list[AgentDecision], append_log]
    warnings: Annotated[list[str], append_log]
