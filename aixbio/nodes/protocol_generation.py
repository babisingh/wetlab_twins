from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

from langchain_core.messages import HumanMessage, SystemMessage

from aixbio.config import (
    ChatOpenRouter,
    LLM_MAX_TOKENS,
    LLM_MODEL,
    OPENROUTER_API_KEY,
    OPENROUTER_BASE_URL,
)
from aixbio.models.audit import AgentDecision
from aixbio.prompts.protocol import PROTOCOL_SYSTEM
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.pubmed import search_expression_literature

logger = logging.getLogger(__name__)


def protocol_generation(state: PipelineState) -> dict:
    protein = state["protein_record"]
    if protein is None:
        return {}

    logger.info("Fetching PubMed literature for %s (%s)", protein.name, protein.uniprot_id)
    literature = search_expression_literature(protein.name, protein.uniprot_id, max_results=5)

    context = _build_context(state, literature)

    llm = ChatOpenRouter(
        model=LLM_MODEL,
        max_tokens=LLM_MAX_TOKENS,
        openai_api_key=OPENROUTER_API_KEY,
        openai_api_base=OPENROUTER_BASE_URL,
        temperature=0.2,
    )

    protocol_text = llm.invoke([
        SystemMessage(content=PROTOCOL_SYSTEM),
        HumanMessage(content=json.dumps(context, default=str)),
    ]).content

    pmid_count = len(literature)
    audit = AgentDecision(
        node="protocol_generation",
        reasoning=(
            f"Generated SOP from {len(state.get('chain_results', []))} chain(s) + "
            f"{pmid_count} PubMed article(s)"
        ),
        action="generate_protocol",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"{protein.name} ({protein.uniprot_id}); {pmid_count} literature refs",
        output_summary=f"Protocol: {len(protocol_text)} chars",
    )

    return {
        "protocol": protocol_text,
        "decision_log": [audit],
    }


def _build_context(state: PipelineState, literature: list[dict]) -> dict:
    protein = state["protein_record"]
    chain_results = state.get("chain_results", [])

    sequence_metrics = []
    for cr in chain_results:
        sequence_metrics.append({
            "chain_id": cr["chain_id"],
            "insert_size_bp": cr["insert_size"],
            "validation_passed": cr["validation_passed"],
            "remediation_rounds": cr["remediation_rounds"],
            "solubility_score": cr.get("solubility_score"),
            "predicted_soluble": (
                cr["solubility_score"] >= 0.45
                if cr.get("solubility_score") is not None
                else None
            ),
            "disulfide_risk": cr.get("disulfide_risk", False),
            "cysteine_count": _count_cys(cr.get("chain_id", ""), protein),
            "tag_recommendation": cr.get("solubility_reasoning", ""),
            "checks": [
                {"name": c.name, "passed": c.passed, "value": c.value}
                for c in cr.get("checks", ())
            ],
        })

    return {
        "protein": {
            "name": protein.name,
            "uniprot_id": protein.uniprot_id,
            "chain_count": len(protein.chains),
            "chains": [
                {"id": c.id, "length": c.length, "aa_sequence": c.aa_sequence}
                for c in protein.chains
            ],
        },
        "expression_system": {
            "host": state.get("host_organism", "Escherichia coli"),
            "vector": state.get("vector", "pET-28a(+)"),
            "tag": state.get("tag_type", "6xHis"),
            "protease_site": state.get("protease_site", "Enterokinase"),
            "cloning_sites": list(state.get("cloning_sites", ("BamHI", "XhoI"))),
        },
        "sequence_metrics": sequence_metrics,
        "overall_validation_passed": state.get("validation_report") is not None
        and state["validation_report"].all_passed,
        "literature": literature,
        "warnings": state.get("warnings", []),
    }


def _count_cys(chain_id: str, protein) -> int:
    for chain in protein.chains:
        if chain.id == chain_id:
            return chain.aa_sequence.upper().count("C")
    return 0
