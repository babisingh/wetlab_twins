from __future__ import annotations

import logging
from datetime import datetime, timezone

from aixbio.models.audit import AgentDecision
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.host_selector import recommend_host

logger = logging.getLogger(__name__)


def host_selection(state: PipelineState) -> dict:
    protein = state["protein_record"]
    if protein is None:
        return {}

    rec = recommend_host(protein.chains, protein.name)

    logger.info(
        "Host selection [%s]: %s (confidence=%s)",
        protein.uniprot_id,
        rec.primary_host,
        rec.confidence,
    )

    audit = AgentDecision(
        node="host_selection",
        reasoning=rec.reasoning,
        action=f"recommend_host:{rec.primary_host}",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=(
            f"{protein.name} — {rec.features.total_chain_length} aa, "
            f"{rec.features.total_cysteine_count} Cys, "
            f"{rec.features.n_glycosylation_sites} N-glyco sites"
        ),
        output_summary=(
            f"{rec.primary_host} [{rec.confidence}]"
            + (f"; alts: {', '.join(rec.alternative_hosts)}" if rec.alternative_hosts else "")
        ),
    )

    out: dict = {
        "host_recommendation": rec,
        "decision_log": [audit],
        "warnings": [f"[host_selection] {c}" for c in rec.caveats],
    }
    return out
