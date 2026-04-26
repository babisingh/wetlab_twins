from __future__ import annotations

import logging

from aixbio.models.solubility import SOLUBILITY_THRESHOLD
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.solubility import predict_solubility

logger = logging.getLogger(__name__)


def solubility_prediction(state: ChainSubgraphState) -> dict:
    chain = state["chain"]
    result = predict_solubility(chain.id, chain.aa_sequence)

    warnings: list[str] = []

    if not result.predicted_soluble:
        warnings.append(
            f"[{chain.id}] Solubility score {result.score:.2f} < {SOLUBILITY_THRESHOLD} — "
            f"inclusion body formation likely in E. coli. "
            + (result.tag_recommendation or "")
        )

    if result.disulfide_risk:
        warnings.append(
            f"[{chain.id}] {result.cysteine_count} cysteine residues detected — "
            f"disulfide bond formation requires oxidative refolding or periplasmic expression."
        )

    logger.info(
        "Solubility [%s]: score=%.2f predicted_soluble=%s disulfide_risk=%s method=%s",
        chain.id,
        result.score,
        result.predicted_soluble,
        result.disulfide_risk,
        result.method,
    )

    out: dict = {"solubility_result": result}
    if warnings:
        out["warnings"] = warnings
    return out
