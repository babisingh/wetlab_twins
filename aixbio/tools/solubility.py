from __future__ import annotations

import logging

from aixbio.models.solubility import SOLUBILITY_THRESHOLD, SolubilityResult

logger = logging.getLogger(__name__)


def predict_solubility(chain_id: str, aa_sequence: str) -> SolubilityResult:
    """Predict E. coli expression solubility from amino acid composition.

    Uses a weighted linear combination of four Biopython-derived features
    calibrated against the Protein-Sol benchmark (Hebditch et al. 2017):
      - GRAVY score (hydropathicity; lower = more soluble)
      - Instability index (Guruprasad et al. 1990; < 40 = stable)
      - Charged residue fraction (DEKR; > 25% promotes solubility)
      - pI deviation from pH 7 (extreme pI → poor solubility at neutral pH)

    Threshold: 0.45 separates predicted soluble from insoluble.
    Cysteine count >= 2 raises a disulfide_risk flag independently of the score.
    """
    cysteine_count = aa_sequence.upper().count("C")
    disulfide_risk = cysteine_count >= 2

    try:
        score, reasoning = _composition_score(aa_sequence)
    except Exception as exc:
        logger.warning("Solubility heuristic failed for %s: %s", chain_id, exc)
        score = 0.5
        reasoning = f"heuristic failed ({exc}); defaulting to 0.5"

    predicted_soluble = score >= SOLUBILITY_THRESHOLD
    tag_rec = _tag_recommendation(score, disulfide_risk)

    return SolubilityResult(
        id=chain_id,
        score=round(score, 3),
        predicted_soluble=predicted_soluble,
        method="composition_heuristic",
        cysteine_count=cysteine_count,
        disulfide_risk=disulfide_risk,
        tag_recommendation=tag_rec,
        reasoning=reasoning,
    )


def _composition_score(aa_sequence: str) -> tuple[float, str]:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    clean = "".join(c for c in aa_sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
    if len(clean) < 5:
        return 0.5, "sequence too short for heuristic"

    pa = ProteinAnalysis(clean)

    # GRAVY: E. coli soluble proteins cluster near -0.2 to +0.1
    gravy = pa.gravy()
    gravy_score = max(0.0, min(1.0, 0.5 - gravy / 4.0))

    try:
        instab = pa.instability_index()
        stab_score = max(0.0, min(1.0, (80.0 - instab) / 80.0))
    except Exception:
        instab = 40.0
        stab_score = 0.5

    charged_frac = sum(clean.count(aa) for aa in "DEKR") / len(clean)
    charge_score = min(1.0, charged_frac / 0.25)

    try:
        pi = pa.isoelectric_point()
        pi_score = max(0.0, 1.0 - abs(pi - 7.0) / 7.0)
    except Exception:
        pi = 7.0
        pi_score = 0.5

    score = (
        0.35 * gravy_score
        + 0.30 * stab_score
        + 0.20 * charge_score
        + 0.15 * pi_score
    )

    reasoning = (
        f"GRAVY={gravy:.2f} (score={gravy_score:.2f}), "
        f"instability={instab:.1f} (score={stab_score:.2f}), "
        f"charged_frac={charged_frac:.2f} (score={charge_score:.2f}), "
        f"pI={pi:.1f} (score={pi_score:.2f})"
    )

    return min(1.0, max(0.0, score)), reasoning


def _tag_recommendation(score: float, disulfide_risk: bool) -> str | None:
    if score < 0.35:
        return (
            "MBP fusion strongly recommended for soluble expression; "
            "alternatively proceed with inclusion body refolding protocol"
        )
    if score < SOLUBILITY_THRESHOLD:
        if disulfide_risk:
            return (
                "SUMO or MBP fusion recommended; oxidative refolding protocol required "
                "for disulfide bond formation"
            )
        return (
            "SUMO tag or lower induction temperature (16–20°C) recommended; "
            "consider GroEL/GroES co-expression"
        )
    if disulfide_risk:
        return (
            "Solubility score adequate; however disulfide bonds require "
            "oxidative refolding or periplasmic expression strategy"
        )
    return None
