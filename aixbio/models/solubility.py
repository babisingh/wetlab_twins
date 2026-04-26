from __future__ import annotations

from dataclasses import dataclass

SOLUBILITY_THRESHOLD = 0.45  # Protein-Sol calibrated boundary (Hebditch et al. 2017)


@dataclass(frozen=True)
class SolubilityResult:
    id: str
    score: float  # 0.0–1.0; >= 0.45 predicted soluble in E. coli
    predicted_soluble: bool
    method: str  # "composition_heuristic" or "protein_sol_api"
    cysteine_count: int
    disulfide_risk: bool  # True if cysteine_count >= 2
    tag_recommendation: str | None
    reasoning: str
