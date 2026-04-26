from __future__ import annotations

import re

from aixbio.models.host import HostFeatures, HostRecommendation


def recommend_host(chains, protein_name: str = "") -> HostRecommendation:
    """Deterministic host-organism recommendation from amino acid sequence features.

    Decision priority:
    1. N-glycosylation sequons (N-X-S/T, X != P) → mammalian (CHO) required.
    2. Very long protein (> 500 aa total) → CHO preferred for proper folding.
    3. Short disulfide-rich protein (< 100 aa, >= 4 Cys, no glyco) → E. coli refolding.
    4. Moderate protein with many Cys (>= 4, no glyco) → P. pastoris secretion.
    5. Highly hydrophobic (GRAVY > 1.5) → Sf9 insect cells.
    6. Some disulfide risk (>= 2 Cys, no glyco) → E. coli with refolding.
    7. Default: E. coli (simplest, fastest, lowest cost).
    """
    features = _extract_features(chains)

    # Rule 1: N-glycosylation requires mammalian expression
    if features.n_glycosylation_sites > 0:
        return HostRecommendation(
            primary_host="CHO",
            confidence="high",
            reasoning=(
                f"{features.n_glycosylation_sites} N-glycosylation sequon(s) detected. "
                "Mammalian CHO cells are required for glycoprotein bioactivity; "
                "E. coli lacks the glycosylation machinery."
            ),
            features=features,
            alternative_hosts=("HEK293",),
            caveats=("Glycan heterogeneity may affect batch-to-batch consistency.",),
        )

    # Rule 2: Large protein — mammalian chaperones and secretory pathway improve folding
    if features.total_chain_length > 500:
        return HostRecommendation(
            primary_host="CHO",
            confidence="medium",
            reasoning=(
                f"Total chain length {features.total_chain_length} aa exceeds 500. "
                "CHO cells provide the folding machinery needed for large proteins. "
                "E. coli often yields insoluble aggregates at this size."
            ),
            features=features,
            alternative_hosts=("HEK293", "Sf9"),
            caveats=("High cell culture cost; consider scale-up planning early.",),
        )

    # Rule 3: Small (< 100 aa), disulfide-rich — E. coli inclusion body route is established
    if features.total_chain_length < 100 and features.total_cysteine_count >= 4:
        return HostRecommendation(
            primary_host="Escherichia coli",
            confidence="medium",
            reasoning=(
                f"Short protein ({features.total_chain_length} aa) with "
                f"{features.total_cysteine_count} cysteines. "
                "E. coli expression as inclusion bodies followed by oxidative refolding "
                "is the established route (e.g. recombinant insulin)."
            ),
            features=features,
            alternative_hosts=("Pichia pastoris",),
            caveats=(
                "Denaturing lysis (urea/guanidinium) required.",
                "Refolding yield is protein-dependent; optimise glutathione redox ratio.",
            ),
        )

    # Rule 4: Moderate size, many Cys — yeast secretory pathway handles disulfides better
    if features.total_cysteine_count >= 4:
        return HostRecommendation(
            primary_host="Pichia pastoris",
            confidence="medium",
            reasoning=(
                f"{features.total_cysteine_count} cysteines detected. "
                "P. pastoris secretory pathway provides an oxidising environment for "
                "disulfide bond formation without inclusion body refolding."
            ),
            features=features,
            alternative_hosts=("Escherichia coli",),
            caveats=(
                "Requires methanol-inducible AOX1 promoter; scale-up is more complex than E. coli.",
                "E. coli refolding is a viable fallback for small proteins.",
            ),
        )

    # Rule 5: Highly hydrophobic — Sf9 insect cells tolerate membrane proteins better
    if features.max_gravy > 1.5:
        return HostRecommendation(
            primary_host="Sf9",
            confidence="medium",
            reasoning=(
                f"High GRAVY score ({features.max_gravy:.2f}) suggests a hydrophobic or "
                "membrane-associated protein. Sf9 insect cells (baculovirus system) "
                "provide lipid bilayers and post-translational processing."
            ),
            features=features,
            alternative_hosts=("CHO",),
            caveats=("Baculovirus construct design adds lead time.",),
        )

    # Rule 6: Some disulfide risk — E. coli refolding or P. pastoris
    if features.disulfide_risk:
        return HostRecommendation(
            primary_host="Escherichia coli",
            confidence="medium",
            reasoning=(
                f"{features.total_cysteine_count} cysteine(s) present. "
                "E. coli expression is feasible with oxidative refolding; "
                "P. pastoris is a lower-risk alternative."
            ),
            features=features,
            alternative_hosts=("Pichia pastoris",),
            caveats=(
                "Disulfide bonds will not form correctly in the E. coli cytoplasm; "
                "use refolding buffer with reduced/oxidised glutathione.",
            ),
        )

    # Rule 7: Default — straightforward E. coli expression
    return HostRecommendation(
        primary_host="Escherichia coli",
        confidence="high",
        reasoning=(
            "No glycosylation sites, low cysteine count, and moderate hydrophobicity. "
            "E. coli BL21(DE3) + pET vector system is the simplest and most cost-effective choice."
        ),
        features=features,
        alternative_hosts=(),
        caveats=(),
    )


def _extract_features(chains) -> HostFeatures:
    total_length = sum(c.length for c in chains)
    total_cys = sum(c.aa_sequence.upper().count("C") for c in chains)
    disulfide_risk = any(c.aa_sequence.upper().count("C") >= 2 for c in chains)
    n_glyco = sum(_count_nglycosylation(c.aa_sequence) for c in chains)
    max_gravy = max((_gravy(c.aa_sequence) for c in chains), default=0.0)
    return HostFeatures(
        total_chain_length=total_length,
        n_glycosylation_sites=n_glyco,
        total_cysteine_count=total_cys,
        disulfide_risk=disulfide_risk,
        max_gravy=round(max_gravy, 3),
        multi_chain=len(list(chains)) > 1,
    )


def _count_nglycosylation(aa_seq: str) -> int:
    return len(re.findall(r"N[^P][ST]", aa_seq.upper()))


def _gravy(aa_seq: str) -> float:
    # Kyte-Doolittle hydropathy scale
    KD = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
        "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    }
    seq = aa_seq.upper()
    scores = [KD.get(aa, 0.0) for aa in seq]
    return sum(scores) / len(scores) if scores else 0.0
