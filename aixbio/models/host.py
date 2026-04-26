from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HostFeatures:
    total_chain_length: int      # sum of all mature chain lengths
    n_glycosylation_sites: int   # N-X-S/T sequons (X != P), summed across all chains
    total_cysteine_count: int    # all Cys in all chains
    disulfide_risk: bool         # any chain has >= 2 Cys
    max_gravy: float             # highest GRAVY score among chains
    multi_chain: bool            # protein has more than one chain


@dataclass(frozen=True)
class HostRecommendation:
    primary_host: str            # e.g. "Escherichia coli"
    confidence: str              # "high" | "medium" | "low"
    reasoning: str
    features: HostFeatures
    alternative_hosts: tuple[str, ...]
    caveats: tuple[str, ...]
