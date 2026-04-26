from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StructureResult:
    id: str
    plddt_mean: float
    rmsd_to_ref: float | None
    perplexity: float | None
    structure_file: str
    # "esmfold" | "afdb" | "afdb_fallback" | "esmfold_failed" | "unknown"
    method: str = "unknown"


@dataclass(frozen=True)
class StructureReport:
    chains: tuple[StructureResult, ...]