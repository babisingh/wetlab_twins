from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class BiosafetyResult:
    safe: bool
    matched_agent: str | None       # human-readable name of the flagged agent
    match_type: Literal["uniprot_id", "keyword", "aa_signature"] | None
    reason: str | None              # one-sentence explanation shown to the user
