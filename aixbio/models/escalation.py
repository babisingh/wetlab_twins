from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from aixbio.models.remediation import PlannedFix


@dataclass(frozen=True)
class EscalationApplyPlan:
    kind: Literal["apply_plan"] = "apply_plan"
    actions: tuple[PlannedFix, ...] = ()
    reasoning: str = ""
    diagnosis: str = ""


@dataclass(frozen=True)
class EscalationIncompatible:
    kind: Literal["incompatible"] = "incompatible"
    reason: str = ""
    suggested_action: str = ""


@dataclass(frozen=True)
class EscalationChangeStrategy:
    kind: Literal["change_strategy"] = "change_strategy"
    field: Literal["tag_type", "protease_site", "vector", "cloning_sites"] = "tag_type"
    new_value: str = ""
    reason: str = ""


@dataclass(frozen=True)
class EscalationGiveUp:
    kind: Literal["give_up"] = "give_up"
    diagnosis: str = ""


EscalationDecision = (
    EscalationApplyPlan | EscalationIncompatible
    | EscalationChangeStrategy | EscalationGiveUp
)
