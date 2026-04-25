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
from aixbio.models.escalation import (
    EscalationApplyPlan,
    EscalationChangeStrategy,
    EscalationGiveUp,
    EscalationIncompatible,
)
from aixbio.models.remediation import PlannedFix, RemediationPlan
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.codon_tables import split_codons, translate_codon

logger = logging.getLogger(__name__)


ESCALATION_SYSTEM = """\
You are an escalation reviewer for a deterministic codon-optimization pipeline. \
The pipeline tried N rounds of synonymous-codon hill-climbing and failed to \
satisfy validation. Choose ONE of four outcomes and return strict JSON.

You may NOT invent DNA. If you propose codon changes, every replacement MUST \
be a synonymous codon for the AA at that position; the harness rejects any \
non-synonymous swap.

Outcomes:
- apply_plan: deterministic loop is stuck in a local minimum; propose a \
small set of swaps the loop's greedy ordering would not have tried.
- incompatible: the compound is fundamentally incompatible with the host \
(e.g. requires post-translational modifications E. coli cannot do).
- change_strategy: a non-sequence config change would unblock (different \
tag, protease, vector, or cloning sites).
- give_up: failure indicates a pipeline bug or a check we cannot fix.

Output JSON shape:
{"kind": "apply_plan", "actions": [{"check_name": "...", "strategy": "...", \
"target_positions": [int], "replacement_codons": ["..."]}], \
"reasoning": "...", "diagnosis": "..."}
{"kind": "incompatible", "reason": "...", "suggested_action": "..."}
{"kind": "change_strategy", "field": "tag_type|protease_site|vector|cloning_sites", \
"new_value": "...", "reason": "..."}
{"kind": "give_up", "diagnosis": "..."}"""


def escalation_agent(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    failed = state["failed_checks"]
    history = state.get("remediation_history", [])
    chain = state["chain"]

    payload = {
        "compound_uniprot": state["protein_record"].uniprot_id,
        "host": state["host_organism"],
        "chain_id": chain.id,
        "aa_sequence": chain.aa_sequence,
        "dna_codons": split_codons(dna_chain.dna_sequence),
        "cai": dna_chain.cai_score,
        "gc": dna_chain.gc_content,
        "failed_checks": [
            {"name": c.name, "value": c.value, "threshold": c.threshold}
            for c in failed
        ],
        "remediation_attempts": state["remediation_attempt"],
        "remediation_history": [
            {
                "check": a.check_name,
                "fix": a.fix_type,
                "pos": list(a.positions_affected),
                "before": list(a.codons_before),
                "after": list(a.codons_after),
            }
            for a in history
        ],
        "tag_type": state["tag_type"],
        "protease_site": state["protease_site"],
        "vector": state["vector"],
    }

    llm = ChatOpenRouter(
        model=LLM_MODEL,
        max_tokens=LLM_MAX_TOKENS,
        openai_api_key=OPENROUTER_API_KEY,
        openai_api_base=OPENROUTER_BASE_URL,
        temperature=0,
    )
    raw = llm.invoke([
        SystemMessage(content=ESCALATION_SYSTEM),
        HumanMessage(content=json.dumps(payload)),
    ]).content

    decision = _parse_and_validate(raw, dna_chain)

    audit = AgentDecision(
        node="escalation_agent",
        reasoning=getattr(decision, "diagnosis", "") or getattr(decision, "reason", ""),
        action=decision.kind,
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"{chain.id}: {len(failed)} failed, {len(history)} prior fixes",
        output_summary=_summarize(decision),
    )

    out: dict = {
        "decision_log": [audit],
        "escalation_decision": decision,
        "escalation_used": True,
    }

    if isinstance(decision, EscalationApplyPlan):
        out["remediation_plan"] = RemediationPlan(
            actions=decision.actions,
            reasoning=f"[escalation] {decision.reasoning}",
            priority_order=tuple({a.check_name for a in decision.actions}),
        )
    elif isinstance(decision, EscalationChangeStrategy):
        out[decision.field] = decision.new_value
        out["warnings"] = [
            f"escalation changed {decision.field} -> {decision.new_value}: {decision.reason}"
        ]
        out["remediation_attempt"] = 0
    elif isinstance(decision, EscalationIncompatible):
        out["warnings"] = [
            f"INCOMPATIBLE: {decision.reason}. Suggest: {decision.suggested_action}"
        ]

    return out


def _parse_and_validate(raw: str, dna_chain):
    obj = json.loads(_strip_codefence(raw))
    kind = obj.get("kind")

    if kind == "apply_plan":
        codons = split_codons(dna_chain.dna_sequence)
        clean: list[PlannedFix] = []
        for a in obj.get("actions", []):
            ok = True
            for pos, new in zip(a["target_positions"], a["replacement_codons"]):
                if pos >= len(codons):
                    ok = False
                    break
                if translate_codon(codons[pos]) != translate_codon(new.upper()):
                    ok = False
                    break
            if ok:
                clean.append(PlannedFix(
                    check_name=a["check_name"],
                    strategy=a.get("strategy", "escalation_swap"),
                    target_positions=tuple(a["target_positions"]),
                    replacement_codons=tuple(c.upper() for c in a["replacement_codons"]),
                ))
        return EscalationApplyPlan(
            actions=tuple(clean),
            reasoning=obj.get("reasoning", ""),
            diagnosis=obj.get("diagnosis", ""),
        )
    if kind == "incompatible":
        return EscalationIncompatible(
            reason=obj.get("reason", ""),
            suggested_action=obj.get("suggested_action", ""),
        )
    if kind == "change_strategy":
        return EscalationChangeStrategy(
            field=obj["field"],
            new_value=obj["new_value"],
            reason=obj.get("reason", ""),
        )
    return EscalationGiveUp(diagnosis=obj.get("diagnosis", "no diagnosis"))


def _strip_codefence(s: str) -> str:
    s = s.strip()
    if s.startswith("```"):
        s = s.split("\n", 1)[1] if "\n" in s else s[3:]
        if s.endswith("```"):
            s = s[:-3]
    return s.strip()


def _summarize(d) -> str:
    if isinstance(d, EscalationApplyPlan):
        return f"apply_plan ({len(d.actions)} swaps)"
    if isinstance(d, EscalationIncompatible):
        return f"incompatible: {d.reason[:80]}"
    if isinstance(d, EscalationChangeStrategy):
        return f"change_strategy: {d.field}={d.new_value}"
    return f"give_up: {d.diagnosis[:80]}"
