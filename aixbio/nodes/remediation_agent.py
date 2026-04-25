from __future__ import annotations

import json
from datetime import datetime, timezone

from langchain_core.messages import HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI

from aixbio.config import LLM_MODEL, OPENROUTER_API_KEY, OPENROUTER_BASE_URL
from aixbio.models.audit import AgentDecision
from aixbio.models.remediation import PlannedFix, RemediationAction, RemediationPlan
from aixbio.prompts.remediation import REMEDIATION_SYSTEM
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.tools.codon_tables import split_codons, translate_codon


def remediation_agent(state: ChainSubgraphState) -> dict:
    dna_chain = state["optimized_dna"]
    failed_checks = state["failed_checks"]
    history = state.get("remediation_history", [])
    attempt = state["remediation_attempt"]

    codons = split_codons(dna_chain.dna_sequence)
    codons_display = " ".join(f"{i}:{c}" for i, c in enumerate(codons))

    failed_summary = json.dumps([
        {"name": c.name, "value": str(c.value), "threshold": c.threshold}
        for c in failed_checks
    ], indent=2)

    history_summary = ""
    if history:
        history_summary = "\nPrevious fix attempts:\n" + json.dumps([
            {"check": h.check_name, "positions": h.positions_affected, "from": h.codons_before, "to": h.codons_after}
            for h in history
        ], indent=2, default=str)

    llm = ChatOpenAI(
        model=LLM_MODEL,
        temperature=0,
        openai_api_key=OPENROUTER_API_KEY,
        openai_api_base=OPENROUTER_BASE_URL,
    )
    response = llm.invoke([
        SystemMessage(content=REMEDIATION_SYSTEM),
        HumanMessage(content=(
            f"Chain: {dna_chain.id}\n"
            f"Attempt: {attempt + 1}\n"
            f"DNA codons ({len(codons)} total):\n{codons_display}\n\n"
            f"Failed checks:\n{failed_summary}\n"
            f"{history_summary}\n\n"
            "Produce targeted fixes."
        )),
    ])

    response_text = response.content
    if isinstance(response_text, list):
        response_text = response_text[0].get("text", "") if response_text else ""

    try:
        json_start = response_text.index("{")
        json_end = response_text.rindex("}") + 1
        parsed = json.loads(response_text[json_start:json_end])
    except (ValueError, json.JSONDecodeError):
        parsed = {"actions": [], "reasoning": "Failed to parse LLM response", "priority_order": []}

    plan = RemediationPlan(
        actions=tuple(
            PlannedFix(
                check_name=a["check_name"],
                strategy=a.get("strategy", "unknown"),
                target_positions=tuple(a["target_positions"]),
                replacement_codons=tuple(a["replacement_codons"]),
            )
            for a in parsed.get("actions", [])
        ),
        reasoning=parsed.get("reasoning", ""),
        priority_order=tuple(parsed.get("priority_order", [])),
    )

    decision = AgentDecision(
        node="remediation_agent",
        reasoning=plan.reasoning,
        action=f"Planned {len(plan.actions)} fixes for {len(failed_checks)} failed checks",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"Chain {dna_chain.id}, attempt {attempt + 1}, {len(failed_checks)} failures",
        output_summary=f"Fixes: {', '.join(a.check_name for a in plan.actions)}",
    )

    return {
        "remediation_plan": plan,
        "decision_log": [decision],
    }


def apply_fixes(state: ChainSubgraphState) -> dict:
    from aixbio.models.dna import DNAChain
    from aixbio.tools.cai import compute_cai
    from aixbio.tools.gc import compute_gc

    dna_chain = state["optimized_dna"]
    plan = state["remediation_plan"]
    codons = split_codons(dna_chain.dna_sequence)

    actions: list[RemediationAction] = []

    for fix in plan.actions:
        for pos, new_codon in zip(fix.target_positions, fix.replacement_codons):
            if pos >= len(codons):
                continue
            old_codon = codons[pos]
            old_aa = translate_codon(old_codon)
            new_aa = translate_codon(new_codon.upper())
            if old_aa != new_aa:
                continue

            actions.append(RemediationAction(
                check_name=fix.check_name,
                fix_type=fix.strategy,
                positions_affected=(pos,),
                codons_before=(old_codon,),
                codons_after=(new_codon.upper(),),
                reasoning=plan.reasoning,
            ))
            codons[pos] = new_codon.upper()

    new_dna = "".join(codons)
    new_dna_chain = DNAChain(
        id=dna_chain.id,
        dna_sequence=new_dna,
        cai_score=compute_cai(new_dna),
        gc_content=compute_gc(new_dna),
    )

    return {
        "optimized_dna": new_dna_chain,
        "remediation_history": actions,
        "remediation_attempt": state["remediation_attempt"] + 1,
    }
