from __future__ import annotations

import asyncio
import json
from datetime import datetime, timezone

from langchain_core.messages import HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI

from aixbio.config import LLM_MODEL, OPENROUTER_API_KEY, OPENROUTER_BASE_URL
from aixbio.models.audit import AgentDecision
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.prompts.sequence_retrieval import SEQUENCE_RETRIEVAL_SYSTEM
from aixbio.state.pipeline_state import PipelineState
from aixbio.tools.uniprot import extract_protein_name, extract_sequence, fetch_uniprot_entry


def sequence_retrieval_agent(state: PipelineState) -> dict:
    compound_id = state["compound_id"]

    entry = asyncio.run(fetch_uniprot_entry(compound_id))
    full_sequence = extract_sequence(entry)
    protein_name = extract_protein_name(entry)

    features_summary = json.dumps(entry.get("features", []), indent=2, default=str)

    llm = ChatOpenAI(
        model=LLM_MODEL,
        temperature=0,
        openai_api_key=OPENROUTER_API_KEY,
        openai_api_base=OPENROUTER_BASE_URL,
    )
    response = llm.invoke([
        SystemMessage(content=SEQUENCE_RETRIEVAL_SYSTEM),
        HumanMessage(content=(
            f"UniProt ID: {compound_id}\n"
            f"Protein name: {protein_name}\n"
            f"Full sequence ({len(full_sequence)} aa): {full_sequence}\n\n"
            f"Feature annotations:\n{features_summary}\n\n"
            "Extract the mature chain(s) for E. coli expression."
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
        parsed = _fallback_single_chain(compound_id, protein_name, full_sequence)

    chains = tuple(
        Chain(
            id=c["id"],
            aa_sequence=c["aa_sequence"],
            length=c["length"],
        )
        for c in parsed["chains"]
    )

    protein_record = ProteinRecord(
        uniprot_id=compound_id,
        name=protein_name,
        chains=chains,
    )

    reasoning = parsed.get("reasoning", "No reasoning provided")

    decision = AgentDecision(
        node="sequence_retrieval_agent",
        reasoning=reasoning,
        action=f"Extracted {len(chains)} chain(s) from {compound_id}",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_summary=f"UniProt {compound_id} ({protein_name}, {len(full_sequence)} aa)",
        output_summary=", ".join(f"{c.id}: {c.length} aa" for c in chains),
    )

    return {
        "protein_record": protein_record,
        "chain_extraction_reasoning": reasoning,
        "decision_log": [decision],
    }


def _fallback_single_chain(
    compound_id: str, name: str, sequence: str
) -> dict:
    return {
        "chains": [{
            "id": f"{name}_{compound_id}",
            "aa_sequence": sequence,
            "length": len(sequence),
        }],
        "reasoning": "Fallback: LLM response parsing failed, using full sequence.",
    }
