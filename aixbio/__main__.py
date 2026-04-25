from __future__ import annotations

import argparse
import json
import sys

from dotenv import load_dotenv

load_dotenv()

from aixbio.config import (
    DEFAULT_AVOID_SITES,
    DEFAULT_CLONING_SITES,
    DEFAULT_HOST,
    DEFAULT_MAX_REMEDIATION_ATTEMPTS,
    DEFAULT_PROTEASE_SITE,
    DEFAULT_TAG_TYPE,
    DEFAULT_VECTOR,
)
from aixbio.graph.main_graph import compile_pipeline


def main():
    parser = argparse.ArgumentParser(description="Biocompound replica pipeline")
    parser.add_argument("compound_id", help="UniProt accession (e.g. P01308)")
    parser.add_argument("--host", default=DEFAULT_HOST, help="Target host organism")
    parser.add_argument("--structural", action="store_true", help="Run structural validation (Step 6)")
    parser.add_argument("--max-retries", type=int, default=DEFAULT_MAX_REMEDIATION_ATTEMPTS)
    parser.add_argument("--auto-approve", action="store_true", help="Skip human checkpoints")
    args = parser.parse_args()

    app = compile_pipeline()
    config = {"configurable": {"thread_id": f"pipeline-{args.compound_id}"}}

    initial_state = {
        "compound_id": args.compound_id,
        "host_organism": args.host,
        "avoid_sites": DEFAULT_AVOID_SITES,
        "tag_type": DEFAULT_TAG_TYPE,
        "protease_site": DEFAULT_PROTEASE_SITE,
        "vector": DEFAULT_VECTOR,
        "cloning_sites": DEFAULT_CLONING_SITES,
        "run_structural_validation": args.structural,
        "protein_record": None,
        "chain_extraction_reasoning": "",
        "chain_results": [],
        "validation_report": None,
        "structure_report": None,
        "pipeline_status": "running",
        "max_remediation_attempts": args.max_retries,
        "decision_log": [],
        "warnings": [],
        "human_approval_chain_extraction": None,
        "human_approval_plasmid": None,
    }

    print(f"Starting pipeline for {args.compound_id} targeting {args.host}")
    print("=" * 60)

    result = app.invoke(initial_state, config)

    while _has_interrupt(result):
        interrupt_data = _get_interrupt_data(result)
        if interrupt_data:
            print(f"\n{'=' * 60}")
            print(f"HUMAN REVIEW: {interrupt_data.get('stage', 'unknown')}")
            print(f"{'=' * 60}")
            print(json.dumps(interrupt_data, indent=2, default=str))
            print()

            if args.auto_approve:
                human_input = "approve"
                print("Auto-approving (--auto-approve flag)")
            else:
                human_input = input("Your decision (approve/reject): ").strip()

            from langgraph.types import Command
            result = app.invoke(Command(resume=human_input), config)
        else:
            break

    _print_results(result)
    return 0


def _has_interrupt(result) -> bool:
    if isinstance(result, dict):
        return "__interrupt__" in result
    return False


def _get_interrupt_data(result) -> dict | None:
    interrupts = result.get("__interrupt__", [])
    if interrupts:
        return interrupts[0].value if hasattr(interrupts[0], "value") else interrupts[0]
    return None


def _print_results(result: dict):
    print(f"\n{'=' * 60}")
    print("PIPELINE RESULTS")
    print(f"{'=' * 60}")

    status = result.get("pipeline_status", "unknown")
    print(f"Status: {status}")

    warnings = result.get("warnings", [])
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for w in warnings:
            print(f"  - {w}")

    validation = result.get("validation_report")
    if validation:
        print(f"\nValidation: {'PASSED' if validation.all_passed else 'FAILED'}")
        for cv in validation.chains:
            print(f"  Chain {cv.id}: {'PASS' if cv.passed else 'FAIL'}")
            for check in cv.checks:
                mark = "+" if check.passed else "X"
                print(f"    [{mark}] {check.name}: {check.value} (threshold: {check.threshold})")

    chain_results = result.get("chain_results", [])
    if chain_results:
        print(f"\nChain Results ({len(chain_results)}):")
        for cr in chain_results:
            print(f"  {cr['chain_id']}:")
            print(f"    Insert size: {cr['insert_size']} bp")
            print(f"    Remediation rounds: {cr['remediation_rounds']}")
            print(f"    Validation passed: {cr['validation_passed']}")

    decisions = result.get("decision_log", [])
    if decisions:
        print(f"\nAgent Decisions ({len(decisions)}):")
        for d in decisions:
            print(f"  [{d.node}] {d.action}")
            print(f"    Reasoning: {d.reasoning}")


if __name__ == "__main__":
    sys.exit(main())
