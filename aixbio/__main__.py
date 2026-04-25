from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

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
    parser.add_argument("--escalation", action="store_true", help="Enable LLM escalation agent on remediation exhaustion")
    parser.add_argument("--output-dir", default="output", help="Directory for output files")
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
        "enable_escalation": args.escalation,
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
    _write_artifacts(result, args.compound_id, args.output_dir)
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


def _write_artifacts(result: dict, compound_id: str, output_dir: str):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    chain_results = result.get("chain_results", [])
    if not chain_results:
        return

    for cr in chain_results:
        chain_id = cr["chain_id"]
        safe_id = chain_id.replace(" ", "_").replace("/", "_")

        if cr.get("optimized_dna_sequence"):
            fasta_path = out / f"{safe_id}.fasta"
            fasta_path.write_text(
                f">{chain_id} optimized for E. coli\n{cr['optimized_dna_sequence']}\n"
            )
            print(f"  Wrote {fasta_path}")

        if cr.get("genbank_file"):
            gb_path = out / f"{safe_id}_plasmid.gb"
            gb_path.write_text(cr["genbank_file"])
            print(f"  Wrote {gb_path}")

    summary = {
        "compound_id": compound_id,
        "status": result.get("pipeline_status", "unknown"),
        "warnings": result.get("warnings", []),
        "chains": [
            {
                "chain_id": cr["chain_id"],
                "insert_size": cr["insert_size"],
                "validation_passed": cr["validation_passed"],
                "remediation_rounds": cr["remediation_rounds"],
                "status": cr["status"],
                "checks": [
                    {"name": c.name, "passed": c.passed, "value": c.value, "threshold": c.threshold}
                    for c in cr.get("checks", ())
                ],
            }
            for cr in chain_results
        ],
    }
    if result.get("protein_record"):
        pr = result["protein_record"]
        summary["protein"] = {
            "name": pr.name,
            "uniprot_id": pr.uniprot_id,
            "chains": [{"id": c.id, "length": c.length} for c in pr.chains],
        }

    json_path = out / f"{compound_id}_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"  Wrote {json_path}")
    print(f"\nAll artifacts written to {out}/")


if __name__ == "__main__":
    sys.exit(main())
