from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

from aixbio.config import (
    DEFAULT_CLONING_SITES,
    DEFAULT_HOST,
    DEFAULT_MAX_REMEDIATION_ATTEMPTS,
    DEFAULT_PROTEASE_SITE,
    DEFAULT_TAG_TYPE,
    DEFAULT_VECTOR,
)
from aixbio.graph.main_graph import compile_pipeline
from aixbio.tools.restriction_sites import get_native_enzymes


def main():
    parser = argparse.ArgumentParser(description="Biocompound replica pipeline")
    parser.add_argument("compound_id", help="UniProt accession (e.g. P01308)")
    parser.add_argument("--host", default=DEFAULT_HOST, help="Target host organism")
    parser.add_argument("--structural", action="store_true", help="Run structural validation (Step 6)")
    parser.add_argument("--max-retries", type=int, default=DEFAULT_MAX_REMEDIATION_ATTEMPTS)
    parser.add_argument("--auto-approve", action="store_true", help="Skip human checkpoints")
    parser.add_argument("--escalation", action="store_true", help="Enable LLM escalation agent on remediation exhaustion")
    parser.add_argument("--protocol", action="store_true", help="Generate literature-backed wet-lab SOP after pipeline completes")
    parser.add_argument("--output-dir", default="output", help="Directory for output files")
    args = parser.parse_args()

    app = compile_pipeline()
    config = {"configurable": {"thread_id": f"pipeline-{args.compound_id}"}}

    initial_state = {
        "compound_id": args.compound_id,
        "host_organism": args.host,
        "avoid_sites": get_native_enzymes(args.host),
        "tag_type": DEFAULT_TAG_TYPE,
        "protease_site": DEFAULT_PROTEASE_SITE,
        "vector": DEFAULT_VECTOR,
        "cloning_sites": DEFAULT_CLONING_SITES,
        "run_structural_validation": args.structural,
        "run_protocol_generation": args.protocol,
        "enable_escalation": args.escalation,
        "protein_record": None,
        "host_recommendation": None,
        "protocol": None,
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

    rec = result.get("host_recommendation")
    if rec:
        print(f"\nHost Recommendation: {rec.primary_host} [{rec.confidence} confidence]")
        print(f"  Reasoning: {rec.reasoning}")
        if rec.alternative_hosts:
            print(f"  Alternatives: {', '.join(rec.alternative_hosts)}")
        for caveat in rec.caveats:
            print(f"  ! {caveat}")

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

    structure = result.get("structure_report")
    if structure:
        print(f"\nStructural Validation (AlphaFold):")
        for sr in structure.chains:
            rmsd_str = f"{sr.rmsd_to_ref:.3f} Å" if sr.rmsd_to_ref is not None else "N/A"
            print(f"  Chain {sr.id}: pLDDT={sr.plddt_mean:.1f}, RMSD={rmsd_str}")
            if sr.structure_file:
                print(f"    Structure file: {sr.structure_file}")

    chain_results = result.get("chain_results", [])
    if chain_results:
        print(f"\nChain Results ({len(chain_results)}):")
        for cr in chain_results:
            print(f"  {cr['chain_id']}:")
            print(f"    Insert size: {cr['insert_size']} bp")
            print(f"    Remediation rounds: {cr['remediation_rounds']}")
            print(f"    Validation passed: {cr['validation_passed']}")
            sol_score = cr.get("solubility_score")
            if sol_score is not None:
                sol_label = "SOLUBLE" if sol_score >= 0.45 else "INCLUSION BODY RISK"
                print(f"    Solubility: {sol_score:.2f} [{sol_label}]")
            if cr.get("disulfide_risk"):
                print(f"    Disulfide risk: YES (refolding required)")

    protocol = result.get("protocol")
    if protocol:
        print(f"\nProtocol generated ({len(protocol)} chars) — written to output/protocol.md")

    # Synthesis feasibility (always computed if chain results are available)
    chain_results = result.get("chain_results", [])
    if chain_results:
        from aixbio.tools.synthesis_feasibility import get_synthesis_quotes
        print(f"\nSynthesis Feasibility:")
        for cr in chain_results:
            if not cr.get("cassette_dna"):
                continue
            q = get_synthesis_quotes(cr["chain_id"], cr["cassette_dna"])
            for vq in q.quotes:
                status = "FEASIBLE" if vq.feasible else "FLAGGED"
                cost = f"~${vq.estimated_cost_usd:.0f}" if vq.estimated_cost_usd else "N/A"
                print(f"  {cr['chain_id']} / {vq.vendor}: [{status}] {cost}")
                for flag in vq.rejection_flags:
                    print(f"    ! {flag}")

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
                "solubility_score": cr.get("solubility_score"),
                "disulfide_risk": cr.get("disulfide_risk", False),
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

    rec = result.get("host_recommendation")
    if rec:
        summary["host_recommendation"] = {
            "primary_host": rec.primary_host,
            "confidence": rec.confidence,
            "reasoning": rec.reasoning,
            "alternative_hosts": list(rec.alternative_hosts),
            "caveats": list(rec.caveats),
            "features": {
                "total_chain_length": rec.features.total_chain_length,
                "n_glycosylation_sites": rec.features.n_glycosylation_sites,
                "total_cysteine_count": rec.features.total_cysteine_count,
                "disulfide_risk": rec.features.disulfide_risk,
                "max_gravy": rec.features.max_gravy,
            },
        }

    json_path = out / f"{compound_id}_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"  Wrote {json_path}")

    protocol = result.get("protocol")
    if protocol:
        protocol_path = out / "protocol.md"
        protocol_path.write_text(protocol)
        print(f"  Wrote {protocol_path}")

    # Synthesis feasibility report
    synthesis_quotes = []
    for cr in chain_results:
        if cr.get("cassette_dna"):
            from aixbio.tools.synthesis_feasibility import format_quotes_text, get_synthesis_quotes
            q = get_synthesis_quotes(cr["chain_id"], cr["cassette_dna"])
            synthesis_quotes.append(q)

    if synthesis_quotes:
        from aixbio.tools.synthesis_feasibility import format_quotes_text
        report_path = out / "synthesis_report.txt"
        report_path.write_text(format_quotes_text(synthesis_quotes))
        print(f"  Wrote {report_path}")

        # Embed compact summary in JSON
        summary["synthesis_quotes"] = [
            {
                "chain_id": q.chain_id,
                "insert_length": q.sequence_length,
                "gc_content": q.gc_content,
                "longest_homopolymer": q.longest_homopolymer,
                "vendors": [
                    {
                        "vendor": vq.vendor,
                        "feasible": vq.feasible,
                        "estimated_cost_usd": vq.estimated_cost_usd,
                        "estimated_turnaround": vq.estimated_turnaround,
                        "rejection_flags": list(vq.rejection_flags),
                        "notes": list(vq.notes),
                    }
                    for vq in q.quotes
                ],
            }
            for q in synthesis_quotes
        ]

    # LC-MS/MS peptide mass tables (one TSV per chain, always written)
    protein_record = result.get("protein_record")
    if protein_record:
        from aixbio.tools.ms_prediction import format_tsv, predict_ms_peptides
        for chain in protein_record.chains:
            safe_id = chain.id.replace(" ", "_").replace("/", "_")
            rows = predict_ms_peptides(chain.id, chain.aa_sequence)
            if rows:
                tsv_path = out / f"{safe_id}_peptides.tsv"
                tsv_path.write_text(format_tsv(rows))
                print(f"  Wrote {tsv_path} ({len(rows)} tryptic peptides)")

    print(f"\nAll artifacts written to {out}/")


if __name__ == "__main__":
    sys.exit(main())
