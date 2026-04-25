from __future__ import annotations

from aixbio.models.validation import ChainValidation, ValidationReport
from aixbio.state.chain_state import ChainSubgraphState
from aixbio.state.pipeline_state import ChainProcessingResult, PipelineState


def package_result(state: ChainSubgraphState) -> dict:
    validation = state["chain_validation"]
    is_passed = validation.passed if validation else False
    return {
        "chain_results": [ChainProcessingResult(
            chain_id=state["chain"].id,
            optimized_dna_sequence=state["optimized_dna"].dna_sequence if state["optimized_dna"] else "",
            cassette_dna=state["cassette"].full_dna if state["cassette"] else "",
            genbank_file=state["plasmid"].genbank_file if state["plasmid"] else "",
            insert_size=state["plasmid"].insert_size if state["plasmid"] else 0,
            validation_passed=is_passed,
            checks=validation.checks if validation else (),
            remediation_rounds=state["remediation_attempt"],
            remediation_history=tuple(state.get("remediation_history", [])),
            status="passed" if is_passed else "failed",
        )],
    }


def package_result_failed(state: ChainSubgraphState) -> dict:
    """Package result for chains that exhausted all remediation attempts."""
    result = package_result(state)
    # Override status to distinguish from a simple validation failure
    result["chain_results"][0]["status"] = "max_retries_exceeded"
    return result


def halt_pipeline(state: ChainSubgraphState) -> dict:
    result = package_result(state)
    result["chain_results"][0]["status"] = "failed"
    result["warnings"] = [
        f"HALT: Back-translation check failed for chain {state['chain'].id}. Pipeline bug detected."
    ]
    return result


def package_result_escalated(state: ChainSubgraphState) -> dict:
    result = package_result(state)
    decision = state.get("escalation_decision")
    if decision is not None and decision.kind == "incompatible":
        result["chain_results"][0]["status"] = "host_incompatible"
    else:
        result["chain_results"][0]["status"] = "escalation_failed"
    return result


def merge_all_chain_results(state: PipelineState) -> dict:
    chain_results = state.get("chain_results", [])
    chain_validations = []
    for cr in chain_results:
        chain_validations.append(ChainValidation(
            id=cr["chain_id"],
            passed=cr["validation_passed"],
            checks=cr["checks"],
        ))

    all_passed = all(cv.passed for cv in chain_validations)

    return {
        "validation_report": ValidationReport(
            chains=tuple(chain_validations),
            all_passed=all_passed,
        ),
        "pipeline_status": "completed" if all_passed else "failed",
    }
