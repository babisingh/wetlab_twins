"""Integration test: run the chain subgraph end-to-end (no LLM calls needed)."""
from aixbio.graph.chain_subgraph import compile_chain_subgraph
from aixbio.models.protein import Chain, ProteinRecord


INSULIN_B = Chain(
    id="Insulin_B",
    aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    length=30,
)

PROTEIN = ProteinRecord(uniprot_id="P01308", name="Insulin", chains=(INSULIN_B,))


def _make_input_state():
    return {
        "chain": INSULIN_B,
        "host_organism": "E. coli K12",
        "avoid_sites": ("BamHI", "XhoI", "EcoRI"),
        "tag_type": "6xHis",
        "protease_site": "Enterokinase",
        "vector": "pET-28a(+)",
        "cloning_sites": ("BamHI", "XhoI"),
        "protein_record": PROTEIN,
        "optimized_dna": None,
        "cassette": None,
        "plasmid": None,
        "chain_validation": None,
        "remediation_attempt": 0,
        "max_remediation_attempts": 0,
        "failed_checks": (),
        "remediation_plan": None,
        "chain_results": [],
        "remediation_history": [],
        "enable_escalation": False,
        "escalation_used": False,
        "escalation_decision": None,
        "decision_log": [],
        "warnings": [],
    }


def test_chain_subgraph_completes():
    """With max_remediation_attempts=0, the subgraph runs Steps 2-5 and terminates."""
    graph = compile_chain_subgraph()
    result = graph.invoke(_make_input_state())

    assert result["chain_validation"] is not None
    validation = result["chain_validation"]
    assert validation.id == "Insulin_B"

    assert result["plasmid"] is not None
    assert result["plasmid"].insert_size > 0
    assert "LOCUS" in result["plasmid"].genbank_file

    back_trans_passed = any(
        c.name == "back_translation" and c.passed
        for c in validation.checks
    )
    assert back_trans_passed, "Back-translation check must pass"

    assert len(result["chain_results"]) == 1
    cr = result["chain_results"][0]
    assert cr["chain_id"] == "Insulin_B"
    assert cr["validation_passed"] == validation.passed
