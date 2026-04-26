"""Test deterministic nodes with known insulin B-chain sequence."""
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.nodes.codon_optimization import codon_optimization
from aixbio.nodes.cassette_assembly import cassette_assembly
from aixbio.nodes.plasmid_assembly import plasmid_assembly
from aixbio.nodes.sequence_validation import sequence_validation
from aixbio.tools.codon_tables import translate_dna
from aixbio.tools.restriction_sites import get_native_enzymes

INSULIN_B_CHAIN = Chain(
    id="Insulin_B",
    aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    length=30,
)

PROTEIN_RECORD = ProteinRecord(
    uniprot_id="P01308",
    name="Insulin",
    chains=(INSULIN_B_CHAIN,),
)


def _make_chain_state(**overrides):
    base = {
        "chain": INSULIN_B_CHAIN,
        "host_organism": "Escherichia coli",
        "avoid_sites": get_native_enzymes("Escherichia coli"),
        "tag_type": "6xHis",
        "protease_site": "Enterokinase",
        "vector": "pET-28a(+)",
        "cloning_sites": ("BamHI", "XhoI"),
        "protein_record": PROTEIN_RECORD,
        "optimized_dna": None,
        "cassette": None,
        "plasmid": None,
        "chain_validation": None,
        "remediation_attempt": 0,
        "max_remediation_attempts": 3,
        "failed_checks": (),
        "remediation_plan": None,
        "remediation_history": [],
        "decision_log": [],
        "warnings": [],
    }
    base.update(overrides)
    return base


def test_codon_optimization():
    state = _make_chain_state()
    result = codon_optimization(state)
    dna_chain = result["optimized_dna"]

    assert dna_chain.id == "Insulin_B"
    assert len(dna_chain.dna_sequence) == 30 * 3
    translated = translate_dna(dna_chain.dna_sequence)
    assert translated == INSULIN_B_CHAIN.aa_sequence
    assert 0.0 < dna_chain.gc_content < 1.0
    assert 0.0 < dna_chain.cai_score <= 1.0


def test_cassette_assembly():
    state = _make_chain_state()
    opt_result = codon_optimization(state)
    state = _make_chain_state(optimized_dna=opt_result["optimized_dna"])
    result = cassette_assembly(state)
    cassette = result["cassette"]

    assert cassette.id == "Insulin_B"
    assert cassette.full_dna.startswith("ATG")
    assert cassette.full_dna.endswith("TAATAA")
    assert cassette.elements.tag == "CACCACCACCACCACCAC"
    assert opt_result["optimized_dna"].dna_sequence in cassette.full_dna


def test_plasmid_assembly():
    state = _make_chain_state()
    opt_result = codon_optimization(state)
    state = _make_chain_state(optimized_dna=opt_result["optimized_dna"])
    cas_result = cassette_assembly(state)
    state = _make_chain_state(
        optimized_dna=opt_result["optimized_dna"],
        cassette=cas_result["cassette"],
    )
    result = plasmid_assembly(state)
    plasmid = result["plasmid"]

    assert plasmid.id == "Insulin_B"
    assert plasmid.vector == "pET-28a(+)"
    assert plasmid.insert_size > 0
    assert "LOCUS" in plasmid.genbank_file


def test_sequence_validation():
    state = _make_chain_state()
    opt_result = codon_optimization(state)
    state = _make_chain_state(optimized_dna=opt_result["optimized_dna"])
    cas_result = cassette_assembly(state)
    state = _make_chain_state(
        optimized_dna=opt_result["optimized_dna"],
        cassette=cas_result["cassette"],
    )
    result = sequence_validation(state)
    validation = result["chain_validation"]

    assert validation.id == "Insulin_B"
    assert len(validation.checks) == 7

    check_names = {c.name for c in validation.checks}
    assert check_names == {
        "gc_content", "cai_score", "restriction_sites",
        "rna_secondary_structure", "back_translation", "rare_codons",
        "direct_repeats",
    }

    back_translation = next(c for c in validation.checks if c.name == "back_translation")
    assert back_translation.passed, f"Back-translation failed: {back_translation.value}"


def test_full_chain_pipeline():
    """Run Steps 2-5 sequentially for insulin B chain."""
    state = _make_chain_state()

    opt_result = codon_optimization(state)
    state = _make_chain_state(optimized_dna=opt_result["optimized_dna"])

    cas_result = cassette_assembly(state)
    state = _make_chain_state(
        optimized_dna=opt_result["optimized_dna"],
        cassette=cas_result["cassette"],
    )

    plas_result = plasmid_assembly(state)
    state = _make_chain_state(
        optimized_dna=opt_result["optimized_dna"],
        cassette=cas_result["cassette"],
        plasmid=plas_result["plasmid"],
    )

    val_result = sequence_validation(state)
    validation = val_result["chain_validation"]

    back_trans = next(c for c in validation.checks if c.name == "back_translation")
    assert back_trans.passed

    print(f"\nInsulin B-chain pipeline results:")
    for check in validation.checks:
        mark = "PASS" if check.passed else "FAIL"
        print(f"  [{mark}] {check.name}: {check.value} (threshold: {check.threshold})")
    print(f"  Overall: {'PASS' if validation.passed else 'FAIL'}")