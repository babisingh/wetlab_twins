from aixbio.tools.codon_tables import (
    best_ecoli_codon,
    split_codons,
    synonymous_alternatives,
    translate_codon,
    translate_dna,
    RARE_CODONS_ECOLI,
)
from aixbio.tools.gc import compute_gc
from aixbio.tools.cai import compute_cai
from aixbio.tools.restriction_sites import (
    find_restriction_sites,
    get_native_enzymes,
    get_recognition_site,
    has_restriction_sites,
)


def test_translate_codon():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"
    assert translate_codon("CGT") == "R"


def test_best_ecoli_codon():
    codon = best_ecoli_codon("R")
    assert translate_codon(codon) == "R"
    assert codon not in RARE_CODONS_ECOLI


def test_split_codons():
    assert split_codons("ATGCGT") == ["ATG", "CGT"]
    assert split_codons("ATGCG") == ["ATG", "CG"]


def test_translate_dna():
    assert translate_dna("ATGCGT") == "MR"


def test_synonymous_alternatives():
    alts = synonymous_alternatives("CGT")
    assert "CGT" not in alts
    assert all(translate_codon(a) == "R" for a in alts)


def test_gc_content():
    assert compute_gc("GCGC") == 1.0
    assert compute_gc("ATAT") == 0.0
    assert abs(compute_gc("GCATAT") - 1 / 3) < 0.01


def test_cai_perfect():
    aa = "MRK"
    dna = "".join(best_ecoli_codon(a) for a in aa)
    cai = compute_cai(dna)
    assert cai > 0.99


def test_restriction_sites():
    dna = "AAAGGATCCAAA"  # contains BamHI (GGATCC)
    sites = find_restriction_sites(dna, ("BamHI", "XhoI"))
    assert len(sites) == 1
    assert sites[0][0] == "BamHI"
    assert has_restriction_sites(dna, ("BamHI",))
    assert not has_restriction_sites("AAAAAA", ("BamHI",))


def test_rna_fold():
    from aixbio.tools.rna_fold import estimate_five_prime_dg
    dg = estimate_five_prime_dg("ATGATGATGATGATG")
    assert isinstance(dg, float)


def test_native_enzymes_ecoli():
    enzymes = get_native_enzymes("Escherichia coli")
    assert len(enzymes) > 0
    assert all(e.startswith("Eco") for e in enzymes)


def test_native_enzymes_unknown_organism():
    enzymes = get_native_enzymes("Unknown organism")
    assert enzymes == ()


def test_recognition_site():
    assert get_recognition_site("EcoRI") == "GAATTC"
    assert get_recognition_site("BamHI") == "GGATCC"
    assert get_recognition_site("NotI") == "GCGGCCGC"


# ---------------------------------------------------------------------------
# Solubility tool
# ---------------------------------------------------------------------------

def test_solubility_insulin_b():
    from aixbio.tools.solubility import predict_solubility
    # Insulin B: 30 aa, 2 Cys — disulfide risk must be flagged
    result = predict_solubility("Insulin_B", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    assert 0.0 <= result.score <= 1.0
    assert result.cysteine_count == 2
    assert result.disulfide_risk is True
    assert result.method == "composition_heuristic"
    assert result.id == "Insulin_B"
    assert isinstance(result.reasoning, str) and result.reasoning


def test_solubility_insulin_a():
    from aixbio.tools.solubility import predict_solubility
    # Insulin A: 21 aa, 4 Cys — higher disulfide risk
    result = predict_solubility("Insulin_A", "GIVEQCCTSICSLYQLENYCN")
    assert result.cysteine_count == 4
    assert result.disulfide_risk is True
    assert result.tag_recommendation is not None  # some recommendation expected


def test_solubility_hydrophobic_protein():
    from aixbio.tools.solubility import predict_solubility
    # Highly hydrophobic sequence — should score low
    result = predict_solubility("hydrophobic", "ILILILILIVVVVVLLLLFFFF")
    assert result.score < 0.45, f"Expected inclusion body risk, got score={result.score}"
    assert not result.predicted_soluble


def test_solubility_charged_protein():
    from aixbio.tools.solubility import predict_solubility
    # Highly charged sequence — should score higher
    result = predict_solubility("charged", "DEKRDEKRDEKRDEKRDEKR")
    assert result.score > 0.45, f"Expected soluble, got score={result.score}"
    assert result.predicted_soluble


def test_solubility_node_returns_warnings_for_disulfide():
    from aixbio.models.protein import Chain, ProteinRecord
    from aixbio.nodes.solubility_prediction import solubility_prediction
    from aixbio.tools.restriction_sites import get_native_enzymes

    chain = Chain(id="Insulin_A", aa_sequence="GIVEQCCTSICSLYQLENYCN", length=21)
    protein = ProteinRecord(uniprot_id="P01308", name="Insulin", chains=(chain,))
    state = {
        "chain": chain,
        "protein_record": protein,
        "host_organism": "Escherichia coli",
        "avoid_sites": get_native_enzymes("Escherichia coli"),
        "tag_type": "6xHis",
        "protease_site": "Enterokinase",
        "vector": "pET-28a(+)",
        "cloning_sites": ("BamHI", "XhoI"),
    }
    result = solubility_prediction(state)

    assert "solubility_result" in result
    sol = result["solubility_result"]
    assert sol.disulfide_risk is True
    assert sol.cysteine_count == 4
    # Node must surface warnings for disulfide risk
    warnings = result.get("warnings", [])
    assert any("disulfide" in w.lower() for w in warnings)
