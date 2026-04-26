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


# ---------------------------------------------------------------------------
# Host selector
# ---------------------------------------------------------------------------

def _make_chain(seq):
    from aixbio.models.protein import Chain
    return Chain(id="test", aa_sequence=seq, length=len(seq))


def test_host_selector_ecoli_default():
    from aixbio.tools.host_selector import recommend_host
    # Simple soluble protein — no special features
    chains = [_make_chain("MAAAKKKDDDEEELLLL")]
    rec = recommend_host(chains, "TestProtein")
    assert rec.primary_host == "Escherichia coli"
    assert rec.confidence == "high"
    assert rec.features.n_glycosylation_sites == 0


def test_host_selector_glycosylated_goes_cho():
    from aixbio.tools.host_selector import recommend_host
    # N-X-S/T sequon present (NAS)
    chains = [_make_chain("MAAANASKKK")]
    rec = recommend_host(chains, "GlycoProtein")
    assert rec.primary_host == "CHO"
    assert rec.confidence == "high"
    assert rec.features.n_glycosylation_sites >= 1


def test_host_selector_insulin_ecoli_refolding():
    from aixbio.tools.host_selector import recommend_host
    from aixbio.models.protein import Chain
    # Insulin A + B: 6 Cys total, 51 aa, no glycosylation
    chain_a = Chain(id="A", aa_sequence="GIVEQCCTSICSLYQLENYCN", length=21)
    chain_b = Chain(id="B", aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT", length=30)
    rec = recommend_host([chain_a, chain_b], "Insulin")
    assert rec.primary_host == "Escherichia coli"
    assert rec.features.total_cysteine_count == 6
    assert rec.features.n_glycosylation_sites == 0
    # Must mention refolding in caveats or reasoning
    full_text = rec.reasoning + " ".join(rec.caveats)
    assert "refolding" in full_text.lower()


def test_host_selector_large_protein_cho():
    from aixbio.tools.host_selector import recommend_host
    # > 500 aa, no glyco, no Cys
    chains = [_make_chain("MAAAKKKDDD" * 55)]  # 550 aa
    rec = recommend_host(chains, "BigProtein")
    assert rec.primary_host == "CHO"


def test_host_selector_no_nglycosylation_sequon():
    from aixbio.tools.host_selector import recommend_host
    # NPS is NOT a sequon (P blocks it)
    chains = [_make_chain("MAAANPSKKK")]
    rec = recommend_host(chains, "NPS")
    assert rec.features.n_glycosylation_sites == 0


# ---------------------------------------------------------------------------
# LC-MS/MS prediction
# ---------------------------------------------------------------------------

def test_ms_insulin_b_peptides():
    from aixbio.tools.ms_prediction import predict_ms_peptides
    peptides = predict_ms_peptides("Insulin_B", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    assert len(peptides) >= 2
    names = [p["peptide"] for p in peptides]
    # Trypsin cuts after R at position 22 and K at position 29
    assert any("FVNQHLCGSHLVEALYLVCGER" in n for n in names)
    assert any("GFFYTPK" in n for n in names)


def test_ms_peptide_mass_positive():
    from aixbio.tools.ms_prediction import predict_ms_peptides
    rows = predict_ms_peptides("test", "MAAAKKKDDD")
    for r in rows:
        assert r["mono_mass"] > 0
        assert r["mz_2"] > 0
        assert r["mz_3"] > 0
        assert r["start"] >= 1
        assert r["end"] >= r["start"]


def test_ms_format_tsv():
    from aixbio.tools.ms_prediction import format_tsv, predict_ms_peptides
    rows = predict_ms_peptides("Insulin_B", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    tsv = format_tsv(rows)
    assert tsv.startswith("chain_id\tpeptide\tstart\tend\tlength\tmono_mass\tmz_2+\tmz_3+")
    assert "Insulin_B" in tsv


# ---------------------------------------------------------------------------
# Synthesis feasibility
# ---------------------------------------------------------------------------

def _make_insulin_b_cassette():
    """Build the real insulin B-chain cassette the pipeline would produce."""
    from aixbio.tools.codon_tables import best_ecoli_codon
    gene = "".join(best_ecoli_codon(aa) for aa in "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    # ATG + 6xHis + EK site + gene + TAATAA  (same as cassette_assembly.py)
    return "ATG" + "CACCACCACCACCACCAC" + "GACGACGACAAAGAC" + gene + "TAATAA"


def test_synthesis_insulin_b_idt_feasible():
    from aixbio.tools.synthesis_feasibility import get_synthesis_quotes
    cassette = _make_insulin_b_cassette()
    q = get_synthesis_quotes("Insulin_B", cassette)
    assert q.sequence_length > 125          # within IDT range
    assert 0.25 <= q.gc_content <= 0.65     # healthy GC
    idt = next(vq for vq in q.quotes if vq.vendor == "IDT")
    assert idt.feasible, f"IDT rejected insulin B: {idt.rejection_flags}"
    assert idt.estimated_cost_usd is not None
    assert idt.estimated_cost_usd > 0


def test_synthesis_insulin_b_twist_short():
    from aixbio.tools.synthesis_feasibility import get_synthesis_quotes
    cassette = _make_insulin_b_cassette()
    q = get_synthesis_quotes("Insulin_B", cassette)
    twist = next(vq for vq in q.quotes if vq.vendor == "Twist")
    # Insert + flanks is ~210 bp — below Twist 300 bp minimum
    if q.sequence_length < 300:
        assert not twist.feasible
        assert any("too short" in f.lower() for f in twist.rejection_flags)
    else:
        assert twist.feasible


def test_synthesis_bad_gc_flagged():
    from aixbio.tools.synthesis_feasibility import get_synthesis_quotes
    # Very low GC sequence
    low_gc_cassette = "ATATATAT" * 30   # ~0% GC
    q = get_synthesis_quotes("low_gc_test", low_gc_cassette)
    for vq in q.quotes:
        if not vq.feasible:
            assert any("GC" in f for f in vq.rejection_flags)


def test_synthesis_homopolymer_flagged():
    from aixbio.tools.synthesis_feasibility import get_synthesis_quotes
    # Embed a 12-bp poly-A run into a valid GC sequence
    base = "ATGCGTAAAGGT" * 5  # normal GC
    bad_insert = base + "AAAAAAAAAAAA" + base  # 12-bp poly-A
    q = get_synthesis_quotes("homopolymer_test", bad_insert)
    assert q.longest_homopolymer >= 12
    idt = next(vq for vq in q.quotes if vq.vendor == "IDT")
    assert not idt.feasible
    assert any("homopolymer" in f.lower() or "A/T" in f for f in idt.rejection_flags)


def test_synthesis_format_report():
    from aixbio.tools.synthesis_feasibility import format_quotes_text, get_synthesis_quotes
    cassette = _make_insulin_b_cassette()
    q = get_synthesis_quotes("Insulin_B", cassette)
    report = format_quotes_text([q])
    assert "Insulin_B" in report
    assert "IDT" in report
    assert "Twist" in report
    assert "bp" in report
