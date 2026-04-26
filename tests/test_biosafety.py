"""Tests for the biosafety screening module."""
from aixbio.tools.biosafety import (
    BLOCKED_UNIPROT_IDS,
    DANGEROUS_KEYWORDS,
    screen_compound_id,
    screen_protein,
)


# ---------------------------------------------------------------------------
# Layer 1 — UniProt ID blocklist
# ---------------------------------------------------------------------------

def test_ricin_blocked_by_id():
    result = screen_compound_id("P02879")
    assert not result.safe
    assert result.match_type == "uniprot_id"
    assert "Ricin" in result.matched_agent
    assert "P02879" in result.reason


def test_botulinum_type_a_blocked():
    result = screen_compound_id("P10844")
    assert not result.safe
    assert "Botulinum" in result.matched_agent


def test_tetanus_toxin_blocked():
    result = screen_compound_id("P04958")
    assert not result.safe


def test_anthrax_lf_blocked():
    result = screen_compound_id("P15917")
    assert not result.safe
    assert "Anthrax" in result.matched_agent


def test_insulin_passes_id_check():
    result = screen_compound_id("P01308")
    assert result.safe
    assert result.matched_agent is None
    assert result.match_type is None


def test_lowercase_id_normalised():
    """Accession lookup must be case-insensitive."""
    result = screen_compound_id("p02879")
    assert not result.safe


def test_unknown_id_passes():
    result = screen_compound_id("XXXXX")
    assert result.safe


def test_all_blocked_ids_have_agent_names():
    """Every blocked accession must have a human-readable name."""
    from aixbio.tools.biosafety import _UNIPROT_AGENT_NAMES
    for uid in BLOCKED_UNIPROT_IDS:
        assert uid in _UNIPROT_AGENT_NAMES, f"Missing agent name for {uid}"


# ---------------------------------------------------------------------------
# Layer 2 — Protein name keyword screen
# ---------------------------------------------------------------------------

def test_ricin_blocked_by_keyword():
    result = screen_protein("Ricin precursor", [])
    assert not result.safe
    assert result.match_type == "keyword"
    assert "ricin" in result.reason.lower()


def test_botulinum_keyword():
    result = screen_protein("Botulinum neurotoxin type A heavy chain", [])
    assert not result.safe
    assert result.match_type == "keyword"


def test_shiga_toxin_keyword():
    result = screen_protein("Shiga toxin 1 A subunit", [])
    assert not result.safe


def test_safe_protein_name():
    result = screen_protein("Insulin", ["MAAAKKKDDDEEELLLL"])
    assert result.safe


def test_keyword_check_case_insensitive():
    result = screen_protein("RICIN A-CHAIN FRAGMENT", [])
    assert not result.safe


# ---------------------------------------------------------------------------
# Layer 3 — Amino-acid signature scan
# ---------------------------------------------------------------------------

def test_ricin_a_chain_signature_detected():
    # Insert the published N-terminal diagnostic signature into a dummy chain.
    seq = "AAAAAAA" + "IFPKQYPIINFT" + "AAAAAAA"
    result = screen_protein("Unknown protein", [seq])
    assert not result.safe
    assert result.match_type == "aa_signature"
    assert "Ricin" in result.matched_agent


def test_botulinum_zinc_motif_detected():
    seq = "MAAAA" + "HELIH" + "KKKKKK"
    result = screen_protein("Uncharacterised protein", [seq])
    assert not result.safe
    assert result.match_type == "aa_signature"
    assert "Botulinum" in result.matched_agent


def test_signature_check_is_case_insensitive():
    seq = "aaaa" + "ifpkqypiinft" + "aaaa"
    result = screen_protein("Unknown", [seq])
    assert not result.safe


def test_safe_sequence_passes_all_layers():
    insulin_b = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    result = screen_protein("Insulin", [insulin_b])
    assert result.safe


# ---------------------------------------------------------------------------
# Node-level integration
# ---------------------------------------------------------------------------

def test_biosafety_node_blocks_ricin():
    from aixbio.models.protein import Chain, ProteinRecord
    from aixbio.nodes.biosafety_screen import biosafety_screen

    chain = Chain(id="A", aa_sequence="IFPKQYPIINFTAAAAAAAA", length=20)
    protein = ProteinRecord(uniprot_id="P02879", name="Ricin precursor", chains=(chain,))
    state = {
        "protein_record": protein,
        "pipeline_status": "running",
        "warnings": [],
    }
    result = biosafety_screen(state)
    assert result["pipeline_status"] == "biosafety_rejected"
    assert result["biosafety_result"].safe is False
    assert any("BIOSAFETY" in w for w in result["warnings"])


def test_biosafety_node_passes_insulin():
    from aixbio.models.protein import Chain, ProteinRecord
    from aixbio.nodes.biosafety_screen import biosafety_screen

    chain = Chain(id="B", aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT", length=30)
    protein = ProteinRecord(uniprot_id="P01308", name="Insulin", chains=(chain,))
    state = {
        "protein_record": protein,
        "pipeline_status": "running",
        "warnings": [],
    }
    result = biosafety_screen(state)
    assert result.get("pipeline_status") != "biosafety_rejected"
    assert result["biosafety_result"].safe is True
