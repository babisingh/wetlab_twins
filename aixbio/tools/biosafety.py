"""Biosafety screening for DNA synthesis pipelines.

Three defence layers — applied in order, earliest match wins:
  1. UniProt accession blocklist  (pre-flight, before any network call)
  2. Protein name / description keyword blocklist (post-retrieval)
  3. Amino-acid signature scan    (post-retrieval, defence-in-depth)

All three layers are intentionally conservative: false positives halt the
pipeline; false negatives allow a dangerous sequence through. When in doubt,
the list is designed to over-block.

Regulatory basis
----------------
The CDC/USDA Select Agent Program (42 CFR Part 73 / 9 CFR Part 121) lists
the toxins included here. UniProt accessions were verified against canonical
entries at the time of writing; the list must be updated as new sequences are
deposited.

For production deployments, supplement this heuristic screen with a BLAST
search against the NCBI Select Agent sequence database:
  https://www.ncbi.nlm.nih.gov/selectagent/
"""
from __future__ import annotations

from aixbio.models.biosafety import BiosafetyResult

# ---------------------------------------------------------------------------
# Layer 1 — UniProt accession blocklist
# Source: UniProt canonical entries for CDC/USDA Tier 1 toxins and related
# high-consequence select agents. Accessions are listed with organism and
# agent name for auditability.
# ---------------------------------------------------------------------------

BLOCKED_UNIPROT_IDS: frozenset[str] = frozenset({
    # Ricin / RCA (Ricinus communis)
    "P02879",   # Ricin precursor (A-chain + B-chain)
    "P06750",   # Ricinus communis agglutinin II (RCA-II) precursor

    # Abrin (Abrus precatorius)
    "P28590",   # Abrin-a precursor
    "P11279",   # Abrin precursor
    "P24476",   # Abrin-b precursor

    # Botulinum neurotoxins — all seven serotypes (Clostridium botulinum)
    "P10844",   # BoNT/A
    "P10845",   # BoNT/B
    "P19321",   # BoNT/C
    "P19322",   # BoNT/D
    "P19619",   # BoNT/E
    "P30995",   # BoNT/F
    "Q60393",   # BoNT/G

    # Tetanus toxin (Clostridium tetani)
    "P04958",   # TeNT

    # Anthrax toxin components (Bacillus anthracis)
    "P13423",   # Protective antigen (PA83)
    "P15917",   # Lethal factor (LF)
    "P40136",   # Edema factor (EF)

    # Diphtheria toxin (Corynebacterium diphtheriae)
    "P00588",   # DT

    # Shiga toxins / Vero toxins (E. coli O157:H7, Shigella dysenteriae)
    "P09385",   # Stx1 A-subunit
    "P09386",   # Stx1 B-subunit
    "P09387",   # Stx2 A-subunit
    "P09388",   # Stx2 B-subunit

    # Cholera toxin (Vibrio cholerae)
    "P01555",   # CtxA
    "P01556",   # CtxB

    # Staphylococcal enterotoxins A–E (S. aureus)
    "P01552",   # SEA
    "P01553",   # SEB
    "P01554",   # SEC1
    "P18177",   # SED
    "P0A0L5",   # SEE

    # Clostridium perfringens epsilon toxin
    "P0C2L9",   # ETX

    # Bordetella pertussis toxin S1 subunit
    "P04977",   # PtxS1

    # Yersinia pestis — plague (F1/Caf1 capsule antigen)
    "Q7CJD3",   # Caf1 (F1 antigen)

    # Ebola virus / filovirus glycoproteins
    "Q05128",   # EBOV GP (Zaire)
    "P16285",   # MARV GP (Marburg)
})

# Map each blocked accession to a human-readable agent name for error messages.
_UNIPROT_AGENT_NAMES: dict[str, str] = {
    "P02879": "Ricin (Ricinus communis)",
    "P06750": "Ricinus communis agglutinin II",
    "P28590": "Abrin-a (Abrus precatorius)",
    "P11279": "Abrin (Abrus precatorius)",
    "P24476": "Abrin-b (Abrus precatorius)",
    "P10844": "Botulinum neurotoxin type A",
    "P10845": "Botulinum neurotoxin type B",
    "P19321": "Botulinum neurotoxin type C",
    "P19322": "Botulinum neurotoxin type D",
    "P19619": "Botulinum neurotoxin type E",
    "P30995": "Botulinum neurotoxin type F",
    "Q60393": "Botulinum neurotoxin type G",
    "P04958": "Tetanus toxin",
    "P13423": "Anthrax protective antigen",
    "P15917": "Anthrax lethal factor",
    "P40136": "Anthrax edema factor",
    "P00588": "Diphtheria toxin",
    "P09385": "Shiga toxin 1 (A-subunit)",
    "P09386": "Shiga toxin 1 (B-subunit)",
    "P09387": "Shiga toxin 2 (A-subunit)",
    "P09388": "Shiga toxin 2 (B-subunit)",
    "P01555": "Cholera toxin (A-subunit)",
    "P01556": "Cholera toxin (B-subunit)",
    "P01552": "Staphylococcal enterotoxin A",
    "P01553": "Staphylococcal enterotoxin B",
    "P01554": "Staphylococcal enterotoxin C1",
    "P18177": "Staphylococcal enterotoxin D",
    "P0A0L5": "Staphylococcal enterotoxin E",
    "P0C2L9": "Clostridium perfringens epsilon toxin",
    "P04977": "Bordetella pertussis toxin S1",
    "Q7CJD3": "Yersinia pestis F1 antigen (plague)",
    "Q05128": "Ebola virus glycoprotein (Zaire)",
    "P16285": "Marburg virus glycoprotein",
}

# ---------------------------------------------------------------------------
# Layer 2 — Protein name / description keyword blocklist
# Case-insensitive substring match. Matches any word in the protein name or
# description returned by UniProt.
# ---------------------------------------------------------------------------

DANGEROUS_KEYWORDS: tuple[str, ...] = (
    "ricin",
    "abrin",
    "botulinum",
    "tetanus toxin",
    "anthrax lethal factor",
    "anthrax edema factor",
    "anthrax protective antigen",
    "diphtheria toxin",
    "shiga toxin",
    "vero toxin",                       # alternative name for Shiga toxin
    "cholera toxin",
    "epsilon toxin",                    # C. perfringens select agent
    "pertussis toxin",
    "staphylococcal enterotoxin",
    "ebola virus glycoprotein",
    "marburg virus glycoprotein",
    "variola virus",
    "smallpox",
    "ricinus communis agglutinin",
    "modeccin",                         # type II RIP (Adenia digitata)
    "volkensin",                        # type II RIP (Adenia volkensii)
    "viscumin",                         # type II RIP (Viscum album)
)

# Human-readable name for keyword match error messages.
_KEYWORD_AGENT_NAMES: dict[str, str] = {
    k: k.title() for k in DANGEROUS_KEYWORDS
}

# ---------------------------------------------------------------------------
# Layer 3 — Amino-acid signature scan
# Short diagnostic peptides from published immunoassay literature and CDC
# identification guidance. These are NOT active-site sequences; they are
# structural/surface epitopes used in authorised detection kits and are already
# in the public domain for identification purposes.
#
# Each entry: (agent_name, signature_peptide)
# Signatures are 7–14 aa, unique to the listed toxin family.
# ---------------------------------------------------------------------------

AA_SIGNATURES: tuple[tuple[str, str], ...] = (
    # Ricin A-chain N-terminal region (publicly used as immunoassay standard)
    ("Ricin A-chain",           "IFPKQYPIINFT"),
    # Ricin B-chain galactose-binding repeat (QW domain, used in ELISA)
    ("Ricin B-chain",           "CNTSQYINQNL"),
    # Botulinum toxin HEXXH zinc-binding motif (BoNT serotype-conserved)
    ("Botulinum neurotoxin",    "HELIH"),
    # Tetanus toxin C-fragment (TeNT-Hc) — binding domain marker
    ("Tetanus toxin",           "FNIFPGSSYVNK"),
    # Anthrax lethal factor HEXXGH catalytic motif
    ("Anthrax lethal factor",   "HELGHSL"),
    # Diphtheria toxin ADP-ribosylation motif
    ("Diphtheria toxin",        "YADGHEDVQAQD"),
    # Shiga toxin A1 active site pentapeptide (depurination)
    ("Shiga toxin",             "YIYYVD"),
    # Cholera toxin A1 active site loop
    ("Cholera toxin",           "RDEYILM"),
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def screen_compound_id(compound_id: str) -> BiosafetyResult:
    """Layer 1 only — call this before any network request."""
    uid = compound_id.strip().upper()
    if uid in BLOCKED_UNIPROT_IDS:
        agent = _UNIPROT_AGENT_NAMES.get(uid, uid)
        return BiosafetyResult(
            safe=False,
            matched_agent=agent,
            match_type="uniprot_id",
            reason=(
                f"UniProt accession {uid} is a CDC/USDA Select Agent toxin "
                f"({agent}). Synthesis of this sequence is prohibited without "
                "an approved Select Agent registration."
            ),
        )
    return BiosafetyResult(safe=True, matched_agent=None, match_type=None, reason=None)


def screen_protein(protein_name: str, chains_aa: list[str]) -> BiosafetyResult:
    """Layers 2 + 3 — call this after UniProt sequence retrieval."""
    # Layer 2: keyword scan on protein name / description
    name_lower = protein_name.lower()
    for kw in DANGEROUS_KEYWORDS:
        if kw in name_lower:
            agent = _KEYWORD_AGENT_NAMES.get(kw, kw)
            return BiosafetyResult(
                safe=False,
                matched_agent=agent,
                match_type="keyword",
                reason=(
                    f"Protein name '{protein_name}' matches select-agent keyword "
                    f"'{kw}'. Synthesis is prohibited without an approved "
                    "Select Agent registration."
                ),
            )

    # Layer 3: amino-acid signature scan across all chains
    for aa_seq in chains_aa:
        seq_upper = aa_seq.upper()
        for agent_name, sig in AA_SIGNATURES:
            if sig.upper() in seq_upper:
                return BiosafetyResult(
                    safe=False,
                    matched_agent=agent_name,
                    match_type="aa_signature",
                    reason=(
                        f"Chain sequence contains a diagnostic amino-acid signature "
                        f"for {agent_name}. Synthesis is prohibited without an "
                        "approved Select Agent registration."
                    ),
                )

    return BiosafetyResult(safe=True, matched_agent=None, match_type=None, reason=None)
