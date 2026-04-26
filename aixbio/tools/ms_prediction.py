from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

_PROTON = 1.007276  # Da


def predict_ms_peptides(chain_id: str, aa_sequence: str) -> list[dict]:
    """Tryptic digest + monoisotopic mass prediction for LC-MS/MS QC.

    Returns a list of dicts, one per peptide, ready for TSV output.
    Columns: chain_id, peptide, start, end, length, mono_mass, mz_2, mz_3.

    Uses Trypsin rule: cleave C-terminal to K or R, not before P.
    Missed cleavages: 0 (standard LC-MS/MS experiment assumption).
    """
    try:
        from pyteomics import mass, parser
    except ImportError:
        logger.error("pyteomics is not installed; run 'uv sync' to add it.")
        return []

    rule = parser.expasy_rules["trypsin"]
    peptides = list(parser.cleave(aa_sequence, rule, missed_cleavages=0))

    # Track start positions
    rows = []
    pos = 1  # 1-based
    for pep in peptides:
        if not pep:
            continue
        pep_len = len(pep)
        try:
            mono = mass.calculate_mass(sequence=pep)
        except Exception:
            mono = 0.0
        rows.append({
            "chain_id": chain_id,
            "peptide": pep,
            "start": pos,
            "end": pos + pep_len - 1,
            "length": pep_len,
            "mono_mass": round(mono, 4),
            "mz_2": round((mono + 2 * _PROTON) / 2, 4),
            "mz_3": round((mono + 3 * _PROTON) / 3, 4),
        })
        pos += pep_len

    logger.info("MS prediction [%s]: %d tryptic peptides", chain_id, len(rows))
    return rows


def format_tsv(rows: list[dict]) -> str:
    if not rows:
        return ""
    header = "chain_id\tpeptide\tstart\tend\tlength\tmono_mass\tmz_2+\tmz_3+\n"
    lines = [
        f"{r['chain_id']}\t{r['peptide']}\t{r['start']}\t{r['end']}\t"
        f"{r['length']}\t{r['mono_mass']}\t{r['mz_2']}\t{r['mz_3']}"
        for r in rows
    ]
    return header + "\n".join(lines) + "\n"
