from __future__ import annotations

import re

from aixbio.models.synthesis import SynthesisQuote, VendorQuote
from aixbio.tools.gc import compute_gc

# Restriction site flanks added by pET-28a(+) cloning
_BAMHI = "GGATCC"
_XHOI = "CTCGAG"

# --------------------------------------------------------------------------
# Public pricing tiers (USD, 2024 list prices — subject to change)
# IDT gBlocks Gene Fragments: https://www.idtdna.com/pages/products/genes-and-gene-fragments/gene-fragments
# Twist Gene Fragments: https://www.twistbioscience.com/products/genes
# --------------------------------------------------------------------------
_IDT_TIERS = [
    (125, 500,  95.0,  "10–15 business days"),
    (500, 1000, 149.0, "10–15 business days"),
    (1000, 2000, 249.0, "10–15 business days"),
    (2000, 3000, 399.0, "10–15 business days"),
]

_TWIST_TIERS = [
    (300, 500,  99.0,  "10–14 business days"),
    (500, 1000, 149.0, "10–14 business days"),
    (1000, 5000, 249.0, "10–14 business days"),
]

# IDT homopolymer thresholds (public complexity guidelines)
_IDT_AT_HOMOPOLYMER_MAX = 9   # run of 10+ A or T → flag
_IDT_GC_HOMOPOLYMER_MAX = 7   # run of  8+ G or C → flag
# Twist uses a single threshold for all bases
_TWIST_HOMOPOLYMER_MAX = 9    # run of 10+ any base → flag

_LOCAL_GC_WINDOW = 50
_LOCAL_GC_MIN = 0.20
_LOCAL_GC_MAX = 0.80


def get_synthesis_quotes(chain_id: str, cassette_dna: str) -> SynthesisQuote:
    """Evaluate a cassette sequence for synthesis feasibility at IDT and Twist.

    The insert submitted to the vendor is the cassette flanked by the
    pET-28a(+) BamHI and XhoI cloning sites.
    """
    insert = _BAMHI + cassette_dna + _XHOI
    gc = compute_gc(insert)
    longest_hp = _longest_homopolymer(insert)

    idt = _idt_quote(insert, gc, longest_hp)
    twist = _twist_quote(insert, gc, longest_hp)

    return SynthesisQuote(
        chain_id=chain_id,
        sequence=insert,
        sequence_length=len(insert),
        gc_content=round(gc, 3),
        longest_homopolymer=longest_hp,
        quotes=(idt, twist),
    )


# ---------------------------------------------------------------------------
# Vendor-specific checkers
# ---------------------------------------------------------------------------

def _idt_quote(seq: str, gc: float, longest_hp: int) -> VendorQuote:
    length = len(seq)
    flags: list[str] = []
    notes: list[str] = []

    if length < 125:
        flags.append(
            f"Sequence too short ({length} bp); IDT gBlocks minimum is 125 bp. "
            "Consider ordering as two overlapping oligos and assembling by PCR."
        )
    elif length > 3000:
        flags.append(f"Sequence too long ({length} bp); IDT gBlocks maximum is 3000 bp.")

    if gc < 0.25:
        flags.append(f"GC content {gc:.0%} is below IDT minimum of 25%.")
    elif gc > 0.65:
        flags.append(f"GC content {gc:.0%} exceeds IDT maximum of 65%.")

    at_run = _longest_homopolymer_bases(seq, "AT")
    gc_run = _longest_homopolymer_bases(seq, "GC")
    if at_run > _IDT_AT_HOMOPOLYMER_MAX:
        flags.append(f"A/T homopolymer run of {at_run} bp exceeds IDT limit of {_IDT_AT_HOMOPOLYMER_MAX + 1}.")
    if gc_run > _IDT_GC_HOMOPOLYMER_MAX:
        flags.append(f"G/C homopolymer run of {gc_run} bp exceeds IDT limit of {_IDT_GC_HOMOPOLYMER_MAX + 1}.")

    local_flags = _local_gc_flags(seq)
    flags.extend(local_flags)

    if 125 <= length <= 500:
        notes.append("Expedited synthesis available (3–5 business days, additional cost).")

    feasible = len(flags) == 0
    cost, turnaround = _lookup_price(_IDT_TIERS, length, feasible)

    return VendorQuote(
        vendor="IDT",
        feasible=feasible,
        rejection_flags=tuple(flags),
        estimated_cost_usd=cost,
        estimated_turnaround=turnaround,
        notes=tuple(notes),
    )


def _twist_quote(seq: str, gc: float, longest_hp: int) -> VendorQuote:
    length = len(seq)
    flags: list[str] = []
    notes: list[str] = []

    if length < 300:
        flags.append(
            f"Sequence too short ({length} bp); Twist minimum is 300 bp. "
            "Use IDT gBlocks (min 125 bp) or oligo assembly for short inserts."
        )
    elif length > 5000:
        flags.append(f"Sequence too long ({length} bp); Twist maximum is 5000 bp.")

    if gc < 0.25:
        flags.append(f"GC content {gc:.0%} is below Twist minimum of 25%.")
    elif gc > 0.65:
        flags.append(f"GC content {gc:.0%} exceeds Twist maximum of 65%.")

    if longest_hp > _TWIST_HOMOPOLYMER_MAX:
        flags.append(
            f"Homopolymer run of {longest_hp} bp exceeds Twist limit of {_TWIST_HOMOPOLYMER_MAX + 1}."
        )

    local_flags = _local_gc_flags(seq)
    flags.extend(local_flags)

    feasible = len(flags) == 0
    cost, turnaround = _lookup_price(_TWIST_TIERS, length, feasible)

    return VendorQuote(
        vendor="Twist",
        feasible=feasible,
        rejection_flags=tuple(flags),
        estimated_cost_usd=cost,
        estimated_turnaround=turnaround,
        notes=tuple(notes),
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _longest_homopolymer(seq: str) -> int:
    return max(
        (max((len(m) for m in re.findall(f"{b}+", seq.upper())), default=0) for b in "ACGT"),
        default=0,
    )


def _longest_homopolymer_bases(seq: str, bases: str) -> int:
    return max(
        (max((len(m) for m in re.findall(f"{b}+", seq.upper())), default=0) for b in bases),
        default=0,
    )


def _local_gc_flags(seq: str) -> list[str]:
    flags = []
    for i in range(len(seq) - _LOCAL_GC_WINDOW + 1):
        window = seq[i: i + _LOCAL_GC_WINDOW]
        gc = sum(window.upper().count(b) for b in "GC") / _LOCAL_GC_WINDOW
        if gc < _LOCAL_GC_MIN:
            flags.append(
                f"Low-GC window ({gc:.0%}) at positions {i + 1}–{i + _LOCAL_GC_WINDOW} "
                f"(threshold: {_LOCAL_GC_MIN:.0%})."
            )
            break
        if gc > _LOCAL_GC_MAX:
            flags.append(
                f"High-GC window ({gc:.0%}) at positions {i + 1}–{i + _LOCAL_GC_WINDOW} "
                f"(threshold: {_LOCAL_GC_MAX:.0%})."
            )
            break
    return flags


def _lookup_price(
    tiers: list[tuple], length: int, feasible: bool
) -> tuple[float | None, str | None]:
    if not feasible:
        return None, None
    for lo, hi, price, turnaround in tiers:
        if lo <= length <= hi:
            return price, turnaround
    return None, None


def format_quotes_text(quotes: list[SynthesisQuote]) -> str:
    lines = ["Synthesis Feasibility Report", "=" * 44, ""]
    for q in quotes:
        lines.append(f"Chain: {q.chain_id}")
        lines.append(f"  Insert length : {q.sequence_length} bp")
        lines.append(f"  GC content    : {q.gc_content:.1%}")
        lines.append(f"  Max homopolymer: {q.longest_homopolymer} bp")
        lines.append("")
        for vq in q.quotes:
            status = "FEASIBLE" if vq.feasible else "FLAGGED"
            cost = f"${vq.estimated_cost_usd:.0f}" if vq.estimated_cost_usd else "N/A"
            eta = vq.estimated_turnaround or "N/A"
            lines.append(f"  {vq.vendor:6s}  [{status}]  ~{cost}  {eta}")
            for flag in vq.rejection_flags:
                lines.append(f"           ✗ {flag}")
            for note in vq.notes:
                lines.append(f"           ℹ {note}")
        lines.append("")
    return "\n".join(lines)
