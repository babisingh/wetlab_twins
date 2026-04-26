from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class VendorQuote:
    vendor: str                          # "IDT" | "Twist"
    feasible: bool
    rejection_flags: tuple[str, ...]     # reasons the vendor would reject the sequence
    estimated_cost_usd: float | None
    estimated_turnaround: str | None     # e.g. "10–15 business days"
    notes: tuple[str, ...]               # informational (not blocking)


@dataclass(frozen=True)
class SynthesisQuote:
    chain_id: str
    sequence: str                        # insert sequence submitted to vendor
    sequence_length: int
    gc_content: float
    longest_homopolymer: int             # longest single-base run in insert
    quotes: tuple[VendorQuote, ...]
