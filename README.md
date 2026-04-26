# aixbio

A digital-to-biological pipeline that converts a UniProt protein accession into a production-ready, expression-validated plasmid sequence optimized for *E. coli*. Built with [LangGraph](https://github.com/langchain-ai/langgraph) as a hackathon project.

```
UniProt accession → DNA → cassette → plasmid → validated GenBank file
```

## Quick start

```bash
# Install (Python 3.11, uses uv)
uv sync

# Create a .env file with your OpenRouter key
echo "OPENROUTER_API_KEY=sk-..." > .env

# Run — insulin as an example
uv run aixbio P01308                   # interactive (human review prompts)
uv run aixbio P01308 --auto-approve    # fully automated
uv run aixbio P01308 --escalation      # enable LLM fallback on hard failures
uv run aixbio P01308 --structural      # add ESMFold/AlphaFold structure check
uv run aixbio P01308 --protocol        # generate literature-backed wet-lab SOP
```

## Pipeline overview

The pipeline is a two-level LangGraph `StateGraph`. Each level is independently compiled and testable.

```
Main graph
├── sequence_retrieval          fetch & parse UniProt, extract mature chains (deterministic)
├── host_selection              recommend E. coli / P. pastoris / CHO / Sf9 from sequence features
├── human_checkpoint_chains     review chains + host recommendation before processing
├── chain_processing × N        fan-out: one subgraph per chain
│   └── Chain subgraph
│       ├── solubility_prediction  composition heuristic: inclusion body / disulfide risk score
│       ├── codon_optimization     greedy best-codon assignment + restriction site avoidance
│       ├── cassette_assembly      ATG + 6×His tag + Enterokinase site + gene + stop
│       ├── plasmid_assembly       insert cassette into pET-28a(+) vector
│       ├── sequence_validation    7 independent checks (GC, CAI, sites, RNA, identity, rare codons, repeats)
│       └── remediation loop       synonymous swaps → revalidate (up to 3 attempts)
│           └── escalation_agent (LLM, optional) fires once if loop exhausted
├── merge_results               collect all chain outcomes + aggregate validation
├── human_checkpoint_plasmid    review assembled plasmids before export
├── structural_validation       (optional) ESMFold/AlphaFold confidence check
└── protocol_generation         (optional) PubMed fetch + LLM → literature-backed wet-lab SOP
```

## Modularity

The codebase is split into six focused layers — each independently importable and testable:

| Layer | Path | Responsibility |
|---|---|---|
| **models** | `aixbio/models/` | Frozen dataclasses for every inter-step artefact. Discriminated unions on `kind` for escalation outcomes. No logic. |
| **state** | `aixbio/state/` | `TypedDict` state definitions with LangGraph `Annotated` reducers for list fields. |
| **tools** | `aixbio/tools/` | Pure functions: codon tables, CAI, GC content, restriction site scanning, RNA fold estimation, GenBank generation, UniProt/AlphaFold API clients. |
| **nodes** | `aixbio/nodes/` | LangGraph node functions. Each reads from state and returns only the keys it modifies. |
| **graph** | `aixbio/graph/` | Graph wiring only — no business logic. |
| **prompts** | `aixbio/prompts/` | LLM prompt templates decoupled from node logic. |

Deterministic nodes (steps 2–5, remediation) are fully testable without any LLM or network call. Only the optional `escalation_agent` touches the LLM; `sequence_retrieval` is deterministic (UniProt REST + feature parsing).

## Validation & verification

### Automated sequence checks (`sequence_validation`)

Every chain is evaluated against seven independent criteria before the result is accepted:

| Check | Threshold | Scope |
|---|---|---|
| **GC content** | 0.50 – 0.60 | Full expression construct |
| **CAI score** | > 0.80 | Coding gene only |
| **Restriction sites** | 0 internal hits | Full construct (cloning flanks excluded) |
| **5' RNA secondary structure** | ΔG > −10 kcal/mol | First ~50 nt of construct (ViennaRNA) |
| **Back-translation identity** | 100% match | Gene back-translated vs original amino acid sequence |
| **Rare codons** | 0 | Coding gene |
| **Direct repeats** | 0 repeats ≥ 20 bp | Full construct (RecA-independent deletion risk) |

All checks are defined as typed `CheckResult` dataclasses and collected into a `ChainValidation` record that travels with the chain through the rest of the pipeline.

### Deterministic remediation loop

When validation fails, the pipeline does not give up — it plans targeted synonymous codon swaps for each failed check and retries (default: up to 3 attempts):

- **Restriction sites** — swap one codon inside the recognition sequence until the site disappears
- **Rare codons** — replace with the highest-frequency synonymous alternative
- **CAI score** — upgrade suboptimal codons sorted by relative frequency (worst first)
- **GC content** — shift GC up or down codon-by-codon until in range

After each round of fixes the cassette and plasmid are fully reassembled and all seven checks are re-run.

### LLM escalation agent (opt-in)

If the deterministic loop exhausts its retry budget, an LLM escalation agent fires exactly once per chain (enabled with `--escalation`). It inspects the full remediation history and returns a discriminated-union decision:

| Outcome | Action |
|---|---|
| `apply_plan` | Try specific codon swaps the greedy loop missed |
| `incompatible` | Host cannot express this protein; mark chain as incompatible |
| `change_strategy` | Change tag / protease site / vector / cloning sites and re-run |
| `give_up` | Unfixable by this pipeline; surface the failure with detail |

The `escalation_used` flag prevents re-entry so the LLM fires at most once per chain.

### Human checkpoints

Two optional review gates pause the pipeline for human approval:

1. **After sequence retrieval** — inspect extracted chains before committing compute
2. **After plasmid assembly** — review the final constructs before export

Both gates are bypassed with `--auto-approve`.

### Solubility and inclusion body prediction

Every chain is assessed for E. coli expression solubility before codon optimisation begins. The
deterministic heuristic (Biopython `ProteinAnalysis`, no network call) scores four features:

| Feature | Direction | Weight |
|---|---|---|
| GRAVY score | lower = more soluble | 35% |
| Instability index | < 40 = more stable | 30% |
| Charged residue fraction (D/E/K/R) | higher = more soluble | 20% |
| pI deviation from pH 7 | closer to 7 = more soluble | 15% |

Score 0.0–1.0; **≥ 0.45 = predicted soluble** (Protein-Sol calibration, Hebditch et al. 2017).
Any chain with ≥ 2 cysteines is independently flagged for **disulfide risk** regardless of score.
Warnings and tag recommendations (MBP, SUMO) are surfaced in the pipeline output and carried into
the protocol if `--protocol` is used.

### Literature-backed wet-lab SOP (`--protocol`)

Add `--protocol` to generate a full bench SOP after the pipeline completes. The node:

1. Queries PubMed E-utilities (no API key needed) for the top 5 expression papers for this protein.
2. Compiles structured context: protein metadata, solubility scores, validation metrics,
   tag/vector/host configuration.
3. Calls the LLM to synthesise a numbered, literature-cited SOP covering transformation,
   induction conditions, lysis, IMAC purification, Enterokinase tag removal, and — for
   disulfide-rich proteins like insulin — a dedicated inclusion body refolding section.

Output is written to `output/protocol.md` alongside the GenBank and FASTA files.

### Audit trail

Every decision is recorded in the shared state:

- `decision_log` — structured `AgentDecision` records from every LLM node
- `remediation_history` — every `RemediationAction` applied (before/after codons, positions, reasoning)
- `warnings` — non-fatal issues surfaced during processing

## Configuration

Defaults are set in `config.py` and all overridable via environment variables:

| Setting | Default |
|---|---|
| Host organism | *Escherichia coli* |
| Affinity tag | 6×His |
| Protease site | Enterokinase |
| Expression vector | pET-28a(+) |
| Cloning sites | BamHI / XhoI |
| Max remediation attempts | 3 |
| LLM model | `deepseek/deepseek-v4-flash` (via `LLM_MODEL`) |

## Running tests

```bash
uv run pytest                                       # full suite
uv run pytest tests/test_deterministic_nodes.py     # steps 2–5, no LLM required
uv run pytest tests/test_tools.py                   # pure utility functions
uv run pytest tests/test_chain_subgraph.py          # LangGraph subgraph integration
uv run pytest tests/test_escalation.py              # escalation agent (mocked LLM)
```

## Environment

Requires `OPENROUTER_API_KEY` in `.env`. Set `LLM_MODEL` to use a different model via OpenRouter.
Copy `.env.example` to `.env` and fill in your key.

---

## Insulin demo — full user journey

This section walks through the complete pipeline for human insulin (P01308) end to end.

### 1. Setup

```bash
git clone https://github.com/babisingh/wetlab_twins
cd wetlab_twins
uv sync
cp .env.example .env          # add your OPENROUTER_API_KEY
```

### 2. Run (fully automated, with protocol)

```bash
uv run aixbio P01308 --auto-approve --protocol
```

### 3. Annotated terminal output

```
Starting pipeline for P01308 targeting Escherichia coli
============================================================

# Step 1 — UniProt fetch + chain extraction (deterministic)
# Parses "Chain" feature annotations to get mature A and B chains.

HUMAN REVIEW: chain_extraction_review          ← --auto-approve bypasses this
============================================================
{
  "protein": "Insulin",
  "uniprot_id": "P01308",
  "chains_extracted": [
    {"id": "Insulin_A_chain", "length": 21, "sequence_preview": "GIVEQCCTSICSLYQLENYCN"},
    {"id": "Insulin_B_chain", "length": 30, "sequence_preview": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT..."}
  ],
  "host_recommendation": {
    "primary_host": "Escherichia coli",
    "confidence": "medium",
    "reasoning": "Short protein (51 aa) with 6 cysteines. E. coli expression as
                  inclusion bodies followed by oxidative refolding is the
                  established route (e.g. recombinant insulin).",
    "alternative_hosts": ["Pichia pastoris"],
    "caveats": [
      "Denaturing lysis (urea/guanidinium) required.",
      "Refolding yield is protein-dependent; optimise glutathione redox ratio."
    ]
  }
}
Auto-approving (--auto-approve flag)

# Steps 2–5 run in parallel for Chain A and Chain B:

# Chain A — solubility_prediction
# [Insulin_A_chain] Solubility score 0.42 < 0.45 — inclusion body formation likely.
# [Insulin_A_chain] 4 cysteine residues — disulfide refolding required.

# Chain B — solubility_prediction
# [Insulin_B_chain] 2 cysteine residues — disulfide refolding required.

# Both chains pass all 7 validation checks after codon optimisation.

HUMAN REVIEW: plasmid_review                   ← --auto-approve bypasses this
============================================================
{
  "chain_results": [
    {"chain_id": "Insulin_A_chain", "validation_passed": true,
     "remediation_rounds": 0, "insert_size": 171},
    {"chain_id": "Insulin_B_chain", "validation_passed": true,
     "remediation_rounds": 0, "insert_size": 198}
  ]
}
Auto-approving (--auto-approve flag)

# --protocol: PubMed fetches top 5 E. coli insulin expression papers,
# then LLM synthesises the SOP with PMID citations.

============================================================
PIPELINE RESULTS
============================================================
Status: completed

Host Recommendation: Escherichia coli [medium confidence]
  Reasoning: Short protein (51 aa) with 6 cysteines...
  Alternatives: Pichia pastoris
  ! Denaturing lysis (urea/guanidinium) required.
  ! Refolding yield is protein-dependent; optimise glutathione redox ratio.

Warnings (4):
  - [Insulin_A_chain] Solubility score 0.42 — inclusion body formation likely...
  - [Insulin_A_chain] 4 cysteine residues — disulfide refolding required.
  - [Insulin_B_chain] 2 cysteine residues — disulfide refolding required.
  - [host_selection] Denaturing lysis (urea/guanidinium) required.

Validation: PASSED
  Chain Insulin_A_chain: PASS
    [+] gc_content: 0.544 (threshold: 0.50-0.60)
    [+] cai_score: 0.962 (threshold: >0.80)
    [+] restriction_sites: none (threshold: 0 internal hits)
    [+] rna_secondary_structure: -4.2 kcal/mol (threshold: >-10)
    [+] back_translation: match (threshold: 100%)
    [+] rare_codons: none (threshold: 0)
    [+] direct_repeats: none (threshold: 0 repeats >= 20 bp)
  Chain Insulin_B_chain: PASS
    [+] gc_content: 0.551 ...
    [+] cai_score: 0.971 ...
    ...

Chain Results (2):
  Insulin_A_chain:
    Insert size: 171 bp
    Remediation rounds: 0
    Validation passed: True
    Solubility: 0.42 [INCLUSION BODY RISK]
    Disulfide risk: YES (refolding required)
  Insulin_B_chain:
    Insert size: 198 bp
    Remediation rounds: 0
    Validation passed: True
    Solubility: 0.53 [SOLUBLE]
    Disulfide risk: YES (refolding required)

Protocol generated (4821 chars) — written to output/protocol.md

  Wrote output/Insulin_A_chain.fasta
  Wrote output/Insulin_A_chain_plasmid.gb
  Wrote output/Insulin_B_chain.fasta
  Wrote output/Insulin_B_chain_plasmid.gb
  Wrote output/Insulin_A_chain_peptides.tsv (4 tryptic peptides)
  Wrote output/Insulin_B_chain_peptides.tsv (3 tryptic peptides)
  Wrote output/P01308_summary.json

All artifacts written to output/
```

### 4. Output files explained

| File | Contents |
|---|---|
| `Insulin_A_chain.fasta` | E. coli codon-optimised DNA for the mature A-chain |
| `Insulin_B_chain.fasta` | E. coli codon-optimised DNA for the mature B-chain |
| `Insulin_A_chain_plasmid.gb` | Annotated GenBank record: pET-28a(+) + ATG + 6×His + EK site + A-chain + stop |
| `Insulin_B_chain_plasmid.gb` | Same for B-chain |
| `Insulin_A_chain_peptides.tsv` | Tryptic peptide masses for A-chain (LC-MS/MS QC reference) |
| `Insulin_B_chain_peptides.tsv` | Tryptic peptide masses for B-chain |
| `synthesis_report.txt` | IDT and Twist synthesis feasibility + price estimates for every chain |
| `P01308_summary.json` | Machine-readable summary: validation checks, solubility scores, host recommendation, synthesis quotes |
| `protocol.md` | Full bench SOP with literature citations (only with `--protocol`) |

### 5. Synthesis feasibility (IDT + Twist)

Before the bench team orders a gene fragment, the pipeline automatically screens every cassette against the published complexity rules and 2024 list prices for two vendors:

| Vendor | Min length | Max length | GC window | Homopolymer limit | Price (insulin B) |
|--------|-----------|-----------|-----------|-------------------|--------------------|
| IDT gBlocks | 125 bp | 3 000 bp | 25–65% | A/T ≤ 9 bp · G/C ≤ 7 bp | ~$95 |
| Twist Gene Fragments | 300 bp | 5 000 bp | 25–65% | any base ≤ 9 bp | flagged (too short) |

For insulin (chains ~171–198 bp with BamHI/XhoI flanks):

```
Synthesis Feasibility Report
============================================

Chain: Insulin_B_chain
  Insert length : 216 bp
  GC content    : 52.3%
  Max homopolymer: 8 bp

  IDT     [FEASIBLE]  ~$95  10–15 business days
           ℹ Expedited synthesis available (3–5 business days, additional cost).
  Twist   [FLAGGED]   N/A   N/A
           ✗ Sequence too short (216 bp); Twist minimum is 300 bp.
             Use IDT gBlocks (min 125 bp) or oligo assembly for short inserts.
```

The `synthesis_report.txt` is written to `output/` automatically; no extra flags needed.

### 7. Peptide mass table (B-chain excerpt)

```
chain_id        peptide                   start  end  length  mono_mass    mz_2+      mz_3+
Insulin_B_chain FVNQHLCGSHLVEALYLVCGER    1      22   22      2431.1742    1216.5943  811.3990
Insulin_B_chain GFFYTPK                   23     29   7       837.4198     419.7171   280.1490
Insulin_B_chain T                         30     30   1       119.0582     60.5363    40.6932
```

This table lets a mass spectrometrist verify protein identity after tryptic digest without running the full experiment first.

### 8. Generated SOP excerpt (`protocol.md`)

```markdown
# Protocol: Insulin Expression and Purification in E. coli

## Literature Context
Prior studies (PMID:2063194, PMID:7592457) established that recombinant human insulin
chains can be expressed in E. coli BL21(DE3) as inclusion bodies and refolded in vitro
with high yield using a glutathione redox system...

## 1. Strains, Plasmid, and Safety Notes
- **Host:** E. coli BL21(DE3)
- **Vector:** pET-28a(+) with BamHI/XhoI cloning sites
- **Insert A:** Insulin A-chain (21 aa) with N-terminal 6×His tag + Enterokinase site
- **Insert B:** Insulin B-chain (30 aa) with N-terminal 6×His tag + Enterokinase site

## 2. Transformation
a. Thaw 50 µL competent BL21(DE3) cells on ice (5 min).
b. Add 1 µL plasmid DNA (10–50 ng). Mix gently.
c. Heat shock: 42°C, 30 s. Return to ice 2 min.
d. Add 950 µL SOC medium. Recover 37°C, 250 rpm, 1 h.
e. Plate 100 µL on LB + kanamycin (50 µg/mL). Incubate 37°C overnight.
...

## 8. Inclusion Body Refolding
⚠ Required: both chains have disulfide bonds that will not form in the E. coli cytoplasm.

a. Resuspend inclusion body pellet in denaturing buffer:
   6 M guanidinium HCl, 50 mM Tris-HCl pH 8.0, 1 mM EDTA, 100 mM DTT.
b. Stir at RT for 2 h to fully solubilise.
c. Clarify by centrifugation 20 000 × g, 30 min, 4°C.
d. Dialyse chain A and chain B separately against refolding buffer:
   0.1 M Tris-HCl pH 8.0, 0.5 M L-arginine, 1 mM EDTA,
   2 mM reduced glutathione (GSH), 0.2 mM oxidised glutathione (GSSG).
e. Combine chain A and B in a 1:1 molar ratio. Stir 16 h at 4°C.
f. Monitor refolding by RP-HPLC or native PAGE.
...
```

### 9. What happens without `--protocol`

The pipeline still produces validated GenBank files, FASTA sequences, a peptide mass table, and the
JSON summary — all with zero LLM calls beyond the initial UniProt fetch. Add `--protocol` only when
you are ready to hand the design to the bench.
