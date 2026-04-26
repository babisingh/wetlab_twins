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
```

## Pipeline overview

The pipeline is a two-level LangGraph `StateGraph`. Each level is independently compiled and testable.

```
Main graph
├── sequence_retrieval          fetch & parse UniProt, extract mature chains (deterministic)
├── human_checkpoint_chains     review extracted chains before processing
├── chain_processing × N        fan-out: one subgraph run per protein chain
│   └── Chain subgraph
│       ├── codon_optimization  greedy best-codon assignment + restriction avoidance
│       ├── cassette_assembly   ATG + 6×His tag + protease site + gene + stop
│       ├── plasmid_assembly    insert cassette into pET-28a(+) vector
│       ├── sequence_validation 7 independent checks (see below)
│       └── remediation loop   synonymous swaps → revalidate (up to 3 attempts)
│           └── escalation_agent (LLM, optional) fires once if loop exhausted
├── merge_results               collect all chain outcomes
├── human_checkpoint_plasmid    review assembled plasmids
└── structural_validation       (optional) ESMFold/AlphaFold confidence check
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
