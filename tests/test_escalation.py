"""Tests for the escalation agent: parser, router, and subgraph integration."""
from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

from aixbio.graph.chain_subgraph import compile_chain_subgraph
from aixbio.models.dna import DNAChain
from aixbio.models.escalation import (
    EscalationApplyPlan,
    EscalationChangeStrategy,
    EscalationGiveUp,
    EscalationIncompatible,
)
from aixbio.models.protein import Chain, ProteinRecord
from aixbio.models.validation import CheckResult
from aixbio.nodes.escalation_agent import _parse_and_validate, _strip_codefence
from aixbio.nodes.routers import escalation_router, validation_router


INSULIN_B = Chain(
    id="Insulin_B",
    aa_sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    length=30,
)

PROTEIN = ProteinRecord(uniprot_id="P01308", name="Insulin", chains=(INSULIN_B,))


def _make_dna_chain() -> DNAChain:
    from aixbio.tools.codon_tables import best_ecoli_codon
    from aixbio.tools.cai import compute_cai
    from aixbio.tools.gc import compute_gc

    dna = "".join(best_ecoli_codon(aa) for aa in INSULIN_B.aa_sequence)
    return DNAChain(
        id="Insulin_B",
        dna_sequence=dna,
        cai_score=compute_cai(dna),
        gc_content=compute_gc(dna),
    )


def _make_input_state(**overrides):
    state = {
        "chain": INSULIN_B,
        "host_organism": "E. coli K12",
        "avoid_sites": ("BamHI", "XhoI", "EcoRI"),
        "tag_type": "6xHis",
        "protease_site": "Enterokinase",
        "vector": "pET-28a(+)",
        "cloning_sites": ("BamHI", "XhoI"),
        "protein_record": PROTEIN,
        "optimized_dna": None,
        "cassette": None,
        "plasmid": None,
        "chain_validation": None,
        "remediation_attempt": 0,
        "max_remediation_attempts": 3,
        "failed_checks": (),
        "remediation_plan": None,
        "chain_results": [],
        "remediation_history": [],
        "enable_escalation": False,
        "escalation_used": False,
        "escalation_decision": None,
        "decision_log": [],
        "warnings": [],
    }
    state.update(overrides)
    return state


# ---------------------------------------------------------------------------
# Unit: _parse_and_validate
# ---------------------------------------------------------------------------

class TestParseAndValidate:
    def test_apply_plan_valid(self):
        dna = _make_dna_chain()
        from aixbio.tools.codon_tables import split_codons, synonymous_alternatives

        codons = split_codons(dna.dna_sequence)
        alt = synonymous_alternatives(codons[0])
        raw = json.dumps({
            "kind": "apply_plan",
            "actions": [{
                "check_name": "cai",
                "strategy": "swap",
                "target_positions": [0],
                "replacement_codons": [alt[0]],
            }],
            "reasoning": "test",
            "diagnosis": "test",
        })
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationApplyPlan)
        assert len(result.actions) == 1
        assert result.actions[0].target_positions == (0,)

    def test_apply_plan_rejects_nonsynonymous(self):
        dna = _make_dna_chain()
        from aixbio.tools.codon_tables import split_codons, translate_codon

        codons = split_codons(dna.dna_sequence)
        original_aa = translate_codon(codons[0])
        bad_codon = "ATG" if original_aa != "M" else "GGG"
        raw = json.dumps({
            "kind": "apply_plan",
            "actions": [{
                "check_name": "cai",
                "target_positions": [0],
                "replacement_codons": [bad_codon],
            }],
            "reasoning": "bad",
            "diagnosis": "bad",
        })
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationApplyPlan)
        assert len(result.actions) == 0

    def test_incompatible(self):
        dna = _make_dna_chain()
        raw = json.dumps({
            "kind": "incompatible",
            "reason": "requires glycosylation",
            "suggested_action": "use CHO cells",
        })
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationIncompatible)
        assert "glycosylation" in result.reason

    def test_change_strategy(self):
        dna = _make_dna_chain()
        raw = json.dumps({
            "kind": "change_strategy",
            "field": "tag_type",
            "new_value": "GST",
            "reason": "6xHis interferes",
        })
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationChangeStrategy)
        assert result.field == "tag_type"
        assert result.new_value == "GST"

    def test_give_up(self):
        dna = _make_dna_chain()
        raw = json.dumps({"kind": "give_up", "diagnosis": "pipeline bug"})
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationGiveUp)
        assert result.diagnosis == "pipeline bug"

    def test_codefence_stripped(self):
        dna = _make_dna_chain()
        raw = '```json\n{"kind": "give_up", "diagnosis": "test"}\n```'
        result = _parse_and_validate(raw, dna)
        assert isinstance(result, EscalationGiveUp)


# ---------------------------------------------------------------------------
# Unit: _strip_codefence
# ---------------------------------------------------------------------------

class TestStripCodefence:
    def test_plain(self):
        assert _strip_codefence('{"a": 1}') == '{"a": 1}'

    def test_json_fence(self):
        assert _strip_codefence('```json\n{"a": 1}\n```') == '{"a": 1}'

    def test_bare_fence(self):
        assert _strip_codefence('```\n{"a": 1}\n```') == '{"a": 1}'


# ---------------------------------------------------------------------------
# Unit: escalation_router
# ---------------------------------------------------------------------------

class TestEscalationRouter:
    def test_give_up(self):
        state = _make_input_state(
            escalation_decision=EscalationGiveUp(diagnosis="bug"),
        )
        assert escalation_router(state) == "halt_pipeline"

    def test_none_decision(self):
        state = _make_input_state(escalation_decision=None)
        assert escalation_router(state) == "halt_pipeline"

    def test_incompatible(self):
        state = _make_input_state(
            escalation_decision=EscalationIncompatible(reason="no glyco"),
        )
        assert escalation_router(state) == "package_result_escalated"

    def test_apply_plan_with_actions(self):
        from aixbio.models.remediation import PlannedFix
        fix = PlannedFix(
            check_name="cai",
            strategy="swap",
            target_positions=(0,),
            replacement_codons=("GGG",),
        )
        state = _make_input_state(
            escalation_decision=EscalationApplyPlan(actions=(fix,)),
        )
        assert escalation_router(state) == "apply_fixes"

    def test_apply_plan_empty_actions(self):
        state = _make_input_state(
            escalation_decision=EscalationApplyPlan(actions=()),
        )
        assert escalation_router(state) == "package_result_escalated"

    def test_change_strategy(self):
        state = _make_input_state(
            escalation_decision=EscalationChangeStrategy(
                field="tag_type", new_value="GST",
            ),
        )
        assert escalation_router(state) == "codon_optimization"


# ---------------------------------------------------------------------------
# Unit: validation_router escalation gating
# ---------------------------------------------------------------------------

class TestValidationRouterEscalation:
    def _failed_state(self, **kw):
        from aixbio.models.validation import ChainValidation
        checks = (CheckResult(name="cai", passed=False, value=0.5, threshold=">0.8"),)
        return _make_input_state(
            chain_validation=ChainValidation(id="B", passed=False, checks=checks),
            remediation_attempt=3,
            max_remediation_attempts=3,
            **kw,
        )

    def test_escalation_disabled_goes_to_failed(self):
        state = self._failed_state(enable_escalation=False)
        assert validation_router(state) == "package_result_failed"

    def test_escalation_enabled_goes_to_escalation(self):
        state = self._failed_state(enable_escalation=True, escalation_used=False)
        assert validation_router(state) == "escalation_agent"

    def test_escalation_used_blocks_reentry(self):
        state = self._failed_state(enable_escalation=True, escalation_used=True)
        assert validation_router(state) == "package_result_failed"


# ---------------------------------------------------------------------------
# Integration: subgraph with enable_escalation=False preserves current behavior
# ---------------------------------------------------------------------------

def test_escalation_disabled_preserves_behavior():
    graph = compile_chain_subgraph()
    state = _make_input_state(max_remediation_attempts=0, enable_escalation=False)
    result = graph.invoke(state)
    assert len(result["chain_results"]) == 1
    cr = result["chain_results"][0]
    assert cr["status"] in ("passed", "max_retries_exceeded", "failed")
    assert cr["status"] != "escalation_failed"
    assert cr["status"] != "host_incompatible"


# ---------------------------------------------------------------------------
# Integration: mocked LLM escalation paths
# ---------------------------------------------------------------------------

def _mock_llm_response(json_obj):
    """Return a mock ChatOpenRouter whose invoke() returns json_obj as content."""
    mock_llm = MagicMock()
    mock_msg = MagicMock()
    mock_msg.content = json.dumps(json_obj)
    mock_llm.invoke.return_value = mock_msg
    return mock_llm


def _failing_validation(state):
    """Sequence validation that always fails with a CAI check."""
    from aixbio.models.validation import ChainValidation
    return {
        "chain_validation": ChainValidation(
            id=state["chain"].id,
            passed=False,
            checks=(CheckResult(name="cai", passed=False, value=0.5, threshold=">0.8"),),
        ),
        "failed_checks": (CheckResult(name="cai", passed=False, value=0.5, threshold=">0.8"),),
    }


@patch("aixbio.nodes.escalation_agent.ChatOpenRouter")
def test_escalation_apply_plan_valid_swaps(mock_ctor):
    from aixbio.tools.codon_tables import best_ecoli_codon, split_codons, synonymous_alternatives

    dna = "".join(best_ecoli_codon(aa) for aa in INSULIN_B.aa_sequence)
    codons = split_codons(dna)
    alt = synonymous_alternatives(codons[0])

    mock_ctor.return_value = _mock_llm_response({
        "kind": "apply_plan",
        "actions": [{
            "check_name": "cai",
            "strategy": "swap",
            "target_positions": [0],
            "replacement_codons": [alt[0]],
        }],
        "reasoning": "break local min",
        "diagnosis": "greedy ordering stuck",
    })

    graph = compile_chain_subgraph()
    state = _make_input_state(
        max_remediation_attempts=0,
        enable_escalation=True,
    )
    with patch(
        "aixbio.graph.chain_subgraph.sequence_validation",
        _failing_validation,
    ):
        graph = compile_chain_subgraph()
        result = graph.invoke(state)
    assert result["escalation_used"] is True


@patch("aixbio.nodes.escalation_agent.ChatOpenRouter")
def test_escalation_incompatible(mock_ctor):
    mock_ctor.return_value = _mock_llm_response({
        "kind": "incompatible",
        "reason": "EPO requires N-glycosylation",
        "suggested_action": "Use Pichia pastoris",
    })

    state = _make_input_state(
        max_remediation_attempts=0,
        enable_escalation=True,
    )
    with patch(
        "aixbio.graph.chain_subgraph.sequence_validation",
        _failing_validation,
    ):
        graph = compile_chain_subgraph()
        result = graph.invoke(state)
    assert result["escalation_used"] is True
    cr = result["chain_results"][-1]
    assert cr["status"] == "host_incompatible"
    assert any("INCOMPATIBLE" in w for w in result["warnings"])


@patch("aixbio.nodes.escalation_agent.ChatOpenRouter")
def test_escalation_nonsynonymous_swaps_dropped(mock_ctor):
    mock_ctor.return_value = _mock_llm_response({
        "kind": "apply_plan",
        "actions": [{
            "check_name": "cai",
            "target_positions": [0],
            "replacement_codons": ["ATG"],
        }],
        "reasoning": "bad swap",
        "diagnosis": "test",
    })

    state = _make_input_state(
        max_remediation_attempts=0,
        enable_escalation=True,
    )
    with patch(
        "aixbio.graph.chain_subgraph.sequence_validation",
        _failing_validation,
    ):
        graph = compile_chain_subgraph()
        result = graph.invoke(state)
    assert result["escalation_used"] is True
    cr = result["chain_results"][-1]
    assert cr["status"] in ("escalation_failed", "host_incompatible", "failed")


@patch("aixbio.nodes.escalation_agent.ChatOpenRouter")
def test_escalation_give_up(mock_ctor):
    mock_ctor.return_value = _mock_llm_response({
        "kind": "give_up",
        "diagnosis": "pipeline bug",
    })

    state = _make_input_state(
        max_remediation_attempts=0,
        enable_escalation=True,
    )
    with patch(
        "aixbio.graph.chain_subgraph.sequence_validation",
        _failing_validation,
    ):
        graph = compile_chain_subgraph()
        result = graph.invoke(state)
    assert result["escalation_used"] is True


@patch("aixbio.nodes.escalation_agent.ChatOpenRouter")
def test_escalation_used_blocks_reentry(mock_ctor):
    """After escalation fires once (apply_plan -> apply_fixes -> revalidate still fails),
    the escalation_used flag prevents a second LLM call."""
    from aixbio.tools.codon_tables import best_ecoli_codon, split_codons, synonymous_alternatives

    dna = "".join(best_ecoli_codon(aa) for aa in INSULIN_B.aa_sequence)
    codons = split_codons(dna)
    alt = synonymous_alternatives(codons[0])

    mock_ctor.return_value = _mock_llm_response({
        "kind": "apply_plan",
        "actions": [{
            "check_name": "cai",
            "strategy": "swap",
            "target_positions": [0],
            "replacement_codons": [alt[0]],
        }],
        "reasoning": "attempt",
        "diagnosis": "stuck",
    })

    state = _make_input_state(
        max_remediation_attempts=0,
        enable_escalation=True,
    )
    with patch(
        "aixbio.graph.chain_subgraph.sequence_validation",
        _failing_validation,
    ):
        graph = compile_chain_subgraph()
        result = graph.invoke(state)

    assert result["escalation_used"] is True
    assert mock_ctor.return_value.invoke.call_count == 1
