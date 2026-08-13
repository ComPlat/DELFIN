"""Reading a spreadsheet used to vouch for every energy in the answer.

The physical-quantity guard asked a turn-level question: did ANY lookup
tool run, or was ANY calculation-output-like file observed. If so,
nothing was flagged. Measured on the shipped scanner:

    3 flags | Die Messung ergab 2.31 eV für S1, 3.90 eV für S2, 412.7 kcal/mol
    0 flags | ... the same sentence, with evidence_tools_used={"read_document"}

One read of one unrelated cost-centre sheet exempted every number. That
is the office ledger's defect on the other side of the house, and it has
the same repair, which is why it was deliberately not fixed as a second
mechanism: a claim is grounded when the FIGURE appears in what a tool
returned, not when some tool ran at all.

The two ledgers stay apart because the evidence has different shapes. The
office one takes structured figures out of result dicts and never parses
prose, so two readings of the same facts cannot disagree. Here there is
no dict: the evidence for "2.31 eV" is the text of an ORCA output the
agent read, and that text is the primary source rather than a rendering
of one.
"""

from __future__ import annotations

import pytest

from delfin.agent import verify_guard as vg


CLAIM = ("Die Messung ergab 2.31 eV für S1, 3.90 eV für S2 und "
         "412.7 kcal/mol Bindungsenergie.")

TDDFT_OUTPUT = """
STATE  1:  E=   0.084912 au      2.310 eV     18631.3 cm**-1
STATE  2:  E=   0.143331 au      3.900 eV     31456.9 cm**-1
Final binding energy   412.7 kcal/mol
"""

SPREADSHEET_OUTPUT = "Kostenstelle 4711  1.234,50 EUR  Betrag 289.90"


@pytest.fixture(autouse=True)
def _fresh_turn():
    vg.reset_observed_numbers()
    yield
    vg.reset_observed_numbers()


# ---------------------------------------------------------------------------
# The defect
# ---------------------------------------------------------------------------

def test_reading_an_unrelated_sheet_does_not_vouch_for_the_energies():
    """The whole point. Four numbers came back and none of them is an
    energy in this answer."""
    vg.record_tool_numbers(SPREADSHEET_OUTPUT)
    flags = vg.scan_for_unsourced_quantities(
        CLAIM,
        observed_files={"/w/calc/tddft.out"},        # would have exempted
        evidence_tools_used={"read_document"},       # would have exempted
        numbers=vg.observed_numbers())
    assert len(flags) == 3
    assert {f.unit for f in flags} == {"eV", "kcal/mol"}


def test_reading_the_output_that_holds_them_does():
    vg.record_tool_numbers(TDDFT_OUTPUT)
    assert vg.scan_for_unsourced_quantities(
        CLAIM, numbers=vg.observed_numbers()) == []


def test_a_number_the_output_does_not_hold_is_still_caught():
    vg.record_tool_numbers(TDDFT_OUTPUT)
    flags = vg.scan_for_unsourced_quantities(
        "Die Anregung liegt bei 7.77 eV.", numbers=vg.observed_numbers())
    assert [f.quantity for f in flags] == ["7.77 eV"]


def test_a_gap_derived_from_two_observed_energies_is_not_an_invention():
    """3.90 - 2.31 = 1.59. An answer that read both and reports the gap
    has done arithmetic; flagging it would teach the model to stop
    reporting gaps, and a guard that punishes honest work gets worked
    around rather than obeyed."""
    vg.record_tool_numbers(TDDFT_OUTPUT)
    assert vg.scan_for_unsourced_quantities(
        "Die Lücke zwischen S1 und S2 beträgt 1.59 eV.",
        numbers=vg.observed_numbers()) == []


# ---------------------------------------------------------------------------
# None is not empty
# ---------------------------------------------------------------------------

def test_a_turn_this_pool_cannot_see_keeps_the_old_shape():
    """None means no result ever reached this pool -- a headless caller, a
    path that does not record. There is nothing to check a value against,
    and flagging an answer whose evidence the function is blind to is a
    guard punishing work it cannot see."""
    assert vg.observed_numbers() is None
    assert vg.scan_for_unsourced_quantities(
        CLAIM, observed_files={"/w/calc/tddft.out"}) == []


def test_results_that_carried_no_number_are_not_the_same_as_no_results():
    """An empty pool says tools ran and returned nothing numeric. A
    quantity in the answer is then unsourced, and saying so is right."""
    vg.record_tool_numbers("no numerals here at all")
    assert vg.observed_numbers() == []
    assert len(vg.scan_for_unsourced_quantities(
        CLAIM, observed_files={"/w/calc/tddft.out"},
        numbers=vg.observed_numbers())) == 3


def test_the_pool_does_not_survive_the_turn():
    """Carried over, a stale energy would ground a later unrelated one --
    the failure this exists to catch, one turn late."""
    vg.record_tool_numbers(TDDFT_OUTPUT)
    assert vg.scan_for_unsourced_quantities(
        CLAIM, numbers=vg.observed_numbers()) == []
    vg.reset_observed_numbers()
    assert vg.observed_numbers() is None


# ---------------------------------------------------------------------------
# Bounded, because an .out file is megabytes
# ---------------------------------------------------------------------------

def test_a_huge_result_neither_hangs_nor_fills_the_pool():
    """This module has already had one incident where a single long line
    made the grounding scan take minutes."""
    import time
    huge = " ".join(str(i) for i in range(400_000))
    t0 = time.monotonic()
    vg.record_tool_numbers(huge)
    elapsed = time.monotonic() - t0
    assert len(vg.observed_numbers()) <= vg.MAX_OBSERVED_NUMBERS
    assert elapsed < 5.0, elapsed


def test_a_result_that_arrived_cut_short_makes_the_pool_say_nothing():
    """A tool result reaches the recorder as a HEAD slice. One truncated
    one leaves a hole of unknown size, and absence from a pool with a
    hole is not evidence -- it would flag a chemist's correctly-quoted
    energy from deep inside a long .out file."""
    vg.record_tool_numbers("STATE 1: 2.310 eV", truncated=True)
    assert vg.observations_are_complete() is False
    assert vg.observed_numbers() is None
    assert vg.scan_for_unsourced_quantities(
        CLAIM, observed_files={"/w/calc/tddft.out"},
        numbers=vg.observed_numbers()) == []


def test_one_cut_result_taints_the_whole_turn():
    """A complete read afterwards does not repair the hole the first one
    left; the numbers past that cut are still unseen."""
    vg.record_tool_numbers(SPREADSHEET_OUTPUT)
    assert vg.observed_numbers() is not None
    vg.record_tool_numbers("more output", truncated=True)
    assert vg.observed_numbers() is None
    vg.record_tool_numbers(TDDFT_OUTPUT)
    assert vg.observed_numbers() is None


def test_recording_never_raises():
    for junk in (None, 12, object(), b"\x00\x01"):
        vg.record_tool_numbers(junk)


# ---------------------------------------------------------------------------
# ... and it is wired, not merely written
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    """The smallest prompt pack the engine will start on."""
    import textwrap
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def test_a_tool_result_reaches_the_pool_through_a_real_turn(agent_tree):
    """A mechanism nothing calls is a comment, and today alone this
    repository found two guards that existed and were skipped by one hop
    in the middle. So this drives an actual turn and looks at the pool
    afterwards, rather than reading the engine's source for the call --
    a source check passes on a comment and breaks on a rename."""
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    fake = MagicMock()
    fake._observed_files_session = set()

    def _stream(*a, **k):
        yield StreamEvent(type="tool_use", tool_name="read_file",
                          tool_input='{"path": "calc/tddft.out"}')
        yield StreamEvent(type="tool_result", tool_name="read_file",
                          tool_output=TDDFT_OUTPUT)
        yield StreamEvent(type="text_delta",
                          text="Die S1-Energie ist 2.31 eV.")
        yield StreamEvent(type="message_delta", output_tokens=5,
                          cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)
    engine.stream_response("was ist die S1-Energie?")

    pool = vg.observed_numbers()
    assert pool is not None, "the tool result never reached the pool"
    assert any(abs(v - 2.310) < 1e-6 for v in pool), pool[:20]
