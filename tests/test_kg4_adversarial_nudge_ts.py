"""Red control for the nudge-without-timestamp crash (nacht-s19 finding,
adversarial pass on 1ebf8a8e)."""
import json

from delfin.agent import longrun_metric as lm


def test_malformed_nudge_ts_does_not_crash_and_is_not_counted():
    """A nudge record with a missing (or non-numeric) ts must neither crash
    the reducer nor be counted as a nudge-while-waiting — the module's own
    contract: what the records cannot attest is None/not counted."""
    inp = lm.LongrunInputs(
        turns=[{"ts": 1000.0}],
        next_turn_ts=2000.0,
        nudges=[
            {"ts": 1500.0, "kind": "auto-continue"},   # inside window: counts
            {"kind": "auto-continue"},                 # ts missing: must not crash
            {"ts": "not-a-number", "kind": "other"},   # ts garbage: must not crash
            {"ts": None},                              # ts None: must not crash
        ],
    )
    res = lm.summarize(inp)
    assert res.nudges_while_waiting == 1  # only the well-formed one


def test_malformed_nudge_ts_only_records_still_reduce_to_none_or_zero():
    """Only malformed nudges: no crash, and nothing attested."""
    inp = lm.LongrunInputs(
        turns=[{"ts": 1000.0}],
        next_turn_ts=2000.0,
        nudges=[{"kind": "auto-continue"}, {"ts": "garbage"}],
    )
    res = lm.summarize(inp)
    assert res.nudges_while_waiting in (0, None)
