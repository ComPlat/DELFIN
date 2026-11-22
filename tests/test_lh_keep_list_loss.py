"""F2: the long-horizon metric carries the keep-list loss.

The compaction record on disk carries the quality report per keep-list
category (compaction_log.record_compaction writes it as ``lost``);
delfin/agent/longhorizon/metrics.py receives those records and must carry
the loss per category, in compaction_quality.category_names() order —
otherwise the yardstick cannot see the loss the product measures.
"""
from __future__ import annotations

from delfin.agent.longhorizon import metrics as m
from delfin.agent import compaction_quality as cq


def test_keep_list_loss_carried_per_category_in_category_order():
    """One compaction record attesting losses in two categories: the
    metric reports every category_names() category, in that exact
    order; the attested categories carry their counts and bounded
    names; an unattested category reports None, never a guessed zero."""
    comps = [{
        "ts": 100.0, "kind": "auto", "messages_compacted": 12,
        "keep_list_loss": {
            "goal": {"total": 1, "lost": ["refine the beta_HRS workflow"]},
            "refusals": {"total": 2, "lost": ["Refused: touch archive",
                                              "Refused: rm results"]},
            "tests": {"total": 0, "lost": []},
        },
    }]
    r = m.measure([], [], comps, planted_facts=[])
    cats = cq.category_names()
    assert [e["category"] for e in r.keep_list_loss] == list(cats), \
        "keep_list_loss must follow category_names() order exactly"
    by_cat = {e["category"]: e for e in r.keep_list_loss}
    # attested: counts and names carried
    assert by_cat["goal"]["lost_total"] == 1
    assert by_cat["goal"]["lost_names"] == ["refine the beta_HRS workflow"]
    assert by_cat["refusals"]["lost_total"] == 2
    assert len(by_cat["refusals"]["lost_names"]) == 2
    # attested zero stays a zero (the record attests 0 lost)
    assert by_cat["tests"]["lost_total"] == 0
    # not attested: None, never 0 (the module invariant)
    for cat in ("instruction", "constraints", "decisions", "waits"):
        assert by_cat[cat]["lost_total"] is None, \
            f"{cat} unattested must be None, not 0"
        assert by_cat[cat]["lost_names"] == []


def test_keep_list_loss_worst_compaction_wins_and_names_stay_bounded():
    """Several compactions: the count is the WORST (max) across them —
    a sum would double-count an item every compaction dropped; the
    names stay bounded (never an unbounded dump from a long session)."""
    comps = [
        {"ts": 10.0, "kind": "auto",
         "keep_list_loss": {"goal": {"total": 1, "lost": ["lost-a"]}}},
        {"ts": 20.0, "kind": "auto",
         "keep_list_loss": {"goal": {"total": 3,
                           "lost": [f"lost-{i}" for i in range(30)]}}},
        {"ts": 30.0, "kind": "auto",
         "keep_list_loss": {"goal": {"total": 2, "lost": ["lost-x"]}}},
    ]
    r = m.measure([], [], comps, planted_facts=[])
    by_cat = {e["category"]: e for e in r.keep_list_loss}
    assert by_cat["goal"]["lost_total"] == 3, \
        "max across compactions, not a sum: got " \
        f"{by_cat['goal']['lost_total']}"
    assert len(by_cat["goal"]["lost_names"]) <= 10, \
        "names must stay bounded"


def test_no_compaction_records_no_loss_is_reported_per_category():
    """No compaction records: every category reports None (unmeasurable),
    not 0 — and no crash."""
    r = m.measure([], [], [], planted_facts=[])
    assert len(r.keep_list_loss) == len(cq.category_names())
    assert all(e["lost_total"] is None for e in r.keep_list_loss)
    assert all(e["lost_names"] == [] for e in r.keep_list_loss)


# --- controls: unchanged behaviour, expected to pass on BOTH sides ------


def test_control_compactions_size_unchanged():
    """Control: the pre-existing compactions_size reduction is untouched."""
    comps = [{"ts": 100.0, "kind": "auto", "messages_compacted": 7,
              "tokens_before": 900, "tokens_after": 300}]
    r = m.measure([], [], comps, planted_facts=[])
    assert r.compactions_size == [{"kind": "auto",
                                   "messages_compacted": 7,
                                   "tokens_before": 900,
                                   "tokens_after": 300}]
    assert r.compactions == 1


def test_control_record_with_no_lost_key_is_tolerated():
    """Control: an old compaction record without ``lost`` changes nothing
    in the pre-existing reduction — no crash, no new size entries."""
    comps = [{"ts": 5.0, "kind": "auto"}]
    r = m.measure([], [], comps, planted_facts=[])
    assert r.compactions == 1
    assert r.compactions_size == [{"kind": "auto", "messages_compacted": None,
                                   "tokens_before": None, "tokens_after": None}]
