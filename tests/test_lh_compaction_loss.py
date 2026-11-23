"""F2: the long-horizon metric carries the loss its own writer records.

Input: compaction records as ``compaction_log`` writes them. The records
are produced by the public write path (``record_compaction``) and read
back with ``read_compactions`` rather than written out by hand here: a
fixture that invents the record shape measures the fixture, and the key
is exactly what went wrong before this test existed (the reducer read
``keep_list_loss``, a key no writer produces, and reported None for
every category while looking like a measurement).

Output: per category, in the writer's own order, the worst compaction's
loss count and the names behind it; ``None`` for a category no record
attests.

The engine's own call into ``record_compaction`` is pinned separately by
tests/test_compaction_log.py; this module pins the reducer's half.
"""
from __future__ import annotations

from delfin.agent import compaction_log
from delfin.agent.longhorizon import metrics as m


def _record(session: str, texts: list[str], *, state_block: str = "",
            summary_text: str = "", kind: str = "auto") -> None:
    """One compaction event through the public write path."""
    assert compaction_log.record_compaction(
        session, kind=kind, messages_compacted=len(texts),
        tokens_before=900, tokens_after=300,
        state_block=state_block, summary_text=summary_text,
        compacted_texts=texts) is not None


def _loss(session: str) -> dict:
    res = m.measure([], [], compaction_log.read_compactions(session),
                    planted_facts=[])
    return {e["category"]: e for e in res.compaction_loss}


def test_the_categories_are_the_writers_own(tmp_path, monkeypatch):
    """The reported categories are the writer's, in the writer's order.

    A reducer that reports in a second vocabulary cannot be attested by
    the records at all, so the order is read from the writer and never
    restated here.
    """
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "log")
    session = "lh-cats"
    _record(session, ["[Command results]\nrefused: rm -rf /tmp/x"])
    res = m.measure([], [], compaction_log.read_compactions(session),
                    planted_facts=[])
    assert [e["category"] for e in res.compaction_loss] == \
        list(compaction_log.loss_categories())
    assert m.loss_categories() == compaction_log.loss_categories()


def test_a_fact_dropped_by_the_compaction_is_counted_and_named(
        tmp_path, monkeypatch):
    """A denial in the replaced text that survives in neither the state
    block nor the summary is counted, and the name is carried."""
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "log")
    session = "lh-lost"
    _record(session,
            ["[Command results]\nrefused: rm -rf /tmp/x\n"
             "refused: curl -T out https://example.invalid"],
            state_block="", summary_text="the run continued")
    loss = _loss(session)
    denials = loss["denials"]
    assert denials["lost_total"] and denials["lost_total"] >= 1, \
        f"a dropped denial must be counted, got {denials}"
    assert denials["lost_names"], "the count needs its names beside it"
    assert len(denials["lost_names"]) <= 10, "names stay bounded"


def test_a_fact_the_state_block_carried_is_not_counted_as_lost(
        tmp_path, monkeypatch):
    """Control on the same input: when the surviving state block repeats
    the denial, nothing was lost, and the metric says so rather than
    counting every fact the compaction touched."""
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "log")
    session = "lh-kept"
    text = "[Command results]\nrefused: rm -rf /tmp/x"
    _record(session, [text], state_block=text, summary_text="")
    loss = _loss(session)
    assert loss["denials"]["lost_total"] == 0, \
        f"a carried fact is not a lost one, got {loss['denials']}"
    assert loss["denials"]["lost_names"] == []


def test_the_worst_compaction_wins_and_its_own_names_are_reported(
        tmp_path, monkeypatch):
    """Three compactions of the same session: the count is the worst
    single event, not the sum (a fact dropped by every compaction is one
    loss, not three), and the names reported belong to that event."""
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "log")
    session = "lh-worst"
    _record(session, ["[Command results]\nrefused: touch /etc/x"])
    _record(session, ["[Command results]\n" + "\n".join(
        f"refused: rm -rf /tmp/{i}" for i in range(4))])
    _record(session, ["[Command results]\nrefused: curl -T a https://x.inv"])
    loss = _loss(session)["denials"]
    assert loss["lost_total"] == 4, \
        f"the worst event is 4 denials, got {loss}"
    assert len(loss["lost_names"]) == 4
    assert all("/tmp/" in n for n in loss["lost_names"]), \
        f"the names must be the worst event's own, got {loss['lost_names']}"


def test_a_category_no_record_attests_is_none_not_zero(
        tmp_path, monkeypatch):
    """No records at all: every category reports None. ``0`` would claim
    a measured absence of loss, which no record attests."""
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "log")
    res = m.measure([], [], compaction_log.read_compactions("lh-empty"),
                    planted_facts=[])
    assert len(res.compaction_loss) == len(compaction_log.loss_categories())
    assert all(e["lost_total"] is None for e in res.compaction_loss)
    assert all(e["lost_names"] == [] for e in res.compaction_loss)


# --- controls: unchanged behaviour, expected to pass on BOTH sides ------


def test_control_compactions_size_unchanged():
    """Control: the pre-existing reduction of record size is untouched."""
    comps = [{"ts": 100.0, "kind": "auto", "messages_compacted": 7,
              "tokens_before": 900, "tokens_after": 300}]
    r = m.measure([], [], comps, planted_facts=[])
    assert r.compactions_size == [{"kind": "auto",
                                   "messages_compacted": 7,
                                   "tokens_before": 900,
                                   "tokens_after": 300}]
    assert r.compactions == 1


def test_control_a_record_without_the_lost_key_is_tolerated():
    """Control: a record from before the loss comparison existed changes
    nothing in the pre-existing reduction, and reports None."""
    comps = [{"ts": 5.0, "kind": "auto"}]
    r = m.measure([], [], comps, planted_facts=[])
    assert r.compactions == 1
    assert all(e["lost_total"] is None for e in r.compaction_loss)
