"""A damaged JSON state file is moved aside, and its loss is discoverable.

Measured on 2026-11-17 in the agent's long-term memory: two remembered
facts, a file truncated the way a full disk or a kill mid-write truncates
it, ONE ``save_memory`` call -- and the file held only the new fact. The
two were gone with no copy anywhere.

The cause is one read answering two questions. ``except
(FileNotFoundError, json.JSONDecodeError): return {"facts": []}`` is a
true answer for a READER (nothing is remembered) and destroys data for a
WRITER, which builds the new file from it. ``agent/state_file.py`` now
separates the two, and the writers quarantine instead of overwriting.

Three further properties are pinned here because each was a real defect
during the work, not a hypothetical:

* the write is atomic -- the plain ``write_text`` it replaced truncated
  first, which is how the damaged file in the first place comes about;
* a second corruption does not overwrite the first rescue (``os.replace``
  overwrites its destination silently, so the name is reserved first);
* the empty store is built fresh per call. As a module-level constant,
  ``dict(CONSTANT)`` copied only the outer mapping, the ``facts`` list
  stayed shared, and a caller appending to it mutated the constant -- so
  every later reader of an "empty" store saw that fact.

And because the one production caller of ``save_memory`` swallows every
exception (``memory_distill``: ``except Exception: continue``), the
notification alone would leave the loss invisible. The doctor therefore
reports both a state file that does not parse and a file kept aside, so
it reaches a person through the report they already read.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import memory_store as M
from delfin.agent import state_file as SF

TRUNCATED = '{"facts": [{"text": "Fact A, kept"}, {"text": "Fact B'


@pytest.fixture
def store(tmp_path):
    p = tmp_path / "agent_memory.json"
    M.save_memory("Fact A, kept", path=p)
    M.save_memory("Fact B, kept too", path=p)
    assert len(M.load_memories(path=p)) == 2
    return p


def _aside(path):
    return sorted(q for q in path.parent.iterdir()
                  if "unreadable" in q.name)


# --- the shared helper ----------------------------------------------------

def test_an_absent_file_is_not_an_error_in_either_mode(tmp_path):
    """Nothing was stored, so there is nothing to lose by writing."""
    missing = tmp_path / "nothing.json"
    assert SF.read_json(missing, default={"a": 1}) == {"a": 1}
    assert SF.read_json(missing, default=None, strict=True) is None


def test_an_unreadable_file_separates_the_two_modes(tmp_path):
    p = tmp_path / "broken.json"
    p.write_text(TRUNCATED)
    assert SF.read_json(p, default=None) is None
    with pytest.raises(SF.StateFileUnreadable):
        SF.read_json(p, default=None, strict=True)


def test_the_quarantine_keeps_the_bytes_and_returns_where(tmp_path):
    p = tmp_path / "broken.json"
    p.write_text(TRUNCATED)
    aside = SF.quarantine(p)
    assert aside is not None and aside.exists()
    assert not p.exists()
    assert aside.read_text() == TRUNCATED, "a rename, not a rewrite"
    assert SF.quarantine(p) is None, "nothing to move is not an error"


def test_two_corruptions_in_one_second_keep_two_files(tmp_path):
    """os.replace overwrites its destination without a word."""
    p = tmp_path / "broken.json"
    p.write_text("first damage")
    SF.quarantine(p)
    p.write_text("second damage")
    SF.quarantine(p)
    kept = sorted(q.read_text() for q in tmp_path.iterdir()
                  if "unreadable" in q.name)
    assert kept == ["first damage", "second damage"]


def test_the_atomic_write_leaves_no_half_file_and_keeps_the_mode(tmp_path):
    import stat

    p = tmp_path / "s.json"
    SF.write_json_atomic(p, {"a": [1, 2]}, mode=0o600)
    assert json.loads(p.read_text()) == {"a": [1, 2]}
    assert stat.S_IMODE(p.stat().st_mode) == 0o600
    assert not (tmp_path / "s.json.tmp").exists(), "no leftover"


# --- the memory store -----------------------------------------------------

def test_a_reader_still_answers_from_a_damaged_store(store):
    """Lenient on purpose: a turn must not raise over a broken file."""
    store.write_text(TRUNCATED)
    assert M.load_memories(path=store) == []


def test_saving_a_fact_keeps_the_damaged_store_and_reports(store):
    store.write_text(TRUNCATED)
    with pytest.raises(M.MemoryStoreMovedAside) as exc:
        M.save_memory("Fact C, new", path=store)

    # The new fact IS stored: the report is not a failed write.
    assert [m["text"] for m in M.load_memories(path=store)] == ["Fact C, new"]
    kept = _aside(store)
    assert len(kept) == 1
    assert "Fact A, kept" in kept[0].read_text(), (
        "what the damaged file held must still be there")
    assert "Fact A" not in str(exc.value), "the message names the file only"


def test_deleting_from_a_damaged_store_refuses_and_changes_nothing(store):
    store.write_text(TRUNCATED)
    before = store.read_text()
    with pytest.raises(SF.StateFileUnreadable):
        M.delete_memory(0, path=store)
    assert store.read_text() == before
    assert not _aside(store), "a delete does not decide the file's fate"


def test_the_empty_store_is_fresh_each_time(tmp_path):
    """A shared `facts` list made every later "empty" store carry a fact."""
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"
    a.write_text(TRUNCATED)
    b.write_text(TRUNCATED)
    first = M._read(a)
    first["facts"].append({"text": "appended by a caller"})
    assert M._read(b)["facts"] == [], (
        "one caller's append must not reach another read")


def test_an_intact_store_is_untouched_by_all_of_this(store):
    """The control."""
    M.save_memory("Fact C, new", path=store)
    assert [m["text"] for m in M.load_memories(path=store)] == [
        "Fact A, kept", "Fact B, kept too", "Fact C, new"]
    assert not _aside(store)
    assert M.delete_memory(0, path=store) is True
    assert [m["text"] for m in M.load_memories(path=store)] == [
        "Fact B, kept too", "Fact C, new"]


# --- the doctor makes it discoverable -------------------------------------

def test_the_doctor_names_a_state_file_that_does_not_parse(
        tmp_path, monkeypatch):
    from delfin.agent import doctor as D

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    (tmp_path / ".delfin_settings.json").write_text('{"agent": {"model')
    root = tmp_path / ".delfin"
    root.mkdir()
    monkeypatch.setattr("delfin.agent.state_paths.state_root",
                        lambda *a, **k: str(root))

    rows = D._check_state_files({})
    bad = [r for r in rows if r["status"] == D.FAIL]
    assert bad, "a file that does not parse must be reported"
    assert ".delfin_settings.json" in bad[0]["detail"]
    assert "settings" in bad[0]["detail"], "name what the file holds"
    assert not bad[0].get("command"), (
        "deliberately prose-only: starting fresh discards the person's own "
        "settings, which no approval-only step may do")


def test_the_doctor_reports_a_file_kept_aside(tmp_path, monkeypatch):
    from delfin.agent import doctor as D

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    root = tmp_path / ".delfin"
    root.mkdir()
    (root / "agent_memory.json.unreadable-20261117-120000").write_text("x")
    monkeypatch.setattr("delfin.agent.state_paths.state_root",
                        lambda *a, **k: str(root))

    rows = D._check_state_files({})
    warned = [r for r in rows if r["status"] == D.WARN]
    assert warned, "a kept file must be reported: its notification is "\
                   "swallowed by the one production caller"
    assert "unreadable" in warned[0]["detail"]


def test_the_doctor_passes_when_everything_parses(tmp_path, monkeypatch):
    from delfin.agent import doctor as D

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    root = tmp_path / ".delfin"
    root.mkdir()
    (root / "credentials.json").write_text('{"A": "stand-in"}')
    monkeypatch.setattr("delfin.agent.state_paths.state_root",
                        lambda *a, **k: str(root))

    rows = D._check_state_files({})
    assert [r["status"] for r in rows] == [D.PASS]
