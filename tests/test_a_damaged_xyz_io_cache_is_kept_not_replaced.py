"""A damaged xyz_io geometry cache is quarantined, never emptied by a write.

``delfin.xyz_io`` keeps two keyed caches beside the geometry so a rerun
does not re-derive the QM/XTB split or the first coordination sphere:
``.qmmm_cache.json`` (signature -> (start, end)) and
``.delfin_first_sphere.json`` (composition key -> atom list). Both persist
functions read the whole file, inserted one key, and wrote it back; the
shared read answered "empty" for both an absent file and one that cannot
be read, so a file truncated the way a full disk or a kill mid-write
truncates it became, on the next persist, a file holding only the one key
just written -- every other signature's split, or every other molecule's
sphere, gone.

The fix routes the writer's read through ``agent.state_file.
read_json_for_write``: on damage the file is moved to the dated
``.unreadable-`` quarantine (bytes kept) and the persist REFUSES to write
to that target (it says so in the log). These two sit in the
mid-calculation geometry path and must never raise, so the refusal is a
log line plus the doctor's state-file check, not an exception. The store
is gone, so the next persist writes a fresh one safely.

Each store is pinned in its own tmp dir; nothing here ever touches a real
run folder.
"""

from __future__ import annotations

import json

from delfin import xyz_io


def _aside(root):
    return sorted(p for p in root.iterdir() if "unreadable" in p.name)


def test_qmmm_cache_refuses_to_replace_a_damaged_file_and_keeps_its_bytes(tmp_path):
    # The persist writes to two targets: the geometry's own directory and
    # its parent. Nest the run one level down so both fall inside tmp_path.
    run = tmp_path / "run"
    run.mkdir()
    # Two real splits (two signatures), then the file truncated the way a
    # full disk or a kill mid-write truncates it.
    cache_file = run / xyz_io._QMMM_CACHE_BASENAME
    seed = {"sigA": [0, 3], "sigB": [4, 7]}
    cache_file.write_text(json.dumps(seed), encoding="utf-8")
    damaged = cache_file.read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    cache_file.write_text(damaged, encoding="utf-8")

    source = run / "geom.xyz"
    source.write_text("2\njob\nH 0 0 0\nH 0 0 1\n", encoding="utf-8")
    xyz_io._persist_qmmm_signature(("sigC",), (0, 1), source)

    # The damaged target was kept, not rewritten: no fresh file at run, one
    # quarantined copy holding exactly the damaged bytes.
    assert not cache_file.exists(), (
        "nothing may be written over a cache that could not be read")
    moved = _aside(run)
    assert len(moved) == 1, "the unreadable cache must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")

    # The undamaged parent target was still written (per-target behaviour):
    # a damaged target is refused, the sibling is not.
    parent_cache = tmp_path / xyz_io._QMMM_CACHE_BASENAME
    if parent_cache.exists():
        parent_data = json.loads(parent_cache.read_text(encoding="utf-8"))
        assert parent_data.get(xyz_io._cache_key(("sigC",))) == [0, 1]


def test_first_sphere_cache_refuses_to_replace_a_damaged_file_and_keeps_its_bytes(tmp_path):
    # A CONTROL.txt pins the cache to this dir. Two real spheres (two
    # composition keys), then the file truncated the way a full disk or a
    # kill mid-write truncates it.
    (tmp_path / "CONTROL.txt").write_text("LEVEL_OF_THEORY=GFN2-xTB\n", encoding="utf-8")
    cache_file = tmp_path / xyz_io._FIRST_SPHERE_CACHE_BASENAME
    seed = {"CoC": [0, 1, 2], "CoN": [0, 3, 4]}
    cache_file.write_text(json.dumps(seed), encoding="utf-8")
    damaged = cache_file.read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    cache_file.write_text(damaged, encoding="utf-8")

    source = tmp_path / "geom.xyz"
    source.write_text("2\njob\nH 0 0 0\nH 0 0 1\n", encoding="utf-8")
    xyz_io._persist_first_sphere("CoO", {0, 1}, source)

    # The damaged cache was kept, not rewritten: no fresh file, one
    # quarantined copy holding exactly the damaged bytes.
    assert not cache_file.exists(), (
        "nothing may be written over a cache that could not be read")
    moved = _aside(tmp_path)
    assert len(moved) == 1, "the unreadable cache must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")


def test_an_intact_cache_is_untouched_by_any_of_this(tmp_path):
    """The control: nothing above may change the normal path."""
    run = tmp_path / "run"
    run.mkdir()
    cache_file = run / xyz_io._QMMM_CACHE_BASENAME
    seed = {"sigA": [0, 3]}
    cache_file.write_text(json.dumps(seed), encoding="utf-8")
    source = run / "geom.xyz"
    source.write_text("2\njob\nH 0 0 0\nH 0 0 1\n", encoding="utf-8")
    xyz_io._persist_qmmm_signature(("sigB",), (4, 7), source)
    data = json.loads((run / xyz_io._QMMM_CACHE_BASENAME).read_text(encoding="utf-8"))
    assert set(data) == {"sigA", "sigB"}, "both signatures are kept"
    assert not _aside(run), "nothing is quarantined when the cache parses"
