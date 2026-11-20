"""A damaged remote-watch store is quarantined, never emptied by a tick.

``delfin.remote.watch`` keeps what the cluster watcher found of every job
DELFIN submitted (``outcome``, ``seen_state``, ``followed`` read it) in
``~/.delfin/remote_watch.json``. ``tick`` is the one writer: it reads the
store, follows the open jobs one step, and writes the whole file back;
the shared read answered "nothing is known" for both an absent file and
one that cannot be read, so a store truncated the way a full disk or a
kill mid-write truncates it became, on the next tick, a store with no job
outcome in it -- every followed job re-probed from scratch, and an ended
job's results re-fetched or lost depending on what the cluster still
kept.

The fix splits the two reads. Readers keep the lenient answer -- a
damaged store must not crash whatever asked. The tick's write goes
through ``_load_for_write``, which on damage moves the file to the dated
``.unreadable-`` quarantine of ``agent.state_file`` (bytes kept): the
tick REFUSES to write (it logs where the store was kept and stores
nothing this round), and the doctor's state-file check tells the person.
The store is gone, so the next tick writes a fresh one safely.

The store is the suite's autouse redirect of ``watch._state_path`` into
the test's own tmp dir; nothing here ever touches the real home.

``cluster_tunnel`` is a host dependency the ``delfin.remote`` import
chain pulls in (``clusters.py`` references it only inside function
bodies, never at import time). This test measures the store's read/write
behaviour of ``watch.tick`` -- cluster probing is behind ``_follow`` /
``_fetch`` and is never reached with no open jobs -- so ``cluster_tunnel``
is irrelevant to what is measured. Where the host lacks it, a minimal
import-only stub shadows the missing module (the real one is used
everywhere it is installed), rather than a ``skipif`` that would leave
the store untested on those machines.
"""

from __future__ import annotations

import json
import sys
import types


def _ensure_cluster_tunnel_importable() -> None:
    """Make ``delfin.remote.watch`` importable even where the host lacks
    ``cluster_tunnel``; the real package is used wherever it is present."""
    try:
        import cluster_tunnel  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    pkg = types.ModuleType("cluster_tunnel")
    pkg.__path__ = []                       # mark it a package
    api = types.ModuleType("cluster_tunnel.api")
    config = types.ModuleType("cluster_tunnel.config")
    pkg.api, pkg.config = api, config
    sys.modules.setdefault("cluster_tunnel", pkg)
    sys.modules.setdefault("cluster_tunnel.api", api)
    sys.modules.setdefault("cluster_tunnel.config", config)


_ensure_cluster_tunnel_importable()

from delfin.remote import watch  # noqa: E402


def _state():
    return watch._state_path()


def _aside():
    return sorted(p for p in _state().parent.iterdir()
                  if "unreadable" in p.name)


def test_a_tick_refuses_to_replace_a_damaged_store_and_keeps_its_bytes():
    # A real store with one judged job outcome, then the file truncated
    # the way a full disk or a kill mid-write truncates it.
    seeded = {"jobs": {"cluster:42": {"state": "COMPLETED", "outcome": "ok"}},
              "probed": {"cluster": 1.0}}
    _state().parent.mkdir(parents=True, exist_ok=True)
    _state().write_text(json.dumps(seeded), encoding="utf-8")
    damaged = _state().read_text(encoding="utf-8")[:30]
    assert not damaged.endswith("}")   # genuinely torn, not valid JSON
    _state().write_text(damaged, encoding="utf-8")

    # Readers stay lenient while the damaged file is still present: a
    # broken store answers as if the job were unknown, it does not crash.
    assert watch.seen_state("cluster", "42") == ""

    out = watch.tick(now=1000.0, probe_every_s=0.0)

    # The tick refused to write: the damaged bytes are not rewritten,
    # truncated further, or silently replaced by a store that drops the
    # job outcomes.
    moved = _aside()
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert moved[0].read_text(encoding="utf-8") == damaged, (
        "the quarantined file must still hold exactly what was in it")
    assert not _state().exists(), (
        "nothing may be written over a store that could not be read")
    assert out.ended == [] and out.following == 0, (
        "a refused tick reports nothing as ended or followed")

    # The store is gone, so the next tick writes a fresh one safely.
    watch.tick(now=1001.0, probe_every_s=0.0)
    assert _state().exists(), "the next tick writes a fresh store"
    assert len(_aside()) == 1, "the second tick quarantines nothing"
    fresh = json.loads(_state().read_text(encoding="utf-8"))
    assert set(fresh) == {"jobs", "probed"}


def test_a_tick_on_an_intact_store_is_untouched_by_any_of_this():
    """The control: nothing above may change the normal path."""
    seeded = {"jobs": {"cluster:42": {"state": "RUNNING"}}, "probed": {}}
    _state().parent.mkdir(parents=True, exist_ok=True)
    _state().write_text(json.dumps(seeded), encoding="utf-8")

    watch.tick(now=1000.0, probe_every_s=0.0)

    assert _state().exists()
    assert not _aside(), "nothing is quarantined when the store parses"
    fresh = json.loads(_state().read_text(encoding="utf-8"))
    assert fresh["jobs"]["cluster:42"] == {"state": "RUNNING"}, (
        "an intact store is read, updated, and written back as before")
