"""pytest plugin: run one shard of the fast suite and report what ran.

Input
    Loaded with ``-p shard_plugin`` (this directory on PYTHONPATH) and the
    options --shard K/N, --shard-report DIR and optionally
    --shard-durations PATH (default: test-durations.json next to this
    file). CI passes all of it through PYTEST_ADDOPTS, so the workflow's
    command stays ``pytest -q -m "not slow" ...``.

Output
    The pytest run over the files ci_shards.plan() assigns to shard K, and
    in DIR: universe.txt (every node id collected after ``-m``
    deselection), selected.txt (this shard's node ids, in run order),
    ran.txt (node ids whose test protocol finished, in order), meta.json
    (shard, shards, commit, pytest exit status, plan digest) and
    durations.json (measured seconds per file: setup + call + teardown).

Semantics
    The selection runs in a wrapper around pytest_collection_modifyitems,
    after every other implementation -- the conftest's slow marking and
    pytest's own ``-m`` deselection included -- so the universe is exactly
    the set an unsharded ``pytest -m "not slow"`` runs. Files outside the
    shard are reported through pytest_deselected, as ``-m`` does.
    universe.txt and selected.txt are written right after collection;
    ran.txt and meta.json only at session end, so a run that dies leaves
    no meta.json and the aggregator fails.

    At configure time, before the first test, the plugin removes
    PYTEST_ADDOPTS and its own directory from PYTHONPATH in os.environ: a
    test that starts pytest or python in a subprocess sees the environment
    an unsharded run gives it. Without --shard the plugin does nothing.

Rejected alternatives
    A wrapper script (``python .github/ci/run_shard.py ...``) as the run
    line: delfin.agent.verify_recipe reads the CI run lines as this
    project's test command, and the first one would no longer be the
    ``pytest -m "not slow"`` command an agent may run on changed files.
    Shard options on the run line: the same command, copied for a local
    run, would fail without the plugin on the path.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import ci_shards  # noqa: E402


def pytest_addoption(parser):
    g = parser.getgroup("shard", "fast-suite sharding (CI)")
    g.addoption("--shard", default=None, metavar="K/N",
                help="run shard K of N (1-based)")
    g.addoption("--shard-report", default=None, metavar="DIR",
                help="where the shard writes its report")
    g.addoption("--shard-durations", default=str(ci_shards.DURATIONS),
                metavar="PATH", help="seconds per test file")


def _scrub_environment() -> None:
    os.environ.pop("PYTEST_ADDOPTS", None)
    parts = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep)
             if p and Path(p).resolve() != HERE]
    if parts:
        os.environ["PYTHONPATH"] = os.pathsep.join(parts)
    else:
        os.environ.pop("PYTHONPATH", None)


class _Shard:
    def __init__(self, k: int, n: int, report: Path, durations: Path):
        self.k, self.n, self.report = k, n, report
        self.durations = ci_shards.load_durations(durations)
        self.digest = ""
        self.ran: list[str] = []
        self.seconds: dict[str, float] = {}

    def _write(self, name: str, text: str) -> None:
        tmp = self.report / (name + ".tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(self.report / name)

    @pytest.hookimpl(hookwrapper=True)
    def pytest_collection_modifyitems(self, session, config, items):
        yield
        universe = [it.nodeid for it in items]
        assignment = ci_shards.plan(
            {ci_shards.file_of(n) for n in universe}, self.durations, self.n)
        self.digest = ci_shards.plan_digest(assignment)
        keep = [it for it in items
                if assignment[ci_shards.file_of(it.nodeid)] == self.k]
        drop = [it for it in items
                if assignment[ci_shards.file_of(it.nodeid)] != self.k]
        if drop:
            config.hook.pytest_deselected(items=drop)
        items[:] = keep
        self.report.mkdir(parents=True, exist_ok=True)
        self._write("universe.txt", "".join(n + "\n" for n in sorted(universe)))
        self._write("selected.txt", "".join(it.nodeid + "\n" for it in keep))

    def pytest_runtest_logreport(self, report):
        f = ci_shards.file_of(report.nodeid)
        self.seconds[f] = self.seconds.get(f, 0.0) + float(report.duration)

    def pytest_runtest_logfinish(self, nodeid, location):
        self.ran.append(nodeid)

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session, exitstatus):
        self.report.mkdir(parents=True, exist_ok=True)
        self._write("ran.txt", "".join(n + "\n" for n in self.ran))
        self._write("durations.json", json.dumps(
            {f: round(s, 3) for f, s in sorted(self.seconds.items())},
            indent=0))
        self._write("meta.json", json.dumps({
            "shard": self.k,
            "shards": self.n,
            "sha": os.environ.get("GITHUB_SHA", ""),
            "exitstatus": int(exitstatus),
            "plan_digest": self.digest,
        }, indent=1))


def pytest_configure(config):
    spec = config.getoption("--shard")
    if not spec:
        return
    k, _, n = spec.partition("/")
    try:
        k, n = int(k), int(n)
    except ValueError:
        raise pytest.UsageError(f"--shard expects K/N, got {spec!r}")
    if not 1 <= k <= n:
        raise pytest.UsageError(f"--shard {spec}: K must be in 1..N")
    report = config.getoption("--shard-report")
    if not report:
        raise pytest.UsageError("--shard needs --shard-report DIR")
    _scrub_environment()
    config.pluginmanager.register(
        _Shard(k, n, Path(report).resolve(),
               Path(config.getoption("--shard-durations"))),
        "ci-shard")
