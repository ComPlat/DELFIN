"""A machine with no cluster tunnel is not a machine with a defect.

``check_clusters`` imports ``delfin.remote``, which imports the cluster
tunnel package. That package is optional: a workstation that never submits
to a cluster is complete without it. Before this, the import raised,
``run_all`` caught it in its defensive handler and produced

    clusters | broken | check crashed: No module named 'cluster_tunnel'
             | report this as a DELFIN bug

so the person was told to file a bug about software they had chosen not to
install, and ``exit_code`` returned 1 -- the doctor's "something is broken"
answer -- on a healthy machine.

Measured on 2026-11-19: on this machine three tests in
``tests/test_delfin_doctor.py`` failed for that reason, and they had been
failing on every machine without the package since the cluster check was
added. They asserted the number of rows and that none is broken, and the
cluster check reads the LIVE machine's ``remote.clusters`` -- so the result
set under test depended on whether whoever ran the tests happened to have a
cluster configured, and logged in. ``run_all`` now takes the same
``settings`` / ``ctun_config`` seam ``check_clusters`` already had, and
those tests pass ``settings={}``.

The import is forced to fail here with ``sys.modules``, so these tests
measure the same thing whether or not the package is installed.
"""

from __future__ import annotations

import sys

import pytest

from delfin.doctor import BROKEN, MISSING, check_clusters, exit_code, run_all


@pytest.fixture
def no_tunnel(monkeypatch):
    """Make ``from delfin.remote import clusters`` raise ImportError.

    A None entry in sys.modules is what the import system treats as "this
    module is not importable", which is exactly the state of a machine
    without the package -- and it is independent of this machine's own.
    """
    monkeypatch.setitem(sys.modules, "delfin.remote.clusters", None)
    monkeypatch.setitem(sys.modules, "delfin.remote.connection", None)


def test_no_cluster_configured_says_nothing_at_all(no_tunnel):
    """Nothing to install, nothing to act on, so no row."""
    assert check_clusters(settings={}) == []


def test_a_configured_cluster_reports_the_missing_package_not_a_crash(
        no_tunnel):
    rows = check_clusters(settings={"remote": {"clusters": {"hpc": {}}}})

    assert len(rows) == 1, rows
    row = rows[0]
    assert row.status == MISSING, (
        "an optional package that is not installed is missing, not broken: "
        "broken is what makes the doctor exit 1 and tells the person to "
        "report a DELFIN bug")
    assert "hpc" in row.detail, "say which cluster needs it"
    assert "not installed" in row.detail
    assert row.fix_hint, "and say what to do about it"


def test_a_healthy_machine_without_the_package_still_exits_zero(
        no_tunnel, tmp_path):
    """The whole point: the doctor's verdict must not turn on a package the
    person chose not to install."""
    results = run_all(scratch_dir=tmp_path, settings={})

    broken = [(r.name, r.detail) for r in results if r.status == BROKEN]
    assert not broken, broken
    assert exit_code(results) == 0


def test_run_all_passes_the_cluster_source_through(tmp_path, monkeypatch):
    """The seam the tests depend on: what run_all is given is what the
    cluster check is asked. Without it no test could keep the live
    machine's configuration out of its result set."""
    seen = {}

    def _spy(*, budget=False, settings=None, ctun_config=None):
        seen.update(budget=budget, settings=settings, ctun_config=ctun_config)
        return []

    monkeypatch.setattr("delfin.doctor.check_clusters", _spy)
    cfg = tmp_path / "ctun.toml"

    run_all(scratch_dir=tmp_path, settings={"remote": {}}, ctun_config=cfg,
            cluster_budget=True)

    assert seen == {"budget": True, "settings": {"remote": {}},
                    "ctun_config": cfg}


def test_an_unrelated_import_error_is_still_a_crash(monkeypatch, tmp_path):
    """The handler must not swallow a real defect. Only the cluster-tunnel
    imports are allowed to be absent; an ImportError from inside the check
    is a broken check and says so."""
    def boom(**_kw):
        raise ImportError("no module named 'something_delfin_needs'")

    monkeypatch.setattr("delfin.doctor.check_clusters", boom)

    results = run_all(scratch_dir=tmp_path, settings={})

    crashed = [r for r in results if r.status == BROKEN]
    assert crashed, "a check that raises must surface as broken"
    assert any("something_delfin_needs" in r.detail for r in crashed)
