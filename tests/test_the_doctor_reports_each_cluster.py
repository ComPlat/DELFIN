"""``delfin doctor`` reports each of DELFIN's clusters: usable, logged in, budget left.

Input: a home under ``tmp_path`` with ctun's config at its default place,
DELFIN's ``remote`` settings as a dict, and a fake ``ssh`` that answers
``ssh -O check`` from a marker file. Output: ``doctor.check_clusters``.

Semantics: no cluster configured, no line; a cluster DELFIN refuses is
"broken" with its problems; usable without a login is "missing" with the
login command; logged in is "ok". Only ``budget=True`` probes the budget
over the tunnel (the doctor otherwise opens no network connection), and a
failed probe is "broken". Inside an agent session the login is not checked.
"""

from __future__ import annotations

import sys

import pytest

from cluster_tunnel import api as ctun_api
from delfin import doctor

_FAKE_SSH = """\
#!{python}
import os, sys
args = sys.argv[1:]
if "-O" in args and os.path.exists(os.path.join(os.environ["FAKE_SSH_STATE"], "master")):
    sys.exit(0)
sys.exit(255)
"""

SETTINGS = {"remote": {"clusters": {"kcist": {}}}}


@pytest.fixture
def home(tmp_path, monkeypatch):
    root = tmp_path / "home"
    ctun_dir = root / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    for key, value in (("HOME", root), ("XDG_CONFIG_HOME", root / ".config"),
                       ("XDG_CACHE_HOME", root / ".cache"), ("FAKE_SSH_STATE", tmp_path)):
        monkeypatch.setenv(key, str(value))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "ssh").write_text(_FAKE_SSH.format(python=sys.executable), encoding="utf-8")
    (bin_dir / "ssh").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:/usr/bin:/bin")
    (ctun_dir / "config.yaml").write_text(
        "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {ctun_dir / 'known_hosts'}\n", encoding="utf-8")
    return tmp_path


def _one(**kwargs) -> doctor.CheckResult:
    (result,) = doctor.check_clusters(settings=SETTINGS, **kwargs)
    assert result.name == "cluster kcist"
    return result


def test_no_configured_cluster_adds_no_line(home):
    assert doctor.check_clusters(settings={}) == []


def test_a_refused_cluster_is_broken_with_its_problems(home):
    cfg = home / "home" / ".config" / "cluster-tunnel" / "config.yaml"
    cfg.write_text("clusters:\n  kcist:\n    host: login.example\n", encoding="utf-8")
    result = _one()
    assert result.status == doctor.BROKEN
    assert "hostkey --pin" in result.detail and "ctun config" in result.fix_hint


def test_usable_without_a_login_is_missing_with_the_login_command(home):
    result = _one()
    assert result.status == doctor.MISSING
    assert result.fix_hint == "ctun login -t kcist"


def test_logged_in_is_ok_and_the_budget_is_probed_only_when_asked(home, monkeypatch):
    (home / "master").touch()
    probed = []

    def budget(name, **kwargs):
        probed.append(name)
        return ctun_api.Decision(allowed=True, reason="", used=2.5, limit=100.0, unit="jobh")

    monkeypatch.setattr(ctun_api, "budget", budget)
    result = _one()
    assert result.status == doctor.OK and "logged in (me@login.example)" in result.detail
    assert probed == []
    result = _one(budget=True)
    assert result.status == doctor.OK and "budget 2.5 of 100.0 jobh used" in result.detail
    assert probed == ["kcist"]


def test_a_failed_budget_probe_is_broken(home, monkeypatch):
    (home / "master").touch()

    def budget(name, **kwargs):
        raise ctun_api.TransportError("sacct: command not found")

    monkeypatch.setattr(ctun_api, "budget", budget)
    result = _one(budget=True)
    assert result.status == doctor.BROKEN and "sacct" in result.detail


def test_inside_an_agent_session_the_login_is_not_checked(home, monkeypatch):
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "1")
    result = _one()
    assert result.status == doctor.OK and "not checked" in result.detail
