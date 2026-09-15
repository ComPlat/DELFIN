"""The CI a push starts is watched, like a cluster job.

After report 20260915-085107 pushed to main, CI went red and nobody looked:
the agent had no way to read it -- ``gh`` is not installed on the cluster
nodes -- and when the user pasted the log it guessed at a licence header.
A pushed commit's GitHub Actions runs are now a watched job: LLM-free,
asked at most once a minute, reported once, with the job and step that
failed and a link to the run.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import job_monitor as jm
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

_SHA = "269429eb"
_FULL = _SHA + "0" * 32


def _run(name, status="completed", conclusion="success", sha=_FULL, rid=1):
    return {"name": name, "status": status, "conclusion": conclusion,
            "head_sha": sha,
            "html_url": f"https://github.com/ComPlat/DELFIN/actions/runs/{rid}",
            "jobs_url": f"https://api.github.com/jobs/{rid}"}


class _Fetch:
    def __init__(self, pages):
        self.pages = pages
        self.urls: list[str] = []

    def __call__(self, url):
        self.urls.append(url)
        for key, page in self.pages.items():
            if key in url:
                return page
        return None


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    jm.register_ci_watch(tmp_path, "ComPlat/DELFIN", _SHA, "main")
    return tmp_path


def _watch_file(ws):
    return ws / ".delfin" / "agent_watched_jobs.json"


def _check(ws, fetch):
    """One check, with the once-a-minute limit taken out of the way."""
    data = jm.load_watched(_watch_file(ws))
    for entry in data["jobs"].values():
        entry["last_checked"] = 0
    jm.save_watched(data, _watch_file(ws))
    return jm.check_agent_jobs(ws, fetch_fn=fetch)


def test_a_ci_id_is_its_own_kind(ws):
    assert jm._classify_job_id(f"ci:ComPlat/DELFIN@{_SHA}", ws) == "ci"
    entry = jm.load_watched(_watch_file(ws))["jobs"][f"ci:ComPlat/DELFIN@{_SHA}"]
    assert entry["kind"] == "ci" and entry["branch"] == "main"


def test_a_running_ci_is_not_reported(ws):
    fetch = _Fetch({"actions/runs": {"workflow_runs": [
        _run("CI", status="in_progress", conclusion=None)]}})
    assert _check(ws, fetch) == []
    assert jm.load_watched(_watch_file(ws))["jobs"]


def test_a_green_ci_is_reported_once(ws):
    fetch = _Fetch({"actions/runs": {"workflow_runs": [
        _run("CI"), _run("CodeQL", rid=2)]}})
    done = _check(ws, fetch)
    assert len(done) == 1
    assert done[0]["kind"] == "ci" and done[0]["ok"] is True
    assert done[0]["state"] == "SUCCESS"
    assert _check(ws, fetch) == []


def test_a_red_ci_names_the_job_and_the_step_that_failed(ws):
    fetch = _Fetch({
        "jobs/7": {"jobs": [
            {"name": "tests (py3.11)", "conclusion": "failure", "steps": [
                {"name": "Checkout", "conclusion": "success"},
                {"name": "Run fast test suite", "conclusion": "failure"}]},
            {"name": "lint (ruff)", "conclusion": "success", "steps": []},
        ]},
        "actions/runs": {"workflow_runs": [
            _run("CI", conclusion="failure", rid=7), _run("CodeQL", rid=8)]},
    })
    done = _check(ws, fetch)
    assert len(done) == 1
    assert done[0]["ok"] is False and done[0]["state"] == "FAILURE"
    assert done[0]["signatures"] == ["CI › tests (py3.11) › Run fast test suite"]
    assert done[0]["url"].endswith("/runs/7")


def test_the_runs_of_another_commit_are_not_this_commits(ws):
    fetch = _Fetch({"actions/runs": {"workflow_runs": [
        _run("CI", conclusion="failure", sha="d507039a" + "0" * 32)]}})
    assert _check(ws, fetch) == []


def test_github_out_of_reach_is_said_once_and_not_as_green(ws):
    fetch = _Fetch({})
    done = _check(ws, fetch)
    assert len(done) == 1
    assert done[0]["state"] == jm.STATE_UNAVAILABLE and done[0]["ok"] is False
    assert "not green" in done[0]["degraded"]
    assert _check(ws, fetch) == []
    assert jm.load_watched(_watch_file(ws))["jobs"], "the watch goes on"


def test_a_push_that_never_got_a_run_is_reported_after_half_an_hour(ws):
    data = jm.load_watched(_watch_file(ws))
    for entry in data["jobs"].values():
        entry["added_at"] = time.time() - 31 * 60
    jm.save_watched(data, _watch_file(ws))
    done = _check(ws, _Fetch({"actions/runs": {"workflow_runs": []}}))
    assert len(done) == 1
    assert done[0]["state"] == "NO CI RUN" and done[0]["ok"] is False


def test_github_is_asked_at_most_once_a_minute(ws):
    fetch = _Fetch({"actions/runs": {"workflow_runs": [
        _run("CI", status="queued", conclusion=None)]}})
    _check(ws, fetch)
    jm.check_agent_jobs(ws, fetch_fn=fetch)
    assert len(fetch.urls) == 1


def test_watch_job_takes_a_ci_id(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    out = json.loads(_DocToolExecutor().execute(
        "watch_job", {"job_id": f"ci:ComPlat/DELFIN@{_SHA}"}, perms))
    assert out.get("status") == "watching", out
    assert out.get("kind") == "ci"
