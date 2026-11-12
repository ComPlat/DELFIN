"""The dashboard's Swarm tab starts a swarm the guided way.

Pick or create a spec, read what a swarm is and what each field means,
Check (problems first, then exactly what will happen), then one Start
button that submits the swarm job the same way as ``delfin-agent swarm
start``. The Refresh button is in the layout (an audit found it built but
never shown). Nothing here submits a job or starts a session.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin import scheduler_client
from delfin.agent.swarm import job
from delfin.dashboard import tab_swarm

GOOD = ("goal: Ship it\n"
        "packages:\n"
        "  - name: core\n"
        "    owned: [\"core/**\"]\n"
        "    task: Build core\n")


def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.org",
                    "-c", "commit.gpgsign=false", *args],
                   cwd=repo, check=True, capture_output=True, text=True)


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    scheduler_client.reset()
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    (r / "README").write_text("x\n")
    _git(r, "add", "README")
    _git(r, "commit", "-q", "-m", "init")
    yield r
    scheduler_client.reset()


# --------------------------------------------------------------------------
# Help text
# --------------------------------------------------------------------------

def test_the_tab_explains_what_a_swarm_is_in_plain_english():
    text = tab_swarm.SWARM_INTRO_HTML
    for words in ("What is a swarm?", "packages", "builder", "reviewer",
                  "own git worktree", "only write the files its package owns",
                  "Nothing is merged or pushed", "one SLURM job on a compute "
                  "node", "at most 8 member sessions", "sbatch script"):
        assert words in text, words


def test_every_spec_field_has_help():
    fields = dict(tab_swarm.FIELD_HELP)
    for field in ("goal", "packages", "packages[].name", "packages[].owned",
                  "packages[].task", "packages[].roles",
                  "limits.max_sessions", "limits.start_spacing_s",
                  "limits.memory_gate"):
        assert fields.get(field), field
    html = tab_swarm.field_help_html()
    assert "<table>" in html and "packages[].owned" in html
    # Every key the example spec uses is explained.
    import yaml
    from delfin.agent.swarm.launch import EXAMPLE_SPEC
    example = yaml.safe_load(EXAMPLE_SPEC)
    for key in example["limits"]:
        assert f"limits.{key}" in fields, key
    for key in example["packages"][0]:
        assert f"packages[].{key}" in fields, key


# --------------------------------------------------------------------------
# Validation before submit
# --------------------------------------------------------------------------

def test_a_broken_spec_is_named_before_anything_starts(repo):
    pre = tab_swarm.check_start("goal: [", repo / "swarm.yaml", route_="job")
    assert not pre["ok"]
    html = tab_swarm.check_html(pre)
    assert "Cannot start yet" in html and "not valid YAML" in html
    assert "What happens" not in html


def test_too_many_members_for_one_job_are_named(repo):
    many = "goal: g\nlimits: {max_sessions: 24}\npackages:\n" + "".join(
        f"  - name: p{i}\n    owned: [\"a{i}/**\"]\n"
        f"    roles: {{builder: true, reviewer: true}}\n" for i in range(5))
    pre = tab_swarm.check_start(many, repo / "swarm.yaml", route_="job")
    assert not pre["ok"]
    assert "at most 8" in tab_swarm.check_html(pre)


def test_no_spec_file_named_is_a_message_not_a_crash(repo):
    pre = tab_swarm.check_start(GOOD, "", route_="job")
    assert not pre["ok"]
    assert "Name the spec file" in " ".join(pre["errors"])


def test_a_good_spec_says_exactly_what_start_will_do(repo):
    path = repo / "swarm.yaml"
    pre = tab_swarm.check_start(GOOD, path, route_="job", settings={})
    assert pre["ok"], pre["errors"]
    html = tab_swarm.check_html(pre)
    assert "What happens when you click Start" in html
    assert "Submit one SLURM job" in html and "core-builder" in html
    assert "Swarm Watch" in html
    # Not saved yet: the tab says so (the board reads the file).
    assert "not saved to" in html
    path.write_text(GOOD)
    pre = tab_swarm.check_start(GOOD, path, route_="job", settings={})
    assert not any("not saved" in w for w in pre["warnings"])


def test_start_submits_only_a_spec_that_passes_its_checks(repo):
    submitted = []

    def fake_submit(text, pre):
        submitted.append((text, pre["route"]))
        return {"job": "4242"}

    bad = tab_swarm.start_checked("goal: [", repo / "swarm.yaml",
                                  route_="job", submit=fake_submit)
    assert not bad["ok"] and submitted == []
    assert bad["message"].startswith("Not started")
    good = tab_swarm.start_checked(GOOD, repo / "swarm.yaml", route_="job",
                                   submit=fake_submit, settings={})
    assert good["ok"] and good["job"] == "4242"
    assert submitted == [(GOOD, "job")]
    assert "Submitted swarm job 4242" in good["message"]
    assert "delfin-agent swarm stop 4242" in good["message"]


def test_start_here_creates_the_worktrees_in_the_specs_repository(repo):
    seen = []

    def fake_start(spec, *, execute=False, repo=None):
        seen.append((execute, repo))
        return {"started": [{"name": "core-builder"}]}

    out = tab_swarm.start_checked(GOOD, repo / "swarm.yaml", route_="here",
                                  start_here=fake_start)
    assert out["ok"] and seen == [(True, str(repo.resolve()))]
    assert "core-builder" in out["message"]


def test_save_keeps_the_older_version_as_a_copy(repo):
    path = repo / "swarm.yaml"
    assert "Saved" in tab_swarm.save_spec_text(GOOD, path)
    newer = GOOD.replace("Ship it", "Ship it twice")
    msg = tab_swarm.save_spec_text(newer, path)
    assert "kept as" in msg
    assert path.read_text() == newer
    copies = list(repo.glob("swarm.yaml.*.bak"))
    assert len(copies) == 1 and copies[0].read_text() == GOOD
    with pytest.raises(Exception):
        tab_swarm.save_spec_text("goal: [", path)
    assert path.read_text() == newer


def test_find_specs_looks_in_one_folder_only(tmp_path):
    (tmp_path / "swarm.yaml").write_text(GOOD)
    (tmp_path / "swarm-b.yml").write_text(GOOD)
    (tmp_path / "deep").mkdir()
    (tmp_path / "deep" / "swarm.yaml").write_text(GOOD)
    found = tab_swarm.find_specs(tmp_path)
    assert [Path(p).name for p in found] == ["swarm-b.yml", "swarm.yaml"]


def test_the_job_status_line_is_plain_words():
    assert "No swarm job runs" in tab_swarm.job_status_html([])
    html = tab_swarm.job_status_html([{
        "job": "7", "goal": "<b>x</b>", "state": "RUNNING", "node": "n1",
        "left": "1:00:00", "members": ["a"], "progress": {"started": ["a"]}}])
    assert "running on n1" in html and "<b>x</b>" not in html


# --------------------------------------------------------------------------
# The widgets
# --------------------------------------------------------------------------

@pytest.fixture()
def tab(repo, monkeypatch):
    pytest.importorskip("ipywidgets")
    monkeypatch.chdir(repo)
    monkeypatch.setattr(job, "route", lambda *a, **k: "here")

    def never(*_a, **_k):
        raise AssertionError("a test must never submit a job")

    monkeypatch.setattr(job, "submit", never)
    ctx = SimpleNamespace(swarm_spec_path=None)
    widget, refs = tab_swarm.create_tab(ctx)
    yield widget, refs, ctx
    refs["state"]["_stopped"] = True
    timer = refs["state"].get("_timer")
    if timer is not None:
        timer.cancel()


def test_the_refresh_button_is_shown(tab):
    widget, refs, _ctx = tab
    assert refs["refresh_btn"] in refs["header"].children
    assert refs["header"] in widget.children


def test_start_stays_off_until_the_check_passes_and_an_edit_turns_it_off(tab, repo):
    _widget, refs, _ctx = tab
    start = refs["start_btn"]
    assert start.disabled
    refs["path_box"].value = str(repo / "swarm.yaml")
    refs["spec_area"].value = "goal: ["
    refs["check_btn"].click()
    assert start.disabled
    assert "Cannot start yet" in refs["check_out"].value
    refs["spec_area"].value = GOOD
    refs["check_btn"].click()
    assert not start.disabled, refs["check_out"].value
    assert "What happens when you click Start" in refs["check_out"].value
    refs["spec_area"].value = GOOD + "# edited\n"
    assert start.disabled


def test_a_pause_verb_the_package_offers_gets_a_button_that_calls_it(repo, monkeypatch):
    """Hunt finding: the Pause/Resume buttons tested the wrong verb names
    and had no handler -- a button that does nothing. A verb the package
    exposes is shown AND wired; one it does not expose is not shown."""
    pytest.importorskip("ipywidgets")
    import delfin.agent.swarm as swarm_pkg

    called = []
    monkeypatch.setattr(swarm_pkg, "pause_all", lambda: called.append("pause") or "2 paused",
                        raising=False)
    monkeypatch.chdir(repo)
    monkeypatch.setattr(job, "route", lambda *a, **k: "here")
    widget, refs = tab_swarm.create_tab(SimpleNamespace(swarm_spec_path=None))
    try:
        buttons = [w for w in _walk(widget) if getattr(w, "description", "") in ("Pause all", "Resume")]
        assert [b.description for b in buttons] == ["Pause all"]
        buttons[0].click()
        assert called == ["pause"]
        assert "2 paused" in refs["flow_status"].value
    finally:
        refs["state"]["_stopped"] = True
        timer = refs["state"].get("_timer")
        if timer is not None:
            timer.cancel()


def _walk(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _walk(child)
