"""The launch directory is a security decision, so it is checked first.

The pairing in here is deliberate and both halves have to stay. The
constructor floor in `KitToolPermissions.__post_init__` is the GUARANTEE —
it holds for the dashboard, the benchmark and every headless caller.
`launch_guard` is the manners: the same refusal, one layer earlier, as a
sentence rather than a traceback out of a constructor the user never
called. A later change that deletes the constructor check because
"launch_guard already does it" would remove the guarantee and keep the
manners, so there is a test here that fails if that happens.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import launch_guard as lg


# --------------------------------------------------------------------------
# The floor
# --------------------------------------------------------------------------

def test_home_is_refused_in_prose_and_not_as_a_traceback(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(fake_home)
    assert report.refused
    text = report.render()
    assert "Traceback" not in text
    assert str(fake_home) in text
    assert "cd " in text, "the refusal has to say what to do instead"


def test_a_system_root_is_refused():
    assert lg.inspect_launch_dir("/etc").refused
    assert lg.inspect_launch_dir("/").refused


def test_the_constructor_still_refuses_home_on_its_own(monkeypatch, tmp_path):
    """The guarantee, asserted where it actually lives.

    If someone deletes the __post_init__ check on the grounds that
    launch_guard covers it, every non-terminal caller loses the floor. This
    is the test that stops that.
    """
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    from delfin.agent.api_client import KitToolPermissions
    with pytest.raises(ValueError):
        KitToolPermissions(workspace=fake_home)


def test_a_project_under_home_is_fine(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(project, check_trust=False)
    assert not report.refused


def test_a_directory_that_is_not_there_is_refused(tmp_path):
    report = lg.inspect_launch_dir(tmp_path / "nope", check_trust=False)
    assert report.refused
    assert "no such directory" in report.render()


# --------------------------------------------------------------------------
# Grants
# --------------------------------------------------------------------------

def test_a_forbidden_add_dir_is_named_not_dropped(monkeypatch, tmp_path):
    """KitToolPermissions drops a bad extra dir silently.

    Silently is the problem: the session then has less reach than the
    command line asked for and nothing says so.
    """
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(
        project, add_dirs=(str(fake_home),), check_trust=False)
    assert report.refused
    assert str(fake_home) in report.render()
    assert not report.granted_dirs


def test_a_credential_directory_is_never_granted(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    secrets = fake_home / ".ssh"
    project.mkdir(parents=True)
    secrets.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    for kwargs in ({"add_dirs": (str(secrets),)}, {"read_dirs": (str(secrets),)}):
        report = lg.inspect_launch_dir(project, check_trust=False, **kwargs)
        assert report.refused, f"{kwargs} was allowed"
        assert "credential" in report.render()


def test_a_parent_of_the_workspace_is_refused_as_a_grant(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    outer = fake_home / "outer"
    project = outer / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(
        project, add_dirs=(str(outer),), check_trust=False)
    assert report.refused
    assert "widen" in report.render()


def test_a_good_grant_survives(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    sibling = fake_home / "library"
    project.mkdir(parents=True)
    sibling.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(
        project, add_dirs=(str(sibling),), check_trust=False)
    assert not report.refused
    assert report.granted_dirs == (sibling.resolve(),)


# --------------------------------------------------------------------------
# Notices
# --------------------------------------------------------------------------

def _git_repo(path: Path) -> None:
    subprocess.run(["git", "init", "-q"], cwd=path, check=True)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=path, check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=path, check=True)


def test_a_folder_with_no_git_is_told_so(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    report = lg.inspect_launch_dir(project, check_trust=False)
    codes = {f.code for f in report.findings}
    assert "not_a_git_repo" in codes
    assert not report.refused, "no git is a notice, not a refusal"


def test_uncommitted_work_is_counted_and_named(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))
    _git_repo(project)
    (project / "already_mine.txt").write_text("mine\n")

    report = lg.inspect_launch_dir(project, check_trust=False)
    assert report.git.is_repo
    assert "already_mine.txt" in " ".join(report.git.dirty)
    assert "already_mine.txt" in report.render()


def test_an_ephemeral_workspace_says_the_write_gate_is_inert(monkeypatch, tmp_path):
    # $HOME has to sit OUTSIDE /tmp here. pytest's tmp_path is itself under
    # /tmp, so a fake home inside it would make /tmp an ancestor of home and
    # the workspace would be refused as a forbidden root — which is correct
    # behaviour for that arrangement and simply not the case under test.
    monkeypatch.setattr(
        Path, "home", classmethod(lambda cls: Path("/nonexistent/home")))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    report = lg.inspect_launch_dir(scratch, check_trust=False)
    codes = {f.code for f in report.findings}
    assert not report.refused
    assert "ephemeral_root" in codes, (
        f"{scratch} is under /tmp, where writes skip the path gate")


def test_untrusted_definitions_are_reported_and_asked_about(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))
    monkeypatch.setattr(
        "delfin.agent.workspace_trust.pending_notices",
        lambda ws: ["2 hook commands are withheld"])

    report = lg.inspect_launch_dir(project)
    assert report.trust_notices == ("2 hook commands are withheld",)
    asked = {f.code for f in report.questions}
    assert "untrusted_workspace" in asked
    assert "withheld" in report.render()


def test_a_broken_trust_store_does_not_stop_the_agent(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    project = fake_home / "project"
    project.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    def _boom(ws):
        raise RuntimeError("trust store is corrupt")

    monkeypatch.setattr(
        "delfin.agent.workspace_trust.pending_notices", _boom)
    report = lg.inspect_launch_dir(project)
    assert not report.refused


# --------------------------------------------------------------------------
# Posture
# --------------------------------------------------------------------------

def test_a_persisted_bypass_is_not_inherited():
    """The test that must never be deleted.

    default_mode: bypassPermissions is written once, often in a container,
    and then makes every future session on the machine unattended with
    nothing on screen saying so.
    """
    mode, why = lg.resolve_posture(
        flag_mode="", persisted_mode="bypassPermissions")
    assert mode == "plan"
    assert "settings.json" in why
    assert "unattended" in why


def test_bypass_on_the_command_line_needs_the_second_flag():
    mode, _ = lg.resolve_posture(
        flag_mode="bypassPermissions", unattended_opt_in=False)
    assert mode == "plan"
    mode, why = lg.resolve_posture(
        flag_mode="bypassPermissions", unattended_opt_in=True)
    assert mode == "bypassPermissions"
    assert why == "--permission-mode"


def test_the_flag_beats_the_settings_file():
    mode, why = lg.resolve_posture(
        flag_mode="acceptEdits", persisted_mode="plan")
    assert mode == "acceptEdits"
    assert why == "--permission-mode"


def test_an_ordinary_persisted_mode_is_honoured():
    for persisted in ("plan", "default", "acceptEdits"):
        mode, why = lg.resolve_posture(persisted_mode=persisted)
        assert mode == persisted
        assert "settings.json" in why


def test_nothing_configured_starts_read_only():
    mode, why = lg.resolve_posture()
    assert mode == "plan"
    assert why == "default"
