"""The dashboard's approval runner is contained on every host.

Audit 2026-09-16: the runner took bubblewrap whenever it was installed,
working or not; without bubblewrap or firejail its allowlist mode ran
commands with nothing around them (python, pip and env are on the list)
and with every key in the environment; a timeout killed only bash.
"""
import pytest

from delfin.agent import api_client as A
from delfin.agent import sandbox as SB


def test_an_installed_but_broken_bwrap_is_not_chosen(monkeypatch):
    monkeypatch.setattr(SB.shutil, "which", lambda n: "/usr/bin/" + n if n == "bwrap" else None)
    monkeypatch.setattr(SB, "_bwrap_works", lambda: False)
    monkeypatch.delenv("DELFIN_AGENT_SANDBOX", raising=False)
    assert SB.detect_config().mode == "allowlist"


def test_the_runner_hands_no_provider_key_to_a_command(tmp_path, monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "kit-key-in-the-environment")
    res = SB.run_agent_command("env", tmp_path,
                               config=SB.SandboxConfig("allowlist", False, 30))
    assert not res.blocked, res.block_reason
    assert "PATH=" in res.stdout
    assert "kit-key-in-the-environment" not in res.stdout


@pytest.mark.skipif(not A._landlock_functional(), reason="no Landlock on this kernel")
def test_the_allowlist_mode_runs_under_landlock(tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "id_ed25519").write_text("PRIVATE-KEY-PROBE\n")
    monkeypatch.setenv("HOME", str(home))
    repo = tmp_path / "repo"
    repo.mkdir()
    argv = SB._unsandboxed_argv(f"cat {home}/.ssh/id_ed25519; echo x > {tmp_path}/out", repo, "allowlist")
    from delfin.agent import contained_run
    out = contained_run.run(argv, cwd=str(repo), timeout=60)
    assert "PRIVATE-KEY-PROBE" not in out.stdout
    assert not (tmp_path / "out").exists()
    assert SB._unsandboxed_argv("true", repo, "off") == ["bash", "-c", "true"]
