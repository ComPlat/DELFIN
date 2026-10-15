"""sandbox._firejail_argv: the whitelist names the repository and nothing else.

Firejail expands file globs (``*``, ``?``, ``[...]``) and macros
(``${HOME}`` and others) in ``--whitelist=`` and ``--chdir=`` values. The
repository path was written there unescaped, so a path holding one of
those characters whitelisted every path it matched. Such a path is
refused; the command is not run.

Input: a repository path. Output: ``ValueError`` from ``_firejail_argv``;
a blocked ``RunResult`` with the reason from ``run_agent_command``.
"""
from __future__ import annotations

import pytest

from delfin.agent import sandbox


def _whitelist(argv):
    return [a.split("=", 1)[1] for a in argv if a.startswith("--whitelist=")]


def test_a_plain_repository_is_whitelisted_as_its_resolved_path(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    argv = sandbox._firejail_argv("ls", repo, allow_network=False)
    assert _whitelist(argv) == [str(repo.resolve())]


@pytest.mark.parametrize("name", ["re*", "re?o", "re[p]o", "x${HOME}", "a$b"])
def test_a_glob_or_macro_character_in_the_repository_is_refused(tmp_path, name):
    repo = tmp_path / name
    repo.mkdir()
    (tmp_path / "repo-other").mkdir()
    # Before: --whitelist=<tmp>/re* matched <tmp>/repo-other as well.
    with pytest.raises(ValueError) as exc:
        sandbox._firejail_argv("ls", repo, allow_network=False)
    assert str(repo.resolve()) in str(exc.value)


def test_a_control_character_in_the_repository_is_refused(tmp_path):
    repo = tmp_path / "a\nb"
    repo.mkdir()
    with pytest.raises(ValueError):
        sandbox._firejail_argv("ls", repo, allow_network=False)


def test_the_run_is_blocked_with_the_reason_and_nothing_runs(tmp_path, monkeypatch):
    repo = tmp_path / "re*"
    repo.mkdir()
    ran = []
    monkeypatch.setattr("delfin.agent.contained_run.run",
                        lambda *a, **k: ran.append(a))
    cfg = sandbox.SandboxConfig(mode="firejail", allow_network=False,
                                timeout_s=10)
    res = sandbox.run_agent_command("ls", repo, config=cfg)
    assert res.blocked is True
    assert "firejail" in (res.block_reason or "")
    assert str(repo.resolve()) in (res.block_reason or "")
    assert ran == []
