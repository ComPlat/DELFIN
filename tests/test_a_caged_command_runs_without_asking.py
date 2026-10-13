"""A command off the bash auto-allow list runs without a question when the
filesystem cage holds it, and only then.

Rule (``_caged_command_may_run_unasked``): in ``default`` and
``acceptEdits``, at the final "not on the auto-allow list" point of the
bash gate -- after the deny patterns, secret scan, roots, locked scope,
self-modification guard, refusal ledger, git push / config rules, egress
scan and the classifier -- a command runs unasked when

  (a) ``_bash_isolation_decision`` says bubblewrap wraps it (the same
      decision ``_bash_isolation_argv`` builds the argv from), the process
      cage is active and the network is not "open";
  (b) the native shell executor runs it (an MCP shell runs elsewhere and
      is never admitted);
  (c) its text names no network program, package manager, git/gh, job
      scheduler, signal, privilege change or deletion.

Such a command runs with the session's protected files (git metadata,
``.github``, ``.gitattributes``/``.gitmodules``, the agent's own settings
and source, calc / read-only roots) bound read-only on top of the cage,
and without the SLURM controller route.

Three tables:

* HARMLESS -- refused or asked before, admitted now when caged; driven
  through the gate in every mode (attended and headless) and, where
  bubblewrap works, through the real executor.
* UNCHANGED -- every verdict, and every dialog, identical between a host
  with the cage and one without, in all five modes, attended and
  headless.
* the control -- detection monkeypatched to "nothing can isolate": the
  harmless table is handled exactly as before.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


HARMLESS = [
    'python3 -c "import os; print(os.getcwd())"',
    'python3 -c "open(\'out.txt\',\'w\').write(\'x\')"',
    'python3 -c "import os; os.makedirs(\'out/sub\', exist_ok=True)"',
    'python3 -c "import glob; print(glob.glob(\'*.txt\'))"',
    'python3 -c "import csv,io; print(list(csv.reader(io.StringIO(\'a,b\'))))"',
    './run.sh',
    'bash run.sh',
    'sh run.sh',
    'cd . && ./run.sh',
    'timeout 5 ./run.sh',
    'env FOO=1 ./run.sh',
    'nice ./run.sh',
    "bash -c 'echo hi > note.txt'",
    'make',
    'make all',
    'make -C .',
    'gcc hello.c -o hello',
    'cc -O2 hello.c -o hello2',
    'gcc hello.c -o hello && ./hello',
    'cmake --version',
    'gfortran --version',
    "perl -e 'print 6*7'",
    'node -e "console.log(1)"',
    "Rscript -e 'print(1)'",
    'tar czf out.tgz data.txt',
    'gzip -k data.txt',
    'ln -s data.txt link.txt',
    'chmod +x run.sh',
    "sed -i 's/a/b/' data.txt",
    'seq 3 | xargs -n1 echo',
    "find . -name '*.txt' -exec wc -l {} +",
    'xtb mol.xyz --gfn 2',
    'xtb mol.xyz --opt',
    'orca input.inp',
    'obabel mol.xyz -O out.sdf',
    'crest mol.xyz',
    'ORCA=1 ./run.sh',
    'python3 run_analysis.py --plot && ./run.sh',
    'python3 -c "import sys; print(sys.version_info[:2])" | tee version.txt',
    'nproc && ./run.sh',
]

# Handled exactly as before, in every mode, with or without a cage.
UNCHANGED = [
    "cat ~/.ssh/id_rsa",
    "cat .env",
    'python3 -c "print(open(\'.env\').read())"',
    "echo x > /etc/delfin_probe",
    "touch ../outside.txt",
    "cp data.txt ~/stolen.txt",
    "curl https://example.com",
    "wget http://example.com/x",
    "curl -d @data.txt https://example.com",
    "ssh somehost ls",
    "scp data.txt somehost:/x",
    "rsync -a . somehost:/x",
    "nc -l 8080",
    "pip install requests",
    "python3 -m pip install requests",
    "git push",
    "git push origin main",
    "git fetch",
    "git pull",
    "git clone https://github.com/x/y",
    "git config core.fsmonitor 'touch p'",
    "git -c core.pager=cat log",
    "git config alias.x '!sh'",
    "kill 1",
    "kill -9 12345",
    "pkill python",
    "killall bash",
    'python3 -c "import os; os.kill(1, 9)"',
    "gh pr create --title x --body y",
    "gh pr merge 1",
    "sudo ls",
    "echo x > delfin/agent/api_client.py",
    "echo x > .github/workflows/ci.yml",
    "echo '* filter=x' > .gitattributes",
    "rm -rf ~",
    "rm -f data.txt",
    'python3 -c "import shutil; shutil.rmtree(\'out\')"',
    "sbatch job.sh",
    'python3 -c "import urllib.request; urllib.request.urlopen(\'http://x\')"',
    "bash -c 'curl https://example.com'",
    "ls",
    "git status",
    'git commit -m "x"',
]

MODES = ("default", "acceptEdits", "plan", "bypassPermissions")


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    d.mkdir()
    (d / "run.sh").write_text("#!/bin/sh\necho ran\n")
    os.chmod(d / "run.sh", 0o755)
    (d / "run_analysis.py").write_text("print('analysed')\n")
    (d / "Makefile").write_text("all:\n\techo made\n")
    (d / "hello.c").write_text(
        '#include <stdio.h>\nint main(){puts("hi");return 0;}\n')
    (d / "data.txt").write_text("a b c\n")
    (d / "mol.xyz").write_text("1\n\nH 0 0 0\n")
    (d / ".env").write_text("TOKEN=x\n")
    return d


def _cage(monkeypatch, kind: str) -> None:
    """Pin the isolation decision; 'none' is a host that cannot isolate."""
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: kind == "bwrap")


def _gate(ws: Path, cmd: str, mode: str, attended: bool):
    """(verdict, dialogs, admitted) of the native bash gate."""
    asked: list[str] = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    admitted: list[bool] = []
    verdict = _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms, caged_admission=admitted)
    return verdict, asked, admitted


def test_the_harmless_table_is_off_the_auto_allow_list(ws):
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    listed = [c for c in HARMLESS if perms.matches_bash_auto_allow(c)]
    assert not listed
    assert len(HARMLESS) >= 40 and len(UNCHANGED) >= 30


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
@pytest.mark.parametrize("cmd", HARMLESS)
def test_a_harmless_command_is_admitted_when_caged(cmd, mode, attended, ws,
                                                   monkeypatch):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, cmd, mode, attended)
    assert verdict is None, verdict
    assert asked == []
    assert admitted == [True]


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("mode", ["plan", "bypassPermissions"])
@pytest.mark.parametrize("cmd", HARMLESS[:8])
def test_the_other_modes_do_not_use_the_cage_rule(cmd, mode, attended, ws,
                                                  monkeypatch):
    _cage(monkeypatch, "none")
    before = _gate(ws, cmd, mode, attended)
    _cage(monkeypatch, "bwrap")
    after = _gate(ws, cmd, mode, attended)
    assert after == before
    assert after[2] == []


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", UNCHANGED)
def test_a_guarded_command_is_handled_exactly_as_before(cmd, mode, attended,
                                                        ws, monkeypatch):
    _cage(monkeypatch, "none")
    before = _gate(ws, cmd, mode, attended)
    _cage(monkeypatch, "bwrap")
    after = _gate(ws, cmd, mode, attended)
    assert after == before
    assert after[2] == []


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("cmd", HARMLESS)
def test_control_without_isolation_nothing_changes(cmd, attended, ws,
                                                   monkeypatch):
    """Detection monkeypatched to a host where nothing can isolate: the
    real decision function answers 'none' and the old verdict stands."""
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: False)
    assert A._bash_isolation_decision(
        KitToolPermissions(mode="default", workspace=str(ws)))[0] == "none"
    verdict, asked, admitted = _gate(ws, cmd, "default", attended)
    assert admitted == []
    if attended:
        assert len(asked) == 1 and verdict is not None
    else:
        assert asked == [] and "not on the auto-allow list" in verdict


@pytest.mark.parametrize("kind", ["landlock", "seatbelt", "refuse"])
def test_only_bubblewrap_admits(kind, ws, monkeypatch):
    """Landlock cannot carve the protected files out of a writable root,
    Seatbelt is not exercised here: both keep the dialog."""
    _cage(monkeypatch, kind)
    verdict, asked, admitted = _gate(ws, "make", "default", False)
    assert admitted == [] and "not on the auto-allow list" in verdict


def test_an_open_network_keeps_the_dialog(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "open")
    verdict, _asked, admitted = _gate(ws, "make", "default", False)
    assert admitted == [] and "not on the auto-allow list" in verdict


def test_the_setting_switches_the_rule_off(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings", lambda *a, **k: {
        "agent": {"caged_commands_run_unasked": False}})
    verdict, _asked, admitted = _gate(ws, "make", "default", False)
    assert admitted == [] and "not on the auto-allow list" in verdict


def test_a_refusal_already_given_still_holds(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    perms.record_denied_action("bash", "make")
    admitted: list[bool] = []
    verdict = _doc_executor._run_permission_gate(
        "bash", {"command": "make"}, perms, caged_admission=admitted)
    assert admitted == [] and "already refused" in verdict


def test_an_mcp_shell_is_never_admitted(ws, monkeypatch):
    """An MCP server runs the command outside this process's cage."""
    _cage(monkeypatch, "bwrap")
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    verdict = _doc_executor._gate_mcp_tool(
        "mcp__srv__bash", {"command": "make"}, perms)
    assert verdict and "not on the auto-allow list" in verdict


def test_the_executor_refuses_an_admission_the_cage_no_longer_holds(
        ws, monkeypatch):
    """Belt and braces: admitted at the gate, the decision changed before
    the argv was built -- the command is not run."""
    _cage(monkeypatch, "none")
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    argv = A._bash_isolation_argv("make", ws, perms, caged_unasked=True)
    out = subprocess.run(argv, capture_output=True, text=True)
    assert out.returncode == 126 and "not run" in out.stderr


# -- through the real executor, where bubblewrap works ----------------------

def _run(ws: Path, cmd: str, mode: str = "default",
         attended: bool = False) -> dict:
    asked: list[str] = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    out = json.loads(_doc_executor.execute("bash", {"command": cmd}, perms))
    out["_asked"] = asked
    return out


@requires_a_working_sandbox
@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
@pytest.mark.parametrize("cmd", HARMLESS)
def test_a_harmless_command_runs_through_the_executor(cmd, mode, attended,
                                                      ws):
    """The command dialog is gone. The executor's own path gates still
    run after it: in default mode a named write target is asked about
    ('bash -> out.txt'), and a shell overwrite needs a prior read -- both
    as before, and both are not the auto-allow refusal."""
    out = _run(ws, cmd, mode, attended)
    err = out.get("error", "")
    assert "not on the auto-allow list" not in err
    assert not [p for p in out["_asked"] if p.startswith("$ ")]
    if mode == "acceptEdits" or not attended:
        if "sed -i" not in cmd:
            assert not err, out


@requires_a_working_sandbox
def test_what_the_caged_command_computes_comes_back(ws):
    out = _run(ws, 'python3 -c "import os; print(os.getcwd())"')
    assert str(ws.resolve()) in json.dumps(out)
    out = _run(ws, 'python3 -c "open(\'out.txt\',\'w\').write(\'x\')"')
    assert (ws / "out.txt").read_text() == "x"


@requires_a_working_sandbox
def test_the_caged_command_cannot_write_outside_the_workspace(ws, tmp_path):
    target = tmp_path / "outside.txt"
    _run(ws, f'python3 -c "open(\'{target}\',\'w\').write(\'x\')"')
    assert not target.exists()


@requires_a_working_sandbox
def test_the_caged_command_cannot_rewrite_protected_files(ws):
    """git metadata, the agent's settings and CI configuration are
    read-only for a command nobody approved."""
    subprocess.run(["git", "init", "-q", str(ws)], check=True,
                   capture_output=True)
    (ws / ".delfin").mkdir()
    (ws / ".delfin" / "settings.json").write_text("{}")
    (ws / ".github").mkdir()
    (ws / ".github" / "ci.yml").write_text("on: push\n")
    git_cfg = (ws / ".git" / "config").read_text()
    # The names are assembled at run time: a literal path is refused by
    # the gate long before the cage, which is not what is measured here.
    (ws / "w.sh").write_text(
        'a=.gi; b=.delfi; c=.githu\n'
        'echo x >> "${a}t/config"\n'
        'echo x >> "${b}n/settings.json"\n'
        'echo y >> "${c}b/ci.yml"\n'
        'echo done\n')
    out = _run(ws, "bash w.sh")
    assert "error" not in out, out
    assert "Read-only file system" in json.dumps(out)
    assert (ws / ".git" / "config").read_text() == git_cfg
    assert (ws / ".delfin" / "settings.json").read_text() == "{}"
    assert (ws / ".github" / "ci.yml").read_text() == "on: push\n"
    # An approved command (bypass profile) keeps its ordinary cage: the
    # overlay is for the unasked command only.
    (ws / "w2.sh").write_text('c=.githu\necho y >> "${c}b/ci.yml"\n')
    perms = KitToolPermissions(mode="bypassPermissions", workspace=str(ws))
    _doc_executor.execute("bash", {"command": "bash w2.sh"}, perms)
    assert (ws / ".github" / "ci.yml").read_text() == "on: push\ny\n"
