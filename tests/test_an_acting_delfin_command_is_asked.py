"""A ``delfin`` / ``delfin-agent`` command that acts is asked; one that reads is not.

Input: a bash command naming DELFIN's own command-line programs, sent
through the real permission gate (``_DocToolExecutor().execute("bash")``)
in a DELFIN workspace -- the only workspace where the DELFIN-specific
auto-allow entries are merged -- in ``default`` and ``acceptEdits``, with
a confirm callback (attended) and without one (headless).

Output under test: whether the command reaches the shell executor, and
whether the user was asked first.

* An acting subcommand (stops processes, writes or deletes credentials,
  moves bug reports, writes files outside the workspace, starts a model
  session, a calculation or a SLURM job, changes a schedule, answers an
  approval) is put to the user when someone can be asked, and refused
  with a message naming the command when nobody can.
* A read-only subcommand (lists, shows, status, help, capabilities, gate
  explain, doctor without --repair) runs without a question.
* A rule the user wrote (``bash_auto_allow_patterns`` from
  ``remember_permission``, or a ``permissions.allow`` entry in the user
  settings file) still allows an acting one.

The entry this replaces, ``^\\s*delfin(?:-\\w+)?\\b``, matched every
subcommand of every DELFIN program, and also any path starting with
``delfin`` (``delfin/scripts/x.sh``).

The shell itself is replaced by a recorder (``_execute_bash``): the
commands under test stop processes and submit jobs, and what is measured
is the gate's decision, which is taken before that call.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import delfin.agent.api_client as A

MODES = ("default", "acceptEdits")

ACTING = [
    # delfin-agent
    "delfin-agent",
    'delfin-agent "summarise the repo"',
    "delfin-agent --mode solo",
    'delfin-agent chat "hello"',
    'delfin-agent run "hello"',
    "delfin-agent init .",
    "delfin-agent init --force .",
    "delfin-agent bench",
    "delfin-agent bench run --repeats 3",
    "delfin-agent bench set-baseline runs/x.json --yes",
    "delfin-agent bench nightly",
    "delfin-agent credentials set KIT_API_KEY",
    "delfin-agent credentials delete KIT_API_KEY",
    "delfin-agent credentials migrate",
    "delfin-agent session export abc123",
    "delfin-agent watchpost run",
    "delfin-agent watchpost diff",
    "delfin-agent watchpost baseline",
    "delfin-agent watchpost enable --every 30",
    "delfin-agent watchpost disable",
    "delfin-agent memory tidy --apply",
    "delfin-agent memory tidy --app",
    "delfin-agent approvals approve r1 --code 1234",
    "delfin-agent approvals deny r1",
    "delfin-agent approvals answer r1 2",
    "delfin-agent skills accept foo",
    "delfin-agent skills reject foo",
    "delfin-agent bug triage",
    "delfin-agent bug solve 20260101-000000",
    "delfin-agent bug watch --once",
    "delfin-agent scheduler start",
    "delfin-agent scheduler stop",
    "delfin-agent scheduler add-bench --every 60",
    "delfin-agent stop-all",
    "delfin-agent stop-all --yes",
    "delfin-agent stop-all --check --yes",
    "delfin-agent pause worker",
    "delfin-agent resume worker",
    "delfin-agent doctor --repair",
    "delfin-agent doctor --rep",
    "delfin-agent where --attach",
    "delfin-agent where --ssh-config --write",
    "delfin-agent jobs watch --on",
    "delfin-agent jobs watch --off",
    "delfin-agent test-on-slurm repo main cpu 10 tests",
    # the chemistry CLI and its siblings
    "delfin",
    "delfin --define",
    "delfin --recalc",
    "delfin cleanup",
    "delfin stop",
    "delfin co2",
    "delfin run_orca job.inp",
    "delfin qm_run xtb -- --gfn 2 x.xyz",
    "delfin tadf_xtb",
    "delfin hyperpol",
    "delfin cluster submit",
    "delfin review",
    "delfin-pipeline run",
    "delfin-step xtb_opt",
    "delfin-build complex.smi",
    "delfin-guppy x.smi",
    "delfin-guppy-batch list.txt",
    "delfin-manta x.smi",
    "delfin-fukui x.xyz",
    "delfin-cluster bench",
    "delfin-voila",
    "delfin-app",
    "delfin-docs-index",
    "delfin-tools-server",
    "delfin-ops-server",
    "delfin-json",
    # a path that merely starts with "delfin"
    "delfin/scripts/run_all.sh",
    # the same programs by module name
    "python -m delfin",
    "python3 -m delfin --recalc",
    "python -m delfin.cli cleanup",
    "python3 -m delfin.agent.cli stop-all --yes",
    "python -m delfin.agent.cli bench run",
    "python3.11 -m delfin.cli_pipeline run",
]

READ_ONLY = [
    "delfin-agent --help",
    "delfin-agent -h",
    "delfin-agent --version",
    "delfin-agent stop-all --help",
    "delfin-agent bench set-baseline --help",
    "delfin-agent bench list",
    "delfin-agent bench latest --limit 5",
    "delfin-agent bench audit runs/x.json",
    "delfin-agent bench compare runs/a.json runs/b.json --json",
    "delfin-agent bench gate runs/x.json --json",
    "delfin-agent session ls --limit 3",
    "delfin-agent session search 'orca timeout'",
    "delfin-agent messages ls --from k1",
    "delfin-agent messages status m1 --from k1",
    "delfin-agent memory",
    "delfin-agent memory tidy",
    "delfin-agent approvals",
    "delfin-agent approvals ls",
    "delfin-agent approvals show r1",
    "delfin-agent approvals watch --seconds 2",
    "delfin-agent skills",
    "delfin-agent skills proposals",
    "delfin-agent skills show foo",
    "delfin-agent bug",
    "delfin-agent bug ls",
    "delfin-agent scheduler",
    "delfin-agent scheduler status",
    "delfin-agent stop-all --check",
    "delfin-agent doctor",
    "delfin-agent doctor --workspace .",
    "delfin-agent capabilities --json",
    "delfin-agent capabilities --permission-mode acceptEdits --json",
    "delfin-agent gate explain 'git push origin main'",
    "delfin-agent gate explain 'ls -la' --json",
    "delfin-agent report --json",
    "delfin-agent watch --denied --seconds 5",
    "delfin-agent where",
    "delfin-agent where --ssh-config",
    "delfin-agent sessions --limit 5",
    "delfin-agent sessions --load",
    "delfin-agent jobs",
    "delfin-agent jobs --json",
    "delfin-agent jobs watch --status",
    "delfin --help",
    "delfin -h",
    "delfin --version",
    "delfin -V",
    "delfin doctor",
    "delfin doctor --json",
    "delfin qm_check",
    "delfin qm_check xtb crest",
    "delfin csp_check",
    "delfin mlp_check",
    "delfin analysis_check",
    "delfin-agent bench list 2>&1",
    "python -m delfin --help",
    "python3 -m delfin.agent.cli capabilities --json",
    "python -m delfin.agent.cli stop-all --check",
]


@pytest.fixture(autouse=True)
def _no_cage_admission(monkeypatch):
    """The decision under test is the auto-allow tier's. On a host with a
    bubblewrap cage the final tier may run an unlisted command caged
    instead of asking (tests/test_a_caged_command_runs_without_asking.py);
    that admission is a separate rule and is switched off here."""
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    d = tmp_path / "ws"
    (d / "delfin" / "agent").mkdir(parents=True)
    (d / "delfin" / "__init__.py").write_text("")
    (d / "delfin" / "agent" / "__init__.py").write_text("")
    (d / "runs").mkdir()
    return d


@pytest.fixture
def ran(monkeypatch):
    calls: list[str] = []

    def _fake(self, arguments, perms, caged_unasked=False):
        calls.append(arguments.get("command", ""))
        return json.dumps({"stdout": "", "stderr": "", "exit_code": 0})

    monkeypatch.setattr(A._DocToolExecutor, "_execute_bash", _fake)
    return calls


_n = [0]


def _perms(ws: Path, mode: str, asked: list | None) -> A.KitToolPermissions:
    _n[0] += 1
    p = A.KitToolPermissions(mode=mode, workspace=str(ws))
    p.task_session_id = f"delfin-cli-{_n[0]}"
    if asked is not None:
        def _cb(name, args, preview):
            asked.append((name, (args or {}).get("command", "")))
            return False
        p.confirm_callback = _cb
    return p


def _bash(perms, cmd: str) -> str:
    return A._DocToolExecutor().execute("bash", {"command": cmd}, perms)


def test_the_workspace_merges_the_delfin_entries(ws):
    p = _perms(ws, "default", None)
    assert p.is_delfin_workspace
    assert set(A._DELFIN_BASH_AUTO_ALLOW) <= set(p.bash_auto_allow_patterns)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ACTING)
def test_an_acting_command_is_asked_when_someone_can_ask(ws, ran, mode, cmd):
    asked: list = []
    out = _bash(_perms(ws, mode, asked), cmd)
    assert ran == [], f"{cmd!r} ran unasked in {mode}: {out}"
    # Asked -- or refused outright by an earlier tier (the secret-path
    # deny for `credentials`, the approval guard for `approvals
    # approve|deny|answer`), which is stricter than asking.
    assert [c for (_n, c) in asked] == [cmd] or (
        asked == [] and '"error"' in out), (cmd, mode, asked, out)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ACTING)
def test_an_acting_command_is_refused_headless_with_its_name(
        ws, ran, mode, cmd):
    out = _bash(_perms(ws, mode, None), cmd)
    assert ran == [], f"{cmd!r} ran headless in {mode}: {out}"
    # The refusal names the command -- or the earlier rule that refused
    # it outright (the secret-path deny names the word it matched).
    assert cmd.split()[0] in out or "secret-deny path ('credentials')" in out, out


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("attended", [True, False],
                         ids=["attended", "headless"])
@pytest.mark.parametrize("cmd", READ_ONLY)
def test_a_read_only_command_runs_unasked(ws, ran, mode, attended, cmd):
    asked: list = []
    out = _bash(_perms(ws, mode, asked if attended else None), cmd)
    assert ran == [cmd], f"{cmd!r} did not run in {mode}: {out}"
    assert asked == []


@pytest.mark.parametrize("mode", MODES)
def test_a_users_pattern_still_allows_an_acting_one(ws, ran, mode):
    p = _perms(ws, mode, None)
    p.bash_auto_allow_patterns = p.bash_auto_allow_patterns + (
        r"^\s*delfin-agent\s+stop-all\s+--yes\s*$",)
    out = _bash(p, "delfin-agent stop-all --yes")
    assert ran == ["delfin-agent stop-all --yes"], out
    # The user's pattern is exactly that one command.
    out = _bash(p, "delfin-agent credentials delete KIT_API_KEY")
    assert ran == ["delfin-agent stop-all --yes"], out


@pytest.mark.parametrize("mode", MODES)
def test_a_users_allow_rule_still_allows_an_acting_one(ws, ran, mode):
    settings = Path.home() / ".delfin_settings.json"
    settings.write_text(json.dumps(
        {"permissions": {"allow": ["bash(delfin-agent bug solve:*)"]}}))
    asked: list = []
    out = _bash(_perms(ws, mode, asked), "delfin-agent bug solve 2026-x")
    assert ran == ["delfin-agent bug solve 2026-x"], out
    assert asked == []


@pytest.mark.parametrize("cmd", [
    "delfin-agent bench list; delfin-agent stop-all --yes",
    "delfin-agent bench list && delfin-agent credentials delete K",
    "delfin-agent bench list $(delfin-agent stop-all --yes)",
    "delfin-agent session search \"$(delfin-agent stop-all --yes)\"",
    "delfin-agent bench list > notes.txt",
    "delfin-agent gate explain `delfin-agent stop-all --yes`",
])
def test_nothing_rides_along_on_a_read_only_entry(ws, ran, cmd):
    out = _bash(_perms(ws, "default", None), cmd)
    assert ran == [], f"{cmd!r} ran: {out}"


@pytest.mark.parametrize("cmd", [
    "delfin-agent credentials", "delfin-agent credentials list",
    "delfin-agent credentials set K", "delfin-agent credentials delete K",
])
def test_every_credentials_form_is_refused_by_the_secret_path_deny(
        ws, ran, cmd):
    """Not an auto-allow decision: the word is on the secret deny-list,
    which refuses before any allow tier, reading forms included."""
    asked: list = []
    out = _bash(_perms(ws, "acceptEdits", asked), cmd)
    assert ran == [] and asked == [], out
    assert "secret-deny" in out, out


def _entry_modules() -> list[str]:
    """The modules behind [project.scripts], read from pyproject.toml."""
    root = Path(A.__file__).resolve().parents[2]
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    block = text.split("[project.scripts]", 1)[1].split("\n[", 1)[0]
    mods = []
    for line in block.splitlines():
        if "=" in line and ":" in line:
            mods.append(line.split("=", 1)[1].strip().strip('"').split(":")[0])
    assert len(mods) >= 15, mods
    return [m[: -len(".__main__")] if m.endswith(".__main__") else m
            for m in mods]


@pytest.mark.parametrize("mode", MODES)
def test_every_console_script_module_is_asked_by_module_name(
        ws, ran, mode):
    """`python -m <entry module>` is the program by another name; a new
    console script whose module the exclusion does not know would be
    free under the general `python -m <module>` entry."""
    free = []
    for mod in _entry_modules() + ["delfin"]:
        cmd = f"python3 -m {mod}"
        ran.clear()
        _bash(_perms(ws, mode, None), cmd)
        if ran:
            free.append(cmd)
    assert free == []


@pytest.mark.parametrize("cmd", [
    "python -m molkit run",
    "python3 -m delfin_tools_x",
    "python -m delfin.agent.bench_tasks",
])
def test_another_module_still_runs_unasked(ws, ran, cmd):
    """Control: the exclusion is the entry modules, not every module whose
    name starts with delfin."""
    out = _bash(_perms(ws, "default", None), cmd)
    assert ran == [cmd], out
