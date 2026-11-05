"""Every act of the guided MLIP set-up is asked; its reads are not.

Input: the commands the plan proposes (environment create, torch and
orca-mlips installs for each backend and each conda-family tool, the
environment removal, the verify job's sbatch), sent through the real
permission gate (``_run_permission_gate``) in default, acceptEdits and
bypassPermissions -- the last with a working cage and the action watcher
down, the case in which ``_caged_still_asked`` decides.

Output under test: each command either reached the user's dialog or was
refused; none ran unasked. The read-only MCP tools are on the gate's
read-only list and registered in the server's read-only block; the two
acting ones are neither, and are host-gated. The set-up text is not in
the always-on prompt or in the native tool schemas.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin import api
from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor
from delfin.mlp_tools import guided_setup as gs

MODES = ("default", "acceptEdits", "bypassPermissions")


def _probe(tool: str) -> dict:
    return {
        "slurm": {"available": True, "gpu_partitions": [
            {"name": "gpu_a100", "gpus": 8, "gpu_type": "a100", "gres": [],
             "cuda_capable": True}]},
        "cuda": {"driver_cuda": "12.4", "gpus": []},
        "orca": {"ready": True, "installations": [], "modules_6_1_or_newer": ["orca/6.1.1"]},
        "conda": {"tool": tool, "path": f"/opt/{tool}/bin/{tool}", "on_path": True,
                  "envs": ["/opt/envs/base"]},
        "mlip_envs": [], "settings": {"effective": {}},
    }


def _acts() -> list:
    out = []
    for tool in ("micromamba", "mamba", "conda"):
        for backend in gs.ORCA_MLIPS_BACKENDS:
            plan = gs.propose(_probe(tool), backend=backend)
            out += [c["command"] for c in plan["commands"] if c["kind"] == "run"]
        out.append(plan["remove_command"]["command"])
    # The same tool by absolute path (not on PATH).
    p = _probe("micromamba")
    p["conda"]["on_path"] = False
    out += [c["command"] for c in gs.propose(p, backend="orb")["commands"]]
    out.append("sbatch /scratch/job/mlp_verify/mlp_verify.sbatch")
    return sorted(set(out))


ACTS = _acts()


class _WatcherDown:
    def __call__(self, prompt, model):
        raise RuntimeError("endpoint error")


@pytest.fixture
def session(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": "m", "rules": [], "broad": False})

    def make(mode):
        ws = tmp_path / "wt"
        ws.mkdir(exist_ok=True)
        perms = KitToolPermissions(workspace=str(ws), mode=mode)
        perms.classifier_runner = _WatcherDown()
        perms.classifier_messages = [{"role": "user", "content": "set up MLIPs"}]
        shown = []
        perms.confirm_callback = lambda name, args, preview: shown.append(
            (args or {}).get("command")) or False
        return perms, shown
    return make


def test_the_plan_yields_the_commands_under_test():
    assert any("python=3.11" in c for c in ACTS)
    assert any("download.pytorch.org/whl/cu124" in c for c in ACTS)
    for b in gs.ORCA_MLIPS_BACKENDS:
        assert any(f"orca-mlips[{b}]" in c for c in ACTS), b
    assert any(" env remove " in c for c in ACTS)
    assert not any(c.startswith("hf ") for c in ACTS)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ACTS)
def test_every_setup_act_is_asked(session, mode, cmd):
    perms, shown = session(mode)
    out = _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)
    # Refused (the dialog answered no, or an earlier rule refused) -- never
    # allowed without the user having been shown the command.
    assert out is not None, (mode, cmd)
    assert shown == [cmd] or "rejected" in out or "refused" in out.lower(), (mode, cmd, out)


def test_the_caged_installer_rule_covers_the_install_steps():
    for cmd in ACTS:
        if cmd.startswith("sbatch"):
            assert A._is_heavy_compute_command(cmd)
        else:
            assert A._caged_still_asked(cmd), cmd


# ------------------------------------------------------------ tool surface


READ_ONLY = ("mlp_setup_probe", "mlp_setup_plan", "mlp_setup_verify_result")
ACTING = ("mlp_setup_configure", "mlp_setup_verify")


def _server_text() -> str:
    return (Path(api.__file__).resolve().parent / "ops_server" / "server.py").read_text()


def test_read_only_tools_are_on_the_gate_list_and_in_the_read_only_block():
    text = _server_text()
    start = text.index("# Read-only — register module functions directly")
    end = text.index("# P1 — statistical plots", start)
    block = re.findall(r'mcp\.tool\(name="([a-z_]+)"\)', text[start:end])
    for name in READ_ONLY:
        assert name in A._MCP_READONLY_TOOL_BASES, name
        assert name in block, name


def test_acting_tools_are_host_gated_and_not_read_only():
    text = _server_text()
    for name in ACTING:
        assert name not in A._MCP_READONLY_TOOL_BASES, name
        assert re.search(rf'mcp\.tool\(name="{name}"\)\(_host_gated\(', text), name
    for name in ("install", "submit", "remove"):
        assert f'name="mlp_setup_{name}"' not in text


def test_every_setup_tool_is_catalogued():
    names = {e["name"] for e in api._TOOL_CATALOG}
    assert set(READ_ONLY + ACTING) <= names


def test_the_acting_tools_refuse_without_the_host_grant(tmp_path, monkeypatch):
    from delfin.ops_server import server

    monkeypatch.setattr(api, "mlp_setup_configure",
                        lambda **k: pytest.fail("reached the settings write"))
    monkeypatch.setattr(api, "mlp_setup_verify",
                        lambda **k: pytest.fail("reached the job writer"))
    out = json.loads(server.tool_mlp_setup_configure(write=True, allow_mutate=False))
    assert out["written"] is False and "grant" in out["error"]
    out = json.loads(server.tool_mlp_setup_verify(job_dir=str(tmp_path), allow_mutate=False))
    assert "grant" in out["error"]


def test_the_setup_text_is_not_in_the_always_on_surface(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    from delfin.agent.prompt_loader import PromptLoader

    prompt = PromptLoader().build_system_prompt(
        "solo_agent", "solo", task_text="rename a function in a parser module")
    assert "mlp_setup" not in prompt
    assert "mlp_setup" not in json.dumps(A._DOC_TOOLS_OPENAI)
    from delfin.agent import skills
    listing = A._skill_listing(skills.discover_skills(None, domain="code"))
    assert "mlp_setup" not in listing
