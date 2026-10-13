"""The ``method-reviewer`` preset: the level of theory, checked before a number is compared.

Measured gap (benchmark task
``science_energies_from_different_methods_are_not_ranked``,
kit.deepseek-v4-flash, 0/5): the session opened every ``.out`` and no
``.inp``, separated xtb from ORCA, never read a functional, and ranked
B3LYP runs below PBE0 runs as a chemical result. A one-sentence prompt rule
did not close it (5/20 without vs 1/20 with) and made answers longer. This
preset is the delegated alternative: a read-only checker whose whole job is
the method, reached through the tools DELFIN already has.

Asserted here:

* it is discovered, selectable by name, and a pack file rather than code;
* its tool list is read-only, and the refusal holds at the executor, not
  only in the advertised schema;
* the existing readers it is pointed at (``parse_orca_output``,
  ``compare_calculations``, the scf/imaginary-frequency extractors) are
  admitted for the child by the session gate and the plan-mode MCP gate;
* those readers see the mismatch in the two fixtures, so no new tool is
  needed for the functional and the basis;
* driven through ``run_subagent`` with the real ``OpenAIClient`` tool loop
  and a scripted endpoint, the preset prompt reaches the child, a
  ``write_file`` call is refused and writes nothing, and ``read_file``
  returns both keyword lines.

Fixtures: ``tests/fixtures/method_reviewer/`` -- two synthetic ORCA
outputs of the same water geometry, B3LYP/def2-SVP vs PBE0-D3BJ/def2-TZVP.
"""

from __future__ import annotations

import json
import shutil
import types
from pathlib import Path

import pytest

from delfin.agent import api_client as ac
from delfin.agent import subagents as SA

_NAME = "method-reviewer"
_PACK = Path(SA.__file__).resolve().parent / "pack" / "agents"
_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "method_reviewer"

_WRITERS = ("write_file", "edit_file", "multi_edit", "apply_patch",
            "notebook_edit", "bash", "run_tests", "subagent", "orchestrate",
            "remember", "publish_report")

# The existing DELFIN readers the preset is pointed at. Server tools are
# listed with their server prefix: ``unknown_tool_names`` cannot check an
# MCP name at startup, and the bare name is not a native tool.
_OPS_READERS = (
    "mcp__delfin-ops__parse_orca_output",
    "mcp__delfin-ops__compare_calculations",
    "mcp__delfin-ops__extract_scf_convergence",
    "mcp__delfin-ops__extract_imaginary_frequencies",
)


@pytest.fixture(autouse=True)
def _isolated(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    SA.reload_subagent_presets()
    yield
    monkeypatch.undo()
    SA.reload_subagent_presets()


def _preset():
    return SA.SUBAGENT_PRESETS[_NAME]


def _workspace(tmp_path: Path) -> Path:
    ws = tmp_path / "ws"
    shutil.copytree(_FIXTURES, ws)
    return ws


def _child(ws: Path, mode: str = "default"):
    parent = ac.KitToolPermissions(workspace=ws, mode=mode)
    preset = _preset()
    return SA._derive_perms(parent, preset.mode, preset_tools=preset.tools)


# ---------------------------------------------------------------------------
# Selectable
# ---------------------------------------------------------------------------

def test_the_preset_is_discovered_and_selectable():
    assert _NAME in SA.SUBAGENT_PRESETS
    assert _NAME in SA.subagent_type_names()


def test_the_subagent_tool_enum_offers_it():
    enum = None
    for tool in ac._DOC_TOOLS_OPENAI:
        fn = tool.get("function", {})
        if fn.get("name") == "subagent":
            enum = fn["parameters"]["properties"]["subagent_type"]["enum"]
    assert enum is not None and _NAME in enum


def test_it_is_a_pack_file_not_code():
    assert (_PACK / f"{_NAME}_subagent.md").is_file()
    assert _NAME not in SA._BUILTIN_PRESETS


def test_the_role_prompt_names_it():
    text = (_PACK / "solo_agent.md").read_text(encoding="utf-8")
    assert f"`{_NAME}`" in text


# ---------------------------------------------------------------------------
# Read-only
# ---------------------------------------------------------------------------

def test_it_runs_in_plan_mode_with_a_narrow_list():
    preset = _preset()
    assert preset.mode == "plan"
    assert preset.tools and len(preset.tools) < 60
    for writer in _WRITERS:
        assert writer not in preset.tools, writer
    assert "subagent_message" in preset.tools
    assert ac.unknown_tool_names(preset.tools, ac._DOC_TOOLS_OPENAI) == []


def test_a_write_is_refused_by_the_preset_list_even_without_plan_mode(tmp_path):
    """bypassPermissions switches every mode refusal off; the list still binds."""
    ws = _workspace(tmp_path)
    parent = ac.KitToolPermissions(workspace=ws, mode="bypassPermissions")
    child = SA._derive_perms(parent, "bypassPermissions",
                             preset_tools=_preset().tools)
    ex = ac._DocToolExecutor()
    out = ex.execute("write_file",
                     {"path": "verdict.txt", "content": "x"}, child)
    assert "--allowed-tools" in json.loads(out).get("error", "")
    assert not (ws / "verdict.txt").exists()


@pytest.mark.parametrize("name", _OPS_READERS)
def test_the_existing_readers_are_admitted_for_the_child(name, tmp_path):
    preset = _preset()
    assert name in preset.tools
    child = _child(_workspace(tmp_path))
    assert ac._tool_denied_for_session(child, name) is None
    assert ac._bare_tool_name(name) in ac._PLAN_READONLY_MCP_TOOLS


# ---------------------------------------------------------------------------
# The method, in the prompt
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("term", [
    "functional", "basis", "dispersion", "solva", "charge", "multiplicity",
    "ecp", "grid", "converged", "imaginary", "spin contamination",
    "not found", "method difference", "file:line", "severity",
])
def test_the_checklist_is_in_the_prompt(term):
    assert term in _preset().system_prompt.lower()


def test_the_prompt_names_the_control_txt_trap():
    """CONTROL.txt ``method`` is the DELFIN workflow (classic/OCCUPIER),
    not the functional (``calc_indexer._extract_from_control_txt``)."""
    body = _preset().system_prompt
    assert "CONTROL.txt" in body and "OCCUPIER" in body


def test_the_description_says_when_to_use_it():
    desc = _preset().description
    assert 40 < len(desc) < 200
    assert "before" in desc.lower()


# ---------------------------------------------------------------------------
# The existing readers see the mismatch in the fixtures
# ---------------------------------------------------------------------------

def test_parse_orca_output_reads_both_levels_of_theory():
    from delfin.api import parse_orca_output
    a = parse_orca_output(str(_FIXTURES / "water_b3lyp" / "water.out"))
    b = parse_orca_output(str(_FIXTURES / "water_pbe0" / "water.out"))
    assert (a.functional, a.basis) == ("B3LYP", "def2-SVP")
    assert (b.functional, b.basis) == ("PBE0", "def2-TZVP")
    assert a.scf_converged and b.scf_converged
    assert a.opt_converged and b.opt_converged


def test_compare_calculations_reports_the_mismatch():
    from delfin.api import compare_calculations
    diff = compare_calculations(str(_FIXTURES / "water_b3lyp"),
                                str(_FIXTURES / "water_pbe0"))
    assert diff.method_match is False
    assert diff.basis_match is False


# ---------------------------------------------------------------------------
# Driven: run_subagent -> real OpenAIClient tool loop -> scripted endpoint
# ---------------------------------------------------------------------------

_FINAL = (
    "water_b3lyp/water.out:11 B3LYP/def2-SVP, no dispersion. "
    "water_pbe0/water.out:11 PBE0-D3BJ/def2-TZVP. "
    "HIGH: functional, basis and dispersion differ -- a method difference, "
    "not chemistry."
)


def _drive(tmp_path):
    ws = _workspace(tmp_path)
    seen: list[dict] = []
    script = [
        [("write_file", {"path": str(ws / "verdict.txt"), "content": "x"})],
        [("read_file", {"path": str(ws / "water_b3lyp" / "water.out")}),
         ("read_file", {"path": str(ws / "water_pbe0" / "water.out")})],
    ]

    class _Stream:
        def __init__(self, i):
            self.i = i

        def __iter__(self):
            if self.i < len(script):
                calls = [types.SimpleNamespace(
                    index=k, id=f"c{self.i}_{k}",
                    function=types.SimpleNamespace(
                        name=name, arguments=json.dumps(args)))
                    for k, (name, args) in enumerate(script[self.i])]
                delta = types.SimpleNamespace(
                    content="", tool_calls=calls, reasoning_content=None)
                fin = "tool_calls"
            else:
                delta = types.SimpleNamespace(
                    content=_FINAL, tool_calls=None, reasoning_content=None)
                fin = "stop"
            return iter((types.SimpleNamespace(usage=None, choices=[
                types.SimpleNamespace(delta=delta, finish_reason=fin)]),))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    seen.append(json.loads(json.dumps(
                        {"messages": kw.get("messages"),
                         "tools": kw.get("tools")}, default=str)))
                    return _Stream(len(seen) - 1)

    perms = ac.KitToolPermissions(workspace=ws, mode="default",
                                  task_session_id="S")
    client = ac.OpenAIClient(
        api_key="x", model="kit.glm-5.3", provider="kit",
        base_url="https://ki-toolbox.scc.kit.edu/api/v1", permissions=perms)
    client.client = _Stub()
    result = SA.run_subagent(
        subagent_type=_NAME,
        description="check the method of two water runs",
        prompt="Compare water_b3lyp/water.out and water_pbe0/water.out.",
        parent_client=client, parent_perms=perms, model="parent",
        inherit_context=False, max_tool_calls=10, max_wall_s=60)
    return result, seen, ws


def _system_of(request: dict) -> str:
    return "\n".join(str(m.get("content") or "") for m in request["messages"]
                     if m.get("role") == "system")


def _tool_results(request: dict) -> dict[str, str]:
    return {m.get("tool_call_id"): str(m.get("content") or "")
            for m in request["messages"] if m.get("role") == "tool"}


def test_the_preset_prompt_reaches_the_child(tmp_path):
    result, seen, _ = _drive(tmp_path)
    assert seen, "the child never sent a request"
    system = _system_of(seen[0])
    assert "method-reviewer" in system
    assert "spin contamination" in system.lower()
    assert result.error == "", result.error
    assert result.final_text.strip() == _FINAL


def test_the_child_is_not_offered_a_writer(tmp_path):
    _, seen, _ = _drive(tmp_path)
    offered = {t["function"]["name"] for t in (seen[0]["tools"] or [])}
    assert "read_file" in offered
    for writer in ("write_file", "edit_file", "bash", "multi_edit"):
        assert writer not in offered, writer


def test_a_write_the_child_calls_anyway_is_refused_and_writes_nothing(tmp_path):
    _, seen, ws = _drive(tmp_path)
    assert len(seen) >= 2
    refusal = _tool_results(seen[1]).get("c0_0", "")
    assert refusal, "no tool result for the write call"
    assert "error" in refusal.lower() or "refused" in refusal.lower()
    assert not (ws / "verdict.txt").exists()


def test_the_child_reads_both_keyword_lines(tmp_path):
    _, seen, _ = _drive(tmp_path)
    assert len(seen) >= 3
    results = _tool_results(seen[2])
    assert "B3LYP def2-SVP" in results.get("c1_0", "")
    assert "PBE0 D3BJ def2-TZVP" in results.get("c1_1", "")


# ---------------------------------------------------------------------------
# The names the prompt and the list point at exist where they are claimed
# ---------------------------------------------------------------------------

def test_every_listed_server_tool_is_a_registered_read_only_ops_tool(tmp_path):
    """Each ``mcp__delfin-ops__*`` entry is a name the ops server registers,
    and the plan-mode MCP gate admits it for the child. The four-name
    subset above left three entries (compare_across_functionals,
    extract_delfin_json, find_orca_errors) unchecked."""
    import re
    from delfin.ops_server import server as _srv_mod_path  # noqa: F401
    src = Path(_srv_mod_path.__file__).read_text(encoding="utf-8")
    registered = set(re.findall(r'mcp\.tool\(name="([a-z_]+)"\)', src))
    listed = [t for t in _preset().tools if t.startswith("mcp__")]
    assert listed
    child = _child(_workspace(tmp_path))
    ex = ac._DocToolExecutor()
    for name in listed:
        assert name.startswith("mcp__delfin-ops__"), name
        assert ac._bare_tool_name(name) in registered, name
        assert ex._gate_mcp_tool(name, {}, child) is None, name


def test_every_control_key_the_prompt_names_is_in_the_shipped_template():
    """The prompt and ``delfin/define.py`` TEMPLATE are two copies of one
    vocabulary; a key the template does not have sends the child to a
    line no real CONTROL.txt carries."""
    import re
    from delfin.define import TEMPLATE
    template_keys = {ln.split("=", 1)[0].strip()
                     for ln in TEMPLATE.splitlines() if "=" in ln}
    body = _preset().system_prompt
    block = body[body.index("CONTROL.txt ("):]
    block = block[:block.index(")")]
    named = set(re.findall(r"`([A-Za-z_0-9]+)`", block))
    for key in ("functional", "disp_corr", "main_basisset", "metal_basisset",
                "aux_jk", "relativity", "ri_jkx", "implicit_solvation_model",
                "solvent", "multiplicity_0"):
        assert key in named, key
    assert named <= template_keys, named - template_keys


@pytest.mark.parametrize("term", [
    # composite methods fix their own basis and dispersion
    "gfn2", "-3c",
    # dispersion built into the functional
    "wb97x-v",
    # def2 carries an ECP beyond Kr without a keyword
    "def2-ecp",
    # the program version decides the defaults
    "version",
    # what kind of energy, and on which geometry
    "gibbs", "//",
    # open-shell reference
    "uks",
    # CONTROL.txt is the request, the .inp what ran
    "_rel",
])
def test_the_pitfalls_that_mislabel_a_mismatch_are_named(term):
    assert term in _preset().system_prompt.lower(), term
