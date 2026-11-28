"""A solo turn on an office request gets the office tools and the office rules.

Measured 2026-11-27 on a field install where people do everyday office work:
a German request to give every collaborator in column E its own worksheet ran
in Code (solo) mode, which is where a person on that install stays. It failed
in three sessions -- empty turns, nine repeated reads, and finally a correct
script the shell could not run. The package's question was why neither of the
two intended routes (work with the full office tool set and rules, or delegate
to the office specialist) happened at all.

Driven through ``AgentEngine.stream_response`` -- the product's own turn
assembly, not a helper -- the base tree answers with:

* **0 of 12 document tools**, because the office tool GROUP was decided from a
  filesystem scan, and the vocabulary that lets a request ask for the group
  had no German worksheet word in it. The prompt module for the same subject
  fired correctly on its own table (``tabelle``, ``spalte``), so the system
  prompt said the module was active while every tool it names had been
  stripped. The empty folder is not an edge case here: it is what a "make me
  a workbook" request starts in.
* **none of the three checking rules** that make a spreadsheet number true.
  They were written only into ``office_agent.md``, which only the office role
  loads, and no mechanism let one rule be stated in two roles.

Both are RED on the base for those reasons and green after. The control below
is a coding request in a code repository: office tools and office rules are
both absent there, on either tree, and that is the point -- a test that only
showed the office request gaining things could not tell "the request decided"
from "everything is always on".

Nothing here depends on what the host has installed. The document-backend and
model-capability probes are stubbed to "available" (they are availability
questions, not behaviour) and the HTTP transport is the only thing faked. The
prompt, the tool group decision and the tool surface are the product's own.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as AC
from delfin.agent import engine as E
from delfin.agent import mcp_client as MCP
from delfin.agent import model_capabilities as MC
from delfin.agent.api_client import _OFFICE_TOOL_NAMES
from delfin.agent.prompt_loader import PromptLoader

#: The request the field install failed on, verbatim.
COLUMN_SPLIT = (
    "lege für jeden collaborator in Spalte E ein Tabellenblatt an, in dem "
    "alle zugehörigen Einträge zu dem Collaborator erfasst sind")

CODING = "fix the failing test in tests/test_login.py"

#: The three rules, by the sentence that carries each. Named as the rule
#: reads, so deleting one from the shared definition fails this test.
RULES = (
    "Read before you write",
    "Verify after you write",
    "Report coverage",
)


class _Delta:
    content, tool_calls, reasoning_content = "ok", None, None


class _Choice:
    index, delta, finish_reason = 0, _Delta(), "stop"


class _Chunk:
    choices, usage = [_Choice()], None


class _Stream:
    def __iter__(self):
        return iter([_Chunk()])

    def close(self):
        pass


def _stub_probes(monkeypatch):
    """Answer the availability questions a real host would answer."""
    class _Registry:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    caps = MC.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(AC._doc_executor, "_ensure_loaded", lambda: True)
    monkeypatch.setattr(AC._doc_executor, "_ensure_calc_loaded", lambda: True)
    monkeypatch.setattr(AC, "_office_backends_available", lambda: True)
    monkeypatch.setattr(AC, "_office_backend_set", lambda: None)
    monkeypatch.setattr(MC, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(MCP, "get_registry", lambda *a, **k: _Registry())


def _turn(tmp_path, request, workspace, monkeypatch):
    """One real turn: what the client SENDS for *request*.

    Returns (offered tool names, system prompt). The prompt and the tool
    surface read here are the ones the product assembles for that
    workspace; only ``chat.completions.create`` is replaced.
    """
    _stub_probes(monkeypatch)
    seen: dict = {}

    def _create(**kwargs):
        seen["tools"] = kwargs.get("tools") or []
        seen["messages"] = kwargs.get("messages") or []
        return _Stream()

    engine = E.AgentEngine(repo_dir=workspace, backend="api",
                           provider="ollama", model="ollama.qwen2.5-coder:32b",
                           mode="solo")
    engine.client.client.chat.completions.create = _create
    list(engine.stream_response(request))

    from conftest import offered_tool_names
    names = set(offered_tool_names(seen.get("tools") or []))
    prompt = str((seen.get("messages") or [{}])[0].get("content") or "")
    return names, prompt


@pytest.fixture
def empty_ws(tmp_path):
    """A plain folder: no document, no chemistry file, not a DELFIN tree.

    This is the state the request starts in -- one has to be able to ask
    for a workbook that does not exist yet.
    """
    ws = tmp_path / "notes"
    ws.mkdir()
    (ws / "notes.md").write_text("# notes\n")
    return ws


@pytest.fixture
def code_ws(tmp_path):
    """A code repository, so the control turn has something to be about."""
    ws = tmp_path / "app"
    (ws / "tests").mkdir(parents=True)
    (ws / "tests" / "test_login.py").write_text("def test_a():\n    assert 1\n")
    (ws / "pyproject.toml").write_text('[project]\nname = "app"\n')
    return ws


# --- the two red-on-base rows -------------------------------------------

def test_an_office_request_in_an_empty_folder_gets_the_document_tools(
        tmp_path, empty_ws, monkeypatch):
    """The request names a workbook; the tools that make one must be there.

    Fails on the base: the folder holds no document, so the group was never
    offered and a solo session saw none of the 12. The input that would make
    this fail for the WRONG reason is a request naming no document -- see the
    control below.
    """
    names, _ = _turn(tmp_path, COLUMN_SPLIT, empty_ws, monkeypatch)
    missing = sorted(_OFFICE_TOOL_NAMES - names)
    assert not missing, (
        "a German worksheet request in an empty folder was offered "
        f"{len(_OFFICE_TOOL_NAMES) - len(missing)} of "
        f"{len(_OFFICE_TOOL_NAMES)} document tools; missing: {missing}")


def test_an_office_request_gets_the_three_checking_rules(
        tmp_path, empty_ws, monkeypatch):
    """Advice for the work arrives WITH the tools for the work."""
    _, prompt = _turn(tmp_path, COLUMN_SPLIT, empty_ws, monkeypatch)
    absent = [r for r in RULES if r not in prompt]
    assert not absent, (
        "a solo prompt answering an office request is missing the rules that "
        f"make its numbers true: {absent}")


# --- the control: coding work pays for none of it ----------------------

def test_a_coding_request_gets_neither_the_tools_nor_the_rules(
        tmp_path, code_ws, monkeypatch):
    """Same code path, different subject. It must not have grown.

    Green on both trees: this is what makes the two rows above a finding
    rather than a surface that is simply always on.
    """
    names, prompt = _turn(tmp_path, CODING, code_ws, monkeypatch)
    assert not (names & _OFFICE_TOOL_NAMES), sorted(names & _OFFICE_TOOL_NAMES)
    present = [r for r in RULES if r in prompt]
    assert not present, (
        f"a coding turn was told the office checking rules: {present}")


# --- the rules have ONE definition, which is why they cannot drift -----

def test_the_two_roles_read_the_same_rule_text():
    """A rule stated twice is a rule that drifts.

    The office role's prompt and the solo module both render the same block
    from ``pack/shared/``. This asserts the block is what each of them
    carries, so deleting the shared file, or pasting a second copy back into
    one role, is caught.
    """
    loader = PromptLoader()
    rendered = loader._fill_generated_blocks(
        "<!-- generated:shared:office_work_rules -->\nstale copy\n"
        "<!-- /generated -->\n")
    for rule in RULES:
        assert rule in rendered, rule
    # A name nobody owns is left alone -- the mechanism's existing contract.
    unknown = ("before\n<!-- generated:shared:not_a_file -->\nkept\n"
               "<!-- /generated -->\nafter\n")
    assert loader._fill_generated_blocks(unknown) == unknown
