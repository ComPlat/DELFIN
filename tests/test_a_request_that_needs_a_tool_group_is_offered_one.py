"""A request that needs a tool group is OFFERED that group.

The defect this closes was measured, not imagined. ``toolsets.py`` decided the
two optional groups (``chemistry``, ``office``) from the workspace's FILES, the
role and a setting -- never from what the person asked for. Driven through the
product's own turn assembly (``AgentEngine.stream_response`` ->
``PromptLoader.compose_sections`` -> the real api client's advertised tool
list), the plain request

    "make me a two-page docx summary of this text: ..."

in a workspace that happens to hold no document got **65 tools, none of them
able to write a .docx**, while the system prompt told the model

    Not loaded here: office tools (...). If the task needs them, ask the user
    for /tools enable office

so its best move on an ordinary request was to ask for a switch instead of
doing the work. Full 11-row matrix: ``scratch/offer_matrix.md``.

RED ON THE BASELINE, for the reason it names: ``test_a_document_request_...``
and ``test_a_chemistry_request_...`` fail on origin/main because the group is
absent from the request the client sends. The three cases that must hold on
BOTH trees -- a code task pays for neither group, a role that could not
execute a tool is not offered it, an explicit ``off`` outranks a request --
are green before and after, which is what makes the two red rows a finding
rather than a broken test.

Nothing here depends on what the host has installed: the document-backend and
model-capability probes are stubbed to "available", and the HTTP transport is
the only thing faked. The prompt, the tool surface and the group decision are
the product's own code.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as AC
from delfin.agent import engine as E
from delfin.agent import mcp_client as MCP
from delfin.agent import model_capabilities as MC
from delfin.agent.api_client import _OFFICE_TOOL_NAMES
from delfin.agent.toolsets import (
    CHEMISTRY,
    OFFICE,
    WorkspaceToolsets,
    request_groups,
)

_CHEM_MCP = "mcp__delfin-ops__parse_orca_output"

# A third-party MCP server: not a DELFIN group, so it must never be dropped
# by a group decision. Its presence also proves the request was really built.
_MCP_NAMES = (_CHEM_MCP, "mcp__delfin-ops__extract_energy_table",
              "mcp__delfin-tools__catalog", "mcp__userserver__lookup")


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
    """Answer the probes a real endpoint and MCP server would answer.

    These are availability QUESTIONS, not behaviour: stubbing them keeps the
    test off the host (no ORCA install, no network, no model) without changing
    any decision the product makes about which groups to offer.
    """
    class _Registry:
        def discover_all(self):
            out = []
            for full in _MCP_NAMES:
                _, server, name = full.split("__")
                out.append(MCP.MCPTool(server=server, name=name,
                                       description=f"{name} tool",
                                       schema={"type": "object"}))
            return out

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
    monkeypatch.setattr(AC.time, "sleep", lambda *a, **k: None)


def _turn(tmp_path, request: str, *, mode="solo", workspace=None,
          settings=None, monkeypatch=None):
    """One real turn; returns what the request OFFERS plus the prompt.

    The engine builds its own api client from ``repo_dir``, so the tool list
    read here is the one the product assembles for that workspace; only
    ``chat.completions.create`` is replaced.
    """
    _stub_probes(monkeypatch)
    ws = workspace or tmp_path
    seen: dict = {}

    def _create(**kwargs):
        seen["tools"] = kwargs.get("tools") or []
        seen["messages"] = kwargs.get("messages") or []
        return _Stream()

    if settings is not None:
        import delfin.user_settings as US

        real = US.load_settings

        def _patched(path=None):
            cfg = real(path) or {}
            agent = dict(cfg.get("agent") or {})
            agent["toolsets"] = settings
            cfg["agent"] = agent
            return cfg

        monkeypatch.setattr(US, "load_settings", _patched)

    engine = E.AgentEngine(repo_dir=ws, backend="api", provider="ollama",
                           model="ollama.qwen2.5-coder:32b", mode=mode)
    engine.client.client.chat.completions.create = _create
    list(engine.stream_response(request))

    from conftest import offered_tool_names
    names = set(offered_tool_names(seen.get("tools") or []))
    prompt = str((seen.get("messages") or [{}])[0].get("content") or "")
    return names, prompt


def _has_chemistry(names: set) -> bool:
    return (_CHEM_MCP in names or "mcp__delfin-tools__catalog" in names
            or "show_molecule" in names)


def _has_office(names: set) -> bool:
    return bool(names & _OFFICE_TOOL_NAMES)


# --- workspaces --------------------------------------------------------

@pytest.fixture
def empty_ws(tmp_path):
    """A plain folder: no document, no chemistry file, not a DELFIN tree."""
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


# --- the two red-on-base rows -----------------------------------------

def test_a_document_request_is_offered_the_office_group(
        tmp_path, empty_ws, monkeypatch):
    """The request names a .docx; the tools that write one must be there.

    Fails on origin/main: the folder holds no document, so no document tool is
    offered and the prompt tells the model to ask for ``/tools enable office``.
    The input that would make this fail for the WRONG reason is a request that
    names no document -- see the control below.
    """
    names, prompt = _turn(
        tmp_path,
        "make me a two-page docx summary of this text: The refinery cut its "
        "emissions by 12 percent last year, mostly by retiring two boilers.",
        workspace=empty_ws, monkeypatch=monkeypatch)
    assert _has_office(names), (
        "a document request in a folder with no document file got no "
        "document tool; the model can only refuse or ask for a switch")
    # The advice for the group arrives WITH the group, never instead of it.
    assert "/tools enable office" not in prompt


def test_a_chemistry_request_is_offered_the_chemistry_group(
        tmp_path, empty_ws, monkeypatch):
    """The request names an ORCA output; the ORCA parsers must be there.

    Fails on origin/main for the same reason as the office row: no ``.out`` in
    the folder means no parser in the request.
    """
    names, _ = _turn(tmp_path,
                     "parse the ORCA output and give me the Gibbs free energy "
                     "and the HOMO-LUMO gap",
                     workspace=empty_ws, monkeypatch=monkeypatch)
    assert _has_chemistry(names), (
        "an ORCA question in a folder with no .out got no parser")


# --- green on both trees: the fix must not just switch everything on ----

def test_a_code_request_pays_for_neither_group(tmp_path, code_ws, monkeypatch):
    """The control. A task that needs neither group is charged neither.

    Green before and after the change -- and it fails the moment the request
    vocabulary fires on ordinary programming words, which is the false
    positive ``toolsets`` already refuses for ``energy`` / ``scan`` /
    ``calc(100% - 20px)``.
    """
    names, prompt = _turn(tmp_path, "fix the failing pytest in tests/"
                                    "test_login.py",
                          workspace=code_ws, monkeypatch=monkeypatch)
    assert not _has_chemistry(names)
    assert not _has_office(names)
    assert "mcp__userserver__lookup" in names, "no request was built at all"
    # And the model is still TOLD what is missing, so the escape hatch stays.
    assert "/tools enable" in prompt


def test_a_chemistry_request_does_not_also_pay_for_documents(
        tmp_path, empty_ws, monkeypatch):
    """Offering is per group, not per turn: one group's words buy one group."""
    names, _ = _turn(tmp_path,
                     "compare the DFT single-point energies across the calc "
                     "folders",
                     workspace=empty_ws, monkeypatch=monkeypatch)
    assert _has_chemistry(names)
    assert not _has_office(names)


def test_a_role_that_could_not_execute_a_tool_is_not_offered_it(
        tmp_path, empty_ws, monkeypatch):
    """Offering stays a subset of what may EXECUTE (rule 5).

    ``dashboard_agent`` carries an execution allow-list with no document tool,
    so offering the office group to it would buy a call that comes straight
    back refused -- the "Unknown tool" failure mode. Green on both trees:
    ``advertisable_tools`` derives the surface from the executor's own gates.
    """
    names, _ = _turn(tmp_path,
                     "merge invoice_1.pdf and invoice_2.pdf into one PDF",
                     mode="dashboard", workspace=empty_ws,
                     monkeypatch=monkeypatch)
    assert not _has_office(names), (
        "the office group was offered to a role that would be refused every "
        "call to it")


def test_an_explicit_off_outranks_the_request(tmp_path, empty_ws, monkeypatch):
    """The person's own switch still wins, whatever they then ask for.

    ``agent.toolsets.office: off`` is somebody having already weighed this
    cost; a word in a request must not overturn it.
    """
    names, _ = _turn(tmp_path,
                     "Add a column of totals to quarterly_report.xlsx",
                     workspace=empty_ws, settings={"office": "off"},
                     monkeypatch=monkeypatch)
    assert not _has_office(names)


# --- the decision itself, at the object the prompt and the tools share --

def test_the_reason_is_said_in_the_group_listing():
    """"/tools groups" can explain WHY a group appeared.

    Fails if ``note_request`` offers a group without a reason: the listing
    prints the reason, and a bare True cannot be shown to a person.
    """
    ts = WorkspaceToolsets(None, settings={"chemistry": "auto",
                                           "office": "auto"})
    assert ts.note_request("open the spreadsheet and add a column") == [OFFICE]
    assert OFFICE in ts.decided()
    assert "request" in ts.decided()[OFFICE]
    # Idempotent: the same words on the next turn add nothing new.
    assert ts.note_request("open the spreadsheet and add a column") == []
    # And a group the request does not name stays hidden.
    assert CHEMISTRY not in ts.decided()


def test_a_group_offered_to_a_solo_session_can_actually_be_called(tmp_path):
    """Offered means callable, not merely named (rule 5).

    Two changes once added tool names to a surface, passed their own tests,
    and every call answered "Unknown tool". So the office row above is not
    enough on its own: this drives one of the offered document tools through
    the executor with a SOLO session's permissions -- the role the request
    arrived in -- and judges the file it wrote.

    ``draft_email`` is the pure-stdlib member of the group (``email.message``),
    so unlike the .docx/.xlsx/.pdf tools it needs no optional package and this
    cannot turn red or green on what the host happens to have installed.
    """
    import json

    from delfin.agent import api_client as A

    ws = tmp_path / "post"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws, agent_role="solo_agent")
    perms.mode = "acceptEdits"
    perms.confirm_callback = lambda n, a, prev: True
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(executor.execute(
        "draft_email",
        {"path": "summary.eml", "to": "billing@example.invalid",
         "subject": "Two-page summary", "body": "The refinery cut..."},
        perms))
    assert out.get("to") == ["billing@example.invalid"], out
    assert (ws / "summary.eml").exists(), out


def test_a_request_that_names_neither_group_offers_nothing():
    """The empty answer, asserted directly so the two rows above are trusted.

    Without this, an implementation that offered both groups to every string
    would pass the red-on-base rows and fail only here.
    """
    ts = WorkspaceToolsets(None, settings={"chemistry": "auto",
                                           "office": "auto"})
    assert ts.note_request("rename the variable and rerun the tests") == []
    assert ts.decided() == {}
    assert ts.note_request("") == []
    assert ts.note_request(None) == []


def test_a_later_turn_cannot_move_the_head_a_prefix_cache_has_served():
    """The request is read for the FIRST turn's decision, and only there.

    This is the constraint that makes the feature safe, and it is asserted
    against the implementation that ships rather than argued: an earlier
    reading folded the request in at every turn, so a second turn naming a
    new group rewrote ``role_prompt`` -- a stable row the provider had
    already cached (measured: 7-12 s warm against 199-266 s cold). The
    engine's held prompt keeps the SENT bytes identical, but the composed
    rows are what ``_prompt_head_check`` reports, and the module's own
    "Semantics" section says the decision is taken once, at the first turn.

    Fails against that earlier reading: the second call there returns the
    two groups it found. Passes here: the decision is already taken, so the
    later request adds nothing and the rows cannot move.
    """
    ts = WorkspaceToolsets(None, settings={"chemistry": "auto",
                                           "office": "auto"})
    # Turn 1: a code request, so no group is offered and the session's
    # decision is taken (this is the turn whose rows the cache keeps).
    assert ts.note_request("fix the failing unit test in src/") == []
    assert ts.decided() == {}
    # Turn 2 would plainly need both groups. The head is composed by now.
    assert ts.note_request(
        "the energy of the molecule looks wrong; also check the column in "
        "the spreadsheet") == []
    assert ts.decided() == {}, (
        "a later turn moved the session's tool-group decision, which rewrites "
        "stable prompt rows the prefix cache has already served")
    # The first turn of a NEW session is a different matter, and is offered.
    fresh = WorkspaceToolsets(None, settings={"chemistry": "auto",
                                              "office": "auto"})
    assert fresh.note_request(
        "the energy of the molecule looks wrong; also check the column in "
        "the spreadsheet") == [CHEMISTRY, OFFICE]


def test_ods_is_a_format_and_not_a_substring_of_an_ordinary_word():
    """A false positive here costs a group for the whole session.

    ``request_groups`` matches bare substrings, and "ods" is inside
    "methods" and "periods", so the plainest code request was offered the
    office group -- the trap the vocabulary's own comment names by refusing
    "excel" for "excellent" and "homo" for "homogeneous".
    """
    ts = WorkspaceToolsets(None, settings={"chemistry": "auto",
                                           "office": "auto"})
    for text in ("please refactor these methods, they are duplicated",
                 "show me the periods of that table"):
        assert ts.note_request(text) == [], (
            f"{text!r} was read as asking for a document format")
        assert ts.decided() == {}, ts.decided()
    # The format is still a term when it IS the word.
    for text in ("convert it to ods", "save as .ods please", "ods file"):
        assert request_groups(text) == {OFFICE: "ods"}, text
