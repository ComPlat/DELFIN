"""Office memory: conventions are remembered, records and other domains
are not.

Two separate guarantees are pinned here, and they fail for different
reasons:

1. **Conventions, never records.** Administrative folders hold names,
   addresses, salaries and case numbers. A store that quietly accumulates
   those out of the documents it works on is worse than no store — it
   persists past the session and outlives the folder. The check runs on
   the way IN and rejects; it does not redact, because a half-scrubbed
   record is still a record.
2. **Domains do not mix.** Office mode already has the chemistry tools
   removed where tools are authorised. Memory was the remaining channel
   through which that subject matter could still reach an office turn, so
   the domain is a property of the stored file and is checked before the
   file is read into a prompt — not an instruction in the role prompt.

``Path.home`` is monkey-patched throughout so no test touches the real
memory store.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as ms
from delfin.agent.prompt_loader import PromptLoader


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    # The registry of office folders is cached per process; each test gets
    # a fresh one so registrations cannot leak between them.
    monkeypatch.setattr(ms, "_office_ws_cache", None, raising=False)
    return tmp_path


@pytest.fixture
def office_dir(fake_home, tmp_path):
    """A registered office folder — writes for it resolve to the office
    domain exactly as they do after an office turn has run."""
    d = tmp_path / "office"
    d.mkdir()
    ms.register_office_workspace(d)
    return d


@pytest.fixture
def code_dir(fake_home, tmp_path):
    d = tmp_path / "checkout"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# The guardrail: what may be written
# ---------------------------------------------------------------------------

CONVENTIONS = [
    "the Kostenstelle is in column D, not column C",
    "dates in this folder are day-first (DD.MM.YYYY) and amounts use a "
    "decimal comma",
    "Mahnung.docx is the template for the Offene-Posten list; its key "
    "column is Belegnr",
    "series files are named <Belegnr>_<Nachname>.pdf, lowercase",
    "the Empfaenger field of the form maps to the Name column of the list",
]

RECORDS = [
    ("an amount", "the open amount for the Meier invoice is 1.234,50"),
    ("a date", "the reporting deadline was 31.07.2026"),
    ("a record key", "invoice R-014 is still unpaid"),
    ("a record number", "the case number is 4711234"),
    ("an e-mail address", "send the list to buchhaltung@example.org"),
    ("an IBAN", "the account is DE02 1203 0000 0000 2020 51"),
    ("a table row", "Belegnr;Name;Betrag|R-1;Meier|R-2;Schmidt"),
]


@pytest.mark.parametrize("text", CONVENTIONS)
def test_a_convention_is_storable(text):
    assert ms.check_office_convention(text) is None, text


@pytest.mark.parametrize("label, text", RECORDS)
def test_a_value_out_of_a_document_is_refused(label, text):
    reason = ms.check_office_convention(text)
    assert reason is not None, text
    # The message has to say what was found, or the agent cannot rewrite
    # the memory into a storable form.
    assert label.split()[-1] in reason or "table row" in reason, reason


def test_an_excerpt_is_refused_by_length():
    reason = ms.check_office_convention("a " * 400)
    assert reason is not None
    assert "sentence" in reason


def test_a_pasted_block_is_refused_by_line_count():
    reason = ms.check_office_convention(
        "\n".join(f"line {chr(97 + i)} of the note" for i in range(9)))
    assert reason is not None
    assert "excerpt" in reason


def test_a_format_may_be_named_even_though_a_value_may_not():
    """The line the check draws: DD.MM.YYYY is a convention, 31.07.2026 is
    a record. Naming the format is how the agent stores the same knowledge
    without carrying the document's content with it."""
    assert ms.check_office_convention("dates are DD.MM.YYYY") is None
    assert ms.check_office_convention("dates are 31.07.2026") is not None


def test_saving_a_record_into_an_office_folder_is_rejected_and_writes_nothing(
        office_dir):
    with pytest.raises(ms.OfficeMemoryRejected):
        ms.save_typed_memory(
            "project: the open amount for Meier is 1.234,50",
            repo_root=office_dir)
    assert ms.list_typed_memories(office_dir) == []
    assert ms.list_typed_memories(office_dir, scope="user") == []


def test_a_rejected_title_does_not_slip_a_value_in(office_dir):
    """The title is stored in frontmatter and in the index hook, so it is
    a second way into the store and gets the same check."""
    with pytest.raises(ms.OfficeMemoryRejected):
        ms.save_typed_memory(
            "project: the key column of the list",
            repo_root=office_dir, title="invoice R-014")


def test_a_convention_is_stored_with_the_office_domain(office_dir):
    fpath, _slug, _mtype = ms.save_typed_memory(
        "project: the Kostenstelle is in column D, not column C",
        repo_root=office_dir)
    assert fpath.exists()
    assert "domain: office" in fpath.read_text(encoding="utf-8")
    rec = ms.list_typed_memories(office_dir)[0]
    assert rec["domain"] == ms.DOMAIN_OFFICE


def test_an_office_memory_cannot_be_filed_in_the_user_wide_store(office_dir):
    """The user-wide store is read from every workspace. An office session
    must not be able to reach it, including through the text prefix that
    normally routes a fact there."""
    fpath, _slug, _mtype = ms.save_typed_memory(
        "global: feedback: series files are named <Belegnr>_<Nachname>.pdf",
        repo_root=office_dir,
        allow_scope_prefix=True)
    assert ms.list_typed_memories(office_dir, scope="user") == []
    assert fpath.parent == ms._delfin_memory_dir(office_dir)


def test_the_guardrail_does_not_apply_outside_office_folders(code_dir):
    """Technical memories legitimately cite line numbers and versions; the
    record check would reject most of them and is deliberately scoped."""
    fpath, _slug, _mtype = ms.save_typed_memory(
        "project: the retry limit is set in delfin/agent/engine.py:1234",
        repo_root=code_dir)
    assert fpath.exists()
    assert ms.list_typed_memories(code_dir)[0]["domain"] == ms.DOMAIN_CODE


# ---------------------------------------------------------------------------
# Domain separation
# ---------------------------------------------------------------------------

def test_domain_for_role_maps_the_office_role_and_mode():
    assert ms.domain_for_role("office_agent", "office") == ms.DOMAIN_OFFICE
    assert ms.domain_for_role("solo_agent", "office") == ms.DOMAIN_OFFICE
    assert ms.domain_for_role("office_agent", "") == ms.DOMAIN_OFFICE
    assert ms.domain_for_role("solo_agent", "solo") == ms.DOMAIN_CODE
    assert ms.domain_for_role("", "") == ms.DOMAIN_CODE


def test_domain_visibility_matrix():
    assert ms.domain_visible(ms.DOMAIN_OFFICE, ms.DOMAIN_OFFICE)
    assert ms.domain_visible(ms.DOMAIN_GENERAL, ms.DOMAIN_OFFICE)
    assert ms.domain_visible(ms.DOMAIN_GENERAL, ms.DOMAIN_CODE)
    assert not ms.domain_visible(ms.DOMAIN_CODE, ms.DOMAIN_OFFICE)
    assert not ms.domain_visible(ms.DOMAIN_OFFICE, ms.DOMAIN_CODE)
    # No domain on the turn = a caller that does not track domains; the
    # store behaves as it did before the split.
    assert ms.domain_visible(ms.DOMAIN_CODE, "")


def test_a_memory_written_before_domains_existed_counts_as_technical():
    """Reading the missing field as neutral would make every pre-existing
    memory visible to office turns — the leak the split closes."""
    assert ms.memory_text_domain("---\nname: x\n---\n\nsome old fact\n") == (
        ms.DOMAIN_CODE)
    assert not ms.domain_visible(
        ms.memory_text_domain("---\nname: x\n---\n\nold\n"),
        ms.DOMAIN_OFFICE)


def test_filter_memory_index_drops_out_of_domain_pointers(fake_home, tmp_path):
    """The index hook alone already names a memory's subject, so filtering
    only the bodies would still put it in the prompt."""
    mem = tmp_path / "store"
    mem.mkdir()
    (mem / "project_orca.md").write_text(
        "---\nname: orca\ndomain: code\n---\n\nORCA runs use def2-TZVP\n",
        encoding="utf-8")
    (mem / "project_list.md").write_text(
        "---\nname: list\ndomain: office\n---\n\nkey column is Belegnr\n",
        encoding="utf-8")
    index = (
        "# DELFIN Project Memory\n\n"
        "## Project (current work)\n"
        "- [orca](project_orca.md) — ORCA basis set default\n"
        "- [list](project_list.md) — the list's key column\n"
    )
    out = ms.filter_memory_index(index, mem, ms.DOMAIN_OFFICE)
    assert "project_list.md" in out
    assert "project_orca.md" not in out
    assert "ORCA" not in out
    # Nothing visible at all means the store contributes nothing.
    assert ms.filter_memory_index(index, mem, "research") == ""


# ---------------------------------------------------------------------------
# Recall: what actually reaches the model
# ---------------------------------------------------------------------------

def _office_prompt(loader: PromptLoader, workspace: Path, task: str) -> str:
    loader.workspace_root = workspace
    return loader.build_system_prompt(
        role_id="office_agent", mode_id="office", route=["office_agent"],
        task_text=task, session_key="office-memory-test")


def test_a_chemistry_memory_does_not_reach_an_office_prompt(
        office_dir, code_dir, monkeypatch):
    """The failure this prevents: the agent is filling forms and starts
    talking about ORCA because a memory from another mode was recalled."""
    loader = PromptLoader()
    # Written the way a technical session writes it: project store of the
    # checkout, and user-wide so it follows the user everywhere.
    ms.save_typed_memory(
        "project: ORCA jobs in this repo default to the def2-TZVP basis set",
        repo_root=code_dir)
    ms.save_typed_memory(
        "global: feedback: always ground ORCA keywords in the manual before "
        "answering",
        repo_root=code_dir,
        allow_scope_prefix=True)
    # The office folder has one convention of its own.
    ms.save_typed_memory(
        "project: the Kostenstelle is in column D, not column C",
        repo_root=office_dir)

    prompt = _office_prompt(loader, office_dir, "fill the ORCA form again")

    # Asserted on the memories' own wording: the shared prompt addenda name
    # chemistry in their own right, and that is not the channel under test.
    assert "def2-TZVP" not in prompt
    assert "ground ORCA keywords" not in prompt
    assert "ORCA jobs in this repo" not in prompt
    assert "column D" in prompt


def test_the_folders_own_conventions_do_reach_the_office_prompt(office_dir):
    loader = PromptLoader()
    ms.save_typed_memory(
        "project: Mahnung.docx is the template for the Offene-Posten list",
        repo_root=office_dir)
    prompt = _office_prompt(loader, office_dir, "prepare the reminders")
    assert "Mahnung.docx" in prompt
    assert "Folder Conventions" in prompt


def test_a_neutral_memory_reaches_both_domains(office_dir, code_dir):
    """``general`` is the explicit opt-in for a fact that is genuinely
    about neither kind of work — stored user-wide, so it is the one thing
    an office turn does pick up from outside its own folder."""
    ms.save_typed_memory(
        "global: user: answers stay short and end with the open points",
        repo_root=code_dir, domain=ms.DOMAIN_GENERAL,
        allow_scope_prefix=True)
    loader = PromptLoader()
    office_prompt = _office_prompt(loader, office_dir, "check the list")
    assert "end with the open points" in office_prompt


def test_an_office_convention_does_not_reach_a_technical_prompt(
        office_dir, code_dir):
    ms.save_typed_memory(
        "project: series files are named <Belegnr>_<Nachname>.pdf",
        repo_root=office_dir)
    loader = PromptLoader()
    loader.workspace_root = office_dir
    out = loader._load_external_memory_context(
        task_text="how are the series files named",
        domain=ms.DOMAIN_CODE)
    assert "Belegnr" not in out


def test_another_folders_session_does_not_see_the_office_store(
        office_dir, code_dir, tmp_path):
    """Scope, independently of domain: the conventions belong to the folder
    they were learned in."""
    ms.save_typed_memory(
        "project: the key column of the recurring list is Belegnr",
        repo_root=office_dir)
    other = tmp_path / "other-office"
    other.mkdir()
    ms.register_office_workspace(other)
    assert ms.list_typed_memories(other) == []
    loader = PromptLoader()
    loader.workspace_root = other
    out = loader._load_external_memory_context(
        task_text="what is the key column", domain=ms.DOMAIN_OFFICE)
    assert "Belegnr" not in out


def test_a_legacy_untagged_memory_is_not_recalled_into_an_office_turn(
        office_dir, monkeypatch):
    """Files written before the domain field existed are the realistic
    population of the user-wide store."""
    global_dir = ms._delfin_global_memory_dir()
    global_dir.mkdir(parents=True, exist_ok=True)
    (global_dir / "user_legacy.md").write_text(
        "---\nname: legacy\ndescription: a standing preference\n---\n\n"
        "prefers the xtb pre-optimisation before every ORCA run\n",
        encoding="utf-8")
    (global_dir / "MEMORY.md").write_text(
        "# DELFIN Global Memory (all projects)\n\n"
        "## User (background)\n"
        "- [legacy](user_legacy.md) — a standing preference\n",
        encoding="utf-8")
    loader = PromptLoader()
    loader.workspace_root = office_dir
    out = loader._load_external_memory_context(
        task_text="prepare the run list", domain=ms.DOMAIN_OFFICE)
    assert "xtb" not in out
    assert "legacy" not in out


def test_a_hand_written_index_bullet_carries_no_domain(fake_home, tmp_path):
    """A fact typed straight into MEMORY.md has no file and no domain
    field. It is judged like any other memory that predates the split:
    technical work sees it, an office turn does not."""
    mem = tmp_path / "store"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "# DELFIN Project Memory\n\n- the cluster queue is called uc3\n",
        encoding="utf-8")
    loader = PromptLoader()
    kept = loader._load_external_memory_context(
        memory_root=mem, domain=ms.DOMAIN_CODE)
    assert "uc3" in kept
    dropped = loader._load_external_memory_context(
        memory_root=mem, domain=ms.DOMAIN_OFFICE)
    assert dropped == ""


def test_the_engine_registers_the_office_folder_for_the_write_path(
        fake_home, tmp_path):
    """Without this the store cannot know that a write for this folder is
    an office write, and neither the domain tag nor the record check would
    apply."""
    from types import SimpleNamespace

    from delfin.agent.engine import AgentEngine

    folder = tmp_path / "buero"
    folder.mkdir()
    engine = SimpleNamespace(
        mode="office", repo_dir=folder,
        kit_permissions=SimpleNamespace(workspace=folder))

    assert not ms.is_office_workspace(folder)
    domain = AgentEngine._scope_memory_domain(engine, "office_agent")
    assert domain == ms.DOMAIN_OFFICE
    assert ms.is_office_workspace(folder)
    # A file inside the folder counts too — the whole tree is office work.
    assert ms.is_office_workspace(folder / "2026" / "invoices")

    other = tmp_path / "src"
    other.mkdir()
    engine_code = SimpleNamespace(
        mode="quick", repo_dir=other,
        kit_permissions=SimpleNamespace(workspace=other))
    assert AgentEngine._scope_memory_domain(
        engine_code, "builder_agent") == ms.DOMAIN_CODE
    assert not ms.is_office_workspace(other)


def test_the_office_folder_registry_survives_a_new_process(office_dir,
                                                           monkeypatch):
    """The domain of a write is derived from the folder, so a later
    process must resolve the same folder the same way."""
    monkeypatch.setattr(ms, "_office_ws_cache", None, raising=False)
    assert ms.is_office_workspace(office_dir)


def test_the_role_prompt_says_when_to_write_and_what_is_refused():
    """Recall is worth nothing if nothing is ever written, so the moment
    to write has to be named. The refusal is enforced in the store; the
    prompt only explains it, so the agent can rephrase instead of retrying
    the same rejected memory."""
    prompt = PromptLoader().load_role_prompt("office_agent")
    assert "`remember`" in prompt
    assert "conventions only" in prompt
    assert "DD.MM.YYYY" in prompt


def test_recall_without_a_domain_is_unchanged(code_dir):
    """Entry points that do not pass a domain keep the previous behaviour,
    so nothing that already worked starts returning less."""
    ms.save_typed_memory(
        "project: the retry limit lives in the engine", repo_root=code_dir)
    loader = PromptLoader()
    monkey_root = ms._delfin_memory_dir(code_dir)
    out = loader._load_external_memory_context(memory_root=monkey_root)
    assert "retry limit" in out
