"""FS1 adversarial: the secret-redaction wrappers must fail CLOSED at the
real call paths, and the repo guard must keep them that way.

Reviewed against the builder's fail_closed change (memory_store,
task_state, session_index) and its disk-side fix (TaskState._as_dict).
Every case here drove the production entry point (save_memory,
TaskState.save/render, state_commands.show, session_index/_export._scrub),
not the inner helper: a wrapper that works when called directly can still
be bypassed by the writer that consumes it.

Red cases found by these probes live with the builder (see the review
message); what is committed here is the green half: the behaviour the
package must keep.
"""
import logging

import pytest

from delfin.agent import (
    memory_store,
    session_export,
    session_index,
    state_commands,
    task_state,
)

# Built by concatenation, as the repo's secret guard requires.
GH = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6Q7r8"
# No vendor prefix, no KEY= anchor: only the exact-value pass can remove it.
# This is the shape a held credential really has, and the reason a redactor
# that fails open is worse than none.
SHAPELESS = "Zq7Wm3Rt9Yb2Nv5Kd8Hs4Lx6Pj1Gc0Fa"
PLACEHOLDER = "[redaction failed -- text withheld]"


@pytest.fixture
def redactor_fails(monkeypatch):
    """The house redactor raises.

    The once-per-process log flag is reset only when the module has it: on a
    tree without the fix the attribute is absent, and the fixture must not
    turn that into a setup error -- the case has to fail on its assertion,
    so the control shows the DEFECT and not a missing symbol.
    """
    def boom(*args, **kwargs):
        raise RuntimeError("redactor exploded")

    monkeypatch.setattr("delfin.agent.output_guard.scrub_secrets", boom)
    if hasattr(memory_store, "_REDACT_LOGGED_ONCE"):
        monkeypatch.setattr(memory_store, "_REDACT_LOGGED_ONCE", False)


# -- the memory store: what actually lands on disk --------------------------

def test_memory_store_disk_holds_the_placeholder_not_the_secret(
        tmp_path, redactor_fails):
    store = tmp_path / "agent_memory.json"
    memory_store.save_memory(f"note {GH}", path=store)
    on_disk = store.read_text(encoding="utf-8")
    assert [GH in on_disk, PLACEHOLDER in on_disk] == [False, True]


def test_memory_store_normal_path_redacts_the_token_and_keeps_the_prose(
        tmp_path):
    store = tmp_path / "agent_memory.json"
    memory_store.save_memory(f"Fe complex note {GH}", path=store)
    on_disk = store.read_text(encoding="utf-8")
    assert GH not in on_disk
    assert "Fe complex note" in on_disk      # a redactor eats secrets, not text


def test_a_shapeless_held_credential_is_removed_by_the_wrapper(monkeypatch):
    """Control: the exact-value pass reaches a credential no pattern sees,
    so the withheld-on-failure cases below blame the failure path."""
    from delfin.agent import web_tools

    monkeypatch.setattr(web_tools, "_known_secret_values", lambda: [SHAPELESS])
    out = memory_store._without_secrets(f"note with {SHAPELESS} inside")
    assert SHAPELESS not in out


# -- the task state: the JSON on disk and the rendered recap ----------------

def test_task_state_save_writes_the_placeholder_and_no_secret(
        tmp_path, redactor_fails):
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    st.commit(task="FS1", commit=f"used {GH} to push")
    st.save()
    on_disk = st.path.read_text(encoding="utf-8")
    assert [GH in on_disk, PLACEHOLDER in on_disk] == [False, True]


def test_task_state_save_redacts_with_a_healthy_redactor(tmp_path):
    """The disk file is a consumer too: a wrapper that scrubs only the
    rendered recap leaves the raw text in the JSON beside it."""
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    st.save()
    on_disk = st.path.read_text(encoding="utf-8")
    assert GH not in on_disk


def test_task_state_render_withholds_the_block_when_redaction_fails(
        tmp_path, redactor_fails):
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    rendered = st.render()
    assert GH not in rendered
    assert PLACEHOLDER in rendered


def test_show_and_render_do_not_raise_when_the_redactor_fails(
        tmp_path, redactor_fails):
    """render() is documented never to raise, and its caller (repl's fresh
    turn context) turns an exception into an empty recap: a redactor
    failure must not silently cost the session its state."""
    st = task_state.open(tmp_path / "task_state.json")
    st.set_person_instruction(f"export the key {GH}")
    body = state_commands.show(st)
    assert isinstance(body, str) and body
    assert GH not in body


# -- the index and the export ----------------------------------------------

def test_index_and_export_wrappers_withhold_on_a_redactor_failure(
        redactor_fails):
    text = f"index text with {GH}"
    assert GH not in session_index._scrub(text)
    assert GH not in session_export._scrub(text)
    # Equality, not just absence: the package promises ONE placeholder
    # constant. `session_index._scrub` spells the literal out by hand
    # (session_index.py:321) instead of importing REDACT_FAILED_PLACEHOLDER,
    # so a later edit to one string would silently split the vocabulary while
    # this absence-only check stayed green.
    assert session_index._scrub(text) == memory_store.REDACT_FAILED_PLACEHOLDER
    assert session_export._scrub(text) == memory_store.REDACT_FAILED_PLACEHOLDER


# -- the failure is logged once, with the exception type only --------------

def test_the_failure_is_logged_once_without_the_text(redactor_fails, caplog):
    caplog.set_level(logging.WARNING, logger="delfin.agent.memory_store")
    memory_store._without_secrets(f"a {GH}")
    memory_store._without_secrets(f"b {GH}")
    records = [r for r in caplog.records if r.name.startswith("delfin")]
    assert len(records) == 1
    assert "RuntimeError" in records[0].getMessage()
    assert GH not in records[0].getMessage()


# -- fail closed must not become fail forever ------------------------------

def test_every_wrapper_redacts_again_after_the_redactor_recovers(monkeypatch):
    """Withholding is for the duration of the failure, not for the process.

    A wrapper that latches its failure state is a different defect with the
    same root: after one transient redactor error every later write is
    withheld, so the memory store and the task state silently stop accepting
    text. The operator found exactly this in service/redact.py (the withhold
    flag was never cleared). These wrappers must recover on their own.
    """
    from delfin.agent import output_guard

    healthy = output_guard.scrub_secrets

    def boom(*args, **kwargs):
        raise RuntimeError("redactor exploded")

    # 1. while the redactor fails: withheld, no secret
    monkeypatch.setattr(output_guard, "scrub_secrets", boom)
    failed = [memory_store._without_secrets(f"note {GH}"),
              task_state._scrub(f"note {GH}"),
              session_index._scrub(f"note {GH}"),
              session_export._scrub(f"note {GH}")]
    assert [GH in t for t in failed] == [False, False, False, False]
    assert failed == [PLACEHOLDER] * 4

    # 2. the redactor works again: redaction resumes, plain text flows
    monkeypatch.setattr(output_guard, "scrub_secrets", healthy)
    recovered = [memory_store._without_secrets(f"note {GH}"),
                 task_state._scrub(f"note {GH}"),
                 session_index._scrub(f"note {GH}"),
                 session_export._scrub(f"note {GH}")]
    assert [GH in t for t in recovered] == [False, False, False, False]
    assert [PLACEHOLDER in t for t in recovered] == [False, False, False, False]


def test_a_wrapper_does_not_withhold_text_that_holds_no_secret(monkeypatch):
    """The recovery control, stated the other way round: a plain note must
    survive a healthy redactor untouched, so 'withhold' is not the normal
    outcome and the case above cannot pass by withholding everything."""
    out = memory_store._without_secrets("the complex was optimised at BP86")
    assert out == "the complex was optimised at BP86"


# -- the guard keeps the next wrapper honest -------------------------------

def test_the_guard_detects_the_real_known_followups():
    """The AST guard's exclusion list must name defects it really finds;
    an exclusion invented to make the scan quiet would hide them."""
    import importlib.util
    from pathlib import Path

    repo = Path(__file__).resolve().parent.parent
    path = repo / "tests" / "test_fs1_no_fail_open_redactors.py"
    assert path.is_file(), f"the guard file is missing: {path.name}"
    spec = importlib.util.spec_from_file_location("fs1_guard", path)
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)

    hits = set(guard.fail_open_redactors())
    found = [("api_client.py", "_redact_tool_result") in hits,
             ("working_state.py", "_scrub") in hits]
    assert found == [True, True], (
        f"the guard does not find the excluded defects: {found}")


def test_no_followup_entry_is_stale(tmp_path, monkeypatch):
    """A KNOWN_FOLLOWUPS entry that is no longer a real defect is a hole:
    the guard would keep tolerating that function, so a later regression in
    the file that was just fixed stays green. Every entry must still be
    found by the scan.

    Also checks the sharp end: an entry whose file has been FIXED fails
    here, which is the reminder to delete it in the same commit.
    """
    guard = _guard()
    hits = set(guard.fail_open_redactors())
    real = [entry for entry in guard.KNOWN_FOLLOWUPS if entry in hits]
    stale = sorted(set(guard.KNOWN_FOLLOWUPS) - set(real))
    assert not stale, (
        "these KNOWN_FOLLOWUPS entries no longer describe a fail-open "
        f"function -- delete them in the fix commit: {stale}")


def test_the_placeholder_is_one_constant_not_a_retold_literal():
    assert getattr(memory_store, "REDACT_FAILED_PLACEHOLDER", "") == PLACEHOLDER
    assert task_state._scrub("plain text") == "plain text"
    assert session_index._scrub("plain text") == "plain text"


# -- the guard must catch every ordinary spelling, not one shape -----------
# Driving fail_open_redactors() over a synthetic tree, as the guard itself
# runs in CI. The first guard walked module-level functions only and
# matched a bare `return text`, so `return str(text)`, an alias and a
# wrapper written as a class method all went unseen (review findings 4-5).

FAIL_OPEN_SPELLINGS = {
    "direct": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return text\n"),
    "str()": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return str(text)\n"),
    "strip()": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return text.strip()\n"),
    "f-string": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return f'{text}'\n"),
    "alias": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        out = text\n"
        "        return out\n"),
    "second positional argument": (
        "def f(a, text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return text\n"),
    "swallowed then returned": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        pass\n"
        "    return text\n"),
}

SAFE_WRAPPERS = {
    "placeholder": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return '[redaction failed -- text withheld]'\n"),
    "empty": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        return ''\n"),
    "detection only": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text) != text\n"
        "    except Exception:\n"
        "        return False\n"),
    "re-raises": (
        "def f(text):\n"
        "    try:\n"
        "        return scrub_secrets(text)\n"
        "    except Exception:\n"
        "        raise\n"),
    "no redactor": (
        "def f(text):\n"
        "    try:\n"
        "        return other(text)\n"
        "    except Exception:\n"
        "        return text\n"),
}

SUBPACKAGE_FAIL_OPEN = (
    "class R:\n"
    "    def text(self, t):\n"
    "        try:\n"
    "            return scrub_secrets(t)\n"
    "        except Exception:\n"
    "            return t\n")


def _guard():
    import importlib.util
    from pathlib import Path

    repo = Path(__file__).resolve().parent.parent
    path = repo / "tests" / "test_fs1_no_fail_open_redactors.py"
    assert path.is_file(), f"the guard file is missing: {path.name}"
    spec = importlib.util.spec_from_file_location("fs1_guard", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _scan(tmp_path, monkeypatch, sources, sub_source=None):
    guard = _guard()
    root = tmp_path / "agent"
    (root / "service").mkdir(parents=True)
    for i, src in enumerate(sources):
        (root / f"mod{i}.py").write_text(src, encoding="utf-8")
    if sub_source is not None:
        (root / "service" / "smod.py").write_text(sub_source, encoding="utf-8")
    monkeypatch.setattr(guard, "AGENT", root)
    return set(guard.fail_open_redactors())


def test_the_guard_finds_every_fail_open_spelling(tmp_path, monkeypatch):
    """One file per spelling, one assertion per spelling: a scan that finds
    ANY of them must not be able to stand in for the ones it misses."""
    sources = list(FAIL_OPEN_SPELLINGS.values())
    hits = _scan(tmp_path, monkeypatch, sources)
    missed = [name for i, name in enumerate(FAIL_OPEN_SPELLINGS)
              if (f"mod{i}.py", "f") not in hits]
    assert not missed, f"the guard misses these spellings: {missed}"


def test_the_guard_flags_no_safe_wrapper(tmp_path, monkeypatch):
    """Nothing new is refused: a wrapper that withholds, re-raises or only
    tests for a secret must stay green."""
    hits = _scan(tmp_path, monkeypatch, list(SAFE_WRAPPERS.values()))
    assert not hits, f"false positives: {sorted(hits)}"


def test_the_guard_scans_a_subpackage_and_its_methods(tmp_path, monkeypatch):
    """delfin/agent/service/ holds a fail-open redactor, and a wrapper there
    is written as a class method: a glob('*.py') scan of the package dir
    with a module-level walk saw neither."""
    hits = _scan(tmp_path, monkeypatch, [], SUBPACKAGE_FAIL_OPEN)
    assert ("smod.py", "text") in hits
