"""Diagnosis and the API key are reachable from the Settings tab.

User rule, 2026-11-17: important settings and repairs go through the
Settings tab as well, not only through a terminal. On a shared machine a
fix that exists only as a command line is a fix most people never make.

`delfin-agent doctor` already answers "what is wrong and what fixes it"
and every row carries the machine-actionable form where there is one, so
the tab renders those rows rather than inventing its own checks.

Two things this file pins because getting them wrong would mislead rather
than merely annoy:

* the key field is cleared BEFORE the value is stored, so a failure in the
  store cannot leave the browser holding it, and no code path echoes it;
* the status line names the SOURCE of each credential, and warns when the
  selected one comes from an environment variable -- the lookup order is
  env > keyring > file, so a value stored from the tab stays unused while
  such a variable is set, and a tab that said "stored" without that would
  be telling the person something untrue.
"""

from __future__ import annotations

import inspect

from delfin.dashboard import tab_settings


def _src() -> str:
    return inspect.getsource(tab_settings.create_tab)


def test_the_tab_has_a_diagnosis_section_driven_by_the_doctor():
    src = _src()
    assert "doctor as _doc" in src and "_doc.run_doctor(" in src, (
        "the tab must render the doctor's rows, not its own checks")
    assert "diagnosis_section" in src
    assert "main_accordion.set_title(9, 'Diagnosis & API key')" in src


def test_a_fix_that_installs_software_is_shown_not_run():
    """Three doctor rows carry a command and all three install or set up
    software. A button that installs is the person's decision taken from
    them."""
    src = _src()
    i = src.index("def _render_diagnosis(")
    block = src[i:i + 2000]
    assert "<code" in block, "a command must be shown to copy"
    for forbidden in ("subprocess", "os.system", "Popen"):
        assert forbidden not in block, (
            f"the diagnosis renderer must not run anything ({forbidden})")


def test_the_key_field_is_cleared_before_the_value_is_stored():
    src = _src()
    i = src.index("def _store_key(")
    block = src[i:i + 1200]
    clear = block.index("key_value.value = ''")
    store = block.index("set_credential(")
    assert clear < store, (
        "the field must be cleared first: a store that raises would "
        "otherwise leave the value in the browser")


def test_no_path_in_the_key_handling_echoes_the_value():
    src = _src()
    i = src.index("def _store_key(")
    block = src[i:i + 1200]
    # The error path may name the exception's TYPE, never its message: a
    # credential error's text can carry the value.
    assert "type(exc).__name__" in block
    assert "str(exc)" not in block, (
        "an exception message in the key path may carry the value")
    assert "value" not in block.split("_refresh_key_status")[-1], (
        "the status refresh must not be handed the value")


def test_the_status_names_the_source_and_warns_about_an_env_variable():
    src = _src()
    i = src.index("def _refresh_key_status(")
    block = src[i:i + 2200]
    assert "source" in block, "the source of each credential is shown"
    assert "takes precedence" in block, (
        "a key shadowed by an environment variable must be called out: a "
        "value stored here stays unused")


def test_the_stricter_route_is_named_in_the_tab():
    """The dashboard route sends the value through the browser and the
    websocket; the terminal's getpass does not. The tab says so instead of
    leaving the person to find out."""
    src = _src()
    assert "delfin-agent credentials set" in src
    assert "stricter" in src.lower()


# ---------------------------------------------------------------------------
# Driven: the tab is built and the button is pressed. Reading the source
# proves a line exists; the first version of this section passed every
# source check above and still rendered every row grey and listed all 24
# passing rows, because the doctor's statuses are UPPER CASE and the
# colour table and the filter were not. Only a press shows that.


def _build_tab(tmp_path):
    import ipywidgets as widgets  # noqa: F401
    from delfin.dashboard.context import DashboardContext
    from pathlib import Path

    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(
        calc_dir=Path(tmp_path / "calc"),
        archive_dir=Path(tmp_path / "archive"),
        office_dir=Path(tmp_path / "office"),
        agent_dir=Path(tmp_path / "ws"),
    )
    ctx.run_js = lambda *a, **k: None
    ctx.agent_workspace = str(tmp_path / "ws")
    tab, _refs = tab_settings.create_tab(ctx)
    return tab


def _walk(w):
    yield w
    for child in getattr(w, "children", ()) or ():
        yield from _walk(child)


def _button(tab, description):
    import ipywidgets as widgets

    for w in _walk(tab):
        if isinstance(w, widgets.Button) and w.description == description:
            return w
    raise AssertionError(f"no button {description!r} in the tab")


def _html_containing(tab, needle):
    import ipywidgets as widgets

    return [w for w in _walk(tab)
            if isinstance(w, widgets.HTML) and needle in str(w.value)]


def test_pressing_the_button_lists_only_what_needs_attention(
        tmp_path, monkeypatch):
    rows = [
        {"check": "fine", "status": "PASS", "detail": "all good", "fix": ""},
        {"check": "missing index", "status": "WARN",
         "detail": "no index", "fix": "build it",
         "command": "delfin-docs-index"},
        {"check": "broken", "status": "FAIL", "detail": "gone", "fix": "x"},
    ]
    from delfin.agent import doctor as _doc
    monkeypatch.setattr(_doc, "run_doctor", lambda *a, **k: list(rows))

    tab = _build_tab(tmp_path)
    assert _html_containing(tab, "Not run yet"), "the section starts unrun"
    _button(tab, "Run diagnosis").click()

    out = _html_containing(tab, "missing index")
    assert out, "the WARN row must be rendered"
    body = out[0].value
    assert "1 PASS" in body and "1 WARN" in body and "1 FAIL" in body
    assert "fine" not in body, (
        "a passing row must not be listed: the heading says only what needs "
        "attention is shown")
    # The status words are upper case; the colours must still be found.
    assert "#ef6c00" in body, "WARN must be orange, not the grey fallback"
    assert "#c62828" in body, "FAIL must be red, not the grey fallback"
    assert "delfin-docs-index" in body, "the command is shown to copy"


def test_pressing_store_with_an_empty_field_changes_nothing(tmp_path):
    tab = _build_tab(tmp_path)
    _button(tab, "Store key").click()
    assert _html_containing(tab, "nothing entered"), (
        "an empty field must say so rather than writing an empty credential")
