"""The Settings tab can apply a repair the diagnosis found.

The engine already existed: ``delfin/agent/repair.py`` turns a doctor row
into a step and applies one step per explicit approval. It was reachable
only as the Agent tab's ``/doctor repair`` chat command, and
``delfin-agent doctor --repair`` refuses outright when stdin is not a
terminal. The repair existed; the person who needed it could not reach
it without a terminal.

These controls are a thin caller of that engine. The steps come from
``repair.doctor_repair_plan`` over the same rows the report rendered, and
applying goes through ``repair.doctor_repair_apply``. A dashboard that
planned its own steps would be a second answer to "what can be fixed",
free to drift from the one that acts.

What is pinned here is the shape of the approval, because that is the
part a UI can quietly get wrong:

* Apply is disabled until the person has seen THIS step's what-and-undo;
* changing the selection disables it again, so a preview cannot be
  inherited by a different step;
* a click on a disabled-by-rule Apply runs nothing even if the widget's
  handler is invoked directly, which is what a stray click or a
  scripted caller does.

Driven, not read: the tab is built and the buttons are pressed. Reading
the source would pass for a version where Apply was never gated at all.
"""

from __future__ import annotations

import time

from delfin.dashboard import tab_settings

ROWS = [
    {"check": "fine", "status": "PASS", "detail": "ok", "fix": ""},
    {"check": "test runner", "status": "WARN", "detail": "no pytest",
     "fix": "install it", "command": "/bin/echo repaired-by-the-tab"},
    {"check": "loader path", "status": "WARN", "detail": "no handle",
     "fix": "a change to the machine's environment"},
]


def _build(tmp_path, rows, monkeypatch):
    import ipywidgets  # noqa: F401
    from pathlib import Path

    from delfin.agent import doctor as D
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    monkeypatch.setattr(D, "run_doctor", lambda *a, **k: list(rows))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()
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
    raise AssertionError(f"no button {description!r}")


def _pick(tab):
    import ipywidgets as widgets

    for w in _walk(tab):
        if isinstance(w, widgets.Dropdown) and any(
                "diagnosis" in str(o[0]) or "nothing here" in str(o[0])
                or "—" in str(o[0]) for o in (w.options or [])):
            return w
    raise AssertionError("the repair selector is missing")


def _html(tab, needle):
    import ipywidgets as widgets

    return [w for w in _walk(tab)
            if isinstance(w, widgets.HTML) and needle in str(w.value)]


def _await(tab, needles, budget_s=12.0):
    deadline = time.time() + budget_s
    while time.time() < deadline:
        for needle in needles:
            found = _html(tab, needle)
            if found:
                return found[0].value
        time.sleep(0.05)
    return ""


def test_only_a_row_the_engine_can_act_on_is_offered(tmp_path, monkeypatch):
    tab = _build(tmp_path, ROWS, monkeypatch)
    _button(tab, "Run diagnosis").click()

    labels = [o[0] for o in _pick(tab).options]
    assert len(labels) == 1, labels
    assert "test runner" in labels[0]
    assert not any("loader path" in q for q in labels), (
        "a prose-only remedy is not fixable by the engine and must not be "
        "offered as if it were")
    assert not any("fine" in q for q in labels), "a passing row is not a fix"


def test_apply_is_shut_until_the_step_has_been_previewed(
        tmp_path, monkeypatch):
    tab = _build(tmp_path, ROWS, monkeypatch)
    apply_btn = _button(tab, "Apply this fix")
    assert apply_btn.disabled, "shut before anything has been run"

    _button(tab, "Run diagnosis").click()
    assert apply_btn.disabled, "a plan is not an approval"

    _button(tab, "What would it do?").click()
    assert not apply_btn.disabled
    assert _html(tab, "undo:"), "the preview must state the undo"


def test_a_click_on_a_shut_apply_runs_nothing(tmp_path, monkeypatch):
    """Directly invoking the handler is what a stray click does."""
    tab = _build(tmp_path, ROWS, monkeypatch)
    _button(tab, "Run diagnosis").click()

    _button(tab, "Apply this fix").click()

    assert not _html(tab, "Applied."), "nothing may have run"
    assert not _html(tab, "Running"), "it must not even start"


def test_changing_the_selection_withdraws_the_approval(
        tmp_path, monkeypatch):
    rows = ROWS + [
        {"check": "second fixable", "status": "WARN", "detail": "x",
         "fix": "y", "command": "/bin/echo second"},
    ]
    tab = _build(tmp_path, rows, monkeypatch)
    _button(tab, "Run diagnosis").click()
    _button(tab, "What would it do?").click()
    apply_btn = _button(tab, "Apply this fix")
    assert not apply_btn.disabled

    pick = _pick(tab)
    pick.value = [o[1] for o in pick.options][1]

    assert apply_btn.disabled, (
        "a preview belongs to the step it was shown for; it must not be "
        "inherited by another")


def test_applying_runs_the_step_and_reports_its_output(
        tmp_path, monkeypatch):
    tab = _build(tmp_path, ROWS, monkeypatch)
    _button(tab, "Run diagnosis").click()
    _button(tab, "What would it do?").click()

    _button(tab, "Apply this fix").click()

    body = _await(tab, ["Applied.", "failed"])
    assert "Applied." in body, body[:300]
    assert "repaired-by-the-tab" in body, (
        "the step's own output is what proves it ran")
    assert "Run the diagnosis again" in body, (
        "a repair that is not re-checked is a claim, not a result")


def test_a_failing_step_is_reported_not_swallowed(tmp_path, monkeypatch):
    rows = [{"check": "broken fix", "status": "FAIL", "detail": "x",
             "fix": "y", "command": "/bin/false"}]
    tab = _build(tmp_path, rows, monkeypatch)
    _button(tab, "Run diagnosis").click()
    _button(tab, "What would it do?").click()

    _button(tab, "Apply this fix").click()

    body = _await(tab, ["Applied.", "failed", "exit"])
    assert body, "the outcome must reach the person either way"


def test_the_section_says_which_rows_it_cannot_apply(tmp_path, monkeypatch):
    """A control that silently covers only some rows invites the reader to
    assume it covers all of them."""
    tab = _build(tmp_path, ROWS, monkeypatch)
    explain = _html(tab, "remedy is a command or a setting")
    assert explain, "the limit of this control must be stated in it"


def test_the_plan_comes_from_the_engine_not_from_the_tab():
    import inspect

    src = inspect.getsource(tab_settings.create_tab)
    assert "doctor_repair_plan" in src and "doctor_repair_apply" in src, (
        "the tab must call the engine's seam, not re-implement planning")
    i = src.index("def _repair_refresh(")
    block = src[i:i + 1200]
    for forbidden in ("subprocess", "os.system", "Popen"):
        assert forbidden not in block, (
            f"the tab must not run anything itself ({forbidden})")
