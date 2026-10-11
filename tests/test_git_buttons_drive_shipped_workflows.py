"""The top row drives git through workflows; the bug report keeps a copy.

Requested by Tilmann (2026-10-09): Merge / Branch / Push / PR buttons at
the top, so a person can tell the agent what to do by pressing one. Each
button sends a slash command, and the shipped skill of that name is the
playbook the model follows -- the same text works typed in the terminal.
The base branch is chosen, not always main. Export and Export as notebook
left the row: sending a bug report now also downloads a copy of it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import skills
from delfin.agent.api_client import _asks_for_push

SRC = Path(__file__).resolve().parents[1].joinpath(
    "delfin", "dashboard", "tab_agent.py").read_text()


@pytest.mark.parametrize("name", ["merge", "new-branch", "push", "pr"])
def test_every_button_has_its_shipped_workflow(name):
    sk = skills.get_skill(name)
    assert sk is not None, f"no shipped skill /{name}"
    assert "pack" in str(sk.source)


def test_only_push_and_pr_ask_for_a_push():
    """A button is the user's request -- for exactly what it says. Merge
    and Branch must not hand the model a push grant on the way."""
    asked = {n: _asks_for_push(skills.render_skill_invocation(
                 skills.get_skill(n), "origin/main"))
             for n in ("merge", "new-branch", "push", "pr")}
    assert asked == {"merge": False, "new-branch": False,
                     "push": True, "pr": True}


def test_merge_and_pr_take_the_chosen_base_branch():
    assert "arguments" in skills.get_skill("merge").body
    assert "--base <base>" in skills.get_skill("pr").body
    assert "`from`" in skills.get_skill("new-branch").body


def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    tab = tab_agent.create_tab(ctx)
    return tab[0] if isinstance(tab, tuple) else tab


def _walk(node):
    stack = [node]
    while stack:
        n = stack.pop()
        yield n
        stack.extend(getattr(n, "children", ()) or ())


def test_the_row_has_the_git_buttons_and_no_export(tmp_path):
    import ipywidgets as widgets
    tab = _tab(tmp_path)
    buttons = {b.description: b for b in _walk(tab)
               if isinstance(b, widgets.Button)}
    for label in ("⤓ Merge", "⑂ Branch", "⤒ Push", "⇄ PR"):
        assert label in buttons, f"missing {label}"
    for label in ("Export", "Export as notebook"):
        if label in buttons:
            assert buttons[label].layout.display == "none", label
    bases = [w for w in _walk(tab) if isinstance(w, widgets.Combobox)
             and w.placeholder == "base branch"]
    assert bases and bases[0].value == "origin/main"


def test_sending_a_bug_report_downloads_a_copy():
    i = SRC.index("def _on_bug_report(button):")
    body = SRC[i:i + 6000]
    assert "_download_report_copy(report_dir)" in body
    j = SRC.index("def _download_report_copy(report_dir)")
    helper = SRC[j:j + 3000]
    assert "zipfile" in helper and "ctx.run_js" in helper
    assert "_REPORT_DOWNLOAD_MAX" in helper
