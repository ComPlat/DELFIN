"""The dashboard's "Skill proposals" review panel.

Package 6, phase 4. A self-written skill is only ever a PROPOSAL (the
learning wave's one non-negotiable difference to Hermes): invisible to
discovery until a human accepts it. The panel is where that human
looks — and it has to show everything the decision needs:

  * the list: name, status, why it qualified, its evidence;
  * the FULL text of the proposal, not a snippet — a review without
    the text under review is a rubber stamp;
  * the security findings (skill_safety.check output) in full;
  * buttons that really save: Accept / Reject call
    skill_proposals.accept / .reject, and a BLOCKED proposal has no
    Accept button at all (nothing unsafe becomes a rule by clicking).

The panel logic lives in delfin/dashboard/skill_proposals_panel.py,
pure functions over the contract types, so it is testable without
ipywidgets; tab_agent wires it into the Agent tab next to the
containment panel. skill_proposals itself (package 1) is imported
lazily — while it is not on this branch the tests stand in for it
with a contract-faithful fake, per the wave's interface agreement.
"""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass, field

import pytest


# ---------------------------------------------------------------------------
# A contract-faithful stand-in for delfin.agent.skill_proposals (Paket 1),
# exactly the shape the wave's interface agreement pins.
# ---------------------------------------------------------------------------

@dataclass
class _Evidence:
    kind: str
    ref: str
    detail: str = ""
    verified_at: str = ""


@dataclass
class _Proposal:
    name: str
    text: str
    evidence: list
    source: str
    status: str
    findings: list
    created: str
    base_version: str = ""


def _install_standin(monkeypatch, *, accept_calls=None, reject_calls=None,
                     proposals=None):
    mod = types.ModuleType("delfin.agent.skill_proposals")
    mod.Evidence = _Evidence
    mod.Proposal = _Proposal

    def list_proposals(status=None):
        return [p for p in proposals or [] if status is None
                or p.status == status]

    def accept(name, *, by):
        accept_calls.append((name, by))
        return _file_of(name)

    def get_proposal(name):
        for p in proposals or []:
            if p.name == name:
                return p
        return None

    def reject(name, *, reason, by):
        reject_calls.append((name, reason, by))
        return _file_of(name)

    def _file_of(name):
        import pathlib
        return pathlib.Path(f"/nonexistent/rejected/{name}.md")

    mod.list_proposals = list_proposals
    mod.get_proposal = get_proposal
    mod.accept = accept
    mod.reject = reject
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_proposals", mod)
    # ``from delfin.agent import skill_proposals`` resolves the attribute
    # on the package first, so patch that too — not only sys.modules.
    import delfin.agent as _pkg
    monkeypatch.setattr(_pkg, "skill_proposals", mod, raising=False)
    return mod


def _proposals_fixture():
    return [
        _Proposal(
            name="run-tests-first", text="# Run tests first\n1. Always\n",
            evidence=[_Evidence("test", "tests/test_x.py::test_y")],
            source="session abc123", status="pending",
            findings=[], created="2026-09-27T08:00:00"),
        _Proposal(
            name="skip-the-gate", text="# Skip the gate\n1. Bash past it\n",
            evidence=[_Evidence("test", "tests/test_g.py::test_z")],
            source="session abc123", status="blocked",
            findings=["proposes bypassing the sandbox gate"],
            created="2026-09-27T08:10:00"),
    ]


# ---------------------------------------------------------------------------
# The rows the panel renders
# ---------------------------------------------------------------------------

from delfin.dashboard import skill_proposals_panel as P  # noqa: E402


def test_the_panel_lists_every_proposal_with_status_and_evidence(
        monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    rows = P.collect_rows()
    by_name = {r["name"]: r for r in rows}
    assert set(by_name) == {"run-tests-first", "skip-the-gate"}
    assert by_name["run-tests-first"]["status"] == "pending"
    # The evidence is shown, not just counted: a review needs the ref.
    assert "tests/test_x.py::test_y" in by_name["run-tests-first"]["evidence"]
    assert by_name["run-tests-first"]["source"] == "session abc123"


def test_the_preview_carries_the_full_text_and_the_full_findings(
        monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    preview = P.render_preview("skip-the-gate")
    # The FULL text, verbatim — a review without the text under review
    # is a rubber stamp.
    assert "# Skip the gate" in preview
    assert "Bash past it" in preview
    # The security finding, in full, not elided.
    assert "proposes bypassing the sandbox gate" in preview


def test_the_buttons_really_save(monkeypatch):
    accepted: list = []
    rejected: list = []
    _install_standin(monkeypatch, accept_calls=accepted,
                     reject_calls=rejected, proposals=_proposals_fixture())
    # The widget path: the panel's actions ARE the contract functions,
    # so a click saves for real.
    P.accept_proposal("run-tests-first", by="dashboard-user")
    P.reject_proposal("skip-the-gate", reason="unsafe",
                      by="dashboard-user")
    assert accepted == [("run-tests-first", "dashboard-user")]
    assert rejected == [("skip-the-gate", "unsafe", "dashboard-user")]


def test_a_blocked_proposal_has_no_accept_button(monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    actions = P.available_actions("skip-the-gate")
    assert "accept" not in actions
    assert "reject" in actions


def test_a_pending_proposal_has_both_buttons(monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    assert set(P.available_actions("run-tests-first")) == {"accept",
                                                           "reject"}


def test_the_panel_html_renders_rows_and_findings(monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    html = P.format_panel_html()
    assert "run-tests-first" in html
    assert "pending" in html
    assert "blocked" in html
    assert "sandbox gate" in html
    # The blocked row must not carry an accept affordance in the HTML.
    import re
    blocked_segment = html.split("skip-the-gate", 1)[1][:2_000]
    assert not re.search(r"Accept", blocked_segment)


def test_no_proposals_is_said_not_blank(monkeypatch):
    _install_standin(monkeypatch, proposals=[])
    html = P.format_panel_html()
    assert "no skill proposals" in html.lower()


def test_a_read_failure_is_reported_not_blank(monkeypatch):
    mod = types.ModuleType("delfin.agent.skill_proposals")
    def boom(status=None):
        raise OSError("proposals dir unreadable")
    mod.list_proposals = boom
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_proposals", mod)
    html = P.format_panel_html()
    assert "could not be read" in html
    assert "proposals dir unreadable" in html


def test_tab_agent_wires_the_panel_in():
    """The panel is reachable from the Agent tab, wired like the
    containment panel (an HTML panel + a refresh hook in the tab)."""
    import inspect
    import delfin.dashboard.tab_agent as T
    src = inspect.getsource(T)
    assert "skill_proposals_panel" in src
    assert "skill_proposals_panel_html" in src
    assert "_refresh_skill_proposals_panel" in src


# ---------------------------------------------------------------------------
# Real buttons (the widget path), and escaping of foreign text
# ---------------------------------------------------------------------------

def _click(btn):
    """Trigger an ipywidgets button's handlers without a browser."""
    for cb in btn._click_handlers.callbacks:
        cb(btn)


def test_accept_button_click_saves_through_the_contract(monkeypatch):
    accepted: list = []
    refreshed: list = []
    _install_standin(monkeypatch, accept_calls=accepted,
                     proposals=_proposals_fixture())
    box = P.build_buttons("run-tests-first", on_change=lambda: refreshed.append(1))
    btns = [c for c in box.children if getattr(c, "description", "") == "Accept"]
    assert btns, "pending proposal must have an Accept button"
    _click(btns[0])
    assert accepted == [("run-tests-first", "dashboard-user")]
    assert refreshed == [1]  # the panel refreshes after the save


def test_blocked_proposal_has_no_accept_widget(monkeypatch):
    _install_standin(monkeypatch, proposals=_proposals_fixture())
    box = P.build_buttons("skip-the-gate")
    descs = [getattr(c, "description", "") for c in box.children]
    assert "Accept" not in descs
    assert "Reject" in descs


def test_reject_requires_a_reason_before_it_saves(monkeypatch):
    rejected: list = []
    _install_standin(monkeypatch, reject_calls=rejected,
                     proposals=_proposals_fixture())
    box = P.build_buttons("run-tests-first")
    reason_box = [c for c in box.children
                  if hasattr(c, "placeholder")][0]
    btn = [c for c in box.children
           if getattr(c, "description", "") == "Reject"][0]
    _click(btn)  # empty reason
    assert rejected == []
    reason_box.value = "not reusable enough"
    _click(btn)
    assert rejected == [("run-tests-first", "not reusable enough",
                         "dashboard-user")]


def test_foreign_text_is_escaped_everywhere(monkeypatch):
    """Proposal text can come from the team archive — other people's
    bytes. Nothing raw may reach an HTML widget."""
    evil = _Proposal(
        name="<script>alert(1)</script>",
        text="# Evil\n<img src=x onerror=alert(2)>\n",
        evidence=[_Evidence("test", "tests/x.py::<b>bold</b>")],
        source="team archive <img src=x onerror=alert(3)>",
        status="pending", findings=["<script>finding(4)</script>"],
        created="2026-09-27T09:00:00")
    _install_standin(monkeypatch, proposals=[evil])
    html = P.format_panel_html()
    # The dangerous forms are LIVE tags; escaping makes every "<" into
    # "&lt;" so none of these can occur. (A bare "onerror=" inside
    # escaped text is inert — there is no tag left to fire it.)
    for raw in ("<script>alert(1)", "<img src=x", "<b>bold</b>",
                "<script>finding"):
        assert raw not in html, f"unescaped foreign text reached HTML: {raw}"
    assert "&lt;script&gt;" in html  # and it IS shown, escaped
