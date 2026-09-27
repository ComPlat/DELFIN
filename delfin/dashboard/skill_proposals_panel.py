"""The dashboard's "Skill proposals" review panel (package 6, phase 4).

A self-written skill is only ever a PROPOSAL: invisible to discovery
until a human accepts it. This module renders everything that human
needs to decide — the list with status, qualification reasons and
evidence, the FULL text, the security findings — and its actions ARE
the contract functions (``skill_proposals.accept`` / ``.reject`` from
package 1), so a click saves for real. A ``blocked`` proposal has no
accept action at all.

Pure functions, deliberately: no ipywidgets import here, so the logic
is testable headless. ``delfin.agent.skill_proposals`` is imported
lazily inside each call — a build without package 1 shows "could not
be read" rather than breaking the Agent tab.
"""

from __future__ import annotations

from html import escape
from typing import Optional


def _proposals_module():
    from delfin.agent import skill_proposals
    return skill_proposals


def collect_rows() -> list[dict]:
    """One row per proposal: name, status, source, evidence refs."""
    proposals = _proposals_module().list_proposals()
    rows = []
    for p in proposals:
        evidence = "; ".join(
            f"{getattr(e, 'kind', '?')}: {getattr(e, 'ref', '')}"
            for e in (p.evidence or []))
        rows.append({
            "name": p.name,
            "status": p.status,
            "source": p.source,
            "evidence": evidence,
            "findings": list(p.findings or []),
            "created": p.created,
        })
    return rows


def render_preview(name: str) -> str:
    """The full proposal text and the full security findings, verbatim."""
    mod = _proposals_module()
    p = mod.get_proposal(name)
    if p is None:
        return f"(proposal '{name}' not found)"
    findings = "\n".join(f"- {f}" for f in (p.findings or []))
    return (
        f"--- Skill proposal: {p.name} (status: {p.status}, "
        f"from {p.source}) ---\n"
        f"{p.text}\n"
        f"--- Security findings ---\n"
        f"{findings if findings else '(none)'}\n"
    )


def available_actions(name: str) -> list[str]:
    """Accept is offered only where accepting is allowed at all."""
    p = _proposals_module().get_proposal(name)
    if p is None:
        return []
    if p.status == "blocked":
        return ["reject"]
    if p.status == "pending":
        return ["accept", "reject"]
    return []  # accepted / rejected: decided already


def accept_proposal(name: str, *, by: str):
    """The Accept button's handler — the contract function itself."""
    return _proposals_module().accept(name, by=by)


def reject_proposal(name: str, *, reason: str, by: str):
    """The Reject button's handler — the contract function itself."""
    return _proposals_module().reject(name, reason=reason, by=by)


def build_buttons(name: str, *, on_change=None):
    """Real ipywidgets Accept / Reject buttons for one proposal.

    The user's standing rule for approval UIs: buttons that actually
    save. Accept calls ``skill_proposals.accept`` directly; Reject
    asks for a reason (an empty reason does NOT save — the archive
    should never say "no" without saying why) and calls ``.reject``.
    A ``blocked`` proposal gets NO accept button. ``on_change`` runs
    after a successful save so the caller can refresh its panels.
    """
    import ipywidgets as _widgets
    actions = set(available_actions(name))
    children = []
    if "accept" in actions:
        def _on_accept(_btn):
            accept_proposal(name, by="dashboard-user")
            if on_change is not None:
                on_change()
        accept_btn = _widgets.Button(
            description="Accept", icon="check",
            button_style="success", tooltip=(
                f"Accept skill proposal '{name}' — it becomes a real, "
                "discoverable skill"))
        accept_btn.on_click(_on_accept)
        children.append(accept_btn)
    if "reject" in actions:
        reason_box = _widgets.Text(
            value="", placeholder="reason (required)",
            layout=_widgets.Layout(width="260px"))

        def _on_reject(_btn):
            reason = (reason_box.value or "").strip()
            if not reason:
                reason_box.placeholder = "reason REQUIRED before rejecting"
                return
            reject_proposal(name, reason=reason, by="dashboard-user")
            if on_change is not None:
                on_change()
        reject_btn = _widgets.Button(
            description="Reject", icon="times",
            button_style="danger", tooltip=(
                f"Reject skill proposal '{name}' — moved to rejected/"))
        reject_btn.on_click(_on_reject)
        children.extend([reason_box, reject_btn])
    return _widgets.HBox(children)


_STATUS_COLOURS = {"pending": "#f59e0b", "blocked": "#d32f2f",
                   "accepted": "#2e7d32", "rejected": "#888"}


def format_panel_html(limit: int = 20) -> str:
    """The panel's HTML: list, full previews, findings, button affordances."""
    try:
        rows = collect_rows()
    except Exception as exc:
        return (
            "<div style='border-left:3px solid #d32f2f;padding:6px 10px;"
            "background:#d32f2f11;border-radius:4px;font-size:12px;"
            "color:#d32f2f;'>🛡 Skill proposals: the list could not be "
            f"read — {escape(str(exc))}</div>"
        )
    if not rows:
        return (
            "<div style='color:#888;font-size:12px;font-style:italic;'>"
            "📚 Skill proposals: no skill proposals awaiting review — "
            "none proposed this session.</div>"
        )
    parts: list[str] = []
    for r in rows[:limit]:
        colour = _STATUS_COLOURS.get(r["status"], "#888")
        actions = available_actions(r["name"])
        buttons = ""
        if "accept" in actions:
            buttons += ("<b style='color:#2e7d32'>[Accept]</b> ")
        if "reject" in actions:
            buttons += "<b style='color:#888'>[Reject]</b>"
        findings = "; ".join(r["findings"])
        findings_html = (
            f" <span style='color:#d32f2f'>⚠ {escape(findings)}</span>"
            if findings else "")
        parts.append(
            f"<div style='font-family:monospace;font-size:12px;"
            f"margin:2px 0;'>📚 <b>{escape(r['name'])}</b> "
            f"<span style='color:{colour}'>[{escape(r['status'])}]</span>"
            f" <span style='color:#888'>from {escape(r['source'])}, "
            f"evidence: {escape(r['evidence'])}</span>{findings_html} "
            f"{buttons}</div>"
        )
        # The full preview belongs to every row: the review needs the
        # text under review, not a snippet.
        parts.append(
            f"<pre style='font-size:11px;color:#555;margin:0 0 6px 18px;"
            f"max-height:220px;overflow-y:auto;'>"
            f"{escape(render_preview(r['name']))}</pre>"
        )
    return (
        "<div style='border-left:3px solid #f59e0b;padding:6px 10px;"
        "background:#f59e0b11;border-radius:4px;'>"
        "<div style='font-size:11px;color:#aaa;margin-bottom:3px;'>"
        f"📚 Skill proposals &nbsp; {len(rows)} awaiting review</div>"
        + "".join(parts) + "</div>"
    )


__all__ = [
    "collect_rows", "render_preview", "available_actions",
    "accept_proposal", "reject_proposal", "format_panel_html",
]
