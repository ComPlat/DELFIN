"""Visual-mode tests for the dashboard activity spinner.

Each mode produces a distinct CSS-class + colour combination so the user
can tell at a glance whether the agent is actively streaming, blocked on
approval, idle-but-queued, or stuck.
"""

from __future__ import annotations

import re

import pytest


@pytest.fixture
def working_html_box():
    """Minimal HTML container that mimics ``widgets.HTML.value`` set semantics."""

    class _Box:
        def __init__(self) -> None:
            self.value = ""

    return _Box()


@pytest.fixture
def status_chip():
    class _Chip:
        def __init__(self) -> None:
            self.value = ""

    return _Chip()


def _build_set_working(working_html, status_chip):
    """Replicate the production ``_set_working`` body in isolation so we can
    drive every mode without spinning up the full Dashboard."""

    import html as _html

    def _set_working(active: bool, label: str = "", mode: str = "streaming") -> None:
        if active:
            text = label or "Working..."
            _esc = _html.escape
            mode_overrides = {
                "gated":  ("⏸️", "gated",  "delfin-agent-working--gated"),
                "queued": ("⏳",       "queued", "delfin-agent-working--queued"),
                "stale":  ("⚠️", "stale",  "delfin-agent-working--stale"),
            }
            if mode in mode_overrides:
                icon, cat, mode_cls = mode_overrides[mode]
            else:
                mode_cls = ""
                icon, cat = "💭", "think"

            classes = "delfin-agent-working" + (f" {mode_cls}" if mode_cls else "")
            working_html.value = (
                f'<div class="{classes}">'
                f'<span class="delfin-spinner">'
                f'<span class="delfin-spinner-dot"></span>'
                f'<span class="delfin-spinner-dot"></span>'
                f'<span class="delfin-spinner-dot"></span>'
                f'</span>'
                f'<span class="delfin-activity-icon">{icon}</span>'
                f'<span class="delfin-activity-label">{_esc(cat)}</span>'
                f'<span class="delfin-activity-text">{_esc(text)}</span>'
                f'</div>'
            )
            header_color = {
                "gated": "#f59e0b", "queued": "#94a3b8", "stale": "#ef4444",
            }.get(mode, "#60a5fa")
            short = text[:40] + ("..." if len(text) > 40 else "")
            status_chip.value = (
                f'<span style="font-family:monospace; color:{header_color}; '
                f'padding:2px 6px; background:#1e293b; '
                f'border-radius:4px; font-size:11px;">'
                f'{icon} {_esc(short)}</span>'
            )
        else:
            working_html.value = '<img style="display:none">'
            status_chip.value = ""

    return _set_working


def test_set_working_streaming_uses_default_blue_color(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "Thinking...")
    assert "delfin-agent-working--" not in working_html_box.value
    assert "delfin-spinner-dot" in working_html_box.value
    assert "color:#60a5fa" in status_chip.value


def test_set_working_gated_marks_orange_paused(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "Waiting: Permission Gate", mode="gated")
    assert "delfin-agent-working--gated" in working_html_box.value
    assert "gated" in working_html_box.value
    assert "#f59e0b" in status_chip.value


def test_set_working_queued_marks_grey_idle(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "Idle - 3 queued", mode="queued")
    assert "delfin-agent-working--queued" in working_html_box.value
    assert "queued" in working_html_box.value
    assert "#94a3b8" in status_chip.value


def test_set_working_stale_marks_red_warning(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "Stream stalled (12 min, no output)", mode="stale")
    assert "delfin-agent-working--stale" in working_html_box.value
    assert "stale" in working_html_box.value
    assert "#ef4444" in status_chip.value


def test_set_working_false_clears_widget_and_header(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "Working", mode="streaming")
    sw(False)
    assert "delfin-agent-working" not in working_html_box.value
    assert status_chip.value == ""


def test_unknown_mode_falls_back_to_streaming(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "anything", mode="bogus")
    # no mode-specific modifier class for unknown modes
    assert "delfin-agent-working--" not in working_html_box.value
    # falls back to blue
    assert "#60a5fa" in status_chip.value


def test_text_label_html_escaped(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "<script>alert(1)</script>")
    assert "<script>alert" not in working_html_box.value
    assert "&lt;script&gt;" in working_html_box.value


def test_long_label_truncated_in_header_with_ellipsis(working_html_box, status_chip):
    sw = _build_set_working(working_html_box, status_chip)
    sw(True, "a" * 100)
    # widget keeps full text but header uses 40-char prefix + ellipsis
    assert "..." in status_chip.value or "…" in status_chip.value


def test_css_classes_are_disjoint_across_modes(working_html_box, status_chip):
    """Modifier classes must not collide so the styles compose cleanly."""
    sw = _build_set_working(working_html_box, status_chip)
    classes_per_mode = {}
    for mode in ("streaming", "gated", "queued", "stale"):
        sw(True, mode, mode=mode)
        m = re.search(r'class="([^"]+)"', working_html_box.value)
        assert m, f"no class attribute found for mode={mode}"
        classes_per_mode[mode] = set(m.group(1).split())
    # All four share the base class
    base = {"delfin-agent-working"}
    for mode, cls in classes_per_mode.items():
        assert base.issubset(cls), mode
    # And only the explicit non-streaming modes carry their modifier
    assert classes_per_mode["streaming"] == base
    assert "delfin-agent-working--gated"  in classes_per_mode["gated"]
    assert "delfin-agent-working--queued" in classes_per_mode["queued"]
    assert "delfin-agent-working--stale"  in classes_per_mode["stale"]
