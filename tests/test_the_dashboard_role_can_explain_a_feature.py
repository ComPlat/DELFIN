"""Asked "wie funktioniert aktuell der co2 coordinator" in dashboard
mode, the role had no way to the feature catalog: it guessed what the
coordinator might be and told the user to switch to Code by hand,
which its own prompt forbids (driven on DeepSeek, 2026-09-11; the same
question on GLM is a field report of the same day). The catalog tools
are read-only prose and now belong to the dashboard role, and its
prompt names them for exactly that question.
"""

from pathlib import Path

from delfin.agent.api_client import _DASHBOARD_AGENT_ALLOWED_TOOLS

_ROLE = (Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "agents"
         / "dashboard_agent.md")


def test_the_dashboard_role_may_read_the_feature_catalog():
    assert {"explain_delfin_feature", "list_delfin_features"} <= _DASHBOARD_AGENT_ALLOWED_TOOLS


def test_the_prompt_names_the_tool_and_the_handoff():
    text = _ROLE.read_text(encoding="utf-8")
    i = text.index("**`explain_delfin_feature`**")
    window = text[i:i + 300]
    assert "wie funktioniert X" in window
    assert "ACTION: /mode solo" in window
