# d2 — steering the model to DELFIN's own 3D viewer.

# The observed defect (2026-11-12): "show me benzene in 3D" got a raw
# <iframe srcdoc> with an NGL viewer from a CDN written into the answer text,
# although the agent tab renders its own interactive viewer card behind the
# show_molecule tool. Three fixes, tested here:
#
#   1. the tool schema must SAY when to use it and what never to do instead;
#   2. the system prompt (solo_agent.md decision table) must name the viewer;
#   3. a finished answer that embeds raw viewer HTML with no show_molecule
#      call this turn gets one visible hint pointing at the tool.

from __future__ import annotations

from delfin.agent import api_client as A
from delfin.agent import engine as E
from delfin.agent import show_molecule_guard as G


# --- 1. The schema description steers --------------------------------------


def _schema_description() -> str:
    for t in A._TOOL_DEFS:
        if t.get("function", {}).get("name") == "show_molecule":
            return t["function"]["description"]
    raise AssertionError("show_molecule not in _TOOL_DEFS")


def test_show_molecule_description_names_when_to_use_it():
    text = _schema_description().lower()
    # When-to-use: the request a user actually writes.
    assert "3d" in text
    assert "molecule" in text or "structure" in text
    # The card: what the tool produces instead of markup.
    assert "card" in text


def test_show_molecule_description_names_what_never_to_do():
    text = _schema_description().lower()
    # The exact observed failure: hand-written viewer HTML in the answer.
    assert "html" in text
    assert "iframe" in text or "script" in text


def test_schema_budget_unchanged():
    # The added steering words are paid for elsewhere: the total must not
    # rise over the committed ceiling (tests/test_tool_schema_budget.py
    # owns the budget; this pins that THIS change did not pay by raising).
    report = A.tool_schema_token_report()
    assert report["total_tokens"] <= 11_675, report["total_tokens"]


# --- 2. The system prompt names the viewer ---------------------------------


def test_solo_prompt_names_the_viewer_tool():
    from pathlib import Path
    text = (Path(E.__file__).parent / "pack" / "agents"
            / "solo_agent.md").read_text(encoding="utf-8").lower()
    assert "show_molecule" in text
    # The wrong way, named so the model can recognise it.
    assert "iframe" in text or "hand-written html" in text


# --- 3. The answer-time hint ------------------------------------------------


def test_scan_finds_raw_viewer_html():
    hit = G.scan_raw_viewer_html(
        "<p>Here you go:</p>\n<iframe srcdoc='&lt;script src=\"https://cdn"
        ".somewhere/3Dmol.js\"&gt;'></iframe>")
    assert hit is not None
    assert "iframe" in hit.lower()


def test_scan_passes_plain_answers():
    for text in (
        "Here is the structure, shown with the viewer card above.",
        "The xyz file holds 12 atoms; the C-C bond is 1.39 A.",
        "x = [i for i in range(3)]  # a code block, not markup",
        "<b>emphasis</b> in an answer is not a viewer",
    ):
        assert G.scan_raw_viewer_html(text) is None


def test_hint_names_the_tool_and_the_card():
    hint = G.raw_viewer_hint()
    assert "show_molecule" in hint
    assert "card" in hint.lower()


def test_engine_appends_hint_when_viewer_html_without_tool_call(monkeypatch):
    eng = E.AgentEngine.__new__(E.AgentEngine)  # no init: only the chain runs
    monkeypatch.setattr(eng, "_last_user_message",
                        "Show me benzene in 3D", raising=False)
    answer = ("<p>The molecule:</p>\n<iframe srcdoc='&lt;script "
              "src=\"https://cdn.example/3Dmol.js\"&gt;'></iframe>")
    out = G.caveat_for(answer, show_molecule_called=False)
    assert "show_molecule" in out
    # The answer itself is preserved: a hint, not a rewrite.
    assert "3Dmol.js" in out or "iframe" in out


def test_engine_quiet_when_tool_was_called():
    answer = "<p>The molecule:</p>\n<iframe srcdoc='x'></iframe>"
    assert G.caveat_for(answer, show_molecule_called=True) == ""


def test_engine_quiet_on_plain_answers():
    assert G.caveat_for("plain text, no markup",
                        show_molecule_called=False) == ""
