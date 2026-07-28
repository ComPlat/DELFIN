"""Rendered deliverables: durable markdown+HTML reports under the
workspace instead of chat-only prose."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.deliverables import md_to_html, publish_report


def test_md_to_html_renders_and_escapes():
    md = ("# Title\n\nSome **bold** and `code`.\n\n"
          "- item one\n- item two\n\n"
          "```\nraw <script>alert(1)</script>\n```\n\n"
          "| a | b |\n|---|---|\n| 1 | 2 |\n")
    html = md_to_html(md, title="Title")
    assert "<h1>Title</h1>" in html
    assert "<strong>bold</strong>" in html
    assert "<code>code</code>" in html
    assert "<li>item one</li>" in html
    assert "<td>1</td>" in html
    assert "<script>" not in html          # escaped, never raw
    assert "&lt;script&gt;" in html


def test_publish_report_writes_both_without_overwrite(tmp_path):
    out1 = publish_report(tmp_path, title="S1 Analysis",
                          markdown="result body", fmt="both")
    assert Path(out1["md_path"]).is_file()
    assert Path(out1["html_path"]).is_file()
    assert "s1-analysis" in out1["md_path"]
    out2 = publish_report(tmp_path, title="S1 Analysis",
                          markdown="second run", fmt="md")
    assert out2["md_path"] != out1["md_path"]      # never clobbered
    assert Path(out1["md_path"]).read_text().count("result body") == 1


def test_publish_report_tool_happy_and_plan_mode(tmp_path):
    from delfin.agent.api_client import KitToolPermissions
    ws = tmp_path / "ws"; ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    out = json.loads(A._doc_executor.execute("publish_report", {
        "title": "Audit", "markdown": "## ok\nfine"}, permissions=perms))
    assert out["status"] == "written"
    assert (ws / "reports").is_dir()
    perms.mode = "plan"
    blocked = json.loads(A._doc_executor.execute("publish_report", {
        "title": "X", "markdown": "y"}, permissions=perms))
    assert "plan mode" in blocked["error"]


def test_publish_report_tool_validation(tmp_path):
    from delfin.agent.api_client import KitToolPermissions
    ws = tmp_path / "ws"; ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    missing = json.loads(A._doc_executor.execute("publish_report", {
        "title": "", "markdown": ""}, permissions=perms))
    assert "error" in missing
    badfmt = json.loads(A._doc_executor.execute("publish_report", {
        "title": "t", "markdown": "m", "format": "pdf"}, permissions=perms))
    assert "format" in badfmt["error"]
