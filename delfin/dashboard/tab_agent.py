"""DELFIN Agent tab: AI-powered chat assistant in the dashboard."""

from __future__ import annotations

import html as _html
import os
import re
import shutil
import threading
import time
from pathlib import Path

import ipywidgets as widgets


# ---------------------------------------------------------------------------
# CSS for the chat interface
# ---------------------------------------------------------------------------

_AGENT_CSS = """\
<style>
.delfin-agent-chat {
    max-height: calc(100vh - 460px);
    overflow-y: auto;
    overflow-anchor: auto;
    padding: 10px;
    border: 1px solid #ddd;
    border-radius: 6px;
    background: #fafafa;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    font-size: 13px;
    line-height: 1.5;
    display: flex;
    flex-direction: column;
}
.delfin-chat-msg {
    overflow-anchor: none;
    margin-bottom: 8px;
    padding: 8px 12px;
    border-radius: 8px;
    max-width: 85%;
    word-wrap: break-word;
}
.delfin-chat-tool + .delfin-chat-tool {
    margin-top: -6px;
}
.delfin-chat-user {
    background: #dbeafe;
    margin-left: auto;
    text-align: left;
    border-bottom-right-radius: 2px;
}
.delfin-chat-agent {
    background: #f3f4f6;
    margin-right: auto;
    border-bottom-left-radius: 2px;
}
.delfin-chat-system {
    background: #f3f4f6;
    margin: 4px 0;
    text-align: left;
    font-size: 11px;
    max-width: 100%;
    padding: 4px 10px;
    color: #6b7280;
    border-left: 3px solid #d1d5db;
    border-radius: 0 4px 4px 0;
    font-family: 'SF Mono', 'Consolas', 'Monaco', monospace;
    white-space: pre-wrap;
    word-break: break-all;
}
.delfin-chat-tool {
    background: transparent;
    color: #9ca3af;
    margin: 0;
    padding: 0 10px;
    font-size: 11px;
    max-width: 100%;
    border-radius: 0;
    font-family: 'SF Mono', 'Cascadia Code', 'Consolas', monospace;
    white-space: pre-wrap;
    word-break: break-all;
    line-height: 1.4;
}
.delfin-chat-tool .tool-name { color: #6b7280; font-weight: 600; }
.delfin-chat-tool .tool-path { color: #6b7280; }
.delfin-chat-tool .tool-param { color: #9ca3af; }
.delfin-chat-tool .tool-diff-old { color: #f7768e; }
.delfin-chat-tool .tool-diff-new { color: #9ece6a; }
.delfin-chat-tool details { cursor: pointer; }
.delfin-chat-tool details summary { display: inline; }
.delfin-chat-tool details pre {
    margin: 2px 0 0 16px;
    padding: 4px 8px;
    background: #f3f4f6;
    border-radius: 3px;
    font-size: 10px;
    color: #6b7280;
    max-height: 120px;
    overflow-y: auto;
}
.delfin-streaming-pre {
    margin: 0;
    padding: 0;
    font-family: 'SF Mono', 'Cascadia Code', 'Consolas', monospace;
    font-size: 13px;
    line-height: 1.5;
    white-space: pre-wrap;
    word-wrap: break-word;
    color: inherit;
    background: transparent;
}
.delfin-chat-approval {
    background: #fef3c7;
    margin: 8px 0;
    text-align: left;
    font-size: 12px;
    max-width: 100%;
    padding: 10px 14px;
    border-left: 4px solid #f59e0b;
    border-radius: 0 8px 8px 0;
    color: #92400e;
    font-weight: 500;
}
.delfin-agent-approval-row {
    background: #fffbeb;
    border: 1px solid #f59e0b;
    border-radius: 8px;
    align-items: center;
    gap: 8px;
    padding: 8px 12px;
    margin: 4px 0;
    position: sticky;
    bottom: 0;
    z-index: 10;
}
.delfin-chat-handoff {
    background: #fef3c7;
    margin: 8px auto;
    text-align: center;
    font-style: italic;
    font-size: 12px;
    max-width: 70%;
    padding: 8px 12px;
    border-left: none;
    border-radius: 8px;
    color: #92400e;
}
.delfin-chat-role {
    display: block;
    font-size: 11px;
    font-weight: 600;
    color: #6b7280;
    margin-bottom: 6px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
.delfin-chat-agent pre {
    background: #1e1e2e;
    color: #cdd6f4;
    padding: 10px 12px;
    border-radius: 6px;
    overflow-x: auto;
    font-size: 12px;
    font-family: 'SF Mono', 'Cascadia Code', 'Fira Code', 'Consolas', monospace;
    line-height: 1.45;
    margin: 6px 0;
    border: 1px solid #313244;
}
.delfin-chat-agent pre .code-lang {
    display: block;
    font-size: 10px;
    color: #6c7086;
    margin-bottom: 4px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
.delfin-chat-agent code {
    background: #e5e7eb;
    padding: 1px 4px;
    border-radius: 3px;
    font-size: 12px;
    font-family: 'SF Mono', 'Cascadia Code', 'Fira Code', 'Consolas', monospace;
}
.delfin-chat-agent pre code {
    background: none;
    padding: 0;
    color: inherit;
}
/* Syntax highlighting (dark theme) */
.delfin-chat-agent pre .kw { color: #cba6f7; }
.delfin-chat-agent pre .str { color: #a6e3a1; }
.delfin-chat-agent pre .num { color: #fab387; }
.delfin-chat-agent pre .cmt { color: #6c7086; font-style: italic; }
.delfin-chat-agent pre .fn { color: #89b4fa; }
.delfin-chat-agent pre .op { color: #89dceb; }
.delfin-chat-agent pre .dec { color: #f9e2af; }
/* Diff highlighting */
.delfin-chat-agent pre .diff-add { color: #a6e3a1; }
.delfin-chat-agent pre .diff-del { color: #f38ba8; }
.delfin-chat-agent pre .diff-hdr { color: #89b4fa; font-weight: 600; }
.delfin-agent-queue {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 10px;
    background: #dbeafe;
    color: #1e40af;
    font-size: 11px;
    font-weight: 600;
    margin: 4px 0;
}
.delfin-agent-status {
    font-size: 12px;
    color: #6b7280;
    padding: 4px 8px;
    background: #f9fafb;
    border-radius: 4px;
    border: 1px solid #e5e7eb;
}
.delfin-agent-status .mode-badge {
    display: inline-block;
    padding: 1px 8px;
    border-radius: 10px;
    background: #dbeafe;
    color: #1e40af;
    font-weight: 600;
    font-size: 11px;
    margin-right: 8px;
}
.delfin-agent-status .role-badge {
    display: inline-block;
    padding: 1px 8px;
    border-radius: 10px;
    background: #d1fae5;
    color: #065f46;
    font-weight: 600;
    font-size: 11px;
    margin-right: 8px;
}
.delfin-agent-status .backend-badge {
    display: inline-block;
    padding: 1px 8px;
    border-radius: 10px;
    background: #ede9fe;
    color: #5b21b6;
    font-weight: 600;
    font-size: 11px;
    margin-right: 8px;
}
.delfin-agent-status .tokens-info {
    color: #9ca3af;
    font-size: 11px;
}
.delfin-agent-status .gate-badge {
    display: inline-block;
    padding: 1px 8px;
    border-radius: 10px;
    font-weight: 700;
    font-size: 11px;
    margin-right: 8px;
}
.delfin-agent-status .gate-text {
    color: #92400e;
    font-size: 11px;
    font-weight: 600;
}
.delfin-cycle-inspector {
    margin: 6px 0 8px 0;
    padding: 10px 12px;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
    background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%);
}
.delfin-cycle-inspector .inspector-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 10px;
    margin-bottom: 10px;
    flex-wrap: wrap;
}
.delfin-cycle-inspector .inspector-title {
    font-size: 12px;
    font-weight: 800;
    color: #111827;
    letter-spacing: 0.4px;
    text-transform: uppercase;
}
.delfin-cycle-inspector .inspector-meta {
    font-size: 11px;
    color: #6b7280;
}
.delfin-cycle-inspector .inspector-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 8px;
}
.delfin-cycle-actions {
    margin: 0 0 8px 0;
    padding: 6px 0 2px 0;
    gap: 8px;
    flex-wrap: wrap;
    align-items: center;
}
.delfin-cycle-actions button {
    min-width: 110px;
}
.delfin-cycle-actions-info {
    font-size: 11px;
    color: #6b7280;
    padding-left: 2px;
}
.delfin-cycle-detail-box {
    margin: 0 0 8px 0;
    padding: 10px 12px;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
    background: #ffffff;
}
.delfin-cycle-detail-head {
    display: flex;
    align-items: center;
    gap: 10px;
    flex-wrap: wrap;
    margin-bottom: 8px;
}
.delfin-cycle-detail-label {
    font-size: 11px;
    font-weight: 800;
    color: #6b7280;
    text-transform: uppercase;
    letter-spacing: 0.35px;
}
.delfin-cycle-detail {
    border-top: 1px solid #e5e7eb;
    padding-top: 10px;
}
.delfin-cycle-detail .detail-title {
    font-size: 13px;
    font-weight: 800;
    color: #111827;
    margin-bottom: 4px;
}
.delfin-cycle-detail .detail-meta {
    font-size: 11px;
    color: #6b7280;
    margin-bottom: 8px;
}
.delfin-cycle-detail .detail-body {
    font-size: 12px;
    color: #1f2937;
    line-height: 1.5;
    white-space: normal;
}
.delfin-cycle-detail .detail-body pre {
    max-height: 260px;
    overflow: auto;
}
.delfin-cycle-inspector .inspector-card {
    min-height: 84px;
    padding: 10px 12px;
    border-radius: 10px;
    border: 1px solid #e5e7eb;
    background: #ffffff;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}
.delfin-cycle-inspector .inspector-card.gate-card {
    background: #fffbeb;
    border-color: #fcd34d;
}
.delfin-cycle-inspector .inspector-card.risk-card {
    background: #fff7ed;
    border-color: #fdba74;
}
.delfin-cycle-inspector .inspector-card.retry-card {
    background: #eff6ff;
    border-color: #93c5fd;
}
.delfin-cycle-inspector .inspector-card.history-card {
    background: #f8fafc;
    border-color: #cbd5e1;
}
.delfin-cycle-inspector .card-label {
    display: block;
    margin-bottom: 6px;
    font-size: 10px;
    font-weight: 800;
    letter-spacing: 0.4px;
    text-transform: uppercase;
    color: #6b7280;
}
.delfin-cycle-inspector .card-value {
    display: block;
    font-size: 13px;
    font-weight: 700;
    color: #111827;
    line-height: 1.35;
}
.delfin-cycle-inspector .card-note {
    display: block;
    margin-top: 6px;
    font-size: 11px;
    color: #6b7280;
    line-height: 1.4;
    white-space: pre-wrap;
}
.delfin-cycle-inspector .card-list {
    margin: 6px 0 0 0;
    padding-left: 16px;
    color: #374151;
    font-size: 11px;
    line-height: 1.45;
}
.delfin-cycle-inspector .card-list li {
    margin: 0 0 4px 0;
}
.delfin-cycle-inspector .history-list {
    margin: 8px 0 0 0;
    padding: 0;
    list-style: none;
}
.delfin-cycle-inspector .history-item {
    padding: 7px 0;
    border-top: 1px solid #e5e7eb;
}
.delfin-cycle-inspector .history-item:first-child {
    border-top: none;
    padding-top: 0;
}
.delfin-cycle-inspector .history-kind {
    display: inline-block;
    margin-right: 6px;
    padding: 1px 6px;
    border-radius: 999px;
    background: #e5e7eb;
    color: #374151;
    font-size: 10px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 0.35px;
}
.delfin-cycle-inspector .history-title {
    font-size: 11px;
    font-weight: 700;
    color: #111827;
}
.delfin-cycle-inspector .history-detail {
    margin-top: 4px;
    font-size: 11px;
    color: #6b7280;
    line-height: 1.4;
    white-space: pre-wrap;
}
.delfin-cycle-inspector .gate-pill {
    display: inline-block;
    margin-bottom: 6px;
    padding: 2px 8px;
    border-radius: 999px;
    font-size: 10px;
    font-weight: 800;
    letter-spacing: 0.4px;
    text-transform: uppercase;
    background: #8b5cf6;
    color: white;
}
.delfin-gate-card {
    margin: 8px 0;
    max-width: 100%;
    padding: 10px 14px;
    border-radius: 10px;
    border-left: 4px solid #9ca3af;
    background: #f9fafb;
    color: #374151;
}
.delfin-gate-card .gate-header {
    display: flex;
    align-items: center;
    gap: 8px;
    margin-bottom: 6px;
    flex-wrap: wrap;
}
.delfin-gate-card .gate-title {
    font-size: 12px;
    font-weight: 700;
}
.delfin-gate-card .gate-role {
    font-size: 11px;
    font-weight: 600;
    color: #6b7280;
    text-transform: uppercase;
    letter-spacing: 0.4px;
}
.delfin-gate-card .gate-body {
    white-space: pre-wrap;
    font-size: 12px;
    line-height: 1.45;
}
.delfin-gate-card .gate-hint {
    margin-top: 8px;
    font-size: 11px;
    color: #6b7280;
}
.delfin-gate-badge {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 999px;
    font-size: 10px;
    font-weight: 800;
    letter-spacing: 0.4px;
    text-transform: uppercase;
}
.delfin-gate-risk {
    background: #fef3c7;
    border-left-color: #f59e0b;
}
.delfin-gate-risk .delfin-gate-badge {
    background: #f59e0b;
    color: #fff7ed;
}
.delfin-gate-schema {
    background: #fee2e2;
    border-left-color: #ef4444;
}
.delfin-gate-schema .delfin-gate-badge {
    background: #ef4444;
    color: white;
}
.delfin-gate-question,
.delfin-gate-plan-approval,
.delfin-gate-findings,
.delfin-gate-approval,
.delfin-gate-confidence,
.delfin-gate-conflict,
.delfin-gate-cost {
    background: #ede9fe;
    border-left-color: #8b5cf6;
}
.delfin-gate-question .delfin-gate-badge,
.delfin-gate-plan-approval .delfin-gate-badge,
.delfin-gate-findings .delfin-gate-badge,
.delfin-gate-approval .delfin-gate-badge,
.delfin-gate-confidence .delfin-gate-badge,
.delfin-gate-conflict .delfin-gate-badge,
.delfin-gate-cost .delfin-gate-badge {
    background: #8b5cf6;
    color: white;
}
.delfin-gate-partial,
.delfin-gate-goal-lock,
.delfin-gate-review {
    background: #dbeafe;
    border-left-color: #2563eb;
}
.delfin-gate-partial .delfin-gate-badge,
.delfin-gate-goal-lock .delfin-gate-badge,
.delfin-gate-review .delfin-gate-badge {
    background: #2563eb;
    color: white;
}
.delfin-agent-nokey {
    padding: 20px;
    text-align: center;
    color: #6b7280;
}
.delfin-agent-nokey code {
    background: #f3f4f6;
    padding: 2px 6px;
    border-radius: 3px;
}
@keyframes delfin-pulse {
    0%, 100% { opacity: 0.4; }
    50% { opacity: 1.0; }
}
.delfin-agent-working {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 6px 12px;
    margin: 4px 0;
    border-radius: 6px;
    background: linear-gradient(90deg, #1e293b, #334155);
    color: #e2e8f0;
    font-family: 'SF Mono', 'Cascadia Code', 'Fira Code', monospace;
    font-size: 12.5px;
    font-weight: 500;
    letter-spacing: 0.02em;
    max-width: 100%;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    box-shadow: 0 1px 3px rgba(0,0,0,0.2);
    border-left: 3px solid #60a5fa;
}
.delfin-agent-working--gated {
    background: linear-gradient(90deg, #3b2a14, #4a3a1e);
    border-left-color: #f59e0b;
}
.delfin-agent-working--gated .delfin-spinner-dot { background: #f59e0b; }
.delfin-agent-working--queued {
    background: linear-gradient(90deg, #1f2937, #374151);
    border-left-color: #94a3b8;
}
.delfin-agent-working--queued .delfin-spinner-dot { background: #94a3b8; }
.delfin-agent-working--stale {
    background: linear-gradient(90deg, #3b1414, #4a1e1e);
    border-left-color: #ef4444;
}
.delfin-agent-working--stale .delfin-spinner-dot { background: #ef4444; }
.delfin-spinner {
    display: inline-flex;
    gap: 3px;
}
.delfin-spinner-dot {
    width: 5px; height: 5px;
    border-radius: 50%;
    background: #60a5fa;
    animation: delfin-bounce 1.2s ease-in-out infinite;
}
.delfin-spinner-dot:nth-child(2) { animation-delay: 0.15s; }
.delfin-spinner-dot:nth-child(3) { animation-delay: 0.3s; }
@keyframes delfin-bounce {
    0%, 80%, 100% { opacity: 0.3; transform: scale(0.8); }
    40% { opacity: 1; transform: scale(1.2); }
}
.delfin-activity-icon {
    color: #60a5fa;
    font-size: 13px;
    min-width: 16px;
    text-align: center;
}
.delfin-activity-label {
    color: #94a3b8;
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}
.delfin-activity-text {
    color: #e2e8f0;
    overflow: hidden;
    text-overflow: ellipsis;
}
.delfin-chat-thinking {
    background: #f5f3ff;
    margin: 4px 0;
    padding: 0;
    max-width: 100%;
    border-left: 3px solid #a78bfa;
    border-radius: 0 4px 4px 0;
    font-size: 11px;
    color: #6d28d9;
}
.delfin-chat-thinking summary {
    cursor: pointer;
    padding: 4px 10px;
    font-weight: 600;
    user-select: none;
}
.delfin-chat-thinking .thinking-content {
    padding: 4px 10px 8px 10px;
    color: #5b21b6;
    font-family: 'SF Mono', 'Consolas', monospace;
    font-size: 11px;
    white-space: pre-wrap;
    word-break: break-word;
    max-height: 200px;
    overflow-y: auto;
}
/* Copy button on code blocks */
.delfin-code-wrap {
    position: relative;
}
.delfin-code-wrap .delfin-copy-btn {
    position: absolute;
    top: 4px;
    right: 4px;
    background: #45475a;
    color: #cdd6f4;
    border: 1px solid #585b70;
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 10px;
    cursor: pointer;
    opacity: 0;
    transition: opacity 0.15s;
    z-index: 10;
    font-family: inherit;
}
.delfin-code-wrap:hover .delfin-copy-btn {
    opacity: 1;
}
.delfin-copy-btn.copied {
    background: #a6e3a1 !important;
    color: #1e1e2e !important;
}
/* Search highlight */
.delfin-search-hl {
    background: #fde68a;
    padding: 0 2px;
    border-radius: 2px;
}
/* Permission approval */
.delfin-chat-approval {
    background: #fef3c7;
    border: 2px solid #f59e0b;
    border-radius: 8px;
    padding: 10px 14px;
    margin: 8px 0;
    max-width: 100%;
}
.delfin-chat-approval .approval-title {
    font-weight: 700;
    color: #92400e;
    font-size: 12px;
    margin-bottom: 6px;
}
.delfin-chat-approval .approval-detail {
    font-family: 'SF Mono', 'Consolas', monospace;
    font-size: 11px;
    color: #78350f;
    white-space: pre-wrap;
    word-break: break-all;
    max-height: 120px;
    overflow-y: auto;
    margin-bottom: 8px;
}
/* Image display in chat */
.delfin-chat-agent img {
    max-width: 100%;
    max-height: 400px;
    border-radius: 6px;
    border: 1px solid #d1d5db;
    margin: 6px 0;
}
/* Agent question / interactive dialog */
.delfin-chat-question {
    background: #eff6ff;
    margin: 8px 0;
    text-align: left;
    font-size: 12px;
    max-width: 100%;
    padding: 10px 14px;
    border-left: 4px solid #3b82f6;
    border-radius: 0 8px 8px 0;
    color: #1e3a5f;
    font-weight: 500;
}
.delfin-chat-question .question-hint {
    font-size: 11px;
    color: #6b7280;
    margin-top: 6px;
    font-weight: 400;
}
.delfin-agent-question-row {
    background: #eff6ff;
    border: 1px solid #3b82f6;
    border-radius: 8px;
    flex-direction: column;
    align-items: flex-start;
    gap: 4px;
    padding: 8px 12px;
    margin: 4px 0;
}
.delfin-agent-question-row .widget-checkbox {
    margin: 0;
}
.delfin-agent-question-row .widget-checkbox label {
    font-size: 12px;
}
/* Image upload sits flush with textarea: same 80px height, content
   stacked vertically (icon over label) so nothing truncates with "...".
   Without this, FileUpload renders icon + label on one line which then
   overflows the 90px width. */
.delfin-agent-upload {
    display: flex !important;
    flex-direction: column !important;
    align-items: center !important;
    justify-content: center !important;
    padding: 4px 2px !important;
    line-height: 1.15 !important;
    font-size: 11px !important;
    white-space: normal !important;
    text-align: center !important;
    margin: 0 !important;
}
.delfin-agent-upload i.fa {
    margin: 0 0 3px 0 !important;
    font-size: 16px;
}
</style>
"""

# ---------------------------------------------------------------------------
# Minimal markdown -> HTML conversion
# ---------------------------------------------------------------------------


def _syntax_highlight(code: str, lang: str) -> str:
    """Basic syntax highlighting for code blocks.

    Uses a stash-based approach: comments and strings are extracted first
    (replaced with placeholders), then keywords/numbers are highlighted
    on the remaining text, then stashed items are restored with their
    own styling.  This prevents regex cross-contamination.
    """
    escaped = _html.escape(code)

    # -- Diff: line-based colouring, no keyword highlighting ---------------
    if lang in ("diff", "patch"):
        lines = []
        for line in escaped.split("\n"):
            if line.startswith("+"):
                lines.append(f'<span class="diff-add">{line}</span>')
            elif line.startswith("-"):
                lines.append(f'<span class="diff-del">{line}</span>')
            elif line.startswith("@@"):
                lines.append(f'<span class="diff-hdr">{line}</span>')
            else:
                lines.append(line)
        return "\n".join(lines)

    _SUPPORTED = {
        "python", "py", "javascript", "js", "typescript", "ts",
        "bash", "sh", "shell", "yaml", "yml", "json", "rust",
        "go", "java", "c", "cpp", "ruby", "rb", "toml",
    }
    if lang not in _SUPPORTED:
        return escaped

    # Resolve aliases
    _ALIAS = {
        "py": "python", "js": "javascript", "ts": "javascript",
        "sh": "bash", "shell": "bash", "yml": "yaml",
        "cpp": "c", "rb": "ruby",
    }
    actual = _ALIAS.get(lang, lang)

    # -- 1. Stash comments and strings (order matters) ---------------------
    stash: list[str] = []

    def _put(span_cls: str):
        def _repl(m):
            stash.append(f'<span class="{span_cls}">{m.group(0)}</span>')
            return f"\x01S{len(stash)-1}\x01"
        return _repl

    # Comments
    if actual in ("python", "bash", "ruby", "yaml", "toml"):
        escaped = re.sub(r"#[^\n]*", _put("cmt"), escaped)
    elif actual in ("javascript", "rust", "go", "java", "c"):
        escaped = re.sub(r"//[^\n]*", _put("cmt"), escaped)

    # Strings (HTML-escaped quotes: &quot; and &#x27;)
    escaped = re.sub(
        r'&quot;(?:[^&]|&(?!quot;))*?&quot;', _put("str"), escaped,
    )
    escaped = re.sub(
        r"&#x27;(?:[^&]|&(?!#x27;))*?&#x27;", _put("str"), escaped,
    )

    # -- 2. Highlight keywords, numbers, decorators on clean text ----------
    _KW = {
        "python": r"\b(def|class|import|from|return|if|elif|else|for|while|try|except|finally|with|as|yield|async|await|raise|pass|break|continue|and|or|not|in|is|None|True|False|self|lambda|global|nonlocal)\b",
        "javascript": r"\b(function|const|let|var|return|if|else|for|while|try|catch|finally|throw|new|class|import|export|from|async|await|yield|true|false|null|undefined|this|typeof|instanceof)\b",
        "bash": r"\b(if|then|else|elif|fi|for|while|do|done|case|esac|function|return|in|echo|exit|export|source|local)\b",
        "yaml": r"\b(true|false|null|yes|no)\b",
        "json": r"\b(true|false|null)\b",
        "rust": r"\b(fn|let|mut|pub|struct|enum|impl|trait|use|mod|if|else|for|while|loop|match|return|self|Self|true|false|None|Some|Ok|Err|async|await|move|unsafe|where|type|const|static)\b",
        "go": r"\b(func|var|const|type|struct|interface|if|else|for|range|return|package|import|go|chan|defer|select|switch|case|true|false|nil|map|make|len|append)\b",
        "java": r"\b(public|private|protected|class|interface|extends|implements|if|else|for|while|return|new|this|super|static|final|void|int|String|boolean|true|false|null|try|catch|throw|import|package)\b",
        "c": r"\b(int|char|float|double|void|if|else|for|while|return|struct|typedef|enum|switch|case|break|continue|sizeof|const|static|unsigned|long|short|NULL|true|false|include|define)\b",
        "ruby": r"\b(def|class|module|if|else|elsif|unless|while|until|for|do|end|return|yield|self|nil|true|false|require|include|attr_accessor|puts|print)\b",
        "toml": r"\b(true|false)\b",
    }

    kw_pattern = _KW.get(actual)
    if kw_pattern:
        escaped = re.sub(kw_pattern, r'<span class="kw">\1</span>', escaped)

    # Numbers (only bare numbers, not inside stash placeholders)
    escaped = re.sub(r"(?<!\x01S)\b(\d+\.?\d*)\b", r'<span class="num">\1</span>', escaped)

    # Decorators (Python)
    if actual == "python":
        escaped = re.sub(r"(@\w+)", r'<span class="dec">\1</span>', escaped)

    # -- 3. Restore stashed items ------------------------------------------
    for i, span in enumerate(stash):
        escaped = escaped.replace(f"\x01S{i}\x01", span)

    return escaped


def _md_to_html(text: str) -> str:
    """Convert markdown subset to HTML for chat display."""
    # Protect code blocks first (extract, replace later)
    code_blocks: list[str] = []

    def _stash_code_block(m):
        lang = m.group(1).strip().lower()
        code = m.group(2).strip()
        highlighted = _syntax_highlight(code, lang)
        lang_label = f'<span class="code-lang">{_html.escape(lang)}</span>' if lang else ""
        # data-code stores raw text for copy button
        raw_escaped = _html.escape(code).replace('"', "&quot;")
        idx = len(code_blocks)
        code_blocks.append(
            f'<div class="delfin-code-wrap">'
            f'<button class="delfin-copy-btn" data-codeidx="{idx}" '
            f'onclick="(function(b){{var p=b.closest(&quot;.delfin-code-wrap&quot;);'
            f"var c=p.querySelector('code').textContent;"
            f"navigator.clipboard.writeText(c).then(function(){{"
            f"b.textContent='Copied!';b.classList.add('copied');"
            f"setTimeout(function(){{b.textContent='Copy';b.classList.remove('copied');}},1500);"
            f'}})}})(this)">Copy</button>'
            f"<pre>{lang_label}<code>{highlighted}</code></pre></div>"
        )
        return f"\x00CB{len(code_blocks) - 1}\x00"

    text = re.sub(r"```(\w*)\n?(.*?)```", _stash_code_block, text, flags=re.DOTALL)

    escaped = _html.escape(text)

    # Inline code
    escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
    # Bold
    escaped = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", escaped)
    # Italic
    escaped = re.sub(r"\*(.+?)\*", r"<i>\1</i>", escaped)
    # Headers (## and ###)
    escaped = re.sub(
        r"^(#{1,3})\s+(.+)$",
        lambda m: f"<b style='font-size:{15 - len(m.group(1))}px'>{m.group(2)}</b>",
        escaped,
        flags=re.MULTILINE,
    )
    # Horizontal rules
    escaped = re.sub(r"^-{3,}$", "<hr style='margin:6px 0;border:none;border-top:1px solid #d1d5db'>", escaped, flags=re.MULTILINE)
    # List items (- or *)
    escaped = re.sub(r"^[\-\*]\s+(.+)$", r"&bull; \1", escaped, flags=re.MULTILINE)
    # Numbered list items
    escaped = re.sub(r"^(\d+)\.\s+(.+)$", r"\1. \2", escaped, flags=re.MULTILINE)
    # Table rows: | a | b | → simple rendering
    escaped = re.sub(
        r"^\|[\s\-:|]+\|$", "", escaped, flags=re.MULTILINE,
    )  # remove separator rows
    escaped = re.sub(
        r"^\|(.+)\|$",
        lambda m: "<code>" + " | ".join(
            c.strip() for c in m.group(1).split("|")
        ) + "</code>",
        escaped,
        flags=re.MULTILINE,
    )
    # Newlines to <br>
    escaped = escaped.replace("\n", "<br>")
    # Restore code blocks
    for i, block in enumerate(code_blocks):
        escaped = escaped.replace(f"\x00CB{i}\x00", block)
    return escaped


# ---------------------------------------------------------------------------
# Slash-command catalog (used by the discoverable command palette)
# ---------------------------------------------------------------------------

# Each tuple: (category, command, summary, has_args)
# - command: literal prefix to insert into the input (no trailing space → user
#   completes; trailing space → ready to type args)
# - has_args: True means the command takes arguments — palette adds a trailing
#   space so the user starts typing immediately
_SLASH_COMMANDS: tuple[tuple[str, str, str, bool], ...] = (
    # Session
    ("Session", "/help", "Show this help", False),
    ("Session", "/clear", "Clear chat history", False),
    ("Session", "/cost", "Show token usage & cost", False),
    ("Session", "/usage", "Detailed token usage & session stats", False),
    ("Session", "/status", "Show engine status", False),
    ("Session", "/compact", "Summarize context to reduce tokens", False),
    ("Session", "/context", "Show context-window usage + compaction status", False),
    ("Session", "/agents", "List subagent presets (or /agents stats for telemetry)", False),
    ("Session", "/skills", "List discovered skills (or /skills <name> for body)", False),
    ("Session", "/undo", "Undo last agent turn (drop from context)", False),
    ("Session", "/retry", "Regenerate last response", False),
    ("Session", "/stop", "Stop current generation", False),
    ("Session", "/reset", "Reset engine for new cycle", False),
    ("Session", "/export", "Export chat as Markdown", False),
    ("Session", "/search", "Search in chat history", True),
    # Engine config
    ("Engine", "/provider", "Switch provider (claude/openai/kit)", True),
    ("Engine", "/model", "Switch model (depends on provider)", True),
    ("Engine", "/effort", "Set effort (low/medium/high/xhigh)", True),
    ("Engine", "/mode", "Switch mode (dashboard/solo/quick/reviewed/tdd/cluster/full)", True),
    ("Engine", "/perms", "Show or set permission profile", True),
    ("Engine", "/perm-cycle", "Cycle through permission profiles", False),
    # Git
    ("Git", "/git status", "Show git status", False),
    ("Git", "/git diff", "Show staged/unstaged changes", False),
    ("Git", "/git log", "Show recent commits", False),
    ("Git", "/git branch", "Show branches", False),
    # Memory
    ("Memory", "/remember", "Save a typed memory ([user|feedback|project|reference:] <text>)", True),
    ("Memory", "/memories", "List all memories", False),
    ("Memory", "/memories verify", "Check stored memories for stale file refs", False),
    ("Memory", "/forget", "Delete a memory by index", True),
    ("Memory", "/plans", "List saved Plan-Mode plans (or /plans <name>)", False),
    ("Plan", "/plan approve", "Approve a pending plan when model forgot ExitPlanMode", False),
    ("Plan", "/plan reject", "Reject a pending plan and exit plan mode", False),
    ("Hooks", "/hooks", "List/add/remove/dry-run settings.json hooks", False),
    ("Session", "/session", "ls/restore/search/fork/tree/handoff/bundle/import/archive", False),
    ("MCP", "/mcp", "List/add/remove/toggle MCP servers (~/.delfin/mcp_servers.json)", False),
    ("Commands", "/commands", "List user-defined slash commands from ~/.delfin/commands/", False),
    ("Project", "/init", "Scan repo + write AGENTS.md / .delfin/settings.json scaffold", False),
    ("Bash", "/bash", "List/watch/kill bash_background jobs", False),
    ("Learning", "/failures", "Recurring tool-error patterns from the failure log", False),
    # Workspace (agent_workspace/)
    ("Workspace", "/workspace ls", "List files in agent workspace", False),
    ("Workspace", "/workspace read", "Read a workspace file", True),
    ("Workspace", "/workspace clean", "Remove all workspace files", False),
    # Dashboard navigation
    ("Dashboard", "/tab", "Switch tab (submit/orca/jobs/calc/settings)", True),
    ("Dashboard", "/jobs", "Switch to Job Status tab", False),
    ("Dashboard", "/jobs check", "Check for new job-state events now", False),
    ("Dashboard", "/ui list", "List all controllable widgets", False),
    ("Dashboard", "/ui", "Manipulate a widget (show/click/value/...)", True),
    # CONTROL
    ("CONTROL", "/control show", "Show CONTROL content from Submit tab", False),
    ("CONTROL", "/control set", "Replace entire CONTROL content", True),
    ("CONTROL", "/control key", "Set one CONTROL key (e.g. functional BP86)", True),
    ("CONTROL", "/control validate", "Validate CONTROL syntax", False),
    # ORCA Builder
    ("ORCA", "/orca show", "Show ORCA Builder settings", False),
    ("ORCA", "/orca set", "Set an ORCA Builder param", True),
    ("ORCA", "/orca submit", "Submit ORCA job", False),
    # Submit
    ("Submit", "/submit", "Submit job from current Submit tab state", False),
    # Calculations / archive browsing
    ("Calc", "/calc ls", "List directories/files", True),
    ("Calc", "/calc cd", "Navigate calc folder (syncs browser)", True),
    ("Calc", "/calc select", "Select file in browser", True),
    ("Calc", "/calc read", "Read a calc file", True),
    ("Calc", "/calc tail", "Read last 8KB (convergence check)", True),
    ("Calc", "/calc info", "Show folder summary & status", True),
    ("Calc", "/calc tree", "Show directory tree", True),
    ("Calc", "/calc search", "Search files by glob pattern", True),
    # Analysis
    ("Analysis", "/analyze", "Full analysis (energy + convergence + errors)", True),
    ("Analysis", "/analyze energy", "Extract Gibbs/ZPE/electronic energies", True),
    ("Analysis", "/analyze convergence", "Check SCF convergence", True),
    ("Analysis", "/analyze errors", "Scan for ORCA error patterns", True),
    ("Analysis", "/analyze status", "Overview of all calc folders", False),
    # Recalc / cancel (destructive — confirmations enforced)
    ("Recalc", "/recalc check", "Check if recalc needed (safe)", True),
    ("Recalc", "/recalc check-all", "Scan all folders (safe)", False),
    ("Recalc", "/recalc", "Submit recalc (confirms first)", True),
    ("Recalc", "/recalc auto", "Recalc all that need it (confirms first)", False),
    ("Recalc", "/cancel", "Cancel a job (confirms first)", True),
    ("Recalc", "/cancel all", "Cancel all active jobs (confirms first)", False),
    ("Recalc", "/cancel running", "Cancel only running (R) jobs", False),
    ("Recalc", "/cancel pending", "Cancel only pending (PD) jobs", False),
    # Batch
    ("Batch", "/batch from-calc", "Build batch from calc folders (optional glob)", True),
    ("Batch", "/batch add", "Add one entry (Name;SMILES;...)", True),
    ("Batch", "/batch show", "Show current batch content", False),
    ("Batch", "/batch clear", "Clear batch field", False),
)


def _filter_slash_commands(
    commands: tuple[tuple[str, str, str, bool], ...],
    query: str,
) -> list[tuple[str, str, str, bool]]:
    """Case-insensitive filter on command + summary + category.

    An empty query returns everything in catalog order.
    """
    q = (query or "").strip().lower()
    if not q:
        return list(commands)
    out: list[tuple[str, str, str, bool]] = []
    for entry in commands:
        category, cmd, summary, _ = entry
        haystack = f"{cmd} {summary} {category}".lower()
        if q in haystack:
            out.append(entry)
    return out


def _render_slash_palette_html(
    commands: list[tuple[str, str, str, bool]],
    query: str = "",
) -> str:
    """Render the slash-command palette as grouped HTML.

    The output is wrapped in ``<div class="delfin-slash-palette">`` and each
    command row carries a ``data-command`` attribute so a click handler can
    insert it into the input.  Inline styling keeps the palette portable.
    """
    if not commands:
        q = _html.escape(query)
        return (
            '<div class="delfin-slash-palette" '
            'style="padding:10px;background:#fff;border:1px solid #e5e7eb;'
            'border-radius:6px;color:#6b7280;font-size:12px;">'
            f'No commands match <code>{q}</code>.</div>'
        )

    grouped: dict[str, list[tuple[str, str, str, bool]]] = {}
    for entry in commands:
        grouped.setdefault(entry[0], []).append(entry)

    sections: list[str] = []
    for category, entries in grouped.items():
        rows: list[str] = []
        for cat, cmd, summary, has_args in entries:
            insert = cmd + (" " if has_args else "")
            rows.append(
                '<div class="slash-row" '
                f'data-command="{_html.escape(insert)}" '
                'style="display:flex;gap:10px;padding:5px 8px;cursor:pointer;'
                'border-radius:4px;">'
                f'<code style="color:#3b82f6;font-size:12px;font-weight:600;'
                'min-width:160px;">'
                f'{_html.escape(cmd)}</code>'
                f'<span style="color:#374151;font-size:12px;">'
                f'{_html.escape(summary)}</span>'
                '</div>'
            )
        sections.append(
            f'<div style="margin-bottom:8px;">'
            f'<div style="font-size:10px;font-weight:700;color:#6b7280;'
            f'text-transform:uppercase;letter-spacing:0.4px;'
            f'padding:3px 8px;">{_html.escape(category)}</div>'
            + "".join(rows) +
            '</div>'
        )

    return (
        '<div class="delfin-slash-palette" '
        'style="max-height:340px;overflow-y:auto;background:#f9fafb;'
        'border:1px solid #e5e7eb;border-radius:6px;padding:6px;'
        'font-family:-apple-system,BlinkMacSystemFont,sans-serif;">'
        + "".join(sections) +
        '</div>'
    )


# ---------------------------------------------------------------------------
# D1 — confirmation extraction helpers
# ---------------------------------------------------------------------------

# German + English phrases that signal the agent is asking the user to
# confirm before it proceeds.  Used to surface Approve/Deny buttons
# instead of waiting for a free-text "yes".
_CONFIRMATION_PATTERNS: tuple[str, ...] = (
    "soll ich", "soll das", "darf ich",
    "should i", "shall i",
    "okay so?", "okay so",
    "willst du", "möchtest du",
    "weitermachen?", "ausführen?", "submitten?", "submit?",
    "proceed?", "confirm?", "go ahead?",
    "shall we", "do you want",
)


def _extract_action_commands(agent_text: str) -> list[str]:
    """Pull every ``ACTION: /command`` line out of an agent response.

    Returns the bare slash-commands (no ``ACTION:`` prefix), in original
    order.  Only single-line ACTIONs are recognised (matches the existing
    dispatcher in ``_dashboard_auto_exec``).
    """
    if not agent_text:
        return []
    out: list[str] = []
    for raw in agent_text.splitlines():
        line = raw.rstrip()
        if line.startswith("ACTION:"):
            cmd = line[len("ACTION:"):].strip()
            if cmd.startswith("/"):
                out.append(cmd)
    return out


# ---------------------------------------------------------------------------
# D2 — context-aware suggestion when the user switches dashboard tabs
# ---------------------------------------------------------------------------

# Tab-name → one-shot suggestion the agent surfaces when the user lands
# on that tab.  Keys match the titles registered in
# ``delfin/dashboard/__init__.py``'s tab_specs.  Empty value means
# "never suggest anything for this tab".
_TAB_SUGGESTIONS: dict[str, str] = {
    "Submit Job":
        "Soll ich die aktuellen CONTROL-Werte validieren oder einen Submit "
        "vorbereiten?",
    "Recalc":
        "Soll ich `/recalc check-all` laufen lassen, um Folder mit "
        "fehlgeschlagenen oder unvollständigen Jobs zu finden?",
    "Job Status":
        "Soll ich die letzten Job-Events prüfen (`/jobs check`) und "
        "fehlgeschlagene Jobs analysieren?",
    "Calculations":
        "Soll ich eine Energie-Tabelle für alle Folder hier erzeugen "
        "(`/skill energy-table`)?",
    "Literature":
        "Soll ich nach einem Thema in der DELFIN-Literatur (search_docs) "
        "suchen — z. B. Funktional/Basis für dein aktuelles System?",
    "ORCA Builder":
        "Soll ich aus dem aktuellen Builder-State eine ORCA-Input-Datei "
        "generieren oder die Parameter prüfen?",
    "Agent Activity": "",
    "DELFIN Agent": "",
    "Archive": "",
    "Remote Archive": "",
    "Settings": "",
}


def _suggestion_for_tab(tab_name: str) -> str | None:
    """Return the one-shot suggestion for ``tab_name`` or ``None``.

    Unknown tabs and tabs explicitly mapped to "" return ``None`` so
    the caller knows there's nothing to show.
    """
    if not tab_name:
        return None
    suggestion = _TAB_SUGGESTIONS.get(tab_name)
    if suggestion:
        return suggestion
    return None


# ---------------------------------------------------------------------------
# Helper: distinguish "interactive approval needed" from "structurally
# blocked" denials so we can show the right message instead of an
# Approve/Deny prompt that won't actually unblock anything.
# ---------------------------------------------------------------------------

def _extract_denied_tool_name(denied_raw: str) -> str:
    """Pull ``tool_name`` from a permission-denied payload.

    The payload is a Python-literal dict produced by the CLI; we accept
    string and dict inputs, return ``""`` if no tool name can be read.
    """
    if not denied_raw:
        return ""
    if isinstance(denied_raw, dict):
        return str(denied_raw.get("tool_name", "") or "")
    import ast as _ast
    try:
        d = _ast.literal_eval(denied_raw) if str(denied_raw).strip().startswith("{") else {}
    except Exception:
        d = {}
    if isinstance(d, dict):
        return str(d.get("tool_name", "") or "")
    return ""


def _tool_in_allowlist(tool_name: str, cli_tools: list[str] | None) -> bool:
    """Return True if ``tool_name`` is permitted by the CLI tool allow-list.

    A ``None`` allow-list means no restriction (CLI default = all tools
    available). An entry like ``"Bash(git *)"`` counts as enabling the
    ``Bash`` tool — patterns are interactive guards on top, not extra
    blocks. Empty list explicitly disallows everything.
    """
    if not tool_name:
        return False
    if cli_tools is None:
        return True
    base_names = {entry.split("(", 1)[0].strip() for entry in cli_tools}
    return tool_name in base_names


def _is_structurally_blocked(denied_raw: str, cli_tools: list[str] | None) -> bool:
    """``True`` when the denied tool is missing from the CLI allow-list.

    In that case Approve/Deny buttons can't help — the CLI subprocess
    rejects the call before any interactive prompt fires. The caller
    should surface a mode-switch hint instead of the normal prompt.
    """
    if cli_tools is None:
        return False
    tool = _extract_denied_tool_name(denied_raw)
    if not tool:
        return False
    return not _tool_in_allowlist(tool, cli_tools)


def _should_show_action_confirmation(agent_text: str) -> bool:
    """Decision rule for the D1 confirm-buttons: agent proposed
    ``ACTION:`` lines AND its prose explicitly invites confirmation.

    Returns False if either signal is missing — auto-exec then runs
    normally and the user types "ja"/"nein" as before.
    """
    if not agent_text:
        return False
    if not _extract_action_commands(agent_text):
        return False
    return _text_requests_confirmation(agent_text)


def _render_action_confirmation_html(commands: list[str]) -> str:
    """Render the proposed ACTION list as a compact preview."""
    if not commands:
        return ""
    rows = "".join(
        f'<li style="font-family:monospace;font-size:11px;color:#92400e;'
        f'padding:2px 0;">{_html.escape(cmd[:160])}</li>'
        for cmd in commands
    )
    return (
        '<div style="flex:1;min-width:0;">'
        '<div style="font-size:11px;font-weight:700;color:#92400e;'
        'margin-bottom:4px;">'
        f'Soll ich diese {len(commands)} Aktion(en) ausführen?</div>'
        '<ul style="list-style:disc;margin:0;padding-left:20px;">'
        + rows + '</ul></div>'
    )


def _text_requests_confirmation(agent_text: str) -> bool:
    """Heuristic: does the agent's text invite explicit user confirmation?

    Trigger phrases match common German + English question patterns that
    follow proposed actions.  Used to decide whether to show
    Approve/Deny buttons instead of streaming auto-execution.
    """
    if not agent_text:
        return False
    haystack = agent_text.lower()
    return any(p in haystack for p in _CONFIRMATION_PATTERNS)


def _render_molecule_to_png_b64(
    mol_or_smiles, *, size: tuple[int, int] = (320, 240)
) -> str | None:
    """Render an RDKit Mol or SMILES string to a base64 PNG.

    Returns ``None`` if RDKit is unavailable or the molecule is invalid.
    Used by the artifact renderer to draw 2D structures inline in the chat.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
    except ImportError:
        return None
    try:
        if isinstance(mol_or_smiles, str):
            mol = Chem.MolFromSmiles(mol_or_smiles)
        else:
            mol = mol_or_smiles
        if mol is None:
            return None
        # Compute 2D coords if missing (XYZ-derived mols often lack them)
        if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D() is False:
            try:
                AllChem.Compute2DCoords(mol)
            except Exception:
                pass
        img = Draw.MolToImage(mol, size=size)
        import base64, io
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        return None


def _render_xyz_summary(path) -> tuple[int, str] | None:
    """Parse an XYZ file and return ``(atom_count, summary_string)``.

    Summary is a compact "C12 H10 O2"-style formula plus first three
    atom lines.  Returns None for unreadable / non-XYZ files.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None
    try:
        n_atoms = int(lines[0].strip().split()[0])
    except (ValueError, IndexError):
        return None
    # Skip comment line (second line) and parse atoms
    atom_lines = lines[2:2 + n_atoms]
    if not atom_lines:
        return None
    counts: dict[str, int] = {}
    for raw in atom_lines:
        parts = raw.strip().split()
        if not parts:
            continue
        sym = parts[0]
        counts[sym] = counts.get(sym, 0) + 1
    formula = " ".join(
        f"{sym}{n}" if n > 1 else sym
        for sym, n in sorted(counts.items())
    )
    return n_atoms, formula


def _render_artifact_inline(path) -> str | None:
    """Render a workspace artifact as inline HTML for the chat.

    Supported types:
        - PNG / JPG / GIF: base64 ``<img>`` (Voila-safe — no file:// link)
        - SVG: inline SVG (kept under 200kB to avoid blowing up the DOM)
        - CSV / TSV: first 12 rows as a compact table
        - JSON: pretty-printed up to 6kB inside ``<pre>``

    Returns ``None`` for unsupported types or unreadable files so the
    caller can simply skip them.
    """
    from pathlib import Path

    p = Path(path)
    try:
        if not p.is_file():
            return None
    except OSError:
        return None

    suffix = p.suffix.lower()
    name_html = _html.escape(p.name)

    # ---- Images via base64 ---------------------------------------------
    if suffix in (".png", ".jpg", ".jpeg", ".gif"):
        try:
            data = p.read_bytes()
        except OSError:
            return None
        # Cap raster images at ~2 MB to keep the chat snappy
        if len(data) > 2_000_000:
            return (
                f'<div class="delfin-artifact" style="font-size:11px;'
                f'color:#9ca3af;padding:6px 0;">'
                f'📷 <code>{name_html}</code> ({len(data) // 1024} KB) — '
                f'too large to embed inline.</div>'
            )
        import base64
        b64 = base64.b64encode(data).decode("ascii")
        mime = "image/jpeg" if suffix in (".jpg", ".jpeg") else f"image/{suffix[1:]}"
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'📷 <code>{name_html}</code></div>'
            f'<img src="data:{mime};base64,{b64}" '
            f'style="max-width:100%;border-radius:4px;display:block;" '
            f'alt="{name_html}"/></div>'
        )

    # ---- Inline SVG ----------------------------------------------------
    if suffix == ".svg":
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None
        if len(text) > 200_000:
            return (
                f'<div class="delfin-artifact" style="font-size:11px;'
                f'color:#9ca3af;padding:6px 0;">'
                f'🎨 <code>{name_html}</code> ({len(text) // 1024} KB) — '
                f'too large to embed inline.</div>'
            )
        # Strip a leading XML declaration for clean inline embedding
        if text.lstrip().startswith("<?xml"):
            text = text.split("?>", 1)[-1]
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'🎨 <code>{name_html}</code></div>'
            f'<div style="max-width:100%;overflow:auto;">{text}</div></div>'
        )

    # ---- CSV / TSV preview ---------------------------------------------
    if suffix in (".csv", ".tsv"):
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None
        sep = "\t" if suffix == ".tsv" else ","
        all_lines = [ln for ln in text.splitlines() if ln.strip()]
        if not all_lines:
            return None
        rows = [ln.split(sep) for ln in all_lines[:12]]
        max_cols = max(len(r) for r in rows)
        header = rows[0]
        body = rows[1:]
        head_html = "".join(
            f'<th style="background:#f3f4f6;padding:4px 8px;'
            f'text-align:left;font-size:11px;color:#374151;'
            f'border-bottom:1px solid #e5e7eb;">{_html.escape(c)}</th>'
            for c in header
        )
        body_html_rows: list[str] = []
        for row in body:
            cells = "".join(
                f'<td style="padding:4px 8px;font-size:11px;color:#1f2937;'
                f'border-bottom:1px solid #f3f4f6;">{_html.escape(c)}</td>'
                for c in row
            )
            body_html_rows.append(f'<tr>{cells}</tr>')
        more = ""
        if len(all_lines) > 12:
            more = (
                f'<div style="font-size:10px;color:#9ca3af;'
                f'padding:4px 0 0 0;">'
                f'… {len(all_lines) - 12} more rows</div>'
            )
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'📊 <code>{name_html}</code> '
            f'({max_cols} cols, {len(all_lines)} rows)</div>'
            f'<div style="overflow-x:auto;"><table style="border-collapse:collapse;'
            f'font-family:-apple-system,BlinkMacSystemFont,sans-serif;">'
            f'<thead><tr>{head_html}</tr></thead>'
            f'<tbody>{"".join(body_html_rows)}</tbody></table></div>'
            f'{more}</div>'
        )

    # ---- SMILES file (single SMILES per line) --------------------------
    if suffix == ".smi":
        try:
            text = p.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            return None
        if not text:
            return None
        first = text.splitlines()[0].split()[0]
        b64 = _render_molecule_to_png_b64(first)
        if b64 is None:
            return (
                f'<div class="delfin-artifact" style="font-size:11px;'
                f'color:#9ca3af;padding:6px 0;">'
                f'🧪 <code>{name_html}</code> SMILES: '
                f'<code>{_html.escape(first[:120])}</code></div>'
            )
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'🧪 <code>{name_html}</code> '
            f'<code style="color:#374151;">{_html.escape(first[:120])}</code></div>'
            f'<img src="data:image/png;base64,{b64}" '
            f'style="display:block;border-radius:4px;" '
            f'alt="{name_html}"/></div>'
        )

    # ---- 2D structure files (.mol, .sdf) --------------------------------
    if suffix in (".mol", ".sdf"):
        try:
            from rdkit import Chem
        except ImportError:
            return None
        try:
            if suffix == ".sdf":
                supplier = Chem.SDMolSupplier(str(p))
                mol = next((m for m in supplier if m is not None), None)
            else:
                mol = Chem.MolFromMolFile(str(p))
        except Exception:
            mol = None
        if mol is None:
            return None
        b64 = _render_molecule_to_png_b64(mol)
        if b64 is None:
            return None
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'🧪 <code>{name_html}</code> '
            f'({mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds)</div>'
            f'<img src="data:image/png;base64,{b64}" '
            f'style="display:block;border-radius:4px;" '
            f'alt="{name_html}"/></div>'
        )

    # ---- XYZ coordinates: formula summary + first atoms -----------------
    if suffix == ".xyz":
        summary = _render_xyz_summary(p)
        if summary is None:
            return None
        n_atoms, formula = summary
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'⚛️ <code>{name_html}</code> '
            f'({n_atoms} atoms)</div>'
            f'<div style="font-size:13px;color:#1f2937;font-family:monospace;'
            f'padding:4px 0;">{_html.escape(formula)}</div></div>'
        )

    # ---- JSON pretty-print ---------------------------------------------
    if suffix == ".json":
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None
        if len(text) > 6_000:
            return (
                f'<div class="delfin-artifact" style="font-size:11px;'
                f'color:#9ca3af;padding:6px 0;">'
                f'📄 <code>{name_html}</code> ({len(text) // 1024} KB JSON) — '
                f'too large to inline; open the file directly.</div>'
            )
        return (
            f'<div class="delfin-artifact" style="margin:6px 0;padding:6px;'
            f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
            f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
            f'📄 <code>{name_html}</code></div>'
            f'<pre style="margin:0;padding:6px;background:#f9fafb;'
            f'border-radius:4px;font-size:11px;max-height:240px;'
            f'overflow:auto;">{_html.escape(text)}</pre></div>'
        )

    return None


def _build_full_transcript_handoff(
    chat_messages: list[dict],
    old_mode: str,
    new_mode: str,
    *,
    max_chars_per_msg: int = 4000,
) -> str:
    """Format the UI conversation history as a single handoff message.

    Used when the user switches modes mid-session — the new engine
    starts with no message history, so we replay the conversation as
    one user-side block at the top of the next turn so the model has
    the full context without re-asking. Per-message text is capped at
    ``max_chars_per_msg`` to keep huge runs from blowing the context.

    Returns an empty string when there is nothing to hand off (no
    chat_messages, or only system banners).

    Roles understood: ``user``, ``assistant``, ``system``, ``tool``.
    Empty content is skipped. ``role`` keys missing or unknown are
    rendered under a ``Note`` heading so nothing silently disappears.
    """
    if not chat_messages:
        return ""

    def _truncate(text: str) -> str:
        if not isinstance(text, str):
            text = str(text or "")
        if len(text) <= max_chars_per_msg:
            return text
        head = max_chars_per_msg - 64
        return text[:head] + f"\n... [truncated {len(text) - head} chars]"

    rendered: list[str] = []
    for msg in chat_messages:
        role = (msg.get("role") or "").strip().lower()
        content = msg.get("content")
        if not content:
            continue
        body = _truncate(content)
        if role == "user":
            heading = "User"
        elif role == "assistant":
            heading = "Assistant"
        elif role == "system":
            heading = "System"
        elif role == "tool":
            tool_name = msg.get("tool_name") or "tool"
            heading = f"Tool ({tool_name})"
        else:
            heading = f"Note ({role or 'unknown'})"
        rendered.append(f"### {heading}\n{body}")

    if not rendered:
        return ""

    header = (
        f"[Mode switch: {old_mode or '?'} → {new_mode}]\n\n"
        f"You are now operating as the **{new_mode}** agent. The user "
        f"kept the same conversation; below is the full prior transcript "
        f"so you can continue without re-asking for context. Pick up "
        f"where the previous mode left off."
    )
    return "\n\n".join(
        [header, "--- Prior conversation (full transcript) ---", *rendered,
         "--- End of prior context ---"]
    )


def _format_session_boot(
    *,
    outcomes: list | None = None,
    jobs: list | None = None,
    commits: list | None = None,
    calc_path: str = "",
    branch: str = "",
    max_chars: int = 2400,
) -> str:
    """Render a one-shot "session boot" snippet for the dashboard agent.

    Goal: on the first user message of a fresh session, give the agent
    enough domain context to answer "what's going on?" without firing a
    dozen tool calls. Cap is tight (~600 tokens) — this is a *primer*,
    not a replacement for live state.

    All inputs are optional: missing data → that block is skipped.
    Empty everything → returns "" so the caller can omit the block.

    ``outcomes`` items may be ``CycleOutcome`` dataclasses or plain
    dicts (so the helper stays mock-friendly in tests).
    """
    parts: list[str] = []

    def _attr(o, name, default=""):
        if isinstance(o, dict):
            return o.get(name, default)
        return getattr(o, name, default)

    if outcomes:
        recent = outcomes[-5:]
        lines = ["Recent outcomes (last 5 cycles):"]
        for o in recent:
            verdict = _attr(o, "verdict", "?") or "?"
            mode = _attr(o, "mode", "")
            task = _attr(o, "task", "") or ""
            cost = _attr(o, "cost_usd", 0.0) or 0.0
            task_short = str(task).replace("\n", " ")[:60]
            lines.append(
                f"  [{verdict}] {mode}: {task_short!r} (${float(cost):.3f})"
            )
        parts.append("\n".join(lines))

    if jobs:
        # Group by status for a compact summary
        from collections import Counter
        statuses = Counter()
        named: list[str] = []
        for j in jobs:
            status = (_attr(j, "status", "") or "").upper()
            statuses[status] += 1
            name = _attr(j, "name", "") or _attr(j, "job_id", "")
            if name and len(named) < 8:
                named.append(f"{name}({status[:1] or '?'})")
        summary = ", ".join(f"{k}:{v}" for k, v in statuses.most_common())
        line = f"Active SLURM jobs: {summary}"
        if named:
            line += f" — {', '.join(named)}"
        parts.append(line)

    if commits:
        head = "Recent commits:"
        body = "\n".join(f"  {c}" for c in commits[:3])
        parts.append(f"{head}\n{body}")

    if branch:
        parts.append(f"Active branch: {branch}")

    if calc_path:
        parts.append(f"Active calc folder: {calc_path}")

    if not parts:
        return ""

    body = "\n\n".join(parts)
    if len(body) > max_chars:
        body = body[: max_chars - 32].rstrip() + "\n... [truncated]"
    return f"[Session boot context]\n{body}"


def _format_solo_domain_state(snapshot: dict) -> str:
    """Render a compact ``[Domain state]`` block for the solo-mode prompt.

    Empty values are skipped.  The block is intentionally short — solo
    mode focuses on code, not UI; this is just the *background* state
    so the agent doesn't have to ask "what folder are you in?".

    ``snapshot`` keys (all optional):
        - ``calc_dir`` (str): root of the active calculations folder
        - ``selected`` (str): currently selected file/folder, relative to calc_dir
        - ``control`` (dict): {functional, basis, pal, charge, mult, ...}
        - ``orca_builder`` (dict): {method, basis, charge, mult, ...}
        - ``job_summary`` (str): one-line job overview ("2 RUNNING, 5 PENDING")
        - ``workspace_files`` (list[str]): file names in agent_workspace/
        - ``active_tab`` (str): currently open dashboard tab
        - ``perm_profile`` (str): active permission profile
    """
    if not snapshot:
        return ""

    lines: list[str] = []

    calc_dir = (snapshot.get("calc_dir") or "").strip()
    if calc_dir:
        lines.append(f"calc_dir: {calc_dir}")

    selected = (snapshot.get("selected") or "").strip()
    if selected:
        lines.append(f"selected: {selected}")

    control = snapshot.get("control") or {}
    if isinstance(control, dict):
        ctl_parts: list[str] = []
        functional = control.get("functional")
        basis = control.get("main_basisset") or control.get("basis")
        if functional and basis:
            ctl_parts.append(f"{functional}/{basis}")
        elif functional:
            ctl_parts.append(str(functional))
        for k in ("PAL", "charge", "multiplicity", "solvent"):
            v = control.get(k)
            if v not in (None, "", 0):
                ctl_parts.append(f"{k}={v}")
        if ctl_parts:
            lines.append("control: " + ", ".join(ctl_parts))

    builder = snapshot.get("orca_builder") or {}
    if isinstance(builder, dict):
        b_parts: list[str] = []
        m = builder.get("method")
        b = builder.get("basis")
        if m and b:
            b_parts.append(f"{m}/{b}")
        elif m:
            b_parts.append(str(m))
        for k in ("charge", "mult", "pal"):
            v = builder.get(k)
            if v not in (None, "", 0):
                b_parts.append(f"{k}={v}")
        if b_parts:
            lines.append("orca_builder: " + ", ".join(b_parts))

    job_summary = (snapshot.get("job_summary") or "").strip()
    if job_summary:
        lines.append(f"jobs: {job_summary}")

    ws = snapshot.get("workspace_files") or []
    if isinstance(ws, list) and ws:
        preview = ", ".join(str(f) for f in ws[:8])
        tail = f" (+{len(ws) - 8} more)" if len(ws) > 8 else ""
        lines.append(f"workspace: {preview}{tail}")

    active_tab = (snapshot.get("active_tab") or "").strip()
    if active_tab:
        lines.append(f"active_tab: {active_tab}")

    perm = (snapshot.get("perm_profile") or "").strip()
    if perm:
        lines.append(f"perms: {perm}")

    if not lines:
        return ""
    return "--- Domain State ---\n" + "\n".join(lines)


def _render_subagent_pane_html(calls: list[dict]) -> str:
    """Render the active/completed subagent calls as a compact panel.

    Each call dict has: ``subagent_type`` (str), ``description`` (str),
    ``prompt`` (str), ``status`` ("running" | "done"), ``output`` (str).
    Returns empty string if there are no calls — caller hides the widget.
    """
    if not calls:
        return ""
    n_running = sum(1 for c in calls if c.get("status") == "running")
    rows: list[str] = []
    for call in calls:
        status = call.get("status", "running")
        if status == "running":
            glyph, color = "↻", "#3b82f6"
        else:
            glyph, color = "✓", "#10b981"
        sub_type = _html.escape(str(call.get("subagent_type") or "agent"))
        desc = _html.escape(str(call.get("description") or "")[:80])
        prompt_preview = _html.escape(
            str(call.get("prompt") or "")[:200].replace("\n", " ").strip()
        )
        output = call.get("output") or ""
        details = ""
        if output:
            output_preview = _html.escape(str(output)[:600])
            details = (
                '<details style="margin-top:4px;">'
                '<summary style="cursor:pointer;font-size:11px;color:#6b7280;">'
                f'Result ({len(str(output))} chars)</summary>'
                '<pre style="margin:4px 0 0 0;padding:6px 8px;'
                'background:#f3f4f6;border-radius:4px;font-size:11px;'
                'white-space:pre-wrap;max-height:160px;overflow-y:auto;">'
                f'{output_preview}</pre></details>'
            )
        rows.append(
            f'<li style="padding:6px 8px;border-bottom:1px solid #f3f4f6;">'
            f'<div style="display:flex;gap:8px;align-items:center;">'
            f'<span style="color:{color};font-weight:700;">{glyph}</span>'
            f'<span style="background:#dbeafe;color:#1e40af;font-size:10px;'
            f'font-weight:700;padding:1px 6px;border-radius:8px;'
            f'text-transform:uppercase;">{sub_type}</span>'
            f'<span style="color:#1f2937;font-size:12px;font-weight:600;">'
            f'{desc}</span></div>'
            + (f'<div style="font-size:11px;color:#6b7280;margin-top:2px;'
               f'padding-left:24px;">{prompt_preview}…</div>' if prompt_preview else '')
            + details
            + '</li>'
        )

    n_total = len(calls)
    header = (
        f'<div style="display:flex;align-items:center;justify-content:space-between;'
        f'gap:8px;margin-bottom:6px;">'
        f'<div style="font-size:11px;font-weight:700;color:#374151;'
        f'text-transform:uppercase;letter-spacing:0.4px;">Sub-Agents</div>'
        f'<div style="font-size:11px;color:#6b7280;">'
        f'{n_total} call(s)'
        + (f' · {n_running} running' if n_running else '')
        + '</div></div>'
    )
    return (
        '<div style="border:1px solid #e5e7eb;border-radius:8px;'
        'background:#fff;padding:10px 12px;'
        'font-family:-apple-system,BlinkMacSystemFont,sans-serif;">'
        + header
        + '<ul style="list-style:none;margin:0;padding:0;">'
        + "".join(rows) + '</ul></div>'
    )


def _render_todo_pane_html(todos: list[dict]) -> str:
    """Render the TodoWrite plan as a persistent sidebar HTML block.

    Returns an empty string if there are no todos — the caller should
    hide the widget.  Each todo carries ``content``, ``status``, and
    ``activeForm`` (shown while ``status == "in_progress"``).
    """
    if not todos:
        return ""
    n_done = sum(1 for t in todos if t.get("status") == "completed")
    n_total = len(todos)
    n_in_progress = sum(1 for t in todos if t.get("status") == "in_progress")

    rows: list[str] = []
    for t in todos:
        status = t.get("status", "pending")
        if status == "completed":
            glyph, color, opacity = "✓", "#10b981", "0.55"
            text = t.get("content", "")
        elif status == "in_progress":
            glyph, color, opacity = "→", "#3b82f6", "1.0"
            text = t.get("activeForm") or t.get("content", "")
        else:
            glyph, color, opacity = "○", "#9ca3af", "0.85"
            text = t.get("content", "")
        text = _html.escape(str(text)[:160])
        rows.append(
            f'<li style="display:flex;gap:8px;padding:4px 0;opacity:{opacity};">'
            f'<span style="color:{color};font-weight:700;flex:0 0 16px;">{glyph}</span>'
            f'<span style="color:#1f2937;font-size:12px;line-height:1.4;">{text}</span>'
            '</li>'
        )

    progress_pct = int((n_done / n_total) * 100) if n_total else 0
    header = (
        f'<div style="display:flex;align-items:center;justify-content:space-between;'
        f'gap:8px;margin-bottom:6px;">'
        f'<div style="font-size:11px;font-weight:700;color:#374151;'
        f'text-transform:uppercase;letter-spacing:0.4px;">Plan</div>'
        f'<div style="font-size:11px;color:#6b7280;">'
        f'{n_done}/{n_total} done'
        + (f' · {n_in_progress} active' if n_in_progress else '')
        + '</div></div>'
        f'<div style="height:3px;background:#e5e7eb;border-radius:2px;'
        f'overflow:hidden;margin-bottom:8px;">'
        f'<div style="height:100%;width:{progress_pct}%;'
        f'background:linear-gradient(90deg,#3b82f6,#10b981);"></div></div>'
    )
    return (
        '<div style="border:1px solid #e5e7eb;border-radius:8px;'
        'background:#f9fafb;padding:10px 12px;'
        'font-family:-apple-system,BlinkMacSystemFont,sans-serif;">'
        + header
        + '<ul style="list-style:none;margin:0;padding:0;">'
        + "".join(rows) + '</ul></div>'
    )


_FAIL_PHRASES = (
    "i'm unable", "i am unable", "i cannot", "can't complete",
    "ich kann nicht", "ich konnte nicht", "konnte das nicht",
    "abgebrochen", "aborting", "i give up", "ich gebe auf",
    "blocked by", "blockiert durch", "permission denied",
)
_PARTIAL_PHRASES = (
    "couldn't fully", "konnte nicht vollständig", "teilweise erfolgreich",
    "partially", "incomplete", "unfinished", "not all of", "nur einen teil",
)


def _compute_turn_verdict(
    *,
    engine,
    response_text: str,
    denied_commands: list[str],
    state: dict,
    error_type: str | None = None,
) -> str:
    """Decide PASS / FAIL / PARTIAL for a finished conversational turn.

    Heuristics are deliberately conservative — a wrong FAIL erodes trust
    in the success-rate metric faster than missing the occasional real
    failure does. The default is PASS.

    FAIL triggers (any):
    - explicit error_type (exception bubbled up from the stream loop)
    - stop_requested AND we have very little assistant output
    - denied commands present AND response contains a give-up phrase

    PARTIAL triggers (when no FAIL):
    - response ends with the QUESTION: tag (agent paused for the user)
    - response contains an "incomplete" phrase
    - denied commands present (work continued, but with restrictions)

    Otherwise: PASS.
    """
    text = (response_text or "").strip()
    text_low = text.lower()

    # ---- FAIL conditions -----------------------------------------------
    if error_type:
        return "FAIL"

    if getattr(engine, "_stop_requested", False) and len(text) < 80:
        # User hit stop and the agent had barely produced anything.
        return "FAIL"

    if denied_commands and any(p in text_low for p in _FAIL_PHRASES):
        # Tool calls were blocked AND the agent surrendered.
        return "FAIL"

    # ---- PARTIAL conditions --------------------------------------------
    # Look at the tail (last 800 chars) — most "I need more info" markers
    # show up at the end of the assistant reply.
    tail = text[-800:].lower()
    if "question:" in tail:
        return "PARTIAL"

    if any(p in text_low for p in _PARTIAL_PHRASES):
        return "PARTIAL"

    if denied_commands:
        # Tool calls were blocked but the agent didn't surrender — counts
        # as PARTIAL (the user got something, but not everything).
        return "PARTIAL"

    return "PASS"


def _build_inline_diff(old: str, new: str, _e, soft_cap: int = 60) -> str:
    """Render a +/- diff for a write/edit operation, always visible in chat.

    The user has explicitly asked to see every code change in the
    dashboard rather than have it tucked behind a `<details>` summary.
    This helper renders all `-` lines first then all `+` lines (mirroring
    `Edit`/`Write` semantics where the substring is replaced rather than
    line-diffed). For long blocks the middle portion is folded into a
    `<details>` so the head + tail remain visible without overwhelming
    the chat — head + tail is what tells you "did the right thing happen
    at the start" and "did the right thing land at the end".
    """
    old_lines = old.split("\n") if old else []
    new_lines = new.split("\n") if new else []
    diff_lines: list[str] = []
    for ln in old_lines:
        diff_lines.append(f'<span class="tool-diff-old">- {_e(ln)}</span>')
    for ln in new_lines:
        diff_lines.append(f'<span class="tool-diff-new">+ {_e(ln)}</span>')

    if not diff_lines:
        return ""

    style = (
        'margin:4px 0; padding:4px 6px; font-size:11px; line-height:1.45; '
        'background:rgba(127,127,127,0.06); border-left:2px solid #94a3b8;'
    )
    if len(diff_lines) <= soft_cap:
        return f'<pre style="{style}">{chr(10).join(diff_lines)}</pre>'

    half = max(soft_cap // 2, 8)
    head = chr(10).join(diff_lines[:half])
    middle = chr(10).join(diff_lines[half:-half])
    tail = chr(10).join(diff_lines[-half:])
    omitted = len(diff_lines) - 2 * half
    return (
        f'<pre style="{style} border-bottom:0;">{head}</pre>'
        f'<details style="margin:0 0 0 0;"><summary style="font-size:10px; '
        f'color:#9ca3af; padding:2px 6px; cursor:pointer;">'
        f'… {omitted} more lines</summary>'
        f'<pre style="{style} margin:0;">{middle}</pre></details>'
        f'<pre style="{style} border-top:0;">{tail}</pre>'
    )


def _record_turn_outcome(
    engine,
    user_task: str,
    response_text: str,
    state: dict,
    start_time: float | None,
    *,
    error_type: str | None = None,
) -> dict[str, str]:
    """Persist one completed turn into outcome history/profile.

    Tracks every conversational mode (solo, dashboard, quick, …) so the
    Activity tab can show real PASS/FAIL/PARTIAL across mode usage.
    Pipeline-internal cycles (session_manager etc.) still go through
    engine.record_cycle_outcome from the cycle gate, not this helper.

    B3: when ``state['_retry_pending_retries']`` > 0 (set by
    ``_retry_last``), the result of this turn is treated as a retry of
    the LAST outcome — we bump that outcome's retries counter and
    update its verdict in place instead of writing a fresh row.
    """
    if not (user_task or "").strip() or not (response_text or "").strip():
        return {}

    denied_commands: list[str] = []
    if isinstance(state, dict):
        denied_commands = list(state.get("_denied_commands", []))

    verdict = _compute_turn_verdict(
        engine=engine,
        response_text=response_text,
        denied_commands=denied_commands,
        state=state,
        error_type=error_type,
    )

    pending_retries = 0
    if isinstance(state, dict):
        pending_retries = int(state.get("_retry_pending_retries", 0) or 0)

    if pending_retries > 0:
        # Update last outcome's retries + verdict in place. Always
        # consume the flag so we don't double-bump on the next normal
        # turn even if the update fails.
        try:
            from delfin.agent.outcome_tracker import update_last_outcome
            ok = update_last_outcome(
                retries=pending_retries,
                verdict=verdict,
                denied_commands=denied_commands,
                error_type=error_type,
            )
        except Exception:
            ok = False
        if isinstance(state, dict):
            state["_retry_pending_retries"] = 0
        if ok:
            return {"retried": str(pending_retries), "verdict": verdict}
        # Fall through to a normal append if the bump failed (file
        # missing, permissions, etc.) — better to record something than
        # nothing.

    try:
        return engine.record_cycle_outcome(
            verdict=verdict,
            user_task=user_task[:200],
            denied_commands=denied_commands,
            error_type=error_type,
            start_time=start_time,
        )
    except Exception:
        return {}


# Backwards-compat alias — old callers might still import this name.
_record_solo_turn_outcome = _record_turn_outcome


# ---------------------------------------------------------------------------
# Tab creation
# ---------------------------------------------------------------------------


def create_tab(ctx):
    """Create the DELFIN Agent tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- check if pyyaml is available (needed for prompt loading) ----------
    try:
        import yaml  # noqa: F401
        _yaml_ok = True
    except ImportError:
        _yaml_ok = False

    _cli_available = bool(shutil.which("claude"))
    _codex_cli_available = bool(shutil.which("codex"))

    # Detect which providers are actually usable
    _available_providers: list[tuple[str, str]] = []
    if _cli_available or os.environ.get("ANTHROPIC_API_KEY", ""):
        _available_providers.append(("Claude", "claude"))
    if _codex_cli_available or os.environ.get("OPENAI_API_KEY", ""):
        _available_providers.append(("OpenAI", "openai"))
    if os.environ.get("KIT_TOOLBOX_API_KEY", ""):
        _available_providers.append(("KIT Toolbox", "kit"))
    # Detect a local OpenAI-compatible server (Ollama / vLLM / LM Studio /
    # llama.cpp). We probe the host's /api/tags (Ollama-native) and fall
    # back to /v1/models (everyone else). A 200 means at least one model
    # is loaded — that's enough to surface the provider in the dropdown.
    _ollama_host = (
        os.environ.get("OLLAMA_HOST")
        or os.environ.get("OLLAMA_BASE_URL")
        or "http://localhost:11434"
    )
    try:
        import urllib.request as _urlreq
        _probe = _urlreq.Request(
            _ollama_host.rstrip("/") + "/api/tags",
            headers={"Accept": "application/json"},
        )
        with _urlreq.urlopen(_probe, timeout=0.6) as _r:
            if _r.status == 200:
                _available_providers.append(("Ollama (local)", "ollama"))
    except Exception:
        # Try the OpenAI-compatible /v1/models endpoint as a fallback
        # so non-Ollama local servers (vLLM, LM Studio) also appear.
        try:
            import urllib.request as _urlreq
            _probe = _urlreq.Request(
                _ollama_host.rstrip("/") + "/v1/models",
                headers={"Accept": "application/json"},
            )
            with _urlreq.urlopen(_probe, timeout=0.6) as _r:
                if _r.status == 200:
                    _available_providers.append(("Ollama (local)", "ollama"))
        except Exception:
            pass
    # Fallback: always show at least Claude (will error with helpful message)
    if not _available_providers:
        _available_providers.append(("Claude", "claude"))

    # -- state -------------------------------------------------------------
    state = {
        "engine": None,
        "chat_messages": [],
        "streaming": False,
        "_generation_id": 0,      # monotonic counter — prevents stale worker cleanup
        "active_session_id": "",  # currently loaded CLI session ID
        "recent_edits": [],       # list of {"file": path, "tool": name} for undo
        "current_todos": [],      # latest TodoWrite payload (list of {content, status, activeForm})
        "subagent_calls": [],     # active/completed Agent tool calls for the subagent panel
        "_ws_known_files": set(),  # snapshot of agent_workspace at turn start (D4)
        "_job_states_seen": {},   # last-seen {job_id: status} for proactive notifications (D5)
        "_job_event_thread": None,  # background watcher (D5)
        "_job_event_thread_stop": False,
        "_seen_tab_suggestions": set(),  # tab names that already had a suggestion (D2)
        "_pending_action_text": "",      # raw agent text awaiting Approve/Deny (D1)
        "_pending_action_commands": [],  # parsed commands awaiting Approve/Deny (D1)
        "message_queue": [],      # queued messages sent while agent is busy
        "session_start_time": None,  # monotonic time of first message
        "_agent_calc_path": "",       # relative path within calc_dir for browsing
        "_pending_dashboard_action": None,  # {action_id, description, callback}
        "_perm_profile": "ask_all",  # permission profile: plan/ask_all/repo_free/all_free
        "_kit_confirm_broker": None,  # KitConfirmBroker (lazy-built when KIT is active)
        "_active_gate": None,        # {type, role, title, detail}
        "_last_stream_activity": 0.0,  # monotonic ts of last stream event (token/thinking/tool)
        "_stale_threshold_s": 600.0,   # >threshold without activity → spinner mode='stale'
        "_stale_timer": None,           # active threading.Timer for stale watch
        "_stale_kill_timer": None,      # cooperative-kill timer (dashboard only)
        "_stale_seen": False,           # already showed stale state this stream
        "_cycle_history": [],        # recent gate / handoff / retry events
        "_inspector_detail_key": "", # selected cycle inspector detail entry
        "_mode_manual_override": False,
        "_mode_change_internal": False,
        "_state_lock": __import__("threading").Lock(),  # protects streaming/gen_id/deny_count
    }

    # -- widgets -----------------------------------------------------------
    css_widget = widgets.HTML(value=_AGENT_CSS)

    # Mode selector
    _MODE_DESCRIPTIONS = {
        "dashboard": "Cheapest mode (Haiku) — operate the dashboard via slash commands: "
                     "set CONTROL keys, configure ORCA Builder, browse & analyze calculations, "
                     "trigger smart recalc, manage jobs. No code changes, session-only.",
        "solo": "Single agent, direct conversation — ask questions about the codebase, "
                "make small edits, explore code, debug issues. No review pipeline, "
                "fastest for simple tasks.",
        "plan": "Read-only research first — the agent explores the codebase, drafts a "
                "step-by-step plan in markdown, and waits for your approval via "
                "ExitPlanMode before any file edits or bash run.",
        "research": "Research agent — literature search, DFT benchmarks, best practices, "
                    "state-of-the-art methods. Web search enabled, read-only, no code changes.",
        "quick": "Session Manager → Builder → Test — lightweight pipeline for bugfixes, "
                 "docs updates, and isolated module changes. SM plans, Builder implements, "
                 "Test Agent verifies with pytest.",
        "reviewed": "Session Manager → Critic → Builder → Reviewer → Test — adds architectural "
                    "review before and code review after implementation. Use for risky refactors, "
                    "API changes, config semantics, or cross-module work.",
        "tdd": "Session Manager → Test → Builder → Reviewer → Test — test-driven development: "
               "Test Agent writes failing tests first, Builder implements until they pass, "
               "Reviewer checks quality, Test Agent verifies. Best for well-defined features.",
        "cluster": "Session Manager → Runtime → Critic → Builder → Reviewer → Test — "
                   "includes HPC/SLURM runtime specialist. Use for submission scripts, "
                   "job scheduling, scratch handling, error recovery, local vs. cluster differences.",
        "full": "Chief → Session Manager → Runtime → Critic → Builder → Reviewer → Test — "
                "maximum oversight with strategic lead (Chief). For releases, milestones, "
                "broad architectural changes. Most expensive, highest confidence.",
    }
    mode_dropdown = widgets.Dropdown(
        # Only solo + dashboard + plan are production-ready; the pipeline
        # modes (research/quick/reviewed/tdd/cluster/full) are gated until
        # they graduate from experimental status. Plan mode is a
        # read-only-first variant that produces a markdown plan and waits
        # for explicit approval before any edits land.
        options=["dashboard", "solo", "plan"],
        value="dashboard",
        description="Mode:",
        layout=widgets.Layout(width="200px"),
        style={"description_width": "45px"},
    )
    mode_desc_html = widgets.HTML(
        value=(
            f'<span style="color:#888; font-size:12px; margin-left:4px;">'
            f'{_MODE_DESCRIPTIONS.get("dashboard", "")}</span>'
        ),
        layout=widgets.Layout(width="auto"),
    )

    # Provider selector — models are fetched dynamically from the API
    # Fallback lists used only when the API is unreachable.
    _PROVIDER_MODELS_FALLBACK = {
        "claude": [("Opus", "opus"), ("Sonnet", "sonnet"), ("Haiku", "haiku")],
        "openai": [
            ("GPT-5.4", "gpt-5.4"),
            ("GPT-5.4-mini", "gpt-5.4-mini"),
            ("GPT-5.3-codex", "gpt-5.3-codex"),
            ("GPT-5.2-codex", "gpt-5.2-codex"),
            ("GPT-5.2", "gpt-5.2"),
            ("GPT-5.1-codex-max", "gpt-5.1-codex-max"),
            ("GPT-5.1-codex-mini", "gpt-5.1-codex-mini"),
            ("GPT-4.1", "gpt-4.1"),
            ("GPT-4.1-mini", "gpt-4.1-mini"),
            ("o4-mini", "o4-mini"),
            ("o3", "o3"),
        ],
        # Mirrors the live KIT-Toolbox /models response as of 2026-05.
        # Kept generous so a failed/empty fetch still gives users every
        # selectable model — the live fetch overwrites this when it
        # succeeds, so any newer additions appear automatically.
        "kit": [
            ("Azure GPT-5.1", "azure.gpt-5.1"),
            ("Azure GPT-5.4", "azure.gpt-5.4"),
            ("Azure GPT-5", "azure.gpt-5"),
            ("Azure GPT-5-mini", "azure.gpt-5-mini"),
            ("Azure GPT-5-nano", "azure.gpt-5-nano"),
            ("Azure o3", "azure.o3"),
            ("Azure o4-mini", "azure.o4-mini"),
            ("Azure GPT-4.1", "azure.gpt-4.1"),
            ("Azure GPT-4.1-mini", "azure.gpt-4.1-mini"),
            ("Azure GPT-4.1-nano", "azure.gpt-4.1-nano"),
            ("KIT gpt-oss 120B", "kit.gpt-oss-120b"),
            ("KIT gemma4 31B-it", "kit.gemma4-31b-it"),
            ("KIT minimax m2.7 229B", "kit.minimax-m2.7-229b"),
            ("KIT mistral-small-4 119B-a8b", "kit.mistral-small-4-119b-a8b"),
            ("KIT qwen3.5 397B-A17b", "kit.qwen3.5-397b-A17b"),
        ],
        # Common tool-calling-capable Ollama / vLLM / LM Studio models
        # as of 2026-05. The live /api/tags fetch overwrites this so
        # the user only sees models they actually have pulled.
        "ollama": [
            ("Qwen 3 — Coder 32B", "qwen3-coder:32b"),
            ("Qwen 2.5 — Coder 32B", "qwen2.5-coder:32b"),
            ("Qwen 2.5 — Coder 14B", "qwen2.5-coder:14b"),
            ("Qwen 2.5 — Coder 7B", "qwen2.5-coder:7b"),
            ("Llama 3.3 70B", "llama3.3:70b"),
            ("Llama 3.2 — 8B", "llama3.2:8b"),
            ("DeepSeek Coder V2 16B", "deepseek-coder-v2:16b"),
            ("DeepSeek Coder V2 — Lite", "deepseek-coder-v2:lite"),
            ("Mistral Nemo 12B", "mistral-nemo:12b"),
            ("Mistral Small 3 24B", "mistral-small:24b"),
        ],
    }
    _PROVIDER_DEFAULTS = {"claude": "sonnet", "openai": "gpt-5.4",
                          "kit": "azure.gpt-5.1",
                          "ollama": "qwen3-coder:32b"}
    _PROVIDER_CHEAP = {"claude": "haiku", "openai": "gpt-5.4-mini",
                       "kit": "azure.gpt-5-nano",
                       "ollama": "qwen2.5-coder:7b"}

    # Skip patterns: models that should not appear in the dropdown.
    # Embedding models can't generate chat completions — they'd error
    # out on the first send, so we hide them from the chat-model picker.
    _KIT_SKIP = {"standard-external", "standard-extern", "standard-local"}
    _MODEL_ID_HIDE_SUBSTRINGS = ("embedding", "embed-")

    def _fetch_models(provider):
        """Fetch model list from API. Returns [(label, id), ...] or None."""
        import json as _json
        try:
            if provider == "kit":
                key = os.environ.get("KIT_TOOLBOX_API_KEY", "")
                if not key:
                    return None
                url = "https://ki-toolbox.scc.kit.edu/api/v1/models"
                import urllib.request
                req = urllib.request.Request(
                    url, headers={"Authorization": f"Bearer {key}"}
                )
                with urllib.request.urlopen(req, timeout=8) as resp:
                    data = _json.loads(resp.read())
            elif provider == "openai":
                key = os.environ.get("OPENAI_API_KEY", "")
                if not key:
                    return None
                url = "https://api.openai.com/v1/models"
                import urllib.request
                req = urllib.request.Request(
                    url, headers={"Authorization": f"Bearer {key}"}
                )
                with urllib.request.urlopen(req, timeout=8) as resp:
                    data = _json.loads(resp.read())
            elif provider == "ollama":
                # Ollama serves a non-OpenAI ``/api/tags`` payload that
                # is friendlier (size, modified_at). Fall back to the
                # OpenAI-compatible ``/v1/models`` for vLLM / LM Studio.
                host = (
                    os.environ.get("OLLAMA_HOST")
                    or os.environ.get("OLLAMA_BASE_URL")
                    or "http://localhost:11434"
                ).rstrip("/")
                import urllib.request
                try:
                    with urllib.request.urlopen(
                        host + "/api/tags", timeout=4,
                    ) as resp:
                        raw = _json.loads(resp.read())
                    models = raw.get("models", [])
                    result = []
                    for m in models:
                        mid = m.get("name") or m.get("model", "")
                        if not mid:
                            continue
                        size_gb = (m.get("size") or 0) / 1_000_000_000
                        label = f"{mid}" + (
                            f"  ({size_gb:.1f} GB)" if size_gb > 0 else ""
                        )
                        result.append((label, mid))
                    if result:
                        result.sort(key=lambda x: x[1])
                        return result
                except Exception:
                    pass
                # OpenAI-compatible fallback (vLLM / LM Studio / llama.cpp)
                try:
                    with urllib.request.urlopen(
                        host + "/v1/models", timeout=4,
                    ) as resp:
                        raw = _json.loads(resp.read())
                    models = raw.get("data", [])
                    result = [(m.get("id", ""), m.get("id", ""))
                              for m in models if m.get("id")]
                    return result or None
                except Exception:
                    return None
            else:
                return None

            models = data.get("data", [])
            result = []
            for m in models:
                mid = m.get("id", "")
                if not mid or mid in _KIT_SKIP:
                    continue
                # Hide embedding / re-ranker models — they can't generate
                # chat completions and would error out at send time.
                if any(sub in mid.lower()
                       for sub in _MODEL_ID_HIDE_SUBSTRINGS):
                    continue
                label = mid.replace("azure.", "Azure ").replace("kit.", "KIT ")
                result.append((label, mid))
            # Sort: azure/cloud first, then local, alphabetically within groups
            result.sort(key=lambda x: (
                0 if x[1].startswith("azure.") else 1 if x[1].startswith("kit.") else 2,
                x[1],
            ))
            return result if result else None
        except Exception:
            return None

    # Use _PROVIDER_MODELS_FALLBACK as the initial value; dynamically updated later
    _PROVIDER_MODELS = dict(_PROVIDER_MODELS_FALLBACK)

    provider_dropdown = widgets.Dropdown(
        options=_available_providers,
        value=_available_providers[0][1],
        description="Provider:",
        layout=widgets.Layout(width="160px"),
        style={"description_width": "55px"},
    )

    # Model selector (options depend on selected provider)
    _init_provider = _available_providers[0][1]
    # Try fetching live models at startup
    _init_fetched = _fetch_models(_init_provider)
    if _init_fetched:
        _PROVIDER_MODELS[_init_provider] = _init_fetched
    # Stash init-fetch result so we can emit a status banner once
    # _append_system_message is wired up further down the build.
    state["_models_init_status"] = {
        "provider": _init_provider,
        "live": bool(_init_fetched),
        "count": (
            len(_init_fetched) if _init_fetched
            else len(_PROVIDER_MODELS[_init_provider])
        ),
        "key_present": (
            bool(os.environ.get("KIT_TOOLBOX_API_KEY", ""))
            if _init_provider == "kit"
            else bool(os.environ.get("OPENAI_API_KEY", ""))
            if _init_provider == "openai"
            else True
        ),
    }
    _init_models = _PROVIDER_MODELS[_init_provider]
    _init_default = _PROVIDER_DEFAULTS.get(_init_provider, _init_models[0][1])
    _init_valid = {v for _, v in _init_models}
    model_dropdown = widgets.Dropdown(
        options=_init_models,
        value=_init_default if _init_default in _init_valid else _init_models[0][1],
        description="Model:",
        layout=widgets.Layout(width="170px"),
        style={"description_width": "45px"},
    )
    model_dropdown.add_class("delfin-agent-model-dropdown")
    # Hidden refresh-trigger button — clicked programmatically from JS
    # when the user opens the model dropdown. Keeps the live fetch
    # off the visible UI while still giving users a fresh list every
    # time they want to pick a model.
    model_refresh_btn = widgets.Button(
        description="↻",
        tooltip="hidden auto-refresh trigger",
        layout=widgets.Layout(width="0px", height="0px",
                              display="none"),
    )
    model_refresh_btn.add_class("delfin-agent-model-refresh")

    def _on_refresh_models(_btn):
        provider = provider_dropdown.value
        fetched = _fetch_models(provider)
        if fetched:
            _PROVIDER_MODELS[provider] = fetched
            model_dropdown.options = fetched
            valid = {v for _, v in fetched}
            if model_dropdown.value not in valid:
                model_dropdown.value = fetched[0][1]
            model_refresh_btn.tooltip = (
                f"Refresh model list. Current: LIVE "
                f"({len(fetched)} models)"
            )
            _append_system_message(
                f"↻ {len(fetched)} models loaded live from {provider}."
            )
        else:
            hint = ""
            if provider == "kit" and not os.environ.get(
                "KIT_TOOLBOX_API_KEY", ""
            ):
                hint = (
                    " → KIT_TOOLBOX_API_KEY is not set in the dashboard "
                    "process. Set it before startup or via "
                    "`os.environ['KIT_TOOLBOX_API_KEY']='…'` in a "
                    "notebook cell and press ↻ again."
                )
            elif provider == "openai" and not os.environ.get(
                "OPENAI_API_KEY", ""
            ):
                hint = " → OPENAI_API_KEY is not set."
            _append_system_message(
                f"↻ Live fetch failed, fallback list active.{hint}"
            )

    model_refresh_btn.on_click(_on_refresh_models)
    # Effort selector (only affects API backend; CLI manages thinking internally)
    effort_dropdown = widgets.Dropdown(
        options=[
            ("Low", "low"),
            ("Medium", "medium"),
            ("High", "high"),
            ("Extra High", "xhigh"),
        ],
        value="medium",
        description="Effort:",
        layout=widgets.Layout(width="155px"),
        style={"description_width": "42px"},
        tooltip="Thinking budget multiplier. xhigh = 3× base (capped at 200k tokens).",
    )

    # Unified permission profile selector
    # Maps to BOTH zone permissions (slash commands) and CLI permission_mode (tools)
    perm_dropdown = widgets.Dropdown(
        options=[
            ("Plan", "plan"),
            ("Ask All", "ask_all"),
            ("Repo Free", "repo_free"),
            ("All Free", "all_free"),
        ],
        value="ask_all",
        description="Perms:",
        layout=widgets.Layout(width="195px"),
        style={"description_width": "42px"},
        tooltip="Permission profile: Plan=read + navigate/UI, Default=ask all changes, "
                "Erlaubt=repo free/calc asks, Full=all free (archive always read-only)",
    )

    # Load saved preferences
    try:
        from delfin.user_settings import load_settings
        _saved = (load_settings().get("agent", {}) or {})
        # Restore provider first (updates model options)
        _saved_provider = _saved.get("provider", "")
        if _saved_provider in ("claude", "openai", "kit"):
            provider_dropdown.value = _saved_provider
            model_dropdown.options = _PROVIDER_MODELS[_saved_provider]
            model_dropdown.value = _PROVIDER_DEFAULTS[_saved_provider]
        _saved_model = _saved.get("model", "")
        _valid_models = {v for _, v in model_dropdown.options}
        if _saved_model in _valid_models:
            model_dropdown.value = _saved_model
        _saved_effort = _saved.get("effort", "")
        if _saved_effort in ("low", "medium", "high", "xhigh"):
            effort_dropdown.value = _saved_effort
        _saved_perm = _saved.get("permission_profile", _saved.get("permission_mode", ""))
        # Migrate old permission names to new profiles
        _perm_migration = {
            "default": "ask_all", "erlaubt": "repo_free", "full": "all_free",
            "acceptEdits": "repo_free", "auto": "all_free",
            "bypassPermissions": "all_free",
        }
        _saved_perm = _perm_migration.get(_saved_perm, _saved_perm)
        if _saved_perm in ("plan", "ask_all", "repo_free", "all_free"):
            perm_dropdown.value = _saved_perm
            state["_perm_profile"] = _saved_perm
    except Exception:
        pass

    # Control buttons
    new_cycle_btn = widgets.Button(
        description="New Session",
        button_style="warning",
        layout=widgets.Layout(width="110px"),
        tooltip="Save current session and start fresh",
    )
    stop_btn = widgets.Button(
        description="Stop",
        button_style="danger",
        layout=widgets.Layout(width="80px"),
    )
    stop_btn.disabled = True

    advance_btn = widgets.Button(
        description="Next Role \u25b8",
        button_style="info",
        layout=widgets.Layout(width="110px"),
        tooltip="Advance to the next agent role in the cycle",
    )
    advance_btn.disabled = True

    # Git commit (agent-assisted) and push
    commit_btn = widgets.Button(
        description="Commit",
        button_style="",
        layout=widgets.Layout(width="90px"),
        tooltip="Ask the agent to commit current changes with a descriptive message",
    )
    push_btn = widgets.Button(
        description="Git Push",
        button_style="",
        layout=widgets.Layout(width="100px"),
        tooltip="Push committed changes (requires confirmation)",
    )
    push_confirm_btn = widgets.Button(
        description="Confirm Push",
        button_style="danger",
        layout=widgets.Layout(width="120px", display="none"),
    )
    push_cancel_btn = widgets.Button(
        description="Cancel",
        button_style="",
        layout=widgets.Layout(width="80px", display="none"),
    )
    push_status_html = widgets.HTML(value="")

    undo_btn = widgets.Button(
        description="Undo Edit",
        button_style="warning",
        layout=widgets.Layout(width="100px"),
        tooltip="Revert the last file edit (git checkout)",
    )
    undo_btn.disabled = True

    # Export button
    export_btn = widgets.Button(
        description="Export",
        button_style="",
        layout=widgets.Layout(width="80px"),
        tooltip="Export chat as Markdown file",
    )

    controls_row = widgets.VBox([
        widgets.HBox(
            [mode_dropdown, provider_dropdown, model_dropdown,
             effort_dropdown, perm_dropdown,
             new_cycle_btn, advance_btn, stop_btn, undo_btn, export_btn,
             commit_btn, push_btn, push_confirm_btn, push_cancel_btn, push_status_html,
             model_refresh_btn],
            layout=widgets.Layout(flex_flow="row wrap"),
        ),
        mode_desc_html,
    ], layout=widgets.Layout(margin="0 0 6px 0"))

    # Session selector
    session_dropdown = widgets.Dropdown(
        options=[("+ New Session", "")],
        value="",
        description="Session:",
        layout=widgets.Layout(width="380px"),
        style={"description_width": "55px"},
    )
    load_session_btn = widgets.Button(
        description="Load",
        button_style="info",
        layout=widgets.Layout(width="70px"),
        tooltip="Load the selected session",
    )
    delete_session_btn = widgets.Button(
        description="Delete",
        button_style="danger",
        layout=widgets.Layout(width="70px"),
        tooltip="Delete the selected session",
    )
    fork_session_btn = widgets.Button(
        description="Fork",
        button_style="warning",
        layout=widgets.Layout(width="70px"),
        tooltip="Duplicate the selected session under a new ID",
    )
    session_row = widgets.HBox(
        [session_dropdown, load_session_btn, fork_session_btn, delete_session_btn],
        layout=widgets.Layout(margin="0 0 6px 0"),
    )

    # Search bar (toggle visibility with Ctrl+K or /search)
    search_input = widgets.Text(
        placeholder="Search in chat...",
        layout=widgets.Layout(width="300px", display="none"),
        continuous_update=True,
    )
    search_count_html = widgets.HTML(value="")
    search_close_btn = widgets.Button(
        description="X",
        layout=widgets.Layout(width="30px", display="none"),
    )
    search_row = widgets.HBox(
        [search_input, search_count_html, search_close_btn],
        layout=widgets.Layout(margin="0 0 4px 0"),
    )

    # Permission approval widgets (shown inline between chat and input)
    approve_btn = widgets.Button(
        description="Approve (Enter)",
        button_style="success",
        layout=widgets.Layout(width="130px", height="34px", display="none"),
    )
    deny_btn = widgets.Button(
        description="Deny (Esc)",
        button_style="danger",
        layout=widgets.Layout(width="110px", height="34px", display="none"),
    )
    approval_info_html = widgets.HTML(value="")
    approval_row = widgets.HBox(
        [approval_info_html, approve_btn, deny_btn],
        layout=widgets.Layout(
            margin="4px 0", padding="6px 10px",
            display="none",
        ),
    )
    approval_row.add_class("delfin-agent-approval-row")

    # D1: Action-confirmation row — shown when the agent proposes ACTION:
    # commands AND ends with a confirmation question.  Lets the user
    # Approve/Deny via clicks instead of typing "ja".
    action_confirm_html = widgets.HTML(value="")
    action_approve_btn = widgets.Button(
        description="Ausführen",
        button_style="success",
        layout=widgets.Layout(width="120px", height="34px"),
    )
    action_deny_btn = widgets.Button(
        description="Abbrechen",
        button_style="danger",
        layout=widgets.Layout(width="110px", height="34px"),
    )
    action_confirm_row = widgets.HBox(
        [action_confirm_html, action_approve_btn, action_deny_btn],
        layout=widgets.Layout(
            margin="4px 0", padding="6px 10px",
            display="none", flex_flow="row wrap", gap="6px",
        ),
    )
    action_confirm_row.add_class("delfin-agent-approval-row")

    # Interactive question widgets (option buttons shown when agent asks)
    question_hint_html = widgets.HTML(value="")
    question_buttons_box = widgets.HBox(
        [],
        layout=widgets.Layout(gap="4px", flex_wrap="wrap"),
    )
    question_row = widgets.HBox(
        [question_hint_html, question_buttons_box],
        layout=widgets.Layout(
            margin="4px 0", padding="6px 10px",
            display="none",
        ),
    )
    question_row.add_class("delfin-agent-question-row")

    # Status bar
    status_html = widgets.HTML(
        value=_render_status("quick", "cli", "session_manager", 0, 3, 0, 0, 0.0),
    )
    cycle_inspector_html = widgets.HTML(
        value=(
            '<div class="delfin-cycle-inspector">'
            '<div class="inspector-header">'
            '<span class="inspector-title">Cycle Inspector</span>'
            '<span class="inspector-meta">Noch kein aktiver Cycle</span>'
            '</div>'
            '<div class="inspector-grid">'
            '<div class="inspector-card"><span class="card-label">Locked Goal</span>'
            '<span class="card-value">Noch kein gelocktes Ziel.</span></div>'
            '<div class="inspector-card gate-card"><span class="card-label">Current Gate</span>'
            '<span class="card-value">Kein aktives Gate</span></div>'
            '<div class="inspector-card"><span class="card-label">Next Role</span>'
            '<span class="card-value">Session Manager</span></div>'
            '<div class="inspector-card risk-card"><span class="card-label">Open Risks</span>'
            '<span class="card-value">Keine offenen Risiken.</span></div>'
            '<div class="inspector-card retry-card"><span class="card-label">Retry Count</span>'
            '<span class="card-value">0</span></div>'
            '<div class="inspector-card history-card"><span class="card-label">Recent Activity</span>'
            '<span class="card-value">Noch keine Aktivität.</span></div>'
            '</div></div>'
        ),
    )
    inspector_primary_btn = widgets.Button(
        description="Continue",
        button_style="success",
        layout=widgets.Layout(width="120px"),
        tooltip="Approve gate or continue the cycle",
    )
    inspector_retry_btn = widgets.Button(
        description="Retry Role",
        button_style="warning",
        layout=widgets.Layout(width="110px"),
        tooltip="Retry the last role output",
    )
    inspector_stop_btn = widgets.Button(
        description="Stop Cycle",
        button_style="danger",
        layout=widgets.Layout(width="110px"),
        tooltip="Stop the running or paused cycle",
    )
    inspector_next_btn = widgets.Button(
        description="Next Role",
        button_style="info",
        layout=widgets.Layout(width="110px"),
        tooltip="Advance to the next role without extra chat input",
    )
    inspector_actions_info = widgets.HTML(
        value='<span class="delfin-cycle-actions-info">Direkte Cycle-Steuerung.</span>',
    )
    inspector_actions_row = widgets.HBox(
        [inspector_primary_btn, inspector_retry_btn, inspector_stop_btn, inspector_next_btn, inspector_actions_info],
        layout=widgets.Layout(margin="0 0 6px 0"),
    )
    inspector_actions_row.add_class("delfin-cycle-actions")
    inspector_detail_dropdown = widgets.Dropdown(
        options=[("Keine Details verfügbar", "")],
        value="",
        layout=widgets.Layout(width="430px"),
    )
    inspector_detail_html = widgets.HTML(
        value=(
            '<div class="delfin-cycle-detail">'
            '<div class="detail-title">Keine Details verfügbar</div>'
            '<div class="detail-meta">Cycle Inspector</div>'
            '<div class="detail-body">Sobald Gates, History oder Rollen-Outputs vorliegen, erscheinen sie hier.</div>'
            '</div>'
        ),
    )
    inspector_detail_box = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<span class="delfin-cycle-detail-label">Role Detail Drawer</span>'),
                    inspector_detail_dropdown,
                ],
                layout=widgets.Layout(margin="0 0 4px 0"),
            ),
            inspector_detail_html,
        ],
        layout=widgets.Layout(margin="0 0 8px 0"),
    )
    inspector_detail_box.add_class("delfin-cycle-detail-box")

    # TodoWrite pane — persistent plan visible while the agent works.
    # Hidden by default; revealed on the first TodoWrite tool call.
    todo_pane_html = widgets.HTML(
        value="",
        layout=widgets.Layout(display="none", margin="0 0 6px 0"),
    )

    # Sub-agent pane — surface ``Agent`` tool calls (general-purpose, Explore,
    # Plan, …) with their status and result.  Hidden until first call.
    subagent_pane_html = widgets.HTML(
        value="",
        layout=widgets.Layout(display="none", margin="0 0 6px 0"),
    )


    # KIT-Toolbox confirmation broker (lazy — only built when KIT becomes
    # the active provider). The container is always present; the broker
    # widget is only mounted on demand so non-KIT users see nothing.
    kit_confirm_container = widgets.VBox(
        [],
        layout=widgets.Layout(display="none", margin="0 0 8px 0"),
    )

    # KIT workspace-roots status: read-only line that shows where the agent
    # can write/edit/bash. Adding new roots happens through the chat — the
    # agent calls remember_permission(kind='extra_dir', value=...) and the
    # Self-Modification Guard panel takes the user's confirmation.
    kit_dirs_status = widgets.HTML(value="")

    def _refresh_kit_dirs_status():
        eng = state.get("engine")
        if eng is None or not hasattr(eng, "list_kit_workspace_dirs"):
            kit_dirs_status.value = ""
            return
        roots = eng.list_kit_workspace_dirs()
        if not roots:
            kit_dirs_status.value = ""
            return

        snapshot: dict = {}
        if hasattr(eng, "kit_settings_snapshot"):
            try:
                snapshot = eng.kit_settings_snapshot() or {}
            except Exception:
                snapshot = {}
        allow_count = len(snapshot.get("allow_patterns", []) or [])
        deny_count = len(snapshot.get("deny_patterns", []) or [])

        # One-line summary: roots inline, comma-separated, repo first.
        roots_html = " · ".join(
            f"<code>{r}</code>" + (" <small>(repo)</small>" if i == 0 else "")
            for i, r in enumerate(roots)
        )
        kit_dirs_status.value = (
            "<small><b>Write access:</b> " + roots_html
            + f" &middot; allow-patterns: {allow_count}"
            + f" &middot; deny-patterns: {deny_count}"
            + " &middot; <i>outside: read-only with confirm</i>"
            + "<br><i>Tip: say <code>'also work in /path'</code> in chat "
            + "&rarr; agent persists it after one confirm click.</i></small>"
        )

    # KIT permission-mode picker (chip + cycle button).
    # Top-of-chat chip = current mode + verbose label.
    # Quick button next to Send = compact cycle button.
    _KIT_MODE_ORDER = ("plan", "default", "acceptEdits", "bypassPermissions")
    _KIT_MODE_VISUAL = {
        "plan":              ("Plan", "#5f6368", "#f1f3f4",
                              "Read-only · Agent schlägt vor, führt nichts aus"),
        "default":           ("Default", "#1a73e8", "#e8f0fe",
                              "Schreiben/Bash bestätigt jeden Schritt"),
        "acceptEdits":       ("Accept Edits", "#e8710a", "#fef7e0",
                              "Schreiben/Edit auto · Bash bestätigt"),
        "bypassPermissions": ("Bypass", "#c5221f", "#fce8e6",
                              "Alles auto · Sandbox + Denylist gelten weiter"),
    }

    kit_mode_chip = widgets.Button(
        description="Plan",
        tooltip="Click cycelt: plan → default → acceptEdits → bypass",
        layout=widgets.Layout(width="auto", height="32px",
                              margin="0 6px 0 0", display="none"),
    )
    kit_mode_label = widgets.HTML(value="", layout=widgets.Layout(display="none"))
    kit_mode_row = widgets.HBox(
        [widgets.HTML("<small><b>KIT-Mode:</b></small>"),
         kit_mode_chip, kit_mode_label],
        layout=widgets.Layout(
            align_items="center",
            margin="2px 0 4px 0",
            display="none",
        ),
    )

    kit_mode_quick_btn = widgets.Button(
        description="Plan",
        tooltip="Cycle KIT mode (next: default)",
        layout=widgets.Layout(width="80px", height="80px",
                              margin="0 0 0 4px", display="none"),
    )

    def _refresh_kit_mode_chip():
        eng = state.get("engine")
        perms = getattr(eng, "kit_permissions", None) if eng is not None else None
        if perms is None:
            for w in (kit_mode_chip, kit_mode_label, kit_mode_quick_btn):
                w.layout.display = "none"
            kit_mode_row.layout.display = "none"
            _refresh_plan_accept_btn()
            return
        mode = perms.mode if perms.mode in _KIT_MODE_VISUAL else "default"
        label, fg, bg, tip = _KIT_MODE_VISUAL[mode]
        # Top chip
        kit_mode_chip.description = label
        kit_mode_chip.tooltip = tip
        kit_mode_chip.style.button_color = bg
        kit_mode_label.value = (
            f"<small style='color:{fg}; margin-left:6px'>{tip}</small>"
        )
        kit_mode_chip.layout.display = ""
        kit_mode_label.layout.display = ""
        kit_mode_row.layout.display = "flex"
        # Quick button next to Send
        idx = _KIT_MODE_ORDER.index(mode)
        next_mode = _KIT_MODE_ORDER[(idx + 1) % len(_KIT_MODE_ORDER)]
        kit_mode_quick_btn.description = label
        kit_mode_quick_btn.tooltip = (
            f"KIT-Mode: {label}\nClick cycelt → {next_mode}"
        )
        kit_mode_quick_btn.style.button_color = bg
        kit_mode_quick_btn.layout.display = ""
        _refresh_plan_accept_btn()

    # Bidirectional mapping between the KIT-Mode chip (4 modes) and the
    # header Perms dropdown (4 profiles). They mean the same thing for KIT;
    # syncing them avoids the "they disagreed and the engine took the
    # dropdown value on recreation" surprise.
    _CHIP_TO_PROFILE = {
        "plan":              "plan",
        "default":           "ask_all",
        "acceptEdits":       "repo_free",
        "bypassPermissions": "all_free",
    }
    _PROFILE_TO_CHIP = {v: k for k, v in _CHIP_TO_PROFILE.items()}

    def _cycle_kit_mode(_btn=None):
        eng = state.get("engine")
        perms = getattr(eng, "kit_permissions", None) if eng is not None else None
        if perms is None:
            return
        cur = perms.mode if perms.mode in _KIT_MODE_ORDER else "default"
        nxt = _KIT_MODE_ORDER[(_KIT_MODE_ORDER.index(cur) + 1) % len(_KIT_MODE_ORDER)]
        ok = bool(eng.set_kit_permission_mode(nxt))
        if ok:
            _append_system_message(f"KIT-Mode → {nxt}")
            # Sync the header Perms dropdown so a future engine-recreate
            # starts in the same mode. Suppress the dropdown's verbose
            # "takes effect on next message" notice since the chip already
            # took effect live.
            target_profile = _CHIP_TO_PROFILE.get(nxt)
            if target_profile and perm_dropdown.value != target_profile:
                state["_chip_syncing_perm"] = True
                try:
                    perm_dropdown.value = target_profile
                finally:
                    state["_chip_syncing_perm"] = False
        _refresh_kit_mode_chip()

    kit_mode_chip.on_click(_cycle_kit_mode)
    kit_mode_quick_btn.on_click(_cycle_kit_mode)

    # Plan-mode accept button: shown only after agent answered while in plan mode.
    plan_accept_btn = widgets.Button(
        description="Accept plan & execute",
        button_style="success",
        tooltip=("Switches to 'acceptEdits' and prompts the agent to "
                 "execute the proposed plan now."),
        layout=widgets.Layout(display="none", margin="4px 0"),
    )

    def _refresh_plan_accept_btn():
        eng = state.get("engine")
        perms = getattr(eng, "kit_permissions", None) if eng is not None else None
        # Show when we're in plan mode AND the agent has produced at least
        # one assistant message in this session.
        has_assistant_turn = bool(state.get("_kit_plan_has_response"))
        if perms is not None and perms.mode == "plan" and has_assistant_turn:
            plan_accept_btn.layout.display = ""
        else:
            plan_accept_btn.layout.display = "none"

    def _on_plan_accept(_btn):
        eng = state.get("engine")
        if eng is None or not hasattr(eng, "set_kit_permission_mode"):
            return
        if not eng.set_kit_permission_mode("acceptEdits"):
            _append_system_message("Mode switch to acceptEdits failed.")
            return
        _append_system_message("Mode → acceptEdits · sending plan-execute command …")
        state["_kit_plan_has_response"] = False
        _refresh_kit_mode_chip()
        # Inject a follow-up user message that triggers execution.
        try:
            input_textarea.value = (
                "Please execute the proposed plan now. "
                "Use the KIT-Toolbox tools (write_file / edit_file / bash) "
                "and briefly report after each step what you did."
            )
            _on_send(None)
        except Exception as exc:
            _append_system_message(f"Auto-send failed: {exc}")

    plan_accept_btn.on_click(_on_plan_accept)

    def _ensure_kit_broker():
        """Build/return the KitConfirmBroker bound to this tab."""
        if state.get("_kit_confirm_broker") is None:
            try:
                from delfin.agent.kit_confirm import KitConfirmBroker
                broker = KitConfirmBroker()
                # Wire the "Erlauben + Dauerhaft" button to the engine's
                # persist hooks. Bound lazily so it always picks up the
                # currently-active engine (provider may switch at runtime).
                # ``kind`` is 'allow' / 'deny' (bash-pattern persistence)
                # or 'extra_dir' (workspace-dir persistence — used when
                # the user permanently grants access to a directory the
                # agent tried to read).
                def _persist(kind: str, value: str) -> tuple[bool, str]:
                    eng = state.get("engine")
                    if eng is None:
                        return False, "KIT-Engine nicht aktiv"
                    if kind in ("allow", "deny"):
                        if not hasattr(eng, "persist_kit_pattern"):
                            return False, "persist_kit_pattern fehlt"
                        return eng.persist_kit_pattern(value, kind=kind)
                    if kind == "extra_dir":
                        if not hasattr(eng, "add_kit_workspace_dir"):
                            return False, "add_kit_workspace_dir fehlt"
                        return eng.add_kit_workspace_dir(value, persist=True)
                    return False, f"unbekannter persist-kind: {kind}"
                broker.set_persist_callback(_persist)
                panel = broker.build_widget()
                kit_confirm_container.children = (panel,)
                state["_kit_confirm_broker"] = broker
            except Exception as exc:
                _append_system_message(f"KIT confirm broker init failed: {exc}")
                state["_kit_confirm_broker"] = None
        _refresh_kit_dirs_status()
        return state.get("_kit_confirm_broker")

    def _show_kit_confirm_panel(visible: bool) -> None:
        kit_confirm_container.layout.display = "flex" if visible else "none"
        if visible:
            _refresh_kit_dirs_status()
        _refresh_kit_mode_chip()

    # Chat display
    chat_html = widgets.HTML(
        value='<div class="delfin-agent-chat"><i>Start a conversation...</i></div>',
        layout=widgets.Layout(min_height="200px"),
    )

    # ----- Slash-command palette (discoverable command browser) -----
    palette_toggle_btn = widgets.Button(
        description="/", icon="terminal",
        layout=widgets.Layout(width="40px", height="32px"),
        tooltip="Show / browse slash commands",
    )
    palette_search = widgets.Text(
        placeholder="Filter commands… (e.g. control, recalc, git)",
        layout=widgets.Layout(width="100%", height="32px", display="none"),
    )
    palette_select = widgets.Select(
        options=[],
        rows=10,
        layout=widgets.Layout(width="100%", display="none", margin="4px 0"),
    )
    palette_row = widgets.HBox(
        [palette_toggle_btn, palette_search],
        layout=widgets.Layout(gap="6px", align_items="center"),
    )

    def _palette_format_options(query: str = "") -> list[tuple[str, str]]:
        """Return Select options: ('/cmd  —  summary', insert-text).

        Built-in slash commands first, then any installed skills as
        ``/skill <name>`` entries (filtered by the same query).
        """
        items = _filter_slash_commands(_SLASH_COMMANDS, query)
        out: list[tuple[str, str]] = [
            (f"{cmd:<22}  {summary}",
             cmd + (" " if has_args else ""))
            for _category, cmd, summary, has_args in items
        ]
        # Append discovered skills, filtered by the same query.
        # The Skill dataclass exposes `.name`, `.description`, `.body`,
        # `.source` only — earlier revisions referenced `.title` and
        # `.slash_command` which never existed; that crashed the palette
        # the moment any skill was discovered.
        try:
            from delfin.agent.skills import discover_skills
            q = (query or "").strip().lower()
            for skill in discover_skills():
                first_line = (skill.description
                              or (skill.body or "").splitlines()[0]
                              if skill.body else "")
                hay = f"/skill {skill.name} {first_line}".lower()
                if not q or q in hay:
                    short = (first_line or "")[:60]
                    label = f"/skill {skill.name:<14}  {short}"
                    out.append((label, f"/skill {skill.name} "))
        except Exception:
            pass
        if not out:
            return [("(no matches)", "")]
        return out

    def _palette_set_visible(visible: bool) -> None:
        palette_search.layout.display = "block" if visible else "none"
        palette_select.layout.display = "block" if visible else "none"
        if visible:
            palette_search.value = ""
            palette_select.options = _palette_format_options("")
        else:
            palette_select.value = None

    def _on_palette_toggle(_btn) -> None:
        _palette_set_visible(palette_select.layout.display == "none")

    def _on_palette_search(change) -> None:
        palette_select.options = _palette_format_options(change.get("new", ""))

    def _on_palette_select(change) -> None:
        insert = change.get("new", "")
        if not insert:
            return
        # Insert at start of input (palette is for starting commands)
        input_textarea.value = insert
        _palette_set_visible(False)
        try:
            input_textarea.focus()
        except Exception:
            pass

    palette_toggle_btn.on_click(_on_palette_toggle)
    palette_search.observe(_on_palette_search, names="value")
    palette_select.observe(_on_palette_select, names="value")

    # Input area — textarea stretches to take all remaining width;
    # send + mode buttons keep fixed footprint at full input height
    # so the row looks like one coherent strip.
    input_textarea = widgets.Textarea(
        placeholder=(
            "Message the agent... (Enter = send, Shift+Enter = newline)\n"
            "KIT tip: say 'also work in /path' to grant write access — "
            "the agent persists it after one confirm click."
        ),
        layout=widgets.Layout(
            flex="1 1 auto", width="auto", height="80px",
        ),
    )
    input_textarea.add_class("delfin-agent-input")
    send_btn = widgets.Button(
        description="Send",
        button_style="primary",
        layout=widgets.Layout(width="80px", height="80px"),
    )
    # image_upload is created later (line ~2240) but inserted here so
    # the whole row [📎 | textarea | Send | quick-mode] is one strip.
    input_row = widgets.HBox(
        [input_textarea, send_btn, kit_mode_quick_btn],
        layout=widgets.Layout(
            margin="0 0 0 4px",
            flex="1 1 auto", width="auto",
            align_items="flex-start",
        ),
    )
    input_row.add_class("delfin-agent-send-row")

    # Working indicator (animated spinner)
    working_html = widgets.HTML(value="")

    # Queue indicator (shown when messages are queued during streaming)
    queue_html = widgets.HTML(value="")

    # Live context-window usage bar — refreshed by _refresh_context_bar()
    # on a periodic timer + after every send/stop/compact. Hidden until the
    # engine has at least one message so the bar doesn't flash at zero on
    # first paint.
    context_bar_html = widgets.HTML(value="")

    # Keyboard shortcuts: dedicated Output widget (never cleared by scroll)
    _enter_js_output = widgets.Output()
    with _enter_js_output:
        from IPython.display import display as _ipyd, Javascript as _JS
        _ipyd(_JS("""
(function() {
    if (window.__delfinAgentKeys) return;
    window.__delfinAgentKeys = true;
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Enter' && !e.shiftKey && !e.ctrlKey && !e.metaKey) {
            if (e.target && e.target.tagName === 'TEXTAREA') {
                var container = e.target.closest
                    ? e.target.closest('.delfin-agent-input') : null;
                if (container) {
                    e.preventDefault();
                    e.stopPropagation();
                    var sendBtn = document.querySelector('.delfin-agent-send-row button');
                    if (sendBtn) sendBtn.click();
                    return;
                }
            }
        }
        if (e.key === 'Escape') {
            var btns = document.querySelectorAll('button');
            for (var i = 0; i < btns.length; i++) {
                if (btns[i].textContent.trim() === 'Stop' && !btns[i].disabled) {
                    btns[i].click(); e.preventDefault(); return;
                }
            }
        }
        if ((e.ctrlKey || e.metaKey) && e.key === 'l') {
            var a = document.querySelector('.delfin-agent-chat');
            if (a) {
                e.preventDefault();
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var ns = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    ns.call(ta, '/clear');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sb = document.querySelector('.delfin-agent-send-row button');
                        if (sb) sb.click();
                    }, 50);
                }
            }
        }
        if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
            var a = document.querySelector('.delfin-agent-chat');
            if (a) {
                e.preventDefault();
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var ns = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    ns.call(ta, '/search');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sb = document.querySelector('.delfin-agent-send-row button');
                        if (sb) sb.click();
                    }, 50);
                }
            }
        }
        // Shift+Tab: cycle permission mode
        if (e.key === 'Tab' && e.shiftKey && !e.ctrlKey && !e.metaKey) {
            var a = document.querySelector('.delfin-agent-chat');
            if (a) {
                e.preventDefault();
                e.stopPropagation();
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var ns = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    ns.call(ta, '/perm-cycle');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sb = document.querySelector('.delfin-agent-send-row button');
                        if (sb) sb.click();
                    }, 50);
                }
            }
        }
    }, true);

    // --- Auto-scroll: poll every 150ms while agent is working ---
    (function() {
        if (window.__delfinChatScroll) return;
        window.__delfinChatScroll = true;
        setInterval(function() {
            // Only auto-scroll when the working indicator is visible
            var working = document.querySelector('.delfin-agent-working');
            if (!working) return;
            var chat = document.querySelector('.delfin-agent-chat');
            if (chat) chat.scrollTop = chat.scrollHeight;
        }, 150);
    })();
})();
"""))

    # -- UI widget registry (for /ui command) --------------------------------
    _ui_widgets: dict[str, widgets.Widget] = {
        "send-btn": send_btn,
        "input": input_textarea,
        "mode": mode_dropdown,
        "perm": perm_dropdown,
    }

    # Widgets that require user confirmation before the agent can click them.
    # The agent CAN click these, but only after the user confirms.
    _CONFIRM_WIDGETS = frozenset({
        "calc-delete-btn", "remote-delete-btn",
        "submit-btn",              # submits real job
        "orca-submit-btn",         # submits real ORCA job
        "calc-recalc-btn",         # triggers recalculation
        "calc-submit-recalc-btn",  # submits recalc
    })
    # Hard-blocked: agent can NEVER click these, even with confirmation
    _BLOCKED_WIDGETS = frozenset({
        "calc-delete-btn", "remote-delete-btn",
    })

    def _build_full_widget_registry() -> dict[str, widgets.Widget]:
        """Build registry including widgets from all dashboard tabs."""
        reg = dict(_ui_widgets)  # start with agent-tab widgets

        # Submit tab
        _submit_map = {
            "job-name": "job_name_widget",
            "control": "control_widget",
            "coords": "coords_widget",
            "submit-btn": "submit_button",
            "batch-smiles": "smiles_batch_widget",
            "time-limit": "job_type_widget",
            "custom-time": "custom_time_widget",
        }
        for alias, ref_key in _submit_map.items():
            w = ctx.submit_refs.get(ref_key)
            if w is not None:
                reg[alias] = w

        # ORCA Builder tab
        _orca_map = {
            "orca-method": "orca_method",
            "orca-basis": "orca_basis",
            "orca-job-type": "orca_job_type",
            "orca-charge": "orca_charge",
            "orca-mult": "orca_multiplicity",
            "orca-pal": "orca_pal",
            "orca-maxcore": "orca_maxcore",
            "orca-coords": "orca_coords",
            "orca-dispersion": "orca_dispersion",
            "orca-solvent": "orca_solvent",
            "orca-preview": "orca_preview",
            # Safe buttons (no real computation, no deletion)
            "orca-convert-btn": "orca_convert_smiles_btn",
            "orca-copy-coords-btn": "orca_copy_coords_btn",
            "orca-check-numbering-btn": "orca_check_numbering_btn",
            "orca-apply-numbering-btn": "orca_apply_numbering_btn",
            "orca-save-btn": "orca_save_btn",
            "orca-mol-prev-btn": "orca_mol_prev_btn",
            "orca-mol-next-btn": "orca_mol_next_btn",
            # Destructive buttons — registered but require user confirmation
            "orca-submit-btn": "orca_submit_btn",
        }
        for alias, ref_key in _orca_map.items():
            w = ctx.orca_builder_refs.get(ref_key)
            if w is not None:
                reg[alias] = w

        # Calculations Browser tab
        _calc_map = {
            # Navigation
            "calc-path": "calc_path_input",
            "calc-sort": "calc_sort_dropdown",
            "calc-filter": "calc_folder_search",
            "calc-search": "calc_search_input",
            # File operations
            "calc-new-folder-btn": "calc_new_folder_btn",
            "calc-new-folder-name": "calc_new_folder_input",
            "calc-rename-btn": "calc_rename_btn",
            "calc-rename-name": "calc_rename_input",
            "calc-duplicate-btn": "calc_duplicate_btn",
            "calc-copy-btn": "calc_copy_btn",
            "calc-copy-path-btn": "calc_copy_path_btn",
            # Transfer / move
            "calc-to-archive-btn": "calc_move_archive_btn",
            "calc-to-calc-btn": "calc_back_to_calculations_btn",
            "calc-ssh-btn": "calc_ssh_transfer_btn",
            # Visualization & report
            "calc-visualize": "calc_view_toggle",
            "calc-png-btn": "calc_view_png_btn",
            "calc-xyz-png-btn": "calc_xyz_png_btn",
            "calc-report-btn": "calc_report_btn",
            # Extract Table
            "calc-table-btn": "calc_table_btn",
            "calc-table-file": "calc_table_file_input",
            "calc-table-scope": "calc_table_scope_dd",
            "calc-table-recursive": "calc_table_recursive_cb",
            "calc-table-decimal": "calc_table_decimal_comma_btn",
            "calc-table-preset": "calc_table_preset_name",
            "calc-table-save-btn": "calc_table_preset_save_btn",
            "calc-table-add-col-btn": "calc_table_add_col_btn",
            "calc-table-run-btn": "calc_table_run_btn",
            "calc-table-csv-btn": "calc_table_csv_btn",
            "calc-table-output": "calc_table_output",
            # Recalc & editor
            "calc-recalc-btn": "calc_recalc_btn",
            "calc-recalc-time": "calc_recalc_time",
            "calc-submit-recalc-btn": "calc_submit_recalc_btn",
            "calc-editor": "calc_edit_area",
            # Options dropdown (Recalc / Smart Recalc / Override / etc.)
            "calc-options": "calc_options_dropdown",
            "calc-override": "calc_override_input",
            "calc-override-time": "calc_override_time",
            "calc-override-btn": "calc_override_btn",
            # Delete (blocked — agent cannot click)
            "calc-delete-btn": "calc_delete_btn",
        }
        for alias, ref_key in _calc_map.items():
            w = ctx.calc_browser_refs.get(ref_key)
            if w is not None:
                reg[alias] = w

        # Remote Archive tab
        _remote_map = {
            "remote-path": "path_input",
            "remote-filter": "filter_input",
            "remote-sort": "sort_dropdown",
            "remote-search": "search_input",
            # File operations
            "remote-new-folder-btn": "new_folder_btn",
            "remote-new-folder-name": "new_folder_input",
            "remote-rename-btn": "rename_btn",
            "remote-rename-name": "rename_input",
            "remote-duplicate-btn": "duplicate_btn",
            "remote-copy-btn": "copy_btn",
            "remote-copy-path-btn": "copy_path_btn",
            "remote-download-btn": "download_btn",
            # Transfer
            "remote-to-calc-btn": "transfer_back_btn",
            "remote-to-archive-btn": "transfer_to_archive_btn",
            "remote-transfers-btn": "transfer_jobs_btn",
            "remote-transfers-refresh-btn": "transfer_jobs_refresh_btn",
            # Visualization
            "remote-visualize": "view_toggle",
            "remote-png-btn": "viewer_png_btn",
            # Extract Table
            "remote-table-btn": "table_btn",
            "remote-table-file": "table_file_input",
            "remote-table-scope": "table_scope_dd",
            "remote-table-recursive": "table_recursive_cb",
            "remote-table-decimal": "table_decimal_comma_btn",
            "remote-table-preset": "table_preset_name",
            "remote-table-save-btn": "table_preset_save_btn",
            "remote-table-add-col-btn": "table_add_col_btn",
            "remote-table-run-btn": "table_run_btn",
            "remote-table-csv-btn": "table_csv_btn",
            "remote-table-output": "table_output_html",
            # Delete (blocked — agent cannot click)
            "remote-delete-btn": "delete_btn",
        }
        _ra_refs = getattr(ctx, 'remote_archive_refs', {})
        for alias, ref_key in _remote_map.items():
            w = _ra_refs.get(ref_key)
            if w is not None:
                reg[alias] = w

        return reg

    # -- Phase 5 UI hookups -----------------------------------------------
    # Live task ticker — reflects TaskStore state, refreshed on tool result.
    task_ticker_html = widgets.HTML(
        value="", layout=widgets.Layout(margin="2px 0 4px 0"),
    )

    def _refresh_task_ticker():
        try:
            from delfin.agent.task_ticker import render_html as _tt_render
            eng = state.get("engine")
            ws = None
            if eng is not None:
                kp = getattr(eng, "kit_permissions", None)
                if kp is not None:
                    ws = kp.workspace
            if ws is None:
                ws = ctx.repo_dir or Path.cwd()
            task_ticker_html.value = _tt_render(
                ws, session_id=str(state.get("active_session_id", "") or "")
            )
        except Exception:
            task_ticker_html.value = ""

    # Status line footer — token / mode / branch summary.
    status_line_html = widgets.HTML(
        value="", layout=widgets.Layout(margin="4px 0 0 0"),
    )

    def _refresh_status_line():
        try:
            from delfin.agent.status_line import (
                StatusContext as _SC, render_status_line as _r_status,
            )
            eng = state.get("engine")
            ws = None
            mode = "default"
            tokens = 0
            cost = 0.0
            model = ""
            if eng is not None:
                kp = getattr(eng, "kit_permissions", None)
                if kp is not None:
                    ws = kp.workspace
                    mode = kp.mode
                tokens = (eng.token_usage.get("input", 0)
                          + eng.token_usage.get("output", 0))
                cost = eng.cost_usd
                model = getattr(eng.client, "model", "")
            ctx_obj = _SC(
                workspace=ws, mode=mode, model=model,
                tokens=tokens, cost_usd=cost,
            )
            line = _r_status(ctx_obj)
            status_line_html.value = (
                f"<div style='font-family:monospace;font-size:11px;"
                f"color:#777;border-top:1px solid #333;padding:3px 6px;'>"
                f"{line}</div>" if line else ""
            )
        except Exception:
            status_line_html.value = ""

    # AskUserQuestion modal — shown only while a tool call is pending.
    ask_user_label = widgets.HTML(value="")
    ask_user_buttons = widgets.HBox([],
        layout=widgets.Layout(flex_flow="row wrap", gap="6px"))
    ask_user_box = widgets.VBox(
        [ask_user_label, ask_user_buttons],
        layout=widgets.Layout(
            display="none",
            border="2px solid #0a84ff",
            padding="8px",
            margin="6px 0",
            background_color="#0a84ff15",
        ),
    )
    state["_ask_user_event"] = None
    state["_ask_user_result"] = None

    def _show_ask_user_modal(args: dict) -> dict:
        """Build buttons for each option and block until the user clicks.

        When at least one option carries a ``preview`` field AND the
        question is single-select, the layout switches to a side-by-side
        comparison: option cards on the left, the currently-focused
        option's markdown preview on the right.
        """
        import threading as _th
        question = args.get("question", "")
        header = args.get("header", "")
        options = args.get("options", []) or []
        multi = bool(args.get("multiSelect", False))
        ev = _th.Event()
        result = {"answers": []}
        state["_ask_user_event"] = ev
        state["_ask_user_result"] = result
        from html import escape as _esc
        header_html = (f"<span style='background:#0a84ff;color:white;"
                       f"padding:2px 6px;border-radius:3px;font-size:11px;"
                       f"margin-right:6px;'>{_esc(header)}</span>"
                       if header else "")
        ask_user_label.value = (
            f"<div style='font-size:14px;'>{header_html}"
            f"<b>{_esc(question)}</b></div>"
        )
        chosen: set[str] = set()
        button_widgets: list = []

        has_previews = (
            not multi
            and any((opt.get("preview") or "").strip() for opt in options)
        )

        # Preview pane state: shared HTML widget; click on an option
        # updates it. The first option's preview is shown on open.
        preview_html = None
        if has_previews:
            preview_html = widgets.HTML(value="")

            def _render_preview(opt: dict) -> str:
                body = (opt.get("preview") or "").strip()
                label = opt.get("label") or ""
                title = (opt.get("description") or "").strip()
                title_block = (
                    f"<div style='color:#888;font-size:11px;margin-bottom:6px;'>"
                    f"{_esc(title)}</div>" if title else ""
                )
                return (
                    f"<div style='border:1px solid #444;border-radius:6px;"
                    f"padding:8px;background:#111;font-family:monospace;"
                    f"font-size:12px;color:#ddd;white-space:pre;"
                    f"overflow-x:auto;max-height:380px;'>"
                    f"<div style='color:#60a5fa;font-weight:600;"
                    f"margin-bottom:4px;'>Preview: {_esc(label)}</div>"
                    f"{title_block}{_esc(body)}</div>"
                )

            if options:
                preview_html.value = _render_preview(options[0])

        def _make_handler(label: str, btn, opt: dict):
            def _on_click(_b):
                if multi:
                    if label in chosen:
                        chosen.discard(label)
                        btn.button_style = ""
                    else:
                        chosen.add(label)
                        btn.button_style = "primary"
                elif has_previews:
                    # Two-stage: first click shows preview + arms confirm.
                    # Already-armed click commits.
                    if chosen == {label}:
                        result["answers"] = [label]
                        ev.set()
                    else:
                        chosen.clear()
                        chosen.add(label)
                        for w in button_widgets:
                            if getattr(w, "_aubu_label", "") == label:
                                w.button_style = "primary"
                                w.description = f"✓ {label[:58]}"
                            elif hasattr(w, "_aubu_label"):
                                w.button_style = ""
                                w.description = w._aubu_label[:60]
                        if preview_html is not None:
                            preview_html.value = _render_preview(opt)
                else:
                    chosen.clear()
                    chosen.add(label)
                    result["answers"] = [label]
                    ev.set()
            return _on_click

        for opt in options:
            label = opt.get("label", "")
            desc = opt.get("description", "")
            btn = widgets.Button(
                description=label[:60],
                tooltip=desc[:200],
                layout=widgets.Layout(
                    width=("100%" if has_previews else "auto"),
                ),
            )
            btn._aubu_label = label  # tag for re-styling on selection
            btn.on_click(_make_handler(label, btn, opt))
            button_widgets.append(btn)
        if multi:
            confirm_btn = widgets.Button(
                description="OK", button_style="success",
                layout=widgets.Layout(width="60px"),
            )
            def _on_confirm(_b):
                result["answers"] = list(chosen)
                ev.set()
            confirm_btn.on_click(_on_confirm)
            button_widgets.append(confirm_btn)

        if has_previews:
            # Vertical option list on the left, preview pane on the right.
            options_col = widgets.VBox(
                button_widgets,
                layout=widgets.Layout(
                    width="40%",
                    gap="4px",
                ),
            )
            hint = widgets.HTML(
                value=(
                    "<div style='color:#888;font-size:11px;margin-top:6px;'>"
                    "Click an option to preview it; click again to confirm."
                    "</div>"
                ),
            )
            preview_col = widgets.VBox(
                [preview_html, hint],
                layout=widgets.Layout(width="60%"),
            )
            ask_user_buttons.layout.flex_flow = "row"
            ask_user_buttons.layout.gap = "12px"
            ask_user_buttons.children = (options_col, preview_col)
        else:
            ask_user_buttons.layout.flex_flow = "row wrap"
            ask_user_buttons.layout.gap = "6px"
            ask_user_buttons.children = tuple(button_widgets)
        ask_user_box.layout.display = ""
        # Block up to 5 minutes
        ev.wait(timeout=300)
        ask_user_box.layout.display = "none"
        ask_user_label.value = ""
        ask_user_buttons.children = ()
        state["_ask_user_event"] = None
        state["_ask_user_result"] = None
        if not result["answers"]:
            return {"answers": [], "timed_out": True}
        return result

    # ExitPlanMode approval — uses existing plan_accept_btn but flips perms
    # via the structured callback so the tool result reflects the choice.
    state["_plan_approval_event"] = None
    state["_plan_approval_result"] = None

    def _show_plan_approval(plan: str) -> dict:
        """Render the plan as a system message, wait for user click."""
        import threading as _th
        ev = _th.Event()
        result: dict = {"approved": False, "new_mode": "default"}
        state["_plan_approval_event"] = ev
        state["_plan_approval_result"] = result
        # Stash the plan body so the accept-handler can persist it.
        state["_pending_plan_body"] = plan
        # Surface the plan so the user sees what they're approving.
        _append_system_message(
            "📋 **Plan zur Freigabe** (klicke 'Plan akzeptieren' "
            "oder wechsle Mode-Chip):\n\n" + plan
        )
        # Reuse existing plan_accept_btn machinery — the callback below
        # observes state changes from _on_plan_accept.
        state["_kit_plan_has_response"] = True
        _refresh_plan_accept_btn()
        ev.wait(timeout=600)   # 10 min
        return result

    # Hook into the existing _on_plan_accept: register a second on_click
    # handler that signals any pending exit_plan_mode tool call. The
    # original handler (registered earlier) still does its UI/mode work
    # — both fire on click.
    def _on_plan_accept_phase5(_btn):
        ev = state.get("_plan_approval_event")
        result = state.get("_plan_approval_result")
        if ev is not None and result is not None:
            result["approved"] = True
            result["new_mode"] = "acceptEdits"
            # Persist the approved plan so the user can re-open it later
            # without going through the tool round-trip.
            plan_body = state.pop("_pending_plan_body", "") or ""
            if plan_body.strip():
                try:
                    from delfin.agent.memory_store import save_plan
                    fpath = save_plan(plan_body, repo_root=ctx.repo_dir or ".")
                    short = str(fpath).replace(str(Path.home()), "~")
                    result["plan_path"] = str(fpath)
                    _append_system_message(f"📝 Plan saved → {short}")
                except Exception as exc:
                    _append_system_message(f"Plan save failed: {exc}")
            ev.set()

    plan_accept_btn.on_click(_on_plan_accept_phase5)

    # File upload widget alongside chat input. Matches the textarea
    # height (80px) so the input row looks like one coherent strip.
    # Accepts any file type; uploads land in <workspace>/.delfin/uploads/
    # and the agent picks them up via read_file / notebook_read.
    image_upload = widgets.FileUpload(
        multiple=True,
        description="📎 Files",
        layout=widgets.Layout(width="90px", height="80px"),
    )
    image_upload.add_class("delfin-agent-upload")
    state["_pending_images"] = []

    # Hard size cap per file to keep dashboard memory bounded. 32 MB is
    # generous for code, configs, PDFs, and small datasets but blocks
    # accidental drops of huge archives.
    _UPLOAD_SIZE_CAP = 32 * 1024 * 1024

    def _on_image_upload(change):
        """Save uploaded files to <workspace>/.delfin/uploads/ for the agent.

        Drops the file under the workspace and lets the agent read it via
        read_file / notebook_read / find_definition / etc. The agent
        learns about the upload via a system message inserted on the
        next send.
        """
        files = image_upload.value
        eng = state.get("engine")
        ws = None
        if eng is not None:
            kp = getattr(eng, "kit_permissions", None)
            if kp is not None:
                ws = kp.workspace
        if ws is None:
            ws = ctx.repo_dir or Path.cwd()
        upload_dir = Path(ws) / ".delfin" / "uploads"
        upload_dir.mkdir(parents=True, exist_ok=True)
        saved: list[Path] = []
        items = (files.items() if isinstance(files, dict)
                 else [(f["name"], f) for f in files])
        for fname, fmeta in items:
            content = fmeta.get("content") or fmeta.get("data")
            if isinstance(content, memoryview):
                content = bytes(content)
            if not content:
                _append_system_message(f"File empty, skipped: {fname}")
                continue
            if len(content) > _UPLOAD_SIZE_CAP:
                _append_system_message(
                    f"File too large (>{_UPLOAD_SIZE_CAP // (1024*1024)} MB): "
                    f"{fname}"
                )
                continue
            try:
                target = upload_dir / Path(fname).name
                target.write_bytes(content)
                saved.append(target)
            except OSError as exc:
                _append_system_message(f"Save failed: {exc}")
        state["_pending_images"] = saved
        if saved:
            paths = "\n".join(f"  - {p}" for p in saved)
            _append_system_message(
                f"📎 {len(saved)} file(s) saved — will be referenced in "
                f"the next message:\n{paths}"
            )

    image_upload.observe(_on_image_upload, names="value")

    # Resume-Last-Session button.
    resume_last_btn = widgets.Button(
        description="↩ Last session",
        tooltip="Load the most recent session from ~/.delfin/agent_sessions",
        layout=widgets.Layout(width="160px"),
    )

    def _on_resume_last(_btn):
        try:
            from delfin.agent.session_store import resume_latest
            data = resume_latest(max_age_s=7 * 86_400)
        except Exception as exc:
            _append_system_message(f"Resume failed: {exc}")
            return
        if data is None:
            _append_system_message(
                "No recent session found (older than 7 days or empty)."
            )
            return
        sid = str(data.get("session_id", ""))
        if not sid:
            _append_system_message("Session file is corrupt (no session_id).")
            return
        try:
            _load_saved_session(sid)
        except Exception as exc:
            _append_system_message(f"Session load failed: {exc}")

    resume_last_btn.on_click(_on_resume_last)

    # Phase-5 wiring helper: bind callbacks on the engine's permissions
    # object every time _ensure_engine produces a new engine.
    def _wire_phase5_callbacks(engine):
        if engine is None:
            return
        kp = getattr(engine, "kit_permissions", None)
        if kp is None:
            return
        try:
            _ensure_task_session_id(engine, create=False)
            kp.ask_user_callback = _show_ask_user_modal
            kp.plan_approval_callback = _show_plan_approval
        except Exception:
            pass
        # Scheduler fire-callback: when a wake-up triggers, drop the prompt
        # into the input box and send. The thread runs in the scheduler's
        # background, so we marshal back to the UI thread via input_textarea.
        try:
            from delfin.agent import scheduler as _sched_mod
            sch = _sched_mod.get_scheduler()

            def _on_wake(entry):
                try:
                    input_textarea.value = (
                        f"[scheduled] {entry.reason or entry.prompt}\n\n"
                        f"{entry.prompt}"
                    )
                    _on_send(None)
                except Exception:
                    pass

            sch.set_fire_callback(_on_wake)
        except Exception:
            pass
        _refresh_task_ticker()
        _refresh_status_line()

    state["_wire_phase5_callbacks"] = _wire_phase5_callbacks

    # -- layout assembly ---------------------------------------------------
    # Augment the existing controls row with the resume button.
    try:
        controls_row.children = tuple(controls_row.children) + (resume_last_btn,)
    except Exception:
        pass

    agent_content = widgets.VBox(
        [css_widget, _enter_js_output, controls_row, session_row, search_row,
         status_html, cycle_inspector_html, inspector_actions_row, inspector_detail_box,
         kit_mode_row, kit_dirs_status, kit_confirm_container,
         task_ticker_html,
         todo_pane_html, subagent_pane_html,
         chat_html,
         plan_accept_btn, ask_user_box,
         working_html, queue_html, context_bar_html,
         approval_row, action_confirm_row, question_row,
         palette_row, palette_select,
         widgets.HBox(
             [image_upload, input_row],
             layout=widgets.Layout(
                 width="100%",
                 align_items="flex-start",
                 margin="6px 0 0 0",
             ),
         ),
         status_line_html],
    )

    if not _yaml_ok:
        missing_html = widgets.HTML(
            value=(
                '<div class="delfin-agent-nokey">'
                "<h3>DELFIN Agent</h3>"
                "<p>Required package <code>pyyaml</code> not installed.</p>"
                "<p>Run: <code>pip install pyyaml</code></p>"
                "</div>"
            ),
        )
        tab_widget = widgets.VBox([css_widget, missing_html])
        return tab_widget, {}

    # No backend block removed — tab always loads, shows hint on send if needed

    # -- session helpers ---------------------------------------------------

    def _refresh_session_dropdown():
        """Rebuild the session dropdown from saved sessions."""
        try:
            from delfin.agent.session_store import list_sessions
            sessions = list_sessions(limit=30)
        except Exception:
            sessions = []

        options = [("+ New Session", "")]
        for s in sessions:
            title = s.get("title", "Untitled") or "Untitled"
            if len(title) > 50:
                title = title[:50] + "..."
            mode = s.get("mode", "")
            n_msgs = s.get("message_count", 0)
            sid = s.get("session_id", "")
            # Format: "title (mode, N msgs)"
            label = f"{title}  [{mode}, {n_msgs} msgs]"
            options.append((label, sid))

        session_dropdown.options = options
        # Keep current selection if still valid
        active = state.get("active_session_id", "")
        valid_ids = [v for _, v in options]
        if active and active in valid_ids:
            session_dropdown.value = active
        else:
            session_dropdown.value = ""

    def _auto_save_session():
        """Save the current session state to disk."""
        engine = state["engine"]
        if not engine:
            return
        # The Claude CLI populates engine.session_id from its stream events,
        # but the OpenAI / KIT-Toolbox path emits no such event so the field
        # stays empty and the session would silently never persist. Mint a
        # UUID ourselves so auto-save works for every provider.
        if not engine.session_id:
            state_sid = str(state.get("active_session_id", "") or "").strip()
            if state_sid:
                engine.session_id = state_sid
            else:
                import uuid as _uuid
                engine.session_id = str(_uuid.uuid4())
        # Skip if there's nothing to save yet (no chat messages).
        if not state.get("chat_messages"):
            return
        try:
            from delfin.agent.session_store import save_session
            estate = engine.export_state()
            # Capture the full live state so resume restores not just chat
            # history but also UI mode, permissions, the in-flight plan
            # body, the subagent panel snapshot — everything a long session
            # depends on to feel continuous after Ctrl-C / reopen.
            save_session(
                session_id=engine.session_id,
                mode=estate["mode"],
                role_index=estate["role_index"],
                route=estate["route"],
                role_outputs=estate["role_outputs"],
                chat_messages=state["chat_messages"],
                cycle_history=state.get("_cycle_history", []),
                engine_messages=estate["engine_messages"],
                token_usage=estate["token_usage"],
                cost_usd=estate["cost_usd"],
                perm_profile=state.get("_perm_profile", ""),
                provider=provider_dropdown.value,
                model=model_dropdown.value,
                effort=effort_dropdown.value,
                active_gate=state.get("_active_gate"),
                last_compaction_info=getattr(engine, "last_compaction_info", None),
                subagent_calls=state.get("subagent_calls") or [],
                pending_plan_body=state.get("_pending_plan_body", ""),
                todo_payload=state.get("current_todos") or [],
            )
            state["active_session_id"] = engine.session_id
            try:
                kp = getattr(engine, "kit_permissions", None)
                if kp is not None:
                    kp.task_session_id = engine.session_id
            except Exception:
                pass
            _refresh_session_dropdown()
        except Exception:
            pass  # non-critical — don't break the chat

    def _load_saved_session(session_id):
        """Load a saved session and restore engine + UI state."""
        try:
            from delfin.agent.session_store import load_session
            data = load_session(session_id)
        except Exception:
            data = None
        if not data:
            _append_system_message(f"Session not found: {session_id[:12]}...")
            return

        # Restore or create engine with the saved mode (migrate legacy names)
        _legacy_map = {"default": "quick", "high_risk": "reviewed",
                       "runtime_cluster": "cluster", "release": "full"}
        saved_mode = data.get("mode", "quick")
        saved_mode = _legacy_map.get(saved_mode, saved_mode)
        _set_mode_programmatically(saved_mode)

        engine = _ensure_engine()
        if not engine:
            return

        # Restore engine state
        engine.restore_state({
            "mode": saved_mode,
            "role_index": data.get("role_index", 0),
            "role_outputs": data.get("role_outputs", {}),
            "engine_messages": data.get("engine_messages", []),
            "token_usage": data.get("token_usage", {"input": 0, "output": 0}),
            "cost_usd": data.get("cost_usd", 0.0),
            "session_id": session_id,
        })

        # Restore chat UI
        state["chat_messages"] = data.get("chat_messages", [])
        state["_cycle_history"] = data.get("cycle_history", [])
        state["_mode_manual_override"] = True
        state["active_session_id"] = session_id
        try:
            kp = getattr(engine, "kit_permissions", None)
            if kp is not None:
                kp.task_session_id = session_id
        except Exception:
            pass

        # Long-session state restore. Each block is wrapped in try/except
        # so a missing dropdown option or a legacy field never breaks the
        # restore — we always end up with a usable engine + chat.
        saved_perm = data.get("perm_profile") or ""
        if saved_perm:
            try:
                state["_perm_profile"] = saved_perm
                if perm_dropdown.value != saved_perm:
                    perm_dropdown.value = saved_perm
            except Exception:
                pass
        saved_provider = data.get("provider") or ""
        if saved_provider:
            try:
                if provider_dropdown.value != saved_provider:
                    provider_dropdown.value = saved_provider
            except Exception:
                pass
        saved_model = data.get("model") or ""
        if saved_model:
            try:
                valid_models = {v for _, v in (model_dropdown.options or [])}
                if saved_model in valid_models and model_dropdown.value != saved_model:
                    model_dropdown.value = saved_model
            except Exception:
                pass
        saved_effort = data.get("effort") or ""
        if saved_effort:
            try:
                if effort_dropdown.value != saved_effort:
                    effort_dropdown.value = saved_effort
            except Exception:
                pass
        # Subagent panel + todo display
        sa_calls = data.get("subagent_calls") or []
        if sa_calls:
            state["subagent_calls"] = sa_calls
        todo_payload = data.get("todo_payload") or []
        if todo_payload:
            state["current_todos"] = todo_payload
        # last_compaction_info goes back on the engine so /context shows
        # accurate "last compaction" info after resume.
        lci = data.get("last_compaction_info")
        if lci and hasattr(engine, "last_compaction_info"):
            engine.last_compaction_info = lci
        # Pending plan body — if the previous session was awaiting plan
        # approval, surface a note so the user knows what's in flight.
        pending_plan = data.get("pending_plan_body") or ""
        if pending_plan:
            state["_pending_plan_body"] = pending_plan
            _append_system_message(
                "📋 A pending plan was restored from the saved session. "
                "Switch to /mode plan to review + approve it."
            )
        # Active gate — restoring lets a paused approval keep its banner.
        saved_gate = data.get("active_gate")
        if isinstance(saved_gate, dict) and saved_gate.get("type"):
            try:
                _set_active_gate(
                    saved_gate.get("type", ""),
                    saved_gate.get("role", ""),
                    saved_gate.get("title", ""),
                    saved_gate.get("detail", ""),
                )
            except Exception:
                _set_active_gate()
        else:
            _set_active_gate()

        _refresh_chat_html()
        _refresh_task_ticker()
        _update_status()
        _update_button_states()

        title = data.get("title", "")[:50] or session_id[:12]
        _record_cycle_event("session", "Session restored", title)
        _append_system_message(f"Session restored: {title}")

    # -- helpers -----------------------------------------------------------

    def _get_agent_settings():
        """Load agent settings from user settings."""
        try:
            from delfin.user_settings import load_settings
            return load_settings().get("agent", {}) or {}
        except Exception:
            return {}

    def _set_mode_programmatically(new_mode: str):
        state["_mode_change_internal"] = True
        try:
            mode_dropdown.value = new_mode
        finally:
            state["_mode_change_internal"] = False

    def _ensure_task_session_id(engine=None, *, create: bool = False) -> str:
        """Return/sync the current task session id for dashboard tasks."""
        sid = str(state.get("active_session_id", "") or "").strip()
        if not sid and engine is not None:
            sid = str(getattr(engine, "session_id", "") or "").strip()
            if sid:
                state["active_session_id"] = sid
        if not sid and create:
            import uuid as _uuid
            sid = str(_uuid.uuid4())
            state["active_session_id"] = sid
        try:
            kp = getattr(engine, "kit_permissions", None)
            if kp is not None:
                kp.task_session_id = sid
        except Exception:
            pass
        return sid

    def _dropdown_values(options) -> list[str]:
        """Normalize ipywidgets dropdown options to a list of values.

        ipywidgets accepts either plain strings like ["dashboard", "solo"]
        or (label, value) tuples. Slash-command handlers must support both.
        """
        vals: list[str] = []
        for opt in options or []:
            if isinstance(opt, (tuple, list)) and len(opt) >= 2:
                vals.append(str(opt[1]))
            else:
                vals.append(str(opt))
        return vals

    def _resolve_backend():
        """Determine which backend to use: cli or api."""
        if provider_dropdown.value == "kit":
            return "api"  # KIT Toolbox is API-only
        if provider_dropdown.value == "openai":
            settings = _get_agent_settings()
            preferred = settings.get("backend", "cli")
            if preferred == "cli" and _codex_cli_available:
                return "cli"
            if preferred == "api" and os.environ.get("OPENAI_API_KEY", ""):
                return "api"
            # Fallback: CLI first, then API
            if _codex_cli_available:
                return "cli"
            if os.environ.get("OPENAI_API_KEY", ""):
                return "api"
            return "cli"  # will error at runtime with helpful message
        settings = _get_agent_settings()
        preferred = settings.get("backend", "cli")
        if preferred == "cli" and _cli_available:
            return "cli"
        if preferred == "api":
            if os.environ.get("ANTHROPIC_API_KEY", ""):
                return "api"
        # Fallback: try CLI first, then API
        if _cli_available:
            return "cli"
        if os.environ.get("ANTHROPIC_API_KEY", ""):
            return "api"
        return "cli"  # will error at runtime

    def _ensure_engine():
        """Create or re-use the engine."""
        if state["engine"] is not None:
            return state["engine"]

        settings = _get_agent_settings()
        provider = provider_dropdown.value
        backend = _resolve_backend()
        if provider == "kit":
            api_key = os.environ.get("KIT_TOOLBOX_API_KEY", "")
        elif provider == "openai":
            api_key = os.environ.get("OPENAI_API_KEY", "")
        else:
            api_key = os.environ.get("ANTHROPIC_API_KEY", "")
        model = model_dropdown.value or settings.get("model", "")

        try:
            from delfin.agent.engine import AgentEngine

            repo_dir = ctx.repo_dir or Path.cwd()
            # MCP config from agent settings (optional)
            _agent_s = _get_agent_settings()
            _mcp_cfg = _agent_s.get("mcp_config", "")

            # Auto-inject doc server MCP config if docs are enabled
            try:
                from delfin.doc_server.config import ensure_mcp_config as _ensure_docs_mcp
                _mcp_cfg = _ensure_docs_mcp(_mcp_cfg)
            except ImportError:
                pass  # doc_server or mcp not installed

            # Auto-inject ops server MCP config (typed DELFIN workflow tools)
            try:
                from delfin.ops_server.config import ensure_mcp_config as _ensure_ops_mcp
                _mcp_cfg = _ensure_ops_mcp(_mcp_cfg, workspace=str(repo_dir))
            except ImportError:
                pass  # ops_server or mcp not installed

            # CLI tool configuration per mode and permission profile
            _cli_tools = None
            _extra_dirs = None
            _ws_dir = str(ctx.agent_dir) if ctx.agent_dir else ""
            if mode_dropdown.value == "dashboard":
                _cli_tools = [
                    "Read", "Grep", "Glob",      # read code + data
                    "Write", "Bash",              # agent_workspace only
                    "WebSearch", "WebFetch",      # literature research
                ]
                if _ws_dir:
                    _extra_dirs = [_ws_dir]
            elif state.get("_perm_profile") == "repo_free":
                # repo_free + non-dashboard: auto-approve safe bash commands
                # so the agent doesn't get stuck in permission-denial loops
                # for routine operations (git, pytest, syntax checks).
                _cli_tools = [
                    "Read", "Grep", "Glob", "Write", "Edit", "Bash",
                    "WebSearch", "WebFetch",
                    # Safe bash auto-approve patterns
                    "Bash(git *)",
                    "Bash(python -m pytest*)", "Bash(python3 -m pytest*)",
                    "Bash(python -c *)", "Bash(python3 -c *)",
                    "Bash(python -m py_compile*)", "Bash(python3 -m py_compile*)",
                ]

            # KIT-Toolbox: build/refresh the confirmation broker BEFORE the
            # engine, so its callback is wired in at construction time. Other
            # providers leave the panel hidden.
            _kit_callback = None
            if provider == "kit":
                broker = _ensure_kit_broker()
                _kit_callback = broker.callback if broker is not None else None
                _show_kit_confirm_panel(True)
            else:
                _show_kit_confirm_panel(False)

            engine = AgentEngine(
                repo_dir=repo_dir,
                backend=backend,
                provider=provider,
                api_key=api_key,
                model=model,
                mode=mode_dropdown.value,
                permission_mode=_active_cli_perm(),
                mcp_config=_mcp_cfg,
                allowed_tools=_cli_tools,
                extra_dirs=_extra_dirs,
                agent_workspace_dir=_ws_dir,
                effort=str(effort_dropdown.value or ""),
                kit_confirm_callback=_kit_callback,
            )

            # Configure calc search directories for OpenAI function calling
            try:
                from delfin.agent.api_client import _doc_executor
                _doc_executor._calc_dirs = {
                    "calc": str(ctx.calc_dir),
                    "archive": str(ctx.archive_dir),
                    "remote_archive": ctx.runtime_settings.get(
                        "remote_archive_dir", ""
                    ),
                }
            except Exception:
                pass

            # Defense in depth: if for any reason the callback wasn't bound
            # at construction (e.g. broker built later), bind at runtime.
            if provider == "kit" and _kit_callback is not None:
                try:
                    engine.set_kit_confirm_callback(_kit_callback)
                except Exception:
                    pass

            state["engine"] = engine
            ctx.agent_engine = engine
            # Phase 5 wiring: bind ask_user / plan-approval / scheduler
            # callbacks now that the engine + perms exist.
            try:
                _wire_p5 = state.get("_wire_phase5_callbacks")
                if _wire_p5 is not None:
                    _wire_p5(engine)
            except Exception:
                pass
            return engine
        except Exception as exc:
            _append_system_message(f"Engine error: {exc}")
            return None

    def _append_chat_message(role, content, role_label="", **meta):
        payload = {"role": role, "content": content, "role_label": role_label}
        payload.update(meta)
        state["chat_messages"].append(payload)
        # Debounce: during streaming, batch tool messages to reduce jitter.
        # Only refresh if >0.2s since last refresh or if not streaming.
        now = time.monotonic()
        if not state.get("streaming") or now - state.get("_last_html_refresh", 0) > 0.2:
            _refresh_chat_html()
            state["_last_html_refresh"] = now
        else:
            # Schedule a delayed refresh so batched messages still appear
            state["_html_refresh_pending"] = True

    def _append_system_message(text):
        _append_chat_message("system", text)

    def _append_tool_message(html_content: str):
        """Append a tool-call message with terminal-style formatting.

        Content is pre-formatted HTML (not markdown) using the
        .tool-name, .tool-path, .tool-param CSS classes.
        """
        _append_chat_message("tool", html_content)

    def _show_action_confirmation(agent_text: str, commands: list[str]) -> None:
        """D1: surface Approve/Deny buttons for a list of pending ACTIONs.

        Stores the originating agent text so Approve can re-run
        ``_dashboard_auto_exec`` without re-prompting the model.
        """
        action_confirm_html.value = _render_action_confirmation_html(commands)
        action_confirm_row.layout.display = "flex"
        state["_pending_action_text"] = agent_text
        state["_pending_action_commands"] = list(commands)

    def _hide_action_confirmation() -> None:
        action_confirm_html.value = ""
        action_confirm_row.layout.display = "none"
        state["_pending_action_text"] = ""
        state["_pending_action_commands"] = []

    def _on_actions_approve(_btn=None) -> None:
        text = state.get("_pending_action_text") or ""
        _hide_action_confirmation()
        if not text:
            return
        try:
            results = _dashboard_auto_exec(text, force_no_confirm=True)
        except Exception as exc:
            _append_system_message(f"Action execution failed: {exc}")
            return
        if results:
            _append_system_message(
                f"User approved → executed {len(results)} action(s)."
            )

    def _on_actions_deny(_btn=None) -> None:
        commands = list(state.get("_pending_action_commands") or [])
        _hide_action_confirmation()
        _append_system_message(
            "User declined → actions NOT executed: "
            + (", ".join(commands[:3]) + ("..." if len(commands) > 3 else ""))
        )
        # Feed the refusal back to the agent so it knows to plan differently.
        try:
            input_textarea.value = "Bitte die Aktionen NICHT ausführen — anderen Vorschlag machen."
        except Exception:
            pass

    def _on_dashboard_tab_change(change=None) -> None:
        """D2: when the user switches tabs, surface a one-shot suggestion
        from the agent perspective.  Triggers at most once per tab per
        session (tracked in ``state["_seen_tab_suggestions"]``)."""
        tabs = getattr(ctx, "tabs_widget", None)
        if tabs is None or not ctx.tab_indices:
            return
        try:
            idx = tabs.selected_index
            tab_name = next(
                (t for t, i in ctx.tab_indices.items() if i == idx),
                "",
            )
        except Exception:
            return
        if not tab_name:
            return
        seen = state.setdefault("_seen_tab_suggestions", set())
        if tab_name in seen:
            return
        suggestion = _suggestion_for_tab(tab_name)
        if not suggestion:
            return
        seen.add(tab_name)
        # Only inject if the engine is initialised (otherwise we'd queue
        # messages that confuse a fresh session).
        if state.get("engine") is None:
            return
        _append_system_message(f"💡 [{tab_name}] {suggestion}")

    def _check_job_events_once() -> int:
        """D5: poll the backend, diff against last-seen statuses, and
        inject system messages for notable transitions.

        Returns the number of system messages emitted (0 when nothing
        changed or no backend).  Safe to call from the UI thread.
        """
        from delfin.dashboard.job_events import (
            diff_job_states, format_job_event_message,
            is_notable_event, snapshot_jobs,
        )
        backend = getattr(ctx, "backend", None)
        if backend is None:
            return 0
        try:
            jobs = backend.list_jobs()
        except Exception:
            return 0
        current = snapshot_jobs(jobs)
        previous = state.get("_job_states_seen") or {}
        events = diff_job_states(previous, current)
        emitted = 0
        for event in events:
            if not is_notable_event(event):
                continue
            try:
                msg = format_job_event_message(event)
            except Exception:
                continue
            _append_system_message(msg)
            emitted += 1
        # Update snapshot to a tiny status-only dict so we don't pin
        # JobInfo references in memory.
        state["_job_states_seen"] = {
            jid: getattr(info, "status", "") or info.get("status", "")
            if isinstance(info, dict) else getattr(info, "status", "")
            for jid, info in current.items()
        }
        return emitted

    def _start_job_event_watcher(interval: float = 30.0) -> None:
        """Launch a daemon thread that calls _check_job_events_once()
        every ``interval`` seconds. Idempotent — a second call is a no-op
        while a watcher is already running."""
        if state.get("_job_event_thread") is not None:
            return
        if getattr(ctx, "backend", None) is None:
            return
        state["_job_event_thread_stop"] = False

        def _loop() -> None:
            while not state.get("_job_event_thread_stop"):
                try:
                    _check_job_events_once()
                except Exception:
                    pass
                # Cooperative sleep so stop is responsive
                slept = 0.0
                while slept < interval and not state.get(
                    "_job_event_thread_stop"
                ):
                    time.sleep(0.5)
                    slept += 0.5
            state["_job_event_thread"] = None

        t = threading.Thread(target=_loop, daemon=True, name="delfin-jobs-watch")
        state["_job_event_thread"] = t
        t.start()

    def _stop_job_event_watcher() -> None:
        """Signal the watcher loop to exit; safe to call multiple times."""
        state["_job_event_thread_stop"] = True

    def _record_cycle_event(kind: str, title: str, detail: str = "", role: str = ""):
        history = state.setdefault("_cycle_history", [])
        entry = {
            "kind": (kind or "event").strip(),
            "title": (title or "").strip(),
            "detail": (detail or "").strip(),
            "role": (role or "").strip(),
        }
        if not entry["title"]:
            return
        if history and all(
            history[-1].get(key, "") == entry.get(key, "")
            for key in ("kind", "title", "detail", "role")
        ):
            return
        history.append(entry)
        if len(history) > 12:
            del history[:-12]

    def _set_active_gate(gate_type="", role="", title="", detail=""):
        if gate_type:
            state["_active_gate"] = {
                "type": gate_type,
                "role": role,
                "title": title,
                "detail": detail,
            }
            # Spinner reflects the wait — orange paused-state so the user
            # always knows the agent is blocked on them, not stuck silently.
            if not state.get("streaming"):
                label = title or _gate_label(gate_type)
                try:
                    _set_working(True, f"Waiting: {label}", mode="gated")
                except Exception:
                    pass
        else:
            state["_active_gate"] = None
            # Don't clear the spinner here — the caller decides; otherwise
            # we'd clobber a streaming spinner when a gate gets dismissed.

    def _gate_label(gate_type: str) -> str:
        labels = {
            "schema": "Schema Gate",
            "risk": "Risk Gate",
            "partial": "Completion Gate",
            "goal-lock": "Goal Lock Gate",
            "review": "Review Gate",
            "question": "Question Gate",
            "plan-approval": "Plan Approval",
            "findings": "Findings Gate",
            "approval": "Permission Gate",
            "confidence": "Confidence Gate",
            "conflict": "Conflict Gate",
            "cost": "Cost Gate",
        }
        return labels.get(gate_type or "", "Gate")

    def _compact_text(text: str, limit: int = 140) -> str:
        clean = " ".join((text or "").strip().split())
        if len(clean) <= limit:
            return clean
        return clean[: max(0, limit - 1)].rstrip() + "…"

    def _extract_open_risks(text: str) -> list[str]:
        risks: list[str] = []
        in_section = False
        for line in (text or "").splitlines():
            stripped = line.strip()
            lowered = stripped.lower()
            if lowered.startswith("### open risks"):
                in_section = True
                continue
            if lowered.startswith("**open risks:**"):
                in_section = True
                inline = stripped.split(":", 1)[1].strip() if ":" in stripped else ""
                if inline and inline.lower() not in {"none", "none.", "n/a", "no", "no risks"}:
                    risks.append(inline)
                continue
            if in_section:
                if stripped.startswith("###") or stripped.startswith("## ") or stripped.startswith("**"):
                    break
                match = re.match(r"^(?:\d+\.\s*|-\s*)(.*)", stripped)
                if match and match.group(1).strip():
                    item = match.group(1).strip()
                    if item.lower() not in {"none", "none.", "n/a", "no", "no risks"}:
                        risks.append(item)
        return risks

    def _render_cycle_inspector() -> str:
        engine = state.get("engine")
        active_gate = state.get("_active_gate") or {}
        goal_text = "Noch kein gelocktes Ziel."
        gate_value = "Kein aktives Gate"
        gate_note = "Die Pipeline kann ohne Benutzerfreigabe weiterlaufen."
        next_role_value = "Session Manager"
        next_role_note = "Startet den nächsten Cycle."
        risks: list[str] = []
        retry_parts: list[str] = []
        history = list(state.get("_cycle_history", []))
        header_meta = "Noch kein aktiver Cycle"

        if engine:
            try:
                from delfin.agent.engine import AgentEngine
            except Exception:
                AgentEngine = None  # type: ignore[assignment]

            sm_output = engine.role_outputs.get("session_manager", "")
            if not sm_output and engine.current_role == "session_manager":
                for msg in reversed(engine.messages):
                    if msg.get("role") == "assistant":
                        sm_output = msg.get("content", "")
                        break

            if AgentEngine and sm_output:
                goal_items = AgentEngine._extract_list_section(sm_output, ("### Goal lock",))
                if goal_items:
                    goal_text = _compact_text(goal_items[0], 180)
                else:
                    task_text = AgentEngine.extract_plan_field(sm_output, "Task")
                    if task_text:
                        goal_text = _compact_text(task_text, 180)

            if active_gate:
                gate_value = active_gate.get("title") or _gate_label(active_gate.get("type", ""))
                gate_note = active_gate.get("detail") or "Wartet auf Benutzerentscheidung."
                next_role_value = "User Input Required"
                paused_role = active_gate.get("role", "")
                next_role_note = (
                    f"Pausiert nach {_format_role_label(paused_role)}."
                    if paused_role else "Pausiert bis zur nächsten Freigabe."
                )
            elif engine.is_cycle_complete:
                next_role_value = "Cycle Complete"
                next_role_note = "Kein weiterer Rollenschritt offen."
            elif engine.current_role:
                next_role_value = _format_role_label(engine.current_role)
                upcoming = ""
                if engine.current_role_index + 1 < len(engine.route):
                    upcoming = _format_role_label(engine.route[engine.current_role_index + 1])
                next_role_note = (
                    f"Danach: {upcoming}" if upcoming else "Danach: Abschluss des Cycles."
                )

            roles_to_scan = list(reversed(engine.route))
            seen_roles: set[str] = set()
            for role_id in roles_to_scan:
                if role_id in seen_roles:
                    continue
                seen_roles.add(role_id)
                out = engine.role_outputs.get(role_id, "")
                if not out:
                    continue
                for item in _extract_open_risks(out):
                    if item not in risks:
                        risks.append(item)
                if len(risks) >= 3:
                    break

            if active_gate.get("type") == "risk" and active_gate.get("detail"):
                gate_risk = active_gate["detail"].strip()
                if gate_risk and gate_risk not in risks:
                    risks.insert(0, gate_risk)

            builder_retries = int(state.get("_builder_retries", 0) or 0)
            if builder_retries:
                retry_parts.append(f"Builder {builder_retries}")
            schema_counts = state.get("_schema_retry_counts", {}) or {}
            schema_total = sum(int(v) for v in schema_counts.values() if str(v).isdigit())
            if schema_total:
                retry_parts.append(f"Schema {schema_total}")

            role_total = len(engine.route)
            if engine.is_cycle_complete:
                header_meta = f"{engine.mode} · cycle complete · ${engine.cost_usd:.2f}"
            else:
                header_meta = (
                    f"{engine.mode} · step {engine.current_role_index + 1}/{max(role_total, 1)}"
                    f" · {_format_role_label(engine.current_role)}"
                )
                if state.get("streaming"):
                    header_meta += " · running"
                header_meta += f" · ${engine.cost_usd:.2f}"

        retry_value = "0"
        retry_note = "Keine aktiven Retries."
        if retry_parts:
            retry_value = " · ".join(retry_parts)
            retry_note = "Aktive Wiederholungen im aktuellen Cycle."

        gate_type = active_gate.get("type", "")
        gate_pill = _gate_label(gate_type) if gate_type else "Open"
        gate_pill_html = (
            f'<span class="gate-pill">{_html.escape(gate_pill)}</span>'
            if active_gate else ""
        )

        risk_body = (
            f'<ul class="card-list">{"".join(f"<li>{_html.escape(_compact_text(risk, 120))}</li>" for risk in risks[:3])}</ul>'
            if risks else
            '<span class="card-note">Keine offenen Risiken aus den letzten Agent-Ausgaben.</span>'
        )
        history_items = history[-4:]
        history_body = (
            '<ul class="history-list">' +
            "".join(
                '<li class="history-item">'
                f'<span class="history-kind">{_html.escape(item.get("kind", "event"))}</span>'
                f'<span class="history-title">{_html.escape(_compact_text(item.get("title", ""), 90))}</span>'
                + (
                    f'<div class="history-detail">{_html.escape(_compact_text(item.get("detail", ""), 140))}</div>'
                    if item.get("detail") else ""
                )
                + '</li>'
                for item in reversed(history_items)
            ) +
            '</ul>'
            if history_items else
            '<span class="card-note">Noch keine Gates, Handoffs oder Retries im aktuellen Cycle.</span>'
        )

        return (
            '<div class="delfin-cycle-inspector">'
            '<div class="inspector-header">'
            '<span class="inspector-title">Cycle Inspector</span>'
            f'<span class="inspector-meta">{_html.escape(header_meta)}</span>'
            '</div>'
            '<div class="inspector-grid">'
            '<div class="inspector-card">'
            '<span class="card-label">Locked Goal</span>'
            f'<span class="card-value">{_html.escape(goal_text)}</span>'
            '<span class="card-note">Gelocktes Ziel aus dem Session-Manager-Plan.</span>'
            '</div>'
            '<div class="inspector-card gate-card">'
            '<span class="card-label">Current Gate</span>'
            f'{gate_pill_html}'
            f'<span class="card-value">{_html.escape(_compact_text(gate_value, 90))}</span>'
            f'<span class="card-note">{_html.escape(_compact_text(gate_note, 160))}</span>'
            '</div>'
            '<div class="inspector-card">'
            '<span class="card-label">Next Role</span>'
            f'<span class="card-value">{_html.escape(next_role_value)}</span>'
            f'<span class="card-note">{_html.escape(next_role_note)}</span>'
            '</div>'
            '<div class="inspector-card risk-card">'
            '<span class="card-label">Open Risks</span>'
            f'<span class="card-value">{_html.escape(str(len(risks)) if risks else "0")}</span>'
            f'{risk_body}'
            '</div>'
            '<div class="inspector-card retry-card">'
            '<span class="card-label">Retry Count</span>'
            f'<span class="card-value">{_html.escape(retry_value)}</span>'
            f'<span class="card-note">{_html.escape(retry_note)}</span>'
            '</div>'
            '<div class="inspector-card history-card">'
            '<span class="card-label">Recent Activity</span>'
            f'<span class="card-value">{_html.escape(str(len(history)) if history else "0")}</span>'
            f'{history_body}'
            '</div>'
            '</div>'
            '</div>'
        )

    def _build_inspector_detail_entries() -> list[dict[str, str]]:
        entries: list[dict[str, str]] = []
        engine = state.get("engine")
        active_gate = state.get("_active_gate") or {}
        history = list(state.get("_cycle_history", []))

        if active_gate:
            entries.append({
                "key": "gate:active",
                "label": f"Active Gate · {_compact_text(active_gate.get('title', ''), 50)}",
                "title": active_gate.get("title") or _gate_label(active_gate.get("type", "")),
                "meta": " · ".join(
                    part for part in [
                        _gate_label(active_gate.get("type", "")),
                        _format_role_label(active_gate.get("role", "")) if active_gate.get("role") else "",
                    ] if part
                ),
                "body": active_gate.get("detail") or "Kein weiteres Detail hinterlegt.",
            })

        for idx, item in enumerate(reversed(history[-6:]), 1):
            entries.append({
                "key": f"history:{idx}",
                "label": f"History · {_compact_text(item.get('title', ''), 54)}",
                "title": item.get("title", "History"),
                "meta": " · ".join(
                    part for part in [
                        item.get("kind", "event").upper(),
                        _format_role_label(item.get("role", "")) if item.get("role") else "",
                    ] if part
                ),
                "body": item.get("detail") or "Kein weiteres Detail hinterlegt.",
            })

        if engine:
            role_order = list(reversed(engine.route))
            seen: set[str] = set()
            for role_id in role_order:
                if role_id in seen:
                    continue
                seen.add(role_id)
                output = engine.role_outputs.get(role_id, "")
                if not output:
                    continue
                entries.append({
                    "key": f"role:{role_id}",
                    "label": f"Role Output · {_format_role_label(role_id)}",
                    "title": f"{_format_role_label(role_id)} output",
                    "meta": "Role output",
                    "body": output[:4000],
                })
                if len(entries) >= 10:
                    break

            if engine.messages:
                for msg in reversed(engine.messages):
                    if msg.get("role") == "assistant":
                        entries.append({
                            "key": "role:current",
                            "label": f"Current Output · {_format_role_label(engine.current_role)}",
                            "title": f"{_format_role_label(engine.current_role)} current output",
                            "meta": "Current assistant message",
                            "body": msg.get("content", "")[:4000],
                        })
                        break

        return entries[:12]

    def _render_inspector_detail_panel(entry: dict[str, str] | None) -> str:
        if not entry:
            return (
                '<div class="delfin-cycle-detail">'
                '<div class="detail-title">Keine Details verfügbar</div>'
                '<div class="detail-meta">Cycle Inspector</div>'
                '<div class="detail-body">Sobald Gates, History oder Rollen-Outputs vorliegen, erscheinen sie hier.</div>'
                '</div>'
            )
        body_html = _md_to_html(entry.get("body", "") or "Kein weiteres Detail hinterlegt.")
        return (
            '<div class="delfin-cycle-detail">'
            f'<div class="detail-title">{_html.escape(entry.get("title", "Detail"))}</div>'
            f'<div class="detail-meta">{_html.escape(entry.get("meta", ""))}</div>'
            f'<div class="detail-body">{body_html}</div>'
            '</div>'
        )

    def _update_cycle_inspector():
        # Hide Cycle Inspector for dashboard/solo — only relevant for pipelines
        _is_pipeline = mode_dropdown.value not in ("dashboard", "solo")
        _disp = "block" if _is_pipeline else "none"
        for _w in (cycle_inspector_html, inspector_actions_row, inspector_detail_box):
            _w.layout.display = _disp
        if not _is_pipeline:
            return
        cycle_inspector_html.value = _render_cycle_inspector()

    def _update_inspector_detail(force_reset: bool = False):
        entries = _build_inspector_detail_entries()
        options = [(entry["label"], entry["key"]) for entry in entries] or [("Keine Details verfügbar", "")]
        selected = state.get("_inspector_detail_key", "")
        valid_keys = {value for _, value in options}
        if force_reset or selected not in valid_keys:
            selected = options[0][1]
        inspector_detail_dropdown.options = options
        inspector_detail_dropdown.value = selected
        state["_inspector_detail_key"] = selected
        selected_entry = next((entry for entry in entries if entry["key"] == selected), None)
        inspector_detail_html.value = _render_inspector_detail_panel(selected_entry)

    def _update_inspector_actions():
        engine = state.get("engine")
        active_gate = state.get("_active_gate") or {}
        is_streaming = bool(state.get("streaming"))
        has_pending_approval = bool(
            state.get("_pending_dashboard_action") or state.get("_pending_approval")
        )
        awaiting_gate = any(
            bool(state.get(key))
            for key in (
                "_awaiting_agent_question",
                "_awaiting_findings_review",
                "_awaiting_gate_review",
                "_awaiting_conflict_resolution",
            )
        )

        inspector_primary_btn.disabled = False
        inspector_retry_btn.disabled = is_streaming or not (engine and engine.messages)
        inspector_stop_btn.disabled = not (is_streaming or engine or active_gate)
        inspector_next_btn.disabled = True

        if has_pending_approval:
            inspector_primary_btn.description = "Approve"
            inspector_primary_btn.tooltip = "Approve the pending permission or dashboard action"
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Eine Freigabe wartet auf Bestätigung.</span>'
            )
        elif active_gate or awaiting_gate:
            inspector_primary_btn.description = "Continue"
            inspector_primary_btn.tooltip = "Continue past the current gate with a go signal"
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Gate aktiv: direkt freigeben, retryen oder stoppen.</span>'
            )
        elif is_streaming:
            inspector_primary_btn.description = "Running"
            inspector_primary_btn.tooltip = "The cycle is currently running"
            inspector_primary_btn.disabled = True
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Agent läuft gerade. Stopptaste beendet den aktuellen Schritt.</span>'
            )
        elif engine and not engine.is_cycle_complete and engine.current_role_index < len(engine.route) - 1 and engine.messages:
            inspector_primary_btn.description = "Continue"
            inspector_primary_btn.tooltip = "Advance to the next role"
            inspector_next_btn.disabled = False
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Cycle bereit für den nächsten Rollenschritt.</span>'
            )
        elif engine and engine.is_cycle_complete:
            inspector_primary_btn.description = "Complete"
            inspector_primary_btn.tooltip = "Cycle is complete"
            inspector_primary_btn.disabled = True
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Cycle abgeschlossen. Neue Session starten oder Follow-up schicken.</span>'
            )
        else:
            inspector_primary_btn.description = "Continue"
            inspector_primary_btn.tooltip = "Start or continue the cycle"
            inspector_primary_btn.disabled = not bool(engine)
            inspector_actions_info.value = (
                '<span class="delfin-cycle-actions-info">'
                'Cycle-Steuerung steht bereit.</span>'
            )

    def _send_control_reply(text: str):
        if state.get("streaming"):
            return
        input_textarea.value = text
        _on_send(None)

    def _on_inspector_detail_change(change):
        if change.get("name") != "value":
            return
        state["_inspector_detail_key"] = change.get("new", "") or ""
        entries = _build_inspector_detail_entries()
        selected_entry = next(
            (entry for entry in entries if entry["key"] == state["_inspector_detail_key"]),
            None,
        )
        inspector_detail_html.value = _render_inspector_detail_panel(selected_entry)

    def _on_inspector_primary(button):
        if state.get("streaming"):
            return
        if state.get("_pending_dashboard_action") or state.get("_pending_approval"):
            _on_approve(button)
            return
        if (state.get("_active_gate")
                or state.get("_awaiting_agent_question")
                or state.get("_awaiting_findings_review")
                or state.get("_awaiting_gate_review")
                or state.get("_awaiting_conflict_resolution")):
            _send_control_reply("go")
            return
        _on_advance_role(button)

    def _on_inspector_retry(button):
        _retry_last()

    def _on_inspector_stop(button):
        if state.get("streaming"):
            _on_stop(button)
            return
        engine = state.get("engine")
        if engine:
            engine.request_stop()
        for key in (
            "_awaiting_agent_question",
            "_awaiting_findings_review",
            "_awaiting_gate_review",
            "_awaiting_conflict_resolution",
        ):
            state.pop(key, None)
        _set_active_gate()
        _record_cycle_event("pause", "Pipeline paused by user")
        _append_system_message("Pipeline paused by user. Use Continue, Next Role or send a new instruction.")
        _update_status()
        _update_button_states()

    def _on_inspector_next(button):
        _on_advance_role(button)

    def _append_gate_message(gate_type, role, title, detail="", hint=""):
        _set_active_gate(gate_type, role, title, detail)
        _record_cycle_event(gate_type or "gate", title, detail, role)
        body = detail.strip()
        if hint:
            body = f"{body}\n\n{hint}".strip()
        _append_chat_message(
            "gate",
            body,
            role_label=_format_role_label(role) if role else "",
            gate_type=gate_type,
            gate_title=title,
            gate_hint=hint,
        )
        _update_status()

    def _update_last_assistant(content, role_label="", finalize=False):
        """Update the last assistant message (for streaming).

        During streaming (*finalize* is False), the message is rendered
        as plain pre-formatted text for speed.  When *finalize* is True,
        the message is converted to full markdown HTML.
        """
        msgs = state["chat_messages"]
        if msgs and msgs[-1]["role"] == "assistant":
            msgs[-1]["content"] = content
            if role_label:
                msgs[-1]["role_label"] = role_label
            msgs[-1]["_streaming"] = not finalize
        else:
            msgs.append(
                {"role": "assistant", "content": content,
                 "role_label": role_label, "_streaming": not finalize}
            )
        _refresh_chat_html(streaming=not finalize)
        if finalize and content.strip():
            # Plan-mode "Akzeptieren"-Button erscheint, sobald der Agent eine
            # Plan-Antwort abgeschlossen hat. Wir merken uns nur den letzten
            # Status — beim Wechsel auf acceptEdits zurücksetzen.
            state["_kit_plan_has_response"] = True
            try:
                _refresh_plan_accept_btn()
            except Exception:
                pass

    # HTML cache: stores rendered HTML for messages that haven't changed.
    # During streaming, only the last (assistant) message is re-rendered.
    _html_cache: dict[str, str] = {"prefix": "", "prefix_count": 0}

    def _render_single_msg(msg):
        """Render a single chat message to HTML."""
        role = msg["role"]
        if role == "user":
            content = _md_to_html(msg["content"])
            return (
                f'<div class="delfin-chat-msg delfin-chat-user">{content}</div>'
            )
        elif role == "assistant":
            label = msg.get("role_label", "Agent")
            # During streaming: plain <pre> — no expensive markdown parsing.
            # On finalize: full _md_to_html with code blocks, formatting etc.
            if msg.get("_streaming"):
                content = (
                    f'<pre class="delfin-streaming-pre">'
                    f'{_html.escape(msg["content"])}</pre>'
                )
            else:
                content = _md_to_html(msg["content"])
            return (
                f'<div class="delfin-chat-msg delfin-chat-agent">'
                f'<div class="delfin-chat-role">{_html.escape(label)}</div>'
                f"{content}</div>"
            )
        elif role == "thinking":
            raw = msg["content"]
            preview = raw[:80].replace("\n", " ").strip()
            if len(raw) > 80:
                preview += "..."
            escaped_full = _html.escape(raw)
            return (
                f'<details class="delfin-chat-msg delfin-chat-thinking">'
                f'<summary>\U0001f9e0 {_html.escape(preview)}</summary>'
                f'<div class="thinking-content">{escaped_full}</div>'
                f'</details>'
            )
        elif role == "approval":
            content = _md_to_html(msg["content"])
            return (
                f'<div class="delfin-chat-msg delfin-chat-approval">'
                f'\u26a0\ufe0f {content}'
                f'<div style="margin-top:4px;font-size:11px;color:#78350f;">'
                f'Press <b>Enter</b> to approve, <b>Esc</b> to deny</div>'
                f'</div>'
            )
        elif role == "gate":
            content = _md_to_html(msg["content"])
            gate_type = msg.get("gate_type", "")
            gate_title = msg.get("gate_title", _gate_label(gate_type))
            gate_hint = msg.get("gate_hint", "")
            gate_role = msg.get("role_label", "")
            card_cls = f"delfin-chat-msg delfin-gate-card delfin-gate-{_html.escape(gate_type or 'review')}"
            badge = _gate_label(gate_type)
            hint_html = (
                f'<div class="gate-hint">{_html.escape(gate_hint)}</div>'
                if gate_hint else ""
            )
            role_html = (
                f'<span class="gate-role">{_html.escape(gate_role)}</span>'
                if gate_role else ""
            )
            return (
                f'<div class="{card_cls}">'
                f'<div class="gate-header">'
                f'<span class="delfin-gate-badge">{_html.escape(badge)}</span>'
                f'<span class="gate-title">{_html.escape(gate_title)}</span>'
                f'{role_html}'
                f'</div>'
                f'<div class="gate-body">{content}</div>'
                f'{hint_html}'
                f'</div>'
            )
        elif role == "tool":
            return (
                f'<div class="delfin-chat-msg delfin-chat-tool">'
                f'{msg["content"]}</div>'
            )
        elif role == "system":
            raw = msg["content"]
            content = _md_to_html(raw)
            if raw.startswith("---") or raw.startswith("Session restored"):
                css_class = "delfin-chat-msg delfin-chat-handoff"
            else:
                css_class = "delfin-chat-msg delfin-chat-system"
            return f'<div class="{css_class}">{content}</div>'
        return ""

    # Always scroll to bottom — simple and reliable.
    _SCROLL_TAG = (
        '<img src="" onerror="'
        "var c=this.closest('.delfin-agent-chat');"
        "if(c)c.scrollTop=c.scrollHeight;"
        "this.remove();"
        '" style="display:none">'
    )

    def _refresh_chat_html(streaming=False):
        """Rebuild the chat display HTML.

        When *streaming* is True, only the last message is re-rendered;
        all earlier messages use a cached HTML prefix.  This avoids
        markdown-converting the entire history on every token chunk.
        """
        msgs = state["chat_messages"]
        if not msgs:
            chat_html.value = (
                '<div class="delfin-agent-chat">'
                "<i>Start a conversation...</i></div>"
            )
            _html_cache["prefix"] = ""
            _html_cache["prefix_count"] = 0
            return

        n = len(msgs)

        # Fast path: during streaming, only re-render the last message.
        # The prefix (all messages before the last) is cached.
        if streaming and n > 1 and _html_cache["prefix_count"] == n - 1:
            last_html = _render_single_msg(msgs[-1])
            chat_html.value = (
                '<div class="delfin-agent-chat">'
                + _html_cache["prefix"]
                + "\n" + last_html
                + "\n" + _SCROLL_TAG
                + "</div>"
            )
            return

        # Full rebuild: render all messages and cache the prefix
        parts = []
        for msg in msgs:
            parts.append(_render_single_msg(msg))

        # Cache everything except the last message
        if n > 1:
            _html_cache["prefix"] = "\n".join(parts[:-1])
            _html_cache["prefix_count"] = n - 1
        else:
            _html_cache["prefix"] = ""
            _html_cache["prefix_count"] = 0

        chat_html.value = (
            '<div class="delfin-agent-chat">'
            + "\n".join(parts)
            + "\n" + _SCROLL_TAG
            + "</div>"
        )

    def _update_status():
        """Update the status bar from engine state."""
        # Live context-window usage bar (cheap, idempotent — fine to call
        # every time the status line refreshes; that already happens after
        # every send/stop/handoff/compact).
        try:
            _refresh_context_bar()
        except Exception:
            pass
        engine = state["engine"]
        active_gate = state.get("_active_gate") or {}
        if engine:
            s = engine.get_status()
            status_html.value = _render_status(
                s["mode"],
                s["backend"],
                s["role"],
                s["role_index"],
                s["role_total"],
                s["input_tokens"],
                s["output_tokens"],
                s["cost_usd"],
                provider=s.get("provider", provider_dropdown.value),
                perm_profile=state.get("_perm_profile", "ask_all"),
                active_gate_type=active_gate.get("type", ""),
                active_gate_text=active_gate.get("title", ""),
            )
        else:
            backend = _resolve_backend() if _cli_available else "api"
            status_html.value = _render_status(
                mode_dropdown.value, backend, "", 0, 0, 0, 0, 0.0,
                provider=provider_dropdown.value,
                perm_profile=state.get("_perm_profile", "ask_all"),
                active_gate_type=active_gate.get("type", ""),
                active_gate_text=active_gate.get("title", ""),
            )
        _update_cycle_inspector()
        # Phase 5e: status-line footer mirrors the same engine snapshot.
        try:
            _refresh_status_line()
        except Exception:
            pass
        _update_inspector_actions()
        _update_inspector_detail()

    def _set_working(active, label="", mode="streaming"):
        """Show or hide the CLI-like activity spinner.

        ``mode`` selects the visual variant:
        - ``streaming`` (default, blue) \u2014 agent is actively streaming
        - ``gated`` (orange) \u2014 paused on user approval / question
        - ``queued`` (grey)  \u2014 idle but messages queued behind a Stop
        - ``stale`` (red)    \u2014 stream silent past the heartbeat threshold
        """
        if active:
            text = label or "Working..."
            _esc = _html.escape

            mode_overrides = {
                "gated":  ("\u23f8\ufe0f", "gated",  "delfin-agent-working--gated"),
                "queued": ("\u23f3",       "queued", "delfin-agent-working--queued"),
                "stale":  ("\u26a0\ufe0f", "stale",  "delfin-agent-working--stale"),
            }
            if mode in mode_overrides:
                icon, cat, mode_cls = mode_overrides[mode]
            else:
                mode_cls = ""
                # Classify streaming activity for icon + label
                _tl = text.lower()
                if _tl.startswith("thinking"):
                    icon, cat = "\U0001f4ad", "think"
                elif "reading" in _tl:
                    icon, cat = "\U0001f4c4", "read"
                elif "editing" in _tl:
                    icon, cat = "\u270f\ufe0f", "edit"
                elif "writing" in _tl:
                    icon, cat = "\u270f\ufe0f", "write"
                elif "searching" in _tl or "finding" in _tl:
                    icon, cat = "\U0001f50d", "grep"
                elif _tl.startswith("$") or _tl.startswith("\u25b8"):
                    icon, cat = "\u25b8", "bash"
                elif "sub-agent" in _tl:
                    icon, cat = "\U0001f916", "agent"
                else:
                    icon, cat = "\u26a1", "work"

            _spinner = (
                '<span class="delfin-spinner">'
                '<span class="delfin-spinner-dot"></span>'
                '<span class="delfin-spinner-dot"></span>'
                '<span class="delfin-spinner-dot"></span>'
                '</span>'
            )

            classes = "delfin-agent-working" + (f" {mode_cls}" if mode_cls else "")
            working_html.value = (
                f'<div class="{classes}">'
                f'{_spinner}'
                f'<span class="delfin-activity-icon">{icon}</span>'
                f'<span class="delfin-activity-label">{_esc(cat)}</span>'
                f'<span class="delfin-activity-text">{_esc(text)}</span>'
                f'</div>'
            )
            # Header colour adapts to mode so the global status chip mirrors
            # the activity-tab spinner exactly.
            header_color = {
                "gated": "#f59e0b", "queued": "#94a3b8", "stale": "#ef4444",
            }.get(mode, "#60a5fa")
            # Global header spinner — visible from all tabs
            short = text[:40] + ("\u2026" if len(text) > 40 else "")
            ctx.agent_status_html.value = (
                f'<span style="font-family:monospace; color:{header_color}; '
                f'padding:2px 6px; background:#1e293b; '
                f'border-radius:4px; font-size:11px;">'
                f'{icon} {_esc(short)}</span>'
            )
        else:
            # Auto-focus input textarea when agent finishes
            working_html.value = (
                '<img src="" onerror="'
                "var ta=document.querySelector('.delfin-agent-input textarea');"
                "if(ta)ta.focus();"
                "this.remove();"
                '" style="display:none">'
            )
            ctx.agent_status_html.value = ""

    def _arm_stale_watcher():
        """Schedule a two-stage watcher that (1) flips the spinner to
        ``mode='stale'`` when the streaming worker hasn't produced
        thinking/tokens/tool-calls for ``state['_stale_threshold_s']``
        seconds and (2) cooperatively stops the stream after a
        configurable kill-threshold so a hung provider doesn't burn
        wall-clock + thinking-tokens indefinitely.

        Re-armed at every send and reset on every activity event.
        Settings:
        - ``agent.stale_threshold_s`` (default 600) — warn-only flip
        - ``agent.stale_kill_after_s`` (default 0 = disabled,
          recommended 120 for dashboard mode) — cooperative
          ``engine.request_stop()`` after this many silent seconds.
        """
        import threading as _threading

        # Cancel any previously armed timer
        prev = state.get("_stale_timer")
        if prev is not None:
            try:
                prev.cancel()
            except Exception:
                pass
        prev_kill = state.get("_stale_kill_timer")
        if prev_kill is not None:
            try:
                prev_kill.cancel()
            except Exception:
                pass

        threshold = float(state.get("_stale_threshold_s") or 600.0)
        # Per-mode default: dashboard is short-action territory, so 120 s
        # is a generous ceiling for "the provider is hung, not thinking".
        # Solo can be long-form reasoning, leave it disabled by default.
        try:
            from delfin.user_settings import load_settings as _ls
            _agent_cfg = (_ls() or {}).get("agent", {}) or {}
            _user_kill = float(_agent_cfg.get("stale_kill_after_s") or 0)
        except Exception:
            _user_kill = 0
        if _user_kill > 0:
            kill_after = _user_kill
        elif state.get("engine") and getattr(state["engine"], "mode", "") == "dashboard":
            kill_after = 120.0
        else:
            kill_after = 0  # disabled

        def _check_stale():
            try:
                if not state.get("streaming"):
                    return
                last = float(state.get("_last_stream_activity") or 0.0)
                if last <= 0.0:
                    return
                elapsed = time.monotonic() - last
                if elapsed >= threshold and not state.get("_stale_seen"):
                    state["_stale_seen"] = True
                    mins = int(elapsed // 60)
                    _set_working(
                        True,
                        f"Stream stalled ({mins} min, no output)",
                        mode="stale",
                    )
            except Exception:
                pass

        def _check_kill():
            """Cooperative kill: send stop signal if still silent."""
            try:
                if not state.get("streaming"):
                    return
                last = float(state.get("_last_stream_activity") or 0.0)
                if last <= 0.0:
                    return
                elapsed = time.monotonic() - last
                if elapsed < kill_after:
                    return
                engine = state.get("engine")
                if engine is None:
                    return
                # Best-effort cooperative stop — same path the user's
                # Stop button uses. Doesn't kill the subprocess, just
                # asks the worker to finish the current turn.
                try:
                    engine.request_stop()
                    if hasattr(engine.client, "signal_stop"):
                        engine.client.signal_stop()
                except Exception:
                    pass
                _append_system_message(
                    f"⏱ Cooperative stop sent — stream silent for "
                    f"{int(elapsed)} s (> kill threshold {int(kill_after)} s). "
                    "No tokens were being produced; the turn is being "
                    "ended to save wall-clock + API cost."
                )
            except Exception:
                pass

        timer = _threading.Timer(threshold + 1.0, _check_stale)
        timer.daemon = True
        timer.start()
        state["_stale_timer"] = timer

        if kill_after > 0:
            kill_timer = _threading.Timer(kill_after + 1.0, _check_kill)
            kill_timer.daemon = True
            kill_timer.start()
            state["_stale_kill_timer"] = kill_timer
        else:
            state["_stale_kill_timer"] = None

    def _update_queue_display():
        """Update the queue indicator."""
        n = len(state["message_queue"])
        if n > 0:
            queue_html.value = (
                f'<span class="delfin-agent-queue">'
                f'{n} message{"s" if n != 1 else ""} queued</span>'
            )
        else:
            queue_html.value = ""

    def _refresh_context_bar() -> None:
        """Repaint the context-window usage bar from the live engine state.

        Reads ``engine._estimate_context_tokens()`` + ``engine.context_window_tokens``
        and draws a 12-px horizontal bar with proportional fill. Colour ramps
        from blue → orange → red as the percentage crosses 60 %/80 %.
        Hidden when the engine is missing or holds zero messages.
        """
        eng = state.get("engine")
        if not eng or not getattr(eng, "messages", None):
            context_bar_html.value = ""
            return
        try:
            tokens = int(eng._estimate_context_tokens())
        except Exception:
            context_bar_html.value = ""
            return
        window = int(getattr(eng, "context_window_tokens", 100_000) or 100_000)
        if window <= 0:
            context_bar_html.value = ""
            return
        pct = min(100.0, tokens / window * 100.0)
        if pct >= 80.0:
            fill = "#ef4444"
        elif pct >= 60.0:
            fill = "#f59e0b"
        else:
            fill = "#60a5fa"
        n_msgs = len(eng.messages)
        compact_pct = float(getattr(eng, "auto_compact_pct", 0.95) or 0.95) * 100.0
        context_bar_html.value = (
            f'<div style="display:flex; align-items:center; gap:8px; '
            f'margin:2px 0 6px 0; font-family:monospace; font-size:11px; '
            f'color:#94a3b8;">'
            f'<span>ctx</span>'
            f'<div style="position:relative; flex:1; height:8px; '
            f'background:#1e293b; border-radius:4px; overflow:hidden; '
            f'box-shadow:inset 0 1px 2px rgba(0,0,0,0.4);">'
            f'<div style="position:absolute; left:0; top:0; bottom:0; '
            f'width:{pct:.2f}%; background:{fill}; '
            f'transition:width 0.4s ease, background 0.4s ease;"></div>'
            f'<div style="position:absolute; left:{compact_pct:.2f}%; top:-2px; '
            f'bottom:-2px; width:1px; background:#475569;" '
            f'title="auto-compact trigger"></div>'
            f'</div>'
            f'<span style="color:{fill};">{tokens:,}</span>'
            f'<span>/ {window:,} ({pct:.1f}%)</span>'
            f'<span style="color:#64748b;">· {n_msgs} msgs</span>'
            f'</div>'
        )

    def _update_button_states():
        engine = state["engine"]
        is_streaming = state["streaming"]
        # Keep textarea + send enabled during streaming for message queuing
        send_btn.disabled = False
        input_textarea.disabled = False
        stop_btn.disabled = not is_streaming
        mode_dropdown.disabled = is_streaming
        model_dropdown.disabled = is_streaming
        effort_dropdown.disabled = is_streaming
        perm_dropdown.disabled = is_streaming

        if engine and not is_streaming:
            can_advance = (
                engine.current_role_index < len(engine.route) - 1
                and len(engine.messages) > 0
            )
            advance_btn.disabled = not can_advance
        else:
            advance_btn.disabled = True

    # -- event handlers ----------------------------------------------------

    def _handle_slash_command(text: str) -> bool:
        """Handle slash commands. Returns True if handled."""
        import subprocess as _sp

        cmd = text.lower().strip()

        # -- Global safety: zone-based permission check (all modes) ----------
        tier = _command_tier(cmd)
        block_msg = _zone_blocks(text, tier)
        if block_msg:
            _append_system_message(block_msg)
            return True

        # /skill <name>  — expand a skill markdown into the input field
        if cmd.startswith("/skill "):
            try:
                from delfin.agent.skills import (
                    get_skill, render_skill_invocation, discover_skills,
                )
            except Exception as exc:  # pragma: no cover
                _append_system_message(f"Skill system unavailable: {exc}")
                return True
            name = cmd[len("/skill "):].strip()
            if not name:
                names = [s.name for s in discover_skills()]
                _append_system_message(
                    "Available skills:\n  " + "\n  ".join(names) if names
                    else "No skills installed in ~/.delfin/skills/."
                )
                return True
            skill = get_skill(name)
            if skill is None:
                _append_system_message(f"Unknown skill '{name}'.")
                return True
            # Inline-expand: load the skill body into the input so the user
            # can review / edit before sending.
            input_textarea.value = render_skill_invocation(skill)
            _append_system_message(
                f"Skill '{skill.name}' loaded into the input. "
                "Review and press Send to run it."
            )
            return True

        if cmd == "/help":
            _append_system_message(
                "Available commands:\n"
                "  /help            — Show this help\n"
                "  /clear           — Clear chat history\n"
                "  /cost            — Show token usage & cost\n"
                "  /compact         — Summarize context (reduce tokens)\n"
                "  /context         — Show context-window usage + compaction state\n"
                "  /agents          — List subagent presets (explore/plan/...)\n"
                "  /skills          — List discovered skills\n"
                "  /undo            — Undo last agent turn (remove from context)\n"
                "  /retry           — Retry last message (undo + re-send)\n"
                "  /stop            — Stop current generation\n"
                "  /status          — Show engine status\n"
                "  /usage           — Detailed token usage, cost & session stats\n"
                "  /export          — Export chat as Markdown file\n"
                "  /search <text>   — Search in chat history\n"
                "  /retry           — Regenerate last response\n"
                "  /git status      — Show git status\n"
                "  /git diff        — Show staged/unstaged changes\n"
                "  /git log         — Show recent commits\n"
                "  /git branch      — Show branches\n"
                "  /provider <name> — Switch provider (claude/openai)\n"
                "  /model <name>    — Switch model (depends on provider)\n"
                "  /effort <lvl>    — Set effort (low/medium/high/xhigh)\n"
                "  /mode <name>     — Switch mode (dashboard/solo/quick/reviewed/tdd/cluster/full)\n"
                "  /perms [profile] — Show/set permission profile (plan/ask_all/repo_free/all_free)\n"
                "  /reset           — Reset engine for new cycle\n"
                "\n"
                "Dashboard control:\n"
                "  /ui list         — List all controllable widgets\n"
                "  /ui <w> show     — Show widget properties & value\n"
                "  /ui <w> click    — Press a button\n"
                "  /ui <w> value <v> — Set widget value (text/number/dropdown)\n"
                "  /ui <w> options  — Show dropdown choices\n"
                "  /ui <w> style <s> — Button color (primary/danger/success/info/warning)\n"
                "  /ui <w> text/disabled/visible/width/height <v>\n"
                "  /tab <name>      — Switch tab (submit/orca/jobs/calc/settings)\n"
                "  /control show    — Show CONTROL content from Submit tab\n"
                "  /control set ... — Set CONTROL content in Submit tab\n"
                "  /control key k v — Change single CONTROL key (e.g. /control key functional BP86)\n"
                "  /control validate — Validate CONTROL syntax\n"
                "  /submit          — Submit job (confirms first)\n"
                "  /orca show       — Show ORCA Builder settings\n"
                "  /orca set <p> <v> — Set ORCA Builder param (method/basis/charge/...)\n"
                "  /orca submit     — Submit ORCA job\n"
                "  /jobs            — Switch to Job Status tab\n"
                "\n"
                "Calculations & analysis:\n"
                "  /calc ls [path]  — List calc directories/files\n"
                "  /calc cd <path>  — Navigate calc folder (syncs browser)\n"
                "  /calc select <f> — Open/select file in browser preview\n"
                "  /calc open <f> — Alias for /calc select\n"
                "  /calc read <file> — Read a calc file\n"
                "  /calc tail <file> — Read last 8KB of output\n"
                "  /calc info <dir> — Show folder summary & status\n"
                "  /calc tree [dir] — Show directory tree\n"
                "  /calc search <p> — Search files by glob pattern\n"
                "  /analyze <dir>   — Full analysis (energy+convergence+errors)\n"
                "  /analyze energy <dir> — Extract energies (Gibbs/ZPE/electronic)\n"
                "  /analyze rank <type> [limit] [root] — Rank energies across folders\n"
                "  /analyze convergence <dir> — Check SCF convergence\n"
                "  /analyze errors <dir> — Scan for ORCA errors\n"
                "  /analyze status  — Overview of all calc folders\n"
                "\n"
                "Recalc & cancel (require confirmation):\n"
                "  /recalc check <dir> — Check if recalc needed\n"
                "  /recalc check-all — Scan all folders\n"
                "  /recalc <dir>    — Submit recalc (confirms first)\n"
                "  /recalc auto     — Recalc all that need it (confirms first)\n"
                "  /cancel <job_id> — Cancel a job (confirms first)\n"
                "  /cancel all      — Cancel all active jobs (confirms first)\n"
                "  /cancel running  — Cancel only running (R) jobs\n"
                "  /cancel pending  — Cancel only pending (PD) jobs\n"
                "\n"
                "Memory:\n"
                "  /remember <text> — Save a persistent memory\n"
                "  /memories        — List all memories\n"
                "  /memories verify — Check memory file-refs against the repo\n"
                "  /forget <index>  — Delete a memory by index\n"
                "\n"
                "Workspace:\n"
                "  /workspace ls    — List files in agent workspace\n"
                "  /workspace read <file> — Read a workspace file\n"
                "  /workspace clean — Remove all workspace files\n"
                "\n"
                "Keyboard shortcuts:\n"
                "  Enter            — Send message\n"
                "  Shift+Enter      — New line\n"
                "  Escape           — Stop generation\n"
                "  Ctrl+L           — Clear chat\n"
                "  Ctrl+K           — Toggle search\n"
                "  Shift+Tab        — Cycle permission mode"
            )
            return True

        if cmd == "/clear":
            state["chat_messages"].clear()
            _refresh_chat_html()
            _append_system_message("Chat cleared.")
            return True

        if cmd == "/cost":
            engine = state["engine"]
            if engine:
                s = engine.get_status()
                inp_t = s["input_tokens"]
                out_t = s["output_tokens"]
                cost = s["cost_usd"]
                cost_str = f"${cost:.4f}" if cost > 0 else _estimate_cost_str(
                    s["backend"], inp_t, out_t,
                    provider=provider_dropdown.value,
                )
                _append_system_message(
                    f"Token usage:\n"
                    f"  Input:  {inp_t:,} tokens\n"
                    f"  Output: {out_t:,} tokens\n"
                    f"  Cost:   {cost_str}\n"
                    f"  Provider: {provider_dropdown.value}\n"
                    f"  Model:  {model_dropdown.value}\n"
                    f"  Effort: {effort_dropdown.value}"
                )
            else:
                _append_system_message("No active engine. Send a message first.")
            return True

        if cmd == "/usage":
            engine = state["engine"]
            if not engine:
                _append_system_message("No active engine. Send a message first.")
                return True
            s = engine.get_status()
            inp_t = s["input_tokens"]
            out_t = s["output_tokens"]
            total_t = inp_t + out_t
            cost = s["cost_usd"]
            model = model_dropdown.value
            backend = s["backend"]
            cost_str = f"${cost:.4f}" if cost > 0 else _estimate_cost_str(
                backend, inp_t, out_t,
                provider=provider_dropdown.value,
            )
            # Per-model pricing estimate
            _PRICING = {
                "opus":   {"input": 15.0, "output": 75.0},
                "sonnet": {"input": 3.0,  "output": 15.0},
                "haiku":  {"input": 0.25, "output": 1.25},
                "gpt-5.4": {"input": 2.0, "output": 8.0},
                "gpt-5.4-mini": {"input": 0.40, "output": 1.60},
                "gpt-5.3-codex": {"input": 2.0, "output": 8.0},
                "gpt-5.2-codex": {"input": 2.0, "output": 8.0},
                "gpt-5.2": {"input": 2.0, "output": 8.0},
                "gpt-5.1-codex-max": {"input": 2.0, "output": 8.0},
                "gpt-5.1-codex-mini": {"input": 0.40, "output": 1.60},
                "gpt-4.1": {"input": 2.0, "output": 8.0},
                "gpt-4.1-mini": {"input": 0.40, "output": 1.60},
                "o4-mini": {"input": 1.10, "output": 4.40},
                "o3": {"input": 2.0, "output": 8.0},
            }
            # Ollama / vLLM / LM Studio are local — no per-token cost.
            if (provider_dropdown.value == "ollama"
                    or model.startswith(("qwen", "llama", "mistral",
                                          "deepseek", "phi", "gemma"))):
                pricing = {"input": 0.0, "output": 0.0}
            else:
                pricing = _PRICING.get(model, _PRICING.get("sonnet", {"input": 3.0, "output": 15.0}))
            est_in = inp_t * pricing["input"] / 1_000_000
            est_out = out_t * pricing["output"] / 1_000_000
            # Message counts
            user_msgs = sum(1 for m in state["chat_messages"] if m["role"] == "user")
            asst_msgs = sum(1 for m in state["chat_messages"] if m["role"] == "assistant")
            total_msgs = user_msgs + asst_msgs
            avg_tok = total_t // max(total_msgs, 1)
            # Session duration
            start = state.get("session_start_time")
            if start:
                elapsed = time.monotonic() - start
                mins, secs = int(elapsed // 60), int(elapsed % 60)
                dur = f"{mins}m {secs}s"
                rate = f"{total_t / max(elapsed, 1):.0f} tok/s"
            else:
                dur = rate = "N/A"
            _append_system_message(
                f"Session Usage:\n"
                f"  Provider:    {provider_dropdown.value}\n"
                f"  Model:       {model} ({backend})\n"
                f"  Permission:  {state.get('_perm_profile', 'default')} (CLI: {_active_cli_perm()})\n"
                f"  Effort:      {effort_dropdown.value}\n"
                f"  Session:     {s['session_id'][:16]}...\n"
                f"\n"
                f"Tokens:\n"
                f"  Input:       {inp_t:,}\n"
                f"  Output:      {out_t:,}\n"
                f"  Total:       {total_t:,}\n"
                f"\n"
                f"Cost:\n"
                f"  Input:       ${est_in:.4f} (${pricing['input']}/MTok)\n"
                f"  Output:      ${est_out:.4f} (${pricing['output']}/MTok)\n"
                f"  Total:       {cost_str}\n"
                f"\n"
                f"Stats:\n"
                f"  Duration:    {dur}\n"
                f"  Messages:    {user_msgs} sent / {asst_msgs} received\n"
                f"  Avg tok/msg: {avg_tok:,}\n"
                f"  Token rate:  {rate}\n"
                f"  Context:     {len(engine.messages)} engine msgs"
            )
            return True

        if cmd == "/status":
            engine = state["engine"]
            if engine:
                s = engine.get_status()
                _append_system_message(
                    f"Engine status:\n"
                    f"  Mode:    {s['mode']}\n"
                    f"  Provider: {s.get('provider', 'claude')}\n"
                    f"  Backend: {s['backend']}\n"
                    f"  Role:    {_format_role_label(s['role'])} "
                    f"({s['role_index']+1}/{s['role_total']})\n"
                    f"  Session: {s['session_id'][:16]}...\n"
                    f"  Cycle:   {'complete' if s['cycle_complete'] else 'in progress'}"
                )
            else:
                _append_system_message("No active engine.")
            return True

        if cmd == "/stop":
            _on_stop(None)
            return True

        if cmd == "/reset":
            _on_new_cycle(None)
            return True

        if cmd == "/compact":
            engine = state["engine"]
            if not engine or not engine.messages:
                _append_system_message("Nothing to compact.")
                return True
            # Keep only the last 4 messages (2 turns) to reduce context
            n_before = len(engine.messages)
            if n_before > 4:
                # Summarize old messages into a single context note
                old_msgs = engine.messages[:-4]
                summary_parts = []
                for m in old_msgs:
                    role = m["role"]
                    content = m["content"][:200]
                    summary_parts.append(f"[{role}]: {content}...")
                summary = "\n".join(summary_parts[-6:])  # last 6 entries
                engine.messages = [
                    {"role": "user", "content":
                     f"[Context summary of {n_before - 4} earlier messages:\n"
                     f"{summary}\n... End of summary]"},
                    {"role": "assistant", "content": "Understood, I have the context."},
                ] + engine.messages[-4:]
                _append_system_message(
                    f"Compacted: {n_before} messages → {len(engine.messages)} "
                    f"(older context summarized)"
                )
            else:
                _append_system_message(
                    f"Only {n_before} messages — too few to compact."
                )
            return True

        if cmd == "/failures" or cmd.startswith("/failures "):
            from delfin.agent import failure_log as _fl
            arg = cmd[len("/failures "):].strip() if cmd.startswith("/failures ") else ""
            if arg in ("", "top"):
                top = _fl.top_recurring(min_count=2)
                if not top:
                    _append_system_message(
                        "No recurring failures recorded yet. "
                        "Failures are auto-logged when tools return errors; "
                        "patterns surface here once they hit 2+ occurrences."
                    )
                    return True
                lines = [
                    "Recurring tool failures (≥2 in the last 7 days):"
                ]
                import time as _t
                for r in top[:10]:
                    when = _t.strftime("%Y-%m-%d %H:%M",
                                       _t.localtime(r["last_seen"] or 0))
                    lines.append(
                        f"  {r['count']:>3}×  {r['tool']:<22}  "
                        f"last {when}\n      err: {r['last_error'][:90]}"
                    )
                lines.append(
                    "\n/failures clear — empty the log\n"
                    "/failures recent — last 20 failures verbatim"
                )
                _append_system_message("\n".join(lines))
                return True
            if arg == "recent":
                recs = _fl.read_failures(last_n=20)
                if not recs:
                    _append_system_message("No failures logged.")
                    return True
                import time as _t
                lines = ["Last 20 failures:"]
                for r in recs:
                    when = _t.strftime("%H:%M:%S",
                                       _t.localtime(r.get("ts") or 0))
                    lines.append(
                        f"  {when}  {r.get('tool',''):<22}  "
                        f"{r.get('error','')[:90]}"
                    )
                _append_system_message("\n".join(lines))
                return True
            if arg == "clear":
                from pathlib import Path as _P
                p = _P.home() / ".delfin" / "failure_log.jsonl"
                try:
                    if p.exists():
                        p.unlink()
                    _append_system_message("Failure log cleared.")
                except Exception as exc:
                    _append_system_message(f"Could not clear log: {exc}")
                return True
            _append_system_message(
                "Usage:\n"
                "  /failures              — top recurring failure shapes\n"
                "  /failures recent       — last 20 failures verbatim\n"
                "  /failures clear        — empty the log"
            )
            return True

        if cmd == "/bash" or cmd.startswith("/bash "):
            from delfin.agent import bash_jobs as _bj
            arg = cmd[len("/bash "):].strip() if cmd.startswith("/bash ") else ""
            reg = _bj.get_registry()
            watchers = state.setdefault("_bg_watchers", {})
            # /bash jobs  — list active + recently-finished jobs
            if arg in ("", "jobs", "ls"):
                jobs = list(reg._jobs.values())
                if not jobs:
                    _append_system_message("No background bash jobs.")
                    return True
                jobs.sort(key=lambda j: j.started_at, reverse=True)
                lines = ["Background bash jobs:"]
                for j in jobs[:20]:
                    st = j.status_dict()
                    flag = "▶" if st["running"] else "■"
                    watched = " (watching)" if j.job_id in watchers else ""
                    lines.append(
                        f"  {flag} {j.job_id}  exit={st['exit_code']}  "
                        f"{st['elapsed_s']:>6.1f}s  {st['command'][:60]}"
                        f"{watched}"
                    )
                lines.append(
                    "\n/bash watch <id> — stream new output to chat\n"
                    "/bash unwatch <id> — stop streaming\n"
                    "/bash kill <id> — SIGTERM then SIGKILL"
                )
                _append_system_message("\n".join(lines))
                return True
            # /bash watch <id>  — stream new lines as system messages
            if arg.startswith("watch "):
                jid = arg[len("watch "):].strip()
                job = reg._jobs.get(jid)
                if job is None:
                    _append_system_message(f"No bash job '{jid}'.")
                    return True
                if jid in watchers:
                    _append_system_message(
                        f"Already watching {jid}. /bash unwatch {jid} to stop."
                    )
                    return True
                import threading as _th
                stop_ev = _th.Event()
                last_size = {"stdout": 0, "stderr": 0}

                def _stream():
                    """Poll the job's stdout/stderr files every second and
                    emit any new bytes as a system message. Auto-stops
                    when the job exits or stop_ev is set."""
                    try:
                        while not stop_ev.is_set():
                            try:
                                for stream_name, path in (
                                    ("stdout", job.stdout_path),
                                    ("stderr", job.stderr_path),
                                ):
                                    if not path.is_file():
                                        continue
                                    sz = path.stat().st_size
                                    if sz > last_size[stream_name]:
                                        with path.open("r", encoding="utf-8", errors="replace") as f:
                                            f.seek(last_size[stream_name])
                                            chunk = f.read()
                                        last_size[stream_name] = sz
                                        if chunk.strip():
                                            tag = "stderr" if stream_name == "stderr" else "stdout"
                                            _append_system_message(
                                                f"[{jid} · {tag}]\n{chunk.rstrip()}"
                                            )
                            except Exception:
                                pass
                            # Stop if the job finished
                            if job.poll() is not None:
                                stop_ev.set()
                                break
                            stop_ev.wait(timeout=1.0)
                    finally:
                        watchers.pop(jid, None)
                        _append_system_message(
                            f"[{jid}] watcher stopped "
                            f"(exit={job.status_dict()['exit_code']})."
                        )

                t = _th.Thread(target=_stream, daemon=True)
                t.start()
                watchers[jid] = (t, stop_ev)
                _append_system_message(
                    f"▶ Watching {jid} — new output will stream to chat. "
                    f"Stop with /bash unwatch {jid}."
                )
                return True
            # /bash unwatch <id>  — stop a watcher
            if arg.startswith("unwatch "):
                jid = arg[len("unwatch "):].strip()
                entry = watchers.pop(jid, None)
                if entry is None:
                    _append_system_message(f"Not watching '{jid}'.")
                else:
                    _, stop_ev = entry
                    stop_ev.set()
                    _append_system_message(f"Stopped watching {jid}.")
                return True
            # /bash kill <id>  — SIGTERM then SIGKILL the underlying process
            if arg.startswith("kill "):
                jid = arg[len("kill "):].strip()
                job = reg._jobs.get(jid)
                if job is None:
                    _append_system_message(f"No bash job '{jid}'.")
                    return True
                try:
                    job.proc.terminate()
                    try:
                        job.proc.wait(timeout=3)
                    except Exception:
                        job.proc.kill()
                except Exception as exc:
                    _append_system_message(f"Kill failed: {exc}")
                    return True
                _append_system_message(
                    f"🛑 {jid} killed (exit={job.proc.returncode})."
                )
                return True
            _append_system_message(
                "Usage:\n"
                "  /bash jobs               — list active background jobs\n"
                "  /bash watch <id>         — stream new output to chat\n"
                "  /bash unwatch <id>       — stop streaming\n"
                "  /bash kill <id>          — kill a job"
            )
            return True

        if cmd == "/init" or cmd.startswith("/init "):
            arg = cmd[len("/init "):].strip() if cmd.startswith("/init ") else ""
            force = arg in ("force", "--force", "-f", "overwrite")
            try:
                from delfin.agent.project_init import init_project
                result = init_project(ctx.repo_dir or ".", overwrite=force)
            except Exception as exc:
                _append_system_message(f"/init failed: {exc}")
                return True
            p = result["profile"]
            lines = [
                f"🚀 /init complete — detected **{p.language}** project "
                f"'{p.name}'."
            ]
            if p.test_command:
                lines.append(f"  test: `{p.test_command}`")
            if p.build_command:
                lines.append(f"  build: `{p.build_command}`")
            if p.lint_command:
                lines.append(f"  lint: `{p.lint_command}`")
            lines.append("")
            if result["files"]:
                lines.append("Created:")
                for f in result["files"]:
                    short = f.replace(str(Path.home()), "~")
                    lines.append(f"  + {short}")
            if result["skipped"]:
                lines.append("Skipped (already exists, use `/init force`):")
                for f in result["skipped"]:
                    short = f.replace(str(Path.home()), "~")
                    lines.append(f"  · {short}")
            lines.append(
                "\nThe new AGENTS.md is auto-loaded into the prompt on the "
                "next turn — edit it to refine the project description."
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/commands" or cmd.startswith("/commands "):
            from delfin.agent import slash_commands as _scm
            arg = cmd[len("/commands "):].strip() if cmd.startswith("/commands ") else ""
            cmds = _scm.discover_commands(ctx.repo_dir or None)
            if arg:
                # Detail view
                tpl = _scm.get_command(arg, ctx.repo_dir or None)
                if tpl is None:
                    _append_system_message(
                        f"No command '{arg}'. Try /commands for the full list."
                    )
                    return True
                src = str(tpl.source).replace(str(Path.home()), "~")
                hint = f"  arg-hint: {tpl.argument_hint}\n" if tpl.argument_hint else ""
                _append_system_message(
                    f"## /{tpl.name}\n"
                    f"_{tpl.description}_\n{hint}\n"
                    f"---\n\n{tpl.body}\n\n---\n"
                    f"Source: `{src}`\n"
                    f"Invoke: `/{tpl.name} <args>` — $ARGUMENTS / $1 / $2 substituted."
                )
                return True
            if not cmds:
                _append_system_message(
                    "No custom commands installed.\n"
                    "Drop a markdown file in ~/.delfin/commands/<name>.md or "
                    "in <workspace>/.delfin/commands/<name>.md.\n"
                    "The file body becomes the prompt; $ARGUMENTS expands to "
                    "the text the user types after the command."
                )
                return True
            lines = ["User-defined slash commands:"]
            for tpl in cmds:
                hint = f" {tpl.argument_hint}" if tpl.argument_hint else ""
                lines.append(
                    f"  /{tpl.name}{hint:<18}  {tpl.description[:60]}"
                )
            lines.append(
                "\nDetail view: /commands <name>"
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/mcp" or cmd.startswith("/mcp "):
            from delfin.agent import mcp_editor as _me
            arg = cmd[len("/mcp "):].strip() if cmd.startswith("/mcp ") else ""

            def _hot_reload_mcp(reason: str = "") -> None:
                """Tear down + re-load the process-wide MCP registry so
                new server entries take effect without a dashboard restart."""
                try:
                    from delfin.agent import mcp_client as _mcp
                    _mcp.reset_registry()
                    new_reg = _mcp.get_registry(ctx.repo_dir or None)
                    n_servers = len(getattr(new_reg, "servers", []) or [])
                    n_tools = len(getattr(new_reg, "tools", {}) or {})
                    _append_system_message(
                        f"♻️ MCP registry reloaded{(' — ' + reason) if reason else ''}: "
                        f"{n_servers} server(s), {n_tools} tool(s) live."
                    )
                except Exception as exc:
                    _append_system_message(f"MCP reload failed: {exc}")

            # /mcp reload  — force a process-wide reload now
            if arg in ("reload", "refresh"):
                _hot_reload_mcp("manual reload")
                return True
            # /mcp add <name> <command> [arg ...]
            if arg.startswith("add "):
                rest = arg[4:].strip()
                parts = rest.split(None, 1)
                if len(parts) < 2:
                    _append_system_message(
                        "Usage: /mcp add <name> <command> [arg ...]"
                    )
                    return True
                name = parts[0]
                # Tokenize the rest via shlex so quoted args survive
                import shlex
                try:
                    tokens = shlex.split(parts[1])
                except ValueError as exc:
                    _append_system_message(f"Bad quoting: {exc}")
                    return True
                command = tokens[0] if tokens else ""
                args_list = tokens[1:]
                try:
                    rec = _me.add_mcp_server(name, command, args_list)
                except Exception as exc:
                    _append_system_message(f"Add failed: {exc}")
                    return True
                _append_system_message(
                    f"✅ MCP server '{rec['name']}' added: "
                    f"{rec['command']} {' '.join(rec['args'])}"
                )
                _hot_reload_mcp(f"added '{rec['name']}'")
                return True
            # /mcp remove <name>
            if arg.startswith("remove ") or arg.startswith("rm "):
                name = arg.split(None, 1)[1].strip() if " " in arg else ""
                if not name:
                    _append_system_message("Usage: /mcp remove <name>")
                    return True
                try:
                    removed = _me.remove_mcp_server(name)
                except Exception as exc:
                    _append_system_message(f"Remove failed: {exc}")
                    return True
                if removed is None:
                    _append_system_message(f"No MCP server named '{name}'.")
                else:
                    _append_system_message(f"🗑 MCP server '{name}' removed.")
                    _hot_reload_mcp(f"removed '{name}'")
                return True
            # /mcp enable <name>  /  /mcp disable <name>
            for verb, flag in (("enable ", True), ("disable ", False)):
                if arg.startswith(verb):
                    name = arg[len(verb):].strip()
                    try:
                        rec = _me.toggle_mcp_server(name, enabled=flag)
                    except Exception as exc:
                        _append_system_message(f"Toggle failed: {exc}")
                        return True
                    if rec is None:
                        _append_system_message(f"No MCP server named '{name}'.")
                    else:
                        _append_system_message(
                            f"{'✅' if flag else '⏸'} MCP server "
                            f"'{name}' → enabled={flag}."
                        )
                        _hot_reload_mcp(f"{'enabled' if flag else 'disabled'} '{name}'")
                    return True
            # /mcp  → list
            try:
                rows = _me.list_mcp_servers()
            except Exception as exc:
                _append_system_message(f"Could not list MCP servers: {exc}")
                return True
            if not rows:
                _append_system_message(
                    "No user-global MCP servers configured.\n"
                    "Add one with: /mcp add <name> <command> [args ...]\n"
                    "Example: /mcp add fs npx -y @modelcontextprotocol/server-filesystem /tmp"
                )
                return True
            lines = ["User-global MCP servers (~/.delfin/mcp_servers.json):"]
            for r in rows:
                badge = "✓" if r["enabled"] else "✗"
                argv = " ".join([r["command"]] + list(r["args"]))[:80]
                lines.append(f"  [{badge}] {r['name']:<14}  {argv}")
            lines.append(
                "\nCommands: /mcp add | remove | enable | disable | reload\n"
                "Note: add/remove/toggle auto-reload the registry; "
                "use /mcp reload to pick up hand-edits to mcp_servers.json."
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/session" or cmd.startswith("/session "):
            from delfin.agent import session_store as _ss
            arg = cmd[len("/session "):].strip() if cmd.startswith("/session ") else ""
            # /session fork [<id>]  — branch a session (default: current)
            if arg == "fork" or arg.startswith("fork "):
                target = arg[len("fork "):].strip() if arg.startswith("fork ") else ""
                if not target:
                    target = state.get("active_session_id", "") or ""
                if not target:
                    _append_system_message(
                        "No active session yet. /session fork <id> or send "
                        "a message first to mint a session ID."
                    )
                    return True
                try:
                    new_sid = _ss.fork_session(target)
                except Exception as exc:
                    _append_system_message(f"Fork failed: {exc}")
                    return True
                if not new_sid:
                    _append_system_message(f"Could not fork session '{target}'.")
                    return True
                _refresh_session_dropdown()
                # Load the fork into the current tab so the user is now
                # editing the branch — the parent stays untouched on disk.
                try:
                    _load_saved_session(new_sid)
                except Exception:
                    pass
                _append_system_message(
                    f"🌿 Forked `{target[:12]}…` → `{new_sid}`. "
                    "You are now editing the fork; parent unchanged."
                )
                return True
            # /session tree [<id>]  — ancestry + immediate children
            if arg == "tree" or arg.startswith("tree "):
                target = arg[len("tree "):].strip() if arg.startswith("tree ") else ""
                if not target:
                    target = state.get("active_session_id", "") or ""
                if not target:
                    _append_system_message("Usage: /session tree <session_id>")
                    return True
                ancestry = _ss.session_lineage(target)
                children = _ss.session_children(target)
                if not ancestry:
                    _append_system_message(f"No session found for '{target}'.")
                    return True
                import time as _t
                lines = [f"Session lineage for {target[:12]}…:"]
                # Ancestry chain, oldest first for readability
                for i, rec in enumerate(reversed(ancestry)):
                    indent = "  " * i
                    when = _t.strftime("%Y-%m-%d %H:%M",
                                       _t.localtime(rec.get("updated_at", 0)))
                    arrow = "→" if i > 0 else "•"
                    title = (rec.get("title") or "")[:50]
                    n = rec.get("message_count", 0)
                    sid_short = rec["session_id"][:12]
                    lines.append(
                        f"{indent}{arrow} {sid_short:<14}  {when}  "
                        f"{n:>3} msgs  {title}"
                    )
                # Immediate children (siblings of leaf)
                if children:
                    lines.append("\nForks from this session:")
                    for child in children:
                        when = _t.strftime("%Y-%m-%d %H:%M",
                                           _t.localtime(child.get("updated_at", 0)))
                        title = (child.get("title") or "")[:50]
                        n = child.get("message_count", 0)
                        lines.append(
                            f"  ⤷ {child['session_id'][:14]:<14}  {when}  "
                            f"{n:>3} msgs  {title}"
                        )
                _append_system_message("\n".join(lines))
                return True
            # /session restore <id>  — load a saved session into this tab
            if arg.startswith("restore "):
                sid = arg[len("restore "):].strip()
                if not sid:
                    _append_system_message("Usage: /session restore <session_id>")
                    return True
                try:
                    _load_saved_session(sid)
                except Exception as exc:
                    _append_system_message(f"Restore failed: {exc}")
                return True
            # /session handoff [new] [<id>]  — produce a fresh-agent brief
            if arg == "handoff" or arg.startswith("handoff "):
                rest = arg[len("handoff"):].strip()
                start_fresh = rest.startswith("new")
                if start_fresh:
                    rest = rest[len("new"):].strip()
                target = rest or state.get("active_session_id", "") or ""
                # Build the brief from the saved session, or — if we're
                # handing off the LIVE session and it hasn't been saved
                # yet — from the in-memory state.
                data = None
                if target:
                    try:
                        data = _ss.load_session(target)
                    except Exception:
                        data = None
                if data is None:
                    # Live session: assemble a minimal data dict from state
                    eng = state.get("engine")
                    estate = eng.export_state() if eng else {}
                    data = {
                        "session_id": target or state.get("active_session_id", "live"),
                        "title": "",
                        "chat_messages": state.get("chat_messages", []),
                        "engine_messages": estate.get("engine_messages", []),
                        "token_usage": estate.get("token_usage", {}),
                        "cost_usd": estate.get("cost_usd", 0.0),
                        "todo_payload": state.get("current_todos", []),
                        "active_gate": state.get("_active_gate"),
                        "pending_plan_body": state.get("_pending_plan_body", ""),
                    }
                try:
                    brief = _ss.build_handoff_brief(data)
                    saved = _ss.save_handoff_brief(
                        data.get("session_id", "session"), brief,
                    )
                except Exception as exc:
                    _append_system_message(f"Handoff generation failed: {exc}")
                    return True
                short = str(saved).replace(str(Path.home()), "~")
                if start_fresh:
                    # Mint a fresh session: new engine, new session_id,
                    # zero history — seeded with the handoff brief as the
                    # opening context. This is the "pass to a new agent".
                    try:
                        _on_new_cycle(None)  # reset engine + chat
                    except Exception:
                        pass
                    _append_chat_message("system", (
                        "🤝 **New agent session** — seeded with a handoff "
                        f"brief from `{data.get('session_id','?')[:16]}`. "
                        "The agent below has no prior conversation history; "
                        "everything it needs is in the brief."
                    ))
                    _append_chat_message("user", brief)
                    # Pre-fill the input so the user can just hit Send to
                    # kick the fresh agent off with the brief as context.
                    input_textarea.value = (
                        "Continue from the handoff brief above. "
                        "Start with the recommended next step."
                    )
                    _append_system_message(
                        f"📝 Handoff brief saved → {short}\n"
                        "Fresh session ready — press Send to start the new "
                        "agent, or edit the prompt first."
                    )
                else:
                    _append_system_message(
                        f"📝 Handoff brief for `{data.get('session_id','?')[:16]}`"
                        f" saved → {short}\n\n{brief}"
                    )
                return True
            # /session bundle <id>  — export a portable .delfin-session zip
            if arg.startswith("bundle"):
                sid = arg[len("bundle"):].strip() or state.get("active_session_id", "")
                if not sid:
                    _append_system_message(
                        "Usage: /session bundle <session_id>  (current "
                        "session has no id yet — send a message first)"
                    )
                    return True
                try:
                    bundle = _ss.export_bundle(sid)
                except Exception as exc:
                    _append_system_message(f"Bundle export failed: {exc}")
                    return True
                if bundle is None:
                    _append_system_message(f"No session '{sid}' to bundle.")
                    return True
                short = str(bundle).replace(str(Path.home()), "~")
                _append_system_message(
                    f"📦 Session bundled → {short}\n"
                    "Contains session.json + transcript archive + handoff "
                    "brief. Hand this file to another machine / agent and "
                    "import it with `/session import <path>`."
                )
                return True
            # /session import <path>  — ingest a .delfin-session bundle
            if arg.startswith("import "):
                path = arg[len("import "):].strip().strip('"').strip("'")
                if not path:
                    _append_system_message("Usage: /session import <bundle-path>")
                    return True
                try:
                    new_sid = _ss.import_bundle(Path(path).expanduser())
                except Exception as exc:
                    _append_system_message(f"Import failed: {exc}")
                    return True
                if new_sid is None:
                    _append_system_message(
                        f"Could not import '{path}' — not a valid "
                        ".delfin-session bundle."
                    )
                    return True
                _refresh_session_dropdown()
                _append_system_message(
                    f"📥 Bundle imported as session `{new_sid}`. "
                    f"Load it with `/session restore {new_sid}`."
                )
                return True
            # /session search <query>  — grep across saved sessions + archives
            if arg.startswith("search "):
                query = arg[len("search "):].strip()
                if not query:
                    _append_system_message("Usage: /session search <query>")
                    return True
                q = query.lower()
                hits: list[str] = []
                # Chat-history of saved sessions (newest first, capped)
                try:
                    for r in _ss.list_sessions(limit=50):
                        sid = r.get("session_id", "")
                        try:
                            data = _ss.load_session(sid) or {}
                        except Exception:
                            continue
                        msgs = data.get("chat_messages") or []
                        for i, m in enumerate(msgs):
                            content = str(m.get("content", ""))
                            if q in content.lower():
                                snippet = content.replace("\n", " ")[:120]
                                hits.append(
                                    f"  session[{sid[:12]}] msg#{i} "
                                    f"({m.get('role','?')}): {snippet}"
                                )
                                if len(hits) >= 30:
                                    break
                        if len(hits) >= 30:
                            break
                except Exception:
                    pass
                # Pre-compaction archives
                try:
                    for archive in _ss.list_transcript_archives()[:50]:
                        sid = archive["session_id"]
                        for rec in _ss.load_transcript_archive(sid):
                            for i, m in enumerate(rec.get("messages") or []):
                                content = str(m.get("content", ""))
                                if q in content.lower():
                                    snippet = content.replace("\n", " ")[:120]
                                    hits.append(
                                        f"  archive[{sid[:12]}] msg#{i} "
                                        f"({m.get('role','?')}): {snippet}"
                                    )
                                    if len(hits) >= 60:
                                        break
                            if len(hits) >= 60:
                                break
                        if len(hits) >= 60:
                            break
                except Exception:
                    pass
                if not hits:
                    _append_system_message(
                        f"No matches for {query!r} in saved sessions or archives."
                    )
                else:
                    _append_system_message(
                        f"Found {len(hits)} match{'es' if len(hits)!=1 else ''} for "
                        f"{query!r}:\n" + "\n".join(hits[:30])
                        + ("\n  ..." if len(hits) > 30 else "")
                    )
                return True
            # /session archive ls
            if arg in ("archive", "archive ls"):
                rows = _ss.list_transcript_archives()
                if not rows:
                    _append_system_message(
                        "No archived transcripts yet. Archives are written "
                        "automatically when the engine compacts the message "
                        "history — long sessions never truly lose context."
                    )
                    return True
                import time as _t
                lines = ["Archived transcripts (pre-compaction snapshots):"]
                for r in rows[:20]:
                    when = _t.strftime("%Y-%m-%d %H:%M",
                                       _t.localtime(r["mtime"]))
                    kb = r["bytes"] / 1024
                    lines.append(
                        f"  {when}  {r['session_id'][:16]:<18}  "
                        f"{r['compactions']:>2} compactions, {kb:>6.1f} kB"
                    )
                lines.append(
                    "\nLoad a transcript: /session archive <session_id>"
                )
                _append_system_message("\n".join(lines))
                return True
            # /session archive <session_id>
            if arg.startswith("archive "):
                sid = arg[len("archive "):].strip()
                if not sid:
                    _append_system_message("Usage: /session archive <session_id>")
                    return True
                records = _ss.load_transcript_archive(sid)
                if not records:
                    _append_system_message(f"No archive for session '{sid}'.")
                    return True
                import time as _t
                lines = [f"Transcript archive for {sid} ({len(records)} compactions):"]
                for i, rec in enumerate(records):
                    when = _t.strftime("%Y-%m-%d %H:%M",
                                       _t.localtime(rec.get("compacted_at", 0)))
                    info = rec.get("info") or {}
                    n_msgs = rec.get("n_messages", 0)
                    saved = info.get("tokens_before", 0)
                    lines.append(
                        f"  [{i}] {when}  {n_msgs:>3} msgs  "
                        f"~{saved:>6} tokens"
                    )
                lines.append(
                    "\nThe full snapshots live at ~/.delfin/transcript_archive/."
                )
                _append_system_message("\n".join(lines))
                return True
            # /session ls
            if arg == "ls" or arg == "list":
                rows = _ss.list_sessions(limit=20)
                if not rows:
                    _append_system_message("No saved sessions.")
                    return True
                import time as _t
                lines = ["Recent sessions (newest first):"]
                for r in rows:
                    when = _t.strftime("%Y-%m-%d %H:%M",
                                       _t.localtime(r.get("updated_at", 0)))
                    lines.append(
                        f"  {when}  {r.get('session_id', '')[:16]:<18}  "
                        f"{(r.get('title') or '')[:50]}"
                    )
                _append_system_message("\n".join(lines))
                return True
            _append_system_message(
                "Usage:\n"
                "  /session ls                  — recent sessions\n"
                "  /session restore <id>        — load a saved session into this tab\n"
                "  /session search <query>      — grep across saved sessions + archives\n"
                "  /session fork [<id>]         — branch a session (default: current)\n"
                "  /session tree [<id>]         — show ancestry + immediate children\n"
                "  /session handoff [<id>]      — generate a fresh-agent handoff brief\n"
                "  /session handoff new [<id>]  — brief + start a fresh seeded session\n"
                "  /session bundle [<id>]       — export a portable .delfin-session zip\n"
                "  /session import <path>       — ingest a .delfin-session bundle\n"
                "  /session archive             — archived pre-compaction transcripts\n"
                "  /session archive <id>        — view archive for a session\n"
            )
            return True

        if cmd == "/hooks" or cmd.startswith("/hooks "):
            from delfin.agent import hooks_editor as _he
            arg = cmd[len("/hooks "):].strip() if cmd.startswith("/hooks ") else ""
            repo = ctx.repo_dir or None
            # /hooks dry-run <event> [tool]
            if arg.startswith("dry-run") or arg.startswith("dry "):
                parts = arg.split(None, 2)
                event = parts[1] if len(parts) > 1 else ""
                tool = parts[2] if len(parts) > 2 else ""
                if event not in _he._VALID_EVENTS:
                    _append_system_message(
                        f"Usage: /hooks dry-run <event> [tool]\n"
                        f"  events: {', '.join(_he._VALID_EVENTS)}"
                    )
                    return True
                try:
                    rows = _he.dry_run_hook(
                        event,
                        tool_name=tool,
                        tool_input={"file_path": "test.py"} if tool else {},
                        workspace=repo,
                    )
                except Exception as exc:
                    _append_system_message(f"Dry-run failed: {exc}")
                    return True
                if not rows:
                    _append_system_message(
                        f"No hooks fired for {event}"
                        + (f" / {tool}" if tool else "")
                        + ". Add one with `/hooks add`."
                    )
                    return True
                lines = [f"Dry-run {event}" + (f" / {tool}" if tool else "") + ":"]
                for r in rows:
                    flag = "MATCH" if r["matched"] else "skip "
                    lines.append(
                        f"  [{flag}] exit={r['exit_code']:<3} "
                        f"{r['duration_s']:>5.2f}s  {r['command'][:80]}"
                    )
                    if r["decision"]:
                        lines.append(f"    → decision={r['decision']} ({r['reason'][:60]})")
                    if r["stderr_tail"].strip():
                        lines.append(f"    stderr: {r['stderr_tail'][:120]}")
                _append_system_message("\n".join(lines))
                return True
            # /hooks add <event> <matcher> <command...>
            if arg.startswith("add "):
                rest = arg[4:].strip()
                parts = rest.split(None, 2)
                if len(parts) < 3:
                    _append_system_message(
                        "Usage: /hooks add <event> <matcher> <command>\n"
                        "  Use \"\" for an empty matcher."
                    )
                    return True
                event, matcher, command = parts
                matcher = matcher.strip('"').strip("'")
                try:
                    rec = _he.add_hook(event, matcher, command)
                except Exception as exc:
                    _append_system_message(f"Add failed: {exc}")
                    return True
                _append_system_message(
                    f"✅ Hook added → {event}: matcher={rec['matcher']!r}, "
                    f"command={command!r}\nSettings: {rec['settings_path']}"
                )
                return True
            # /hooks remove <event> <index>
            if arg.startswith("remove ") or arg.startswith("rm "):
                rest = arg.split(None, 1)[1] if " " in arg else ""
                parts = rest.split()
                if len(parts) < 2 or not parts[1].isdigit():
                    _append_system_message("Usage: /hooks remove <event> <index>")
                    return True
                event, idx_s = parts[0], parts[1]
                try:
                    removed = _he.remove_hook(event, int(idx_s))
                except Exception as exc:
                    _append_system_message(f"Remove failed: {exc}")
                    return True
                if removed is None:
                    _append_system_message(
                        f"No hook at {event}[{idx_s}] — check /hooks first."
                    )
                else:
                    _append_system_message(
                        f"🗑 Removed {event}[{idx_s}]: {removed}"
                    )
                return True
            # /hooks  → list
            try:
                rows = _he.list_hooks(repo)
            except Exception as exc:
                _append_system_message(f"Could not load hooks: {exc}")
                return True
            if not rows:
                _append_system_message(
                    "No hooks registered.\n\n"
                    "Add one with: /hooks add <event> <matcher> <command>\n"
                    "  events: PreToolUse, PostToolUse, UserPromptSubmit, Stop\n"
                    "  example: /hooks add PreToolUse Edit|Write 'ruff check ${file}'\n"
                )
                return True
            lines = ["Registered hooks (from settings.json layers):"]
            cur_event = ""
            for r in rows:
                if r["event"] != cur_event:
                    lines.append(f"\n{r['event']}:")
                    cur_event = r["event"]
                lines.append(
                    f"  [{r['index']}] matcher={r['matcher']!r:<24}  "
                    f"{r['command'][:80]}"
                )
            lines.append(
                "\nDry-run a hook: /hooks dry-run <event> [tool_name]"
                "\nRemove: /hooks remove <event> <index>"
            )
            _append_system_message("\n".join(lines))
            return True

        # /plan approve | /plan reject — manual fallback for plan-mode
        # when the model forgets to call exit_plan_mode. Especially
        # important for weak models that hallucinate or never emit the
        # tool call. The user can read the plan in chat and unblock the
        # session without restarting.
        if cmd in ("/plan approve", "/plan reject"):
            ev = state.get("_plan_approval_event")
            result = state.get("_plan_approval_result")
            if ev is None or result is None:
                _append_system_message(
                    "No plan is currently awaiting approval. "
                    "Switch to /mode plan and ask the agent to draft a plan."
                )
                return True
            if cmd == "/plan approve":
                result["approved"] = True
                result["new_mode"] = "acceptEdits"
                # Persist the pending plan body the same way the click
                # handler does, so /plan approve and the button are
                # functionally identical.
                plan_body = state.pop("_pending_plan_body", "") or ""
                if plan_body.strip():
                    try:
                        from delfin.agent.memory_store import save_plan
                        fpath = save_plan(
                            plan_body, repo_root=ctx.repo_dir or "."
                        )
                        short = str(fpath).replace(str(Path.home()), "~")
                        result["plan_path"] = str(fpath)
                        _append_system_message(
                            f"📝 Plan saved → {short}"
                        )
                    except Exception as exc:
                        _append_system_message(
                            f"Plan save failed: {exc}"
                        )
                _append_system_message(
                    "✅ Plan approved via /plan approve fallback. "
                    "Mode flipped to acceptEdits — the agent can now "
                    "execute the plan."
                )
            else:
                result["approved"] = False
                result["new_mode"] = "default"
                _append_system_message(
                    "🚫 Plan rejected via /plan reject fallback."
                )
            ev.set()
            return True

        if cmd == "/plans" or cmd.startswith("/plans "):
            from delfin.agent.memory_store import (
                list_plans, get_plan, delete_plan,
            )
            arg = cmd[len("/plans "):].strip() if cmd.startswith("/plans ") else ""
            repo = ctx.repo_dir or "."
            # /plans delete <name>
            if arg.startswith("delete "):
                target = arg[len("delete "):].strip()
                removed = delete_plan(repo, target)
                if removed is None:
                    _append_system_message(
                        f"No plan matching '{target}'. /plans for the list."
                    )
                else:
                    short = str(removed).replace(str(Path.home()), "~")
                    _append_system_message(f"🗑 Plan deleted → {short}")
                return True
            # /plans <name>: view body
            if arg:
                rec = get_plan(repo, arg)
                if rec is None:
                    _append_system_message(
                        f"No plan matching '{arg}'. /plans for the list."
                    )
                    return True
                short = str(rec["path"]).replace(str(Path.home()), "~")
                _append_system_message(
                    f"## Plan: {rec['name']}\n"
                    f"_{rec['description']}_\n\n"
                    f"---\n\n{rec['body']}\n\n---\n"
                    f"Source: `{short}`"
                )
                return True
            # /plans: list mode
            plans = list_plans(repo)
            if not plans:
                _append_system_message(
                    "No saved plans yet. Plans are written here when you "
                    "approve a Plan-Mode plan via ExitPlanMode."
                )
                return True
            import time as _time
            lines = ["Saved plans (newest first):"]
            for p in plans[:25]:
                ts = _time.strftime("%Y-%m-%d %H:%M",
                                    _time.localtime(p["created_at"]))
                desc = (p["description"] or "")[:60]
                lines.append(f"  {ts}  {p['name']:<32}  {desc}")
            lines.append(
                "\nUse `/plans <name>` to view a plan body, "
                "`/plans delete <name>` to remove one."
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/agents stats":
            try:
                from delfin.agent.subagents import (
                    read_telemetry, summarize_telemetry,
                )
                records = read_telemetry()
                stats = summarize_telemetry(records)
            except Exception as exc:
                _append_system_message(f"Could not read telemetry: {exc}")
                return True
            if not stats:
                _append_system_message(
                    "No subagent runs recorded yet. Invoke a subagent via "
                    "the `subagent` tool and check back."
                )
                return True
            lines = ["Subagent telemetry (lifetime totals):"]
            lines.append(
                f"  {'type':<18} {'n':>4} {'avg s':>7} {'in_tok':>8} "
                f"{'out_tok':>8} {'err':>4} {'trunc':>6}"
            )
            for t in sorted(stats):
                b = stats[t]
                avg_s = (b["elapsed_s_total"] / b["count"]) if b["count"] else 0.0
                lines.append(
                    f"  {t:<18} {b['count']:>4} {avg_s:>7.1f} "
                    f"{b['input_tokens_total']:>8} {b['output_tokens_total']:>8} "
                    f"{b['errors']:>4} {b['truncated']:>6}"
                )
            lines.append(
                f"\nLog: ~/.delfin/subagent_telemetry.jsonl "
                f"({len(records)} records)"
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/agents":
            try:
                from delfin.agent.subagents import list_subagents
                presets = list_subagents()
            except Exception as exc:
                _append_system_message(f"Could not list subagents: {exc}")
                return True
            if not presets:
                _append_system_message("No subagent presets registered.")
                return True
            lines = ["Subagent presets (use via the `subagent` tool):"]
            for p in presets:
                name = p.get("subagent_type") or p.get("name") or "?"
                desc = p.get("description") or ""
                lines.append(f"  {name:<18} {desc}")
            lines.append(
                "\nType `/agents stats` for lifetime cost/duration aggregates."
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/skills" or cmd.startswith("/skills "):
            try:
                from delfin.agent.skills import discover_skills, get_skill
                skills = discover_skills(ctx.repo_dir or None)
            except Exception as exc:
                _append_system_message(f"Could not list skills: {exc}")
                return True
            arg = cmd[len("/skills "):].strip() if cmd.startswith("/skills ") else ""
            # Detail view: /skills <name> prints the full body so the user
            # can read what the skill does without loading it into the input.
            if arg:
                sk = get_skill(arg, ctx.repo_dir or None)
                if sk is None:
                    _append_system_message(
                        f"Unknown skill '{arg}'. "
                        "Type /skills (no arg) for the full list."
                    )
                    return True
                src = str(getattr(sk, "source", "")).replace(
                    str(Path.home()), "~"
                )
                _append_system_message(
                    f"## /skill {sk.name}\n"
                    f"_{(sk.description or '').strip()}_\n\n"
                    f"---\n\n{sk.body.strip()}\n\n---\n"
                    f"Source: `{src}`\n"
                    f"Usage: `/skill {sk.name}` loads the body into your "
                    "input so you can edit it before sending."
                )
                return True
            # No arg → list mode
            if not skills:
                _append_system_message(
                    "No skills installed. Drop a SKILL.md under "
                    "~/.delfin/skills/<name>/ or use the built-ins in "
                    "delfin/agent/pack/skills/."
                )
                return True
            lines = ["Discovered skills (use `/skills <name>` to view body):"]
            for sk in skills:
                slash = f"/skill {sk.name}"
                desc = (sk.description or "").strip().splitlines()[0][:80] if sk.description else ""
                lines.append(f"  {slash:<28} {desc}")
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/tools" or cmd.startswith("/tools "):
            arg = cmd[7:].strip() if len(cmd) > 7 else ""
            # Categorise the native function-calling tools the agent has
            # access to. MCP tools live in their own namespace and are
            # huge in number, so we just report the count + how to
            # enumerate them per category via list_tools().
            _CATS = {
                "Filesystem": [
                    ("read_file", "Read a file"),
                    ("write_file", "Create / overwrite a file"),
                    ("edit_file", "Replace text in an existing file"),
                    ("multi_edit", "Multiple edits to one file in a single call"),
                    ("apply_patch", "Apply a unified diff"),
                    ("list_files", "List files matching a glob"),
                    ("grep_file", "Regex search inside files"),
                ],
                "Bash & jobs": [
                    ("bash", "Run a shell command (foreground, ≤120 s default)"),
                    ("bash_background", "Start a long-running command, returns job_id"),
                    ("bash_status", "Poll background job"),
                    ("bash_output", "Live head/tail of a background job's stdout/stderr"),
                    ("bash_kill", "SIGTERM then SIGKILL a background job"),
                    ("run_tests", "Run pytest on a target module/path"),
                ],
                "Web research": [
                    ("web_search", "DuckDuckGo search, returns top hits"),
                    ("web_fetch", "Fetch a single URL as plain text"),
                ],
                "Notebooks": [
                    ("notebook_read", "Read .ipynb cells (text+output summaries)"),
                    ("notebook_edit", "Insert / replace / delete cells safely"),
                ],
                "Sub-agent & planning": [
                    ("subagent", "Delegate to explore / plan / code-reviewer / general-purpose (+ md-defined)"),
                    ("exit_plan_mode", "Submit a plan for user approval (plan-mode only)"),
                    ("enter_worktree", "Create an isolated git worktree"),
                    ("exit_worktree", "Tear down a worktree"),
                ],
                "Task / Skill": [
                    ("task_create", "Add a todo (subject, description, active_form)"),
                    ("task_update", "Change a todo status (pending / in_progress / completed)"),
                    ("task_list", "List todos grouped by status"),
                    ("skill", "Invoke a discovered skill by name"),
                ],
                "User interaction": [
                    ("ask_user_question", "Render multi-choice question + collect reply"),
                    ("push_notification", "Send a desktop notification"),
                ],
                "Scheduling": [
                    ("schedule_wakeup", "Re-invoke this agent after N seconds"),
                    ("cron_create", "Persistent cron-style schedule"),
                    ("cron_list", "List active cron entries"),
                    ("cron_delete", "Remove a cron entry"),
                    ("remote_trigger", "Generic webhook-style trigger"),
                ],
                "Code navigation (DELFIN)": [
                    ("find_definition", "Jump to a symbol's definition"),
                    ("find_references", "Find all callers/usages of a symbol"),
                    ("project_introspect", "High-level repo layout dump"),
                ],
                "Docs & calcs (DELFIN)": [
                    ("search_docs", "Search indexed PDFs (ORCA manual, xTB, papers)"),
                    ("list_docs", "List available docs"),
                    ("list_sections", "List sections in a doc"),
                    ("read_section", "Read a specific section"),
                    ("search_calcs", "Find calculations in calc/ + archive/"),
                    ("calc_summary", "High-level overview of all calculations"),
                    ("get_calc_info", "Detailed info for one calculation"),
                ],
                "Permissions": [
                    ("remember_permission", "Persist a single allow/deny pattern"),
                    ("remember_permission_bundle", "Persist a profile (e.g. project_dev)"),
                ],
            }
            # Filter by category if arg passed
            if arg:
                filt = arg.lower()
                shown = {
                    cat: tools for cat, tools in _CATS.items()
                    if filt in cat.lower()
                    or any(filt in name for name, _ in tools)
                }
                if not shown:
                    _append_system_message(
                        f"No tools match {arg!r}. Try `/tools` for the full list."
                    )
                    return True
            else:
                shown = _CATS
            # Count subagent presets + MCP tools for the footer
            try:
                from delfin.agent.subagents import subagent_type_names
                sa_count = len(subagent_type_names())
            except Exception:
                sa_count = 0
            lines = ["Native tools (function-calling surface):"]
            total = 0
            for cat, tools in shown.items():
                lines.append(f"\n{cat}:")
                for name, desc in tools:
                    lines.append(f"  {name:<28} — {desc}")
                    total += 1
            lines.append(
                f"\n{total} native tools shown. "
                f"+ {sa_count} subagent presets (`/agents`). "
                "+ ~74 DELFIN ops MCP tools and ~10 doc-search MCP tools "
                "(`mcp__delfin-ops__list_tools(category='parsing')` etc.)."
            )
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/context":
            engine = state["engine"]
            if not engine:
                _append_system_message("No active engine. Send a message first.")
                return True
            n_msgs = len(getattr(engine, "messages", []) or [])
            try:
                tokens = int(engine._estimate_context_tokens())
            except Exception:
                tokens = 0
            window = int(getattr(engine, "context_window_tokens", 100_000) or 100_000)
            pct = (tokens / window * 100.0) if window else 0.0
            compact_pct = float(getattr(engine, "auto_compact_pct", 0.95) or 0.95)
            last_compact = getattr(engine, "last_compaction_info", None)
            if last_compact:
                lc = (
                    f"  Last compaction: {last_compact.get('messages_compacted', '?')} msgs, "
                    f"saved ~{last_compact.get('tokens_saved', 0)} tokens"
                )
            else:
                lc = "  Last compaction: (none this session)"
            _append_system_message(
                "Context status:\n"
                f"  Messages:    {n_msgs}\n"
                f"  Est. tokens: {tokens:,} ({pct:.1f}% of {window:,} window)\n"
                f"  Auto-compact trigger: ≥12 msgs OR ≥{compact_pct*100:.0f}% of window\n"
                f"{lc}\n"
                f"  Mode:        {engine.mode}\n"
                f"  Model:       {model_dropdown.value}\n"
                f"  Provider:    {provider_dropdown.value}\n"
                f"  Profile:     {state.get('_perm_profile', 'ask_all')}"
            )
            return True

        if cmd == "/undo":
            engine = state["engine"]
            if not engine or len(engine.messages) < 2:
                _append_system_message("Nothing to undo.")
                return True
            # Remove the last assistant + user message pair
            removed = 0
            while engine.messages and removed < 2:
                last = engine.messages[-1]
                engine.messages.pop()
                removed += 1
                if last["role"] == "user":
                    break
            # Also remove last chat messages from display
            msgs = state["chat_messages"]
            to_remove = 0
            for m in reversed(msgs):
                if m["role"] in ("assistant", "tool", "system"):
                    to_remove += 1
                elif m["role"] == "user":
                    to_remove += 1
                    break
                else:
                    break
            if to_remove:
                state["chat_messages"] = msgs[:-to_remove]
                _refresh_chat_html()
            _append_system_message(
                f"Undone last turn ({removed} engine messages removed)."
            )
            return True

        # /git commands
        _repo_dir = str(ctx.repo_dir or ".")
        if cmd == "/git status":
            try:
                r = _sp.run(
                    ["git", "status", "--short"],
                    capture_output=True, text=True, cwd=_repo_dir, timeout=10,
                )
                output = r.stdout.strip() or "Working tree clean."
                _append_system_message(f"git status:\n{output}")
            except Exception as e:
                _append_system_message(f"git status error: {e}")
            return True

        if cmd == "/git diff":
            try:
                r = _sp.run(
                    ["git", "diff", "--stat"],
                    capture_output=True, text=True, cwd=_repo_dir, timeout=10,
                )
                output = r.stdout.strip() or "No unstaged changes."
                # Also check staged
                r2 = _sp.run(
                    ["git", "diff", "--cached", "--stat"],
                    capture_output=True, text=True, cwd=_repo_dir, timeout=10,
                )
                staged = r2.stdout.strip()
                full = output
                if staged:
                    full += f"\n\nStaged:\n{staged}"
                _append_system_message(f"git diff:\n{full}")
            except Exception as e:
                _append_system_message(f"git diff error: {e}")
            return True

        if cmd == "/git log":
            try:
                r = _sp.run(
                    ["git", "log", "--oneline", "-15"],
                    capture_output=True, text=True, cwd=_repo_dir, timeout=10,
                )
                output = r.stdout.strip() or "No commits."
                _append_system_message(f"git log (last 15):\n{output}")
            except Exception as e:
                _append_system_message(f"git log error: {e}")
            return True

        if cmd == "/git branch":
            try:
                r = _sp.run(
                    ["git", "branch", "-a", "--no-color"],
                    capture_output=True, text=True, cwd=_repo_dir, timeout=10,
                )
                output = r.stdout.strip() or "No branches."
                _append_system_message(f"git branch:\n{output}")
            except Exception as e:
                _append_system_message(f"git branch error: {e}")
            return True

        # /provider <name>
        if cmd.startswith("/provider "):
            name = cmd[10:].strip().lower()
            if name in ("claude", "openai", "kit"):
                provider_dropdown.value = name
                _append_system_message(f"Provider switched to {name}.")
            else:
                _append_system_message(
                    f"Unknown provider '{name}'. Options: claude, openai, kit"
                )
            return True

        # /model <name>
        if cmd.startswith("/model "):
            name = cmd[7:].strip()
            valid = {v for _, v in model_dropdown.options}
            if name in valid:
                model_dropdown.value = name
                _append_system_message(f"Model switched to {name}.")
            else:
                _append_system_message(
                    f"Unknown model '{name}'. Options: {', '.join(sorted(valid))}"
                )
            return True

        # /effort <level>  — accepts low/medium/high/xhigh
        if cmd.startswith("/effort "):
            level = cmd[8:].strip()
            valid = {"low", "medium", "high", "xhigh"}
            if level in valid:
                effort_dropdown.value = level
                _append_system_message(f"Effort set to {level}.")
            else:
                _append_system_message(
                    f"Unknown effort '{level}'. Options: {', '.join(sorted(valid))}"
                )
            return True

        # /mode <name>
        if cmd.startswith("/mode "):
            name = cmd[6:].strip()
            opts = _dropdown_values(
                mode_dropdown.options if hasattr(mode_dropdown, "options") else []
            )
            if not opts:
                opts = [mode_dropdown.value]
            if name in opts:
                mode_dropdown.value = name
                _append_system_message(f"Mode switched to {name}.")
            else:
                _append_system_message(
                    f"Unknown mode '{name}'. Options: {', '.join(opts)}"
                )
            return True

        # /perms — switch permission profile
        if cmd == "/perms" or cmd.startswith("/perms "):
            arg = cmd[6:].strip() if len(cmd) > 6 else ""
            valid = list(_PERM_PROFILES.keys())
            if not arg:
                cur = state.get("_perm_profile", "ask_all")
                lines = [f"Permission profile: **{cur}**\n"]
                _labels = {
                    "plan": "Read-only — agent can only browse and analyze",
                    "ask_all": "Ask for everything that changes",
                    "repo_free": "Repo free, calc needs confirmation",
                    "all_free": "All free (archive + remote archive always read-only)",
                }
                for pid in valid:
                    marker = " ◀" if pid == cur else ""
                    lines.append(f"  {pid:10s} — {_labels.get(pid, '')}{marker}")
                lines.append(f"\nUsage: /perms <{'/'.join(valid)}>")
                _append_system_message("\n".join(lines))
            elif arg in valid:
                state["_perm_profile"] = arg
                perm_dropdown.value = arg  # sync dropdown → triggers _on_perm_change
            else:
                _append_system_message(
                    f"Unknown profile '{arg}'. Options: {', '.join(valid)}"
                )
            return True

        # /export — save chat as Markdown file
        if cmd == "/export":
            _export_chat()
            return True

        # /search <text>
        if cmd.startswith("/search"):
            query = text[7:].strip()  # use original case
            if query:
                _do_search(query)
            else:
                # Toggle search bar visibility
                _toggle_search()
            return True

        # /retry — regenerate last response
        if cmd == "/retry":
            _retry_last()
            return True

        # Internal: Shift+Tab permission cycling (not shown in /help)
        if cmd == "/perm-cycle":
            if state["streaming"]:
                _append_system_message("Cannot change permissions while streaming.")
                return True
            perm_values = [v for _, v in perm_dropdown.options]
            perm_labels = {v: label for label, v in perm_dropdown.options}
            idx = perm_values.index(perm_dropdown.value) if perm_dropdown.value in perm_values else 0
            next_idx = (idx + 1) % len(perm_values)
            perm_dropdown.value = perm_values[next_idx]
            _append_system_message(
                f"Permission: {perm_labels[perm_values[next_idx]]} (Shift+Tab to cycle)"
            )
            return True

        # -- Memory commands ---------------------------------------------------

        if cmd == "/memories verify" or cmd == "/memories check":
            try:
                from delfin.agent.memory_store import verify_typed_memories
                stale = verify_typed_memories(ctx.repo_dir or ".")
            except Exception as exc:
                _append_system_message(f"Memory verify failed: {exc}")
                return True
            if not stale:
                _append_system_message(
                    "All memory references resolve — no stale paths found."
                )
            else:
                lines = ["Stale references found in memory:"]
                for entry in stale:
                    lines.append(f"  {entry['file']}:")
                    for ref in entry["stale_refs"][:8]:
                        lines.append(f"    - {ref}")
                lines.append(
                    "Consider /forget <idx> for outdated entries, or edit "
                    "the .md files in ~/.claude/projects/<slug>/memory/."
                )
                _append_system_message("\n".join(lines))
            return True

        if cmd == "/memories":
            from delfin.agent.memory_store import load_memories
            facts = load_memories()
            if not facts:
                _append_system_message("No memories stored. Use /remember <text> to add one.")
            else:
                lines = []
                for i, f in enumerate(facts):
                    lines.append(f"  [{i}] ({f.get('source', '?')}) {f['text']}")
                _append_system_message("Agent memories:\n" + "\n".join(lines))
            return True

        if cmd.startswith("/remember "):
            text_to_save = text[len("/remember "):].strip()
            if text_to_save:
                from delfin.agent.memory_store import (
                    save_memory,
                    save_typed_memory,
                    parse_memory_type,
                )
                mem_type, body = parse_memory_type(text_to_save)
                # Always keep a flat JSON record for legacy retrieval.
                idx = save_memory(body or text_to_save, source=mem_type)
                # Mirror into the .delfin typed-memory layout so the
                # next session picks it up via the prompt loader.
                try:
                    fpath, slug, _t = save_typed_memory(
                        text_to_save,
                        repo_root=ctx.repo_dir or ".",
                    )
                    short = str(fpath).replace(str(Path.home()), "~")
                    _append_system_message(
                        f"Memory [{idx}] saved as {mem_type} → {short}"
                    )
                except Exception as exc:
                    _append_system_message(
                        f"Memory [{idx}] saved (flat); typed-mirror failed: {exc}"
                    )
            else:
                _append_system_message(
                    "Usage: /remember [user|feedback|project|reference:] <text>"
                )
            return True

        if cmd.startswith("/forget"):
            arg = text[7:].strip()
            if arg.isdigit():
                from delfin.agent.memory_store import delete_memory
                if delete_memory(int(arg)):
                    _append_system_message(f"Memory [{arg}] deleted.")
                else:
                    _append_system_message(f"Invalid index: {arg}")
            else:
                _append_system_message("Usage: /forget <index>")
            return True

        # -- Workspace commands ------------------------------------------------

        if cmd.startswith("/workspace"):
            ws_cmd = text[10:].strip().lower()
            ws_dir = ctx.agent_dir
            if ws_cmd == "ls" or ws_cmd == "":
                try:
                    items = sorted(ws_dir.iterdir())
                    if items:
                        lines = []
                        for p in items[:50]:
                            kind = "dir" if p.is_dir() else f"{p.stat().st_size:,}B"
                            lines.append(f"  {p.name}  ({kind})")
                        _append_system_message(
                            f"Agent workspace ({ws_dir}):\n" + "\n".join(lines)
                        )
                    else:
                        _append_system_message(f"Agent workspace is empty: {ws_dir}")
                except Exception as e:
                    _append_system_message(f"Workspace error: {e}")
            elif ws_cmd == "clean":
                import shutil
                for p in ws_dir.iterdir():
                    if p.is_dir():
                        shutil.rmtree(p, ignore_errors=True)
                    else:
                        p.unlink(missing_ok=True)
                _append_system_message("Agent workspace cleaned.")
            elif ws_cmd.startswith("read "):
                fname = text[10:].strip()[5:].strip()
                fpath = (ws_dir / fname).resolve()
                if not str(fpath).startswith(str(ws_dir)):
                    _append_system_message("Path outside workspace.")
                elif fpath.is_file():
                    content = fpath.read_text(encoding="utf-8", errors="replace")[:8000]
                    _append_system_message(f"{fname}:\n{content}")
                else:
                    _append_system_message(f"File not found: {fname}")
            else:
                _append_system_message(
                    "Workspace commands: /workspace ls, /workspace read <file>, /workspace clean"
                )
            return True

        # -- Dashboard control commands ----------------------------------------

        # /tab <name> — navigate to a dashboard tab
        if cmd.startswith("/tab "):
            tab_name = text[5:].strip()
            _tab_map = {
                "submit": "Submit Job",
                "recalc": "Recalc",
                "jobs": "Job Status",
                "orca": "ORCA Builder",
                "calc": "Calculations",
                # Common aliases — accept English plural/full forms and
                # German equivalents so the agent doesn't have to guess.
                "calcs": "Calculations",
                "calculation": "Calculations",
                "calculations": "Calculations",
                "berechnung": "Calculations",
                "berechnungen": "Calculations",
                "archive": "Archive",
                "archiv": "Archive",
                "literature": "Literature",
                "literatur": "Literature",
                "agent": "Agent",
                "settings": "Settings",
                "einstellungen": "Settings",
                "job": "Job Status",
                "joblist": "Job Status",
                "submit_job": "Submit Job",
                "submitjob": "Submit Job",
                "orca_builder": "ORCA Builder",
                "orcabuilder": "ORCA Builder",
            }
            # Canonical short keys to surface in error messages — keep the
            # alias map fat but the help short.
            _CANONICAL_TAB_KEYS = (
                "submit", "recalc", "jobs", "orca", "calc",
                "archive", "literature", "agent", "settings",
            )
            tab_key = tab_name.lower()
            title = _tab_map.get(tab_key)
            if title is None:
                # Fuzzy match: typos like "claultions" / "submmt" / "joobs"
                # should still land on the right tab instead of hard-erroring
                # and burning a roundtrip for the agent to retry. Uses
                # SequenceMatcher ratio with a generous 0.6 cut-off (matches
                # 2-character transpositions / drops on 8-12 letter names).
                import difflib as _difflib
                candidates = list(_tab_map.keys()) + list(_tab_map.values())
                scored = sorted(
                    (
                        (
                            _difflib.SequenceMatcher(
                                None, tab_key, c.lower()
                            ).ratio(),
                            c,
                        )
                        for c in candidates
                    ),
                    reverse=True,
                )
                if scored and scored[0][0] >= 0.6:
                    best = scored[0][1]
                    title = _tab_map.get(best.lower(), best)
                    _append_system_message(
                        f"(Fuzzy-matched '{tab_name}' → '{title}')"
                    )
                else:
                    title = tab_name  # falls through to the error path below
            if ctx.select_tab(title):
                _append_system_message(f"Switched to tab: {title}")
            else:
                avail = ", ".join(_CANONICAL_TAB_KEYS)
                _append_system_message(
                    f"Tab not found: {tab_name}\nAvailable: {avail}"
                )
            return True

        # /ui <widget> <property> [value] — modify dashboard widgets (session-only)
        # /ui list — show all available widgets
        if cmd == "/ui list":
            _all_w = _build_full_widget_registry()
            lines = []
            for name in sorted(_all_w):
                w = _all_w[name]
                wtype = type(w).__name__
                val_preview = ""
                if hasattr(w, "value"):
                    v = str(w.value)
                    val_preview = f" = {v[:40]}{'...' if len(v) > 40 else ''}"
                lines.append(f"  {name} ({wtype}){val_preview}")
            _append_system_message("Available widgets:\n" + "\n".join(lines))
            return True

        if cmd.startswith("/ui "):
            _all_w = _build_full_widget_registry()
            _ui_parts = text[4:].strip().split(None, 2)
            if len(_ui_parts) < 2:
                _append_system_message(
                    "Usage: /ui <widget> <property> [value]\n"
                    "       /ui list — show all widgets\n"
                    f"Widgets: {', '.join(sorted(_all_w))}\n"
                    "Properties:\n"
                    "  click          — press a button\n"
                    "  value <v>      — set widget value (text, number, dropdown selection)\n"
                    "  replace <old> <new> — find & replace text in widget value\n"
                    "  show           — show current properties & value\n"
                    "  options        — show dropdown options\n"
                    "  style <s>      — button style (primary|danger|success|info|warning)\n"
                    "  text <v>       — description / placeholder text\n"
                    "  disabled true|false\n"
                    "  visible true|false\n"
                    "  width <css>    — e.g. 120px\n"
                    "  height <css>   — e.g. 40px"
                )
                return True
            _ui_wname = _ui_parts[0].lower()
            _ui_prop = _ui_parts[1].lower()
            _ui_val = _ui_parts[2] if len(_ui_parts) > 2 else ""
            _ui_w = _all_w.get(_ui_wname)
            if _ui_w is None:
                _append_system_message(
                    f"Unknown widget: {_ui_wname}\n"
                    f"Available: {', '.join(sorted(_all_w))}"
                )
                return True

            # show — display current properties
            if _ui_prop == "show":
                _props = []
                if hasattr(_ui_w, "value"):
                    v = str(_ui_w.value)
                    _props.append(f"value: {v[:200]}{'...' if len(v) > 200 else ''}")
                if hasattr(_ui_w, "description"):
                    _props.append(f"text: {_ui_w.description}")
                if hasattr(_ui_w, "button_style"):
                    _props.append(f"style: {_ui_w.button_style}")
                if hasattr(_ui_w, "options"):
                    opts = list(_ui_w.options) if _ui_w.options else []
                    _props.append(f"options: {opts[:20]}{'...' if len(opts) > 20 else ''}")
                if hasattr(_ui_w, "disabled"):
                    _props.append(f"disabled: {_ui_w.disabled}")
                if hasattr(_ui_w, "layout"):
                    _props.append(f"width: {_ui_w.layout.width}")
                    _props.append(f"height: {_ui_w.layout.height}")
                _props.append(f"visible: {getattr(_ui_w.layout, 'display', '') != 'none'}")
                _props.append(f"type: {type(_ui_w).__name__}")
                _append_system_message(
                    f"Widget '{_ui_wname}':\n" + "\n".join(f"  {p}" for p in _props)
                )
                return True

            # options — show dropdown options
            if _ui_prop == "options":
                if hasattr(_ui_w, "options"):
                    opts = list(_ui_w.options)
                    _append_system_message(
                        f"'{_ui_wname}' options ({len(opts)}):\n"
                        + "\n".join(f"  {o}" for o in opts)
                    )
                else:
                    _append_system_message(f"Widget '{_ui_wname}' has no options.")
                return True

            # click — press a button
            if _ui_prop == "click":
                if not hasattr(_ui_w, "click"):
                    _append_system_message(f"Widget '{_ui_wname}' is not a button.")
                    return True
                # Hard-blocked buttons: agent can never click these
                if _ui_wname in _BLOCKED_WIDGETS:
                    _append_system_message(f"⛔ Blocked: '{_ui_wname}' cannot be clicked by the agent.")
                    return True
                # Confirmation-required buttons: ask user before executing
                if _ui_wname in _CONFIRM_WIDGETS:
                    _readable = _ui_wname.replace("-", " ").replace("btn", "").strip()
                    _request_confirmation(
                        f"ui-click-{_ui_wname}",
                        f"Agent wants to click **{_readable}**",
                        lambda: (_ui_w.click(), _append_system_message(f"'{_ui_wname}' clicked (confirmed).")),
                    )
                    return True
                _ui_w.click()
                _append_system_message(f"'{_ui_wname}' clicked.")
                return True

            # value — set the widget's value (text, int, dropdown selection)
            if _ui_prop == "value":
                if not hasattr(_ui_w, "value"):
                    _append_system_message(f"Widget '{_ui_wname}' has no value property.")
                    return True
                # Dropdown: validate against options
                if hasattr(_ui_w, "options") and _ui_w.options:
                    opts = list(_ui_w.options)
                    # Check for tuple options like (label, value)
                    opt_values = [o[1] if isinstance(o, tuple) else o for o in opts]
                    opt_labels = [o[0] if isinstance(o, tuple) else o for o in opts]
                    # Try match by value or label (case-insensitive)
                    matched = None
                    for ov, ol in zip(opt_values, opt_labels):
                        if str(ov).lower() == _ui_val.lower() or str(ol).lower() == _ui_val.lower():
                            matched = ov
                            break
                    if matched is None:
                        _append_system_message(
                            f"Invalid option: {_ui_val}\n"
                            f"Available: {', '.join(str(o) for o in opt_labels)}"
                        )
                        return True
                    _ui_w.value = matched
                    _append_system_message(f"'{_ui_wname}' -> {matched}")
                # IntText / BoundedIntText
                elif isinstance(_ui_w, (widgets.IntText, widgets.BoundedIntText)):
                    try:
                        _ui_w.value = int(_ui_val)
                        _append_system_message(f"'{_ui_wname}' -> {_ui_w.value}")
                    except ValueError:
                        _append_system_message(f"'{_ui_wname}' expects an integer, got: {_ui_val}")
                # FloatText
                elif isinstance(_ui_w, (widgets.FloatText, widgets.BoundedFloatText)):
                    try:
                        _ui_w.value = float(_ui_val)
                        _append_system_message(f"'{_ui_wname}' -> {_ui_w.value}")
                    except ValueError:
                        _append_system_message(f"'{_ui_wname}' expects a number, got: {_ui_val}")
                # Text / Textarea — convert escaped \n to real newlines
                else:
                    _final_val = _ui_val.replace("\\n", "\n") if isinstance(_ui_w, widgets.Textarea) else _ui_val
                    _ui_w.value = _final_val
                    _preview = _final_val.replace("\n", " ")[:80]
                    _append_system_message(f"'{_ui_wname}' -> {_preview}{'...' if len(_final_val) > 80 else ''}")
                return True

            # replace <old> <new> — find & replace text in a widget's value
            if _ui_prop == "replace":
                if not hasattr(_ui_w, "value") or not isinstance(_ui_w.value, str):
                    _append_system_message(f"Widget '{_ui_wname}' has no text value to replace in.")
                    return True
                # _ui_val contains "old_text new_text" — split on first whitespace
                _r_parts = _ui_val.split(None, 1)
                if len(_r_parts) < 2:
                    _append_system_message("Usage: /ui <widget> replace <old> <new>")
                    return True
                _r_old, _r_new = _r_parts
                if _r_old not in _ui_w.value:
                    _append_system_message(f"'{_r_old}' not found in '{_ui_wname}'.")
                    return True
                _ui_w.value = _ui_w.value.replace(_r_old, _r_new, 1)
                _append_system_message(f"'{_ui_wname}': replaced '{_r_old}' with '{_r_new}'")
                return True

            # style — button color
            if _ui_prop == "style":
                if not hasattr(_ui_w, "button_style"):
                    _append_system_message(f"Widget '{_ui_wname}' has no button_style.")
                    return True
                _valid_styles = ("primary", "danger", "success", "info", "warning", "")
                if _ui_val not in _valid_styles:
                    _append_system_message(
                        f"Invalid style: {_ui_val}\n"
                        f"Valid: {', '.join(s or '(empty)' for s in _valid_styles)}"
                    )
                    return True
                _ui_w.button_style = _ui_val
                _append_system_message(f"'{_ui_wname}' style -> {_ui_val or '(default)'}")
                return True

            # text — description / placeholder
            if _ui_prop == "text":
                if hasattr(_ui_w, "description"):
                    _ui_w.description = _ui_val
                    _append_system_message(f"'{_ui_wname}' text -> {_ui_val}")
                elif hasattr(_ui_w, "placeholder"):
                    _ui_w.placeholder = _ui_val
                    _append_system_message(f"'{_ui_wname}' placeholder -> {_ui_val}")
                else:
                    _append_system_message(f"Widget '{_ui_wname}' has no text property.")
                return True

            # disabled
            if _ui_prop == "disabled":
                _ui_w.disabled = _ui_val.lower() in ("true", "1", "yes")
                _append_system_message(f"'{_ui_wname}' disabled -> {_ui_w.disabled}")
                return True

            # visible
            if _ui_prop == "visible":
                _show = _ui_val.lower() in ("true", "1", "yes")
                _ui_w.layout.display = "" if _show else "none"
                _append_system_message(f"'{_ui_wname}' visible -> {_show}")
                return True

            # width / height
            if _ui_prop in ("width", "height"):
                setattr(_ui_w.layout, _ui_prop, _ui_val)
                _append_system_message(f"'{_ui_wname}' {_ui_prop} -> {_ui_val}")
                return True

            _append_system_message(
                f"Unknown property: {_ui_prop}\n"
                "Valid: click, value, replace, show, options, style, text, disabled, visible, width, height"
            )
            return True

        # /control show — show current CONTROL content from Submit tab
        if cmd == "/control show":
            cw = ctx.submit_refs.get("control_widget")
            if cw:
                val = cw.value.strip()
                if val:
                    _append_system_message(f"Current CONTROL content:\n```\n{val}\n```")
                else:
                    _append_system_message("CONTROL widget is empty.")
            else:
                _append_system_message("Submit tab not available.")
            return True

        # /control set <content> — set CONTROL content in Submit tab
        if cmd.startswith("/control set "):
            content = text[len("/control set "):].strip()
            cw = ctx.submit_refs.get("control_widget")
            if cw is None:
                _append_system_message("Submit tab not available.")
                return True
            old_value = (cw.value or "").strip()
            # Guard: if the agent passed a single 'key=value' line while the
            # existing CONTROL is multi-line, this is almost always a mis-call
            # of '/control set' instead of '/control key K V'. Refuse and
            # explain so the agent learns the difference instead of silently
            # wiping the file.
            looks_single_kv = (
                "\n" not in content
                and "=" in content
                and len(content) < 80
            )
            existing_is_multiline = old_value.count("\n") >= 1
            if looks_single_kv and existing_is_multiline:
                hint_key = content.split("=", 1)[0].strip()
                hint_val = content.split("=", 1)[1].strip()
                _append_system_message(
                    f"`/control set` was called with a single "
                    f"`{hint_key}={hint_val}` and would have overwritten "
                    f"the existing multi-line CONTROL. Refusing — use "
                    f"`/control key {hint_key} {hint_val}` to change one "
                    "keyword without clobbering the rest. If you really "
                    "want to replace everything, send the full multi-line "
                    "CONTROL via `/control set`."
                )
                return True
            cw.value = content
            _append_system_message(
                f"CONTROL updated ({len(content)} chars). "
                "Switch to Submit tab to review."
            )
            return True

        # /control key <key> <value> — change a single key in CONTROL
        if cmd.startswith("/control key "):
            parts = text[len("/control key "):].strip().split(None, 1)
            if len(parts) < 2:
                _append_system_message("Usage: /control key <key> <value>")
                return True
            key, value = parts[0], parts[1]
            cw = ctx.submit_refs.get("control_widget")
            if not cw:
                _append_system_message("Submit tab not available.")
                return True
            import re as _re
            old = cw.value
            # Fuzzy-match the key against keys already present in CONTROL,
            # so typos like "fucntional" land on "functional" instead of
            # silently appending a new bogus key. Only kicks in when the
            # literal key is NOT already in the file — preserves the
            # "add new key" capability for genuinely new keys.
            existing_keys = _re.findall(r"^([A-Za-z_][\w-]*)\s*=", old,
                                         flags=_re.MULTILINE)
            if existing_keys and key.lower() not in {k.lower() for k in existing_keys}:
                import difflib as _difflib
                scored = sorted(
                    (
                        (
                            _difflib.SequenceMatcher(
                                None, key.lower(), ek.lower()
                            ).ratio(),
                            ek,
                        )
                        for ek in existing_keys
                    ),
                    reverse=True,
                )
                if scored and scored[0][0] >= 0.75:
                    best = scored[0][1]
                    _append_system_message(
                        f"(Fuzzy-matched CONTROL key '{key}' → '{best}')"
                    )
                    key = best
            # Match key=... (case-insensitive key match)
            pattern = _re.compile(
                r"^(" + _re.escape(key) + r")\s*=.*$",
                _re.MULTILINE | _re.IGNORECASE,
            )
            if pattern.search(old):
                new = pattern.sub(f"{key}={value}", old)
                cw.value = new
                _append_system_message(f"CONTROL: {key}={value}")
            else:
                # Key not found — append it
                cw.value = old.rstrip() + f"\n{key}={value}\n"
                _append_system_message(f"CONTROL: added {key}={value}")
            return True

        # /control validate — validate CONTROL content
        if cmd == "/control validate":
            cw = ctx.submit_refs.get("control_widget")
            if not cw:
                _append_system_message("Submit tab not available.")
                return True
            try:
                from delfin.dashboard.tab_submit import validate_control_text, get_esd_hints
                errors = validate_control_text(cw.value)
                lines = []
                if errors:
                    lines.append("CONTROL validation failed:")
                    for err in errors:
                        lines.append(f"  - {err}")
                else:
                    lines.append("CONTROL looks valid.")
                hints = get_esd_hints(cw.value)
                if hints:
                    lines.append("ESD hints:")
                    for h in hints:
                        lines.append(f"  i {h}")
                _append_system_message("\n".join(lines))
            except Exception as exc:
                _append_system_message(f"Validation error: {exc}")
            return True

        # /submit — trigger job submission from Submit tab
        if cmd == "/submit":
            submit_fn = ctx.submit_refs.get("handle_submit")
            if submit_fn:
                jn = ctx.submit_refs.get("job_name_widget")
                job_name = jn.value.strip() if jn else "?"
                def _do_submit():
                    try:
                        submit_fn(None)
                        _append_system_message("Job submitted. Check Job Status tab.")
                        ctx.select_tab("Job Status")
                    except Exception as exc:
                        _append_system_message(f"Submit error: {exc}")
                _confirm_or_exec("submit_job", f"Submit job '{job_name}' from Submit tab?", _do_submit, cmd_for_zone="/submit")
            else:
                _append_system_message("Submit tab not available.")
            return True

        # /orca show — show ORCA Builder settings
        if cmd == "/orca show":
            refs = ctx.orca_builder_refs
            if refs:
                method = refs.get("orca_method")
                basis = refs.get("orca_basis")
                job_type = refs.get("orca_job_type")
                charge = refs.get("orca_charge")
                mult = refs.get("orca_multiplicity")
                disp = refs.get("orca_dispersion")
                solv = refs.get("orca_solvent")
                pal = refs.get("orca_pal")
                maxcore = refs.get("orca_maxcore")
                preview = refs.get("orca_preview")
                info = (
                    f"ORCA Builder settings:\n"
                    f"  Method:       {method.value if method else '?'}\n"
                    f"  Job type:     {job_type.value if job_type else '?'}\n"
                    f"  Basis:        {basis.value if basis else '?'}\n"
                    f"  Charge:       {charge.value if charge else '?'}\n"
                    f"  Multiplicity: {mult.value if mult else '?'}\n"
                    f"  Dispersion:   {disp.value if disp else '?'}\n"
                    f"  Solvent:      {solv.value if solv else '?'}\n"
                    f"  PAL:          {pal.value if pal else '?'}\n"
                    f"  Maxcore:      {maxcore.value if maxcore else '?'}\n"
                )
                if preview and preview.value.strip():
                    info += f"\nInput preview:\n```\n{preview.value.strip()}\n```"
                _append_system_message(info)
            else:
                _append_system_message("ORCA Builder not available.")
            return True

        # /orca set <param> <value> — set ORCA Builder parameter
        if cmd.startswith("/orca set "):
            parts = text[len("/orca set "):].strip().split(None, 1)
            if len(parts) < 2:
                _append_system_message("Usage: /orca set <param> <value>\nParams: method, basis, job_type, charge, mult, dispersion, solvent, pal, maxcore, coords")
                return True
            param, value = parts[0].lower(), parts[1]
            refs = ctx.orca_builder_refs
            if not refs:
                _append_system_message("ORCA Builder not available.")
                return True
            _orca_param_map = {
                "method": "orca_method",
                "basis": "orca_basis",
                "job_type": "orca_job_type",
                "charge": "orca_charge",
                "mult": "orca_multiplicity",
                "multiplicity": "orca_multiplicity",
                "dispersion": "orca_dispersion",
                "solvent": "orca_solvent",
                "pal": "orca_pal",
                "maxcore": "orca_maxcore",
                "coords": "orca_coords",
            }
            widget_key = _orca_param_map.get(param)
            if not widget_key:
                # Fuzzy-match: typos like "metod", "basiss", "soltent",
                # "dispersio", "multipliciti" should still land on the
                # right ORCA Builder parameter instead of erroring.
                import difflib as _difflib
                scored = sorted(
                    (
                        (
                            _difflib.SequenceMatcher(
                                None, param, alias
                            ).ratio(),
                            alias,
                        )
                        for alias in _orca_param_map.keys()
                    ),
                    reverse=True,
                )
                if scored and scored[0][0] >= 0.7:
                    best = scored[0][1]
                    _append_system_message(
                        f"(Fuzzy-matched ORCA param '{param}' → '{best}')"
                    )
                    param = best
                    widget_key = _orca_param_map[param]
                else:
                    _append_system_message(
                        f"Unknown param: {param}\n"
                        f"Available: {', '.join(_orca_param_map.keys())}"
                    )
                    return True
            w = refs.get(widget_key)
            if w:
                try:
                    if hasattr(w, "options") and isinstance(value, str):
                        # Dropdown — try exact match first, then case-insensitive
                        opt_values = [v for _, v in w.options] if isinstance(w.options[0], tuple) else list(w.options)
                        if value in opt_values:
                            w.value = value
                        else:
                            match = next((v for v in opt_values if v.lower() == value.lower()), None)
                            if match:
                                w.value = match
                            else:
                                w.value = value
                    elif isinstance(w.value, int):
                        w.value = int(value)
                    else:
                        w.value = value
                    _append_system_message(f"ORCA Builder: {param} = {w.value}")
                    # Refresh preview
                    update_fn = refs.get("update_orca_preview")
                    if update_fn:
                        try:
                            update_fn()
                        except Exception:
                            pass
                except Exception as exc:
                    _append_system_message(f"Error setting {param}: {exc}")
            else:
                _append_system_message(f"Widget not available: {param}")
            return True

        # /orca submit — submit ORCA job (requires confirmation)
        if cmd == "/orca submit":
            btn = ctx.orca_builder_refs.get("orca_submit_btn")
            if btn:
                def _do_orca_submit():
                    btn.click()
                    _append_system_message("ORCA job submitted. Check Job Status tab.")
                _confirm_or_exec(
                    "orca_submit", "Submit ORCA job?", _do_orca_submit,
                    cmd_for_zone="/orca submit",
                )
            else:
                _append_system_message("ORCA Builder not available.")
            return True

        # /jobs — show job status
        if cmd == "/jobs":
            refresh = ctx.job_status_refs.get("refresh_job_list")
            if refresh:
                try:
                    refresh()
                except Exception:
                    pass
            ctx.select_tab("Job Status")
            _append_system_message("Switched to Job Status tab.")
            return True

        # /jobs check — D5: manually trigger the proactive watcher once
        if cmd == "/jobs check":
            try:
                emitted = _check_job_events_once()
            except Exception as exc:
                _append_system_message(f"Job event check failed: {exc}")
                return True
            if emitted == 0:
                _append_system_message(
                    "No job-state changes since the last check."
                )
            return True

        # -- Calculations browsing (all SAFE) --------------------------------

        def _resolve_calc_path(subpath: str) -> Path:
            """Resolve a path relative to agent's current calc dir."""
            base = ctx.calc_dir
            agent_rel = state.get("_agent_calc_path", "")
            if agent_rel:
                base = base / agent_rel
            if subpath:
                target = base / subpath
            else:
                target = base
            # Security: don't escape calc_dir
            try:
                target.resolve().relative_to(ctx.calc_dir.resolve())
            except ValueError:
                return ctx.calc_dir
            return target

        if cmd == "/calc ls" or cmd.startswith("/calc ls "):
            subpath = text[len("/calc ls"):].strip()
            target = _resolve_calc_path(subpath)
            if not target.exists():
                _append_system_message(f"Not found: {target}")
                return True
            if target.is_file():
                _append_system_message(f"{target.name} ({target.stat().st_size:,} bytes)")
                return True
            items = sorted(target.iterdir(), key=lambda p: (p.is_file(), p.name))
            lines = []
            for p in items[:100]:
                if p.is_dir():
                    lines.append(f"  [DIR]  {p.name}/")
                else:
                    sz = p.stat().st_size
                    lines.append(f"  {sz:>10,}  {p.name}")
            if len(items) > 100:
                lines.append(f"  ... and {len(items) - 100} more")
            rel = target.relative_to(ctx.calc_dir) if target != ctx.calc_dir else Path(".")
            listing = "\n".join(lines) or "(empty)"
            _append_system_message(f"calc/{rel}:\n{listing}")
            return True

        if cmd.startswith("/calc cd "):
            subpath = text[len("/calc cd"):].strip()
            if subpath in ("..", ".."):
                cur = state.get("_agent_calc_path", "")
                state["_agent_calc_path"] = str(Path(cur).parent) if cur else ""
            elif subpath == "/":
                state["_agent_calc_path"] = ""
            else:
                target = _resolve_calc_path(subpath)
                if target.is_dir():
                    try:
                        rel = str(target.resolve().relative_to(ctx.calc_dir.resolve()))
                        state["_agent_calc_path"] = rel if rel != "." else ""
                    except ValueError:
                        state["_agent_calc_path"] = ""
                else:
                    _append_system_message(f"Not a directory: {subpath}")
                    return True
            cur = state.get("_agent_calc_path", "") or "."
            # Sync the actual browser widget to the same directory
            _path_w = ctx.calc_browser_refs.get("calc_path_input")
            if _path_w:
                try:
                    _full = str(ctx.calc_dir / cur) if cur != "." else str(ctx.calc_dir)
                    _path_w.value = _full
                except Exception:
                    pass
            _append_system_message(f"calc dir: calc/{cur}")
            return True

        # /calc select|open <filename> — open/select a file in the browser
        # so it appears in the Calculations preview pane.
        if cmd.startswith("/calc select ") or cmd.startswith("/calc open "):
            filename = (
                text[len("/calc select"):].strip()
                if cmd.startswith("/calc select ")
                else text[len("/calc open"):].strip()
            )
            _path_w = ctx.calc_browser_refs.get("calc_path_input")
            if not _path_w:
                _append_system_message("Calc browser not available.")
                return True
            # Navigate the browser directly to the file — this triggers
            # the browser's own path handler which selects the file,
            # opens the preview pane, and also refreshes the options
            # dropdown (Recalc/Smart Recalc/etc.).
            target = _resolve_calc_path(filename)
            if not target.exists():
                _append_system_message(f"Not found: {filename}")
                return True
            _path_w.value = str(target)
            _append_system_message(f"Opened in browser: {target.name}")
            return True

        if cmd.startswith("/calc read "):
            filepath = text[len("/calc read"):].strip()
            target = _resolve_calc_path(filepath)
            if not target.is_file():
                _append_system_message(f"Not a file: {filepath}")
                return True
            size = target.stat().st_size
            limit = 8192 if target.suffix in (".out", ".log") else 32768
            try:
                content = target.read_text(encoding="utf-8", errors="replace")[:limit]
                if size > limit:
                    content += f"\n... [truncated, {size:,} bytes total]"
                _append_system_message(f"```\n{content}\n```")
            except Exception as exc:
                _append_system_message(f"Error reading: {exc}")
            return True

        if cmd.startswith("/calc tail "):
            filepath = text[len("/calc tail"):].strip()
            target = _resolve_calc_path(filepath)
            if not target.is_file():
                _append_system_message(f"Not a file: {filepath}")
                return True
            try:
                size = target.stat().st_size
                read_size = min(size, 8192)
                with open(target, "rb") as f:
                    f.seek(max(0, size - read_size))
                    tail = f.read().decode("utf-8", errors="replace")
                _append_system_message(f"Last {read_size:,} bytes of {target.name}:\n```\n{tail}\n```")
            except Exception as exc:
                _append_system_message(f"Error: {exc}")
            return True

        if cmd.startswith("/calc info "):
            folder = text[len("/calc info"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            try:
                from delfin.smart_recalc import has_ok_marker
            except ImportError:
                has_ok_marker = None
            lines = [f"Folder: {target.name}"]
            for f in sorted(target.iterdir()):
                if f.is_file():
                    sz = f.stat().st_size
                    status = ""
                    if f.suffix == ".out" and has_ok_marker:
                        ok = has_ok_marker(f)
                        status = " [OK]" if ok else " [INCOMPLETE/ERROR]"
                    lines.append(f"  {sz:>10,}  {f.name}{status}")
            control = target / "CONTROL.txt"
            if control.exists():
                lines.append(f"\n  CONTROL.txt present ({control.stat().st_size:,} bytes)")
            _append_system_message("\n".join(lines))
            return True

        if cmd.startswith("/calc tree"):
            folder = text[len("/calc tree"):].strip()
            target = _resolve_calc_path(folder) if folder else _resolve_calc_path("")
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            lines = [f"{target.name}/"]
            def _tree(d, prefix, depth):
                if depth > 2:
                    return
                items = sorted(d.iterdir(), key=lambda p: (p.is_file(), p.name))[:50]
                for i, p in enumerate(items):
                    connector = "\u2514\u2500 " if i == len(items) - 1 else "\u251c\u2500 "
                    if p.is_dir():
                        lines.append(f"{prefix}{connector}{p.name}/")
                        ext = "   " if i == len(items) - 1 else "\u2502  "
                        _tree(p, prefix + ext, depth + 1)
                    else:
                        lines.append(f"{prefix}{connector}{p.name}")
            _tree(target, "", 0)
            _append_system_message("\n".join(lines[:200]))
            return True

        if cmd.startswith("/calc search "):
            pattern = text[len("/calc search"):].strip()
            target = ctx.calc_dir
            matches = sorted(target.glob(pattern))[:50]
            if matches:
                lines = [str(m.relative_to(target)) for m in matches]
                _append_system_message(f"Found {len(matches)} matches:\n" + "\n".join(lines))
            else:
                _append_system_message(f"No matches for: {pattern}")
            return True

        # -- Analysis (all SAFE) ---------------------------------------------

        if cmd.startswith("/analyze energy "):
            folder = text[len("/analyze energy"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            try:
                from delfin.energies import find_gibbs_energy, find_ZPE, find_electronic_energy
            except ImportError:
                _append_system_message("Energy parsing not available.")
                return True
            lines = [f"Energies for {target.name}:"]
            for out_file in sorted(target.glob("*.out")):
                e = find_electronic_energy(str(out_file))
                g = find_gibbs_energy(str(out_file))
                z = find_ZPE(str(out_file))
                parts = []
                if e is not None:
                    parts.append(f"E={e:.6f}")
                if g is not None:
                    parts.append(f"G={g:.6f}")
                if z is not None:
                    parts.append(f"ZPE={z:.6f}")
                if parts:
                    lines.append(f"  {out_file.name}: {', '.join(parts)} Eh")
                else:
                    lines.append(f"  {out_file.name}: no energies found")
            if len(lines) == 1:
                lines.append("  (no .out files)")
            _append_system_message("\n".join(lines))
            return True

        # /analyze rank <type> [limit] [root]
        # Cross-folder energy ranking: scans every *.out in calc/,
        # archive/, and remote_archive/, extracts the requested energy
        # type, sorts ascending (lower = better for gibbs/electronic/zpe).
        # One tool call instead of ~50 cd/ls/tail rounds.
        if cmd.startswith("/analyze rank"):
            rest = text[len("/analyze rank"):].strip().split()
            if not rest:
                _append_system_message(
                    "Usage: /analyze rank <type> [limit] [root]\n"
                    "  type   : gibbs | electronic | spe | zpe\n"
                    "  limit  : how many top hits to show (default 10)\n"
                    "  root   : calc | archive | remote_archive | all "
                    "(default all)"
                )
                return True
            etype_raw = rest[0].lower()
            _ETYPE_ALIASES = {
                "gibbs": "gibbs", "g": "gibbs", "gfe": "gibbs",
                "electronic": "electronic", "e": "electronic",
                "spe": "electronic", "single": "electronic",
                "zpe": "zpe", "z": "zpe",
            }
            etype = _ETYPE_ALIASES.get(etype_raw)
            if etype is None:
                _append_system_message(
                    f"Unknown energy type: {etype_raw}. "
                    f"Choose from: gibbs, electronic, spe, zpe."
                )
                return True
            limit = 10
            root_filter = "all"
            for arg in rest[1:]:
                if arg.isdigit():
                    limit = max(1, min(int(arg), 100))
                elif arg.lower() in ("calc", "archive", "remote_archive",
                                      "remote", "all"):
                    root_filter = (
                        "remote_archive" if arg.lower() == "remote"
                        else arg.lower()
                    )
            try:
                from delfin.energies import (
                    find_gibbs_energy, find_ZPE, find_electronic_energy,
                )
            except ImportError:
                _append_system_message("Energy parsing not available.")
                return True
            _extractor = {
                "gibbs": find_gibbs_energy,
                "electronic": find_electronic_energy,
                "zpe": find_ZPE,
            }[etype]
            roots: list[tuple[str, Path]] = []
            if root_filter in ("calc", "all"):
                roots.append(("calc", ctx.calc_dir))
            if root_filter in ("archive", "all"):
                roots.append(("archive", ctx.archive_dir))
            if root_filter in ("remote_archive", "all"):
                _ra = ctx.runtime_settings.get("remote_archive_dir", "")
                if _ra:
                    roots.append(("remote_archive", Path(_ra)))
            import time as _time
            t0 = _time.monotonic()
            _BUDGET_S = 15.0   # hard cap to keep dashboard responsive
            results: list[tuple[float, str, str]] = []
            scanned = 0
            for root_name, root in roots:
                if not root.is_dir():
                    continue
                # Limit recursion depth and skip obvious cache dirs.
                for out_file in root.rglob("*.out"):
                    if _time.monotonic() - t0 > _BUDGET_S:
                        break
                    parts = out_file.parts
                    if any(p.startswith(".") or p == "__pycache__"
                           for p in parts):
                        continue
                    scanned += 1
                    val = _extractor(str(out_file))
                    if val is None:
                        continue
                    try:
                        rel = str(out_file.relative_to(root))
                    except ValueError:
                        rel = str(out_file)
                    results.append((float(val), root_name, rel))
                if _time.monotonic() - t0 > _BUDGET_S:
                    break
            results.sort(key=lambda r: r[0])
            elapsed = _time.monotonic() - t0
            header = (
                f"Energy ranking — {etype} (ascending, lowest first). "
                f"Scanned {scanned} .out files in "
                f"{', '.join(r[0] for r in roots)} "
                f"in {elapsed:.1f}s. Showing top {min(limit, len(results))} "
                f"of {len(results)} with extractable {etype}:"
            )
            if not results:
                _append_system_message(
                    f"{header}\n  (no {etype} values found)"
                )
                return True
            lines = [header]
            for rank, (val, root_name, rel) in enumerate(results[:limit], 1):
                lines.append(f"  {rank:2d}. {val:>16.6f} Eh  "
                              f"[{root_name}] {rel}")
            if elapsed > _BUDGET_S:
                lines.append(
                    f"  (time budget reached at {_BUDGET_S}s — "
                    f"more .out files were not scanned)"
                )
            _append_system_message("\n".join(lines))
            return True

        if cmd.startswith("/analyze convergence "):
            folder = text[len("/analyze convergence"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            lines = [f"Convergence check for {target.name}:"]
            for out_file in sorted(target.glob("*.out")):
                try:
                    size = out_file.stat().st_size
                    with open(out_file, "rb") as f:
                        f.seek(max(0, size - 20480))
                        tail = f.read().decode("utf-8", errors="replace")
                    if "ORCA TERMINATED NORMALLY" in tail:
                        lines.append(f"  {out_file.name}: CONVERGED OK")
                    elif "SCF NOT CONVERGED" in tail or "CONVERGENCE" in tail.upper() and "FAILURE" in tail.upper():
                        lines.append(f"  {out_file.name}: SCF NOT CONVERGED")
                    elif "ABORTING" in tail or "ERROR" in tail:
                        lines.append(f"  {out_file.name}: ERROR/ABORTED")
                    else:
                        lines.append(f"  {out_file.name}: INCOMPLETE (no termination marker)")
                except Exception as exc:
                    lines.append(f"  {out_file.name}: read error ({exc})")
            if len(lines) == 1:
                lines.append("  (no .out files)")
            _append_system_message("\n".join(lines))
            return True

        if cmd.startswith("/analyze errors "):
            folder = text[len("/analyze errors"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            _ERROR_PATTERNS = [
                "ABORTING THE RUN", "FATAL ERROR", "SCF NOT CONVERGED",
                "Insufficient memory", "ran out of disk space",
                "mpirun detected that one or more processes exited",
                "TRAH STEP ABORTING", "ORCA finished by error termination",
            ]
            lines = [f"Error scan for {target.name}:"]
            found_any = False
            for out_file in sorted(target.glob("*.out")):
                try:
                    size = out_file.stat().st_size
                    with open(out_file, "rb") as f:
                        f.seek(max(0, size - 20480))
                        tail = f.read().decode("utf-8", errors="replace")
                    errors = [p for p in _ERROR_PATTERNS if p in tail]
                    if errors:
                        found_any = True
                        lines.append(f"  {out_file.name}: {', '.join(errors)}")
                except Exception:
                    pass
            if not found_any:
                lines.append("  No errors found in .out files.")
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/analyze status":
            lines = ["Calculation status overview:"]
            try:
                from delfin.smart_recalc import has_ok_marker
            except ImportError:
                has_ok_marker = None
            # Get running jobs
            running_dirs = set()
            try:
                jobs = ctx.backend.list_jobs()
                for j in jobs:
                    if j.status.upper() in ("RUNNING", "R", "PENDING", "PD"):
                        running_dirs.add(j.job_dir)
            except Exception:
                pass
            ok_count, fail_count, running_count, empty_count = 0, 0, 0, 0
            for d in sorted(ctx.calc_dir.iterdir()):
                if not d.is_dir():
                    continue
                out_files = list(d.glob("*.out"))
                if str(d) in running_dirs or str(d.resolve()) in running_dirs:
                    lines.append(f"  RUNNING   {d.name}")
                    running_count += 1
                elif not out_files:
                    empty_count += 1
                else:
                    all_ok = all(has_ok_marker(f) for f in out_files) if has_ok_marker else False
                    if all_ok:
                        ok_count += 1
                    else:
                        fail_count += 1
                        lines.append(f"  FAILED    {d.name}")
            lines.insert(1, f"  Completed: {ok_count}  |  Failed: {fail_count}  |  Running: {running_count}  |  Empty: {empty_count}")
            _append_system_message("\n".join(lines))
            return True

        if (cmd.startswith("/analyze ")
                and not cmd.startswith("/analyze energy")
                and not cmd.startswith("/analyze rank")
                and not cmd.startswith("/analyze convergence")
                and not cmd.startswith("/analyze errors")
                and cmd != "/analyze status"):
            # /analyze <folder> — full analysis
            folder = text[len("/analyze"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            # Run all three: energy + convergence + errors
            for sub_cmd in [f"/analyze energy {folder}", f"/analyze convergence {folder}", f"/analyze errors {folder}"]:
                _handle_slash_command(sub_cmd)
            return True

        # -- Recalc (check = SAFE, submit = CONFIRMATION REQUIRED) -----------

        if cmd.startswith("/recalc check-all"):
            lines = ["Recalc check (all folders):"]
            try:
                from delfin.smart_recalc import has_ok_marker, fingerprint_unchanged
            except ImportError:
                _append_system_message("smart_recalc not available.")
                return True
            needs = []
            for d in sorted(ctx.calc_dir.iterdir()):
                if not d.is_dir():
                    continue
                for inp in d.glob("*.inp"):
                    out = inp.with_suffix(".out")
                    ok = has_ok_marker(out) if out.exists() else False
                    fp = fingerprint_unchanged(inp) if ok else False
                    if not ok or not fp:
                        needs.append(d.name)
                        reason = "no output" if not out.exists() else ("incomplete" if not ok else "input changed")
                        lines.append(f"  NEEDS RECALC  {d.name} ({reason})")
                        break
            if not needs:
                lines.append("  All calculations are up to date.")
            else:
                lines.append(f"\n  Total: {len(needs)} folders need recalc")
            _append_system_message("\n".join(lines))
            return True

        if cmd.startswith("/recalc check "):
            folder = text[len("/recalc check"):].strip()
            target = _resolve_calc_path(folder)
            if not target.is_dir():
                _append_system_message(f"Not a directory: {folder}")
                return True
            try:
                from delfin.smart_recalc import has_ok_marker, fingerprint_unchanged
            except ImportError:
                _append_system_message("smart_recalc not available.")
                return True
            lines = [f"Recalc check for {target.name}:"]
            for inp in sorted(target.glob("*.inp")):
                out = inp.with_suffix(".out")
                ok = has_ok_marker(out) if out.exists() else False
                fp = fingerprint_unchanged(inp) if ok else False
                needs = not ok or not fp
                lines.append(
                    f"  {inp.name}: ok={ok}, fingerprint_match={fp}, "
                    f"{'NEEDS RECALC' if needs else 'up to date'}"
                )
            if len(lines) == 1:
                lines.append("  (no .inp files)")
            _append_system_message("\n".join(lines))
            return True

        if cmd == "/recalc auto":
            try:
                from delfin.smart_recalc import has_ok_marker, fingerprint_unchanged
            except ImportError:
                _append_system_message("smart_recalc not available.")
                return True
            needs = []
            for d in sorted(ctx.calc_dir.iterdir()):
                if not d.is_dir():
                    continue
                for inp in d.glob("*.inp"):
                    out = inp.with_suffix(".out")
                    ok = has_ok_marker(out) if out.exists() else False
                    fp = fingerprint_unchanged(inp) if ok else False
                    if not ok or not fp:
                        needs.append(d)
                        break
            if not needs:
                _append_system_message("All calculations are up to date. Nothing to recalc.")
                return True
            names = ", ".join(d.name for d in needs[:10])
            if len(needs) > 10:
                names += f" ... (+{len(needs) - 10} more)"
            def _do_recalc_auto():
                submitted = 0
                for job_dir in needs:
                    control = job_dir / "CONTROL.txt"
                    if not control.exists():
                        _append_system_message(f"  Skip {job_dir.name}: no CONTROL.txt")
                        continue
                    try:
                        result = ctx.backend.submit_delfin(
                            job_dir=job_dir, job_name=job_dir.name,
                            mode="delfin-recalc-classic",
                        )
                        if result.returncode == 0:
                            submitted += 1
                        else:
                            _append_system_message(f"  Failed {job_dir.name}: {result.stderr[:200]}")
                    except Exception as exc:
                        _append_system_message(f"  Error {job_dir.name}: {exc}")
                _append_system_message(f"Submitted {submitted}/{len(needs)} recalc jobs.")
                refresh = ctx.job_status_refs.get("refresh_job_list")
                if refresh:
                    try:
                        refresh()
                    except Exception:
                        pass
            _confirm_or_exec(
                "recalc_auto",
                f"Submit recalc for {len(needs)} folders: {names}",
                _do_recalc_auto,
                cmd_for_zone="/recalc auto",
            )
            return True

        if cmd.startswith("/recalc ") and not cmd.startswith("/recalc check"):
            raw_arg = text[len("/recalc"):].strip()
            if raw_arg == "auto":
                pass  # handled above
            else:
                # Support multiple folders: /recalc folder1 folder2 folder3
                folder_names = raw_arg.split()
                targets: list[Path] = []
                for fname in folder_names:
                    t = _resolve_calc_path(fname)
                    if not t.is_dir():
                        _append_system_message(f"Not a directory: {fname}")
                        continue
                    control = t / "CONTROL.txt"
                    if not control.exists():
                        _append_system_message(f"No CONTROL.txt in {fname}")
                        continue
                    targets.append(t)
                if not targets:
                    return True

                def _do_recalc_batch(job_dirs=targets):
                    submitted = 0
                    for job_dir in job_dirs:
                        try:
                            result = ctx.backend.submit_delfin(
                                job_dir=job_dir, job_name=job_dir.name,
                                mode="delfin-recalc-classic",
                            )
                            if result.returncode == 0:
                                submitted += 1
                                _append_system_message(f"Recalc submitted: {job_dir.name}")
                            else:
                                _append_system_message(f"Recalc failed {job_dir.name}: {result.stderr[:200]}")
                        except Exception as exc:
                            _append_system_message(f"Recalc error {job_dir.name}: {exc}")
                    _append_system_message(f"Done: {submitted}/{len(job_dirs)} recalc jobs submitted.")
                    refresh = ctx.recalc_refs.get("refresh_recalc_folders")
                    if refresh:
                        try:
                            refresh()
                        except Exception:
                            pass
                    refresh2 = ctx.job_status_refs.get("refresh_job_list")
                    if refresh2:
                        try:
                            refresh2()
                        except Exception:
                            pass
                    ctx.select_tab("Job Status")

                names = ", ".join(t.name for t in targets)
                _confirm_or_exec(
                    "recalc_batch",
                    f"Submit recalc for {len(targets)} folder(s): {names}?",
                    _do_recalc_batch,
                    cmd_for_zone="/recalc batch",
                )
            return True

        # -- Cancel jobs (CONFIRMATION REQUIRED) -----------------------------

        if cmd.startswith("/cancel "):
            target = text[len("/cancel"):].strip()
            if target in ("all", "running", "pending"):
                # Filter by SLURM state. "all" = both pending+running,
                # "running" = R only, "pending" = PD only.
                if target == "running":
                    state_match = {"RUNNING", "R"}
                    label_human = "running"
                elif target == "pending":
                    state_match = {"PENDING", "PD"}
                    label_human = "pending"
                else:
                    state_match = {"RUNNING", "R", "PENDING", "PD"}
                    label_human = "active"
                try:
                    jobs = ctx.backend.list_jobs()
                    active = [
                        j for j in jobs
                        if j.status.upper() in state_match
                    ]
                except Exception:
                    active = []
                if not active:
                    # Helpful hint when the user has 0 jobs but expected
                    # to see some (often: they're looking at someone
                    # else's jobs in squeue, or $USER points elsewhere).
                    try:
                        import os as _os
                        my_user = _os.environ.get("USER", "?")
                    except Exception:
                        my_user = "?"
                    _append_system_message(
                        f"No {label_human} jobs to cancel for user "
                        f"`{my_user}`. (DELFIN only sees your own jobs; "
                        f"jobs from other users in the cluster queue are "
                        f"visible in `squeue` but not cancellable from "
                        f"this account.)"
                    )
                    return True
                names = ", ".join(f"{j.name}({j.job_id})" for j in active[:5])
                if len(active) > 5:
                    names += f" ... (+{len(active) - 5} more)"
                def _do_cancel_all():
                    cancelled = 0
                    for j in active:
                        try:
                            ok, msg = ctx.backend.cancel_job(j.job_id)
                            if ok:
                                cancelled += 1
                            else:
                                _append_system_message(f"  Failed to cancel {j.name}: {msg}")
                        except Exception as exc:
                            _append_system_message(f"  Error cancelling {j.name}: {exc}")
                    _append_system_message(f"Cancelled {cancelled}/{len(active)} jobs.")
                    refresh = ctx.job_status_refs.get("refresh_job_list")
                    if refresh:
                        try:
                            refresh()
                        except Exception:
                            pass
                _confirm_or_exec(
                    "cancel_all" if target == "all" else f"cancel_{target}",
                    f"Cancel {len(active)} {label_human} jobs: {names}",
                    _do_cancel_all,
                    cmd_for_zone=f"/cancel {target}",
                )
            else:
                job_id = target
                def _do_cancel(jid=job_id):
                    try:
                        ok, msg = ctx.backend.cancel_job(jid)
                        if ok:
                            _append_system_message(f"Job {jid} cancelled.")
                        else:
                            _append_system_message(f"Cancel failed: {msg}")
                        refresh = ctx.job_status_refs.get("refresh_job_list")
                        if refresh:
                            try:
                                refresh()
                            except Exception:
                                pass
                        ctx.select_tab("Job Status")
                    except Exception as exc:
                        _append_system_message(f"Cancel error: {exc}")
                _confirm_or_exec("cancel_single", f"Cancel job {job_id}?", _do_cancel, cmd_for_zone=f"/cancel {job_id}")
            return True

        # -- /batch commands ---------------------------------------------------

        if cmd.startswith("/batch "):
            sub = text[len("/batch"):].strip()

            if sub == "from-calc" or sub.startswith("from-calc"):
                # Use the Calculations Browser's Build Batch from XYZ system
                # 1) Select all folders, 2) prepare export, 3) fill Submit batch field
                calc_refs = ctx.calc_browser_refs or {}
                prepare_fn = calc_refs.get("xyz_batch_prepare_export")
                batch_state = calc_refs.get("xyz_batch_state")
                refresh_fn = calc_refs.get("xyz_batch_refresh")
                show_panel_fn = calc_refs.get("xyz_batch_show_panel")

                if not prepare_fn or not batch_state:
                    _append_system_message("Batch export function not available.")
                    return True

                # Optional filter: /batch from-calc Casagrande*
                parts = sub.split(None, 1)
                name_filter = parts[1].strip() if len(parts) > 1 else ""

                # Mark all folders (or filtered ones) in calc root
                calc_root = ctx.calc_dir
                all_folders = []
                for folder in sorted(calc_root.iterdir()):
                    if not folder.is_dir() or folder.name.startswith("."):
                        continue
                    if name_filter:
                        import fnmatch as _fnmatch
                        if not _fnmatch.fnmatch(folder.name, name_filter):
                            continue
                    all_folders.append(folder.name)

                if not all_folders:
                    if name_filter and (
                        name_filter.lower().endswith(".xyz")
                        or "." in name_filter
                    ):
                        # Common mistake: agent passed a file pattern as
                        # if the filter applied to file contents. Spell
                        # out the actual semantics.
                        _append_system_message(
                            f"No calculation folders matching folder-name "
                            f"glob '{name_filter}'.\n"
                            f"Tip: /batch from-calc (no filter) already "
                            f"walks every calc folder and collects each "
                            f"folder's initial.xyz / input.txt / coords.xyz "
                            f"automatically. The optional filter is for "
                            f"folder names (e.g. /batch from-calc Casagrande*), "
                            f"NOT for files inside the folders."
                        )
                    else:
                        hint = (
                            f" matching folder-name glob '{name_filter}'"
                            if name_filter else ""
                        )
                        _append_system_message(
                            f"No calculation folders{hint} found."
                        )
                    return True

                # Set marked paths in batch state
                batch_state["xyz_batch_marked_paths"] = list(all_folders)

                # Set batch filename
                batch_filename_widget = calc_refs.get("xyz_batch_filename")
                if batch_filename_widget:
                    batch_filename_widget.value = "batch"

                # Prepare export using the existing function
                export_data, skipped, error = prepare_fn()
                batch_state["xyz_batch_marked_paths"] = []  # clean up

                if error:
                    skip_info = f" (skipped: {', '.join(skipped[:5])})" if skipped else ""
                    _append_system_message(f"Batch build failed: {error}{skip_info}")
                    return True

                # Fill the Submit tab batch textarea
                batch_text = export_data["payload_text"]
                entry_count = export_data["entry_count"]
                batch_widget = ctx.submit_refs.get("smiles_batch_widget")
                if batch_widget:
                    batch_widget.value = batch_text
                    ctx.select_tab("Submit Job")
                    skip_info = ""
                    if skipped:
                        skip_info = f"\nSkipped {len(skipped)}: {', '.join(skipped[:5])}"
                    _append_system_message(
                        f"Batch filled with {entry_count} entries from calculations.{skip_info}\n"
                        f"Switched to Submit tab."
                    )
                else:
                    _append_system_message(
                        f"Built {entry_count} batch entries but batch widget not found."
                    )
                return True

            if sub == "clear":
                batch_widget = ctx.submit_refs.get("smiles_batch_widget")
                if batch_widget:
                    batch_widget.value = ""
                    _append_system_message("Batch field cleared.")
                else:
                    _append_system_message("Batch widget not found.")
                return True

            if sub == "show":
                batch_widget = ctx.submit_refs.get("smiles_batch_widget")
                if batch_widget:
                    val = batch_widget.value.strip()
                    if val:
                        line_count = len([l for l in val.splitlines() if l.strip()])
                        _append_system_message(f"Batch content ({line_count} lines):\n{val[:1000]}")
                    else:
                        _append_system_message("Batch field is empty.")
                else:
                    _append_system_message("Batch widget not found.")
                return True

            if sub.startswith("add "):
                entry = sub[4:].strip()
                if not entry or ";" not in entry:
                    _append_system_message(
                        "Format: /batch add Name;SMILES;key=value\n"
                        "    or: /batch add Name;charge=0;\\nXYZ\\ncoords\\n*"
                    )
                    return True
                batch_widget = ctx.submit_refs.get("smiles_batch_widget")
                if batch_widget:
                    current = batch_widget.value
                    if current and not current.endswith("\n"):
                        current += "\n"
                    batch_widget.value = current + entry
                    _append_system_message(f"Added to batch: {entry[:80]}")
                else:
                    _append_system_message("Batch widget not found.")
                return True

            _append_system_message(
                "Batch commands:\n"
                "  /batch from-calc          — Fill batch from calculations (input.txt, then *.xyz)\n"
                "  /batch from-calc *.xyz    — Fill batch from *.xyz files only\n"
                "  /batch from-calc input.txt — Fill batch from input.txt only\n"
                "  /batch add Name;SMILES;key=value  — Add one entry\n"
                "  /batch show  — Show current batch content\n"
                "  /batch clear  — Clear batch field"
            )
            return True

        return False

    # -- dashboard mode helpers --------------------------------------------

    def _collect_solo_domain_snapshot() -> dict:
        """Collect a compact dict of dashboard/calc state for solo mode.

        Pure data — handed to the module-level :func:`_format_solo_domain_state`
        for rendering.  Each lookup is wrapped so a missing widget never
        breaks the snapshot.
        """
        snap: dict = {"calc_dir": str(ctx.calc_dir)}

        # Selected file / folder in the calc browser
        try:
            cb = ctx.calc_browser_refs or {}
            cur_path = state.get("_agent_calc_path", "") or ""
            sel_widget = cb.get("calc_selected") or cb.get("calc_search")
            sel = ""
            if sel_widget is not None and getattr(sel_widget, "value", ""):
                sel = str(sel_widget.value)
            if cur_path:
                snap["selected"] = (
                    f"{cur_path}/{sel}".strip("/") if sel else cur_path
                )
        except Exception:
            pass

        # CONTROL keys parsed from the Submit tab textarea
        try:
            cw = ctx.submit_refs.get("control_widget")
            text = (cw.value if cw else "") or ""
            ctl: dict[str, str] = {}
            for raw_line in text.splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, _, v = line.partition("=")
                ctl[k.strip()] = v.strip()
            if ctl:
                snap["control"] = ctl
        except Exception:
            pass

        # ORCA Builder (dropdown values)
        try:
            refs = ctx.orca_builder_refs or {}
            builder: dict = {}
            for key, snap_key in (
                ("orca_method", "method"), ("orca_basis", "basis"),
                ("orca_charge", "charge"), ("orca_mult", "mult"),
                ("orca_pal", "pal"),
            ):
                w = refs.get(key)
                if w is not None and getattr(w, "value", "") not in ("", None):
                    builder[snap_key] = w.value
            if builder:
                snap["orca_builder"] = builder
        except Exception:
            pass

        # Job summary (compact "N RUNNING, M PENDING" string)
        try:
            jrefs = ctx.job_status_refs or {}
            counts = jrefs.get("status_counts")
            if isinstance(counts, dict) and counts:
                parts = [
                    f"{n} {st}" for st, n in sorted(counts.items()) if n
                ]
                if parts:
                    snap["job_summary"] = ", ".join(parts)
        except Exception:
            pass

        # Workspace inventory (small)
        try:
            files = sorted(
                p.name for p in ctx.agent_dir.iterdir() if p.is_file()
            )
            if files:
                snap["workspace_files"] = files[:20]
        except Exception:
            pass

        # Active tab title for cross-tab awareness
        try:
            tabs = ctx.tabs_widget
            if tabs is not None and ctx.tab_indices:
                idx = tabs.selected_index
                title = next(
                    (t for t, i in ctx.tab_indices.items() if i == idx),
                    "",
                )
                if title:
                    snap["active_tab"] = title
        except Exception:
            pass

        # Permission profile
        perm = state.get("_perm_profile", "")
        if perm:
            snap["perm_profile"] = perm

        return snap

    def _build_dashboard_context() -> str:
        """Build current dashboard state as context for the dashboard agent."""
        parts = []
        # CONTROL content
        cw = ctx.submit_refs.get("control_widget")
        if cw and cw.value.strip():
            parts.append(f"Current CONTROL (Submit tab):\n```\n{cw.value.strip()}\n```")
        # ORCA Builder settings
        refs = ctx.orca_builder_refs
        if refs:
            method = refs.get("orca_method")
            basis = refs.get("orca_basis")
            charge = refs.get("orca_charge")
            if method:
                parts.append(
                    f"ORCA Builder: method={method.value}, "
                    f"basis={basis.value if basis else '?'}, "
                    f"charge={charge.value if charge else '?'}"
                )
        # Job name
        jn = ctx.submit_refs.get("job_name_widget")
        if jn and jn.value.strip():
            parts.append(f"Job name: {jn.value.strip()}")
        # User identity (so agent knows who it's talking to)
        try:
            import os as _os
            _user = _os.environ.get("USER", "") or _os.getlogin()
            if _user:
                parts.append(f"Current user: {_user}")
        except Exception:
            pass
        # calc_dir + archive dirs + workspace + permissions
        parts.append(f"Calculations dir: {ctx.calc_dir}")
        parts.append(f"Archive dir: {ctx.archive_dir} (READ-ONLY: you can read/browse)")
        _remote_arch = ctx.runtime_settings.get("remote_archive_dir", "")
        if _remote_arch:
            parts.append(f"Remote archive dir: {_remote_arch} (READ-ONLY: you can read/browse)")
        parts.append(f"Agent workspace: {ctx.agent_dir}")
        perm = state.get("_perm_profile", "ask_all")
        _perm_desc = {
            "plan": "READ-ONLY everywhere",
            "ask_all": "ask for all changes",
            "repo_free": "repo free, calc asks",
            "all_free": "full access (archive + remote archive read-only)",
        }
        parts.append(f"Permissions: {perm} — {_perm_desc.get(perm, perm)}")
        # List workspace files if any
        try:
            ws_files = [p.name for p in ctx.agent_dir.iterdir() if p.is_file()]
            if ws_files:
                parts.append(f"Workspace files: {', '.join(ws_files[:20])}")
        except Exception:
            pass
        return "\n".join(parts)

    def _build_dashboard_session_boot() -> str:
        """Collect once-per-session domain snapshot and format it.

        Used on the FIRST user-send of a fresh session so the agent
        starts with the same domain awareness a returning developer
        would have: recent outcomes, active jobs, recent commits,
        active branch, current calc folder.
        """
        outcomes: list = []
        try:
            from delfin.agent.outcome_tracker import load_outcomes
            outcomes = load_outcomes(max_entries=5)
        except Exception:
            outcomes = []

        jobs: list = []
        try:
            backend = getattr(ctx, "backend", None)
            if backend is not None and hasattr(backend, "list_jobs"):
                jobs = list(backend.list_jobs() or [])
        except Exception:
            jobs = []

        commits: list[str] = []
        branch_name = ""
        try:
            import subprocess as _sp
            _r = _sp.run(
                ["git", "log", "-3", "--oneline"],
                cwd=str(ctx.repo_dir or Path.cwd()),
                capture_output=True, text=True, timeout=2,
            )
            if _r.returncode == 0:
                commits = [
                    line for line in _r.stdout.splitlines() if line.strip()
                ]
            _b = _sp.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                cwd=str(ctx.repo_dir or Path.cwd()),
                capture_output=True, text=True, timeout=2,
            )
            if _b.returncode == 0:
                branch_name = _b.stdout.strip()
        except Exception:
            pass

        calc_path = ""
        try:
            _w = ctx.calc_browser_refs.get("calc_path_input")
            if _w and _w.value:
                calc_path = str(_w.value)
        except Exception:
            calc_path = ""

        return _format_session_boot(
            outcomes=outcomes,
            jobs=jobs,
            commits=commits,
            calc_path=calc_path,
            branch=branch_name,
        )

    def _build_solo_session_boot() -> str:
        """One-shot domain primer for SOLO-mode first-send.

        Solo agents do code work, not workflow orchestration, so the
        primer skips active SLURM jobs and active calc-folder (those
        belong to dashboard). Keeps recent outcomes (last sessions'
        successes/failures inform the next task), active branch, and
        recent commits — exactly the orient-yourself info a returning
        developer would want.
        """
        outcomes: list = []
        try:
            from delfin.agent.outcome_tracker import load_outcomes
            outcomes = load_outcomes(max_entries=5)
        except Exception:
            outcomes = []

        commits: list[str] = []
        branch_name = ""
        try:
            import subprocess as _sp
            _r = _sp.run(
                ["git", "log", "-3", "--oneline"],
                cwd=str(ctx.repo_dir or Path.cwd()),
                capture_output=True, text=True, timeout=2,
            )
            if _r.returncode == 0:
                commits = [
                    line for line in _r.stdout.splitlines() if line.strip()
                ]
            _b = _sp.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                cwd=str(ctx.repo_dir or Path.cwd()),
                capture_output=True, text=True, timeout=2,
            )
            if _b.returncode == 0:
                branch_name = _b.stdout.strip()
        except Exception:
            pass

        return _format_session_boot(
            outcomes=outcomes,
            commits=commits,
            branch=branch_name,
            # jobs + calc_path intentionally omitted for solo
        )

    # -- command safety tiers (enforced at CODE level, not prompt) ----------

    _TIER3_EXACT = {
        "/submit", "/orca submit", "/recalc auto",
        "/cancel all", "/cancel running", "/cancel pending",
    }
    _TIER3_PREFIX = ("/recalc ", "/cancel ")
    _TIER3_SAFE_PREFIX = ("/recalc check",)  # these stay tier 0

    def _command_tier(cmd: str) -> int:
        """Classify a slash command: 0=read, 1=navigate, 2=configure, 3=mutate."""
        cl = cmd.lower().strip()
        # Tier 3: destructive / irreversible
        if cl in _TIER3_EXACT:
            return 3
        for pfx in _TIER3_PREFIX:
            if cl.startswith(pfx):
                if any(cl.startswith(sp) for sp in _TIER3_SAFE_PREFIX):
                    return 0
                return 3
        # Tier 2: session-only config changes
        if cl.startswith(("/control set", "/control key", "/orca set")):
            return 2
        # /ui: tier depends on property and target widget
        if cl.startswith("/ui "):
            _ui_sub = cl[4:].strip().split(None, 2)
            if cl == "/ui list":
                return 0
            if len(_ui_sub) >= 2 and _ui_sub[1] in ("show", "options"):
                return 0
            if len(_ui_sub) >= 2 and _ui_sub[1] == "click":
                # Clicking submit/recalc/transfer buttons = Tier 3
                _btn = _ui_sub[0] if _ui_sub else ""
                if any(k in _btn for k in ("submit", "recalc", "ssh", "archive", "transfer", "cancel", "override")):
                    return 3
                return 2  # other button clicks = Tier 2
            if len(_ui_sub) >= 2 and _ui_sub[1] in ("value", "replace"):
                return 2
            return 1  # style, text, visible, width, height, disabled
        # Tier 1: navigation
        if cl.startswith(("/tab ", "/jobs", "/mode ", "/model ", "/provider ")):
            return 1
        # Tier 0: read-only
        return 0

    # -- permission profiles -------------------------------------------------
    # Each profile maps zone → (max_tier, confirm_tier3).
    #   max_tier:       highest tier allowed (-1=blocked, 0=read, 3=all)
    #   confirm_tier3:  True = tier-3 commands need Approve/Deny dialog
    #
    # remote_archive is ALWAYS read-only — no profile can override this.
    _PERM_PROFILES: dict[str, dict[str, tuple[int, bool]]] = {
        "plan": {
            # Browse + dashboard control — agent can navigate, style UI, fill fields
            # but cannot submit jobs, cancel, recalc, or write files
            "workspace":      (2, True),
            "calc":           (2, True),
            "repo":           (2, True),
            "archive":        (1, True),
            "remote_archive": (0, True),
            "unknown":        (-1, True),
        },
        "ask_all": {
            # Ask for everything that changes
            "workspace":      (3, False),  # agent sandbox — free
            "calc":           (3, True),   # mutate with confirmation
            "repo":           (3, True),   # mutate with confirmation
            "archive":        (0, True),   # read-only ALWAYS
            "remote_archive": (0, True),   # read-only ALWAYS
            "unknown":        (-1, True),
        },
        "repo_free": {
            # Repo free, calc asks
            "workspace":      (3, False),
            "calc":           (3, True),   # still needs confirmation
            "repo":           (3, False),  # auto-approved
            "archive":        (0, True),   # read-only ALWAYS
            "remote_archive": (0, True),   # read-only ALWAYS
            "unknown":        (-1, True),
        },
        "all_free": {
            # Everything works except archive & remote archive
            "workspace":      (3, False),
            "calc":           (3, False),  # auto-approved
            "repo":           (3, False),  # auto-approved
            "archive":        (0, True),   # read-only ALWAYS
            "remote_archive": (0, True),   # read-only ALWAYS
            "unknown":        (-1, True),
        },
    }

    # Map DELFIN profile → Claude CLI permission_mode.
    # Safety: never use bypassPermissions for ANY mode — all modes cap at
    # 'acceptEdits', which auto-approves file edits but still asks for Bash.
    # The DELFIN zone system is the primary safety layer; CLI permissions
    # are defense-in-depth.
    _PROFILE_TO_CLI_PERM: dict[str, str] = {
        "plan":      "plan",           # read-only
        "ask_all":   "default",        # CLI asks before every write/bash
        "repo_free": "acceptEdits",    # CLI auto-approves edits, asks for bash
        "all_free":  "auto",           # user chose "all free" → skip permission prompts
    }

    def _active_perms() -> dict[str, tuple[int, bool]]:
        """Return the zone permissions for the active profile."""
        profile = state.get("_perm_profile", "ask_all")
        return _PERM_PROFILES.get(profile, _PERM_PROFILES["ask_all"])

    def _active_cli_perm() -> str:
        """Return the Claude CLI permission_mode for the active profile.

        Dashboard mode always uses 'default' (asks before write tools).
        All other modes map through _PROFILE_TO_CLI_PERM which caps at
        'acceptEdits' — never bypassPermissions.
        """
        cur_mode = mode_dropdown.value
        if cur_mode == "dashboard":
            return "default"
        profile = state.get("_perm_profile", "ask_all")
        return _PROFILE_TO_CLI_PERM.get(profile, "default")

    # -- path zone classification -------------------------------------------

    def _path_zone(cmd: str) -> str:
        """Classify a /calc command's target into a permission zone.

        Returns one of: 'workspace', 'calc', 'archive', 'remote_archive',
        'repo', 'unknown'.  Commands without a path argument return 'calc'
        (the default browsing root).
        """
        parts = cmd.split(None, 2)
        if len(parts) < 3:
            return "calc"  # no path arg → default calc_dir

        raw = parts[2]

        # Resolve to absolute path for reliable comparison
        agent_rel = state.get("_agent_calc_path", "")
        base = ctx.calc_dir / agent_rel if agent_rel else ctx.calc_dir
        target = (base / raw).resolve()

        # Check zones in specificity order (most specific first)
        try:
            target.relative_to(ctx.agent_dir.resolve())
            return "workspace"
        except (ValueError, RuntimeError):
            pass

        # Remote archive (check before regular archive — may be a subdir)
        _remote = ctx.runtime_settings.get("remote_archive_dir", "")
        if _remote:
            try:
                target.relative_to(Path(_remote).resolve())
                return "remote_archive"
            except (ValueError, RuntimeError):
                pass

        # Archive
        try:
            target.relative_to(ctx.archive_dir.resolve())
            return "archive"
        except (ValueError, RuntimeError):
            pass

        # Fallback: check keywords in the raw path for commands like
        # "/calc ls archive/" that may not resolve cleanly
        p = raw.lower()
        if "remote" in p and "archive" in p:
            return "remote_archive"
        if "archive" in p:
            archive_str = str(ctx.archive_dir).lower()
            if "archive" in p or p.startswith(archive_str):
                return "archive"

        # Calculations dir
        try:
            target.relative_to(ctx.calc_dir.resolve())
            return "calc"
        except (ValueError, RuntimeError):
            pass

        # Repo dir (for code agents)
        if ctx.repo_dir:
            try:
                target.relative_to(Path(ctx.repo_dir).resolve())
                return "repo"
            except (ValueError, RuntimeError):
                pass

        return "unknown"

    def _zone_blocks(cmd: str, tier: int) -> str | None:
        """Return a block message if *tier* exceeds zone permissions, else None.

        Messages are intentionally generic to avoid leaking zone
        classification to the agent (security: no directory enumeration).
        """
        zone = _path_zone(cmd)
        perms = _active_perms()
        max_tier, _ = perms.get(zone, (-1, True))
        if tier > max_tier:
            profile = state.get("_perm_profile", "ask_all")
            if profile == "plan":
                return "⛔ Blocked: Current permission profile is read-only."
            # Generic message — don't reveal zone name or directory structure
            return "⛔ Blocked: This action is not permitted. Ask the user for approval."
        return None

    def _zone_needs_confirm(cmd: str) -> bool:
        """True if a tier-3 command in this zone needs the confirmation dialog."""
        zone = _path_zone(cmd)
        perms = _active_perms()
        _, confirm = perms.get(zone, (-1, True))
        return confirm

    # Bulk operation intent detection.
    # Requires BOTH an action keyword AND a scope keyword in the user's
    # recent messages (last 3, to allow "kill alle" → "yes" follow-up).
    # This prevents the agent from self-triggering bulk ops by mentioning
    # just one keyword in its own output, while letting the user confirm
    # with a short "ja" / "yes" / "ok" without re-typing the full intent.
    _BULK_ACTION_KW = (
        "recalc", "neuberechn", "submit", "absend", "abschick",
        "cancel", "abbrech", "stopp", "stop",
        # killing / terminating phrasings (added 2026-05-11 after the
        # "kill alle running jobs" user-request was blocked because
        # "kill" was missing from the list).
        "kill", "killen", "terminate", "terminier", "beenden",
        "halt", "shutdown", "zerstör",
    )
    _BULK_SCOPE_KW = (
        "all", "alle", "auto", "alles", "every", "sämtliche", "bulk",
        "komplett", "gesamt", "running", "laufend", "aktiv",
        "pending", "wartend", "queued", "queue",
    )
    # How many of the most recent user messages to scan for intent.
    # Set to 3 so "kill alle running jobs" → agent asks "soll ich?" →
    # user says "ja" still finds the original action+scope intent.
    _BULK_INTENT_LOOKBACK = 3

    def _dashboard_auto_exec(agent_text: str, force_no_confirm: bool = False):
        """Scan agent output for ACTION: /command lines and execute them.

        Safety enforcement (code-level, not prompt-level):
        - Zone-based permissions: workspace=free, calc=confirm,
          archive/remote_archive=read-only (writes blocked, reads allowed), unknown=blocked
        - Tier 3: max 1 per response, bulk ops need explicit user intent
        - Workspace zone: tier 3 skips confirmation gate

        D1 — when the agent ends its response with a confirmation question
        and proposes ACTIONs, surface Approve/Deny buttons instead of
        running auto-exec.  ``force_no_confirm=True`` bypasses this gate
        (used by the Approve handler).
        """
        # D1 short-circuit: agent is asking for confirmation → show
        # Approve/Deny buttons and skip auto-exec until the user clicks.
        if (
            not force_no_confirm
            and _should_show_action_confirmation(agent_text)
        ):
            commands_preview = _extract_action_commands(agent_text)
            if commands_preview:
                _show_action_confirmation(agent_text, commands_preview)
                return []  # nothing executed — wait for the user

        import re as _re

        # ACTION: /done is a sentinel — the agent uses it to signal that
        # all requested actions are now complete and the auto-continue
        # loop should NOT re-prompt the model for a wrap-up commentary.
        # Saves the 30-120 s "(commands executed)" turn that otherwise
        # costs $0.02-0.05 for zero user value. The marker is added to
        # the END of the results list AFTER real ACTIONs have executed,
        # so a response with both real commands + /done both runs the
        # commands AND skips the wrap-up turn. A response with /done
        # alone returns only ["__DONE__"], which the continuation loop
        # treats as "agent had no actions to execute".
        lines = agent_text.split("\n")
        _done_seen = any(_re.match(r"^ACTION:\s*/done\b", ln) for ln in lines)

        commands: list[str] = []
        i = 0
        while i < len(lines):
            m = _re.match(r"^ACTION:\s*(/\S+.*)$", lines[i])
            if m:
                cmd = m.group(1)
                # /done is the wrap-up sentinel — never executed as a real
                # command, just signals to the outer loop. Skip the
                # dispatch path entirely.
                if cmd.strip().lower().startswith("/done"):
                    i += 1
                    continue
                # For /control set: collect continuation lines as content
                if cmd.lower().startswith("/control set "):
                    content_lines = [cmd[len("/control set "):]]
                    i += 1
                    while i < len(lines):
                        if _re.match(r"^ACTION:\s*/", lines[i]):
                            break
                        if lines[i].strip() in ("```", ""):
                            i += 1
                            continue
                        content_lines.append(lines[i])
                        i += 1
                    commands.append("/control set " + "\n".join(content_lines))
                else:
                    # Convert escaped \n to real newlines (agents often
                    # output literal \n for multiline values like coords)
                    commands.append(cmd.replace("\\n", "\n"))
                    i += 1
            else:
                i += 1

        results: list[str] = []
        mutate_count = 0

        for cmd_line in commands[:10]:  # safety limit
            tier = _command_tier(cmd_line)
            short = cmd_line[:80] + ("..." if len(cmd_line) > 80 else "")

            # --- Zone-based permission check ---
            block_msg = _zone_blocks(cmd_line, tier)
            if block_msg:
                _append_system_message(block_msg)
                results.append(f"BLOCKED: {block_msg}")
                continue

            zone = _path_zone(cmd_line)

            # --- Tier 3: max 1 per response ---
            if tier == 3:
                mutate_count += 1
                if mutate_count > 1:
                    _append_system_message(
                        "\u26d4 Blocked: Only one destructive action per response. "
                        "Ask the user for the next step."
                    )
                    results.append("BLOCKED: max 1 destructive action per response")
                    continue

                # Bulk ops need explicit user intent (action + scope) in
                # the last N user messages \u2014 not just the immediate one,
                # so "kill alle running jobs" \u2192 "ja" follow-up still
                # finds the action+scope keywords in recent history.
                cl = cmd_line.lower().strip()
                if cl in (
                    "/recalc auto",
                    "/cancel all", "/cancel running", "/cancel pending",
                ):
                    recent_user_msgs: list[str] = []
                    for m in reversed(state.get("chat_messages", []) or []):
                        if m.get("role") == "user":
                            recent_user_msgs.append(
                                str(m.get("content", "")).lower()
                            )
                            if len(recent_user_msgs) >= _BULK_INTENT_LOOKBACK:
                                break
                    # Backward-compat fallback: include _last_user_message
                    # in case chat_messages wasn't populated yet.
                    fallback = state.get("_last_user_message", "").lower()
                    if fallback and fallback not in recent_user_msgs:
                        recent_user_msgs.append(fallback)
                    haystack = " ".join(recent_user_msgs)
                    has_action = any(kw in haystack for kw in _BULK_ACTION_KW)
                    has_scope = any(kw in haystack for kw in _BULK_SCOPE_KW)
                    if not (has_action and has_scope):
                        _append_system_message(
                            "\u26d4 Blocked: Bulk operation requires explicit "
                            "user request (e.g. 'kill alle', 'cancel all', "
                            "'recalc alle'). Report findings and ask the user."
                        )
                        results.append("BLOCKED: bulk op without user intent")
                        continue

            # --- Execute ---
            _append_system_message(f"\u25b6 Executing: {short}")
            n_before = len(state["chat_messages"])
            try:
                handled = _handle_slash_command(cmd_line)
                if not handled:
                    # Friendly hint: if the agent tried `/calc`, `/jobs`,
                    # `/archive` etc. as a plain slash-command, suggest the
                    # proper `/tab <name>` form. Saves a round-trip when
                    # the agent meant "open this tab".
                    _TAB_NAMES = {
                        "submit", "recalc", "jobs", "orca", "calc",
                        "archive", "literature", "agent", "settings",
                        "calculations", "calculation", "calcs",
                        "berechnung", "berechnungen", "literatur",
                        "archiv", "einstellungen", "job",
                    }
                    _bare = cmd_line.strip().lstrip("/").split()[0].lower()
                    if _bare in _TAB_NAMES:
                        _append_system_message(
                            f"Unknown command: {short}. "
                            f"Did you mean: /tab {_bare} ?"
                        )
                    else:
                        _append_system_message(f"Unknown command: {short}")
            except Exception as exc:
                _append_system_message(f"Error executing {short}: {exc}")
            # Collect new system messages as feedback
            for msg in state["chat_messages"][n_before:]:
                if msg["role"] == "system":
                    results.append(msg["content"])
        # Append the done sentinel AFTER real-command results. Callers
        # check ``"__DONE__" in results`` for the early-break signal;
        # ``results == ["__DONE__"]`` specifically means /done was the
        # ONLY thing the agent emitted (no real commands fired).
        if _done_seen:
            results.append("__DONE__")
        return results

    # -- feature implementations -------------------------------------------

    def _export_chat():
        """Export the chat history as a Markdown file."""
        if not state["chat_messages"]:
            _append_system_message("Nothing to export.")
            return
        lines = [f"# DELFIN Agent Chat Export\n"]
        engine = state["engine"]
        if engine:
            s = engine.get_status()
            lines.append(f"**Mode:** {s['mode']} | **Model:** {model_dropdown.value} | **Tokens:** {s['input_tokens']:,} in / {s['output_tokens']:,} out\n")
        lines.append("---\n")
        for msg in state["chat_messages"]:
            role = msg["role"]
            content = msg["content"]
            label = msg.get("role_label", "")
            if role == "user":
                lines.append(f"### User\n\n{content}\n")
            elif role == "assistant":
                rl = label or "Agent"
                lines.append(f"### {rl}\n\n{content}\n")
            elif role == "thinking":
                lines.append(f"<details><summary>Thinking</summary>\n\n{content}\n\n</details>\n")
            elif role == "gate":
                gate_title = msg.get("gate_title", "Gate")
                lines.append(f"> [GATE] {gate_title}\n>\n> {content}\n")
            elif role == "system":
                lines.append(f"> {content}\n")
        md_text = "\n".join(lines)
        # Write to file
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        export_dir = ctx.agent_dir / "exports"
        export_dir.mkdir(parents=True, exist_ok=True)
        export_path = export_dir / f"chat_{ts}.md"
        export_path.write_text(md_text, encoding="utf-8")
        _append_system_message(f"Chat exported to: {export_path}")

    def _toggle_search():
        """Toggle the search bar visibility."""
        visible = search_input.layout.display != "none"
        if visible:
            search_input.layout.display = "none"
            search_close_btn.layout.display = "none"
            search_count_html.value = ""
            search_input.value = ""
            # Restore original chat HTML (remove highlights)
            _refresh_chat_html()
        else:
            search_input.layout.display = "inline-flex"
            search_close_btn.layout.display = "inline-flex"
            # Focus will be handled by keyboard shortcut JS

    def _do_search(query):
        """Search chat messages and highlight matches."""
        if not query:
            _refresh_chat_html()
            search_count_html.value = ""
            return
        count = 0
        for msg in state["chat_messages"]:
            if query.lower() in msg["content"].lower():
                count += 1
        search_count_html.value = (
            f'<span style="font-size:11px;color:#6b7280;margin-left:6px;">'
            f'{count} match{"es" if count != 1 else ""}</span>'
        )
        # Rebuild chat HTML with highlighted matches
        if not state["chat_messages"]:
            return
        parts = ['<div class="delfin-agent-chat">']
        q_lower = query.lower()
        q_escaped = _html.escape(query)
        for msg in state["chat_messages"]:
            role = msg["role"]
            content = _md_to_html(msg["content"])
            # Highlight matches (case-insensitive)
            if q_lower in msg["content"].lower():
                # Highlight in the HTML output (rough but effective)
                content = re.sub(
                    re.escape(q_escaped),
                    f'<span class="delfin-search-hl">{q_escaped}</span>',
                    content,
                    flags=re.IGNORECASE,
                )
            if role == "user":
                parts.append(f'<div class="delfin-chat-msg delfin-chat-user">{content}</div>')
            elif role == "assistant":
                label = msg.get("role_label", "Agent")
                parts.append(
                    f'<div class="delfin-chat-msg delfin-chat-agent">'
                    f'<div class="delfin-chat-role">{_html.escape(label)}</div>'
                    f'{content}</div>'
                )
            elif role == "thinking":
                raw = msg["content"]
                preview = raw[:80].replace("\n", " ").strip()
                if len(raw) > 80:
                    preview += "..."
                escaped_full = _html.escape(raw)
                parts.append(
                    f'<details class="delfin-chat-msg delfin-chat-thinking">'
                    f'<summary>\U0001f9e0 {_html.escape(preview)}</summary>'
                    f'<div class="thinking-content">{escaped_full}</div></details>'
                )
            elif role == "gate":
                gate_type = msg.get("gate_type", "")
                gate_title = msg.get("gate_title", _gate_label(gate_type))
                gate_hint = msg.get("gate_hint", "")
                gate_role = msg.get("role_label", "")
                card_cls = f"delfin-chat-msg delfin-gate-card delfin-gate-{_html.escape(gate_type or 'review')}"
                badge = _gate_label(gate_type)
                hint_html = (
                    f'<div class="gate-hint">{_html.escape(gate_hint)}</div>'
                    if gate_hint else ""
                )
                role_html = (
                    f'<span class="gate-role">{_html.escape(gate_role)}</span>'
                    if gate_role else ""
                )
                parts.append(
                    f'<div class="{card_cls}">'
                    f'<div class="gate-header">'
                    f'<span class="delfin-gate-badge">{_html.escape(badge)}</span>'
                    f'<span class="gate-title">{_html.escape(gate_title)}</span>'
                    f'{role_html}'
                    f'</div>'
                    f'<div class="gate-body">{content}</div>'
                    f'{hint_html}'
                    f'</div>'
                )
            elif role == "system":
                raw = msg["content"]
                css_class = "delfin-chat-msg delfin-chat-handoff" if (
                    raw.startswith("---") or raw.startswith("Session restored")
                ) else "delfin-chat-msg delfin-chat-system"
                parts.append(f'<div class="{css_class}">{content}</div>')
        parts.append(
            '<img src="" onerror="'
            "var c=this.closest('.delfin-agent-chat');"
            "if(c)c.scrollTop=c.scrollHeight;"
            "this.remove();"
            '" style="display:none">'
        )
        parts.append("</div>")
        chat_html.value = "\n".join(parts)

    def _retry_last():
        """Remove the last assistant response and re-send the last user message."""
        if state["streaming"]:
            _append_system_message("Cannot retry while streaming.")
            return
        engine = state["engine"]
        if not engine or not engine.messages:
            _append_system_message("Nothing to retry.")
            return
        # Find the last user message
        last_user_text = ""
        # Remove trailing assistant + thinking messages from chat
        while state["chat_messages"] and state["chat_messages"][-1]["role"] in ("assistant", "thinking", "system", "gate", "approval"):
            state["chat_messages"].pop()
        _set_active_gate()
        if state["chat_messages"] and state["chat_messages"][-1]["role"] == "user":
            last_user_text = state["chat_messages"][-1]["content"]
            state["chat_messages"].pop()  # will be re-added by _on_send flow
        # Remove from engine messages too
        while engine.messages and engine.messages[-1]["role"] == "assistant":
            engine.messages.pop()
        if engine.messages and engine.messages[-1]["role"] == "user":
            engine.messages.pop()
        _refresh_chat_html()
        if last_user_text:
            # B3 — flag the next outcome write as a retry attempt so it
            # bumps the previous outcome's retries counter instead of
            # appending a fresh entry. Counter persists across the send
            # cycle and is consumed (cleared) by _record_turn_outcome.
            _prev_retries = int(state.get("_retry_pending_retries", 0))
            state["_retry_pending_retries"] = _prev_retries + 1
            # Re-send via the send flow
            input_textarea.value = last_user_text
            _on_send(None)
        else:
            _append_system_message("Could not find a user message to retry.")

    def _check_auto_compact():
        """Check if context is getting large and auto-compact or warn."""
        engine = state["engine"]
        if not engine:
            return
        total_tokens = engine.token_usage.get("input", 0)
        n_msgs = len(engine.messages)
        # Claude CLI auto-compacts internally, so we only do a silent
        # fallback compact on our engine messages if they get very large.
        if n_msgs > 30:
            old_count = n_msgs
            old_msgs = engine.messages[:-6]
            summary_parts = []
            for m in old_msgs[-8:]:
                content = m["content"][:200]
                summary_parts.append(f"[{m['role']}]: {content}...")
            summary = "\n".join(summary_parts)
            engine.messages = [
                {"role": "user", "content":
                 f"[Auto-compacted {old_count - 6} earlier messages:\n"
                 f"{summary}\n... End]"},
                {"role": "assistant", "content": "Understood, continuing with context."},
            ] + engine.messages[-6:]

    def _format_tool_description(raw):
        """Parse a raw permission denial string into a readable description."""
        import ast as _ast
        try:
            d = _ast.literal_eval(raw) if raw.strip().startswith("{") else {}
        except Exception:
            d = {}
        if not d:
            return raw[:200]
        tool = d.get("tool_name", "Unknown")
        inp = d.get("tool_input", {})
        if not isinstance(inp, dict):
            return f"{tool}"
        # Build a human-readable summary per tool type
        if tool == "Edit":
            fp = inp.get("file_path", "")
            old = inp.get("old_string", "")
            new = inp.get("new_string", "")
            fname = fp.rsplit("/", 1)[-1] if "/" in fp else fp
            old_short = old[:60].replace("\n", " ")
            new_short = new[:60].replace("\n", " ")
            return f"Edit {fname}: \"{old_short}\" → \"{new_short}\""
        if tool == "Write":
            fp = inp.get("file_path", "")
            fname = fp.rsplit("/", 1)[-1] if "/" in fp else fp
            return f"Write file: {fname}"
        if tool == "Bash":
            cmd = inp.get("command", "")
            return f"Run: {cmd[:120]}"
        if tool in ("Read", "Glob", "Grep"):
            target = inp.get("file_path", "") or inp.get("pattern", "") or inp.get("path", "")
            return f"{tool}: {target[:120]}"
        # Generic: show tool + first key
        parts = [f"{tool}:"]
        for k, v in list(inp.items())[:2]:
            sv = str(v)[:80]
            parts.append(f" {k}={sv}")
        return "".join(parts)

    def _request_confirmation(action_id, description, callback):
        """Show confirmation prompt for destructive dashboard operations.

        Stores the pending action and shows Approve/Deny buttons.
        When user approves, executes the callback.
        """
        state["_pending_dashboard_action"] = {
            "action_id": action_id,
            "description": description,
            "callback": callback,
        }
        _set_active_gate("approval", "", "Dashboard approval required", description)
        _append_chat_message(
            "approval",
            f"\u26a0\ufe0f {description}",
        )
        _update_status()
        approval_info_html.value = ""
        approve_btn.layout.display = "inline-flex"
        deny_btn.layout.display = "inline-flex"
        approval_row.layout.display = "flex"

    def _confirm_or_exec(action_id, description, callback, cmd_for_zone=""):
        """Either request confirmation or execute directly, based on permission profile.

        If the active profile says this zone needs confirmation for tier-3,
        show the Approve/Deny dialog.  Otherwise execute immediately.
        """
        if cmd_for_zone and not _zone_needs_confirm(cmd_for_zone):
            # Profile allows auto-execution in this zone
            callback()
            return
        # Needs confirmation
        _request_confirmation(action_id, description, callback)

    def _compute_approval_diff(detail: str) -> str:
        """Render a unified-diff preview when the pending tool is Edit /
        Write / multi_edit. Returns "" when the tool isn't file-mutating
        or the diff can't be computed (binary file, missing path, etc.).

        Reads the actual file from disk; pure-stdlib (difflib), capped to
        ~120 lines so the chat doesn't drown in a huge rewrite preview.
        """
        import ast as _ast
        import difflib as _difflib
        try:
            d = _ast.literal_eval(detail) if detail and detail.lstrip().startswith("{") else {}
        except Exception:
            return ""
        tool = (d.get("tool_name") or "").strip()
        inp = d.get("tool_input") or {}
        if not isinstance(inp, dict):
            return ""
        file_path = (inp.get("file_path") or inp.get("path") or "").strip()
        if not file_path or tool not in ("Edit", "Write", "multi_edit"):
            return ""
        try:
            target = Path(file_path).expanduser()
            current = target.read_text(encoding="utf-8") if target.is_file() else ""
        except OSError:
            return ""

        if tool == "Write":
            proposed = inp.get("content", "")
        elif tool == "Edit":
            old = inp.get("old_string", "")
            new = inp.get("new_string", "")
            if not old:
                return ""
            replace_all = bool(inp.get("replace_all", False))
            if old not in current:
                return f"(diff preview unavailable — old_string not found in {target.name})"
            proposed = (current.replace(old, new) if replace_all
                        else current.replace(old, new, 1))
        elif tool == "multi_edit":
            edits = inp.get("edits") or []
            proposed = current
            for e in edits:
                old = e.get("old_string", "")
                new = e.get("new_string", "")
                if old and old in proposed:
                    if e.get("replace_all"):
                        proposed = proposed.replace(old, new)
                    else:
                        proposed = proposed.replace(old, new, 1)
        else:
            return ""

        if proposed == current:
            return f"(no change to {target.name})"
        diff_lines = list(_difflib.unified_diff(
            current.splitlines(keepends=False),
            proposed.splitlines(keepends=False),
            fromfile=str(target),
            tofile=f"{target} (proposed)",
            lineterm="",
            n=3,
        ))
        if not diff_lines:
            return ""
        if len(diff_lines) > 120:
            diff_lines = diff_lines[:60] + [
                f"... ({len(diff_lines) - 120} lines truncated) ...",
            ] + diff_lines[-60:]
        return "```diff\n" + "\n".join(diff_lines) + "\n```"

    def _show_approval_prompt(tool_name, detail):
        """Show approval request inline in chat + approval buttons."""
        state["_pending_approval"] = {"tool": tool_name, "detail": detail}
        readable = _format_tool_description(detail)
        _set_active_gate("approval", "", "Tool approval required", readable)
        # Compute + surface a diff preview for file-mutating tools so the
        # user can see the actual change instead of just the tool name.
        try:
            diff_block = _compute_approval_diff(detail)
        except Exception:
            diff_block = ""
        if diff_block:
            _append_chat_message(
                "system", f"**Proposed change:**\n\n{diff_block}"
            )
        # Show in chat as a special approval message
        _append_chat_message(
            "approval",
            readable,
        )
        _update_status()
        # Show only the buttons between chat and input (description is in chat)
        approval_info_html.value = ""
        approve_btn.layout.display = "inline-flex"
        deny_btn.layout.display = "inline-flex"
        approval_row.layout.display = "flex"

    def _hide_approval():
        """Hide approval UI."""
        approve_btn.layout.display = "none"
        deny_btn.layout.display = "none"
        approval_row.layout.display = "none"
        approval_info_html.value = ""
        state.pop("_pending_approval", None)
        if (state.get("_active_gate") or {}).get("type") == "approval":
            _set_active_gate()
            _update_status()

    def _on_approve(button):
        """User approves a blocked operation.

        For dashboard actions: execute the stored callback.
        For Bash commands: execute directly via subprocess (reliable).
        For file ops: upgrade permission mode and let the agent retry.
        """
        import ast as _ast
        import subprocess as _sp

        # Dashboard action confirmations (from /recalc, /cancel, /submit)
        dashboard_action = state.get("_pending_dashboard_action")
        if dashboard_action:
            _hide_approval()
            state["_pending_dashboard_action"] = None
            desc = dashboard_action["description"]
            _append_system_message(f"\u2705 Approved: {desc}")
            try:
                dashboard_action["callback"]()
            except Exception as exc:
                _append_system_message(f"\u274c Error: {exc}")
            return

        pending = state.get("_pending_approval")
        if not pending:
            _hide_approval()
            return
        raw = pending.get("tool", "")
        readable = _format_tool_description(raw)
        _hide_approval()

        # Parse the blocked tool call
        try:
            d = _ast.literal_eval(raw) if raw.strip().startswith("{") else {}
        except Exception:
            d = {}
        tool = d.get("tool_name", "")
        inp = d.get("tool_input", {}) if isinstance(d.get("tool_input"), dict) else {}

        # --- Bash commands: route through delfin.agent.sandbox -----------
        # Defense in depth: allow-list + bwrap/firejail sandbox + audit log.
        # See delfin/agent/sandbox.py for the layer details.  The full,
        # untruncated command is shown to the user above this branch already
        # (Layer 4 \u2014 approval); we re-print it here so the audit trail in
        # the chat matches the audit JSONL.
        if tool == "Bash" and inp.get("command"):
            cmd = inp["command"]
            from delfin.agent import sandbox as _sandbox
            cfg = _sandbox.detect_config()
            _append_system_message(
                f"\u2705 Approved & executing (mode={cfg.mode}, "
                f"net={'on' if cfg.allow_network else 'off'}):\n  $ {cmd}"
            )
            try:
                res = _sandbox.run_agent_command(
                    cmd, Path(ctx.repo_dir or "."), config=cfg,
                )
            except Exception as exc:
                _append_system_message(f"\u26a0 Sandbox error: {exc}")
                input_textarea.value = f"Sandbox error: {exc}"
                _on_send(None)
                return

            if res.blocked:
                _append_system_message(
                    f"\u26d4 BLOCKED by allow-list: {res.block_reason}\n"
                    f"  $ {cmd}\n"
                    f"Run it manually in a terminal if you really need it."
                )
                input_textarea.value = (
                    f"That command was blocked by the allow-list "
                    f"({res.block_reason}). Please suggest a safer alternative."
                )
                _on_send(None)
                return

            if res.timed_out:
                _append_system_message(
                    f"\u26a0 Command timed out after {cfg.timeout_s}s."
                )
                input_textarea.value = (
                    "The command timed out. Please try a different approach."
                )
                _on_send(None)
                return

            output = (res.stdout.strip() + "\n" + res.stderr.strip()).strip()
            if res.returncode == 0:
                _append_system_message(
                    f"\u2714 Command succeeded ({res.elapsed_s:.1f}s):\n"
                    f"{output[:500]}"
                )
                input_textarea.value = (
                    f"I ran the command for you. Result:\n```\n{output[:1000]}\n```\n"
                    f"Continue with the next step."
                )
            else:
                _append_system_message(
                    f"\u2718 Command failed (exit {res.returncode}, "
                    f"{res.elapsed_s:.1f}s):\n{output[:500]}"
                )
                input_textarea.value = (
                    f"The command failed (exit {res.returncode}):\n"
                    f"```\n{output[:1000]}\n```\n"
                    f"Please suggest an alternative approach."
                )
            _on_send(None)
            return

        # --- File operations: upgrade permission profile if needed ---
        current_profile = state.get("_perm_profile", "ask_all")
        _PROFILE_RANK = {"plan": 0, "ask_all": 1, "repo_free": 2, "all_free": 3}
        current_rank = _PROFILE_RANK.get(current_profile, 1)
        # Edit/Write need "repo_free", Bash needs "all_free"
        needed_rank = 2 if tool in ("Edit", "Write", "Read", "Glob", "Grep", "") else 3
        need_upgrade = current_rank < needed_rank

        if need_upgrade:
            new_profile = "repo_free" if needed_rank == 2 else "all_free"
            _append_system_message(
                f"\u2705 Approved: {readable}\n"
                f"\u2191 Upgrading permissions: {current_profile} \u2192 {new_profile}"
            )
            old_engine = state["engine"]
            session_id = ""
            if old_engine:
                session_id = old_engine.session_id
            perm_dropdown.value = new_profile  # triggers _on_perm_change → syncs state
            engine = _ensure_engine()
            if engine and session_id:
                engine.session_id = session_id
            input_textarea.value = f"Please retry: {readable}"
            _on_send(None)
        else:
            _append_system_message(f"\u2705 Approved: {readable}")
            input_textarea.value = f"Yes, proceed with: {readable}"
            _on_send(None)

    def _on_deny(button):
        """User denies a blocked operation."""
        # Dashboard action denial
        dashboard_action = state.get("_pending_dashboard_action")
        if dashboard_action:
            _hide_approval()
            state["_pending_dashboard_action"] = None
            _append_system_message(f"\u274c Cancelled: {dashboard_action['description']}")
            return
        pending = state.get("_pending_approval")
        _hide_approval()
        raw = pending.get("tool", "") if pending else ""
        readable = _format_tool_description(raw) if raw else "operation"
        _append_system_message(f"\u274c Denied: {readable}")
        input_textarea.value = f"No, do NOT do that. Find an alternative approach."
        _on_send(None)

    # -- Interactive question detection & UI ------------------------------------

    def _detect_question(text: str) -> dict | None:
        """Detect if the agent's response ends with a question requiring user input.

        Returns a dict with 'type' and 'options' if a question is detected,
        or None if the response is a normal statement.

        Types:
        - 'numbered': Options like "1) foo  2) bar  3) baz"
        - 'yesno': Yes/no confirmation question
        - 'open': Open-ended question (ends with ?)
        """
        if not text or len(text) < 10:
            return None
        # Only look at the last ~2000 chars (the tail of the response)
        tail = text[-2000:].strip()
        # Skip if the response ended with a code block (likely not a question)
        if tail.rstrip().endswith("```"):
            return None

        # --- Numbered options: 1) / 1. / (1) patterns ---
        # Look for 2+ numbered items in the tail
        # Patterns: "1) text", "1. text", "(1) text", "**1.** text"
        option_patterns = [
            # "1) description" or "1. description" at line start
            re.compile(r'^\s*\*?\*?(\d+)[).]\*?\*?\s+(.+)', re.MULTILINE),
            # "(1) description" at line start
            re.compile(r'^\s*\((\d+)\)\s+(.+)', re.MULTILINE),
            # "- **Option 1**: description"
            re.compile(r'^\s*[-*]\s+\*?\*?(?:Option\s+)?(\d+)\*?\*?[.:]\s*(.+)', re.MULTILINE),
        ]
        for pat in option_patterns:
            matches = pat.findall(tail)
            if len(matches) >= 2:
                # Distinguish choosable options from numbered explanation steps.
                # Heuristic: real options are short (< 50 chars avg) and appear
                # near a question mark. Long numbered items are instructions.
                avg_len = sum(len(d.strip()) for _, d in matches) / len(matches)
                has_question = "?" in tail
                if avg_len > 120 or not has_question:
                    continue  # very long steps or no question context — skip
                options = []
                for num, desc in matches:
                    label = desc.strip().rstrip("*").strip()
                    options.append((num, label))
                return {"type": "numbered", "options": options}

        # --- QUESTION: tag (from solo_agent.md) ---
        if "QUESTION:" in tail:
            return {"type": "open", "options": []}

        # --- Yes/No questions ---
        last_lines = tail.split("\n")[-3:]
        last_text = " ".join(last_lines).strip().lower()
        yesno_indicators = [
            "shall i", "should i", "do you want me to", "would you like me to",
            "soll ich", "möchtest du dass ich", "willst du dass ich",
            "proceed?", "continue?", "go ahead?",
            "fortfahren?", "weitermachen?",
            "(yes/no)", "(y/n)", "(ja/nein)",
        ]
        if any(ind in last_text for ind in yesno_indicators) and "?" in last_text:
            return {"type": "yesno", "options": []}

        # --- Open question (ends with ?) ---
        # Only trigger if the very last meaningful line ends with ?
        for line in reversed(last_lines):
            line = line.strip()
            if line:
                if line.endswith("?"):
                    return {"type": "open", "options": []}
                break

        return None

    def _show_question_ui(question_info: dict):
        """Show interactive widgets based on detected question type.

        - numbered: Each option as a row with checkbox + quick-pick button
        - yesno: Yes / No buttons
        - open: Focus input area with hint
        """
        state["_pending_question"] = question_info
        qtype = question_info["type"]
        children = []

        if qtype == "numbered":
            options = question_info["options"]
            # Each option = one row: [checkbox] [quick-pick button]
            checkboxes = []
            rows = []
            for num, label in options:
                cb = widgets.Checkbox(
                    value=False,
                    description=label,
                    indent=False,
                    layout=widgets.Layout(width="auto", min_width="200px"),
                    style={"description_width": "0px"},
                )
                cb._option_num = num
                checkboxes.append(cb)
                btn = widgets.Button(
                    description=f"{num}",
                    button_style="info",
                    tooltip=f"Quick select {num}",
                    layout=widgets.Layout(width="32px", min_width="32px", height="26px"),
                )
                btn._option_value = num
                btn.on_click(_on_question_option)
                row = widgets.HBox(
                    [btn, cb],
                    layout=widgets.Layout(gap="6px", align_items="center"),
                )
                rows.append(row)
            state["_question_checkboxes"] = checkboxes

            # Submit button for multi-select
            submit_btn = widgets.Button(
                description="Submit selection",
                button_style="primary",
                layout=widgets.Layout(width="130px", height="30px"),
            )
            submit_btn.on_click(_on_question_submit)

            question_hint_html.value = (
                '<span style="font-size:11px;color:#3b82f6;font-weight:600;">'
                'Pick one or select multiple + Submit selection:</span>'
            )
            option_list = widgets.VBox(rows, layout=widgets.Layout(gap="0px"))
            children = [widgets.VBox([option_list, submit_btn], layout=widgets.Layout(gap="6px"))]

        elif qtype == "yesno":
            question_hint_html.value = (
                '<span style="font-size:11px;color:#3b82f6;font-weight:600;">'
                'Confirm:</span>'
            )
            yes_btn = widgets.Button(
                description="Yes",
                button_style="success",
                layout=widgets.Layout(width="70px", height="30px"),
            )
            yes_btn._option_value = "yes"
            yes_btn.on_click(_on_question_option)
            no_btn = widgets.Button(
                description="No",
                button_style="danger",
                layout=widgets.Layout(width="70px", height="30px"),
            )
            no_btn._option_value = "no"
            no_btn.on_click(_on_question_option)
            children = [yes_btn, no_btn]

        elif qtype == "open":
            question_hint_html.value = (
                '<span style="font-size:11px;color:#3b82f6;font-weight:600;">'
                'Agent awaits your answer</span>'
            )
            input_textarea.placeholder = "Type your answer..."

        question_buttons_box.children = children
        question_row.layout.display = "flex"

    def _on_question_submit(button):
        """Submit multi-select: gather checked options, send as answer."""
        parts = []
        checkboxes = state.get("_question_checkboxes", [])
        for cb in checkboxes:
            if cb.value:
                parts.append(getattr(cb, "_option_num", ""))
        if not parts:
            return  # nothing selected
        answer = ", ".join(parts)
        _hide_question_ui()
        input_textarea.value = answer
        _on_send(None)

    def _hide_question_ui():
        """Hide question UI and reset state."""
        question_row.layout.display = "none"
        question_hint_html.value = ""
        question_buttons_box.children = []
        state.pop("_pending_question", None)
        state.pop("_question_checkboxes", None)
        input_textarea.placeholder = "Message the agent... (Enter to send, Shift+Enter for newline)"

    def _on_question_option(button):
        """User clicked an option button — send as answer."""
        value = getattr(button, "_option_value", "")
        if not value:
            return
        _hide_question_ui()
        input_textarea.value = value
        _on_send(None)

    def _on_export(button):
        """Export button handler."""
        _export_chat()

    def _on_search_change(change):
        """Live search as user types."""
        _do_search(change["new"])

    def _on_search_close(button):
        """Close search bar."""
        _toggle_search()

    # -- pipeline display helpers ----------------------------------------------

    def _update_pipeline_display(eng):
        """Show visual pipeline progress: SM ✓ → Critic ✓ → Builder ⏳ → Test ○"""
        if not eng:
            return
        # Solo/Dashboard: no pipeline noise
        if mode_dropdown.value in ("solo", "dashboard"):
            return
        steps = eng.pipeline_status()
        _icons = {"done": "\u2705", "active": "\u23f3", "pending": "\u25cb"}
        parts = []
        for s in steps:
            icon = _icons.get(s["status"], "?")
            label = _format_role_label(s["role"])
            parts.append(f"{icon} {label}")
        pipeline_str = " \u2192 ".join(parts)
        _append_system_message(f"Pipeline: {pipeline_str}")

    def _extract_retry_context(agent_output: str, source: str) -> str:
        """Extract specific error context from agent output for targeted retry.

        Parses test failures (tracebacks, assertion errors) or reviewer
        findings (CRITICAL/MAJOR items) into a concise summary that helps
        the Builder fix the exact issues.
        """
        lines = agent_output.split("\n")
        findings = []

        if source == "test":
            # Extract pytest failures: FAILED lines, tracebacks, assertions
            in_failure = False
            failure_block: list[str] = []
            for line in lines:
                if "FAILED" in line or "FAIL" in line:
                    findings.append(line.strip())
                elif "Error" in line or "assert" in line.lower():
                    findings.append(line.strip())
                elif "Traceback" in line:
                    in_failure = True
                    failure_block = [line.strip()]
                elif in_failure:
                    failure_block.append(line.strip())
                    if line.strip() and not line.startswith(" "):
                        findings.append("\n".join(failure_block[-5:]))
                        in_failure = False
                        failure_block = []
            # Also grab lines with file:line references
            for line in lines:
                if re.match(r".*\.py:\d+", line) and "FAIL" in line.upper():
                    if line.strip() not in findings:
                        findings.append(line.strip())

        elif source == "reviewer":
            # Extract CRITICAL and MAJOR findings
            for line in lines:
                stripped = line.strip()
                if stripped.startswith(("1.", "2.", "3.", "4.", "5.")):
                    if any(kw in stripped.upper() for kw in
                           ("CRITICAL", "MAJOR", "BUG", "FIX:")):
                        findings.append(stripped)
                elif "CRITICAL" in stripped.upper() or "MAJOR" in stripped.upper():
                    findings.append(stripped)

        if not findings:
            # Fallback: last 500 chars of output
            return agent_output[-500:]

        return "\n".join(findings[:15])  # max 15 findings

    def _check_acceptance_gate(eng):
        """Check test agent output for acceptance criteria results."""
        test_out = eng.role_outputs.get("test_agent", "")
        if not test_out:
            return ""

        # Count PASS / FAIL / UNTESTED
        pass_count = len(re.findall(r"\bPASS\b", test_out))
        fail_count = len(re.findall(r"\bFAIL\b", test_out))
        untested = len(re.findall(r"\bUNTESTED\b", test_out))

        # Check for approve/reject verdict
        has_approve = bool(re.search(r"\*\*status:\*\*\s*approve", test_out, re.I))
        has_reject = bool(re.search(r"\*\*status:\*\*\s*reject", test_out, re.I))

        parts = []
        if pass_count:
            parts.append(f"{pass_count} PASS")
        if fail_count:
            parts.append(f"{fail_count} FAIL")
        if untested:
            parts.append(f"{untested} UNTESTED")

        summary = f"({', '.join(parts)})" if parts else ""

        if has_reject or fail_count > 0:
            return f"\u274c {summary}"
        if has_approve:
            return f"\u2705 {summary}"
        return summary

    # -- main event handlers -----------------------------------------------

    def _on_send(button):
        user_text = input_textarea.value.strip()
        if not user_text:
            return

        # Hide any pending question UI when user sends a message
        _hide_question_ui()

        # Slash commands execute immediately, even during streaming.
        # But file paths like /home/... are NOT slash commands — only
        # short tokens like /help, /calc, /orca etc. count.
        _first_token = user_text.split()[0].lower() if user_text else ""
        _SLASH_PREFIXES = {
            "/help", "/clear", "/cost", "/compact", "/stop", "/status",
            "/usage", "/export", "/search", "/retry", "/undo", "/git", "/provider",
            "/model", "/effort", "/mode", "/perms", "/perm-cycle", "/reset",
            "/memories", "/remember", "/forget", "/plans", "/plan", "/hooks",
            "/session", "/mcp", "/commands", "/init", "/bash", "/failures",
            "/workspace", "/tab", "/ui",
            "/control", "/submit", "/orca", "/jobs", "/calc", "/analyze",
            "/recalc", "/cancel", "/context", "/agents", "/skills",
        }
        # User-defined slash commands and skill expansion: when /<name>
        # doesn't match a built-in slash command, first look in the
        # commands/ markdown templates (lightweight $ARGUMENTS-style
        # substitution), then fall back to skills/ markdown (richer
        # Skill object). Either way the user's chat input gets rewritten
        # to the expanded body before reaching the agent.
        if (user_text.startswith("/")
                and _first_token not in _SLASH_PREFIXES
                and "/" not in _first_token[1:]
                and len(_first_token) > 1
                and not _first_token.startswith("/.")):
            _engine = state.get("engine")
            _ws = None
            if _engine and getattr(_engine, "kit_permissions", None):
                _ws = _engine.kit_permissions.workspace
            _name = _first_token[1:]
            _rest = user_text[len(_first_token):].strip()
            # 1) custom command template (commands/<name>.md)
            try:
                from delfin.agent import slash_commands as _scm
                _tpl = _scm.get_command(_name, _ws)
            except Exception:
                _tpl = None
            if _tpl is not None:
                user_text = _scm.expand_template(_tpl.body, _rest)
            else:
                # 2) fall back to full Skill (skills/<name>/SKILL.md)
                try:
                    from delfin.agent import skills as _skills_mod
                    _sk = _skills_mod.get_skill(_name, _ws)
                except Exception:
                    _sk = None
                if _sk is not None:
                    expanded = _skills_mod.render_skill_invocation(_sk, _rest)
                    user_text = expanded + ("\n\n" + _rest if _rest else "")
        elif _first_token in _SLASH_PREFIXES:
            input_textarea.value = ""
            _append_chat_message("user", user_text)
            _handle_slash_command(user_text)
            return

        # If streaming, queue the message instead of interrupting.
        # The agent will receive it after finishing the current response.
        # Use /stop to force-interrupt if needed.
        if state["streaming"]:
            # Safety: detect stale streaming state (worker crashed/finished
            # but streaming flag wasn't reset). If no active CLI process,
            # reset the flag and process normally.
            _engine = state.get("engine")
            _stale = False
            if _engine and hasattr(_engine, "client"):
                _proc = getattr(_engine.client, "_proc", None)
                if _proc is not None and _proc.poll() is not None:
                    _stale = True  # process is dead but streaming still True
                elif _proc is None:
                    _stale = True  # no process at all
            if _stale:
                with state["_state_lock"]:
                    state["streaming"] = False
                _set_working(False)
                # Fall through to send normally
            else:
                state["message_queue"].append(user_text)
                input_textarea.value = ""
                _append_chat_message("user", user_text)
                _append_system_message(
                    f"\U0001f4e8 Queued — agent will receive this after finishing. "
                    f"Type /stop to interrupt."
                )
                _update_queue_display()
                return

        # Track user message for safety intent-checking
        state["_last_user_message"] = user_text

        engine = _ensure_engine()
        if engine is None:
            return

        # No auto-mode-switching — the user's chosen mode is always respected.

        # After cycle complete: enter follow-up mode (keep context alive).
        # The user can continue chatting with the builder/solo agent.
        # Use /reset to start a truly fresh cycle.
        if engine.is_cycle_complete and not state.get("_follow_up"):
            # Find best follow-up role: builder > solo > last in route
            _fu_role = "builder_agent"
            if _fu_role not in engine.route:
                _fu_role = "solo_agent" if "solo_agent" in engine.route else engine.route[-1]
            _fu_idx = engine.route.index(_fu_role) if _fu_role in engine.route else 0
            engine.current_role_index = _fu_idx
            state["_follow_up"] = True
            # Solo/Dashboard: silent follow-up (no pipeline noise)
            if mode_dropdown.value not in ("solo", "dashboard"):
                _append_system_message(
                    f"--- Follow-up mode ({_fu_role.replace('_', ' ').title()}) "
                    f"— continue chatting or /reset for new task ---"
                )

        # Detect user approval to start pipeline from Session Manager
        # If SM has already responded and user sends a short confirmation,
        # advance past SM and auto-run the remaining pipeline.
        _APPROVAL_WORDS = {
            "ja", "yes", "go", "start", "los", "mach", "weiter", "ok",
            "okay", "passt", "approved", "proceed", "ja bitte", "los gehts",
            "do it", "sieht gut aus", "einverstanden", "starten", "anfangen",
            "beginne", "lgtm", "ship it", "machen",
        }
        _sm_approval = False
        if (engine.current_role == "session_manager"
                and len(engine.messages) >= 2
                and len(user_text) < 80):
            _lower = user_text.lower().strip().rstrip("!.?")
            if _lower in _APPROVAL_WORDS:
                _sm_approval = True

        # Handle QUESTION: protocol — agent asked user, user replied
        _question_role = state.pop("_awaiting_agent_question", None)
        if _question_role:
            # Inject user answer back to same agent (don't advance)
            current_msg = f"[User answer to your question]\n{user_text}"
            state["_question_answer"] = current_msg
            _set_active_gate()
            _update_status()
            _role_lbl = _format_role_label(_question_role)
            _append_system_message(f"Antwort an {_role_lbl} gesendet.")

        # Handle conflict resolution
        if state.pop("_awaiting_conflict_resolution", False):
            _lower_msg = user_text.lower().strip().rstrip("!.?")
            if "reject" in _lower_msg or "stop" in _lower_msg or "nein" in _lower_msg:
                _append_system_message("Pipeline stopped by user. Use /reset for new task.")
                if engine:
                    engine.request_stop()
                input_textarea.value = ""
                return
            # else: continue to test
            _set_active_gate()
            _update_status()

        # Handle findings review response
        _findings_review_role = state.pop("_awaiting_findings_review", None)
        if _findings_review_role:
            state["_findings_approved"] = True
            _lower_msg = user_text.lower().strip().rstrip("!.?")
            # If user wants to skip/filter findings, modify the role output
            if _lower_msg not in _APPROVAL_WORDS and "skip" in _lower_msg:
                # Append user's filter instructions to the role output
                # so Builder sees them alongside the findings
                old_output = engine.role_outputs.get(_findings_review_role, "")
                engine.role_outputs[_findings_review_role] = (
                    old_output + f"\n\n--- USER FILTER ---\n{user_text}\n"
                    f"Address only the findings NOT mentioned above as skipped."
                )
            _set_active_gate()
            _update_status()

        # Handle generic gate review response
        _gate_review_role = state.pop("_awaiting_gate_review", None)
        if _gate_review_role:
            _lower_msg = user_text.lower().strip().rstrip("!.?")
            if _lower_msg not in _APPROVAL_WORDS:
                # Not an approval — re-send to same agent with user guidance
                _append_system_message(
                    f"Pipeline paused after {_format_role_label(_gate_review_role)}. "
                    f"Reply 'go' to continue, or /reset to stop this cycle."
                )
                input_textarea.value = ""
                return
            # Gate approved — advance to next role instead of re-sending
            # "go" to the same agent (which would re-review unchanged code).
            _set_active_gate()
            prev_label = _format_role_label(engine.current_role)
            engine.compact_for_next_role()
            has_next = engine.advance_role()
            if has_next:
                next_label = _format_role_label(engine.current_role)
                _record_cycle_event("handoff", f"{prev_label} -> {next_label}")
                _append_system_message(
                    f"--- Gate approved: {prev_label} \u2192 {next_label} ---"
                )
                # Replace user_text with handoff message for the next agent
                original_task = ""
                for msg in state.get("_chat_messages", []):
                    if msg.get("role") == "user":
                        original_task = msg.get("content", "")
                        break
                user_text = engine.build_handoff_message(original_task)
            else:
                _record_cycle_event("cycle", "Cycle complete")
                _append_system_message("--- Cycle complete ---")
                _update_status()
                _update_pipeline_display(engine)
                input_textarea.value = ""
                return
            _update_status()
            _update_pipeline_display(engine)

        input_textarea.value = ""
        _append_chat_message("user", user_text)

        # Fire UserPromptSubmit hooks. Block reasons surface as a system
        # message and abort the send; non-blocking hook output (stderr,
        # non-zero exits without "block" decision) also surfaces as an
        # info note so the user actually sees their hooks running.
        try:
            from delfin.agent import hooks as _hooks_mod
            _cfg = _hooks_mod.load_hooks(ctx.repo_dir or None)
            if not _cfg.is_empty():
                _ups = _hooks_mod.run_hooks(
                    "UserPromptSubmit", _cfg,
                    user_prompt=user_text,
                    workspace=ctx.repo_dir or None,
                )
                _blk = _hooks_mod.first_block(_ups)
                if _blk is not None:
                    _append_system_message(
                        "⛔ UserPromptSubmit hook blocked the message:\n"
                        f"  reason: {(_blk.reason or _blk.stderr or '').strip()[:200]}"
                    )
                    return
                # Non-blocking, but report if any hook produced visible output
                _noisy = [r for r in _ups
                          if r.matched and (r.stderr.strip() or r.exit_code)]
                if _noisy:
                    _lines = ["🎣 UserPromptSubmit hooks fired:"]
                    for r in _noisy:
                        _lines.append(
                            f"  {r.command[:60]}  exit={r.exit_code}  "
                            f"{r.duration_s:.1f}s"
                        )
                        if r.stderr.strip():
                            _lines.append(f"    stderr: {r.stderr.strip()[:120]}")
                    _append_system_message("\n".join(_lines))
        except Exception:
            pass

        with state["_state_lock"]:
            state["streaming"] = True
            state["_generation_id"] = state.get("_generation_id", 0) + 1
            _my_gen_id = state["_generation_id"]
            state["_deny_count"] = 0  # Reset retry counter for new message
        if state["session_start_time"] is None:
            state["session_start_time"] = time.monotonic()
        _ensure_task_session_id(engine, create=True)
        _update_button_states()
        # Mark turn start + arm the stale-stream watcher.
        state["_last_stream_activity"] = time.monotonic()
        state["_stale_seen"] = False
        _arm_stale_watcher()
        # Snapshot the engine's pre-turn cost/tokens so the worker can
        # emit a per-turn delta footer once the turn finishes.
        try:
            _pre = engine.get_status()
            state["_turn_pre_cost"] = float(_pre.get("cost_usd") or 0.0)
            state["_turn_pre_in"] = int(_pre.get("input_tokens") or 0)
            state["_turn_pre_out"] = int(_pre.get("output_tokens") or 0)
            state["_turn_started_monotonic"] = time.monotonic()
        except Exception:
            state["_turn_pre_cost"] = 0.0
            state["_turn_pre_in"] = 0
            state["_turn_pre_out"] = 0
            state["_turn_started_monotonic"] = time.monotonic()
        _set_working(True, "Thinking...")

        role_label = _format_role_label(engine.current_role)

        def _worker():
            nonlocal _sm_approval
            chunks = []
            thinking_chunks = []
            tool_count = [0]  # mutable counter for tool calls in this turn
            last_update = 0.0
            try:

                def _on_thinking(text):
                    nonlocal last_update
                    state["_last_stream_activity"] = time.monotonic()
                    if state.get("_stale_seen"):
                        # Recovered from stale — clear the warning state
                        state["_stale_seen"] = False
                    thinking_chunks.append(text)
                    now = time.monotonic()
                    if now - last_update > 0.15:
                        # Show last ~80 chars of thinking in spinner
                        full = "".join(thinking_chunks)
                        snippet = full[-80:].replace("\n", " ").strip()
                        if len(full) > 80:
                            snippet = "..." + snippet
                        _set_working(True, f"Thinking: {snippet}")
                        last_update = now

                def _on_token(text):
                    nonlocal last_update
                    state["_last_stream_activity"] = time.monotonic()
                    if state.get("_stale_seen"):
                        state["_stale_seen"] = False
                    # When first text arrives, flush thinking as collapsed block
                    if thinking_chunks and not chunks:
                        full_thinking = "".join(thinking_chunks)
                        if full_thinking.strip():
                            # Show as collapsed thinking block in chat
                            _append_chat_message(
                                "thinking",
                                full_thinking.strip(),
                            )
                        thinking_chunks.clear()
                    chunks.append(text)
                    now = time.monotonic()
                    # Flush any pending tool message refreshes
                    if state.get("_html_refresh_pending"):
                        state["_html_refresh_pending"] = False
                        _refresh_chat_html()
                        state["_last_html_refresh"] = now
                    if now - last_update > 0.25:
                        _set_working(True, "Writing...")
                        _update_last_assistant("".join(chunks), role_label)
                        _update_status()
                        last_update = now

                def _on_tool_use(tool_name, tool_input):
                    state["_last_stream_activity"] = time.monotonic()
                    if state.get("_stale_seen"):
                        state["_stale_seen"] = False
                    # Flush thinking if no text came before tool use
                    if thinking_chunks:
                        full_thinking = "".join(thinking_chunks)
                        if full_thinking.strip():
                            _append_chat_message("thinking", full_thinking.strip())
                        thinking_chunks.clear()
                    # KIT-Toolbox emits MCP-style names ("mcp__kit-coding__
                    # edit_file"); strip the prefix so the renderer cascade
                    # below can match on the bare tool name and produce the
                    # same diff-block UX as Claude CLI's Edit/Write tools.
                    if tool_name and tool_name.startswith("mcp__"):
                        parts = tool_name.split("__")
                        if len(parts) >= 3:
                            tool_name = parts[-1]
                    # Parse tool input
                    try:
                        import json as _j
                        parsed = _j.loads(tool_input)
                    except Exception:
                        parsed = {}

                    # Spinner detail with tool counter
                    tool_count[0] += 1
                    _tc = tool_count[0]
                    _fname = (parsed.get("file_path") or parsed.get("path") or "")
                    if _fname:
                        _fname = _fname.rsplit("/", 1)[-1]
                    _detail = {
                        "Read":  f"Reading {_fname}..." if _fname else "Reading...",
                        "Edit":  f"Editing {_fname}..." if _fname else "Editing...",
                        "Write": f"Writing {_fname}..." if _fname else "Writing...",
                        "Grep":  f"Searching: {parsed.get('pattern', '')[:40]}...",
                        "Glob":  f"Finding: {parsed.get('pattern', '')[:40]}...",
                        "Bash":  f"$ {(parsed.get('command') or '')[:50]}...",
                        "Agent": f"Sub-agent: {(parsed.get('description') or '')[:40]}...",
                    }.get(tool_name, f"Running {tool_name}...")
                    _set_working(True, f"[{_tc}] {_detail}")

                    # Flush pending text as a finalized assistant message
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label, finalize=True)
                        chunks.clear()

                    # --- Build terminal-style tool display ---
                    _e = _html.escape

                    def _short_path(p, segments=3):
                        """Last N path segments for readability."""
                        parts = p.rstrip("/").split("/")
                        return "/".join(parts[-segments:]) if len(parts) > segments else p

                    if tool_name == "Read":
                        fpath = parsed.get("file_path", "")
                        sp = _short_path(fpath)
                        line_info = ""
                        if parsed.get("offset"):
                            lim = parsed.get("limit", "")
                            line_info = f":{parsed['offset']}"
                            if lim:
                                line_info += f"-{parsed['offset'] + lim}"
                        _append_tool_message(
                            f'<span class="tool-name">Read</span>  '
                            f'<span class="tool-path">{_e(sp)}{_e(line_info)}</span>'
                        )

                    elif tool_name == "Grep":
                        pat = parsed.get("pattern", "")
                        path = _short_path(parsed.get("path", "."))
                        glb = parsed.get("glob", "")
                        typ = parsed.get("type", "")
                        scope = f" in {_e(path)}" if path != "." else ""
                        filt = ""
                        if glb:
                            filt = f"  ({_e(glb)})"
                        elif typ:
                            filt = f"  (*.{_e(typ)})"
                        _append_tool_message(
                            f'<span class="tool-name">Grep</span>  '
                            f'<span class="tool-param">"{_e(pat[:60])}"</span>'
                            f'{scope}{filt}'
                        )

                    elif tool_name == "Glob":
                        pat = parsed.get("pattern", "")
                        path = _short_path(parsed.get("path", "."))
                        scope = f"  in {_e(path)}" if path != "." else ""
                        _append_tool_message(
                            f'<span class="tool-name">Glob</span>  '
                            f'<span class="tool-param">{_e(pat)}</span>{scope}'
                        )

                    elif tool_name == "Bash":
                        cmd = parsed.get("command", "")
                        # Show full command up to 200 chars
                        if len(cmd) > 200:
                            cmd = cmd[:200] + "..."
                        _append_tool_message(
                            f'<span class="tool-name">$</span> {_e(cmd)}'
                        )

                    elif tool_name in (
                        "Edit", "Write",
                        "edit_file", "write_file", "multi_edit",
                    ):
                        # KIT-Toolbox uses 'path'; Claude CLI uses 'file_path'.
                        fpath = parsed.get("file_path") or parsed.get("path") or ""
                        sp = _short_path(fpath)

                        # Display label normalises capitalised + lowercase forms
                        # so the user sees a consistent "Edit"/"Write"/"MultiEdit"
                        # banner regardless of which provider produced the call.
                        label = {
                            "Edit": "Edit", "edit_file": "Edit",
                            "Write": "Write", "write_file": "Write",
                            "multi_edit": "MultiEdit",
                        }.get(tool_name, tool_name)

                        if tool_name in ("Edit", "edit_file"):
                            old = parsed.get("old_string", "") or ""
                            new = parsed.get("new_string", "") or ""
                            diff_html = _build_inline_diff(old, new, _e)
                            _append_tool_message(
                                f'<span class="tool-name">{label}</span>  '
                                f'<span class="tool-path">{_e(sp)}</span>\n'
                                f'{diff_html}'
                            )
                        elif tool_name == "multi_edit":
                            edits = parsed.get("edits", []) or []
                            blocks = []
                            for i, ed in enumerate(edits):
                                if not isinstance(ed, dict):
                                    continue
                                old = ed.get("old_string", "") or ""
                                new = ed.get("new_string", "") or ""
                                blocks.append(
                                    f'<div style="font-size:10px; '
                                    f'color:#9ca3af; margin-top:4px;">'
                                    f'edit {i+1}/{len(edits)}</div>'
                                    f'{_build_inline_diff(old, new, _e)}'
                                )
                            body = "\n".join(blocks) if blocks else ""
                            _append_tool_message(
                                f'<span class="tool-name">{label}</span>  '
                                f'<span class="tool-path">{_e(sp)}</span>  '
                                f'<span class="tool-param">'
                                f'({len(edits)} edits)</span>\n{body}'
                            )
                        else:
                            # Write or write_file — render full content as
                            # a single + diff so the user sees what's being
                            # written, not just a chars count.
                            content = parsed.get("content", "") or ""
                            diff_html = _build_inline_diff("", content, _e)
                            _append_tool_message(
                                f'<span class="tool-name">{label}</span>  '
                                f'<span class="tool-path">{_e(sp)}</span>  '
                                f'<span class="tool-param">'
                                f'({len(content)} chars, '
                                f'{content.count(chr(10)) + 1} lines)</span>\n'
                                f'{diff_html}'
                            )

                        state["recent_edits"].append({
                            "file": fpath, "tool": tool_name,
                        })
                        undo_btn.disabled = False

                    elif tool_name == "bash":
                        cmd = parsed.get("command", "") or ""
                        if len(cmd) > 200:
                            cmd = cmd[:200] + "..."
                        _append_tool_message(
                            f'<span class="tool-name">$</span> {_e(cmd)}'
                        )

                    elif tool_name == "Agent":
                        desc = parsed.get("description", "")
                        prompt = parsed.get("prompt", "")[:80]
                        sub_type = parsed.get("subagent_type", "general-purpose")
                        # Append to subagent panel state and refresh widget.
                        state["subagent_calls"].append({
                            "subagent_type": sub_type,
                            "description": desc,
                            "prompt": parsed.get("prompt", ""),
                            "status": "running",
                            "output": "",
                        })
                        _pane_html = _render_subagent_pane_html(state["subagent_calls"])
                        if _pane_html:
                            subagent_pane_html.value = _pane_html
                            subagent_pane_html.layout.display = "block"
                        # Keep the inline tool message for chat-history continuity
                        _append_tool_message(
                            f'<span class="tool-name">Agent</span>  '
                            f'<span class="tool-param">[{_e(sub_type)}] {_e(desc)}</span>'
                            + (f'\n  {_e(prompt)}...' if prompt else '')
                        )

                    elif tool_name == "WebSearch":
                        query = parsed.get("query", "")
                        _append_tool_message(
                            f'<span class="tool-name">WebSearch</span>  '
                            f'<span class="tool-param">"{_e(query[:80])}"</span>'
                        )

                    elif tool_name == "WebFetch":
                        url = parsed.get("url", "")
                        _append_tool_message(
                            f'<span class="tool-name">WebFetch</span>  '
                            f'<span class="tool-path">{_e(url[:100])}</span>'
                        )

                    elif tool_name == "TodoWrite":
                        todos = parsed.get("todos") or []
                        # Persist current todo plan for inspectors / status line
                        state["current_todos"] = todos
                        # Update the persistent plan pane (visible while agent works)
                        _pane_html = _render_todo_pane_html(todos)
                        if _pane_html:
                            todo_pane_html.value = _pane_html
                            todo_pane_html.layout.display = "block"
                        else:
                            todo_pane_html.value = ""
                            todo_pane_html.layout.display = "none"
                        if not todos:
                            _append_tool_message(
                                '<span class="tool-name">TodoWrite</span>  '
                                '<span class="tool-param">cleared</span>'
                            )
                        else:
                            _STATUS_GLYPH = {
                                "completed":   "[x]",
                                "in_progress": "[~]",
                                "pending":     "[ ]",
                            }
                            _STATUS_COLOR = {
                                "completed":   "tool-diff-new",
                                "in_progress": "tool-name",
                                "pending":     "tool-param",
                            }
                            n_done = sum(1 for t in todos if t.get("status") == "completed")
                            n_total = len(todos)
                            lines = [
                                f'<span class="tool-name">TodoWrite</span>  '
                                f'<span class="tool-param">{n_done}/{n_total} done</span>'
                            ]
                            for t in todos:
                                status = t.get("status", "pending")
                                glyph = _STATUS_GLYPH.get(status, "[ ]")
                                cls = _STATUS_COLOR.get(status, "tool-param")
                                # Show activeForm while running, content otherwise
                                if status == "in_progress":
                                    text = t.get("activeForm") or t.get("content", "")
                                else:
                                    text = t.get("content", "")
                                text = text[:120]
                                lines.append(
                                    f'  <span class="{cls}">{_e(glyph)}</span> '
                                    f'<span class="tool-param">{_e(text)}</span>'
                                )
                            _append_tool_message("\n".join(lines))

                    else:
                        # Generic: show first meaningful param
                        param = ""
                        for key in ("file_path", "path", "pattern",
                                    "prompt", "description", "query", "command"):
                            if key in parsed:
                                val = str(parsed[key])[:80]
                                param = f'<span class="tool-param">{_e(val)}</span>'
                                break
                        if not param:
                            param = _e(tool_input[:100])
                        _append_tool_message(
                            f'<span class="tool-name">{_e(tool_name)}</span>  {param}'
                        )

                def _emit_new_workspace_artifacts() -> None:
                    """D4: scan agent_workspace for newly created files and
                    append inline previews (PNG/SVG/CSV/JSON) to the chat."""
                    try:
                        before = state.get("_ws_known_files") or set()
                        current = {
                            p.name: p for p in ctx.agent_dir.iterdir()
                            if p.is_file()
                        }
                    except Exception:
                        return
                    new_names = set(current.keys()) - set(before)
                    if not new_names:
                        return
                    for name in sorted(new_names):
                        path = current.get(name)
                        if path is None:
                            continue
                        try:
                            html_block = _render_artifact_inline(path)
                        except Exception:
                            html_block = None
                        if html_block:
                            _append_tool_message(html_block)
                    # Update the snapshot so we don't re-emit the same artifact.
                    state["_ws_known_files"] = set(current.keys())

                def _on_tool_result(tool_name, tool_output):
                    """Append tool result as collapsible detail to the last tool message."""
                    if not tool_output:
                        return
                    if tool_name and tool_name.startswith("mcp__"):
                        parts = tool_name.split("__")
                        if len(parts) >= 3:
                            tool_name = parts[-1]
                    # Subagent panel: mark the next running call as done.
                    if tool_name == "Agent":
                        for entry in state["subagent_calls"]:
                            if entry.get("status") == "running":
                                entry["status"] = "done"
                                entry["output"] = tool_output
                                break
                        _pane_html = _render_subagent_pane_html(state["subagent_calls"])
                        if _pane_html:
                            subagent_pane_html.value = _pane_html
                            subagent_pane_html.layout.display = "block"
                    # Phase 5d: refresh the live task ticker whenever a
                    # task_* tool fires so the panel reflects new state
                    # without waiting for the next user turn.
                    if tool_name in ("task_create", "task_update", "task_list"):
                        try:
                            _refresh_task_ticker()
                        except Exception:
                            pass
                    # Truncate for display
                    output = tool_output
                    _MAX_LINES = 8
                    _MAX_CHARS = 600
                    lines = output.split("\n")
                    truncated = False
                    if len(lines) > _MAX_LINES:
                        output = "\n".join(lines[:_MAX_LINES])
                        truncated = True
                    if len(output) > _MAX_CHARS:
                        output = output[:_MAX_CHARS]
                        truncated = True
                    suffix = ""
                    if truncated:
                        total = len(tool_output)
                        suffix = f"\n... ({total:,} chars total)"
                    # Integrate into the last tool message as collapsible detail
                    msgs = state["chat_messages"]
                    if msgs and msgs[-1]["role"] == "tool":
                        first_line = output.split("\n")[0][:80]
                        detail_html = (
                            f'<details><summary> → {_html.escape(first_line)}'
                            f'{"..." if truncated or len(output.split(chr(10))) > 1 else ""}'
                            f'</summary>'
                            f'<pre>{_html.escape(output)}{_html.escape(suffix)}</pre>'
                            f'</details>'
                        )
                        msgs[-1]["content"] += detail_html
                        _refresh_chat_html()
                    else:
                        _append_tool_message(
                            f'<span style="color:#9ca3af;">'
                            f'{_html.escape(output)}</span>'
                        )

                    # D4: After file-creating tools, surface any new artifacts
                    # in the chat so the user sees plots/CSVs without
                    # leaving the conversation. The MCP plot tool writes
                    # PNGs straight into agent_workspace/, so it qualifies
                    # too; check the tool-name suffix because MCP tools
                    # are namespaced as ``mcp__<server>__<name>``.
                    _is_file_creating_tool = (
                        tool_name in ("Write", "Edit", "Bash", "NotebookEdit")
                        or tool_name.endswith("__plot_energy_distribution")
                        or tool_name.endswith("__plot_energy_correlation")
                        or tool_name.endswith("__plot_orbital_diagram")
                        or tool_name.endswith("__plot_optimization_convergence")
                        or tool_name.endswith("__plot_uvvis_spectrum")
                        or tool_name.endswith("__plot_scf_convergence")
                        or tool_name.endswith("__plot_population_charges")
                        or tool_name.endswith("__plot_vibrational_spectrum")
                    )
                    if _is_file_creating_tool:
                        _emit_new_workspace_artifacts()

                def _on_permission_denied(description):
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label)
                        chunks.clear()
                    denied_raw = str(description)
                    readable = _format_tool_description(denied_raw)
                    state["_last_denied"] = denied_raw
                    # Step 5 — distinguish "interactive approval" from
                    # "structurally blocked": if the denied tool isn't in
                    # the CLI allowlist for the current mode, Approve/Deny
                    # can't help. Tell the user what to switch instead.
                    _engine_obj = state.get("engine")
                    _allowed = getattr(
                        getattr(_engine_obj, "client", None), "allowed_tools", None,
                    )
                    if _is_structurally_blocked(denied_raw, _allowed):
                        _tool = _extract_denied_tool_name(denied_raw) or "this tool"
                        _append_system_message(
                            f"⛔ {_tool} is blocked by the {mode_dropdown.value} "
                            f"mode CLI allowlist — Approve won't help. "
                            f"Switch the Mode dropdown to **solo** "
                            f"(live, conversation context is preserved) "
                            f"or run the action manually."
                        )
                        engine.request_stop()
                        return
                    # Accumulate denied commands for context injection
                    _denied_cmds = state.setdefault("_denied_commands", [])
                    if readable not in _denied_cmds:
                        _denied_cmds.append(readable)
                    # Track denials to stop retry loops (1 denial = stop)
                    with state["_state_lock"]:
                        deny_count = state.get("_deny_count", 0) + 1
                        state["_deny_count"] = deny_count
                    if deny_count >= 2:
                        _append_system_message(
                            f"\u26d4 Blocked {deny_count}x — stopping. "
                            f"Change permission mode or do it manually."
                        )
                        # Inject denial context so the model knows on the
                        # next turn that commands were blocked.
                        _denied_list = ", ".join(_denied_cmds[-5:])
                        engine.messages.append({
                            "role": "user",
                            "content": (
                                f"[System] Bash commands were BLOCKED by the "
                                f"permission system: {_denied_list}. "
                                f"Do NOT retry these or any variation. "
                                f"Tell the user what commands to run manually."
                            ),
                        })
                        engine.request_stop()
                        return
                    _append_system_message(
                        f"\u26d4 Blocked ({deny_count}/2): {readable}"
                    )
                    # Show approval buttons for interactive approval
                    _show_approval_prompt(denied_raw, denied_raw)

                # Effort multiplier: scales the role-specific budget.
                # xhigh applies an aggressive multiplier and a higher cap.
                _effort_mult = {"low": 0.5, "medium": 1.0, "high": 2.0, "xhigh": 3.0}
                # Dashboard mode auto-clamps effort: tab-open / control-key /
                # submit-job dispatch never needs deep reasoning, and high-
                # effort against a reasoning model (Azure GPT-5.x, o-series)
                # routinely takes 2+ minutes for what should be a 100 ms
                # routing decision. The user can still override via the
                # dropdown but the default is "low" so the typical case is
                # fast + cheap.
                _effective_effort = effort_dropdown.value
                if (mode_dropdown.value == "dashboard"
                        and _effective_effort in ("high", "xhigh")):
                    _effective_effort = "low"
                _mult = _effort_mult.get(_effective_effort, 1.0)
                _budget_cap = 200_000 if _effective_effort == "xhigh" else 128_000

                # Store original user task for handoff messages
                original_task = user_text
                # S1 — Per-turn live state goes into the SYSTEM prompt via
                # engine.set_live_state(), not into the user message body.
                # That keeps engine.messages history small and cache-friendly:
                # old turns no longer carry their stale dashboard state.
                _live_state_text = ""
                if mode_dropdown.value == "dashboard":
                    try:
                        _live_state_text = _build_dashboard_context()
                    except Exception:
                        _live_state_text = ""
                elif mode_dropdown.value == "solo":
                    try:
                        _solo_snap = _collect_solo_domain_snapshot()
                        _live_state_text = _format_solo_domain_state(_solo_snap)
                    except Exception:
                        _live_state_text = ""
                if hasattr(engine, "set_live_state"):
                    engine.set_live_state(_live_state_text)
                current_msg = user_text

                # S4 — Session boot context: ONCE on the first user-send of
                # a fresh dashboard session. Gives the agent a 600-token
                # primer (recent outcomes, SLURM jobs, commits, calc-folder)
                # so it doesn't have to fire a dozen tool calls just to
                # know "where are we and what just failed?".
                if mode_dropdown.value == "dashboard" and not engine.messages:
                    try:
                        _boot = _build_dashboard_session_boot()
                    except Exception:
                        _boot = ""
                    if _boot:
                        current_msg = f"{_boot}\n\n{current_msg}"
                elif mode_dropdown.value == "solo" and not engine.messages:
                    try:
                        _boot = _build_solo_session_boot()
                    except Exception:
                        _boot = ""
                    if _boot:
                        current_msg = f"{_boot}\n\n{current_msg}"

                # Live mode-switch handoff: if the user just changed mode,
                # _on_mode_change stashed a full-transcript block here.
                # Prepend it to the first user message of the new engine so
                # the model picks up where the previous mode left off.
                _handoff = state.pop("_pending_mode_handoff", "")
                if _handoff:
                    current_msg = f"{_handoff}\n\n{current_msg}"

                # D4: snapshot agent_workspace before the turn so we can
                # diff and inline-render any new artifacts the agent
                # creates (PNG/SVG/CSV/JSON).
                try:
                    state["_ws_known_files"] = {
                        p.name for p in ctx.agent_dir.iterdir() if p.is_file()
                    }
                except Exception:
                    state["_ws_known_files"] = set()

                _turn_start_time = time.monotonic()
                max_auto_steps = len(engine.route) + 1  # safety limit

                for _step in range(max_auto_steps):
                    if engine._stop_requested or engine.is_cycle_complete:
                        break

                    # Role-specific thinking budget and model routing
                    from delfin.agent.engine import AgentEngine as _AE
                    _cur_role = engine.current_role
                    _task_class = ""
                    try:
                        _route_info = _AE.recommend_task_route(
                            original_task or current_msg,
                            engine.mode,
                            is_delfin_workspace=getattr(
                                engine, "_is_delfin_workspace", True,
                            ),
                        )
                        _task_class = _route_info.get("task_class", "")
                    except Exception:
                        pass
                    # Solo: adaptive thinking (0 = CLI manages it)
                    if _cur_role == "solo_agent":
                        _budget = 0
                    else:
                        _base_budget = _AE.thinking_budget_for_role(
                            _cur_role,
                            task_class=_task_class,
                        )
                        _budget = min(int(_base_budget * _mult), _budget_cap)

                    # Per-role model: switch to optimal model (Claude only)
                    _effective_model = model_dropdown.value
                    if provider_dropdown.value == "claude":
                        _role_model = _AE.model_for_role(_cur_role)
                        _agent_settings = _get_agent_settings()
                        _role_models_cfg = _agent_settings.get("role_models", {})
                        if _cur_role in _role_models_cfg:
                            _role_model = _role_models_cfg[_cur_role]
                        _user_model = model_dropdown.value
                        _effective_model = _user_model if _role_model == "auto" else _role_model
                        if (hasattr(engine.client, "switch_model")
                                and _effective_model != getattr(engine.client, "model", "")):
                            engine.client.switch_model(_effective_model)

                    # Track per-role costs
                    _cost_before = engine.cost_usd
                    _in_before = engine.token_usage["input"]
                    _out_before = engine.token_usage["output"]

                    role_label = _format_role_label(engine.current_role)
                    chunks.clear()
                    thinking_chunks.clear()

                    # Load persistent memory for system prompt
                    from delfin.agent.memory_store import format_memory_context
                    from delfin.agent.project_memory import load_project_memory
                    _memory = format_memory_context(task_text=original_task)
                    # Auto-load CLAUDE.md / AGENTS.md / DELFIN.md from cwd up.
                    try:
                        _kperms = engine.kit_permissions
                        _extra = list(_kperms.extra_workspace_dirs) if _kperms else []
                        _ws = _kperms.workspace if _kperms else None
                    except Exception:
                        _extra = []
                        _ws = None
                    _proj_mem = load_project_memory(
                        cwd=_ws, extra_roots=_extra,
                    )
                    if _proj_mem:
                        _memory = (_memory + "\n\n" + _proj_mem).strip() if _memory else _proj_mem

                    # Provider profile is injected by PromptLoader.
                    # Keep memory_context reserved for session memory + transient state
                    # so we do not pay twice for the same profile tokens.
                    # Phase 5f: note pending file uploads so the agent
                    # knows where to find them. Cleared after dispatch.
                    _pending_imgs = state.get("_pending_images") or []
                    if _pending_imgs:
                        _img_lines = "\n".join(
                            f"  - {p}" for p in _pending_imgs
                        )
                        current_msg = (
                            f"{current_msg}\n\n"
                            f"[The user attached these files — "
                            f"read them with read_file (for text/code/configs) "
                            f"or notebook_read (for .ipynb). Images can be "
                            f"described with whatever vision tools you have:\n"
                            f"{_img_lines}]"
                        )
                        state["_pending_images"] = []
                        try:
                            image_upload.value = ()
                        except Exception:
                            pass
                    _denied = state.get("_denied_commands", [])
                    if _denied:
                        _denial_ctx = (
                            "\n[System] Previously BLOCKED commands in this session: "
                            + ", ".join(_denied[-5:])
                            + ". Do NOT retry these."
                        )
                        current_msg = current_msg + _denial_ctx

                    engine.stream_response(
                        user_message=current_msg,
                        on_token=_on_token,
                        on_tool_use=_on_tool_use,
                        on_tool_result=_on_tool_result,
                        on_permission_denied=_on_permission_denied,
                        on_thinking=_on_thinking,
                        thinking_budget=_budget,
                        memory_context=_memory,
                    )
                    # Final update: finalize=True triggers full markdown rendering
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label, finalize=True)

                    # Auto-execute slash commands from agent output (all modes).
                    # Dashboard, Solo, Builder — any agent can control the UI
                    # via ACTION: /command lines. Safety tiers still enforced.
                    # Cap at 3 continuation rounds (down from 4): empirical
                    # data on multi-step dashboard requests shows >3 rounds
                    # is always model-confusion, not user intent. The
                    # explicit `ACTION: /done` sentinel below skips even
                    # the post-execute commentary turn so a 1- or 2-action
                    # sequence doesn't burn an extra 30-120 s + $0.02-0.05
                    # on the model emitting "(commands executed)".
                    _MAX_ACTION_CONT = 3
                    _cont_turn = 0
                    while chunks and _cont_turn < _MAX_ACTION_CONT:
                        _cont_turn += 1
                        raw = "".join(chunks)
                        exec_results = _dashboard_auto_exec(raw)
                        # Split the done-sentinel from real results so the
                        # placeholder accurately reflects what happened.
                        done_seen = "__DONE__" in exec_results
                        real_results = [r for r in exec_results if r != "__DONE__"]
                        # Remove ACTION: lines from visible output
                        import re as _re_strip
                        cleaned = _re_strip.sub(
                            r"^ACTION:\s*/.*$", "", raw, flags=_re_strip.MULTILINE
                        ).strip()
                        if cleaned != raw.strip():
                            if cleaned:
                                placeholder = cleaned
                            elif real_results:
                                placeholder = "(commands executed)"
                            elif done_seen:
                                # /done alone — agent emitted no real
                                # action. Make the empty bubble explicit
                                # so the user doesn't think something
                                # ran silently.
                                placeholder = (
                                    "(agent had no action to execute — "
                                    "please clarify or rephrase)"
                                )
                            else:
                                placeholder = "(commands executed)"
                            _update_last_assistant(placeholder, role_label)
                        # Explicit done-sentinel from the agent — skip
                        # the wrap-up continuation entirely. Triggers
                        # whether /done was emitted alone or alongside
                        # real ACTIONs.
                        if done_seen:
                            break
                        # No results → nothing to continue on
                        if not real_results:
                            break
                        # Re-bind exec_results to the real slice so the
                        # feedback that goes back to the model doesn't
                        # contain the sentinel.
                        exec_results = real_results
                        # Inject results and run the agent again
                        feedback = "[Command results]\n" + "\n".join(exec_results)
                        engine.messages.append(
                            {"role": "user", "content": feedback}
                        )
                        if engine._stop_requested:
                            break
                        # Continuation turn
                        chunks.clear()
                        thinking_chunks.clear()
                        engine.stream_response(
                            user_message=feedback,
                            on_token=_on_token,
                            on_tool_use=_on_tool_use,
                            on_tool_result=_on_tool_result,
                            on_permission_denied=_on_permission_denied,
                            on_thinking=_on_thinking,
                            thinking_budget=_budget,
                            memory_context=_memory,
                        )
                        if chunks:
                            _update_last_assistant("".join(chunks), role_label, finalize=True)
                    # Cap-hit notification — only when the loop actually
                    # exhausted all rounds without finishing on its own.
                    if _cont_turn >= _MAX_ACTION_CONT and chunks:
                        _post_raw = "".join(chunks)
                        if _extract_action_commands(_post_raw):
                            _append_system_message(
                                f"⏸ Stopped after {_MAX_ACTION_CONT} ACTION rounds — "
                                f"agent kept emitting commands. Send a follow-up to continue."
                            )

                    # -- Interactive question detection (solo/dashboard) --
                    # After the agent finishes a turn, check if the response
                    # ends with a question and show option buttons if applicable.
                    _hide_question_ui()  # always reset first
                    if chunks and engine.mode in ("solo", "dashboard"):
                        _full_text = "".join(chunks)
                        _q_info = _detect_question(_full_text)
                        if _q_info:
                            _show_question_ui(_q_info)

                    # B1 — Outcome tracking for ALL conversational modes,
                    # not just solo. Pipeline modes still go through the
                    # cycle gate's record_cycle_outcome; this captures the
                    # interactive dashboard / solo / quick / reviewed / etc.
                    # turns the user actually drives.
                    if chunks and engine.mode in (
                        "solo", "dashboard", "quick", "reviewed",
                        "tdd", "cluster", "full",
                    ):
                        _record_turn_outcome(
                            engine,
                            user_task=original_task,
                            response_text="".join(chunks),
                            state=state,
                            start_time=_turn_start_time,
                        )

                    # S7 — Show per-turn cost in ALL modes including solo +
                    # dashboard. Pipeline rows stay verbose (they include the
                    # role label so the user can attribute cost per agent);
                    # solo/dashboard get a single compact line.
                    _role_cost = engine.cost_usd - _cost_before
                    _role_in = engine.token_usage["input"] - _in_before
                    _role_out = engine.token_usage["output"] - _out_before
                    if _role_in > 0 or _role_out > 0:
                        _cost_str = f"${_role_cost:.3f}" if _role_cost > 0 else ""
                        _is_pipeline_mode = engine.mode not in ("solo", "dashboard")
                        if _is_pipeline_mode:
                            _append_system_message(
                                f"{role_label}: {_role_in:,} in / {_role_out:,} out"
                                f"{' · ' + _cost_str if _cost_str else ''}"
                                f" [{_effective_model}]"
                            )
                        elif _cost_str:
                            _append_system_message(
                                f"Turn cost: {_cost_str} "
                                f"({_role_in:,} in / {_role_out:,} out · {_effective_model})"
                            )

                    if engine._stop_requested:
                        break

                    # --- Cost Governor ---
                    # No hard stop — pausing + restarting wastes more money
                    # than letting the pipeline finish. Show milestones so the
                    # user can see how much is being spent.
                    _cost_milestones = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
                    _passed = state.setdefault("_cost_milestones_passed", set())
                    for _ms in _cost_milestones:
                        if engine.cost_usd >= _ms and _ms not in _passed:
                            _passed.add(_ms)
                            _append_system_message(
                                f"Cost milestone: ${engine.cost_usd:.2f} (passed ${_ms:.0f})"
                            )

                    # S7 — Soft-limit banner: show ONCE per session when the
                    # configured threshold is crossed. Doesn't stop anything;
                    # the user decides whether to /stop or continue.
                    if not state.get("_cost_soft_limit_warned"):
                        try:
                            _soft = float(
                                _get_agent_settings().get("cost_soft_limit_usd", 5.0)
                            )
                        except (TypeError, ValueError):
                            _soft = 5.0
                        if _soft > 0 and engine.cost_usd >= _soft:
                            state["_cost_soft_limit_warned"] = True
                            _append_system_message(
                                f"💸 Cost soft-limit reached: ${engine.cost_usd:.2f} "
                                f"≥ ${_soft:.2f}. Continue or /stop — limit "
                                f"configurable via agent.cost_soft_limit_usd."
                            )

                    # --- Auto-advance logic ---
                    prev_role_id = engine.current_role

                    # --- Get last assistant output for all gate checks ---
                    last_out = ""
                    for msg in reversed(engine.messages):
                        if msg["role"] == "assistant":
                            last_out = msg["content"]
                            break

                    # --- Structured output validation gate ---
                    _schema_errors = engine.validate_role_output(prev_role_id, last_out)
                    if _schema_errors:
                        _schema_retry_counts = state.setdefault("_schema_retry_counts", {})
                        _schema_retry = int(_schema_retry_counts.get(prev_role_id, 0))
                        if _schema_retry < 2:
                            _schema_retry_counts[prev_role_id] = _schema_retry + 1
                            _err_lines = "\n".join(f"- {err}" for err in _schema_errors)
                            _append_gate_message(
                                "schema",
                                prev_role_id,
                                f"{_format_role_label(prev_role_id)} output invalid",
                                _err_lines,
                                f"Automatic schema retry {_schema_retry + 1}/2 is being attempted.",
                            )
                            current_msg = (
                                "Your previous output did not match the required structured format.\n"
                                "Fix the schema violations below and resend the full output in the "
                                "mandatory format for your role.\n\n"
                                f"Schema violations:\n{_err_lines}\n\n"
                                f"Original task:\n{original_task}"
                            )
                            _update_status()
                            _update_pipeline_display(engine)
                            continue
                        _append_gate_message(
                            "schema",
                            prev_role_id,
                            f"{_format_role_label(prev_role_id)} output still invalid",
                            "\n".join(f"- {err}" for err in _schema_errors),
                            "Give guidance or let the same role try again.",
                        )
                        state["_awaiting_agent_question"] = prev_role_id
                        _update_status()
                        _update_pipeline_display(engine)
                        break
                    else:
                        _schema_retry_counts = state.setdefault("_schema_retry_counts", {})
                        _schema_retry_counts.pop(prev_role_id, None)

                    # --- Communication gate: explicit review for risky/partial handoffs ---
                    _gate_action, _gate_type, _gate_message = engine.evaluate_role_gate(
                        prev_role_id, last_out
                    )
                    if _gate_action == "pause":
                        _append_gate_message(
                            _gate_type or "review",
                            prev_role_id,
                            f"{_format_role_label(prev_role_id)} gate triggered",
                            _gate_message,
                            "Reply 'go' to continue anyway, or give guidance.",
                        )
                        state["_awaiting_gate_review"] = prev_role_id
                        _update_status()
                        _update_pipeline_display(engine)
                        break

                    if "SKIP" in last_out[:200].upper() and prev_role_id not in (
                        "session_manager", "builder_agent", "test_agent",
                    ):
                        _append_system_message(
                            f"--- {_format_role_label(prev_role_id)}: "
                            f"nothing to do, skipping ---"
                        )
                        engine.advance_role()
                        _update_status()
                        _update_pipeline_display(engine)
                        if engine.is_cycle_complete:
                            _append_system_message("--- Cycle complete ---")
                            break
                        # Build handoff for next role
                        current_msg = engine.build_handoff_message(original_task)
                        continue

                    # --- QUESTION: Protocol — any agent can ask the user ---
                    if "QUESTION:" in last_out:
                        import re as _re_q
                        _q_match = _re_q.search(
                            r"QUESTION:\s*(.+?)(?:\n\n|\Z)", last_out, _re_q.DOTALL
                        )
                        if _q_match:
                            _q_text = _q_match.group(1).strip()
                            _role_label = _format_role_label(prev_role_id)
                            _append_gate_message(
                                "question",
                                prev_role_id,
                                f"{_role_label} fragt dich:",
                                _q_text,
                                f"Antworte dem {_role_label}.",
                            )
                            state["_awaiting_agent_question"] = prev_role_id
                            _update_status()
                            _update_pipeline_display(engine)
                            break

                    # --- Confidence Gate — pause on low confidence ---
                    if "**confidence:** low" in last_out.lower():
                        _append_gate_message(
                            "confidence",
                            prev_role_id,
                            f"{_format_role_label(prev_role_id)} has low confidence",
                            "The agent reported low confidence for this step.",
                            "Review the output and reply 'go' to continue or give guidance.",
                        )
                        state["_awaiting_agent_question"] = prev_role_id
                        _update_status()
                        _update_pipeline_display(engine)
                        break

                    # Session Manager: STOP and wait for user approval.
                    # If the SM responded conversationally (greeting, clarification),
                    # just wait for the next user message — no gate needed.
                    if prev_role_id == "session_manager" and not _sm_approval:
                        if _AE.is_conversational(prev_role_id, last_out):
                            # Conversational response — stay in chat mode,
                            # next user message will loop back here.
                            _update_status()
                            _update_pipeline_display(engine)
                            break

                        # Always show plan for user approval — no auto-approve.
                        # The user should always be involved in the planning process.
                        _append_gate_message(
                            "plan-approval",
                            prev_role_id,
                            "Plan approval required",
                            (
                                "The Session Manager produced a plan. "
                                "Review it and approve via the Continue button "
                                "or type corrections."
                            ),
                            "",
                        )
                        _update_status()
                        _update_pipeline_display(engine)
                        break

                    # Dynamic routing: parse SM's routing directives.
                    if prev_role_id == "session_manager" and _sm_approval:
                        # --- Skip agents ---
                        _SKIPPABLE = {"critic_agent", "reviewer_agent", "runtime_agent", "research_agent"}
                        skip_section = re.search(
                            r"### Skip agents\s*\n(.*?)(?:\n###|\n\*\*|$)",
                            last_out, re.DOTALL
                        )
                        if skip_section:
                            skip_text = skip_section.group(1)
                            skipped = []
                            for role_name in _SKIPPABLE:
                                if role_name in skip_text:
                                    skipped.append(role_name)
                            if skipped:
                                engine.route = [
                                    r for r in engine.route if r not in skipped
                                ]
                                labels = [_format_role_label(s) for s in skipped]
                                _append_system_message(
                                    f"--- Dynamic routing: skipping "
                                    f"{', '.join(labels)} (SM recommendation) ---"
                                )

                        # --- SKIP_RESEARCH shorthand ---
                        if "SKIP_RESEARCH" in last_out and "research_agent" in engine.route:
                            engine.route = [r for r in engine.route if r != "research_agent"]
                            _append_system_message(
                                "--- Dynamic routing: skipping Research Agent "
                                "(no external info needed) ---"
                            )

                        # --- Add agents (NEW) ---
                        _ADDABLE = {"research_agent", "runtime_agent", "critic_agent", "reviewer_agent"}
                        add_section = re.search(
                            r"### Add agents\s*\n(.*?)(?:\n###|\n\*\*|$)",
                            last_out, re.DOTALL
                        )
                        if add_section:
                            add_text = add_section.group(1)
                            added = []
                            for role_name in _ADDABLE:
                                if role_name in add_text and role_name not in engine.route:
                                    # Insert before builder_agent
                                    _b_idx = (engine.route.index("builder_agent")
                                              if "builder_agent" in engine.route
                                              else len(engine.route))
                                    engine.route.insert(_b_idx, role_name)
                                    added.append(role_name)
                            if added:
                                labels = [_format_role_label(a) for a in added]
                                _append_system_message(
                                    f"--- Dynamic routing: adding "
                                    f"{', '.join(labels)} (SM recommendation) ---"
                                )

                        # --- RESEARCH_NEEDED shorthand ---
                        if "RESEARCH_NEEDED:" in last_out and "research_agent" not in engine.route:
                            _b_idx = (engine.route.index("builder_agent")
                                      if "builder_agent" in engine.route
                                      else len(engine.route))
                            engine.route.insert(_b_idx, "research_agent")
                            _append_system_message(
                                "--- Dynamic routing: adding Research Agent "
                                "(SM requested research) ---"
                            )

                    # --- Critic/Runtime findings gate ---
                    # Only pause on REJECT — approve/approve_with_risks auto-continues.
                    # This avoids wasting tokens by pausing on every review.
                    _REVIEW_ROLES = {"critic_agent", "runtime_agent"}
                    if (prev_role_id in _REVIEW_ROLES
                            and not state.get("_findings_approved")
                            and last_out.strip()):
                        _lower500 = last_out[:500].lower()
                        _is_reject = ("reject" in _lower500
                                      and "approve" not in _lower500)
                        if _is_reject:
                            _append_gate_message(
                                "findings",
                                prev_role_id,
                                f"{_format_role_label(prev_role_id)} rejected",
                                f"The review found critical issues. Review and decide.",
                                "Reply 'go' to continue anyway, or /reset to stop.",
                            )
                            state["_awaiting_findings_review"] = prev_role_id
                            _update_status()
                            _update_pipeline_display(engine)
                            break

                    # --- Reviewer → Builder loop (max retries) ---
                    if prev_role_id == "reviewer_agent":
                        has_issues = "ISSUES" in last_out[:500].upper()
                        _retries = state.get("_builder_retries", 0)
                        if has_issues and _retries < 3:
                            state["_builder_retries"] = _retries + 1
                            # Extract specific findings for Builder
                            _findings = _extract_retry_context(
                                last_out, "reviewer"
                            )
                            engine.retry_from_builder()
                            _record_cycle_event(
                                "retry",
                                f"Builder retry {_retries + 1}/3 from Reviewer",
                                _compact_text(_findings, 160),
                                "reviewer_agent",
                            )
                            _append_system_message(
                                f"--- Reviewer found issues. "
                                f"Builder retry {_retries + 1}/3 ---"
                            )
                            engine.compact_for_next_role()
                            current_msg = (
                                f"RETRY from Reviewer (attempt {_retries + 1}/3).\n"
                                f"Fix these specific issues:\n{_findings}\n\n"
                                + engine.build_handoff_message(original_task)
                            )
                            _update_status()
                            _update_pipeline_display(engine)
                            continue

                    # --- Test → Builder loop (max retries shared) ---
                    if prev_role_id == "test_agent":
                        _fail_kw = ("FAIL", "FAILED", "ERROR", "error:",
                                    "Exception", "Traceback", "test failed")
                        has_fail = any(k.lower() in last_out.lower() for k in _fail_kw)
                        _retries = state.get("_builder_retries", 0)
                        if has_fail and _retries < 3:
                            state["_builder_retries"] = _retries + 1
                            # Extract specific test failures for Builder
                            _findings = _extract_retry_context(
                                last_out, "test"
                            )
                            engine.retry_from_builder()
                            _record_cycle_event(
                                "retry",
                                f"Builder retry {_retries + 1}/3 from Test",
                                _compact_text(_findings, 160),
                                "test_agent",
                            )
                            _append_system_message(
                                f"--- Test failures detected. "
                                f"Builder retry {_retries + 1}/3 ---"
                            )
                            engine.compact_for_next_role()
                            current_msg = (
                                f"RETRY from Test Agent (attempt {_retries + 1}/3).\n"
                                f"Fix these specific failures:\n{_findings}\n\n"
                                + engine.build_handoff_message(original_task)
                            )
                            _update_status()
                            _update_pipeline_display(engine)
                            continue

                    # --- Conflict Resolution: Critic rejected but Builder approved ---
                    if prev_role_id == "builder_agent":
                        _critic_out = engine.role_outputs.get("critic_agent", "")
                        if ("reject" in _critic_out[:500].lower()
                                and "approve" in last_out[:500].lower()
                                and not state.get("_conflict_resolved")):
                            _append_gate_message(
                                "conflict",
                                prev_role_id,
                                "Critic / Builder conflict",
                                "Critic rejected the plan but Builder approved the implementation.",
                                "Reply 'go' to continue testing, or 'reject' to stop.",
                            )
                            state["_awaiting_conflict_resolution"] = True
                            _update_status()
                            _update_pipeline_display(engine)
                            break
                        state["_conflict_resolved"] = True

                    # Dashboard / Solo / follow-up: never advance — keep conversation open
                    # No pipeline messages for conversational modes.
                    if mode_dropdown.value in ("dashboard", "solo") or state.get("_follow_up"):
                        _set_active_gate()
                        _update_status()
                        break

                    # Advance to next role
                    has_next = engine.advance_role()
                    _set_active_gate()
                    if not has_next:
                        state.pop("_retry_used", None)
                        state.pop("_builder_retries", None)
                        state.pop("_conflict_resolved", None)
                        _record_cycle_event("cycle", "Cycle complete")
                        # Acceptance gate: check if test agent approved
                        _cycle_verdict = _check_acceptance_gate(engine)
                        _append_system_message(
                            f"--- Cycle complete {_cycle_verdict} ---"
                        )
                        _update_pipeline_display(engine)

                        # --- Persistent Cycle Memory + Provider Profile ---
                        try:
                            import json as _json_mem
                            from datetime import datetime as _dt_mem
                            _mem_path = ctx.agent_dir / ".cycle_memory.jsonl"
                            _cycle_mem = {
                                "task": original_task[:200],
                                "mode": mode_dropdown.value,
                                "route": engine.route,
                                "verdict": _cycle_verdict,
                                "retries": state.get("_builder_retries", 0),
                                "cost_usd": round(engine.cost_usd, 4),
                                "timestamp": _dt_mem.now().isoformat(),
                            }
                            with open(_mem_path, "a", encoding="utf-8") as _mf:
                                _mf.write(_json_mem.dumps(_cycle_mem) + "\n")
                        except Exception:
                            pass
                        # Self-optimization: record outcome + update provider profile
                        try:
                            _denied = state.get("_denied_commands", [])
                            _opt_changes = engine.record_cycle_outcome(
                                verdict=_cycle_verdict,
                                user_task=original_task,
                                denied_commands=_denied,
                                start_time=state.get("session_start_time"),
                            )
                            if _opt_changes:
                                _opt_str = "; ".join(
                                    f"{k}: {v}" for k, v in _opt_changes.items()
                                )
                                _append_system_message(
                                    f"\U0001f504 Self-optimization: {_opt_str}"
                                )
                        except Exception:
                            pass

                        # Enter follow-up mode: route to builder for continued work
                        _fu_role = "builder_agent"
                        if _fu_role not in engine.route:
                            _fu_role = engine.route[-1]
                        _fu_idx = engine.route.index(_fu_role) if _fu_role in engine.route else 0
                        engine.current_role_index = _fu_idx
                        state["_follow_up"] = True
                        _append_system_message(
                            f"Continue chatting or /reset for new task."
                        )
                        break

                    # Show pipeline progress + handoff in chat
                    next_role = _format_role_label(engine.current_role)
                    _record_cycle_event(
                        "handoff",
                        f"{_format_role_label(prev_role_id)} -> {next_role}",
                        role=prev_role_id,
                    )
                    _append_system_message(
                        f"--- Auto-handoff: "
                        f"{_format_role_label(prev_role_id)} \u2192 {next_role} ---"
                    )

                    # Independent analysis: if two consecutive review roles
                    # (e.g., Critic → Runtime), the second should NOT see
                    # the first's output — they analyze independently.
                    _ANALYSIS_ROLES = {"critic_agent", "runtime_agent"}
                    next_id = engine.current_role
                    if (prev_role_id in _ANALYSIS_ROLES
                            and next_id in _ANALYSIS_ROLES):
                        # Temporarily hide previous reviewer's output
                        _saved_output = engine.role_outputs.pop(prev_role_id, None)
                        _append_system_message(
                            f"--- Independent analysis: {next_role} "
                            f"reviews without seeing {_format_role_label(prev_role_id)}'s output ---"
                        )

                    # Compact context between roles to save tokens
                    engine.compact_for_next_role()

                    _update_status()
                    _update_pipeline_display(engine)

                    # Build context-rich handoff message for next agent
                    current_msg = engine.build_handoff_message(original_task)

                    # Restore hidden output after handoff message is built
                    if (prev_role_id in _ANALYSIS_ROLES
                            and next_id in _ANALYSIS_ROLES
                            and _saved_output is not None):
                        engine.role_outputs[prev_role_id] = _saved_output

            except Exception as exc:
                error_text = str(exc)
                is_crash = (
                    "CLI error" in error_text
                    or "Not logged in" in error_text
                    or "Broken pipe" in error_text
                )
                if is_crash:
                    # Auto-recover: invalidate engine so next send recreates it
                    state["engine"] = None
                    _append_system_message(
                        f"CLI process crashed: {error_text[:200]}\n"
                        f"Engine will auto-restart on next message."
                    )
                elif chunks:
                    _update_last_assistant(
                        "".join(chunks) + f"\n\nError: {error_text}", role_label
                    )
                else:
                    _append_system_message(f"Error: {error_text}")
            finally:
                # Only clean up UI state if no newer generation has started
                # (prevents stale worker from clobbering a fresh send after Stop)
                with state["_state_lock"]:
                    _is_current = state.get("_generation_id") == _my_gen_id
                    if _is_current:
                        state["streaming"] = False
                # Disarm the stale + kill watchers — the worker is done
                # one way or another, no need to flag stale-ness or send
                # a cooperative stop on an already-closed turn.
                try:
                    for _slot in ("_stale_timer", "_stale_kill_timer"):
                        _t = state.get(_slot)
                        if _t is not None:
                            try:
                                _t.cancel()
                            except Exception:
                                pass
                        state[_slot] = None
                    state["_stale_seen"] = False
                except Exception:
                    pass
                if _is_current:
                    # Silent-exit detection: if the worker finished without
                    # producing any tokens AND no exception fired, surface a
                    # system note so the user isn't left wondering why the
                    # spinner just disappeared. This is the "agent stops
                    # mid-task" failure mode reported in production.
                    if not chunks and not thinking_chunks:
                        _append_system_message(
                            "Agent returned no output. The CLI ended the turn "
                            "without producing text — try /retry or /status."
                        )
                    _set_working(False)
                    _update_status()
                    _update_button_states()
                    _auto_save_session()
                    _check_auto_compact()
                    # Stalled-task detector — once per turn check
                    # whether any in-progress task has gone untouched
                    # for too long and suggest /mode plan.  Best-effort,
                    # never raises.
                    try:
                        kp = getattr(engine, "kit_permissions", None)
                        ws = (kp.workspace if kp is not None
                              else (ctx.repo_dir or Path.cwd()))
                        from delfin.agent.agent_tasks import get_store
                        store = get_store(Path(ws))
                        stalled = store.find_stalled(
                            max_age_s=600.0,
                            session_id=getattr(engine, "session_id", "") or None,
                        )
                        if stalled and not state.get("_stalled_warning_shown"):
                            top = stalled[0]
                            _append_system_message(
                                f"⚠️ Task #{top['id']} '{top.get('subject','')[:60]}' "
                                f"has been in_progress for {top['age_s']//60} min "
                                "without an update. Consider /mode plan to "
                                "re-think the approach, or /forget it."
                            )
                            state["_stalled_warning_shown"] = True
                        elif not stalled:
                            state["_stalled_warning_shown"] = False
                    except Exception:
                        pass
                    # Live per-turn cost footer: post a one-line system
                    # message summarising the delta in tokens/cost/
                    # duration for the turn we just finished. Best-effort,
                    # silently skips if get_status() raises.
                    try:
                        _post = engine.get_status()
                        d_in = int(_post.get("input_tokens") or 0) - int(state.get("_turn_pre_in") or 0)
                        d_out = int(_post.get("output_tokens") or 0) - int(state.get("_turn_pre_out") or 0)
                        d_cost = float(_post.get("cost_usd") or 0.0) - float(state.get("_turn_pre_cost") or 0.0)
                        dur = time.monotonic() - float(state.get("_turn_started_monotonic") or time.monotonic())
                        if d_in or d_out or d_cost > 0.0001 or tool_count[0]:
                            cost_s = f"${d_cost:.4f}" if d_cost > 0.0001 else "<$0.0001"
                            _append_system_message(
                                f"⏱ turn  {dur:.1f}s  ·  "
                                f"{d_in:,}↓ / {d_out:,}↑ tokens  ·  "
                                f"{tool_count[0]} tool call{'s' if tool_count[0] != 1 else ''}  ·  "
                                f"{cost_s}"
                            )
                    except Exception:
                        pass
                    # Fire Stop hooks at end-of-turn so user-defined linters
                    # / test runs / notifications can react. Surface noisy
                    # output (non-zero exit or stderr) as a system note so
                    # the user sees hook activity instead of it being silent.
                    try:
                        from delfin.agent import hooks as _hooks_mod
                        _cfg = _hooks_mod.load_hooks(ctx.repo_dir or None)
                        if not _cfg.is_empty():
                            _stops = _hooks_mod.run_hooks(
                                "Stop", _cfg,
                                workspace=ctx.repo_dir or None,
                            )
                            _noisy = [r for r in _stops
                                      if r.matched and (r.stderr.strip() or r.exit_code)]
                            if _noisy:
                                _lines = ["🎣 Stop hooks fired:"]
                                for r in _noisy:
                                    _lines.append(
                                        f"  {r.command[:60]}  exit={r.exit_code}  "
                                        f"{r.duration_s:.1f}s"
                                    )
                                    if r.stderr.strip():
                                        _lines.append(f"    stderr: {r.stderr.strip()[:120]}")
                                _append_system_message("\n".join(_lines))
                    except Exception:
                        pass
                    # Process next queued message if any
                    _process_queue()
                else:
                    # Stale worker — just save session, don't touch UI
                    _auto_save_session()

        threading.Thread(target=_worker, daemon=True).start()

    def _process_queue():
        """Send the next queued message, if any."""
        if state["message_queue"] and not state["streaming"]:
            next_msg = state["message_queue"].pop(0)
            _update_queue_display()
            input_textarea.value = next_msg
            _on_send(None)

    def _on_stop(button):
        # Soft-stop: cooperative end of the running turn so the next Send
        # can continue the same conversation (CLI session_id preserved).
        # Queue is intentionally NOT cleared — user may want to keep the
        # follow-up prompts. Old behaviour killed the subprocess and wiped
        # the queue, which forced a full re-init on every ESC.
        engine = state["engine"]
        if engine:
            engine.request_stop()
            if hasattr(engine.client, "signal_stop"):
                engine.client.signal_stop()
        # Bump generation so the old worker's finally block won't touch UI
        state["_generation_id"] = state.get("_generation_id", 0) + 1
        state["streaming"] = False
        # Finalize any in-progress streaming message
        msgs = state["chat_messages"]
        if msgs and msgs[-1].get("_streaming"):
            msgs[-1]["_streaming"] = False
            _refresh_chat_html()
        _set_active_gate()
        # If the user stopped while messages are queued, show a queued
        # spinner so they know the agent is idle-but-pending, not dead.
        queued = len(state.get("message_queue") or [])
        if queued:
            _set_working(True, f"Idle - {queued} queued", mode="queued")
        else:
            _set_working(False)
        _update_button_states()
        _record_cycle_event("pause", "Generation stopped by user")
        if queued:
            _append_system_message(
                f"Stop. {queued} queued message(s) preserved — Send to continue."
            )
        else:
            _append_system_message("Stop. Send a message to continue.")

    def _on_new_cycle(button):
        # Save current session before clearing so it can be restored later
        if state["chat_messages"]:
            _auto_save_session()
            _refresh_session_dropdown()
        engine = state["engine"]
        if engine:
            engine.reset_cycle(mode=mode_dropdown.value)
        else:
            state["engine"] = None
        state["chat_messages"].clear()
        state["streaming"] = False
        state["active_session_id"] = ""
        try:
            kp = getattr(engine, "kit_permissions", None) if engine else None
            if kp is not None:
                kp.task_session_id = ""
        except Exception:
            pass
        state["recent_edits"].clear()
        state.pop("_mode_suggested", None)
        state.pop("_pending_mode_msg", None)
        state.pop("_retry_used", None)
        state.pop("_builder_retries", None)
        state.pop("_follow_up", None)
        state["_cycle_history"] = []
        state["_mode_manual_override"] = False
        _set_active_gate()
        state["message_queue"].clear()
        state["session_start_time"] = None
        queue_html.value = ""
        undo_btn.disabled = True
        session_dropdown.value = ""
        _refresh_chat_html()
        _refresh_task_ticker()
        _update_status()
        _update_button_states()

    def _on_advance_role(button):
        """Manual advance — fallback when auto-advance was stopped."""
        engine = state["engine"]
        if not engine:
            return
        prev_role = _format_role_label(engine.current_role)
        has_next = engine.advance_role()
        _set_active_gate()
        if has_next:
            next_role = _format_role_label(engine.current_role)
            _record_cycle_event("handoff", f"{prev_role} -> {next_role}")
            _append_system_message(
                f"--- Manual handoff: {prev_role} \u2192 {next_role} ---"
            )
        else:
            state.pop("_retry_used", None)
            _record_cycle_event("cycle", "Cycle complete")
            _append_system_message("--- Cycle complete ---")
        _update_status()
        _update_button_states()
        _auto_save_session()

    def _on_mode_change(change):
        new_mode = change["new"]
        old_mode = (change.get("old") or "")
        if not state.get("_mode_change_internal"):
            state["_mode_manual_override"] = True
        # Update mode description label
        desc = _MODE_DESCRIPTIONS.get(new_mode, "")
        mode_desc_html.value = (
            f'<span style="color:#888; font-size:12px; margin-left:4px;">'
            f'{desc}</span>'
        )
        engine = state["engine"]
        if engine and not state["streaming"]:
            if not engine.messages:
                engine.reset_cycle(mode=new_mode)
                _update_status()
            elif getattr(engine, "mode", "") != new_mode:
                # Mode-switch mid-conversation: preserve message history
                # so the new agent sees the user's prompt and doesn't
                # have to ask for a re-paste. Role-tracking is reset
                # because the new mode has a different role route, but
                # the conversation thread is kept intact. Save the old
                # session first so the original mode's run remains
                # browsable in the session list.
                if state.get("chat_messages"):
                    _auto_save_session()
                    _refresh_session_dropdown()
                engine.reset_cycle(mode=new_mode, preserve_messages=True)
                state.pop("_follow_up", None)
                state["_cycle_history"] = []
                _set_active_gate()
                _refresh_task_ticker()
                _append_system_message(
                    f"Mode switched to {new_mode}. Continuing with the "
                    "existing conversation — no need to re-send your "
                    "prompt."
                )
                _update_status()
        # Dashboard mode: lock permission to default, model to cheapest
        if new_mode == "dashboard":
            state["_perm_before_dashboard"] = perm_dropdown.value
            perm_dropdown.value = "ask_all"
            perm_dropdown.disabled = True
            state["_model_before_dashboard"] = model_dropdown.value
            _cheap = _PROVIDER_CHEAP.get(provider_dropdown.value, "haiku")
            model_dropdown.value = _cheap
        elif new_mode == "plan":
            # Plan mode locks the permission profile to "plan" (read-only)
            # so no Edit/Write/Bash can fire until the user clicks "Plan
            # akzeptieren" (which flips back to the saved profile).
            state["_perm_before_plan"] = perm_dropdown.value
            perm_dropdown.value = "plan"
            perm_dropdown.disabled = True
            _append_system_message(
                "🧭 Plan mode active — read-only. The agent will draft a "
                "plan and call ExitPlanMode for your approval before any "
                "edits or bash runs."
            )
        else:
            perm_dropdown.disabled = state.get("streaming", False)
            saved = state.pop("_perm_before_dashboard", None)
            if saved and perm_dropdown.value == "ask_all":
                perm_dropdown.value = saved
            saved_plan = state.pop("_perm_before_plan", None)
            if saved_plan and perm_dropdown.value == "plan":
                perm_dropdown.value = saved_plan
            saved_model = state.pop("_model_before_dashboard", None)
            _cheap = _PROVIDER_CHEAP.get(provider_dropdown.value, "haiku")
            if saved_model and model_dropdown.value == _cheap:
                model_dropdown.value = saved_model
        # Show/hide Cycle Inspector based on mode
        _update_cycle_inspector()

        # Solo-minimal UI: hide pipeline overhead for solo/dashboard
        _is_minimal = new_mode in ("solo", "dashboard", "plan")
        # Hide mode description (saves vertical space)
        mode_desc_html.layout.display = "none" if _is_minimal else "block"
        # Hide pipeline-only buttons (but keep commit/push visible)
        advance_btn.layout.display = "none" if _is_minimal else "inline-flex"

    def _on_provider_change(change):
        """Switch provider (Claude / OpenAI / KIT), update model options."""
        if state["streaming"]:
            return
        provider = change["new"]
        # Try fetching models dynamically from API
        fetched = _fetch_models(provider)
        if fetched:
            _PROVIDER_MODELS[provider] = fetched
        models = _PROVIDER_MODELS.get(provider, _PROVIDER_MODELS_FALLBACK.get(
            provider, _PROVIDER_MODELS_FALLBACK["claude"]))
        model_dropdown.options = models
        default = _PROVIDER_DEFAULTS.get(provider, models[0][1])
        valid_values = {v for _, v in models}
        model_dropdown.value = default if default in valid_values else models[0][1]
        # Invalidate engine
        engine = state["engine"]
        if engine:
            if hasattr(engine.client, "kill"):
                engine.client.kill()
            state["engine"] = None
        # Show/hide the KIT confirmation panel based on the new provider.
        try:
            _show_kit_confirm_panel(provider == "kit")
        except Exception:
            pass
        n_models = len(models)
        src = "live" if fetched else "fallback"
        _append_system_message(
            f"Provider switched to {provider}. {n_models} models loaded ({src})."
        )
        _update_status()
        # Persist
        try:
            from delfin.user_settings import load_settings, save_settings
            s = load_settings()
            s.setdefault("agent", {})
            s["agent"]["provider"] = provider
            save_settings(s)
        except Exception:
            pass

    def _on_model_change(change):
        """Recreate engine with new model on next send."""
        if state["streaming"]:
            return
        engine = state["engine"]
        if engine:
            # Force engine recreation with new model
            state["engine"] = None
            _append_system_message(
                f"Model switched to {change['new']}. Next message uses new model."
            )
        # Persist the choice
        try:
            from delfin.user_settings import load_settings, save_settings
            s = load_settings()
            s.setdefault("agent", {})
            s["agent"]["model"] = change["new"]
            save_settings(s)
        except Exception:
            pass

    def _on_effort_change(change):
        """Persist effort preference."""
        if state["streaming"]:
            return
        try:
            from delfin.user_settings import load_settings, save_settings
            s = load_settings()
            s.setdefault("agent", {})
            s["agent"]["effort"] = change["new"]
            save_settings(s)
        except Exception:
            pass

    def _on_perm_change(change):
        """Sync permission profile from dropdown to state, recreate engine."""
        if state["streaming"]:
            return
        new_profile = change["new"]
        state["_perm_profile"] = new_profile
        # If the change came from the KIT-Mode chip, the live perms.mode
        # already reflects it and we don't need to recreate the engine.
        # Just sync silently and return.
        if state.get("_chip_syncing_perm"):
            return
        engine = state["engine"]
        if engine:
            # Apply to the live KIT engine immediately if we can — avoids
            # the "next message" lag and keeps the chip + dropdown in sync.
            chip_target = _PROFILE_TO_CHIP.get(new_profile)
            if chip_target and hasattr(engine, "set_kit_permission_mode"):
                if engine.set_kit_permission_mode(chip_target):
                    _refresh_kit_mode_chip()
                    _append_system_message(
                        f"Permissions → **{new_profile}** (KIT mode: {chip_target})."
                    )
                    return
            state["engine"] = None
            cli_perm = _PROFILE_TO_CLI_PERM.get(new_profile, "default")
            _append_system_message(
                f"Permissions → **{new_profile}** (CLI: {cli_perm}). "
                f"Takes effect on next message."
            )
        # Warn on full mode
        if new_profile == "all_free":
            _append_system_message(
                "⚠ WARNING: **all_free** mode gives the agent unrestricted "
                "access to files, shell commands, and all directories "
                "(except archive & remote archive). Only use if you trust the setup."
            )
        try:
            from delfin.user_settings import load_settings, save_settings
            s = load_settings()
            s.setdefault("agent", {})
            s["agent"]["permission_profile"] = new_profile
            save_settings(s)
        except Exception:
            pass
        _update_status()

    def _on_commit(button):
        """Ask the agent to stage and commit current changes."""
        if state["streaming"]:
            return
        input_textarea.value = (
            "Please commit the current changes. "
            "Run `git diff --stat` to see what changed, then `git add` the relevant files "
            "and `git commit` with a concise, descriptive commit message in English. "
            "Do NOT push."
        )
        _on_send(None)

    def _on_push(button):
        """Show confirmation for git push."""
        import subprocess as _sp

        # Show current branch and unpushed commits
        try:
            branch = _sp.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                capture_output=True, text=True, cwd=str(ctx.repo_dir or "."),
            ).stdout.strip()
            unpushed = _sp.run(
                ["git", "log", "--oneline", "@{upstream}..HEAD"],
                capture_output=True, text=True, cwd=str(ctx.repo_dir or "."),
            ).stdout.strip()
            if not unpushed:
                push_status_html.value = (
                    '<span style="color:#757575;">No unpushed commits.</span>'
                )
                return
            push_status_html.value = (
                f'<span style="color:#ef6c00;">'
                f'Push <b>{branch}</b>? '
                f'({len(unpushed.splitlines())} commit(s))'
                f'</span>'
            )
        except Exception:
            push_status_html.value = (
                '<span style="color:#ef6c00;">Push to remote?</span>'
            )
        push_btn.layout.display = "none"
        push_confirm_btn.layout.display = "inline-flex"
        push_cancel_btn.layout.display = "inline-flex"

    def _on_push_confirm(button):
        """Execute the push after user confirmation."""
        import subprocess as _sp

        push_confirm_btn.layout.display = "none"
        push_cancel_btn.layout.display = "none"
        push_btn.layout.display = "inline-flex"
        try:
            result = _sp.run(
                ["git", "push"],
                capture_output=True, text=True, cwd=str(ctx.repo_dir or "."),
                timeout=30,
            )
            if result.returncode == 0:
                output = result.stdout.strip() or result.stderr.strip()
                push_status_html.value = (
                    f'<span style="color:#2e7d32;">'
                    f'\u2714 Pushed successfully. {_html.escape(output[:100])}'
                    f'</span>'
                )
                _append_system_message(f"Git push completed: {output[:200]}")
            else:
                push_status_html.value = (
                    f'<span style="color:#d32f2f;">'
                    f'Push failed: {_html.escape(result.stderr.strip()[:150])}'
                    f'</span>'
                )
        except Exception as exc:
            push_status_html.value = (
                f'<span style="color:#d32f2f;">Error: {_html.escape(str(exc))}</span>'
            )

    def _on_push_cancel(button):
        """Cancel the push."""
        push_confirm_btn.layout.display = "none"
        push_cancel_btn.layout.display = "none"
        push_btn.layout.display = "inline-flex"
        push_status_html.value = ""

    def _on_undo(button):
        """Revert the last file edit.

        The repo to operate on is determined by walking up from the
        edited file until a ``.git`` directory is found \u2014 necessary
        because the agent now writes into extra_workspace_dirs that are
        unrelated to ``ctx.repo_dir`` (e.g. Jerome's TestOpt project
        while the dashboard runs from the DELFIN repo).

        If the edit was a fresh ``write_file`` (untracked file in git),
        the file is unlinked instead of checked out \u2014 ``git checkout``
        on an untracked path fails silently and would leave the file in
        place.
        """
        if not state["recent_edits"]:
            return
        import subprocess as _sp
        from pathlib import Path as _Path
        last = state["recent_edits"].pop()
        fpath = _Path(last["file"]).expanduser().resolve()
        if not fpath.exists() and last.get("tool") not in ("Write", "write_file"):
            _append_system_message(f"Undo: file no longer exists: {fpath}")
            undo_btn.disabled = len(state["recent_edits"]) == 0
            return

        # Find the enclosing git repo by walking up.
        repo_root = None
        cur = fpath if fpath.is_dir() else fpath.parent
        while True:
            if (cur / ".git").exists():
                repo_root = cur
                break
            if cur == cur.parent:
                break
            cur = cur.parent

        short = str(fpath)
        if repo_root is not None:
            try:
                short = str(fpath.relative_to(repo_root))
            except ValueError:
                pass

        try:
            if repo_root is None:
                # Not under git \u2014 best we can do is unlink fresh writes;
                # for edits we have no baseline to roll back to.
                if last.get("tool") in ("Write", "write_file"):
                    fpath.unlink(missing_ok=True)
                    _append_system_message(
                        f"\u21a9 Removed (no git): {fpath}"
                    )
                else:
                    _append_system_message(
                        f"Undo unavailable for {fpath}: not under git, "
                        "and the original content was not saved before "
                        "the edit. Tip: keep the project under git for "
                        "reliable undo."
                    )
                undo_btn.disabled = len(state["recent_edits"]) == 0
                return

            # File under git: check if tracked or fresh.
            ls = _sp.run(
                ["git", "ls-files", "--error-unmatch", str(fpath)],
                capture_output=True, text=True,
                cwd=str(repo_root), timeout=5,
            )
            tracked = ls.returncode == 0

            if not tracked and last.get("tool") in ("Write", "write_file"):
                # Fresh untracked file the agent created \u2014 just unlink.
                fpath.unlink(missing_ok=True)
                _append_system_message(
                    f"\u21a9 Removed new file: {short}"
                )
            elif tracked:
                result = _sp.run(
                    ["git", "checkout", "--", str(fpath)],
                    capture_output=True, text=True,
                    cwd=str(repo_root), timeout=10,
                )
                if result.returncode == 0:
                    _append_system_message(f"\u21a9 Reverted: {short}")
                else:
                    _append_system_message(
                        f"Undo failed for {short}: "
                        f"{result.stderr.strip()[:120]}"
                    )
            else:
                _append_system_message(
                    f"Undo unavailable for untracked edit on {short}. "
                    "Add the file to git first if you want this kind of "
                    "undo to work."
                )
        except Exception as exc:
            _append_system_message(f"Undo error: {exc}")
        undo_btn.disabled = len(state["recent_edits"]) == 0

    def _on_load_session(button):
        sid = session_dropdown.value
        if not sid:
            # "New Session" selected — just reset
            _on_new_cycle(button)
            return
        if state["streaming"]:
            return
        _load_saved_session(sid)

    def _on_delete_session(button):
        sid = session_dropdown.value
        if not sid or state["streaming"]:
            return
        try:
            from delfin.agent.session_store import delete_session
            delete_session(sid)
        except Exception:
            pass
        # If we're viewing this session, clear the UI
        if state.get("active_session_id") == sid:
            _on_new_cycle(button)
        _refresh_session_dropdown()

    def _on_fork_session(button):
        sid = session_dropdown.value
        if not sid or state["streaming"]:
            return
        try:
            from delfin.agent.session_store import fork_session
            new_sid = fork_session(sid)
        except Exception as exc:
            _append_system_message(f"Fork failed: {exc}")
            return
        if not new_sid:
            _append_system_message(f"Could not fork session '{sid}'.")
            return
        _refresh_session_dropdown()
        session_dropdown.value = new_sid
        _load_saved_session(new_sid)
        _append_system_message(
            f"Forked session '{sid[:8]}…' → '{new_sid}'. "
            "You are now editing the fork; the original is unchanged."
        )

    # -- wire events -------------------------------------------------------
    send_btn.on_click(_on_send)
    stop_btn.on_click(_on_stop)
    new_cycle_btn.on_click(_on_new_cycle)
    advance_btn.on_click(_on_advance_role)
    inspector_primary_btn.on_click(_on_inspector_primary)
    inspector_retry_btn.on_click(_on_inspector_retry)
    inspector_stop_btn.on_click(_on_inspector_stop)
    inspector_next_btn.on_click(_on_inspector_next)
    inspector_detail_dropdown.observe(_on_inspector_detail_change, names="value")
    load_session_btn.on_click(_on_load_session)
    delete_session_btn.on_click(_on_delete_session)
    fork_session_btn.on_click(_on_fork_session)
    action_approve_btn.on_click(_on_actions_approve)
    action_deny_btn.on_click(_on_actions_deny)
    undo_btn.on_click(_on_undo)
    commit_btn.on_click(_on_commit)
    export_btn.on_click(_on_export)
    approve_btn.on_click(_on_approve)
    deny_btn.on_click(_on_deny)
    search_input.observe(_on_search_change, names="value")
    search_close_btn.on_click(_on_search_close)
    push_btn.on_click(_on_push)
    push_confirm_btn.on_click(_on_push_confirm)
    push_cancel_btn.on_click(_on_push_cancel)
    mode_dropdown.observe(_on_mode_change, names="value")
    provider_dropdown.observe(_on_provider_change, names="value")
    model_dropdown.observe(_on_model_change, names="value")
    effort_dropdown.observe(_on_effort_change, names="value")
    perm_dropdown.observe(_on_perm_change, names="value")

    # Initial KIT confirm-panel visibility based on saved provider.
    try:
        if provider_dropdown.value == "kit":
            _ensure_kit_broker()
            _show_kit_confirm_panel(True)
        else:
            _show_kit_confirm_panel(False)
    except Exception:
        pass

    # -- initial state -----------------------------------------------------
    _update_status()
    _update_button_states()
    _refresh_session_dropdown()

    # D5: prime the job-state snapshot so the first real transition we
    # observe doesn't get reported as a "first time seen" non-event.
    try:
        _check_job_events_once()
    except Exception:
        pass
    # Start the periodic watcher so users see job transitions without
    # actively polling /jobs.  The thread is a daemon and exits with
    # the Python session.
    try:
        _start_job_event_watcher()
    except Exception:
        pass

    # D2: observe dashboard tab changes for one-shot context suggestions.
    # ctx.tabs_widget is created later in delfin.dashboard.__init__ — it
    # may already be set or arrive shortly after.  We re-attempt on first
    # send if it's not ready yet.
    def _try_attach_tab_observer() -> None:
        tabs = getattr(ctx, "tabs_widget", None)
        if tabs is None:
            return
        try:
            tabs.observe(_on_dashboard_tab_change, names="selected_index")
        except Exception:
            pass

    _try_attach_tab_observer()

    # Apply solo-minimal UI at startup (default mode is dashboard)
    _init_mode = mode_dropdown.value
    if _init_mode in ("solo", "dashboard"):
        mode_desc_html.layout.display = "none"
        advance_btn.layout.display = "none"

    _enter_key_init_js = """
(function() {
    if (window.__delfinAgentKeys) return;
    window.__delfinAgentKeys = true;
    document.addEventListener('keydown', function(e) {
        /* Enter = Approve (if approval pending) or Send */
        if (e.key === 'Enter' && !e.shiftKey && !e.ctrlKey && !e.metaKey) {
            /* Check if approval buttons are visible */
            var approveBtn = document.querySelector('.delfin-agent-approval-row button');
            if (approveBtn && approveBtn.offsetParent !== null) {
                e.preventDefault();
                e.stopPropagation();
                approveBtn.click();
                return;
            }
            if (e.target && e.target.tagName === 'TEXTAREA') {
                var container = e.target.closest
                    ? e.target.closest('.delfin-agent-input') : null;
                if (container) {
                    e.preventDefault();
                    e.stopPropagation();
                    var sendBtn = document.querySelector('.delfin-agent-send-row button');
                    if (sendBtn) sendBtn.click();
                    return;
                }
            }
        }
        /* Escape = Deny (if approval pending) or Stop generation */
        if (e.key === 'Escape') {
            /* Check if deny button is visible */
            var approvalRow = document.querySelector('.delfin-agent-approval-row');
            if (approvalRow && approvalRow.offsetParent !== null) {
                var btns = approvalRow.querySelectorAll('button');
                if (btns.length >= 2) {
                    btns[1].click();
                    e.preventDefault();
                    return;
                }
            }
            var stopBtns = document.querySelectorAll('button');
            for (var i = 0; i < stopBtns.length; i++) {
                if (stopBtns[i].textContent.trim() === 'Stop' && !stopBtns[i].disabled) {
                    stopBtns[i].click();
                    e.preventDefault();
                    return;
                }
            }
        }
        /* Ctrl+L = Clear chat */
        if ((e.ctrlKey || e.metaKey) && e.key === 'l') {
            /* Only if focus is in the agent area */
            var agentArea = document.querySelector('.delfin-agent-chat');
            if (agentArea) {
                e.preventDefault();
                /* Trigger /clear by setting textarea and clicking send */
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var nativeSet = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    nativeSet.call(ta, '/clear');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sendBtn = document.querySelector('.delfin-agent-send-row button');
                        if (sendBtn) sendBtn.click();
                    }, 50);
                }
            }
        }
        /* Ctrl+K = Toggle search */
        if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
            var agentArea = document.querySelector('.delfin-agent-chat');
            if (agentArea) {
                e.preventDefault();
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var nativeSet = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    nativeSet.call(ta, '/search');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sendBtn = document.querySelector('.delfin-agent-send-row button');
                        if (sendBtn) sendBtn.click();
                    }, 50);
                }
            }
        }
        /* Shift+Tab = Cycle permission mode */
        if (e.key === 'Tab' && e.shiftKey && !e.ctrlKey && !e.metaKey) {
            var agentArea = document.querySelector('.delfin-agent-chat');
            if (agentArea) {
                e.preventDefault();
                e.stopPropagation();
                var ta = document.querySelector('.delfin-agent-input textarea');
                if (ta) {
                    var nativeSet = Object.getOwnPropertyDescriptor(
                        window.HTMLTextAreaElement.prototype, 'value').set;
                    nativeSet.call(ta, '/perm-cycle');
                    ta.dispatchEvent(new Event('input', {bubbles: true}));
                    setTimeout(function() {
                        var sendBtn = document.querySelector('.delfin-agent-send-row button');
                        if (sendBtn) sendBtn.click();
                    }, 50);
                }
            }
        }
    }, true);

    /* Auto-refresh model list when the model dropdown is opened.
       We listen on mousedown — that fires *before* the native select
       expands its options, so by the time the user picks a value the
       live list is already populated. We throttle to avoid hammering
       the API: max one fetch every 5 seconds per dropdown. */
    var _lastModelRefresh = 0;
    document.addEventListener('mousedown', function(e) {
        var wrapper = e.target.closest
            ? e.target.closest('.delfin-agent-model-dropdown')
            : null;
        if (!wrapper) return;
        var now = Date.now();
        if (now - _lastModelRefresh < 5000) return;
        _lastModelRefresh = now;
        var refresh = document.querySelector(
            '.delfin-agent-model-refresh button'
        );
        if (refresh) refresh.click();
    }, true);
})();
"""

    # Emit init-status banner for the initial /models fetch so the user
    # immediately knows whether the dropdown is showing live or fallback
    # data. Surfaces the API-key hint when the live fetch silently fell
    # back — the most common cause of "I only see 9 models".
    # Note: claude has no /models endpoint — its list is intentionally
    # hard-coded (CLI manages the actual model selection), so we don't
    # call that a "failure".
    _LIVE_FETCH_PROVIDERS = {"kit", "openai"}
    try:
        _mst = state.get("_models_init_status") or {}
        _mp = _mst.get("provider", "")
        _mc = _mst.get("count", 0)
        if _mp not in _LIVE_FETCH_PROVIDERS:
            # No live fetch for this provider by design (claude has no
            # /models endpoint, list is hard-coded). Stay quiet — the
            # dropdown already shows what's available.
            pass
        elif _mst.get("live"):
            # Stay quiet on success too — no news is good news.
            pass
        else:
            if not _mst.get("key_present"):
                _key_name = (
                    "KIT_TOOLBOX_API_KEY" if _mp == "kit"
                    else "OPENAI_API_KEY"
                )
                _hint = (
                    f" — {_key_name} is not set in the dashboard "
                    f"process. Set it before startup or via a notebook "
                    f"cell (`import os; "
                    f"os.environ['{_key_name}']='…'`) and press ↻."
                )
            else:
                _hint = (
                    " — Live fetch failed (network/auth). "
                    "Press ↻ next to the model dropdown to retry."
                )
            _append_system_message(
                f"⚠ Models: fallback active for {_mp} "
                f"({_mc} models){_hint}"
            )
    except Exception:
        pass

    tab_widget = agent_content
    # A2b — expose the agent-tab state to other tabs (Activity tab's live
    # pane reads engine status, session cost, streaming flag, perms).
    # Read-only contract: other tabs MUST NOT mutate this dict.
    try:
        ctx.agent_state = state
    except Exception:
        pass

    # --resume <sid> boot flag: when ctx.initial_session_id is set or the
    # DELFIN_RESUME_SESSION env var points at a saved id (or "latest"),
    # auto-load that session immediately after the widgets render. The
    # user lands directly inside the previous conversation instead of an
    # empty tab.
    try:
        import os as _os
        boot_sid = getattr(ctx, "initial_session_id", "") or _os.environ.get(
            "DELFIN_RESUME_SESSION", ""
        )
        boot_sid = (boot_sid or "").strip()
        if boot_sid:
            if boot_sid == "latest":
                try:
                    from delfin.agent.session_store import resume_latest
                    _data = resume_latest()
                    if _data:
                        boot_sid = str(_data.get("session_id") or "")
                    else:
                        boot_sid = ""
                except Exception:
                    boot_sid = ""
            if boot_sid:
                try:
                    _load_saved_session(boot_sid)
                except Exception:
                    pass
    except Exception:
        pass

    return tab_widget, {"init_js": _enter_key_init_js}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _format_role_label(role_id: str) -> str:
    """Convert role_id to a readable label."""
    if not role_id:
        return "Agent"
    return role_id.replace("_", " ").title()


def _render_status(
    mode: str,
    backend: str,
    role: str,
    role_index: int,
    role_total: int,
    input_tokens: int,
    output_tokens: int,
    cost_usd: float,
    provider: str = "claude",
    perm_profile: str = "ask_all",
    active_gate_type: str = "",
    active_gate_text: str = "",
) -> str:
    """Render the status bar HTML."""
    role_label = _format_role_label(role)
    role_info = ""
    if role_total > 0:
        role_info = (
            f'<span class="role-badge">{_html.escape(role_label)} '
            f"({role_index + 1}/{role_total})</span>"
        )

    if provider == "kit":
        backend_label = "KIT Toolbox"
    elif provider == "openai":
        backend_label = "Codex CLI" if backend == "cli" else "OpenAI API"
    elif backend == "cli":
        backend_label = "CLI (OAuth)"
    else:
        backend_label = "API"
    backend_info = f'<span class="backend-badge">{backend_label}</span>'

    if cost_usd > 0:
        cost_str = f"${cost_usd:.3f}"
    else:
        cost_str = _estimate_cost_str(backend, input_tokens, output_tokens,
                                      provider=provider)

    tokens_str = f"{input_tokens:,} in / {output_tokens:,} out"

    # Permission profile badge (color-coded)
    _perm_colors = {
        "plan": "#6c757d",      # gray
        "ask_all": "#0d6efd",   # blue
        "repo_free": "#198754", # green
        "all_free": "#dc3545",  # red
    }
    perm_color = _perm_colors.get(perm_profile, "#6c757d")
    perm_badge = (
        f'<span class="backend-badge" style="background:{perm_color}">'
        f'{_html.escape(perm_profile)}</span>'
    )

    gate_info = ""
    if active_gate_type:
        _gate_colors = {
            "schema": ("#ef4444", "white"),
            "risk": ("#f59e0b", "#fff7ed"),
            "partial": ("#2563eb", "white"),
            "goal-lock": ("#2563eb", "white"),
            "review": ("#2563eb", "white"),
            "question": ("#8b5cf6", "white"),
            "plan-approval": ("#8b5cf6", "white"),
            "findings": ("#8b5cf6", "white"),
            "approval": ("#8b5cf6", "white"),
            "confidence": ("#8b5cf6", "white"),
            "conflict": ("#8b5cf6", "white"),
            "cost": ("#8b5cf6", "white"),
        }
        bg, fg = _gate_colors.get(active_gate_type, ("#6b7280", "white"))
        label = active_gate_type.replace("-", " ")
        label = " ".join(part.capitalize() for part in label.split())
        gate_text = _html.escape(active_gate_text[:80]) if active_gate_text else ""
        gate_info = (
            f'<span class="gate-badge" style="background:{bg};color:{fg}">{_html.escape(label)}</span>'
            f'<span class="gate-text">{gate_text}</span>'
        )

    return (
        f'<div class="delfin-agent-status">'
        f'<span class="mode-badge">{_html.escape(mode)}</span>'
        f"{role_info}"
        f"{backend_info}"
        f"{perm_badge}"
        f"{gate_info}"
        f'<span class="tokens-info">{tokens_str} · {cost_str}</span>'
        f"</div>"
    )


def _estimate_cost_str(
    backend: str, input_tokens: int, output_tokens: int,
    provider: str = "claude",
) -> str:
    """Rough cost string for the dashboard status row."""
    if provider == "kit":
        # KIT-Toolbox is provided by KIT — no per-call USD cost, but
        # there's a quota and the user explicitly asked to see the
        # burn rate. Showing tokens spent makes runaway agents
        # visible without faking a price.
        total = input_tokens + output_tokens
        if total <= 0:
            return "KIT (no usage yet)"
        return f"KIT quota: {total:,} tokens"
    if backend == "cli":
        return "included in subscription"
    if provider == "openai":
        cost = (input_tokens * 2.0 + output_tokens * 8.0) / 1_000_000
    else:
        cost = (input_tokens * 3.0 + output_tokens * 15.0) / 1_000_000
    return f"~${cost:.3f}"
