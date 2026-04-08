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
    max-height: calc(100vh - 420px);
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
    margin-bottom: 12px;
    padding: 8px 12px;
    border-radius: 8px;
    max-width: 85%;
    word-wrap: break-word;
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
    border-radius: 6px;
    align-items: center;
    gap: 8px;
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
    display: inline-block;
    padding: 3px 10px;
    border-radius: 12px;
    background: #fef3c7;
    color: #92400e;
    font-size: 12px;
    font-weight: 600;
    animation: delfin-pulse 1.5s ease-in-out infinite;
    max-width: 90%;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
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
    # Fallback: always show at least Claude (will error with helpful message)
    if not _available_providers:
        _available_providers.append(("Claude", "claude"))

    # -- state -------------------------------------------------------------
    state = {
        "engine": None,
        "chat_messages": [],
        "streaming": False,
        "active_session_id": "",  # currently loaded CLI session ID
        "recent_edits": [],       # list of {"file": path, "tool": name} for undo
        "message_queue": [],      # queued messages sent while agent is busy
        "session_start_time": None,  # monotonic time of first message
        "_agent_calc_path": "",       # relative path within calc_dir for browsing
        "_pending_dashboard_action": None,  # {action_id, description, callback}
        "_perm_profile": "ask_all",  # permission profile: plan/ask_all/repo_free/all_free
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
        options=["dashboard", "solo", "research", "quick", "reviewed", "tdd", "cluster", "full"],
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
        "kit": [
            ("Azure GPT-4.1-mini", "azure.gpt-4.1-mini"),
            ("Azure GPT-4.1", "azure.gpt-4.1"),
            ("Azure GPT-4.1-nano", "azure.gpt-4.1-nano"),
        ],
    }
    _PROVIDER_DEFAULTS = {"claude": "sonnet", "openai": "gpt-5.4",
                          "kit": "azure.gpt-4.1-mini"}
    _PROVIDER_CHEAP = {"claude": "haiku", "openai": "gpt-5.4-mini",
                       "kit": "azure.gpt-4.1-nano"}

    # Skip patterns: models that should not appear in the dropdown
    _KIT_SKIP = {"standard-external", "standard-local"}

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
            else:
                return None

            models = data.get("data", [])
            result = []
            for m in models:
                mid = m.get("id", "")
                if not mid or mid in _KIT_SKIP:
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
    # Effort selector (only affects API backend; CLI manages thinking internally)
    effort_dropdown = widgets.Dropdown(
        options=[
            ("Low", "low"),
            ("Medium", "medium"),
            ("High", "high"),
        ],
        value="medium",
        description="Effort:",
        layout=widgets.Layout(width="155px"),
        style={"description_width": "42px"},
        tooltip="Thinking budget (only works with API backend; CLI manages thinking internally)",
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
        tooltip="Permission profile: Plan=read-only, Default=ask all changes, "
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
        if _saved_effort in ("low", "medium", "high"):
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
            [mode_dropdown, provider_dropdown, model_dropdown, effort_dropdown, perm_dropdown,
             new_cycle_btn, advance_btn, stop_btn, undo_btn, export_btn,
             commit_btn, push_btn, push_confirm_btn, push_cancel_btn, push_status_html],
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
    session_row = widgets.HBox(
        [session_dropdown, load_session_btn, delete_session_btn],
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

    # Status bar
    status_html = widgets.HTML(
        value=_render_status("quick", "cli", "session_manager", 0, 3, 0, 0, 0.0),
    )

    # Chat display
    chat_html = widgets.HTML(
        value='<div class="delfin-agent-chat"><i>Start a conversation...</i></div>',
        layout=widgets.Layout(min_height="200px"),
    )

    # Input area
    input_textarea = widgets.Textarea(
        placeholder="Type your message here... (Enter = Send, Shift+Enter = new line, queues while busy)",
        layout=widgets.Layout(width="100%", height="70px"),
    )
    input_textarea.add_class("delfin-agent-input")
    send_btn = widgets.Button(
        description="Send",
        button_style="primary",
        layout=widgets.Layout(width="80px", height="70px"),
    )
    input_row = widgets.HBox(
        [input_textarea, send_btn],
        layout=widgets.Layout(margin="6px 0 0 0"),
    )
    input_row.add_class("delfin-agent-send-row")

    # Working indicator (animated spinner)
    working_html = widgets.HTML(value="")

    # Queue indicator (shown when messages are queued during streaming)
    queue_html = widgets.HTML(value="")

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

    # Widgets the agent must NEVER click (hard block)
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
            "orca-submit-btn": "orca_submit_btn",
            "orca-preview": "orca_preview",
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

    # -- layout assembly ---------------------------------------------------
    agent_content = widgets.VBox(
        [css_widget, _enter_js_output, controls_row, session_row, search_row,
         status_html, chat_html, working_html, queue_html, approval_row, input_row],
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
        if not engine or not engine.session_id:
            return
        try:
            from delfin.agent.session_store import save_session
            estate = engine.export_state()
            save_session(
                session_id=engine.session_id,
                mode=estate["mode"],
                role_index=estate["role_index"],
                route=estate["route"],
                role_outputs=estate["role_outputs"],
                chat_messages=state["chat_messages"],
                engine_messages=estate["engine_messages"],
                token_usage=estate["token_usage"],
                cost_usd=estate["cost_usd"],
            )
            state["active_session_id"] = engine.session_id
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
        mode_dropdown.value = saved_mode

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
        state["active_session_id"] = session_id
        _refresh_chat_html()
        _update_status()
        _update_button_states()

        title = data.get("title", "")[:50] or session_id[:12]
        _append_system_message(f"Session restored: {title}")

    # -- helpers -----------------------------------------------------------

    def _get_agent_settings():
        """Load agent settings from user settings."""
        try:
            from delfin.user_settings import load_settings
            return load_settings().get("agent", {}) or {}
        except Exception:
            return {}

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

            engine = AgentEngine(
                repo_dir=repo_dir,
                backend=backend,
                provider=provider,
                api_key=api_key,
                model=model,
                mode=mode_dropdown.value,
                permission_mode=_active_cli_perm(),
                mcp_config=_mcp_cfg,
            )
            state["engine"] = engine
            ctx.agent_engine = engine
            return engine
        except Exception as exc:
            _append_system_message(f"Engine error: {exc}")
            return None

    def _append_chat_message(role, content, role_label=""):
        state["chat_messages"].append(
            {"role": role, "content": content, "role_label": role_label}
        )
        _refresh_chat_html()

    def _append_system_message(text):
        _append_chat_message("system", text)

    def _update_last_assistant(content, role_label=""):
        """Update the last assistant message (for streaming)."""
        msgs = state["chat_messages"]
        if msgs and msgs[-1]["role"] == "assistant":
            msgs[-1]["content"] = content
            if role_label:
                msgs[-1]["role_label"] = role_label
        else:
            msgs.append(
                {"role": "assistant", "content": content, "role_label": role_label}
            )
        _refresh_chat_html()

    def _refresh_chat_html():
        """Rebuild the chat display HTML."""
        if not state["chat_messages"]:
            chat_html.value = (
                '<div class="delfin-agent-chat">'
                "<i>Start a conversation...</i></div>"
            )
            return
        parts = ['<div class="delfin-agent-chat">']
        for msg in state["chat_messages"]:
            role = msg["role"]
            content = _md_to_html(msg["content"])
            if role == "user":
                parts.append(
                    f'<div class="delfin-chat-msg delfin-chat-user">{content}</div>'
                )
            elif role == "assistant":
                label = msg.get("role_label", "Agent")
                parts.append(
                    f'<div class="delfin-chat-msg delfin-chat-agent">'
                    f'<div class="delfin-chat-role">{_html.escape(label)}</div>'
                    f"{content}</div>"
                )
            elif role == "thinking":
                raw = msg["content"]
                # Truncate for display, keep first 500 chars
                preview = raw[:80].replace("\n", " ").strip()
                if len(raw) > 80:
                    preview += "..."
                escaped_full = _html.escape(raw)
                parts.append(
                    f'<details class="delfin-chat-msg delfin-chat-thinking">'
                    f'<summary>\U0001f9e0 {_html.escape(preview)}</summary>'
                    f'<div class="thinking-content">{escaped_full}</div>'
                    f'</details>'
                )
            elif role == "approval":
                parts.append(
                    f'<div class="delfin-chat-msg delfin-chat-approval">'
                    f'\u26a0\ufe0f {content}'
                    f'<div style="margin-top:4px;font-size:11px;color:#78350f;">'
                    f'Press <b>Enter</b> to approve, <b>Esc</b> to deny</div>'
                    f'</div>'
                )
            elif role == "system":
                raw = msg["content"]
                # Handoff/cycle messages get special styling
                if raw.startswith("---") or raw.startswith("Session restored"):
                    css_class = "delfin-chat-msg delfin-chat-handoff"
                else:
                    css_class = "delfin-chat-msg delfin-chat-system"
                parts.append(f'<div class="{css_class}">{content}</div>')
        # Scroll to bottom: <img onerror> is the only way to run JS inside
        # ipywidgets.HTML (script tags are sanitized, but img onerror works).
        parts.append(
            '<img src="" onerror="'
            "var c=this.closest('.delfin-agent-chat');"
            "if(c)c.scrollTop=c.scrollHeight;"
            "this.remove();"
            '" style="display:none">'
        )
        parts.append("</div>")
        chat_html.value = "\n".join(parts)

    def _update_status():
        """Update the status bar from engine state."""
        engine = state["engine"]
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
            )
        else:
            backend = _resolve_backend() if _cli_available else "api"
            status_html.value = _render_status(
                mode_dropdown.value, backend, "", 0, 0, 0, 0, 0.0,
                provider=provider_dropdown.value,
                perm_profile=state.get("_perm_profile", "ask_all"),
            )

    def _set_working(active, label=""):
        """Show or hide the working indicator."""
        if active:
            text = label or "Working..."
            working_html.value = (
                f'<span class="delfin-agent-working">'
                f'\u23f3 {_html.escape(text)}</span>'
            )
            # Global header spinner — visible from all tabs
            short = text[:40] + ("…" if len(text) > 40 else "")
            ctx.agent_status_html.value = (
                f'<span style="font-family:monospace; color:#1565c0; '
                f'padding:2px 6px; border:1px solid #90caf9; '
                f'border-radius:4px; animation:pulse 1.5s infinite;">'
                f'\u23f3 Agent: {_html.escape(short)}</span>'
                f'<style>@keyframes pulse{{0%,100%{{opacity:1}}50%{{opacity:.5}}}}</style>'
            )
        else:
            working_html.value = ""
            ctx.agent_status_html.value = ""

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

        if cmd == "/help":
            _append_system_message(
                "Available commands:\n"
                "  /help            — Show this help\n"
                "  /clear           — Clear chat history\n"
                "  /cost            — Show token usage & cost\n"
                "  /compact         — Summarize context (reduce tokens)\n"
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
                "  /effort <lvl>    — Set effort (low/medium/high)\n"
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
                "  /calc select <f> — Select file in browser (opens options dropdown)\n"
                "  /calc read <file> — Read a calc file\n"
                "  /calc tail <file> — Read last 8KB of output\n"
                "  /calc info <dir> — Show folder summary & status\n"
                "  /calc tree [dir] — Show directory tree\n"
                "  /calc search <p> — Search files by glob pattern\n"
                "  /analyze <dir>   — Full analysis (energy+convergence+errors)\n"
                "  /analyze energy <dir> — Extract energies (Gibbs/ZPE/electronic)\n"
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
                "\n"
                "Memory:\n"
                "  /remember <text> — Save a persistent memory\n"
                "  /memories        — List all memories\n"
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

        # /effort <level>
        if cmd.startswith("/effort "):
            level = cmd[8:].strip()
            valid = {"low", "medium", "high"}
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
            opts = [v for _, v in mode_dropdown.options] if hasattr(mode_dropdown, 'options') else []
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
                from delfin.agent.memory_store import save_memory
                idx = save_memory(text_to_save, source="user")
                _append_system_message(f"Memory [{idx}] saved: {text_to_save}")
            else:
                _append_system_message("Usage: /remember <text>")
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
                "archive": "Archive",
                "agent": "Agent",
                "settings": "Settings",
            }
            title = _tab_map.get(tab_name.lower(), tab_name)
            if ctx.select_tab(title):
                _append_system_message(f"Switched to tab: {title}")
            else:
                avail = ", ".join(_tab_map.keys())
                _append_system_message(f"Tab not found: {tab_name}\nAvailable: {avail}")
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
                if _ui_wname in _BLOCKED_WIDGETS:
                    _append_system_message(
                        "Blocked: agent cannot use delete buttons. "
                        "Use the dashboard UI directly."
                    )
                    return True
                if not hasattr(_ui_w, "click"):
                    _append_system_message(f"Widget '{_ui_wname}' is not a button.")
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
            if cw:
                cw.value = content
                _append_system_message(
                    f"CONTROL updated ({len(content)} chars). "
                    "Switch to Submit tab to review."
                )
            else:
                _append_system_message("Submit tab not available.")
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
                _append_system_message(f"Unknown param: {param}\nAvailable: {', '.join(_orca_param_map.keys())}")
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

        # /calc select <filename> — select a file in the browser (triggers options dropdown)
        if cmd.startswith("/calc select "):
            filename = text[len("/calc select"):].strip()
            _path_w = ctx.calc_browser_refs.get("calc_path_input")
            if not _path_w:
                _append_system_message("Calc browser not available.")
                return True
            # Navigate the browser directly to the file — this triggers
            # the browser's own path handler which opens the file and
            # populates the options dropdown (Recalc/Smart Recalc/etc.)
            target = _resolve_calc_path(filename)
            if not target.exists():
                _append_system_message(f"Not found: {filename}")
                return True
            _path_w.value = str(target)
            _append_system_message(f"Selected: {target.name}")
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

        if cmd.startswith("/analyze ") and not cmd.startswith("/analyze energy") and not cmd.startswith("/analyze convergence") and not cmd.startswith("/analyze errors") and cmd != "/analyze status":
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
            folder = text[len("/recalc"):].strip()
            if folder == "auto":
                pass  # handled above
            else:
                target = _resolve_calc_path(folder)
                if not target.is_dir():
                    _append_system_message(f"Not a directory: {folder}")
                    return True
                control = target / "CONTROL.txt"
                if not control.exists():
                    _append_system_message(f"No CONTROL.txt in {folder}")
                    return True
                def _do_recalc(job_dir=target):
                    try:
                        result = ctx.backend.submit_delfin(
                            job_dir=job_dir, job_name=job_dir.name,
                            mode="delfin-recalc-classic",
                        )
                        if result.returncode == 0:
                            _append_system_message(f"Recalc submitted for {job_dir.name}.")
                        else:
                            _append_system_message(f"Recalc failed: {result.stderr[:300]}")
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
                    except Exception as exc:
                        _append_system_message(f"Recalc error: {exc}")
                _confirm_or_exec(
                    "recalc_single",
                    f"Submit recalc for '{folder}'?",
                    _do_recalc,
                    cmd_for_zone=f"/recalc {folder}",
                )
            return True

        # -- Cancel jobs (CONFIRMATION REQUIRED) -----------------------------

        if cmd.startswith("/cancel "):
            target = text[len("/cancel"):].strip()
            if target == "all":
                try:
                    jobs = ctx.backend.list_jobs()
                    active = [j for j in jobs if j.status.upper() in ("RUNNING", "R", "PENDING", "PD")]
                except Exception:
                    active = []
                if not active:
                    _append_system_message("No active jobs to cancel.")
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
                    "cancel_all",
                    f"Cancel {len(active)} active jobs: {names}",
                    _do_cancel_all,
                    cmd_for_zone="/cancel all",
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

        return False

    # -- dashboard mode helpers --------------------------------------------

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
        # calc_dir + workspace + permissions
        parts.append(f"Calculations dir: {ctx.calc_dir}")
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

    # -- command safety tiers (enforced at CODE level, not prompt) ----------

    _TIER3_EXACT = {"/submit", "/orca submit", "/recalc auto", "/cancel all"}
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
            # Read-only everywhere — agent can only browse and analyze
            "workspace":      (0, True),
            "calc":           (0, True),
            "repo":           (0, True),
            "archive":        (0, True),
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

    # Map DELFIN profile → Claude CLI permission_mode
    _PROFILE_TO_CLI_PERM: dict[str, str] = {
        "plan":      "plan",
        "ask_all":   "default",
        "repo_free": "bypassPermissions",  # CLI runs non-interactive; DELFIN zone system enforces safety
        "all_free":  "bypassPermissions",  # CLI runs non-interactive; DELFIN zone system enforces safety
    }

    def _active_perms() -> dict[str, tuple[int, bool]]:
        """Return the zone permissions for the active profile."""
        profile = state.get("_perm_profile", "ask_all")
        return _PERM_PROFILES.get(profile, _PERM_PROFILES["ask_all"])

    def _active_cli_perm() -> str:
        """Return the Claude CLI permission_mode for the active profile."""
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
    # last message.  This prevents the agent from self-triggering bulk ops
    # by mentioning just one keyword in its own output.
    _BULK_ACTION_KW = (
        "recalc", "neuberechn", "submit", "absend", "abschick",
        "cancel", "abbrech", "stopp",
    )
    _BULK_SCOPE_KW = (
        "all", "alle", "auto", "alles", "every", "sämtliche", "bulk",
        "komplett", "gesamt",
    )

    def _dashboard_auto_exec(agent_text: str):
        """Scan agent output for ACTION: /command lines and execute them.

        Safety enforcement (code-level, not prompt-level):
        - Zone-based permissions: workspace=free, calc=confirm,
          archive/remote_archive=read-only (HARD BLOCK), unknown=blocked
        - Tier 3: max 1 per response, bulk ops need explicit user intent
        - Workspace zone: tier 3 skips confirmation gate
        """
        import re as _re

        lines = agent_text.split("\n")
        commands: list[str] = []
        i = 0
        while i < len(lines):
            m = _re.match(r"^ACTION:\s*(/\S+.*)$", lines[i])
            if m:
                cmd = m.group(1)
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

                # Bulk ops need explicit user intent (action + scope)
                cl = cmd_line.lower().strip()
                if cl in ("/recalc auto", "/cancel all"):
                    user_msg = state.get("_last_user_message", "").lower()
                    has_action = any(kw in user_msg for kw in _BULK_ACTION_KW)
                    has_scope = any(kw in user_msg for kw in _BULK_SCOPE_KW)
                    if not (has_action and has_scope):
                        _append_system_message(
                            "\u26d4 Blocked: Bulk operation requires explicit "
                            "user request (e.g. 'recalc alle'). "
                            "Report findings and ask the user."
                        )
                        results.append("BLOCKED: bulk op without user intent")
                        continue

            # --- Execute ---
            _append_system_message(f"\u25b6 Executing: {short}")
            n_before = len(state["chat_messages"])
            try:
                handled = _handle_slash_command(cmd_line)
                if not handled:
                    _append_system_message(f"Unknown command: {short}")
            except Exception as exc:
                _append_system_message(f"Error executing {short}: {exc}")
            # Collect new system messages as feedback
            for msg in state["chat_messages"][n_before:]:
                if msg["role"] == "system":
                    results.append(msg["content"])
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
            elif role == "system":
                lines.append(f"> {content}\n")
        md_text = "\n".join(lines)
        # Write to file
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        export_dir = Path.home() / ".delfin" / "exports"
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
        while state["chat_messages"] and state["chat_messages"][-1]["role"] in ("assistant", "thinking", "system"):
            state["chat_messages"].pop()
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
        _append_chat_message(
            "approval",
            f"\u26a0\ufe0f {description}",
        )
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

    def _show_approval_prompt(tool_name, detail):
        """Show approval request inline in chat + approval buttons."""
        state["_pending_approval"] = {"tool": tool_name, "detail": detail}
        readable = _format_tool_description(detail)
        # Show in chat as a special approval message
        _append_chat_message(
            "approval",
            readable,
        )
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

        # --- Bash commands: run directly, report result to agent ---
        if tool == "Bash" and inp.get("command"):
            cmd = inp["command"]
            _append_system_message(f"\u2705 Approved & executing: $ {cmd[:200]}")
            try:
                result = _sp.run(
                    cmd, shell=True, capture_output=True, text=True,
                    cwd=str(ctx.repo_dir or "."), timeout=60,
                )
                output = (result.stdout.strip() + "\n" + result.stderr.strip()).strip()
                if result.returncode == 0:
                    _append_system_message(
                        f"\u2714 Command succeeded:\n{output[:500]}"
                    )
                    # Tell agent it worked
                    input_textarea.value = (
                        f"I ran the command for you. Result:\n```\n{output[:1000]}\n```\n"
                        f"Continue with the next step."
                    )
                else:
                    _append_system_message(
                        f"\u2718 Command failed (exit {result.returncode}):\n{output[:500]}"
                    )
                    input_textarea.value = (
                        f"The command failed (exit {result.returncode}):\n"
                        f"```\n{output[:1000]}\n```\n"
                        f"Please suggest an alternative approach."
                    )
            except _sp.TimeoutExpired:
                _append_system_message("\u26a0 Command timed out after 60s.")
                input_textarea.value = "The command timed out. Please try a different approach."
            except Exception as exc:
                _append_system_message(f"\u26a0 Error: {exc}")
                input_textarea.value = f"Error running command: {exc}"
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

        # Slash commands execute immediately, even during streaming
        if user_text.startswith("/"):
            input_textarea.value = ""
            _append_chat_message("user", user_text)
            _handle_slash_command(user_text)
            return

        # If streaming, queue the message for later
        if state["streaming"]:
            state["message_queue"].append(user_text)
            input_textarea.value = ""
            _append_system_message(
                f"Message queued ({len(state['message_queue'])} in queue). "
                f"Will send when current response completes."
            )
            _update_queue_display()
            return

        # Track user message for safety intent-checking
        state["_last_user_message"] = user_text

        engine = _ensure_engine()
        if engine is None:
            return

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
            _append_system_message(
                f"--- Follow-up mode ({_fu_role.replace('_', ' ').title()}) "
                f"— continue chatting or /reset for new task ---"
            )

        # Auto-mode suggestion: only on first message of a cycle, and only
        # if user hasn't already been asked.  User must accept the switch.
        if not engine.messages and not state.get("_mode_suggested"):
            from delfin.agent.engine import AgentEngine as _AE
            suggested = _AE.suggest_mode(user_text, mode_dropdown.value)
            if suggested:
                state["_mode_suggested"] = True
                state["_pending_mode_msg"] = user_text
                _append_chat_message("user", user_text)
                _append_system_message(
                    f"Mode suggestion: your message mentions files that "
                    f"match **{suggested}** mode (current: {mode_dropdown.value}).\n"
                    f"Type `/mode {suggested}` to switch, or just send "
                    f"your next message to keep **{mode_dropdown.value}**."
                )
                input_textarea.value = ""
                return

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

        input_textarea.value = ""
        _append_chat_message("user", user_text)

        state["streaming"] = True
        state["_deny_count"] = 0  # Reset retry counter for new message
        if state["session_start_time"] is None:
            state["session_start_time"] = time.monotonic()
        _update_button_states()
        _set_working(True, "Thinking...")

        role_label = _format_role_label(engine.current_role)

        def _worker():
            chunks = []
            thinking_chunks = []
            last_update = 0.0
            try:

                def _on_thinking(text):
                    nonlocal last_update
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
                    if now - last_update > 0.1:
                        _set_working(True, "Writing...")
                        _update_last_assistant("".join(chunks), role_label)
                        _update_status()
                        last_update = now

                def _on_tool_use(tool_name, tool_input):
                    # Flush thinking if no text came before tool use
                    if thinking_chunks:
                        full_thinking = "".join(thinking_chunks)
                        if full_thinking.strip():
                            _append_chat_message("thinking", full_thinking.strip())
                        thinking_chunks.clear()
                    # Show detailed activity like Claude CLI
                    try:
                        import json as _j
                        _p = _j.loads(tool_input)
                    except Exception:
                        _p = {}
                    _fname = (_p.get("file_path") or _p.get("path") or "")
                    if _fname:
                        _fname = _fname.rsplit("/", 1)[-1]  # basename only
                    _detail = {
                        "Read":  f"Reading {_fname}..." if _fname else "Reading...",
                        "Edit":  f"Editing {_fname}..." if _fname else "Editing...",
                        "Write": f"Writing {_fname}..." if _fname else "Writing...",
                        "Grep":  f"Searching: {_p.get('pattern', '')[:40]}...",
                        "Glob":  f"Finding: {_p.get('pattern', '')[:40]}...",
                        "Bash":  f"$ {(_p.get('command') or '')[:50]}...",
                        "Agent": f"Sub-agent: {(_p.get('description') or '')[:40]}...",
                    }.get(tool_name, f"Running {tool_name}...")
                    _set_working(True, _detail)
                    # Flush pending text as a finalized assistant message
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label)
                        chunks.clear()  # reset — next text starts fresh

                    try:
                        import json as _json
                        parsed = _json.loads(tool_input)
                    except Exception:
                        parsed = {}

                    # --- Edit/Write: show diff preview ---
                    if tool_name in ("Edit", "Write") and parsed.get("file_path"):
                        fpath = parsed["file_path"]
                        short_path = fpath.replace(
                            str(ctx.repo_dir or ""), ""
                        ).lstrip("/")
                        if tool_name == "Edit":
                            old = parsed.get("old_string", "")
                            new = parsed.get("new_string", "")
                            old_preview = old[:150] + ("..." if len(old) > 150 else "")
                            new_preview = new[:150] + ("..." if len(new) > 150 else "")
                            _append_system_message(
                                f"\u270f Edit: {short_path}\n"
                                f"  - {old_preview}\n"
                                f"  + {new_preview}"
                            )
                        else:
                            content = parsed.get("content", "")
                            _append_system_message(
                                f"\u270f Write: {short_path} "
                                f"({len(content)} chars)"
                            )
                        state["recent_edits"].append({
                            "file": fpath, "tool": tool_name,
                        })
                        undo_btn.disabled = False
                        return

                    # --- Bash: show command ---
                    if tool_name == "Bash" and parsed.get("command"):
                        cmd = parsed["command"]
                        if len(cmd) > 120:
                            cmd = cmd[:120] + "..."
                        _append_system_message(f"\u25b8 $ {cmd}")
                        return

                    # --- Other tools: show key param ---
                    short = ""
                    for key in ("file_path", "path", "pattern",
                                "prompt", "description", "query"):
                        if key in parsed:
                            val = str(parsed[key])
                            if len(val) > 80:
                                val = val[:80] + "..."
                            short = val
                            break
                    if not short:
                        short = tool_input[:120]
                    _append_system_message(
                        f"\u25b8 {tool_name}({short})"
                    )

                def _on_permission_denied(description):
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label)
                        chunks.clear()
                    denied_raw = str(description)
                    readable = _format_tool_description(denied_raw)
                    state["_last_denied"] = denied_raw
                    # Track denials to stop retry loops
                    deny_count = state.get("_deny_count", 0) + 1
                    state["_deny_count"] = deny_count
                    if deny_count >= 3:
                        _append_system_message(
                            f"\u26d4 Blocked 3x — stopping retries. "
                            f"Change permission mode or do it manually."
                        )
                        engine.request_stop()
                        return
                    _append_system_message(
                        f"\u26d4 Blocked ({deny_count}/3): {readable}"
                    )
                    # Show approval buttons for interactive approval
                    _show_approval_prompt(denied_raw, denied_raw)

                # Effort multiplier: scales the role-specific budget
                _effort_mult = {"low": 0.5, "medium": 1.0, "high": 2.0}
                _mult = _effort_mult.get(effort_dropdown.value, 1.0)

                # Store original user task for handoff messages
                original_task = user_text
                # Dashboard mode: inject current widget state so the agent
                # can read CONTROL, ORCA settings, etc. without tool calls.
                if mode_dropdown.value == "dashboard":
                    _ctx_text = _build_dashboard_context()
                    current_msg = (
                        f"[Dashboard state]\n{_ctx_text}\n\n"
                        f"[User request]\n{user_text}"
                    )
                else:
                    current_msg = user_text
                max_auto_steps = len(engine.route) + 1  # safety limit

                for _step in range(max_auto_steps):
                    if engine._stop_requested or engine.is_cycle_complete:
                        break

                    # Role-specific thinking budget and model routing
                    from delfin.agent.engine import AgentEngine as _AE
                    _cur_role = engine.current_role
                    _base_budget = _AE.thinking_budget_for_role(_cur_role)
                    _budget = min(int(_base_budget * _mult), 128000)

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
                    _memory = format_memory_context()

                    engine.stream_response(
                        user_message=current_msg,
                        on_token=_on_token,
                        on_tool_use=_on_tool_use,
                        on_permission_denied=_on_permission_denied,
                        on_thinking=_on_thinking,
                        thinking_budget=_budget,
                        memory_context=_memory,
                    )
                    # Final update for this role
                    if chunks:
                        _update_last_assistant("".join(chunks), role_label)

                    # Dashboard mode: auto-execute slash commands from agent output
                    # and strip ACTION: lines from displayed chat.
                    # Continuation loop: if commands returned results, let the
                    # agent process them (up to _MAX_CONT turns to prevent runaway).
                    _MAX_DASHBOARD_CONT = 10
                    _cont_turn = 0
                    while mode_dropdown.value == "dashboard" and chunks and _cont_turn < _MAX_DASHBOARD_CONT:
                        _cont_turn += 1
                        raw = "".join(chunks)
                        exec_results = _dashboard_auto_exec(raw)
                        # Remove ACTION: lines from visible output
                        import re as _re_strip
                        cleaned = _re_strip.sub(
                            r"^ACTION:\s*/.*$", "", raw, flags=_re_strip.MULTILINE
                        ).strip()
                        if cleaned != raw.strip():
                            _update_last_assistant(
                                cleaned or "(commands executed)",
                                role_label,
                            )
                        # No results → nothing to continue on
                        if not exec_results:
                            break
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
                            on_permission_denied=_on_permission_denied,
                            on_thinking=_on_thinking,
                            thinking_budget=_budget,
                            memory_context=_memory,
                        )
                        if chunks:
                            _update_last_assistant("".join(chunks), role_label)

                    # Show per-role cost
                    _role_cost = engine.cost_usd - _cost_before
                    _role_in = engine.token_usage["input"] - _in_before
                    _role_out = engine.token_usage["output"] - _out_before
                    if _role_in > 0 or _role_out > 0:
                        _cost_str = f"${_role_cost:.3f}" if _role_cost > 0 else ""
                        _append_system_message(
                            f"{role_label}: {_role_in:,} in / {_role_out:,} out"
                            f"{' · ' + _cost_str if _cost_str else ''}"
                            f" [{_effective_model}]"
                        )

                    if engine._stop_requested:
                        break

                    # --- Auto-advance logic ---
                    prev_role_id = engine.current_role

                    # --- Conditional skip: if agent says SKIP, advance ---
                    last_out = ""
                    for msg in reversed(engine.messages):
                        if msg["role"] == "assistant":
                            last_out = msg["content"]
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

                    # Session Manager: STOP and wait for user approval
                    # The SM is conversational — user must review the plan
                    # and explicitly approve before the pipeline continues.
                    if prev_role_id == "session_manager" and not _sm_approval:
                        _update_status()
                        _update_pipeline_display(engine)
                        break

                    # Dynamic routing: parse SM's "Skip agents" section
                    # and remove skipped roles from the route.
                    if prev_role_id == "session_manager" and _sm_approval:
                        _SKIPPABLE = {"critic_agent", "reviewer_agent", "runtime_agent"}
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

                    # --- Critic/Runtime findings gate ---
                    # Pause after Critic or Runtime Agent so the user can
                    # review findings and optionally skip some before Builder.
                    _REVIEW_ROLES = {"critic_agent", "runtime_agent"}
                    if (prev_role_id in _REVIEW_ROLES
                            and not state.get("_findings_approved")
                            and last_out.strip()):
                        # Check if there are actual findings (not just "no issues")
                        _no_issue_kw = ("no issues", "no findings", "lgtm",
                                        "looks good", "no critical", "SKIP")
                        _has_findings = not any(
                            k.lower() in last_out[:500].lower() for k in _no_issue_kw
                        )
                        if _has_findings:
                            _append_system_message(
                                f"--- {_format_role_label(prev_role_id)} review done. "
                                f"Reply 'go' to accept all findings, or describe "
                                f"which to skip (e.g. 'skip finding about X') ---"
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

                    # Dashboard / follow-up: never advance — keep conversation open
                    if mode_dropdown.value == "dashboard" or state.get("_follow_up"):
                        _update_status()
                        break

                    # Advance to next role
                    has_next = engine.advance_role()
                    if not has_next:
                        state.pop("_retry_used", None)
                        state.pop("_builder_retries", None)
                        # Acceptance gate: check if test agent approved
                        _cycle_verdict = _check_acceptance_gate(engine)
                        _append_system_message(
                            f"--- Cycle complete {_cycle_verdict} ---"
                        )
                        _update_pipeline_display(engine)
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
                state["streaming"] = False
                _set_working(False)
                _update_status()
                _update_button_states()
                _auto_save_session()
                _check_auto_compact()
                # Process next queued message if any
                _process_queue()

        threading.Thread(target=_worker, daemon=True).start()

    def _process_queue():
        """Send the next queued message, if any."""
        if state["message_queue"] and not state["streaming"]:
            next_msg = state["message_queue"].pop(0)
            _update_queue_display()
            input_textarea.value = next_msg
            _on_send(None)

    def _on_stop(button):
        engine = state["engine"]
        if engine:
            engine.request_stop()
            # Kill the persistent CLI process
            if hasattr(engine.client, "kill"):
                engine.client.kill()
        state["streaming"] = False
        _set_working(False)
        _update_button_states()
        _append_system_message("Generation stopped by user.")
        if state["message_queue"]:
            n = len(state["message_queue"])
            state["message_queue"].clear()
            _update_queue_display()
            _append_system_message(f"Cleared {n} queued message(s).")

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
        state["recent_edits"].clear()
        state.pop("_mode_suggested", None)
        state.pop("_pending_mode_msg", None)
        state.pop("_retry_used", None)
        state.pop("_builder_retries", None)
        state.pop("_follow_up", None)
        state["message_queue"].clear()
        state["session_start_time"] = None
        queue_html.value = ""
        undo_btn.disabled = True
        session_dropdown.value = ""
        _refresh_chat_html()
        _update_status()
        _update_button_states()

    def _on_advance_role(button):
        """Manual advance — fallback when auto-advance was stopped."""
        engine = state["engine"]
        if not engine:
            return
        prev_role = _format_role_label(engine.current_role)
        has_next = engine.advance_role()
        if has_next:
            next_role = _format_role_label(engine.current_role)
            _append_system_message(
                f"--- Manual handoff: {prev_role} \u2192 {next_role} ---"
            )
        else:
            state.pop("_retry_used", None)
            _append_system_message("--- Cycle complete ---")
        _update_status()
        _update_button_states()
        _auto_save_session()

    def _on_mode_change(change):
        new_mode = change["new"]
        # Update mode description label
        desc = _MODE_DESCRIPTIONS.get(new_mode, "")
        mode_desc_html.value = (
            f'<span style="color:#888; font-size:12px; margin-left:4px;">'
            f'{desc}</span>'
        )
        engine = state["engine"]
        if engine and not state["streaming"] and not engine.messages:
            engine.reset_cycle(mode=new_mode)
            _update_status()
        # Dashboard mode: lock permission to default, model to cheapest
        if new_mode == "dashboard":
            state["_perm_before_dashboard"] = perm_dropdown.value
            perm_dropdown.value = "ask_all"
            perm_dropdown.disabled = True
            state["_model_before_dashboard"] = model_dropdown.value
            _cheap = _PROVIDER_CHEAP.get(provider_dropdown.value, "haiku")
            model_dropdown.value = _cheap
        else:
            perm_dropdown.disabled = state.get("streaming", False)
            saved = state.pop("_perm_before_dashboard", None)
            if saved and perm_dropdown.value == "ask_all":
                perm_dropdown.value = saved
            saved_model = state.pop("_model_before_dashboard", None)
            _cheap = _PROVIDER_CHEAP.get(provider_dropdown.value, "haiku")
            if saved_model and model_dropdown.value == _cheap:
                model_dropdown.value = saved_model

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
        engine = state["engine"]
        if engine:
            state["engine"] = None
            cli_perm = _PROFILE_TO_CLI_PERM.get(new_profile, "default")  # CLI default
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
        """Revert the last file edit via git checkout."""
        if not state["recent_edits"]:
            return
        import subprocess as _sp
        last = state["recent_edits"].pop()
        fpath = last["file"]
        short = fpath.replace(str(ctx.repo_dir or ""), "").lstrip("/")
        try:
            result = _sp.run(
                ["git", "checkout", "--", fpath],
                capture_output=True, text=True,
                cwd=str(ctx.repo_dir or "."), timeout=10,
            )
            if result.returncode == 0:
                _append_system_message(f"\u21a9 Reverted: {short}")
            else:
                _append_system_message(
                    f"Undo failed: {result.stderr.strip()[:100]}"
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

    # -- wire events -------------------------------------------------------
    send_btn.on_click(_on_send)
    stop_btn.on_click(_on_stop)
    new_cycle_btn.on_click(_on_new_cycle)
    advance_btn.on_click(_on_advance_role)
    load_session_btn.on_click(_on_load_session)
    delete_session_btn.on_click(_on_delete_session)
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

    # -- initial state -----------------------------------------------------
    _update_status()
    _update_button_states()
    _refresh_session_dropdown()

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
})();
"""

    tab_widget = agent_content
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

    return (
        f'<div class="delfin-agent-status">'
        f'<span class="mode-badge">{_html.escape(mode)}</span>'
        f"{role_info}"
        f"{backend_info}"
        f"{perm_badge}"
        f'<span class="tokens-info">{tokens_str} · {cost_str}</span>'
        f"</div>"
    )


def _estimate_cost_str(
    backend: str, input_tokens: int, output_tokens: int,
    provider: str = "claude",
) -> str:
    """Rough cost string."""
    if provider == "kit":
        return "free (KIT)"
    if backend == "cli":
        return "included in subscription"
    if provider == "openai":
        cost = (input_tokens * 2.0 + output_tokens * 8.0) / 1_000_000
    else:
        cost = (input_tokens * 3.0 + output_tokens * 15.0) / 1_000_000
    return f"~${cost:.3f}"
