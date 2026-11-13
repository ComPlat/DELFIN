"""Whether the agent has a language model it can reach, and how to connect one.

A fresh install has none: the shipped default is the CLI backend, whose
program is either not installed or, as measured on a fresh install
(2026-11-13), an older release that answers DELFIN's command line with
"unknown option '--include-hook-events'". The terminal then failed the turn
with a backend error and the dashboard hid the agent tab. This module is the
one answer both ask instead:

* :func:`cli_program_problem` -- is the CLI backend's program there and new
  enough (its ``--help`` names the flag DELFIN passes)? Asked once per
  program file (path and modification time).
* :func:`selection_problem` -- can a (backend, provider) selection run: the
  key it needs, or its program? Mirrors the branches of
  ``api_client.create_client``. A local server (``ollama``) is not probed
  here; the dashboard probes it when it lists providers.
* :func:`connected_providers` -- the providers with a key or a usable program.
* :func:`no_model_message` / :func:`no_model_card_html` -- the text for the
  terminal and the card for the agent tab.

Keys are read through :func:`delfin.agent.credentials.load_credential`
(environment, keyring, ~/.delfin/credentials.json); values never leave it.

Rejected: falling back to whichever provider has a key. The provider is the
person's choice (cost, data leaving the machine); the message names what is
stored and how to choose it.
"""

from __future__ import annotations

import html
import os
import shutil
import subprocess

#: The flag of the CLI backend's command line (api_client.CLIClient) that an
#: older program does not know.
CLI_REQUIRED_FLAG = "--include-hook-events"
CLI_PROGRAM = "claude"

#: Key name, provider id, label -- for the providers a key connects.
KEYED_PROVIDERS = (
    ("KIT_TOOLBOX_API_KEY", "kit", "KIT Toolbox"),
    ("OPENAI_API_KEY", "openai", "OpenAI"),
    ("ANTHROPIC_API_KEY", "claude", "API"),
)

_HELP_CACHE: dict[tuple[str, int], str] = {}


def _key(name: str) -> str:
    try:
        from .credentials import load_credential
        return load_credential(name) or ""
    except Exception:
        return os.environ.get(name, "")


def cli_program_problem(path: str | None = None) -> str:
    """'' when the CLI backend's program can run DELFIN's command line, else why not."""
    path = path or shutil.which(CLI_PROGRAM)
    if not path:
        return f"the CLI backend is selected, but its program ({CLI_PROGRAM}) is not on PATH"
    try:
        stamp = (os.path.realpath(path), os.stat(path).st_mtime_ns)
    except OSError:
        return f"the CLI backend's program at {path} cannot be read"
    text = _HELP_CACHE.get(stamp)
    if text is None:
        try:
            done = subprocess.run([path, "--help"], capture_output=True, text=True,
                                  timeout=20, stdin=subprocess.DEVNULL)
            text = (done.stdout or "") + (done.stderr or "")
        except (OSError, subprocess.SubprocessError):
            text = ""
        _HELP_CACHE[stamp] = text
    if CLI_REQUIRED_FLAG not in text:
        return (f"the CLI backend's program at {path} is too old for DELFIN "
                f"(its --help has no {CLI_REQUIRED_FLAG}); update it, or connect a model "
                "another way")
    return ""


def selection_problem(backend: str, provider: str, model: str = "") -> str:
    """'' when this backend/provider can run, else why not (no secrets in it)."""
    backend = (backend or "").strip() or "api"
    provider = (provider or "").strip()
    if provider == "ollama":
        return ""
    if provider == "kit":
        return "" if _key("KIT_TOOLBOX_API_KEY") else \
            "the KIT Toolbox is selected, but no KIT_TOOLBOX_API_KEY is stored"
    if provider == "openai":
        if backend == "cli":
            return "" if shutil.which("codex") else \
                "the OpenAI CLI backend is selected, but its program (codex) is not on PATH"
        return "" if _key("OPENAI_API_KEY") else \
            "OpenAI is selected, but no OPENAI_API_KEY is stored"
    if provider not in ("", "claude", "anthropic"):
        return ""
    if backend == "api":
        return "" if _key("ANTHROPIC_API_KEY") else \
            "the API backend is selected, but no key for it is stored"
    return cli_program_problem()


def connected_providers() -> list[str]:
    """Provider ids with a stored key or a usable program, in KEYED_PROVIDERS order."""
    found = [provider for name, provider, _ in KEYED_PROVIDERS if _key(name)]
    if "claude" not in found and shutil.which(CLI_PROGRAM) and not cli_program_problem():
        found.append("claude")
    if "openai" not in found and shutil.which("codex"):
        found.append("openai")
    return found


def _stored_note() -> str:
    stored = [(provider, label) for name, provider, label in KEYED_PROVIDERS[:2] if _key(name)]
    if not stored:
        return ""
    provider, _ = stored[0]
    example = ("--provider kit --model kit.MODEL" if provider == "kit"
               else "--provider openai --model MODEL")
    return (f"A key for {' and '.join(label for _, label in stored)} is stored: choose that "
            "provider and a model in the dashboard's DELFIN Agent tab (the terminal uses "
            f"the same choice), or start delfin-agent {example}.")


def no_model_message(problem: str) -> str:
    """The terminal's answer when no model can be reached."""
    lines = [f"No language model is set up for the agent: {problem}." if problem
             else "No language model is set up for the agent.",
             "",
             "Connect one, then start delfin-agent again:",
             "  KIT Toolbox   delfin-agent credentials set KIT_TOOLBOX_API_KEY",
             "                (the key from https://ki-toolbox.scc.kit.edu)",
             "  OpenAI        delfin-agent credentials set OPENAI_API_KEY",
             "  local server  an OpenAI-compatible server (Ollama, vLLM, LM Studio);",
             "                OLLAMA_HOST if it is not at http://localhost:11434,",
             "                then delfin-agent --provider ollama --model NAME",
             "The key is read without echo and stored in the system keyring or in",
             "~/.delfin/credentials.json (mode 0600). Then choose the provider and the",
             "model in the dashboard's DELFIN Agent tab; the terminal uses that choice."]
    note = _stored_note()
    if note:
        lines += ["", note]
    return "\n".join(lines)


def no_model_card_html(problem: str) -> str:
    """The agent tab's card when no provider is connected."""
    esc = html.escape
    detail = (f'<div style="margin-top:6px;opacity:.8">Found: {esc(problem)}.</div>'
              if problem else "")
    return (
        '<div class="delfin-agent-no-model" style="border:1px solid '
        'var(--jp-border-color1,#c8cdd3);border-radius:8px;padding:12px 16px;'
        'margin:8px 0;line-height:1.5">'
        '<div style="font-weight:600;margin-bottom:4px">No language model is connected yet</div>'
        '<div>The agent needs one. Choose one way, then reload this page:</div>'
        '<ul style="margin:6px 0 6px 18px;padding:0">'
        '<li><b>KIT Toolbox</b>: in a terminal on this machine, '
        '<code>delfin-agent credentials set KIT_TOOLBOX_API_KEY</code> and paste the key '
        'from <a href="https://ki-toolbox.scc.kit.edu" target="_blank" rel="noopener">'
        'ki-toolbox.scc.kit.edu</a>.</li>'
        '<li><b>OpenAI</b>: <code>delfin-agent credentials set OPENAI_API_KEY</code>.</li>'
        '<li><b>A local model server</b> (Ollama, vLLM, LM Studio): start it; set '
        '<code>OLLAMA_HOST</code> if it is not at <code>http://localhost:11434</code>.</li>'
        '</ul>'
        '<div>The key is read without echo and kept in the system keyring or in '
        '<code>~/.delfin/credentials.json</code> (mode 0600). After the reload, choose the '
        'provider and the model above the chat.</div>'
        f'{detail}</div>'
    )
