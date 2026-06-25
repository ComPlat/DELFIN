"""The model dropdown must only offer models that work for the current auth.

OpenAI via the Codex CLI on a ChatGPT-account login (no OPENAI_API_KEY) cannot
use the API-only '*-codex' models — they 400 (bug 2026-06-25, ka_xn0397/Jerome:
gpt-5.3-codex → "not supported when using Codex with a ChatGPT account"). When
an API key IS present (API mode), all models are offered.
"""

from pathlib import Path

_SRC = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
        / "tab_agent.py").read_text(encoding="utf-8")


def test_filter_is_wired_for_codex_chatgpt_mode():
    assert '_codex_cli_available and not _provider_key("OPENAI_API_KEY")' in _SRC
    assert '"-codex" not in str(mid)' in _SRC


def test_filter_predicate_drops_only_codex_models():
    models = [("GPT-5.4", "gpt-5.4"), ("GPT-5.3-codex", "gpt-5.3-codex"),
              ("GPT-5.1-codex-mini", "gpt-5.1-codex-mini"), ("GPT-4.1", "gpt-4.1")]
    filtered = [(l, m) for (l, m) in models if "-codex" not in str(m)]
    ids = {m for _, m in filtered}
    assert "gpt-5.3-codex" not in ids        # API-only model dropped
    assert "gpt-5.1-codex-mini" not in ids
    assert "gpt-5.4" in ids                   # ChatGPT-compatible kept
    assert "gpt-4.1" in ids
