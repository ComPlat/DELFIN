"""The prompt says one thing about language, not three.

A live GLM turn answered the greeting "Hallo" after 190 seconds, and its
own reasoning names why it had to think about it at all:

    "The critical rules say communicate in German, but the session
     language section says English."

Three parts of the built prompt spoke about language and two of them were
hardcoded: the critical anchor ("Communicate with the user in German") and
the shared answer-language rule ("answer in the language of the user's
LATEST message"). The session pin -- set from the user's own first message
-- was the only one that knew the session, and it was outranked by an
anchor that claims to override any conflicting prior instruction.

These build the real prompt and read what is in it.
"""

import pytest

from delfin.agent.prompt_loader import PromptLoader

_ROLES = ("solo_agent", "dashboard_agent", "office_agent", "builder_agent")


def _sections(role: str, lang: str = "") -> dict[str, str]:
    secs = PromptLoader().compose_sections(
        role_id=role, mode_id="default", task_text="read calc.py",
        session_language=lang)
    return {s.name: s.content for s in secs}


@pytest.mark.parametrize("role", _ROLES)
def test_an_english_session_is_told_english_everywhere(role):
    secs = _sections(role, "en")
    for name in ("critical_anchor", "honesty_addendum"):
        text = secs[name]
        assert "in English" in text, (role, name)
        assert "in German" not in text, (role, name, text)


@pytest.mark.parametrize("role", _ROLES)
def test_a_german_session_is_told_german_everywhere(role):
    secs = _sections(role, "de")
    for name in ("critical_anchor", "honesty_addendum"):
        assert "in German" in secs[name], (role, name)


@pytest.mark.parametrize("role", _ROLES)
def test_a_pinned_session_carries_no_latest_message_rule(role):
    """The two rules disagree exactly where it matters -- one German
    sentence inside an English session -- so the pinned session must not
    carry both."""
    joined = "\n".join(_sections(role, "en").values())
    assert "language of the user's LATEST message" not in joined, role


@pytest.mark.parametrize("role", _ROLES)
@pytest.mark.parametrize("lang", ("", "en", "de"))
def test_code_stays_english_in_every_session(role, lang):
    joined = "\n".join(_sections(role, lang).values())
    assert "INTO code stays English" in joined, (role, lang)
    assert "docstrings" in joined, (role, lang)


@pytest.mark.parametrize("role", _ROLES)
def test_an_unpinned_session_reads_exactly_what_it_read_before(role):
    """No session language, no rewrite: the default wording is the file's."""
    secs = _sections(role, "")
    assert "language of the user's LATEST message" in secs["honesty_addendum"]
    assert ("Communicate with the user in German."
            in secs["critical_anchor"]), role


def test_the_rewrite_leaves_other_text_alone():
    loader = PromptLoader()
    text = "- **Report faithfully.** Say what failed.\n"
    assert loader._pin_language_rule(text, "en") == text
    assert loader._pin_language_rule("", "en") == ""


def test_the_rewrite_keeps_the_bullets_around_it():
    loader = PromptLoader()
    text = (
        "- **Before.** Keep me.\n"
        "- **Answer in the language of the user's LATEST message** -- German\n"
        "  in, German out; English in, English out.\n"
        "- **After.** Keep me too.\n"
    )
    out = loader._pin_language_rule(text, "en")
    assert "- **Before.** Keep me.\n" in out
    assert "- **After.** Keep me too.\n" in out
    assert "LATEST message" not in out
    assert "Answer in English" in out


def test_an_unknown_language_changes_nothing():
    loader = PromptLoader()
    text = "- **Answer in the language of the user's LATEST message** -- x\n"
    assert loader._pin_language_rule(text, "fr") == text
    assert loader._pin_language_rule(text, "") == text


def test_the_session_block_names_its_own_rank(tmp_path):
    """A remembered preference sits in the memory section saying German is
    the default. The block that pins the session says which one wins."""
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng._session_language = "en"
    block = eng._session_language_block()
    assert "English" in block
    assert "outranks" in block
    assert "remembered language preference" in block


def test_no_session_language_means_no_block():
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng._session_language = ""
    assert eng._session_language_block() == ""


# --- the learned rules -------------------------------------------------

def test_the_shipped_language_rules_are_the_ones_dropped():
    """Read the file that ships, not a synthetic string. Three of its
    rules name a conversation language; those are the ones that compete
    with the session pin."""
    import json
    from pathlib import Path

    from delfin.agent.provider_profile import _dictates_a_language

    path = (Path(__file__).resolve().parents[1]
            / "delfin" / "agent" / "learned_profiles.json")
    data = json.loads(path.read_text())
    dropped, kept = [], []
    for prov, entry in data.items():
        if not isinstance(entry, dict):
            continue
        for rule in (entry.get("communication") or {}).get("rules", []):
            (dropped if _dictates_a_language(rule) else kept).append(
                (prov, rule))
    assert [p for p, _ in dropped] == ["shared", "kit", "ollama"], dropped
    for _, rule in dropped:
        assert "German for user-facing conversation" in rule
    # Everything else survives: exactness, units, tool discipline.
    assert len(kept) >= 15, len(kept)
    for _, rule in kept:
        assert "German for user-facing" not in rule


def test_a_rule_about_code_language_is_kept():
    from delfin.agent.provider_profile import _dictates_a_language
    assert not _dictates_a_language("Code, commits and artifacts in English")
    assert not _dictates_a_language("Docstrings in English")
    assert _dictates_a_language("Answer the user in German")


def test_a_rule_with_no_language_is_kept():
    from delfin.agent.provider_profile import _dictates_a_language
    assert not _dictates_a_language("Answer exactly, with numbers and units")
    assert not _dictates_a_language("")
    assert not _dictates_a_language(None)


def test_the_pinned_language_does_not_move_the_cacheable_prefix():
    """Measured on the KIT deployment, 2026-09-07: the same 15k-token
    DELFIN prompt answers in 2.8 s when the endpoint can serve the prefix
    from its cache and 341 s when it cannot. The session language is set
    from the first message and never changes, so every turn of a session
    must produce the same stable head."""
    loader = PromptLoader()
    kwargs = dict(role_id="solo_agent", mode_id="default",
                  task_text="read calc.py", session_language="en")
    first = loader.stable_prefix(**kwargs)
    second = loader.stable_prefix(**kwargs)
    assert first == second
    # A different session language is a different session, not a
    # different turn -- these may and must differ.
    assert first != loader.stable_prefix(**{**kwargs,
                                            "session_language": "de"})
