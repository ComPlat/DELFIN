"""What a workspace ships and what is actually in force.

A repository can carry hook commands and MCP servers under ``.delfin/``,
and neither is loaded until the user has trusted that directory. The
refusal is right and thoroughly covered
(``test_a_checked_out_repository_cannot_run_commands.py``). What was not
covered is that nobody is told inside a turn: the withholding is rendered
by ``/hooks`` and ``/trust``, which a user opens only once they already
suspect it, and the agent has no tool that can see it at all.

Measured 2026-09-08, asked why a PreToolUse hook in ``.delfin/settings.json``
never fires: the answer was that there is no hook mechanism and the schema
does not fit — reasoned confidently from the one thing the prompt does say
about that file, which is that it holds permission rules. A hook that
stopped running looks exactly like a hook that had nothing to say, and so
does a hook that never started.

So the prompt names it, in the short form, and only while there is
something to name -- in its own section, because the session environment
block is composed for ``solo_agent`` alone and a workspace ships what it
ships whichever role is reading it.
"""

import json

import pytest

from delfin.agent import workspace_trust as WT
from delfin.agent.prompt_loader import PromptLoader


def _ships_hooks(tmp_path):
    ws = tmp_path / "checked-out-repo"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(json.dumps({
        "hooks": {"PreToolUse": [{
            "matcher": "write_file",
            "hooks": [{"type": "command", "command": "echo hi"}],
        }]}
    }), encoding="utf-8")
    return ws


def _env(ws) -> str:
    loader = PromptLoader()
    loader.workspace_root = ws
    return loader._withheld_config_block()


def test_the_withheld_hook_is_named_with_the_way_out(tmp_path):
    env = _env(_ships_hooks(tmp_path))
    assert "hook" in env.lower(), env
    assert "/hooks trust" in env, env


def test_a_workspace_that_ships_nothing_pays_nothing(tmp_path):
    """The line has to be free in the ordinary case, or it is a tax on
    every session for a situation almost none of them are in."""
    plain = tmp_path / "plain"
    plain.mkdir()
    assert _env(plain) == ""


def test_a_trusted_workspace_says_nothing_either(tmp_path):
    """Nothing is being withheld, so there is nothing to report — and a
    line that stays after the user has acted teaches them to ignore it."""
    ws = _ships_hooks(tmp_path)
    assert _env(ws) != ""
    WT.trust_workspace(ws, actor=WT.ACTOR_USER)
    assert _env(ws) == ""


def test_it_is_the_short_form(tmp_path):
    """The full notice is three sentences of explanation. This one rides
    along in every prompt while it lasts, so it carries the count, the
    reason and the command, and stops."""
    lines = _env(_ships_hooks(tmp_path)).splitlines()
    assert len(lines) == 1, lines
    assert len(lines[0]) < 160, lines[0]


def test_a_trust_layer_that_cannot_answer_does_not_take_the_prompt_down(
        tmp_path, monkeypatch):
    """Orientation is not worth a failed turn."""
    def _boom(*a, **k):
        raise RuntimeError("store unreadable")
    monkeypatch.setattr(WT, "pending_notices", _boom)
    ws = _ships_hooks(tmp_path)
    loader = PromptLoader()
    loader.workspace_root = ws
    assert loader._withheld_config_block() == ""
    assert "cwd:" in loader._build_session_env_block()


def test_the_short_form_still_reaches_the_reader_who_asked_for_all_of_it():
    """`short=` is an option, not a replacement: /trust shows the whole
    explanation because the user opened it on purpose."""
    import inspect
    sig = inspect.signature(WT.pending_notices)
    assert sig.parameters["short"].default is False

def test_it_survives_into_the_prompt_the_model_is_given(tmp_path):
    """A rule is a rule when it is IN the composed prompt. The helper
    returning the line proves only that the helper works."""
    ws = _ships_hooks(tmp_path)
    loader = PromptLoader()
    loader.workspace_root = ws
    for role, mode in (("solo_agent", "solo"),
                       ("dashboard_agent", "dashboard"),
                       ("office_agent", "office")):
        prompt = loader.build_system_prompt(
            role, mode, task_text="Warum feuert mein Hook nicht?")
        assert "/hooks trust" in prompt, (
            f"{role}: the withheld line never reached the model")


def test_a_tracker_that_drops_everything_cannot_drop_this(tmp_path):
    """Sections earn their place by being cited, and one that is never
    cited gets dropped. This is not that kind of section: it states a
    fact about the workspace the model has no other way to learn, and it
    is worth the most in the turns where nobody quotes it back."""
    class _DropsEverything:
        def should_skip(self, *a, **k):
            return True

    loader = PromptLoader()
    loader.workspace_root = _ships_hooks(tmp_path)
    loader._context_tracker = _DropsEverything()
    prompt = loader.build_system_prompt(
        "solo_agent", "solo", task_text="Warum feuert mein Hook nicht?")
    assert "/hooks trust" in prompt
