"""A permission the user granted can be withdrawn through the product.

Input:  a stand-in settings file (``kit_settings.USER_SETTINGS_PATH``
        redirected into tmp_path, because that constant is resolved at
        import and a late HOME would arrive too late), a real
        ``KitToolPermissions``, and the real ``AgentEngine`` methods bound
        to it the way the engine binds them. The typed route is driven
        through ``repl_commands.dispatch``, the line a person types.
Output: whether a grant leaves BOTH the stored settings and the running
        session, whether it is visible before it is withdrawn, and what
        happens to a rule nobody granted or that nobody chose.

Semantics: ``persist_kit_pattern`` writes a bash rule to
``~/.delfin/settings.json`` AND appends it to
``perms.bash_auto_allow_patterns``; ``add_kit_read_dir`` writes a
directory to the same file AND appends it to ``perms.session_read_dirs``.
The gate reads the live tuples, not the file. A revoke that rewrote only
the file would therefore report success while the session kept the grant
until it ended -- which is why every case here asserts both halves.

Measured on 2026-11-22: ``kit_settings.remove_pattern`` and
``remove_read_dir`` existed, were correct, and had no caller anywhere,
while ``persist_pattern`` had 11 and ``persist_read_dir`` 2. The engine
offered ``remove_kit_workspace_dir`` for extra dirs and no counterpart
for the other two, and no command showed the granted entries at all. A
rule approved once with "always allow" could be taken back only by
editing the JSON by hand.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_a_granted_permission_can_be_taken_back.py -q
"""

from __future__ import annotations

import json
import types
from pathlib import Path

import pytest

PATTERN = r"^\s*pytest\b"


@pytest.fixture
def store(tmp_path, monkeypatch):
    """Redirect the user settings file; never the real ~/.delfin."""
    from delfin.agent import kit_settings
    path = tmp_path / "settings.json"
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH", path)
    return path


def _stored(path: Path, key: str) -> list:
    if not path.exists():
        return []
    raw = json.loads(path.read_text(encoding="utf-8"))
    for block in raw.values():
        if isinstance(block, dict) and key in block:
            return list(block[key])
    return list(raw.get(key, []))


@pytest.fixture
def engine(tmp_path):
    """The real engine methods, bound as the engine binds them.

    ``AgentEngine.__new__`` without ``__init__`` is how
    tests/test_allow_for_good_covers_what_was_asked.py drives these same
    methods: the constructor builds a client and loads the prompt pack,
    neither of which decides anything here. The permissions object and
    the settings file are real.
    """
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.engine import AgentEngine

    workspace = tmp_path / "ws"
    workspace.mkdir()
    perms = KitToolPermissions(workspace=str(workspace), mode="default")
    eng = AgentEngine.__new__(AgentEngine)
    eng.client = types.SimpleNamespace(_permissions=perms,
                                       permission_mode="default")
    eng.repo_dir = str(workspace)
    eng.told: list[str] = []
    eng._tell_the_running_turn = eng.told.append
    return eng


def _grantable_dir(tmp_path) -> Path:
    d = tmp_path / "readable"
    d.mkdir()
    return d.resolve()


# -- the pattern --------------------------------------------------------


def test_a_withdrawn_pattern_leaves_the_file_and_the_session(engine, store):
    """The decisive case: both halves, in both directions."""
    from delfin.agent.engine import AgentEngine

    ok, _ = AgentEngine.persist_kit_pattern(engine, PATTERN)
    assert ok
    assert PATTERN in _stored(store, "allow_patterns")
    assert PATTERN in engine.kit_permissions.bash_auto_allow_patterns

    ok, msg = AgentEngine.forget_kit_pattern(engine, PATTERN)
    assert ok, msg
    assert PATTERN not in _stored(store, "allow_patterns"), "stored rule kept"
    assert PATTERN not in engine.kit_permissions.bash_auto_allow_patterns, (
        "the session still allows it: a revoke that reports success and "
        "leaves the gate open for the rest of the run")


def test_the_running_turn_is_told_the_rule_is_gone(engine, store):
    """The grant tells the turn; so must the revoke.

    An agent that was told "it no longer needs a confirmation" will keep
    using the command. If the withdrawal is silent, its next call is
    refused for a reason it was never given.
    """
    from delfin.agent.engine import AgentEngine

    AgentEngine.persist_kit_pattern(engine, PATTERN)
    engine.told.clear()
    AgentEngine.forget_kit_pattern(engine, PATTERN)
    assert engine.told, "the running turn was not told the rule was withdrawn"
    said = " ".join(engine.told)
    # The notice quotes the pattern with !r, as the grant's notice does,
    # so comparing against the raw string would fail over the escaped
    # backslashes rather than over a missing name.
    assert repr(PATTERN)[1:-1] in said, said
    low = said.lower()
    assert "withdrew" in low or "no longer" in low or "removed" in low


def test_a_rule_nobody_granted_is_not_a_success(engine, store):
    """Not silently 'removed': there was nothing to remove."""
    from delfin.agent.engine import AgentEngine

    ok, msg = AgentEngine.forget_kit_pattern(engine, r"^\s*never-granted\b")
    assert not ok
    assert "not" in msg.lower()


def test_a_pattern_delfin_ships_is_refused_not_dropped(engine, store):
    """A shipped default is in no settings file, so its removal could not
    be recorded: the next session would allow it again. Refusing says so;
    dropping it from the live tuple would be a revoke that does not hold.
    """
    from delfin.agent.api_client import _BUILTIN_BASH_AUTO_ALLOW
    from delfin.agent.engine import AgentEngine

    builtin = sorted(_BUILTIN_BASH_AUTO_ALLOW)[0]
    engine.kit_permissions.bash_auto_allow_patterns = (
        tuple(engine.kit_permissions.bash_auto_allow_patterns) + (builtin,)
        if builtin not in engine.kit_permissions.bash_auto_allow_patterns
        else engine.kit_permissions.bash_auto_allow_patterns)

    ok, msg = AgentEngine.forget_kit_pattern(engine, builtin)
    assert not ok
    assert builtin in engine.kit_permissions.bash_auto_allow_patterns
    assert "ship" in msg.lower() or "default" in msg.lower()


def test_a_deny_the_user_wrote_can_also_be_withdrawn(engine, store):
    """Symmetry with persist_kit_pattern, which takes the same kind."""
    from delfin.agent.engine import AgentEngine

    deny = r"^\s*curl\b"
    assert AgentEngine.persist_kit_pattern(engine, deny, kind="deny")[0]
    assert deny in engine.kit_permissions.bash_deny_patterns

    ok, msg = AgentEngine.forget_kit_pattern(engine, deny, kind="deny")
    assert ok, msg
    assert deny not in _stored(store, "deny_patterns")
    assert deny not in engine.kit_permissions.bash_deny_patterns


def _repo_records(engine, key: str, value: str) -> None:
    """Put an entry in the REPO settings, which load() merges in."""
    path = Path(engine.repo_dir) / ".delfin" / "settings.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"kit": {key: [value]}}), encoding="utf-8")


def test_a_deny_the_repo_records_is_not_called_withdrawn_for_good(
        engine, store):
    """The half-withdrawal this method exists to prevent.

    ``load()`` is the merged view and ``remove_pattern`` at user scope
    rewrites only the user file. ``kit_settings._merge`` lets a repository
    TIGHTEN and never widen, so this can only happen to the deny list --
    which does take the repository's entries. The deny therefore stays in
    force and returns on the next start; reporting that as a withdrawal
    full stop would be the same kind of half-truth as a save that
    deletes.
    """
    from delfin.agent.engine import AgentEngine

    deny = r"^\s*curl\b"
    _repo_records(engine, "deny_patterns", deny)
    engine.kit_permissions.bash_deny_patterns = (
        tuple(engine.kit_permissions.bash_deny_patterns) + (deny,))

    ok, msg = AgentEngine.forget_kit_pattern(engine, deny, kind="deny")
    assert ok, msg
    assert deny not in engine.kit_permissions.bash_deny_patterns
    assert "this session only" in msg.lower(), msg
    assert "next start" in msg.lower(), msg
    assert ".delfin/settings.json" in msg, msg


def test_an_allow_rule_a_repo_ships_was_never_in_force(engine, store):
    """The other side of the same merge rule, which is why the message
    above names only the deny list: an ``allow_patterns`` entry in a
    repository's settings widens, so it is never read at all. A test that
    asserted the half-withdrawal message for an allow rule would be
    asserting a state the product cannot reach."""
    from delfin.agent import kit_settings

    _repo_records(engine, "allow_patterns", PATTERN)
    merged = kit_settings.load(repo_dir=engine.repo_dir)
    assert PATTERN not in merged.allow_patterns, merged.allow_patterns


# -- the read directory -------------------------------------------------


def test_a_withdrawn_read_dir_leaves_the_file_and_the_session(
        engine, store, tmp_path):
    from delfin.agent.engine import AgentEngine

    target = _grantable_dir(tmp_path)
    ok, msg = AgentEngine.add_kit_read_dir(engine, target)
    assert ok, msg
    assert str(target) in _stored(store, "read_dirs")
    assert target in engine.kit_permissions.session_read_dirs

    ok, msg = AgentEngine.remove_kit_read_dir(engine, target)
    assert ok, msg
    assert str(target) not in _stored(store, "read_dirs")
    assert target not in engine.kit_permissions.session_read_dirs, (
        "the session can still read there")


def test_a_read_dir_nobody_granted_is_not_a_success(engine, store, tmp_path):
    from delfin.agent.engine import AgentEngine

    ok, msg = AgentEngine.remove_kit_read_dir(engine, _grantable_dir(tmp_path))
    assert not ok
    assert "not" in msg.lower()


# -- the typed route ----------------------------------------------------


def _ctx(engine, workspace):
    from delfin.agent.repl_commands import ReplContext
    return ReplContext(engine=engine, workspace=Path(workspace))


def test_the_permissions_command_shows_what_was_granted(
        engine, store, tmp_path):
    """Invisibility was half the defect: a rule nobody can see is a rule
    nobody can decide to withdraw."""
    from delfin.agent import repl_commands
    from delfin.agent.engine import AgentEngine

    target = _grantable_dir(tmp_path)
    AgentEngine.persist_kit_pattern(engine, PATTERN)
    AgentEngine.add_kit_read_dir(engine, target)

    out = repl_commands.dispatch("/permissions", _ctx(engine, tmp_path)).output
    assert PATTERN in out, out
    assert str(target) in out, out


def test_the_permissions_command_takes_a_pattern_back(
        engine, store, tmp_path):
    from delfin.agent import repl_commands
    from delfin.agent.engine import AgentEngine

    AgentEngine.persist_kit_pattern(engine, PATTERN)
    res = repl_commands.dispatch(f"/permissions forget {PATTERN}",
                                 _ctx(engine, tmp_path))
    assert res.handled
    assert PATTERN not in _stored(store, "allow_patterns"), res.output
    assert PATTERN not in engine.kit_permissions.bash_auto_allow_patterns


def test_the_permissions_command_takes_a_read_dir_back(
        engine, store, tmp_path):
    from delfin.agent import repl_commands
    from delfin.agent.engine import AgentEngine

    target = _grantable_dir(tmp_path)
    AgentEngine.add_kit_read_dir(engine, target)
    res = repl_commands.dispatch(f"/permissions unread {target}",
                                 _ctx(engine, tmp_path))
    assert res.handled
    assert str(target) not in _stored(store, "read_dirs"), res.output
    assert target not in engine.kit_permissions.session_read_dirs


def test_a_withdrawal_word_is_never_read_as_a_mode(engine, store, tmp_path):
    """``/permissions <word>`` switches the approval mode. The new words
    have to be taken first, or `forget` would be handed to
    set_kit_permission_mode as a mode name."""
    from delfin.agent import repl_commands

    switched = []
    engine.set_kit_permission_mode = switched.append
    for line in ("/permissions forget " + PATTERN,
                 "/permissions unread /tmp",
                 "/permissions granted"):
        repl_commands.dispatch(line, _ctx(engine, tmp_path))
    assert switched == [], switched
