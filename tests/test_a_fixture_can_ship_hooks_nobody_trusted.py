"""A workspace that ships hook definitions nobody has trusted.

That is what a freshly cloned repository is, and it is the situation a
user is in when they write a hook, restart, and nothing happens. The
mechanism side is covered by
``test_a_checked_out_repository_cannot_run_commands.py``; this file
covers the fixture the behaviour task needs, and the two properties that
make it a fixture rather than a leak.

The definition cannot be committed inside the fixture: ``.delfin/`` is
ignored checkout-wide, so a settings.json written there would exist only
in the working copy that wrote it and be missing from every clone — a
task that passes locally and measures nothing in CI. It is therefore
carried under tests/fixtures/hooks_seed/ and installed at run time,
inside the pristine guard.
"""

import json
from pathlib import Path

from delfin.agent import hooks as H
from delfin.agent.benchmark_runner import (_PristineWorkspace,
                                           _seed_fixture_hooks,
                                           workspace_for)

_REPO = Path(__file__).resolve().parents[1]
_SEEDS = _REPO / "tests" / "fixtures" / "hooks_seed"


def test_the_seed_directory_is_named_after_its_fixture():
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    assert ws is not None
    assert (_SEEDS / ws.name).is_dir(), sorted(p.name for p in _SEEDS.iterdir())


def test_the_shipped_hook_runs_nothing_that_matters():
    """A directory someone trusts by hand must not then execute something
    with an effect. The fixture is about the trust decision, not about
    what a hook does."""
    for f in sorted(_SEEDS.rglob("*.json")):
        data = json.loads(f.read_text(encoding="utf-8"))
        for entries in (data.get("hooks") or {}).values():
            for entry in entries:
                for hook in entry.get("hooks") or []:
                    cmd = str(hook.get("command", ""))
                    assert cmd.startswith("echo "), f"{f.name}: {cmd!r}"


def test_the_definition_is_not_ignored_by_git():
    """The whole reason the seed exists. If this path ever lands under a
    `.delfin/` directory it becomes invisible to a clone again."""
    assert ".delfin" not in _SEEDS.parts
    for f in _SEEDS.rglob("*.json"):
        assert ".delfin" not in f.relative_to(_REPO).parts


def test_the_seeded_workspace_offers_hooks_and_loads_none():
    """Installed, seen, and withheld — with the user told how many and
    what to type."""
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    with _PristineWorkspace(_REPO):
        assert _seed_fixture_hooks(_REPO, ws) >= 1
        assert (ws / ".delfin" / "settings.json").is_file()
        cfg = H.load_hooks(ws)
        assert cfg.is_empty(), "an untrusted workspace supplied a hook"
        notice = " ".join(cfg.warnings)
        assert "hook" in notice.lower()
        assert "/hooks trust" in notice


def test_the_seed_does_not_survive_the_run():
    """Inside the guard, like the memory seed. A benchmark may write what
    it likes; it does not get to leave it in the checkout."""
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    with _PristineWorkspace(_REPO):
        _seed_fixture_hooks(_REPO, ws)
        assert (ws / ".delfin" / "settings.json").is_file()
    assert not (ws / ".delfin" / "settings.json").exists()
