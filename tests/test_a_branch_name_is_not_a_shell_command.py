"""A branch name the wave is told to run against must be a branch name.

``compose`` fills the spec's ``base`` into shell commands that every member is
instructed to run -- "BEFORE ANYTHING ELSE: ``git log --oneline {base}..HEAD |
wc -l``" in the reviewer brief, and "``git log --oneline {base}..<branch>``" in
the integrator's control step. Nothing validated ``base``: a spec carrying
``origin/main && echo INJECTED`` produced exactly that line in the mandated
text, so a string the wave did not author became a command all four reviewers
were told to execute first.

Task ids are already validated this way (``wave._ID``, ``compose._validate_packages``);
``base`` is the one interpolation into a mandated command that was not. The
fix is the same shape: check it once, before any brief is built, and refuse
with a message that names the offending value.

Input: a spec whose ``base`` carries shell text. Output: ``WaveError`` rather
than a brief that instructs the shell to run it.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm import compose, wave


def _spec(**over) -> dict:
    """The smallest spec that composes, so only ``base`` varies per case."""
    base = dict(
        run="w48",
        packages=[{"id": "p1", "title": "one", "body": "do the thing"}],
        arms=[{"id": "a1", "title": "arm", "body": "work"}],
        integrator={"id": "p1-int", "title": "int", "body": "integrate"},
        lander={"id": "zz-land", "title": "land", "body": "land it"},
        reviewer={"id": "zz-rev", "title": "review", "body": "review it"},
        fixer={"id": "zz-fix", "title": "fix", "body": "fix it"},
    )
    base.update(over)
    return base


def _briefs(tasks) -> str:
    return "\n".join(str(t.get("brief") or "") for t in tasks)


def test_a_base_carrying_shell_text_is_refused():
    """The documented defect: ``base`` used to reach the mandated command."""
    with pytest.raises(wave.WaveError) as exc:
        compose.compose(**_spec(base="origin/main && echo INJECTED"))
    assert "base" in str(exc.value).lower(), (
        "the refusal must name what was refused, not just fail")


@pytest.mark.parametrize("bad", [
    "origin/main; echo pwned",
    "origin/main | tee /tmp/x",
    "origin/main`id`",
    "origin/main$(id)",
    "origin/main\nrm -rf build",
    "origin main",
])
def test_every_shell_metacharacter_is_refused(bad):
    """One case per metacharacter shape, so a partial fix is still red."""
    with pytest.raises(wave.WaveError):
        compose.compose(**_spec(base=bad))


@pytest.mark.parametrize("good", [
    "origin/main",
    "origin/main",
    "main",
    "HEAD",
    "41e15fe9f",
    "origin/wave/w48/zz-land",
    "wave/w48/h1-int",
    "release-2.0",
    "origin/feature/a_b.c",
])
def test_a_real_branch_name_still_composes(good):
    """The never-worse case: every form the product actually uses still works."""
    tasks = compose.compose(**_spec(base=good))
    assert tasks, "a valid base must still produce the task list"
    assert good in _briefs(tasks), (
        f"the validated base must still reach the briefs; {good!r} vanished")


def test_the_refusal_happens_before_any_brief_is_built():
    """Refuse once, up front -- not midway through composing a partial wave."""
    with pytest.raises(wave.WaveError):
        compose.compose(**_spec(base="origin/main; rm -rf /"))
