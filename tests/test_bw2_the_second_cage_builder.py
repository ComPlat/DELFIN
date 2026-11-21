"""BW2: the BW1 guarantee reaches every builder that hands an argv to bwrap.

Input:  the two argv builders in delfin/agent (``sandbox._bwrap_argv`` and
        the one in ``contained_run``), the mask list in
        ``sandbox.home_secret_paths``, and the tool-level path deny globs.
Output: whether a mask/bind target named in a finished argv is guaranteed
        to exist when bwrap mounts it, and whether any home secret is
        protected by the cage mask alone.
Semantics: BW1 closed this for ``contained_run`` and left
        ``ensure_argv_bind_targets`` with exactly ONE caller;
        ``sandbox._bwrap_argv`` handed its argv straight to bwrap. The
        leaf of a ``--tmpfs``/``--ro-bind`` whose parent chain sits under
        the read-only ``--ro-bind / /`` root cannot be created by bwrap,
        so a target that is missing at exec kills the command before it
        runs. Nothing here reads or writes the real HOME: every path is
        inside tmp_path, and the real ``~/.delfin`` is never touched.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_bw2_the_second_cage_builder.py -q
"""

from __future__ import annotations

import ast
import fnmatch
import shutil
import subprocess
from pathlib import Path

import pytest

from delfin.agent import cage_binds

REPO_ROOT = Path(__file__).resolve().parent.parent
AGENT = REPO_ROOT / "delfin" / "agent"

# No --tmpfs /tmp: with it, a target under /tmp has a WRITABLE tmpfs parent
# and bwrap creates the leaf itself, so the probe passes while measuring
# nothing. The first version of this probe did exactly that.
BASE = ["bwrap", "--ro-bind", "/", "/", "--dev", "/dev", "--proc", "/proc"]

needs_bwrap = pytest.mark.skipif(
    not shutil.which("bwrap"),
    reason="bwrap is not installed; the argv-level checks below still run")


def _callers_of(name: str) -> set[str]:
    """Files in delfin/agent whose AST contains a call to ``name``."""
    out = set()
    for path in sorted(AGENT.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                fn = node.func
                if (getattr(fn, "id", None) == name
                        or getattr(fn, "attr", None) == name):
                    out.add(path.name)
    return out


def test_every_bwrap_builder_applies_the_guarantee():
    """The reachability check, which is what was actually wrong: the step
    existed, was tested, and one of the two builders never called it."""
    callers = _callers_of("ensure_argv_bind_targets")
    assert "sandbox.py" in callers, (
        "sandbox._bwrap_argv's argv goes to bwrap unprepared")
    assert "contained_run.py" in callers, "BW1's caller is still there"
    assert len(callers) >= 2, callers


def test_the_two_single_layer_home_secrets_are_denied_at_the_tool_level():
    """A home secret whose only protection is the cage mask is protected
    only when the mask is emitted, and the mask is emitted only for a path
    that already exists when the argv is built. These two had no second
    layer; the other home secrets are denied by name."""
    from delfin.agent.api_client import _DEFAULT_PATH_DENY_GLOBS as globs

    for rel in (".claude/settings.json", ".claude",
                ".local/share/delfin-messenger/db"):
        assert any(fnmatch.fnmatch(rel, g) for g in globs), rel
    # The deny list must not grow into ordinary files. .config/git is named
    # in sandbox.py as deliberately NOT masked (hiding it made git fatal).
    for rel in (".config/git/ignore", "src/main.py", "README.md"):
        assert not any(fnmatch.fnmatch(rel, g) for g in globs), rel


def test_an_existing_home_secret_always_gets_a_mask(tmp_path):
    """What the mask list does guarantee today, pinned so a change to the
    existence test is visible: every secret path that EXISTS is masked."""
    from delfin.agent import sandbox

    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".netrc").write_text("stand-in\n")
    args = sandbox._bwrap_home_mask_args(home=home)
    assert str(home / ".ssh") in args, args
    assert str(home / ".netrc") in args, args
    # A path that does NOT exist gets no mask at all. That is the half of
    # BW2 this commit does not close: closing it needs each entry's SHAPE
    # declared (tmpfs for a directory, /dev/null for a file), and the shape
    # of .delfin/approval_codes and .delfin/trust_store cannot be
    # established without reading the operator's real store. Asserted, not
    # assumed, so the brief that asks for it cannot go stale.
    assert str(home / ".aws") not in args


@needs_bwrap
def test_a_missing_mask_target_under_the_readonly_root_kills_the_command(tmp_path):
    target = tmp_path / "secret_dir"
    argv = BASE + ["--tmpfs", str(target), "true"]
    first = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    if first.returncode == 0:          # a host where bwrap can create it
        pytest.skip("this host's bwrap created the leaf; nothing to measure")
    assert "mkdir" in first.stderr or "Read-only" in first.stderr, first.stderr

    cage_binds.ensure_argv_bind_targets(argv)
    assert target.exists(), "the step creates what the argv names"
    second = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    assert second.returncode == 0, second.stderr


@needs_bwrap
def test_the_mask_still_masks_after_the_target_is_created(tmp_path):
    """Protection both ways: creating the target must not turn a mask into
    a readable directory. Without this, the fix above could be 'make it
    run' at the price of the thing it was hiding."""
    target = tmp_path / "secret_dir"
    target.mkdir()
    (target / "key.txt").write_text("stand-in secret\n")
    argv = BASE + ["--tmpfs", str(target), "cat", str(target / "key.txt")]
    cage_binds.ensure_argv_bind_targets(argv)
    run = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    assert run.returncode != 0
    assert "stand-in secret" not in run.stdout
