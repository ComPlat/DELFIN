"""Adversarial review tests for SWP target 1 (delfin/agent/swarm/policy.py).

Reviewer nacht-s15. These pin the cases a plain deny list misses, so a
deny pattern that only matches its own red control's exact spelling fails
here. They went RED on eeb76a05a (reviewer's adversarial probe found the
home-traversal gap) and are GREEN on the builder's hardening b9d508f88 --
the traversal fix never over-blocks the allow cases.

The cases:
- abbreviated / attached values in dangerous commands,
- env prefixes (X=y cmd),
- numbered and >&word redirections, unicode look-alikes,
- git apply --check MUST be allowed (patch application except --check),
- direct pytest with a harmless name must be refused (tests only via the
  project's test entry),
- home derivation must handle ~, $HOME and ${HOME} style references to
  OTHER users AND refuse a `../` walk OUT of the member's own home, while
  allowing the member's own paths.
"""

import re
from pathlib import Path

import pytest

import delfin.agent.swarm.policy as policy


def _refuses(patterns, cmd):
    return any(re.search(p, cmd) for p in patterns)


HOME = Path("/home/alice")


def pats(role="builder"):
    return policy.deny_patterns_for(role, home=HOME)


# ---------------------------------------------------------------------------
# These MUST be refused. Each one is a shape the builder's red control
# misses; a deny pattern that only matches the exact spelling in that
# control fails here.
# ---------------------------------------------------------------------------

MUST_REFUSE = [
    # attached value: no space between -o and the file
    "python -m pytest -o cache_dir=/tmp/x tests/test_x.py",
    # abbreviated long options, alone and combined
    "git apply --exclude=utils x.patch",
    "git apply --exclude=utils --include=only.py x.patch",
    # env prefix hiding the real command
    "PYTHONPATH=. python -m pytest tests/test_sw1_actions.py",
    "BUILDER_MODE=1 rm -rf results/",
    "X=y git push origin main",
    # numbered and >&word redirections (patch application via redirection)
    "patch -p1 < change.patch 2> err.log",
    "patch 2>err.log -p1 < x.diff",
    "patch -p1 < x.diff > applied.log",
    "patch -p1 < x.diff >&2",
    "patch -p1 < x.diff 1>&2",
    # unicode look-alikes: FULLWIDTH SEMICOLON as a command separator,
    # non-breaking space inside the command word
    "rm -rf results/；cat ~/.ssh/id_ed25519",
    "rm\u00a0-rf\u00a0results/",
    # direct pytest (not the project's test entry) with a harmless name
    "pytest -q tests/test_x.py",
    "python -m pytest tests/test_x.py",
    # patch applied via git am with attached short flag
    "git am -3series.patch",
    # home references in ~ and env style for OTHER users
    "cat ~bob/report.txt",
    # a `../` walk OUT of the member's own home into another user's area
    "cat $HOME/../bob/report.txt",
    "cat ${HOME}/../bob/report.txt",
    "cat ~/../bob/report.txt",
    "cat ~carol/.ssh/id_ed25519",
]


@pytest.mark.parametrize("cmd", MUST_REFUSE)
def test_adversarial_must_refuse(cmd):
    assert _refuses(pats(), cmd), f"COMMON must refuse: {cmd!r}"


# ---------------------------------------------------------------------------
# These MUST still be allowed (the deny list must not over-block).
# ---------------------------------------------------------------------------

MUST_ALLOW = [
    # the sanctioned patch CHECK (task: "no patch application except --check")
    "git apply --check change.patch",
    "git apply --check --3way change.patch",
    # pytest named in the allowed test entry (the gate wrapper)
    "gate -q tests/test_swp_policy.py",
    "/abs/path/gate-tools/gate tests/test_a.py tests/test_b.py::test_x",
    # ordinary reads inside the member's own home, in env style too
    "cat /home/alice/project/README.md",
    "cat $HOME/project/README.md",
    "cat ${HOME}/project/README.md",
    "cat ~/project/README.md",
    # a word that merely CONTAINS a dangerous substring
    "grep rm_from_plan delfin/agent/swarm/launch.py",
    "cat notes/patch_review.md",
]


@pytest.mark.parametrize("cmd", MUST_ALLOW)
def test_adversarial_must_allow(cmd):
    assert not _refuses(pats(), cmd), f"COMMON must allow: {cmd!r}"


# ---------------------------------------------------------------------------
# No hard-coded identity in the SOURCE beyond the banned literals the
# builder already tests (ka_ibcs, /home/ka, hostnames): also check the
# compiled patterns never contain an absolute home literal.
# ---------------------------------------------------------------------------

def test_no_role_pattern_contains_a_literal_home():
    # The STATIC data (COMMON, NO_PUSH, role extras) must carry no home
    # literal. The home-derived patterns are generated at call time from
    # the caller's home by design, so they are not checked here.
    static = (policy.COMMON_DENY_PATTERNS + policy.NO_PUSH_PATTERNS
              + policy.REVIEWER_NO_CODE_PATTERNS
              + policy.QS_TEST_SCOPE_PATTERNS)
    for p in static:
        assert str(HOME) not in p, (
            f"static pattern hard-codes the member home: {p}")
        assert re.search(r"/home/[\w.-]+/", p) is None, (
            f"static pattern hard-codes a home path: {p}")
