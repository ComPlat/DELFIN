"""Package SWP phase 1: role policies as DATA (delfin/agent/swarm/policy.py).

The role deny patterns live in policy.py as data, generalised so no home
path, user name or host is hard-coded: the "other people's homes" patterns
are derived from Path.home() at runtime. A pattern list refuses when at
least one pattern matches a command (must-refuse); a harmless command must
match none (must-allow).

Red control: policy.py does not exist on main, so this module fails to
import there (ModuleNotFoundError) -- that is the red on the previous
commit.
"""

import re
from pathlib import Path

import pytest

import delfin.agent.swarm.policy as policy


class FakeHome:
    """A real-ish home for the runtime derivation (no files touched)."""

    def __init__(self, path):
        self._p = Path(path)

    @property
    def name(self):
        return self._p.name

    def __str__(self):
        return str(self._p)


def _refuses(patterns, cmd):
    """True when at least one deny pattern matches the command."""
    return any(re.search(p, cmd) for p in patterns)


def _allows(patterns, cmd):
    return not _refuses(patterns, cmd)


HOME = Path("/home/alice")


@pytest.fixture(scope="module")
def home():
    return HOME


# ---------------------------------------------------------------------------
# Every role refuses the member-wide COMMON set.
# ---------------------------------------------------------------------------

COMMON_REFUSE = [
    "rm -rf results/",
    "unlink want.txt",
    "shred report.pdf",
    "git worktree add /tmp/x main",
    "git fetch origin",
    "git pull origin main",
    "git apply change.patch",
    "git am series.patch",
    "patch -p1 < change.patch",
    "python -m pytest tests/test_sw1_actions.py",
    "pytest tests/test_x.py && ls",
    "gate tests/test_a.py || echo failed",
    "squeue -u me",
    "sbatch run.sh",
    "scancel 12345",
    "scontrol show job 7",
    "printenv PATH",
    "cat credentials.json",
    "cat ~/.ssh/id_ed25519",
    "cat project/.env",
    "cat /home/bob/report.txt",
    "cat /pfs/data3/home/carol/data.txt",
    "cat .delfin/settings.json",
]


@pytest.mark.parametrize("cmd", COMMON_REFUSE)
def test_common_deny_patterns_refuse(cmd, home):
    pats = policy.deny_patterns_for("builder", home=home)
    assert _refuses(pats, cmd), f"COMMON must refuse: {cmd}"


COMMON_ALLOW = [
    "git log --oneline -5",
    "git show HEAD",
    "git diff --stat",
    "git status",
    "ls tests/",
    "cat delfin/agent/swarm/launch.py",
    "gate -q tests/test_swp_policy.py",
    "/abs/path/to/gate-tools/gate tests/test_swp_policy.py",
    "cat /home/alice/project/README.md",   # own home is allowed
    "cat /tools/readme.txt",
]


@pytest.mark.parametrize("cmd", COMMON_ALLOW)
def test_common_deny_patterns_allow(cmd, home):
    pats = policy.deny_patterns_for("builder", home=home)
    assert _allows(pats, cmd), f"COMMON must allow: {cmd}"


@pytest.fixture(scope="module")
def builder_pats(home):
    return policy.deny_patterns_for("builder", home=home)


@pytest.fixture(scope="module")
def reviewer_pats(home):
    return policy.deny_patterns_for("reviewer", home=home)


@pytest.fixture(scope="module")
def qs_pats(home):
    return policy.deny_patterns_for("qs", home=home)


# ---------------------------------------------------------------------------
# NO_PUSH for every role -- the push is done by the user/bot.
# ---------------------------------------------------------------------------

PUSH_REFUSE = [
    "git push origin main",
    "git push -u origin agent/s14-swpb",
    "git push --force origin head",
]


@pytest.mark.parametrize("role", ["builder", "reviewer", "qs"])
@pytest.mark.parametrize("cmd", PUSH_REFUSE)
def test_no_push_for_every_role(role, cmd, home):
    pats = policy.deny_patterns_for(role, home=home)
    assert _refuses(pats, cmd), f"role {role} must refuse push: {cmd}"


# ---------------------------------------------------------------------------
# Reviewer: no commit of product code (only its adversarial tests).
# ---------------------------------------------------------------------------

REVIEWER_REFUSE = [
    "git add delfin/agent/swarm/launch.py",
    "git commit -am shipped the feature",
    "git add -A",
]


@pytest.mark.parametrize("cmd", REVIEWER_REFUSE)
def test_reviewer_refuses_product_commits(cmd):
    assert _refuses(policy.deny_patterns_for("reviewer"),
                    cmd), f"reviewer must refuse: {cmd}"


def test_reviewer_may_stage_its_own_adversarial_tests():
    pats = policy.deny_patterns_for("reviewer")
    # The reviewer's own test glob is NOT product code.
    assert _allows(pats, "git add tests/test_swp_adversarial.py"), (
        "reviewer may stage its own tests")


def test_builder_may_stage_its_own_files():
    pats = policy.deny_patterns_for("builder")
    assert _allows(pats, "git add delfin/agent/swarm/policy.py")
    assert _allows(pats, "git add tests/test_swp_policy.py")


# ---------------------------------------------------------------------------
# Generalisation: home path and user name are DERIVED, never hard-coded.
# ---------------------------------------------------------------------------

def test_no_literal_home_or_user_in_the_pattern_source():
    src = Path(policy.__file__).read_text(encoding="utf-8")
    # The class of values the task forbids: a concrete user or host.
    for banned in ("ka_ibcs", "/home/ka/", "/pfs/data6/home/ka",
                   "@localhost", "hostname"):
        assert banned not in src, f"policy.py must not hard-code {banned!r}"


def test_another_home_is_derived_at_runtime_from_the_member_home():
    # The patterns for a different member home must refuse THAT member's
    # other-user homes but allow the member's own.
    other_home = FakeHome("/home/bob")
    pats = policy.deny_patterns_for("builder", home=other_home)
    assert _refuses(pats, "cat /home/carol/file.txt")
    assert _allows(pats, "cat /home/bob/project/README.md")


def test_role_policy_is_exposed_as_data():
    # "Role policies as DATA": the caller can read the policy dict.
    data = policy.ROLE_DENY
    assert set(data) == {"builder", "reviewer", "qs"}
    for role, extra in data.items():
        assert isinstance(extra, (tuple, list))


def test_valid_roles_match_the_spec_vocabulary():
    # The role vocabulary agrees with spec.VALID_ROLES.
    from delfin.agent.swarm import spec as _spec
    assert set(policy.ROLE_DENY) == set(_spec.VALID_ROLES)
