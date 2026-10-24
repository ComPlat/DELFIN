# FR1 reviewer adversarial tests (nacht-s17): the four routine shapes in
# delfin/agent/routine_actions.py. Attacks beyond the builder's
# tests/test_fr1_routine_shapes.py: flag look-alikes on merge --abort,
# --orphan, env prefixes, unicode look-alikes, quoted tokens,
# bash_background parity, timeout_s of the wrong type, userinfo URLs.
# Both directions: must-NOT-be-routine and must-still-be-routine.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from delfin.agent.routine_actions import routine_verdict  # noqa: E402


def _bash(cmd):
    return routine_verdict({"tool": "bash", "command": cmd}, {})


def _fetch(**kw):
    action = {"tool": "web_fetch"}
    action.update(kw)
    return routine_verdict(action, {})


def _none_on(cmd):
    assert _bash(cmd) is None, f"unexpected routine verdict: {cmd}"


def test_merge_abort_with_any_extra_argument_is_not_routine():
    _none_on("git merge --abort --quiet")
    _none_on("git merge --no-edit --abort")
    _none_on("git merge --continue")
    _none_on("git merge")
    _none_on("git merge --abort extra-arg")


def test_checkout_b_flag_lookalikes_still_judged():
    _none_on("git checkout --orphan newb")
    _none_on("git checkout -b -f newb")
    _none_on("git checkout -b main -- extrafile")
    _none_on("git switch -C newb")
    _none_on("git switch -c newb --track origin/main")  # corrected: F2


def test_env_prefix_and_unicode_commands_are_not_routine():
    _none_on("X=y git checkout -b newb")
    _none_on("git -C . checkout -b newb")
    _none_on("\u0261it checkout -b newb")  # unicode g look-alike


def test_quoted_routine_tokens_match_the_same_shape():
    # Quotes are not parsed by the shape matcher. A quoted --abort does
    # NOT match the exact shape (token '"--abort"' != '--abort') -- the
    # safe direction, it stays judged. A quoted branch name matches and
    # the allow is harmless (creates a branch literally named 'newb').
    assert _bash('git merge "--abort"') is None
    assert _bash("git checkout -b 'newb'") is not None


def test_bash_background_parity_with_bash():
    # Same tool surface, same judgement: routine on the plain shape,
    # judged on the look-alike.
    assert _bash("git merge --abort") is not None
    v = routine_verdict({"tool": "bash_background",
                         "command": "git merge --abort"}, {})
    assert v is not None
    assert routine_verdict({"tool": "bash_background",
                            "command": "git push"}, {}) is None


def test_checkout_b_onto_existing_branch_verdict_is_harmless():
    # `git checkout -b existing` makes GIT fail (branch exists), not the
    # classifier: the allow verdict may stand -- it writes nothing.
    assert _bash("git checkout -b main") is not None
    assert _bash("git switch -c main") is not None


def test_git_mv_positional_and_unicode_paths_match_shape():
    # A plain rename with no flag or wildcard is the shape; the write
    # gate holds the paths. Empty-ish and flag-like names still judged.
    assert _bash("git mv a b") is not None
    _none_on("git mv -")
    _none_on("git mv -- a b")
    _none_on("git mv a b c")


def test_web_fetch_timeout_type_and_boundary():
    assert _fetch(url="https://example.org/x", timeout_s=True) is None
    assert _fetch(url="https://example.org/x", timeout_s="60") is None
    assert _fetch(url="https://example.org/x", timeout_s=0) is None
    assert _fetch(url="https://example.org/x", timeout_s=61) is None
    assert _fetch(url="https://example.org/x", timeout_s=60) is not None
    assert _fetch(url="https://example.org/x", timeout_s=0.5) is not None


def test_web_fetch_userinfo_and_fragment_are_still_get_only():
    # Credentials in a URL are a secret-handling question the model must
    # see. As committed (F3 open): userinfo IS currently routine -- pinned
    # until the builder rules on F3.
    assert _fetch(url="https://user@host/x") is not None
    assert _fetch(url="http://host/x#frag") is not None


def test_web_fetch_missing_url_and_empty_are_not_routine():
    assert _fetch() is None
    assert _fetch(url="   ") is None
