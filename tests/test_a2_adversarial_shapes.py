"""Package A2 reviewer (nacht-s25): adversarial part 2 -- branch shape and
version-shape escapes. Every test pins a REFUSAL (classify -> None).
"""

import delfin.agent.workspace_selfservice as ws


def _ws(tmp_path):
    return str(tmp_path)


# --- F3: branch-name shape escapes -------------------------------------------
# _SAFE_BRANCH_RE (workspace_selfservice.py:222) accepts [A-Za-z0-9][A-Za-z0-9._-]{0,200}.

def test_branch_unicode_lookalike_is_none(tmp_path):
    # Cyrillic homoglyph of "main" must not pass -- and must not be read as a
    # foreign name either; it is simply unrated -> None.
    assert ws.classify("git switch -c fаature", _ws(tmp_path)) is None, \
        "F3: unicode lookalike branch"


def test_branch_dotted_escape_is_none(tmp_path):
    assert ws.classify("git switch -c a..b", _ws(tmp_path)) is None, \
        "F3: a..b branch (dotted)"


def test_branch_hidden_dot_name_is_none(tmp_path):
    assert ws.classify("git switch -c .hidden", _ws(tmp_path)) is None, \
        "F3: .hidden branch"


def test_branch_master_hyphen_suffix_is_none(tmp_path):
    # "master-2" is NOT main/master itself but a plain word the shape matches;
    # an own-feature name is acceptable, this pins the RE does not over-match
    # whitespace/control chars, not that master-2 is refused.
    assert ws.classify("git switch -c ok name", _ws(tmp_path)) is None, \
        "F3: two-word branch args"


def test_switch_dash_c_attached_value_is_none(tmp_path):
    # `git switch -cfeature` is a DIFFERENT argv than `git switch -c feature`:
    # git would parse -c as its own flag and "feature" as the start point.
    # The pinned form is 4 words ["git","switch","-c",<branch>]; the attached
    # form is 3 words with words[2]=="-cfeature" -- it must NOT auto.
    assert ws.classify("git switch -cfeature", _ws(tmp_path)) is None, \
        "F3: -c attached value"


def test_switch_start_point_form_is_none(tmp_path):
    # `git switch <name> <start>` (5 words) is outside the pinned form.
    assert ws.classify("git switch feature origin/main", _ws(tmp_path)) is None, \
        "F3: switch with start point"


def test_checkout_without_b_is_none(tmp_path):
    # `git checkout <name>` switches branches too but is NOT in the pinned
    # forms (only checkout -b) -> None.
    assert ws.classify("git checkout feature", _ws(tmp_path)) is None, \
        "F3: checkout without -b"


def test_checkout_b_with_start_point_is_none(tmp_path):
    assert ws.classify("git checkout -b feature origin/main", _ws(tmp_path)) is None, \
        "F3: checkout -b with start point"


# --- F4: version-shape escapes ------------------------------------------------

def test_version_dash_v_is_none(tmp_path):
    # -v/-V are not universally version flags (tar -v = verbose); only
    # --version and -V and "version" are pinned; -v must be None.
    assert ws.classify("gh -v", _ws(tmp_path)) is None, "F4: -v not pinned"


def test_version_lowercase_v_uppercase_mismatch_is_none(tmp_path):
    # `bash -V` is NOT a version flag (bash uses --version); -V must only be
    # trusted per known tool, and the builder pinned -V globally. This test
    # pins the decision: if -V is rated for all listed tools, an unknown-to-
    # -V tool like bash must still not be auto -- the builder's list includes
    # bash, so this asserts None as the SAFE verdict.
    assert ws.classify("bash -V", _ws(tmp_path)) is None, "F4: bash -V"


def test_version_aborted_flag_form_is_none(tmp_path):
    assert ws.classify("gh version", _ws(tmp_path)) is None, "F4: 'gh version'"


def test_version_path_to_tool_is_none(tmp_path):
    # A venv-prefixed interpreter --version would still be a pinned tool by
    # words[0] matching; the pinned list matches BARE names only, so a path
    # head must not auto.
    assert ws.classify(".venv/bin/python --version", _ws(tmp_path)) is None, \
        "F4: path-prefixed --version"
