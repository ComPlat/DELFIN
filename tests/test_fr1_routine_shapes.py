"""FR1 phase 1: routine shapes run without a dialog.

Red control for delfin/agent/routine_actions.routine_verdict, wired into
action_classifier.screen(). Each shape is an allow ONLY when the
code-checked conditions hold; the look-alikes stay judged as today.

The operator's additions (wave 20 stage 1d, FR1 package text): under
bypass mode curl, wget, pip install, ssh to other nodes and
git commit -am run in wave trees, so NONE of those shapes may become
routine in the classifier. Each is red-tested here as an ask/refuse,
never an allow.
"""

from delfin.agent import action_classifier as ac


def _verdict(action):
    return ac.screen(action, {"task": "FR1: routine actions"})


# --- (a) git merge --abort -------------------------------------------------

def test_git_merge_abort_is_routine():
    v = _verdict({"tool": "bash", "command": "git merge --abort"})
    assert v is not None and v.decision == "allow"


def test_git_merge_abort_lookalike_still_judged():
    v = _verdict({"tool": "bash",
                  "command": "git merge --abort; git push"})
    assert v is None


# --- (b) new local branch ---------------------------------------------------

def test_checkout_b_is_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout -b agent/s10-topic origin/main"})
    assert v is not None and v.decision == "allow"


def test_switch_c_is_routine():
    v = _verdict({"tool": "bash", "command": "git switch -c new-thing"})
    assert v is not None and v.decision == "allow"


def test_checkout_b_force_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git checkout -B existing"})
    assert v is None


def test_checkout_b_f_force_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git checkout -b -f thing"})
    assert v is None


def test_checkout_onto_existing_branch_with_paths_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout main -- delfin/agent/api_client.py"})
    assert v is None


def test_switch_c_force_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git switch -C thing"})
    assert v is None


# --- (c) git mv inside scope ------------------------------------------------

def test_git_mv_two_paths_is_routine():
    v = _verdict({"tool": "bash",
                  "command": "git mv delfin/agent/routine_actions.py"
                             " .gate/moved.py"})
    assert v is not None and v.decision == "allow"


def test_git_mv_force_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git mv -f a b"})
    assert v is None


def test_git_mv_k_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git mv -k a b"})
    assert v is None


def test_git_mv_directory_wildcard_not_routine():
    v = _verdict({"tool": "bash", "command": "git mv srcdir/*.py dst/"})
    assert v is None


# --- (d) read-only web_fetch -------------------------------------------------

def test_web_fetch_get_is_routine():
    v = _verdict({"tool": "web_fetch",
                  "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC123/"})
    assert v is not None and v.decision == "allow"


def test_web_fetch_file_url_not_routine():
    v = _verdict({"tool": "web_fetch", "url": "file:///etc/passwd"})
    assert v is None


def test_web_fetch_ftp_not_routine():
    v = _verdict({"tool": "web_fetch", "url": "ftp://example.org/x"})
    assert v is None


def test_web_fetch_post_body_not_routine():
    v = _verdict({"tool": "web_fetch",
                  "url": "https://example.org/api",
                  "payload": {"write": "true"}})
    assert v is None


def test_web_fetch_extra_args_not_routine():
    v = _verdict({"tool": "web_fetch",
                  "url": "https://example.org/", "method": "POST"})
    assert v is None


# --- operator's never-routine shapes -----------------------------------------

def test_curl_is_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "curl https://example.org/upload -d x=1"})
    assert v is None


def test_wget_is_not_routine():
    v = _verdict({"tool": "bash", "command": "wget http://example.org/f"})
    assert v is None


def test_pip_install_is_not_routine():
    v = _verdict({"tool": "bash", "command": "pip install ruff"})
    assert v is None


def test_ssh_is_not_routine():
    v = _verdict({"tool": "bash", "command": "ssh login.cluster2"})
    assert v is None


def test_git_commit_am_is_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git commit -am 'quick fix in wave tree'"})
    assert v is None


def test_git_push_is_not_routine():
    v = _verdict({"tool": "bash", "command": "git push origin HEAD"})
    assert v is None


# --- never-rules and injection guard keep precedence --------------------------

def test_user_never_rule_beats_routine():
    rules = ("never allow git mv",)
    v = ac.screen({"tool": "bash", "command": "git mv a b"},
                  {"task": "FR1"}, rules)
    assert v is not None and v.decision == "refuse"


def test_routine_only_for_bash_and_native_web_fetch():
    # An MCP web variant is a different surface, judged as today.
    v = _verdict({"tool": "mcp__kit__web_fetch", "url": "https://example.org/"})
    assert v is None
