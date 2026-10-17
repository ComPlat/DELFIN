"""Package A2 reviewer (nacht-s25): adversarial part 3 -- write-target
containment, quoting, symlink and .. escapes. Every test pins a REFUSAL.
"""

import os

import delfin.agent.workspace_selfservice as ws


def _ws(tmp_path):
    return str(tmp_path)


def _venv(tmp_path):
    venv = tmp_path / ".venv"
    (venv / "bin").mkdir(parents=True)
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
    return venv


# --- F5: containment escapes on the approved shapes --------------------------

def test_redirect_inside_is_none_not_auto(tmp_path):
    # Even an INSIDE redirect on a pinned auto shape is unrated (the pinned
    # forms have no redirects); the gate decides.
    _venv(tmp_path)
    cmd = "git switch -c feature > switch.log"
    assert ws.classify(cmd, _ws(tmp_path)) is None, "F5: redirect on auto shape"


def test_pip_venv_outside_via_dotdot(tmp_path):
    # venv prefix outside the workspace: the workspace is a SUBFOLDER of
    # tmp_path and the venv sits one level above it.
    outside = tmp_path / "venv"
    (outside / "bin").mkdir(parents=True)
    (outside / "reqs.txt").write_text("numpy==1.26.4\n")
    ws_root = tmp_path / "ws"
    ws_root.mkdir()
    (ws_root / "requirements.txt").write_text("numpy==1.26.4\n")
    venv_pip = str(outside / "bin" / "pip")
    cmd = f"{venv_pip} install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(ws_root)) is None, "F5: venv outside"


def test_pip_req_dotdot_escape(tmp_path):
    _venv(tmp_path)
    cmd = ".venv/bin/pip install -r ../reqs.txt --require-hashes"
    assert ws.classify(cmd, _ws(tmp_path)) is None, "F5: ../ requirements"


def test_pip_req_symlink_escape(tmp_path):
    # requirements file inside the workspace but a symlink to outside.
    ws_root = tmp_path / "ws"
    ws_root.mkdir()
    (ws_root / ".venv").mkdir()
    (ws_root / ".venv" / "bin").mkdir()
    (ws_root / "reqs_inside.txt").write_text("numpy==1.26.4\n")
    outside_dir = tmp_path / "outside"
    outside_dir.mkdir()
    (outside_dir / "evil.txt").write_text("numpy==1.26.4\n")
    link = ws_root / "requirements.txt"
    os.symlink(str(outside_dir / "evil.txt"), str(link))
    cmd = ".venv/bin/pip install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(ws_root)) is None, "F5: symlinked requirements"


def test_pip_multiple_r_flags_last_wins_is_none(tmp_path):
    # Two -r flags: pip accepts the LAST one. The matcher records req on each
    # iteration so the last one is checked -- but BOTH must be inside for the
    # pinned form; this pins the last-wins reading stays refused when either
    # is outside. The workspace is a SUBFOLDER of tmp_path so the evil
    # requirements file genuinely sits outside it (it did not, in the first
    # version: outside_dir was a sibling of the workspace root, i.e. inside).
    ws_root = tmp_path / "ws"
    ws_root.mkdir()
    venv = ws_root / ".venv"
    (venv / "bin").mkdir(parents=True)
    (ws_root / "requirements.txt").write_text("numpy==1.26.4\n")
    outside_dir = tmp_path / "outside"
    outside_dir.mkdir()
    (outside_dir / "evil.txt").write_text("numpy==1.26.4\n")
    cmd = (".venv/bin/pip install -r requirements.txt -r "
           f"{outside_dir / 'evil.txt'} --require-hashes")
    assert ws.classify(cmd, str(ws_root)) is None, "F5: double -r escape"


def test_pip_aborted_dash_r_is_none(tmp_path):
    _venv(tmp_path)
    cmd = ".venv/bin/pip install -r --require-hashes"
    assert ws.classify(cmd, _ws(tmp_path)) is None, "F5: -r swallows the flag"


def test_pip_quoted_req_with_space_is_none(tmp_path):
    # A quoted path with a space: containment must resolve the WHOLE word,
    # not split on whitespace.
    _venv(tmp_path)
    (tmp_path / "my reqs.txt").write_text("numpy==1.26.4\n")
    cmd = (".venv/bin/pip install -r 'my reqs.txt' --require-hashes "
           "--only-binary=:all:")
    assert ws.classify(cmd, _ws(tmp_path)) == "auto", \
        "F5: quoted path with space inside workspace should stay auto"


def test_pip_tilde_escape(tmp_path):
    # ~ expands to the user's home -- outside the workspace.
    _venv(tmp_path)
    cmd = ".venv/bin/pip install -r ~/evil/reqs.txt --require-hashes"
    assert ws.classify(cmd, _ws(tmp_path)) is None, "F5: ~ expansion"
