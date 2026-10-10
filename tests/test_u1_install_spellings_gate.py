"""Install-outside-the-session-venv refusal, measured at the executor gate.

Input: a bash command string, driven through
``_doc_executor._run_permission_gate`` in ``bypassPermissions`` mode (the
mode in which nothing else refuses an install). Output: the gate's error.

Each case below is a spelling of an install into a location other than the
session venv that pip, uv, python or micromamba accept as written (pip
resolves long options by unambiguous prefix and short options in clusters;
python accepts interpreter flags before ``-m`` and ``-mpip``). The first cut
of the rule matched fixed substrings and let all of them run. The second
group are installs into the session venv itself, which must stay allowed.
The dashboard approval path (``sandbox.is_allowed``) applies the same rule.
"""

from __future__ import annotations

import pytest

from delfin import installer
from delfin.agent import sandbox
from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _gate(cmd: str, tmp_path) -> str:
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = KitToolPermissions(workspace=str(ws), mode="bypassPermissions")
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms) or ""


OUTSIDE = [
    "pip install --targ /tmp/x somepkg",          # optparse abbreviation
    "pip install --prefi /tmp/p somepkg",
    "pip install -qt /tmp/x somepkg",             # short-option cluster
    "pip install --us''er somepkg",               # shell quote removal
    "python3 -I -m pip install somepkg",          # interpreter flag first
    "python3 -mpip install somepkg",
    "/usr/bin/pip3 install somepkg",              # another interpreter's pip
    "pip --python /usr/bin/python3 install somepkg",
    "uv pip install --python /usr/bin/python3 somepkg",
    "uv pip install -p /usr/bin/python3 somepkg",
    "micromamba install -n base somepkg",
    "conda install --name=base somepkg",
    "conda install -nbase somepkg",
    "pipx install somepkg",
    "uv tool install somepkg",
    "PIP_ROOT=/tmp/r pip install somepkg",
    "pip config set global.user true",
    "true;pip install --user somepkg",
    "bash -c 'pip install --user somepkg'",
    "env -i PATH=/usr/bin pip install --user somepkg",
    'echo "$(pip install --user somepkg)"',
    "echo `pip install --user somepkg`",
]


@pytest.mark.parametrize("cmd", OUTSIDE)
def test_every_spelling_of_an_outside_install_is_refused(cmd, tmp_path):
    err = _gate(cmd, tmp_path)
    assert "unsafe install" in err and "session venv" in err.lower(), (cmd, err)


@pytest.mark.parametrize("cmd", OUTSIDE)
def test_the_dashboard_approval_path_refuses_the_same(cmd):
    assert not sandbox.is_allowed(cmd).allowed, cmd


def _inside() -> list[str]:
    py = installer.session_python()
    return [
        f"{py} -m pip install pytest",
        f"{py} -m pip install --root-user-action=ignore pytest",
        f"{py} -m pip install --upgrade --pre pytest",
        f"{py} -m pip install -U -r requirements-test.txt",
        f"{py} -m pip install --python {py} pytest",
        f"{py} -m pip uninstall pytest pytest-timeout --yes",
        "pip install -r requirements.txt",
        "pip install -e .",
        # Text that names an outside install without running one.
        'git commit -m "Refuse pip install --user in every mode"',
        "grep -n 'pip install --user' delfin/agent/sandbox.py",
    ]


def test_an_install_into_the_session_venv_is_not_refused(tmp_path):
    for cmd in _inside():
        err = _gate(cmd, tmp_path)
        assert "unsafe install" not in err, (cmd, err)


def test_the_catalog_command_passes_the_gate(tmp_path):
    tool = installer.find("pytest")
    cmd = installer.python_tools_install_command(tool)
    assert cmd
    assert "unsafe install" not in _gate(cmd, tmp_path), cmd


def test_without_a_known_session_interpreter_python_m_pip_is_refused():
    from delfin.agent.api_client import _refuse_unsafe_install
    assert _refuse_unsafe_install("/venv/bin/python -m pip install x", "")


@pytest.mark.parametrize("cmd", [
    ".venv-demo/bin/pip install -r requirements.txt",
    ".venv/bin/python -m pip install somepkg",
])
def test_a_venv_inside_the_workspace_is_not_refused(cmd, tmp_path):
    assert "unsafe install" not in _gate(cmd, tmp_path), cmd


def test_a_workspace_venv_is_still_held_to_the_options(tmp_path):
    err = _gate(".venv/bin/pip install --user somepkg", tmp_path)
    assert "unsafe install" in err, err
