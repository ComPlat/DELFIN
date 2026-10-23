"""A dashboard a script starts opens no browser; one a person starts still may.

2026-09-16: every delfin-voila launched from an agent's shell -- browser
drives, the dashboard tests of a gate -- opened a tab on the developer's
desktop, because the shell had inherited the VS Code terminal variables
the launcher reads. A person typing the command in that terminal is the
case the auto-open was built for; a child process of that terminal is not.
"""
import argparse

from delfin import cli_voila


def _args(open_browser=None, no_browser=None):
    return argparse.Namespace(open_browser=open_browser, no_browser=no_browser)


VSCODE = {"TERM_PROGRAM": "vscode", "VSCODE_IPC_HOOK_CLI": "/run/user/1/vscode-ipc.sock"}


def test_a_person_in_a_vscode_terminal_gets_the_tab():
    assert cli_voila._decide_open_browser(_args(), VSCODE, interactive=True) is True


def test_a_script_in_the_same_environment_gets_none():
    assert cli_voila._decide_open_browser(_args(), VSCODE, interactive=False) is False


def test_the_environment_can_switch_it_off_for_everyone():
    env = dict(VSCODE, DELFIN_NO_BROWSER="1")
    assert cli_voila._decide_open_browser(_args(), env, interactive=True) is False


def test_the_explicit_switches_win_over_everything():
    env = dict(VSCODE, DELFIN_NO_BROWSER="1")
    assert cli_voila._decide_open_browser(_args(open_browser=True), env, interactive=False) is True
    assert cli_voila._decide_open_browser(_args(no_browser=True), VSCODE, interactive=True) is False


def test_outside_vscode_nothing_opens_by_itself():
    assert cli_voila._decide_open_browser(_args(), {"TERM_PROGRAM": "tmux"}, interactive=True) is False


def test_main_asks_the_decision_not_the_terminal_alone():
    import inspect
    src = inspect.getsource(cli_voila.main)
    assert "open_browser = _decide_open_browser(args)" in src
    assert "open_browser = _is_vscode_session()" not in src
