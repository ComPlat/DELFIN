"""A fresh install has no language model: the agent says so and how to connect one.

Measured on a fresh install (2026-11-13): with nothing configured the shipped
default (CLI backend) ran whatever program of that name was on PATH; an old
one answered "unknown option '--include-hook-events'" and the turn failed. On
a clean PATH with no key the dashboard hid the DELFIN Agent tab altogether.

Here: an empty HOME, a PATH whose CLI program is either absent or one that
does not know the flag DELFIN passes, no keys. Nothing reaches the network.
"""

import ast
import os
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]


def _program(tmp_path, *, knows_flag):
    bindir = tmp_path / "bin"
    bindir.mkdir(parents=True, exist_ok=True)
    text = "  --include-hook-events  stream hook events\n" if knows_flag else ""
    prog = bindir / "claude"
    prog.write_text(f"#!/bin/sh\nprintf 'Usage: claude [options]\\n  -p, --print\\n{text}'\nexit 0\n")
    prog.chmod(0o755)
    return bindir


def _env(tmp_path, bindir=None):
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    path = [str(pathlib.Path(sys.executable).parent), "/usr/bin", "/bin"]
    if bindir is not None:
        path.insert(0, str(bindir))
    return {"HOME": str(home), "PATH": os.pathsep.join(path),
            "PYTHONPATH": str(REPO), "LANG": "C.UTF-8"}


def test_an_old_cli_program_is_not_a_model(tmp_path, monkeypatch):
    from delfin.agent import model_setup

    old = _program(tmp_path / "old", knows_flag=False)
    monkeypatch.setenv("PATH", str(old))
    problem = model_setup.cli_program_problem()
    assert problem and "too old" in problem and str(old / "claude") in problem

    new = _program(tmp_path / "new", knows_flag=True)
    monkeypatch.setenv("PATH", str(new))
    assert model_setup.cli_program_problem() == ""

    monkeypatch.setenv("PATH", str(tmp_path / "nothing-here"))
    assert "not on PATH" in model_setup.cli_program_problem()


def test_the_dashboard_judges_the_program_from_its_files_without_starting_it(tmp_path):
    """Building the agent tab starts no program (the suite's agent-CLI
    tripwire caught a --help probe there)."""
    from delfin.agent import model_setup

    marker = tmp_path / "was-run"
    npm = tmp_path / "lib" / "node_modules" / "cli"
    npm.mkdir(parents=True)
    (npm / "package.json").write_text('{"version": "1.0.120"}')
    (npm / "cli.js").write_text(f"#!/bin/sh\ntouch {marker}\n")
    (npm / "cli.js").chmod(0o755)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "claude").symlink_to(npm / "cli.js")
    problem = model_setup.cli_program_problem(str(tmp_path / "bin" / "claude"), run_help=False)
    assert "1.0.120" in problem and "too old" in problem

    native = tmp_path / "share" / "versions" / "2.1.322"
    native.parent.mkdir(parents=True)
    native.write_text(f"#!/bin/sh\ntouch {marker}\n")
    native.chmod(0o755)
    assert model_setup.cli_program_problem(str(native), run_help=False) == ""
    assert not marker.exists(), "a program was started"


def test_the_selection_needs_its_own_key_or_program(tmp_path, monkeypatch):
    from delfin.agent import model_setup

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    for name in ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("PATH", str(_program(tmp_path, knows_flag=False)))

    assert model_setup.selection_problem("cli", "", "sonnet")
    assert "KIT_TOOLBOX_API_KEY" in model_setup.selection_problem("api", "kit", "kit.x")
    assert model_setup.selection_problem("api", "ollama", "q:7b") == ""
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "k")
    assert model_setup.selection_problem("api", "kit", "kit.x") == ""
    assert model_setup.connected_providers() == ["kit"]


def test_the_terminal_says_no_model_is_set_up_instead_of_a_backend_error(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    done = subprocess.run(
        [sys.executable, "-c", "import sys; from delfin.agent.cli import main; sys.exit(main())",
         "chat", "-p", "say hi"],
        cwd=work, env=_env(tmp_path, _program(tmp_path, knows_flag=False)),
        capture_output=True, text=True, timeout=180, stdin=subprocess.DEVNULL)

    out = done.stdout + done.stderr
    assert done.returncode == 2, out
    assert "No language model is set up" in out
    assert "delfin-agent credentials set KIT_TOOLBOX_API_KEY" in out
    assert "unknown option" not in out


def test_the_card_says_how_to_connect_a_model():
    from delfin.agent import model_setup

    card = model_setup.no_model_card_html("")
    assert "delfin-agent-no-model" in card
    for needed in ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "OLLAMA_HOST",
                   "delfin-agent credentials set", "reload"):
        assert needed in card, needed
    assert "too old" in model_setup.no_model_card_html("the program is too old")


def _agent_spec():
    tree = ast.parse((REPO / "delfin" / "dashboard" / "__init__.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict):
            keys = [k.value for k in node.keys if isinstance(k, ast.Constant)]
            values = dict(zip(keys, node.values))
            if isinstance(values.get("id"), ast.Constant) and values["id"].value == "agent" \
                    and "default_visible" in values:
                return values
    raise AssertionError("no tab spec with id 'agent'")


def test_the_agent_tab_is_shown_without_a_model():
    spec = _agent_spec()
    visible = spec["default_visible"]
    assert isinstance(visible, ast.Constant) and visible.value is True, \
        "the DELFIN Agent tab depends on a detected backend"


def test_the_agent_tab_shows_the_card_when_nothing_is_connected():
    source = (REPO / "delfin" / "dashboard" / "tab_agent.py").read_text()
    assert "no_model_card_html(" in source


def test_the_card_takes_the_place_of_the_empty_chat_not_a_row_of_its_own():
    """Measured 2026-11-13 at 1366x768: a card above the chat pushed the
    input row below the fold. It is the empty chat's content instead: both
    places that draw an empty chat draw it, and the tab's root rows do not
    grow by one."""
    source = (REPO / "delfin" / "dashboard" / "tab_agent.py").read_text()
    assert source.count("_empty_chat_html()") >= 2, "both empty-chat sites"
    assert "*([_no_model_card]" not in source, "the card is still a row of the tab"


def test_the_card_is_compact():
    from delfin.agent import model_setup

    card = model_setup.no_model_card_html("x")
    assert "max-height" in card and "overflow-y:auto" in card
