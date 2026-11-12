"""`delfin-agent --version` and `delfin-voila --version` print the package
version and exit 0, as git, gh and kubectl do (bug hunt 11.11.2026, #14).

The delfin-agent case needs the cli.py hunk of the copy pair
(.gate/px/new_cli.py).
"""

from __future__ import annotations

import pytest

from delfin import __version__


def test_delfin_agent_version(capsys):
    from delfin.agent import cli
    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"delfin-agent {__version__}"


def test_delfin_agent_version_is_not_routed_into_a_chat():
    from delfin.agent import cli
    parser = cli.build_parser()
    assert cli._route_argv(["--version"], cli._subcommand_names(parser)) == [
        "--version"]


def test_delfin_voila_version(capsys):
    from delfin import cli_voila
    with pytest.raises(SystemExit) as exc:
        cli_voila.main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"delfin-voila {__version__}"
