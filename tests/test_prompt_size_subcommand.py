"""The prompt-size table is reachable from the CLI entry point.

``format_prompt_size_report`` rendered a dict the product computed
(``PromptLoader.prompt_size_report``) and nothing in the product called
-- it sat in the unmounted-capabilities baseline as test-only. The
mount is the ``delfin-agent prompt-size`` subcommand; these controls
drive ``delfin.agent.cli.main(argv)`` -- parser, dispatch, exit code,
stdout -- exactly as the command runs it, per the rule that a test must
go through the real entry point.
"""

from __future__ import annotations

import pytest

from delfin.agent import cli


def _main(argv, capsys):
    rc = cli.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


def test_prompt_size_prints_the_table(capsys):
    rc, out, err = _main(["prompt-size"], capsys)
    assert rc == 0, err
    # The renderer's own header lines, from the dict it renders.
    assert "role=" in out and "mode=" in out
    assert "cacheable prefix (stable)" in out
    assert "per-turn (volatile)" in out
    assert "TOTAL" in out


def test_prompt_size_reports_the_requested_role(capsys):
    rc, out, err = _main(["prompt-size", "--role", "solo_agent",
                          "--mode", "solo"], capsys)
    assert rc == 0, err
    assert "role=solo_agent" in out
    assert "mode=solo" in out


def test_prompt_size_unknown_role_fails_loudly(capsys):
    rc, out, err = _main(["prompt-size", "--role", "no_such_role"], capsys)
    assert rc != 0
    assert "ERROR" in err
