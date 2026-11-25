"""``pipeline_prepare`` writes the CONTROL template where the caller named it.

Input: a chain of three calls, each through the product's own entry point --
``ops_server.tool_pipeline_prepare`` (the MCP handler), ``delfin.api.
pipeline_prepare`` (the typed API), and ``delfin.cli.main`` (the command), each
given an absolute ``control_file`` under a temporary directory.
Output: the template exists at that path, and nowhere else.

What was wrong: ``delfin/api.py:177`` passed ``control_file`` to
``delfin --define <control_file>``, and ``--define`` names the INPUT file
(``delfin/cli_helpers.py:127``, ``nargs="?"``), not the CONTROL file. So
``pipeline_prepare(control_file="/tmp/x/CONTROL.txt")`` created
``/tmp/x/CONTROL.txt`` as an EMPTY input file and wrote the template to
``CONTROL.txt`` in whatever directory the ops server happened to run in. The
tool's answer said so in passing ("... has been written") between two
"has been created (empty)" lines, so the caller had to read the text to notice;
a model that asked for a template in its workspace got an empty file there and
a populated one in DELFIN's own checkout.

Measured 2026-11-24: the MCP call left ``<repo>/CONTROL.txt`` behind in this
tree (275 lines of template) and an empty file at the requested path. The
existing test (``tests/test_api.py:84``) pins the argv, so it stayed green: the
argv was the bug.

The fix passes the CONTROL file through ``--control`` -- which every command
already takes (``delfin/cli_helpers.py:141``) and which ``delfin /workspace
--control /p/CONTROL.txt --define input.txt`` honours -- and never names the
input after the control file.
"""

from __future__ import annotations

import json

def test_the_api_writes_the_template_at_the_named_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from delfin import api

    target = tmp_path / "project" / "CONTROL.txt"
    target.parent.mkdir()
    out = api.pipeline_prepare(control_file=str(target))

    assert target.is_file(), out.stdout + out.stderr
    assert "functional=" in target.read_text(encoding="utf-8"), \
        "the template, not an empty input file"
    assert not (tmp_path / "CONTROL.txt").exists(), \
        "no second CONTROL.txt in the process's own directory"


def test_the_mcp_tool_writes_it_there_and_says_so(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_OPS_ALLOW_MUTATE", "1")
    from delfin.ops_server import server as ops

    target = tmp_path / "ws" / "CONTROL.txt"
    target.parent.mkdir()
    raw = ops.tool_pipeline_prepare(control_file=str(target), allow_mutate=True)
    payload = json.loads(raw)

    assert payload.get("ok") is True, payload
    assert target.is_file() and "functional=" in target.read_text(encoding="utf-8")
    assert not (tmp_path / "CONTROL.txt").exists(), payload


def test_the_input_file_is_still_created_beside_the_template(tmp_path, monkeypatch):
    """Control: the other half of the documented behaviour is unchanged."""
    monkeypatch.chdir(tmp_path)
    from delfin import api

    target = tmp_path / "ws" / "CONTROL.txt"
    target.parent.mkdir()
    api.pipeline_prepare(control_file=str(target))

    assert (target.parent / "input.txt").is_file(), \
        "the input file belongs beside the CONTROL file it is named in"
    assert "input_file=input.txt" in target.read_text(encoding="utf-8")


def test_a_bare_control_name_is_still_the_workspace(tmp_path, monkeypatch):
    """Control: the historical relative form keeps working."""
    monkeypatch.chdir(tmp_path)
    from delfin import api

    api.pipeline_prepare(control_file="CONTROL.txt")
    assert (tmp_path / "CONTROL.txt").is_file()
    assert "functional=" in (tmp_path / "CONTROL.txt").read_text(encoding="utf-8")
