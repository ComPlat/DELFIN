"""IO1: the python -m delfin.hpc_io_check CLI — the same check a person
runs by hand. Prints findings; exit 1 on a refusal, 0 on warn/clean,
2 on usage/IO error. main(argv) is called directly so no subprocess is
needed; the printed text comes back via capsys.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data/workspace/\n"
WARN = "#!/bin/bash\nfor f in a b; do python /pfs/env/x.py; done\n"
CLEAN = "#!/bin/bash\nsleep 1\n"


def _run(monkeypatch, tmp_path, script_text):
    path = tmp_path / "job.sh"
    path.write_text(script_text, encoding="utf-8")
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    return hpc_io_check.main([str(path)])


def test_refuse_exits_1(monkeypatch, tmp_path, capsys):
    rc = _run(monkeypatch, tmp_path, WALK)
    out = capsys.readouterr().out
    assert rc == 1
    assert "tree_walk_on_lustre" in out
    assert "refuse" in out
    assert "instead:" in out


def test_warn_exits_0(monkeypatch, tmp_path, capsys):
    rc = _run(monkeypatch, tmp_path, WARN)
    out = capsys.readouterr().out
    assert rc == 0
    assert "lustre_interpreter_in_loop" in out


def test_clean_exits_0(monkeypatch, tmp_path, capsys):
    rc = _run(monkeypatch, tmp_path, CLEAN)
    assert rc == 0
    assert "no Lustre metadata-load findings" in capsys.readouterr().out


def test_missing_file_exits_2(monkeypatch, tmp_path, capsys):
    rc = hpc_io_check.main([str(tmp_path / "nope.sh")])
    assert rc == 2
    assert "cannot read" in capsys.readouterr().err


def test_no_args_exits_2_with_usage(monkeypatch, capsys):
    rc = hpc_io_check.main([])
    assert rc == 2
    assert "usage" in capsys.readouterr().out
