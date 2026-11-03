"""``python -m delfin.orca_knowledge probe KEYWORD``: ask ORCA's own parser
whether a word is a simple-input keyword. Tests use a fake ``orca`` script
(ORCA_BINARY); the real ORCA is never started here."""
from __future__ import annotations

import os
import textwrap

import pytest

from delfin.orca_knowledge import __main__ as okm
from delfin.orca_knowledge import probe


@pytest.fixture(autouse=True)
def no_real_orca(monkeypatch):
    """No test may reach the real ORCA: without an explicit fake the probe
    must answer 'unavailable', never start /opt/orca*/orca."""
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", None)
    monkeypatch.delenv("ORCA_BINARY", raising=False)
    monkeypatch.setattr(probe.shutil, "which", lambda name: None)
    monkeypatch.setattr(probe.glob, "glob", lambda pat: [])
    yield


def _fake_orca(tmp_path, body):
    """A script pretending to be ORCA; it reads the input file path as argv[1]."""
    script = tmp_path / "fakebin" / "orca"
    script.parent.mkdir(exist_ok=True)
    script.write_text(textwrap.dedent(body))
    script.chmod(0o755)
    return str(script)


def _input_written(run_dir):
    inp = run_dir / "probe.inp"
    return inp.read_text()


@pytest.fixture
def run_dir(tmp_path):
    return tmp_path


def test_rejected_when_orca_prints_unrecognised(run_dir, tmp_path, capsys):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        cat "$1"
        echo
        echo "UNRECOGNIZED OR DUPLICATED KEYWORD"
        echo "The first unknown keyword was $2"
        exit 1
    """)
    result = probe.probe_one("RKS", run_dir)
    assert result == "rejected"
    text = _input_written(run_dir)
    assert "! RKS" in text
    assert "* xyz" in text and "H" in text  # geometry on H2


def test_reject_message_printed_once(tmp_path, capsys):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "UNRECOGNIZED OR DUPLICATED KEYWORD"
    """)
    assert okm.main(["probe", "RKS"]) == 0
    out, _ = capsys.readouterr()
    assert out == "RKS: rejected\n"


def test_accepted_on_normal_header(run_dir, tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        cat "$1"
        echo
        echo "  * Program     ORCA"
        echo "  * Version     6.1.1"
        echo "  # of atoms 2"
        exit 0
    """)
    assert probe.probe_one("RKS", run_dir) == "accepted"
    assert probe.probe_one("B3LYP", run_dir) == "accepted"


def test_cli_prints_one_line_per_word(tmp_path, capsys):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        grep -q '! NOTAKEYWORD' "$1" \\
          && { echo "UNRECOGNIZED OR DUPLICATED KEYWORD"; exit 1; } \\
          || { echo "  * Program     ORCA"; exit 0; }
    """)
    rc = okm.main(["probe", "RI", "NOTAKEYWORD"])
    out, _ = capsys.readouterr()
    assert out == "RI: accepted\nNOTAKEYWORD: rejected\n"
    assert rc == 0


def test_no_orca_binary(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("ORCA_BINARY", raising=False)
    monkeypatch.setattr(probe, "find_orca", lambda: None)
    rc = okm.main(["probe", "RKS"])
    out, _ = capsys.readouterr()
    assert out == "RKS: unavailable\n"
    assert rc == 2


def test_env_variable_selects_fake_orca(run_dir, tmp_path, monkeypatch):
    script = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "UNRECOGNIZED OR DUPLICATED KEYWORD"
    """)
    monkeypatch.setenv("ORCA_BINARY", script)
    probe.ORCA_BINARY_CACHE = None
    assert probe.find_orca() == script
    assert probe.probe_one("XYZQ", run_dir) == "rejected"


def test_run_in_fresh_temp_dir(run_dir, tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "cwd: $(pwd)" > orca_probe_seen.txt
        echo "  * Program     ORCA"
    """)
    assert probe.probe_one("SCF", run_dir) == "accepted"
    assert (run_dir / "orca_probe_seen.txt").read_text() != ""
    # and the seen cwd was the fresh run dir, not the caller's
    assert str(run_dir) in (run_dir / "orca_probe_seen.txt").read_text()


def test_timeout_kills_orca_process_group(run_dir, tmp_path):
    """A hanging fake ORCA is killed within the budget, together with its group."""
    pidfile = tmp_path / "sleeper.pid"
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, f"""
        #!/bin/sh
        echo $$ > {pidfile}
        echo "  * Program     ORCA"   # header first ...
        sleep 60                       # ... then hang instead of finishing
    """)
    result = probe.probe_one("RKS", run_dir, timeout_s=2)
    assert result == "unavailable"  # no parse verdict within the budget
    pid = int(pidfile.read_text().split()[-1])
    with pytest.raises((ProcessLookupError, OSError)):
        os.kill(pid, 0)  # it is gone


def test_probe_is_serial_not_parallel(run_dir, tmp_path, monkeypatch):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        sleep 0.05
        echo "  * Program     ORCA"
    """)
    order = []
    real = probe.probe_one
    def spy(word, run_dir, timeout_s=probe.DEFAULT_TIMEOUT_S):
        result = real(word, run_dir, timeout_s)
        order.append(word)
        return result
    monkeypatch.setattr(probe, "probe_one", spy)
    probe.probe_all(["A", "B", "C"], run_dir)
    assert order == ["A", "B", "C"]


def test_a_word_that_could_add_input_lines_is_refused_unrun(run_dir, tmp_path):
    marker = tmp_path / "ran"
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, f"""
        #!/bin/sh
        touch {marker}
    """)
    assert probe.probe_one("B3LYP\n%pal nprocs 384 end", run_dir) == "rejected"
    assert probe.probe_one("", run_dir) == "rejected"
    assert not marker.exists()


def test_a_valid_method_that_runs_past_the_budget_is_accepted(run_dir, tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                       INPUT FILE"
        sleep 60
    """)
    assert probe.probe_one("DLPNO-CCSD(T)", run_dir, timeout_s=2) == "accepted"
