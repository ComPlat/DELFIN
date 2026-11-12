"""``probe-input`` and ``probe-queue``: ORCA decides on a full input.

All tests use a fake ``orca`` script (ORCA_BINARY_CACHE); the real ORCA is
never started here.
"""
from __future__ import annotations

import json
import textwrap

import pytest

from delfin.orca_knowledge import __main__ as okm
from delfin.orca_knowledge import probe


@pytest.fixture(autouse=True)
def no_real_orca(monkeypatch):
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", None)
    monkeypatch.delenv("ORCA_BINARY", raising=False)
    monkeypatch.setattr(probe.shutil, "which", lambda name: None)
    # find_orca falls back to glob over /opt/orca*/orca; patch it directly
    # (patching glob.glob would hit the real glob module, breaking probe-queue).
    monkeypatch.setattr(probe, "find_orca", lambda: None)
    yield


def _fake_orca(tmp_path, body):
    script = tmp_path / "fakebin" / "orca"
    script.parent.mkdir(exist_ok=True)
    script.write_text(textwrap.dedent(body))
    script.chmod(0o755)
    return str(script)


GOOD_INP = """\
! B3LYP def2-SVP
* xyz 0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.74
*
"""


def test_finished_verdict(tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        echo "        ORCA TERMINATED NORMALLY"
    """)
    result = probe.probe_input(GOOD_INP, run_dir=tmp_path / "run")
    assert result["verdict"] == "finished"
    assert result["messages"] == []


@pytest.mark.parametrize("marker,expected_lines", [
    ("UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE\n  BOGUSWORD",
     None),
    ("Unknown identifier in SCF block line    2  :\n Last token : BOGUSOPT",
     None),
    ("Unrecognized symbol in SCF block line    3  :\n Last token : *", None),
    ("ERROR: expect a '$', '!', '%', '*' or '[' in the input", None),
])
def test_reject_markers(tmp_path, marker, expected_lines):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, f"""
        #!/bin/sh
        echo "                               INPUT FILE"
        echo "{marker}"
    """)
    result = probe.probe_input(GOOD_INP)
    assert result["verdict"] == "rejected"
    assert result["messages"]


def test_parsed_but_not_finished(tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        echo "SCF running..."
    """)
    result = probe.probe_input(GOOD_INP)
    assert result["verdict"] == "parsed"


def test_timeout_before_parse(tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        sleep 60
    """)
    result = probe.probe_input(GOOD_INP, timeout_s=2)
    assert result["verdict"] == "timeout"
    assert result["timed_out"] is True


def test_timeout_after_parse_is_parsed(tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        sleep 60
    """)
    result = probe.probe_input(GOOD_INP, timeout_s=2)
    assert result["verdict"] == "parsed"
    assert result["timed_out"] is True


def test_unavailable_without_orca():
    result = probe.probe_input(GOOD_INP)
    assert result["verdict"] == "unavailable"


def test_guard_refuses_oversized_inputs(tmp_path, capsys):
    ran = tmp_path / "ran"
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, f"#!/bin/sh\ntouch {ran}")
    seven = GOOD_INP.replace(
        "H 0.0 0.0 0.74\n*\n",
        "H 0.0 0.0 0.74\nH 0.0 0.0 0.9\nH 0.0 0.0 1.1\nH 0.0 0.0 1.3\n"
        "H 0.0 0.0 1.5\nH 0.0 0.0 1.7\n*\n")
    assert len([line for line in seven.splitlines()
                if line.startswith("H ")]) == 7
    result = probe.probe_input(seven)
    assert result["verdict"] == "refused"
    assert "probe limit" in result["messages"][0]
    for bad in ("! B3LYP\n%pal nprocs 4 end\n" + GOOD_INP.split("\n", 1)[1],
                "! B3LYP\n%maxcore 4000\n" + GOOD_INP.split("\n", 1)[1]):
        result = probe.probe_input(bad)
        assert result["verdict"] == "refused", bad
    assert not ran.exists()  # ORCA never started


def test_guard_accepts_within_limits(tmp_path):
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "        ORCA TERMINATED NORMALLY"
    """)
    small = "! B3LYP\n%pal nprocs 2\n  nprocs 2\nend\n%maxcore 1000\n" + GOOD_INP
    assert probe.probe_input(small)["verdict"] == "finished"


def test_look_for_reports_hits_with_context(tmp_path):
    out = ("  * Program     ORCA\n"
           "TCutPNO  =   33.833573141\n"
           "some more text\n"
           "TCutPNO  =   1.0d-2\n")
    probe.ORCA_BINARY_CACHE = _fake_orca(tmp_path, f"""
        #!/bin/sh
        cat <<'EOF'
{out}
EOF
    """)
    result = probe.probe_input(GOOD_INP, look_for=["TCutPNO", "Nope"])
    hits = {h["text"]: h for h in result["look_for"]}
    assert hits["TCutPNO"]["found"] is True
    assert hits["TCutPNO"]["count"] == 2
    assert any("TCutPNO" in line for line in hits["TCutPNO"]["context"])
    assert hits["Nope"]["found"] is False
    assert hits["Nope"]["count"] == 0


# ---------------- CLI ----------------

def test_cli_probe_input_rejected_exit_code(tmp_path, monkeypatch, capsys):
    inp = tmp_path / "bad.inp"
    inp.write_text(GOOD_INP)
    script = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        echo "UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE"
    """)
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", script)
    monkeypatch.chdir(tmp_path)
    rc = okm.main(["probe-input", str(inp)])
    out, _ = capsys.readouterr()
    assert rc == 1
    assert "verdict: rejected" in out
    assert "UNRECOGNIZED" in out


def test_cli_probe_input_look_for_flag(tmp_path, monkeypatch, capsys):
    inp = tmp_path / "good.inp"
    inp.write_text(GOOD_INP)
    script = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        echo "        ORCA TERMINATED NORMALLY"
        echo "TCutPNO  =   33.833573141"
    """)
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", script)
    rc = okm.main(["probe-input", str(inp), "--look-for", "TCutPNO"])
    out, _ = capsys.readouterr()
    assert rc == 0
    assert "look-for 'TCutPNO': FOUND (1x)" in out
    assert "TCutPNO  =   33.833573141" in out


def test_cli_probe_input_missing_file(capsys):
    rc = okm.main(["probe-input", "does-not-exist.inp"])
    out, err = capsys.readouterr()
    assert rc == 2
    assert "cannot read" in err


def test_cli_probe_queue_runs_requests_serially(tmp_path, monkeypatch, capsys):
    script = _fake_orca(tmp_path, """
        #!/bin/sh
        echo "                               INPUT FILE"
        case "$(cat "$1")" in
          *TCutPNO*) echo "TCutPNO  =   33.833573141"; echo "ORCA TERMINATED NORMALLY";;
          *) echo "ORCA TERMINATED NORMALLY";;
        esac
    """)
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", script)
    d = tmp_path / "queue"
    d.mkdir()
    (d / "req_a.json").write_text(json.dumps({
        "card_id": "properties.tcutpno", "question": "is TCutPNO printed?",
        "input": GOOD_INP + "TCutPNO\n", "look_for": ["TCutPNO"],
        "max_seconds": 10}))
    (d / "req_b.json").write_text(json.dumps({
        "card_id": "other", "question": "does it finish?",
        "input": GOOD_INP, "look_for": []}))
    rc = okm.main(["probe-queue", str(d)])
    out, _ = capsys.readouterr()
    assert rc == 0
    res_a = json.loads((d / "req_a.result.json").read_text())
    res_b = json.loads((d / "req_b.result.json").read_text())
    assert res_a["verdict"] == "finished"
    assert res_a["card_id"] == "properties.tcutpno"
    assert res_a["look_for"][0]["found"] is True
    assert res_b["verdict"] == "finished"
    assert "req_a.json: finished" in out
    assert "req_b.json: finished" in out


def test_cli_probe_queue_result_files_not_reprobed(tmp_path, monkeypatch,
                                                   capsys):
    """A second run must not re-run the .result.json files it wrote."""
    script = _fake_orca(tmp_path, "#!/bin/sh\necho ok")
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", script)
    d = tmp_path / "queue"
    d.mkdir()
    (d / "req.json").write_text(json.dumps({
        "card_id": "x", "input": GOOD_INP}))
    okm.main(["probe-queue", str(d)])
    capsys.readouterr()
    # result file exists now; a re-run reads only req.json
    rc = okm.main(["probe-queue", str(d)])
    out, _ = capsys.readouterr()
    assert rc == 0
    assert "req.result.json: " not in out


def test_cli_probe_queue_empty_dir(tmp_path, capsys):
    d = tmp_path / "empty"
    d.mkdir()
    rc = okm.main(["probe-queue", str(d)])
    out, err = capsys.readouterr()
    assert rc == 2
    assert "no probe requests" in err


def test_cli_probe_queue_guards_oversized(tmp_path, monkeypatch, capsys):
    script = _fake_orca(tmp_path, "#!/bin/sh\necho hi")
    monkeypatch.setattr(probe, "ORCA_BINARY_CACHE", script)
    d = tmp_path / "q"
    d.mkdir()
    big = GOOD_INP + "\n".join(f"H 0.0 0.0 {1.0 + i * 0.1}" for i in range(7))
    (d / "big.json").write_text(json.dumps({"card_id": "big",
                                            "input": "! B3LYP\n" + big}))
    rc = okm.main(["probe-queue", str(d)])
    res = json.loads((d / "big.result.json").read_text())
    assert res["verdict"] == "refused"
    assert rc == 0  # a refused probe is an answer, not a failure


# -- inputs that could run programs or read files are never probed -----------

@pytest.mark.parametrize("text", [
    "! ExtOpt\n%method\n ProgExt \"/bin/sh\"\nend\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
    "%compound\n Sys_cmd \"rm -rf ~\"\nend\n",
    "! HF\n%compound \"script.cmp\"\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
    "! HF def2-SVP\n* xyzfile 0 1 /etc/passwd\n",
    "! HF def2-SVP\n* gzmtfile 0 1 /etc/shadow\n",
    "! HF def2-SVP MORead\n%moinp \"/home/x/secret.gbw\"\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
    "! HF def2-SVP\n%geom\n InHessName \"/tmp/x.hess\"\nend\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
    "! HF def2-SVP\n%basis\n GTOName = \"/etc/passwd\"\nend\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
    "! HF def2-SVP\n%nosuchblock x 1 end\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n",
])
def test_an_input_that_could_run_a_program_or_read_a_file_is_refused(text):
    with pytest.raises(probe.ProbeGuardError):
        probe.guard_input(text)


def test_an_ordinary_small_probe_passes_the_guard():
    probe.guard_input("! HF def2-SVP TightSCF\n%scf\n MaxIter 50\nend\n"
                      "%mdci\n TCutPNO 3.33e-7\nend\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
