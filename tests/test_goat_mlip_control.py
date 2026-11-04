"""The GOAT method choice in CONTROL: an orca-mlips MLIP instead of xTB.

Keys (opt-in; absent = the classic xTB GOAT, byte for byte,
tests/test_goat_default_input_parity.py):

* ``global_opt_method`` -- empty (the xTB_method) | mlip:uma | mlip:mace |
  mlip:orb | mlip:aimnet2 (the orca-mlips plugins);
* ``global_opt_mlip_model`` -- optional model name passed as ``--model``;
* ``global_opt_mlip_device`` -- cuda (default) | cpu;
* ``refuse_noncommercial_weights`` -- yes refuses mlip:mace (MACE weights
  under the Academic Software License, LICENSE_GUARD.md family 4).

With an MLIP method XTB_GOAT writes GOAT + ExtOpt through the orca-mlips
builder, runs ORCA where the workflow runs, keeps the file names later
steps read (XTB_GOAT.inp, output_XTB_GOAT.out, XTB_GOAT.globalminimum.xyz
in <xTB_method>_GOAT), and records the method in XTB_GOAT.method.json.
device=cuda without a visible GPU stops the step. CREST stays xTB-only.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin import smart_recalc, xtb_crest
from delfin.common.control_validator import (
    CONTROL_FIELD_SPECS,
    _as_global_opt_method,
    _as_global_opt_mlip_device,
    _as_global_opt_mlip_model,
)
from delfin.config import validate_control_text

CONTROL = Path(__file__).resolve().parent / "fixtures" / "control_templates"
SAMPLE = sorted(CONTROL.glob("*.txt"))[-1]
NEW_KEYS = ("global_opt_method", "global_opt_mlip_model", "global_opt_mlip_device",
            "refuse_noncommercial_weights")
BODY = "O 0.000000 0.000000 0.000000\nH 0.757000 0.586000 0.000000\nH -0.757000 0.586000 0.000000\n"
MINIMUM = "3\nE -76.4\nO 0.0 0.0 0.0\nH 0.75 0.58 0.0\nH -0.75 0.58 0.0\n"


# ----------------------------------------------------------------- validator


@pytest.mark.parametrize("raw,want", [
    ("", ""), (None, ""), ("none", ""), ("xtb", ""), ("XTB", ""),
    ("mlip:uma", "mlip:uma"), ("MLIP:UMA", "mlip:uma"), ("mlip: mace", "mlip:mace"),
    ("mlip:orb", "mlip:orb"), ("mlip:aimnet2", "mlip:aimnet2"),
])
def test_method_values(raw, want):
    assert _as_global_opt_method(raw) == want


@pytest.mark.parametrize("raw", ["mlip:ani2x", "uma", "mlip:", "dft", "mlip:uma;rm"])
def test_method_refusals(raw):
    with pytest.raises(ValueError, match="mlip:uma"):
        _as_global_opt_method(raw)


def test_device_and_model():
    assert _as_global_opt_mlip_device("") == "cuda"
    assert _as_global_opt_mlip_device("CPU") == "cpu"
    with pytest.raises(ValueError, match="cuda"):
        _as_global_opt_mlip_device("gpu")
    assert _as_global_opt_mlip_model("") == ""
    assert _as_global_opt_mlip_model("uma-s-1p1") == "uma-s-1p1"
    for bad in ("a b", 'x"y', "m\nend", "x" * 200):
        with pytest.raises(ValueError, match="model"):
            _as_global_opt_mlip_model(bad)


def test_keys_are_in_the_schema_and_documented():
    names = {s.name for s in CONTROL_FIELD_SPECS}
    notes = json.loads((Path(__file__).resolve().parents[1] / "delfin" / "common"
                        / "control_key_notes.json").read_text())["keys"]
    for key in NEW_KEYS:
        assert key in names
        assert key in notes and notes[key]["valid_values"]


def test_template_stays_opt_in():
    from delfin import define

    for key in NEW_KEYS:
        assert key not in getattr(define, "CONTROL_TEMPLATE", "")


def _errors(*lines):
    return validate_control_text(SAMPLE.read_text(encoding="utf-8") + "\n" + "\n".join(lines) + "\n")


def _new_errors(*lines):
    base = set(_errors())
    return [e for e in _errors(*lines) if e not in base]


def test_mlip_with_goat_validates():
    assert _new_errors("global_optimizer=GOAT", "global_opt_method=mlip:uma",
                       "global_opt_mlip_device=cuda", "global_opt_mlip_model=uma-s-1p1") == []


def test_mlip_with_crest_is_refused():
    errs = _new_errors("global_optimizer=CREST", "global_opt_method=mlip:uma")
    assert any("CREST" in e and "global_opt_method" in e for e in errs)


def test_bad_values_are_errors():
    assert any("global_opt_method" in e for e in _new_errors("global_opt_method=mlip:foo"))
    assert any("global_opt_mlip_device" in e for e in _new_errors("global_opt_mlip_device=tpu"))


def test_refuse_noncommercial_refuses_mace_only():
    errs = _new_errors("global_optimizer=GOAT", "global_opt_method=mlip:mace",
                       "refuse_noncommercial_weights=yes")
    assert any("non-commercial" in e.lower() for e in errs)
    assert _new_errors("global_optimizer=GOAT", "global_opt_method=mlip:orb",
                       "refuse_noncommercial_weights=yes") == []


# ----------------------------------------------------------------- XTB_GOAT


@pytest.fixture()
def run_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("DELFIN_RECALC", raising=False)
    (tmp_path / "input.txt").write_text(BODY, encoding="utf-8")
    monkeypatch.setattr(xtb_crest, "_resolve_global_job_limits",
                        lambda cores, maxcore=None: (int(cores), 1000))
    return tmp_path


def _fake_orca(calls):
    def run(inp, out, **kw):
        calls.append((Path(inp).name, Path(out).name, Path(kw["working_dir"]).name))
        wd = Path(kw["working_dir"])
        (wd / "XTB_GOAT.globalminimum.xyz").write_text(MINIMUM)
        Path(out).write_text("****ORCA TERMINATED NORMALLY****\n")
        return True
    return run


def _config(**kw):
    base = {"xTB_method": "XTB2", "PAL": 8, "input_file": "input.txt"}
    base.update(kw)
    return base


def test_mlip_goat_writes_extopt_input_and_keeps_file_names(run_dir, monkeypatch):
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: True)
    xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:uma",
                                     global_opt_mlip_model="uma-s-1p1"))
    folder = run_dir / "XTB2_GOAT"
    assert calls == [("XTB_GOAT.inp", "output_XTB_GOAT.out", "XTB2_GOAT")]
    inp = (folder / "XTB_GOAT.inp").read_text()
    low = inp.lower()
    assert "goat" in low and "extopt" in low and "xtb2" not in low
    assert 'progext "orca-mlips-uma"' in low
    assert "--backend uma --model uma-s-1p1 --device cuda" in inp
    assert "nprocs 8" in low
    assert "*xyzfile 0 1 XTB_GOAT_start.xyz" in inp
    start = (folder / "XTB_GOAT_start.xyz").read_text().splitlines()
    assert start[0] == "3" and start[2].split()[0] == "O"
    # downstream: the global minimum replaced the start geometry as before
    assert (run_dir / "input.txt").read_text() == "".join(MINIMUM.splitlines(True)[2:])
    meta = json.loads((folder / "XTB_GOAT.method.json").read_text())
    assert meta["global_opt_method"] == "mlip:uma"
    assert meta["progext"] == "orca-mlips-uma"
    assert meta["device"] == "cuda" and meta["model"] == "uma-s-1p1"
    assert meta["energies"] == "mlip-screen"


def test_cuda_without_gpu_stops_before_orca(run_dir, monkeypatch):
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: False)
    with pytest.raises(RuntimeError, match="GPU"):
        xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:orb"))
    assert calls == []
    assert (run_dir / "input.txt").read_text() == BODY


def test_cpu_device_needs_no_gpu(run_dir, monkeypatch):
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: False)
    xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:aimnet2",
                                     global_opt_mlip_device="cpu"))
    assert len(calls) == 1
    assert "--device cpu" in (run_dir / "XTB2_GOAT" / "XTB_GOAT.inp").read_text()


def test_noncommercial_refusal_holds_at_run_time(run_dir, monkeypatch):
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: True)
    with pytest.raises(RuntimeError, match="non-commercial"):
        xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:mace",
                                         refuse_noncommercial_weights="yes"))
    assert calls == []


def test_gpu_visibility_reads_cuda_visible_devices(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    assert xtb_crest._gpu_visible() is False
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "-1")
    assert xtb_crest._gpu_visible() is False


def test_recalc_reruns_when_the_method_changes(run_dir, monkeypatch):
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: True)
    calls = []
    fake = _fake_orca(calls)

    def run_and_fingerprint(inp, out, **kw):
        fake(inp, out, **kw)
        smart_recalc.store_fingerprint(Path(inp))
        return True

    monkeypatch.setattr(xtb_crest, "run_orca", run_and_fingerprint)
    xtb_crest.XTB_GOAT(1, 0, _config())
    (run_dir / "input.txt").write_text(BODY)
    xtb_crest.XTB_GOAT(1, 0, _config())
    assert len(calls) == 1  # unchanged xTB GOAT: skipped
    (run_dir / "input.txt").write_text(BODY)
    xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:uma"))
    assert len(calls) == 2  # the method changed: GOAT runs again
    (run_dir / "input.txt").write_text(BODY)
    xtb_crest.XTB_GOAT(1, 0, _config(global_opt_method="mlip:uma"))
    assert len(calls) == 2  # the same MLIP run: skipped


def test_report_names_the_mlip():
    from delfin.reporting.delfin_docx_report import _global_opt_method_label

    assert _global_opt_method_label({"xTB_method": "XTB2"}) == "GFN2-xTB"
    label = _global_opt_method_label({"xTB_method": "XTB2", "global_opt_method": "mlip:uma"})
    assert "UMA" in label and "orca-mlips" in label
