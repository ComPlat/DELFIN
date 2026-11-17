"""XTB_GOAT with the shipped CONTROL samples: input bytes and call sequence.

Pins what the classic workflow's GOAT step writes and runs when no
method choice is made: for every CONTROL sample in
tests/fixtures/control_templates (GOAT switched on, everything else as
shipped) the input is ``!<xTB_method> GOAT \\n%pal nprocs <PAL> end\\n
*xyz <charge> <mult>\\n`` + the start geometry + ``\\n*\\n``, written to
``<xTB_method>_GOAT/XTB_GOAT.inp``, and ORCA is called exactly once with
``output_XTB_GOAT.out`` in that folder.  Produced by xtb_crest.XTB_GOAT at
origin/main 3315360cf; a GOAT method choice must leave it unchanged.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin import xtb_crest
from delfin.config import parse_control_text, set_control_value

CONTROL = Path(__file__).resolve().parent / "fixtures" / "control_templates"
BODY = "O 0.000000 0.000000 0.000000\nH 0.757000 0.586000 0.000000\nH -0.757000 0.586000 0.000000\n"
MINIMUM = "3\nE -76.4\nO 0.0 0.0 0.0\nH 0.75 0.58 0.0\nH -0.75 0.58 0.0\n"


def _samples():
    return sorted(CONTROL.glob("*.txt"))


@pytest.mark.parametrize("sample", _samples(), ids=lambda p: p.name)
def test_default_goat_input_and_call_sequence(sample, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("DELFIN_RECALC", raising=False)
    text = set_control_value(sample.read_text(encoding="utf-8"), "global_optimizer", "GOAT")
    config = parse_control_text(text)
    config["input_file"] = "input.txt"
    (tmp_path / "input.txt").write_text(BODY, encoding="utf-8")
    monkeypatch.setattr(xtb_crest, "_resolve_global_job_limits",
                        lambda cores, maxcore=None: (int(cores), 1000))
    calls = []

    def fake_run_orca(inp, out, **kw):
        calls.append((Path(inp).name, Path(out).name, Path(kw["working_dir"]).name))
        (Path(kw["working_dir"]) / "XTB_GOAT.globalminimum.xyz").write_text(MINIMUM)
        return True

    monkeypatch.setattr(xtb_crest, "run_orca", fake_run_orca)
    xtb_crest.XTB_GOAT(1, 0, config)

    method = config["xTB_method"]
    folder = tmp_path / f"{method}_GOAT"
    assert calls == [("XTB_GOAT.inp", "output_XTB_GOAT.out", f"{method}_GOAT")]
    assert (folder / "XTB_GOAT.inp").read_text() == (
        f"!{method} GOAT \n%pal nprocs {int(config['PAL'])} end\n*xyz 0 1\n" + BODY + "\n*\n"
    )
    assert sorted(p.name for p in folder.iterdir()) == ["XTB_GOAT.globalminimum.xyz", "XTB_GOAT.inp"]
    # the global minimum replaces the start geometry, header stripped
    assert (tmp_path / "input.txt").read_text() == "".join(MINIMUM.splitlines(True)[2:])
