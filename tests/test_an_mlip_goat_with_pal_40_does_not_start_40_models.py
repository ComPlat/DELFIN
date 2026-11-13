"""A PAL-40 CONTROL.txt with global_opt_method=mlip:aimnet2 must not start 40 models.

Measured (docs/MLIP_GOAT_GPU.md): GOAT runs one worker per process, and with
orca-mlips every worker loaded its own model server -- 40 of them ran the
GPU out of memory (CUDA allocation failures, the job OOM-killed) while 8 ran.
The orca-mlips GOAT input now carries ``%goat nworkers 8`` and ``%pal nprocs 8``;
the job keeps its 40 cores.
Imports only what the old code had, so this file is red before the fix.
"""

from __future__ import annotations

import json
from pathlib import Path

from delfin import xtb_crest

BODY = "O 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"
MINIMUM = "3\nE -76.4\nO 0.0 0.0 0.0\nH 0.75 0.58 0.0\nH -0.75 0.58 0.0\n"


def test_a_pal_40_orca_mlips_goat_writes_8_workers(tmp_path, monkeypatch):
    import delfin.mlp_tools as mlp_tools
    import delfin.mlp_tools.orca_mlips as om

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("DELFIN_RECALC", raising=False)
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path / "ai"))
    monkeypatch.delenv("DELFIN_OET_PYTHON", raising=False)
    (tmp_path / "input.txt").write_text(BODY)
    monkeypatch.setattr(xtb_crest, "_resolve_global_job_limits",
                        lambda cores, maxcore=None: (int(cores), 1000))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: True)
    monkeypatch.setattr(om, "detect_orca_mlips", lambda *a, **k: {"available": True})
    goat = getattr(mlp_tools, "mlip_goat", None)
    if goat is None:
        try:
            import delfin.mlp_tools.mlip_goat as goat
        except ImportError:
            goat = None
    if goat is not None:
        monkeypatch.setattr(goat, "gpu_memory_gib", lambda: 80.0)
    seen = []

    def fake_orca(inp, out, **kw):
        seen.append(Path(inp).read_text().lower())
        (Path(kw["working_dir"]) / "XTB_GOAT.globalminimum.xyz").write_text(MINIMUM)
        Path(out).write_text("****ORCA TERMINATED NORMALLY****\n")
        return True

    monkeypatch.setattr(xtb_crest, "run_orca", fake_orca)
    xtb_crest.XTB_GOAT(1, 0, {"xTB_method": "XTB2", "PAL": 40, "input_file": "input.txt",
                              "global_opt_method": "mlip:aimnet2"})
    assert "nprocs 8" in seen[0]  # one process per capped worker; PAL stays 40
    assert "nworkers 8" in seen[0]
    meta = json.loads((tmp_path / "XTB2_GOAT" / "XTB_GOAT.method.json").read_text())
    assert meta["goat_workers"] == 8
