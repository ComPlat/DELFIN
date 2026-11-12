"""MLIP GOAT: which provider runs it, how many workers, one model server.

Measured on an 81-atom chromophore (docs/MLIP_GOAT_GPU.md): with orca-mlips
every GOAT worker loaded a model of its own -- 40 workers (PAL 40) ran the GPU
out of memory, 8 ran -- and every ExtOpt call started a Python process with
torch (1.0-1.5 s per step).  So:

* orca-mlips GOAT is capped at 8 workers on a large GPU (fewer on a small
  one) as ``%goat nworkers``; PAL is not changed; ``goat_max_workers``
  overrides;
* orca-external-tools (OET) runs every worker against ONE server per job
  whose model processes the GPU memory has to hold; ``auto`` prefers it when
  it is installed with the backend; ``global_opt_mlip_provider`` forces one;
* the OET server is started once before ORCA and stopped after it, and the
  ORCA input names a stable client wrapper (no port in the input).
"""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import pytest

from delfin import xtb_crest
from delfin.common.control_validator import (
    CONTROL_FIELD_SPECS,
    _as_global_opt_mlip_provider,
    _as_goat_max_workers,
)
from delfin.mlp_tools import mlip_goat, oet
from delfin.mlp_tools.mlip_goat import plan_goat

NO_OET = {"available": False, "backends": ()}
OET = {"available": True, "backends": ("aimnet2",)}
MLIPS = {"available": True}
NO_MLIPS = {"available": False}
BODY = "O 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"
MINIMUM = "3\nE -76.4\nO 0.0 0.0 0.0\nH 0.75 0.58 0.0\nH -0.75 0.58 0.0\n"


# -- the worker cap ---------------------------------------------------------------

def test_orca_mlips_goat_with_pal_40_is_capped_to_8_workers():
    plan = plan_goat("aimnet2", 40, gpu_gib=80, oet_info=NO_OET, mlips_info=MLIPS)
    assert plan.provider == "orca-mlips"
    assert plan.nworkers == 8 and plan.workers == 8
    assert "one model per GOAT worker" in plan.reason


@pytest.mark.parametrize("gib,cap", [(80, 8), (40, 8), (24, 4), (16, 4), (8, 4)])
def test_the_cap_follows_the_gpu_memory_and_stays_a_multiple_of_4(gib, cap):
    assert mlip_goat.orca_mlips_worker_cap(gib) == cap


def test_few_workers_are_not_touched():
    plan = plan_goat("aimnet2", 8, gpu_gib=80, oet_info=NO_OET, mlips_info=MLIPS)
    assert plan.nworkers is None and plan.workers == 8


def test_goat_max_workers_overrides_and_is_checked():
    plan = plan_goat("aimnet2", 40, gpu_gib=80, max_workers=16, oet_info=NO_OET,
                     mlips_info=MLIPS)
    assert plan.nworkers == 16
    with pytest.raises(ValueError, match="multiple of 4"):
        plan_goat("aimnet2", 40, gpu_gib=80, max_workers=6, oet_info=NO_OET, mlips_info=MLIPS)


# -- the provider -------------------------------------------------------------------

def test_auto_prefers_oet_with_the_backend():
    plan = plan_goat("aimnet2", 40, gpu_gib=80, oet_info=OET, mlips_info=MLIPS)
    assert plan.provider == "oet"
    # every worker queues on one server; the GPU holds its model processes
    assert plan.nworkers is None and plan.workers == 40 and plan.server_processes == 8


def test_auto_falls_back_to_orca_mlips_for_a_backend_oet_lacks():
    plan = plan_goat("orb", 40, gpu_gib=80, oet_info=OET, mlips_info=MLIPS)
    assert plan.provider == "orca-mlips" and plan.nworkers == 8


def test_forcing_a_provider():
    assert plan_goat("aimnet2", 8, gpu_gib=80, forced="orca-mlips", oet_info=OET,
                     mlips_info=MLIPS).provider == "orca-mlips"
    with pytest.raises(ValueError, match="global_opt_mlip_provider=oet"):
        plan_goat("aimnet2", 8, gpu_gib=80, forced="oet", oet_info=NO_OET, mlips_info=MLIPS)
    with pytest.raises(ValueError, match="no calculator"):
        plan_goat("orb", 8, gpu_gib=80, forced="oet", oet_info=OET, mlips_info=MLIPS)


def test_oet_server_covers_one_request_when_workers_are_capped():
    # 40 cores over 4 workers: each request asks for 10 cores (NCores), which
    # the server's -n (its core budget) must hold.
    plan = plan_goat("aimnet2", 40, gpu_gib=80, max_workers=4, oet_info=OET, mlips_info=MLIPS)
    assert plan.server_processes == 10


def test_the_control_keys():
    names = {s.name for s in CONTROL_FIELD_SPECS}
    assert {"global_opt_mlip_provider", "goat_max_workers"} <= names
    assert _as_global_opt_mlip_provider("") == "auto"
    assert _as_global_opt_mlip_provider("orca-external-tools") == "oet"
    assert _as_global_opt_mlip_provider("ORCA_MLIPS") == "orca-mlips"
    with pytest.raises(ValueError):
        _as_global_opt_mlip_provider("xtb")
    assert _as_goat_max_workers("") is None and _as_goat_max_workers("8") == 8
    for bad in ("6", "0", "-4", "many"):
        with pytest.raises(ValueError):
            _as_goat_max_workers(bad)
    notes = json.loads((Path(__file__).resolve().parents[1] / "delfin" / "common"
                        / "control_key_notes.json").read_text())["keys"]
    assert notes["global_opt_mlip_provider"]["valid_values"]
    assert notes["goat_max_workers"]["valid_values"]


# -- orca-external-tools: detection, client wrapper, one server --------------------

def _fake_oet_env(root: Path, *, server_script: str = "") -> Path:
    prefix = root / ".mamba_env" / "oet"
    (prefix / "bin").mkdir(parents=True)
    py = prefix / "bin" / "python"
    py.write_text("#!/bin/sh\n")
    py.chmod(0o755)
    site = prefix / "lib" / "python3.12" / "site-packages"
    for name, version in (("oet", "2.1.0"), ("aimnet", "0.2.6")):
        meta = site / f"{name}-{version}.dist-info"
        meta.mkdir(parents=True)
        meta.joinpath("METADATA").write_text(
            f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")
    client = prefix / "bin" / "oet_client"
    client.write_text("#!/bin/sh\n")
    client.chmod(0o755)
    server = prefix / "bin" / "oet_server"
    server.write_text(server_script or "#!/bin/sh\n")
    server.chmod(0o755)
    return prefix


@pytest.fixture()
def ai_root(tmp_path, monkeypatch):
    root = tmp_path / "ai"
    root.mkdir()
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(root))
    monkeypatch.delenv("DELFIN_OET_PYTHON", raising=False)
    monkeypatch.delenv("DELFIN_ORCA_MLIPS_PYTHON", raising=False)
    return root


def test_oet_is_found_by_its_distribution_with_its_backends(ai_root):
    assert oet.detect_oet()["available"] is False
    prefix = _fake_oet_env(ai_root)
    info = oet.detect_oet()
    assert info["available"] and info["version"] == "2.1.0"
    assert info["backends"] == ("aimnet2",)
    assert info["client"] == str(prefix / "bin" / "oet_client")


def test_the_progext_is_a_stable_wrapper_around_the_documented_client(ai_root):
    prefix = _fake_oet_env(ai_root)
    path = oet.client_wrapper()
    assert path == str(prefix / "bin" / "delfin-oet-client")
    text = Path(path).read_text()
    assert f'exec "{prefix / "bin" / "oet_client"}" -b "$addr" "$@"' in text
    assert oet.ADDRESS_FILE in text
    assert oet.client_wrapper() == path  # unchanged, not rewritten


def test_ext_params_carry_model_and_device_per_request():
    assert oet.ext_params("aimnet2", device="cuda") == "-d cuda"
    assert oet.ext_params("aimnet2", model="aimnet2-2025", device="cpu") == "-m aimnet2-2025 -d cpu"
    with pytest.raises(ValueError, match="not an orca-external-tools method"):
        oet.ext_params("orb")


_FAKE_SERVER = textwrap.dedent("""\
    #!{python}
    import http.server, sys
    method, _, address = sys.argv[1], sys.argv[2], sys.argv[3]
    host, port = address.split(":")
    open("server_args.txt", "w").write(" ".join(sys.argv[1:]))
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200); self.end_headers(); self.wfile.write(b"ok")
        def log_message(self, *a):
            pass
    http.server.HTTPServer((host, int(port)), H).serve_forever()
    """)


def test_one_server_runs_for_the_block_and_is_stopped(ai_root, tmp_path):
    _fake_oet_env(ai_root, server_script=_FAKE_SERVER.format(python=sys.executable))
    work = tmp_path / "run"
    work.mkdir()
    with oet.oet_server("aimnet2", work, processes=8, startup_timeout=30) as address:
        assert (work / oet.ADDRESS_FILE).read_text().strip() == address
        assert oet._healthy(address)
        args = (work / "server_args.txt").read_text().split()
        assert args[0] == "aimnet2" and args[1] == "-b" and args[3:] == ["-n", "8"]
    assert not oet._healthy(address, timeout=1.0)
    assert not (work / oet.ADDRESS_FILE).exists()


# -- through CONTROL.txt ---------------------------------------------------------------

@pytest.fixture()
def run_dir(tmp_path, monkeypatch, ai_root):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("DELFIN_RECALC", raising=False)
    (tmp_path / "input.txt").write_text(BODY)
    monkeypatch.setattr(xtb_crest, "_resolve_global_job_limits",
                        lambda cores, maxcore=None: (int(cores), 1000))
    monkeypatch.setattr(xtb_crest, "_gpu_visible", lambda: True)
    monkeypatch.setattr(mlip_goat, "gpu_memory_gib", lambda: 80.0)
    return tmp_path


def _fake_orca(calls):
    def run(inp, out, **kw):
        wd = Path(kw["working_dir"])
        calls.append((Path(inp).read_text(), (wd / oet.ADDRESS_FILE).exists()))
        (wd / "XTB_GOAT.globalminimum.xyz").write_text(MINIMUM)
        Path(out).write_text("****ORCA TERMINATED NORMALLY****\n")
        return True
    return run


def test_an_oet_goat_runs_inside_one_server(run_dir, monkeypatch, ai_root):
    prefix = _fake_oet_env(ai_root, server_script=_FAKE_SERVER.format(python=sys.executable))
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    xtb_crest.XTB_GOAT(1, 0, {"xTB_method": "XTB2", "PAL": 40, "input_file": "input.txt",
                              "global_opt_method": "mlip:aimnet2"})
    text, server_up = calls[0]
    assert server_up, "ORCA ran while the server was up"
    assert f'"{prefix / "bin" / "delfin-oet-client"}"' in text
    assert '"-d cuda"' in text and "nworkers" not in text.lower()
    meta = json.loads((run_dir / "XTB2_GOAT" / "XTB_GOAT.method.json").read_text())
    assert meta["provider"] == "oet" and meta["server_processes"] == 8
    assert not (run_dir / "XTB2_GOAT" / oet.ADDRESS_FILE).exists()


def test_the_input_does_not_change_between_runs(run_dir, monkeypatch, ai_root):
    """No port in the input: smart recalc keeps a finished OET GOAT."""
    _fake_oet_env(ai_root, server_script=_FAKE_SERVER.format(python=sys.executable))
    calls = []
    monkeypatch.setattr(xtb_crest, "run_orca", _fake_orca(calls))
    cfg = {"xTB_method": "XTB2", "PAL": 8, "input_file": "input.txt",
           "global_opt_method": "mlip:aimnet2"}
    xtb_crest.XTB_GOAT(1, 0, cfg)
    (run_dir / "input.txt").write_text(BODY)
    xtb_crest.XTB_GOAT(1, 0, cfg)
    assert calls[0][0] == calls[1][0]
