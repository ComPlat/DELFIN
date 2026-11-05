"""Guided set-up of an MLIP environment for ORCA ExtOpt via orca-mlips.

Five stages, each a function with data in and data out:

* :func:`probe` -- read-only inventory of the host: Slurm GPU partitions
  and their gres (``scontrol``/``sinfo`` through
  :mod:`delfin.scheduler_client`), the NVIDIA driver (``nvidia-smi`` only
  where it exists), the CUDA toolkit and environment modules, ORCA
  installations, a conda-family tool (micromamba/mamba/conda), existing
  environments holding orca-mlips / fairchem / mace / orb / aimnet and
  their ``orca-mlips-*`` ProgExt entry points, and the current
  ``runtime.slurm`` / ``runtime.mlp`` settings. Every part degrades to a
  stated absence; nothing raises.
* :func:`propose` -- partition, backend (with the code and weight licence
  records of :mod:`delfin.mlp_tools.registry`), and the exact commands:
  environment with Python 3.11, torch from the official PyTorch index for
  the probed CUDA version, ``pip install "orca-mlips[<backend>]"``. The
  commands are data; nothing here runs them.
* :func:`run_install` -- runs the proposed commands one by one, each only
  after ``approve(command)`` returned True. There is no default approver.
* :func:`configure` -- the settings diff for ``runtime.slurm.gpu_partitions``
  and ``runtime.mlp.*``; written through
  :func:`delfin.user_settings.save_settings` only with ``write=True``.
* :func:`prepare_verify_job` / :func:`submit_verify_job` /
  :func:`read_verify_result` -- a water energy+gradient through the
  orca-mlips command (one cold and two warm calls, which separates process
  start + model load from one gradient) and a short ExtOpt Opt, on the GPU
  partition. Preparation writes the job directory only; submission goes
  through :func:`delfin.scheduler_client.query` and the MLP job cap.

Hugging Face tokens: this module never reads, prints or stores one, and
names no token file or token variable. A gated checkpoint (UMA) is the
user's own ``hf auth login`` plus their acceptance of the model terms,
done outside DELFIN; the plan says so and the user re-probes afterwards.

Rejected: a venv route when no conda-family tool exists. The agent's
shell gate refuses ``<venv>/bin/pip install`` outside the workspace
("unsafe install target") and lets ``python -m venv <path>`` through
unasked, so that route is neither usable by the agent nor asked for each
step; micromamba/mamba/conda commands are asked in every permission mode
(``_caged_still_asked``), which is the property the install stage needs.
"""

from __future__ import annotations

import copy
import json
import math
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

from delfin import scheduler_client
from delfin.common.logging import get_logger

logger = get_logger(__name__)

#: orca-mlips plugins (console scripts orca-mlips-<name>).
ORCA_MLIPS_BACKENDS = ("uma", "orb", "mace", "aimnet2")

#: orca-mlips plugin -> registry backend that carries its licence records.
#: mace -> mace_omol: the orca-mlips MACE plugin loads MACE-OMOL weights
#: (mlp-method skill, ExtOpt providers section).
REGISTRY_BACKEND = {
    "uma": "fairchem_uma",
    "orb": "orb",
    "mace": "mace_omol",
    "aimnet2": "aimnet2",
}

DEFAULT_ENV_NAME = "delfin-mlip"
PYTHON_VERSION = "3.11"
TORCH_INDEX_URL = "https://download.pytorch.org/whl/{tag}"

#: CUDA builds the PyTorch index has published wheels for, newest first.
#: A cu12x wheel runs on any CUDA 12 driver (minor-version compatibility),
#: so the rule is: the newest build whose CUDA version does not exceed the
#: driver's. The wheel set changes with each torch release; pip takes the
#: newest torch the chosen index still carries.
_TORCH_CUDA_BUILDS = ((12, 8), (12, 6), (12, 4), (12, 1), (11, 8))
#: Used when no driver is visible on this node (a login node usually has
#: no GPU); the verify job reports the GPU node's driver.
_DEFAULT_CUDA_BUILD = (12, 6)
#: GPU generations that need at least this CUDA build.
_MIN_BUILD_FOR_GPU = (
    (re.compile(r"\b(?:gb?200|b100|b200|rtx\s*50\d\d)\b", re.I), (12, 8)),
)

#: Order of magnitude for the torch CUDA wheel and its nvidia-* runtime
#: wheels; not measured by DELFIN, stated as such in the plan.
TORCH_DOWNLOAD_ESTIMATE = "about 3 GB (torch CUDA wheel + nvidia-* runtime wheels; order of magnitude)"

#: Distribution names (normalised) the env probe looks for.
_MLIP_DISTS = ("orca_mlips", "fairchem_core", "mace_torch", "orb_models",
               "aimnet", "aimnet2calc", "torch")

#: Partition/gres names of non-NVIDIA accelerators a CUDA torch cannot use.
_NON_CUDA_RE = re.compile(r"mi\d{3}|amd|rocm|gaudi|xpu", re.I)

# Scheduler argv as module constants: the direct-call guard
# (tests/test_r3_direct_scheduler_call_guard.py) flags a scheduler command
# spelled at a call site; the runs go through scheduler_client.query.
_SCONTROL_PARTITIONS = ("scontrol", "show", "partition", "--oneliner")
_SINFO_GRES = ("sinfo", "-h", "-o", "%P|%G")
_SBATCH_PROG = "sbatch"

_PARTITION_RE = re.compile(r"^[A-Za-z0-9_.-]+(?:,[A-Za-z0-9_.-]+)*$")
_GRES_RE = re.compile(r"^[A-Za-z0-9_.]+(?::[A-Za-z0-9_.]+){0,2}$")
_TIME_RE = re.compile(r"^(?:\d+-)?\d{1,3}(?::\d{2}){0,2}$")
_ENV_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_SAFE_PATH_RE = re.compile(r"^[^\s'\"`$;&|<>\\]+$")

VERIFY_JOB_NAME = "mlp_verify"
VERIFY_TIME_LIMIT = "00:20:00"
VERIFY_RESULT = "verify_result.json"

#: Water, Angstrom (experimental geometry, r(OH) 0.9572, HOH 104.5 deg).
_WATER = (
    ("O", (0.0, 0.0, 0.1173)),
    ("H", (0.0, 0.7572, -0.4692)),
    ("H", (0.0, -0.7572, -0.4692)),
)


# --------------------------------------------------------------------- host


def _default_run(argv: Sequence[str], timeout: int = 20):
    """(returncode, stdout, stderr); (None, '', reason) when not runnable."""
    try:
        done = subprocess.run(list(argv), capture_output=True, text=True,
                              timeout=timeout)
        return done.returncode, done.stdout or "", done.stderr or ""
    except Exception as exc:  # noqa: BLE001 - reported as an absence
        return None, "", f"{type(exc).__name__}: {exc}"


def _default_modules(patterns: Sequence[str]) -> list:
    from delfin.system_tools import matching_modules
    return list(matching_modules(patterns))


def _default_orca_installs() -> list:
    from delfin.runtime_setup import discover_orca_installations
    return list(discover_orca_installations())


def _default_settings() -> dict:
    from delfin.user_settings import load_settings
    return load_settings()


def _default_current_interpreter() -> dict:
    from delfin.mlp_tools.orca_mlips import detect_orca_mlips
    return detect_orca_mlips()


@dataclass
class Host:
    """What the probe may ask the machine. Every field is replaceable, so
    a test describes a host instead of needing one."""

    which: Callable[[str], Optional[str]] = shutil.which
    run: Callable[..., tuple] = _default_run
    scheduler: Callable[[Sequence[str]], object] = field(
        default=lambda argv: scheduler_client.query(
            scheduler_client.kind_of_argv(argv), list(argv)))
    modules: Callable[[Sequence[str]], list] = _default_modules
    orca_installs: Callable[[], list] = _default_orca_installs
    settings: Callable[[], dict] = _default_settings
    current_interpreter: Callable[[], dict] = _default_current_interpreter
    env: dict = field(default_factory=lambda: dict(os.environ))


def _safe(fn, default, what: str, notes: list):
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001 - one broken part, not the probe
        notes.append(f"{what}: not probed ({type(exc).__name__}: {exc})")
        return default


# -------------------------------------------------------------------- slurm


def _parse_scontrol_partitions(text: str) -> list:
    rows = []
    for line in (text or "").splitlines():
        fields = dict(p.split("=", 1) for p in line.split() if "=" in p)
        name = fields.get("PartitionName", "")
        if not name:
            continue
        tres = fields.get("TRES", "")
        m = re.search(r"gres/gpu(?::([A-Za-z0-9_.-]+))?=(\d+)", tres)
        rows.append({
            "name": name,
            "state": fields.get("State", ""),
            "gpus": int(m.group(2)) if m else 0,
            "gpu_type": (m.group(1) or "") if m else "",
            "max_time": fields.get("MaxTime", ""),
        })
    return rows


def _parse_sinfo_gres(text: str) -> dict:
    out: dict = {}
    for line in (text or "").splitlines():
        if "|" not in line:
            continue
        part, gres = line.split("|", 1)
        part = part.strip().rstrip("*")
        gres = gres.strip()
        if part and gres and gres != "(null)":
            out.setdefault(part, set()).add(gres)
    return {k: sorted(v) for k, v in out.items()}


def _probe_slurm(host: Host, notes: list) -> dict:
    present = {c: bool(host.which(c)) for c in ("sbatch", "squeue", "scontrol", "sinfo")}
    out = {"available": present["sbatch"] and present["squeue"],
           "commands": present, "gpu_partitions": [], "source": "", "reason": ""}
    if not out["available"]:
        out["reason"] = "no Slurm on this host (sbatch/squeue not found)"
        return out
    from delfin.slurm_submit import discover_gpu_partitions

    rows: list = []
    if present["scontrol"]:
        text = host.scheduler(_SCONTROL_PARTITIONS)
        if text is scheduler_client.THROTTLED:
            out["reason"] = "scontrol throttled by the scheduler client; probe again in 25 s"
        elif text:
            names = set(discover_gpu_partitions(str(text)))
            rows = [r for r in _parse_scontrol_partitions(str(text)) if r["name"] in names]
            out["source"] = "scontrol"
        else:
            out["reason"] = "scontrol show partition gave no answer"
    gres_by_part: dict = {}
    if present["sinfo"]:
        text = host.scheduler(_SINFO_GRES)
        if text and text is not scheduler_client.THROTTLED:
            gres_by_part = _parse_sinfo_gres(str(text))
            if not rows:
                # scontrol refused or absent: GPU partitions from sinfo gres.
                rows = [{"name": p, "state": "", "gpus": 0, "gpu_type": "",
                         "max_time": ""}
                        for p, g in gres_by_part.items()
                        if any(x.startswith("gpu") for x in g)
                        and not p.startswith("dev_")]
                out["source"] = out["source"] or ("sinfo" if rows else "")
    for r in rows:
        r["gres"] = [g for g in gres_by_part.get(r["name"], []) if g.startswith("gpu")]
        label = " ".join([r["name"], r["gpu_type"], *r["gres"]])
        r["cuda_capable"] = not _NON_CUDA_RE.search(label)
    out["gpu_partitions"] = rows
    if not rows and not out["reason"]:
        out["reason"] = "Slurm reachable, but no partition reports GPUs"
    return out


# --------------------------------------------------------------------- cuda


def _version_tuple(text: str) -> Optional[tuple]:
    m = re.search(r"(\d+)\.(\d+)", text or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def _probe_cuda(host: Host, notes: list) -> dict:
    out = {"nvidia_smi": bool(host.which("nvidia-smi")), "driver": "",
           "driver_cuda": "", "gpus": [], "toolkit": "", "modules": []}
    if out["nvidia_smi"]:
        rc, text, _err = host.run(["nvidia-smi"])
        if rc == 0:
            m = re.search(r"Driver Version:\s*([\d.]+)", text)
            c = re.search(r"CUDA Version:\s*([\d.]+)", text)
            out["driver"] = m.group(1) if m else ""
            out["driver_cuda"] = c.group(1) if c else ""
        rc, text, _err = host.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"])
        if rc == 0:
            out["gpus"] = [g.strip() for g in text.splitlines() if g.strip()]
    if host.which("nvcc"):
        rc, text, _err = host.run(["nvcc", "--version"])
        m = re.search(r"release\s+([\d.]+)", text or "")
        out["toolkit"] = m.group(1) if (rc == 0 and m) else ""
    out["modules"] = _safe(lambda: host.modules(["cuda"]), [], "cuda modules", notes)
    return out


# --------------------------------------------------------------------- orca


def _probe_orca(host: Host, notes: list) -> dict:
    from delfin.runtime_setup import describe_orca_installation

    installs = []
    for path in _safe(host.orca_installs, [], "ORCA installations", notes):
        label = describe_orca_installation(path)
        ver = _version_tuple(label.replace("ORCA", ""))
        installs.append({"path": str(path), "label": label,
                         "version": ".".join(map(str, ver)) if ver else "",
                         "is_6_1_or_newer": bool(ver and ver >= (6, 1))})
    mods = _safe(lambda: host.modules(["orca"]), [], "orca modules", notes)
    mod_ok = [m for m in mods if (_version_tuple(m.split("/", 1)[-1]) or (0, 0)) >= (6, 1)]
    return {"installations": installs, "modules": mods,
            "ready": any(i["is_6_1_or_newer"] for i in installs) or bool(mod_ok),
            "modules_6_1_or_newer": mod_ok}


# -------------------------------------------------------------------- conda


_CONDA_TOOLS = ("micromamba", "mamba", "conda")


def _probe_conda(host: Host, notes: list) -> dict:
    tool, path, on_path = "", "", False
    for name in _CONDA_TOOLS:
        found = host.which(name)
        if found:
            tool, path, on_path = name, found, True
            break
    if not tool:
        for var, name in (("MAMBA_EXE", "micromamba"), ("CONDA_EXE", "conda")):
            cand = host.env.get(var, "")
            if cand and Path(cand).is_file():
                tool, path = name, cand
                break
    envs: list = []
    if tool:
        rc, text, err = host.run([path, "env", "list", "--json"])
        if rc == 0:
            try:
                envs = [str(e) for e in json.loads(text).get("envs", [])]
            except Exception:
                notes.append(f"{tool} env list: unreadable JSON")
        else:
            notes.append(f"{tool} env list failed: {(err or '').strip()[:200]}")
    return {"tool": tool, "path": path, "on_path": on_path, "envs": envs}


def _norm_dist(name: str) -> str:
    return re.sub(r"[-_.]+", "_", name).lower()


def inspect_env(prefix: str | Path) -> dict:
    """Packages and ProgExt entry points of one environment, from files.

    Reads ``site-packages/*.dist-info`` directory names and ``bin/``; runs
    no interpreter, imports nothing.
    """
    prefix = Path(prefix)
    out = {"prefix": str(prefix), "name": prefix.name, "python": "",
           "packages": {}, "progext": [], "torch_build": ""}
    lib = prefix / "lib"
    try:
        pydirs = sorted(d for d in lib.iterdir() if d.name.startswith("python3"))
    except OSError:
        pydirs = []
    for pyd in pydirs:
        out["python"] = pyd.name.replace("python", "")
        site = pyd / "site-packages"
        try:
            entries = list(site.iterdir())
        except OSError:
            continue
        for e in entries:
            if not e.name.endswith(".dist-info"):
                continue
            stem = e.name[: -len(".dist-info")]
            name, _, version = stem.partition("-")
            key = _norm_dist(name)
            if key in _MLIP_DISTS:
                out["packages"][key] = version
    torch_v = out["packages"].get("torch", "")
    if torch_v:
        m = re.search(r"\+(cu\d+|cpu|rocm[\d.]+)", torch_v)
        out["torch_build"] = m.group(1) if m else "default"
    try:
        out["progext"] = sorted(p.name for p in (prefix / "bin").iterdir()
                                if p.name.startswith("orca-mlips-"))
    except OSError:
        pass
    return out


def _probe_envs(conda: dict, configured_prefix: str) -> list:
    prefixes = list(conda.get("envs") or [])
    if configured_prefix and configured_prefix not in prefixes:
        prefixes.append(configured_prefix)
    found = []
    for p in prefixes:
        info = inspect_env(p)
        if set(info["packages"]) - {"torch"} or info["progext"]:
            info["configured"] = (p == configured_prefix)
            found.append(info)
    return found


# ----------------------------------------------------------------- settings


def _env_prefix_from_activation(activation: str) -> str:
    m = re.search(r"""PATH=["']?([^"':$]+)/bin""", activation or "")
    return m.group(1) if m else ""


def _probe_settings(host: Host, notes: list) -> dict:
    settings = _safe(host.settings, {}, "settings", notes)
    runtime = (settings.get("runtime") or {}) if isinstance(settings, dict) else {}
    slurm = runtime.get("slurm") or {}
    mlp = runtime.get("mlp") or {}
    overrides = sorted(k for k in host.env if k.startswith("DELFIN_MLP_CLUSTER_"))
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig

    cfg = MLPClusterConfig.from_sources(settings=runtime, env=host.env)
    return {
        "gpu_partitions": str(slurm.get("gpu_partitions", "") or ""),
        "mlp": dict(mlp),
        "env_overrides": overrides,
        "effective": {
            "gpu_partitions": cfg.gpu_partitions, "gres": cfg.gres,
            "time_limit": cfg.time_limit, "mem_mb": cfg.mem_mb,
            "cpus_per_task": cfg.cpus_per_task,
            "env_activation": cfg.env_activation,
        },
        "env_prefix": str(mlp.get("env_prefix", "") or "")
        or _env_prefix_from_activation(cfg.env_activation),
    }


# -------------------------------------------------------------------- probe


def probe(host: Optional[Host] = None) -> dict:
    """Read-only inventory of this host for an orca-mlips set-up.

    Output: a dict with the parts ``slurm``, ``cuda``, ``orca``, ``conda``,
    ``mlip_envs``, ``current_interpreter``, ``settings``, plus ``notes``
    (parts that could not be probed, with the reason) and ``summary``
    (short lines for a person). Installs nothing, writes nothing, submits
    nothing; a missing scheduler, GPU or tool is a stated absence.
    """
    host = host or Host()
    notes: list = []
    slurm = _safe(lambda: _probe_slurm(host, notes),
                  {"available": False, "gpu_partitions": [], "reason": "probe failed"},
                  "slurm", notes)
    cuda = _safe(lambda: _probe_cuda(host, notes), {}, "cuda", notes)
    orca = _safe(lambda: _probe_orca(host, notes), {"installations": [], "modules": [],
                                                     "ready": False}, "orca", notes)
    conda = _safe(lambda: _probe_conda(host, notes), {"tool": "", "path": "", "envs": []},
                  "conda", notes)
    settings = _safe(lambda: _probe_settings(host, notes), {}, "settings", notes)
    envs = _safe(lambda: _probe_envs(conda, settings.get("env_prefix", "")), [],
                 "MLIP environments", notes)
    current = _safe(host.current_interpreter, {}, "current interpreter", notes)
    result = {"slurm": slurm, "cuda": cuda, "orca": orca, "conda": conda,
              "mlip_envs": envs, "current_interpreter": current,
              "settings": settings, "notes": notes}
    result["summary"] = summarize_probe(result)
    return result


def summarize_probe(p: dict) -> list:
    lines = []
    s = p.get("slurm") or {}
    if not s.get("available"):
        lines.append(f"Slurm: {s.get('reason') or 'not available'}")
    else:
        parts = s.get("gpu_partitions") or []
        if parts:
            desc = ", ".join(
                r["name"] + (f" ({r['gpus']} {r['gpu_type'] or 'GPU'})" if r.get("gpus") else "")
                + ("" if r.get("cuda_capable", True) else " [not CUDA]")
                for r in parts)
            lines.append(f"Slurm GPU partitions: {desc}")
        else:
            lines.append(f"Slurm: {s.get('reason') or 'no GPU partitions'}")
    c = p.get("cuda") or {}
    if c.get("driver"):
        lines.append(f"NVIDIA driver {c['driver']} (CUDA {c.get('driver_cuda') or '?'}) "
                     f"on this node: {', '.join(c.get('gpus') or []) or 'GPU names unknown'}")
    else:
        lines.append("No NVIDIA GPU visible on this node"
                     + (" (nvidia-smi missing)" if not c.get("nvidia_smi") else ""))
    if c.get("modules"):
        lines.append("CUDA modules: " + ", ".join(c["modules"][:6]))
    o = p.get("orca") or {}
    if o.get("ready"):
        found = [i["label"] for i in o.get("installations", []) if i.get("is_6_1_or_newer")]
        found += o.get("modules_6_1_or_newer", [])
        lines.append("ORCA >= 6.1: " + ", ".join(found))
    else:
        seen = [i["label"] for i in o.get("installations", [])] + list(o.get("modules", []))
        lines.append("ORCA 6.1 not found" + (f" (found: {', '.join(seen)})" if seen else ""))
    k = p.get("conda") or {}
    lines.append(f"Environment tool: {k['tool']} ({len(k.get('envs') or [])} envs)"
                 if k.get("tool") else "No micromamba/mamba/conda found")
    envs = p.get("mlip_envs") or []
    for e in envs:
        pk = ", ".join(f"{n} {v}" for n, v in sorted(e["packages"].items()))
        lines.append(f"MLIP env {e['name']}: {pk or 'no MLIP packages'}; "
                     f"ProgExt: {', '.join(e['progext']) or 'none'}")
    if not envs:
        lines.append("No MLIP environment found")
    st = p.get("settings") or {}
    eff = st.get("effective") or {}
    lines.append("Settings: gpu_partitions="
                 + (eff.get("gpu_partitions") or "(empty)")
                 + f", gres={eff.get('gres', '')}, env_activation="
                 + ("set" if eff.get("env_activation") else "(empty)")
                 + (f"; env overrides: {', '.join(st['env_overrides'])}"
                    if st.get("env_overrides") else ""))
    return lines


# ----------------------------------------------------------------- proposal


def cuda_build_for(driver_cuda: str = "", gpu_names: Sequence[str] = ()) -> tuple:
    """(tag, reason) of the torch CUDA build to install, e.g. ('cu126', ...)."""
    drv = _version_tuple(driver_cuda)
    minimum = (0, 0)
    for rx, need in _MIN_BUILD_FOR_GPU:
        if any(rx.search(g or "") for g in gpu_names):
            minimum = max(minimum, need)
    if drv is None:
        build = max(_DEFAULT_CUDA_BUILD, minimum)
        reason = ("no NVIDIA driver visible on this node; CUDA "
                  f"{build[0]}.{build[1]} build assumed (runs on any CUDA "
                  f"{build[0]} driver); the verify job reports the GPU node's driver")
    else:
        fits = [b for b in _TORCH_CUDA_BUILDS if b <= drv]
        build = fits[0] if fits else _TORCH_CUDA_BUILDS[-1]
        if minimum > build:
            build = minimum
        reason = (f"driver supports CUDA {drv[0]}.{drv[1]}; newest torch build "
                  f"not above it: CUDA {build[0]}.{build[1]}")
    return f"cu{build[0]}{build[1]}", reason


def licence_notes(backend: str) -> dict:
    """Code and weight licence records of an orca-mlips backend."""
    from delfin.mlp_tools import registry
    from delfin.mlp_tools.registry import _ORCA_MLIPS_CODE

    b = backend.strip().lower()
    if b not in REGISTRY_BACKEND:
        raise ValueError(f"backend {backend!r} is not an orca-mlips plugin: "
                         f"use one of {', '.join(ORCA_MLIPS_BACKENDS)}")
    spec = registry.get_backend(REGISTRY_BACKEND[b])
    weight = dict(spec.weight_licence)
    lines = [f"orca-mlips code: {_ORCA_MLIPS_CODE['name']} ({_ORCA_MLIPS_CODE['url']})",
             f"{spec.name} code: {spec.code_licence['name']} ({spec.code_licence['url']})",
             f"weights: {weight['name']} ({weight['url']})"
             + ("" if weight.get("verified") else " -- record not verified against "
                "the upstream text; check the model card")]
    if weight.get("gated"):
        lines.append("gated checkpoint: accept the model terms on its Hugging Face "
                     "page and log in yourself with `hf auth login`; DELFIN never "
                     "handles the token. Re-probe afterwards.")
    if weight.get("non_commercial"):
        lines.append("non-commercial weights: academic / non-commercial use only")
    return {"backend": b, "registry_backend": spec.name,
            "package_code": dict(_ORCA_MLIPS_CODE),
            "code": dict(spec.code_licence), "weights": weight,
            "gated": bool(weight.get("gated")),
            "non_commercial": bool(weight.get("non_commercial")),
            "lines": lines}


def default_backend() -> tuple:
    """The first orca-mlips backend whose weights are neither gated nor
    non-commercial and whose weight record is verified."""
    for b in ORCA_MLIPS_BACKENDS:
        lic = licence_notes(b)
        if not lic["gated"] and not lic["non_commercial"] and lic["weights"].get("verified"):
            return b, ("first orca-mlips backend with verified, ungated, "
                       "commercially usable weights")
    return ORCA_MLIPS_BACKENDS[0], "no ungated backend recorded"


def _choose_partition(p: dict, requested: str) -> tuple:
    slurm = p.get("slurm") or {}
    rows = slurm.get("gpu_partitions") or []
    names = [r["name"] for r in rows]
    if requested:
        if not _PARTITION_RE.match(requested):
            raise ValueError(f"partition {requested!r} is not a partition name list")
        unknown = [n for n in requested.split(",") if names and n not in names]
        why = "requested" + (f"; not among the probed GPU partitions: {', '.join(unknown)}"
                             if unknown else "")
        return requested, why
    configured = ((p.get("settings") or {}).get("effective") or {}).get("gpu_partitions", "")
    if configured:
        return configured, "already configured"
    usable = [r["name"] for r in rows if r.get("cuda_capable", True)]
    if usable:
        return usable[0], "first CUDA-capable GPU partition reported by Slurm"
    return "", slurm.get("reason") or "no GPU partition found"


@dataclass
class Command:
    """One step of a plan. ``kind`` 'run' is executable after approval;
    'manual' is for the user to do themselves (never run by DELFIN)."""

    argv: list
    purpose: str
    kind: str = "run"

    @property
    def display(self) -> str:
        return " ".join(self.argv) if self.kind == "manual" else shlex.join(self.argv)

    def to_dict(self) -> dict:
        return {"command": self.display, "argv": list(self.argv),
                "purpose": self.purpose, "kind": self.kind}


def _tool_cmd(conda: dict) -> tuple:
    tool = conda.get("tool") or ""
    path = conda.get("path") or ""
    # The bare name when the probe found it on PATH (what the user types);
    # the absolute path otherwise, so the command runs as shown.
    exe = tool if (conda.get("on_path") or not path) else path
    return tool, exe


def predicted_prefix(conda: dict, env_name: str, env: Optional[dict] = None) -> str:
    env = env if env is not None else os.environ
    for e in conda.get("envs") or []:
        pe = Path(e)
        if pe.name == env_name:
            return str(pe)
    for e in conda.get("envs") or []:
        pe = Path(e)
        if pe.parent.name == "envs":
            return str(pe.parent / env_name)
    root = env.get("MAMBA_ROOT_PREFIX") or env.get("CONDA_PREFIX_1") or ""
    return str(Path(root) / "envs" / env_name) if root else ""


def install_commands(conda: dict, *, backend: str, cuda_tag: str,
                     env_name: str = DEFAULT_ENV_NAME) -> list:
    """The install steps for one backend. Empty when no conda-family tool."""
    if not _ENV_NAME_RE.match(env_name):
        raise ValueError(f"environment name {env_name!r} is not a plain name")
    if not re.fullmatch(r"cu\d{3,4}", cuda_tag):
        raise ValueError(f"CUDA build tag {cuda_tag!r} is not of the form cuXYZ")
    b = backend.strip().lower()
    if b not in ORCA_MLIPS_BACKENDS:
        raise ValueError(f"backend {backend!r} is not an orca-mlips plugin")
    tool, exe = _tool_cmd(conda)
    if not tool:
        return []
    create = [exe, "create", "-y", "-n", env_name]
    if tool == "micromamba":
        create += ["-c", "conda-forge"]
    create.append(f"python={PYTHON_VERSION}")
    run = [exe, "run", "-n", env_name, "python", "-m", "pip", "install"]
    cmds = [
        Command(create, f"create the environment {env_name} with Python {PYTHON_VERSION}"),
        Command(run + ["torch", "--index-url", TORCH_INDEX_URL.format(tag=cuda_tag)],
                f"PyTorch built for CUDA ({cuda_tag}) from the official PyTorch index; "
                + TORCH_DOWNLOAD_ESTIMATE),
        Command(run + [f"orca-mlips[{b}]"],
                f"orca-mlips with the {b} plugin (installs the orca-mlips-{b} ProgExt command)"),
    ]
    if licence_notes(b)["gated"]:
        cmds.append(Command(["hf", "auth", "login"],
                            "run this yourself in your terminal (it asks for your token; "
                            "DELFIN does not see it), after accepting the model terms "
                            f"at {licence_notes(b)['weights']['url']}; then re-probe",
                            kind="manual"))
    return cmds


def remove_command(conda: dict, env_name: str = DEFAULT_ENV_NAME) -> Optional[Command]:
    if not _ENV_NAME_RE.match(env_name):
        raise ValueError(f"environment name {env_name!r} is not a plain name")
    tool, exe = _tool_cmd(conda)
    if not tool:
        return None
    return Command([exe, "env", "remove", "-y", "-n", env_name],
                   f"delete the environment {env_name} and everything installed in it")


def settings_values(*, partition: str, prefix: str, backend: str,
                    orca_module: str = "", gres: str = "gpu:1",
                    cpus_per_task: int = 2, time_limit: str = "04:00:00",
                    mem_mb: int = 16000) -> dict:
    """The settings the plan would write (dotted key -> value)."""
    activation = []
    if orca_module:
        activation.append(f"module load {orca_module}")
    if prefix:
        if not _SAFE_PATH_RE.match(prefix):
            raise ValueError(f"environment prefix {prefix!r} holds shell characters")
        activation.append(f'export PATH="{prefix}/bin:$PATH"')
    return {
        "runtime.slurm.gpu_partitions": partition,
        "runtime.mlp.gres": gres,
        "runtime.mlp.env_activation": "; ".join(activation),
        "runtime.mlp.cpus_per_task": int(cpus_per_task),
        "runtime.mlp.time_limit": time_limit,
        "runtime.mlp.mem_mb": int(mem_mb),
        "runtime.mlp.backend": backend,
        "runtime.mlp.env_prefix": prefix,
    }


def propose(p: dict, *, backend: str = "", partition: str = "",
            env_name: str = DEFAULT_ENV_NAME, cuda: str = "") -> dict:
    """A set-up plan from a probe result. Runs nothing, writes nothing."""
    blockers: list = []
    if backend:
        b = backend.strip().lower()
        if b not in ORCA_MLIPS_BACKENDS:
            raise ValueError(f"backend {backend!r} is not an orca-mlips plugin: "
                             f"use one of {', '.join(ORCA_MLIPS_BACKENDS)}")
        why_b = "requested"
    else:
        b, why_b = default_backend()
    part, why_p = _choose_partition(p, partition)
    if not part:
        blockers.append(f"GPU partition: {why_p}")
    c = p.get("cuda") or {}
    if cuda:
        tag = cuda if cuda.startswith("cu") else "cu" + cuda.replace(".", "")
        why_c = "requested"
    else:
        gpu_names = list(c.get("gpus") or [])
        for r in (p.get("slurm") or {}).get("gpu_partitions") or []:
            if r["name"] in part.split(","):
                gpu_names += [r.get("gpu_type", "")] + list(r.get("gres") or [])
        tag, why_c = cuda_build_for(c.get("driver_cuda", ""), gpu_names)
    conda = p.get("conda") or {}
    cmds = install_commands(conda, backend=b, cuda_tag=tag, env_name=env_name)
    if not conda.get("tool"):
        blockers.append("no micromamba/mamba/conda on this host: install micromamba "
                        "yourself (https://mamba.readthedocs.io/en/latest/installation/"
                        "micromamba-installation.html), then re-probe")
    orca = p.get("orca") or {}
    if not orca.get("ready"):
        blockers.append("ORCA 6.1 or newer not found (installation or module); "
                        "ExtOpt via orca-mlips needs it on the GPU node")
    existing = next((e for e in p.get("mlip_envs") or []
                     if f"orca-mlips-{b}" in e.get("progext", [])), None)
    prefix = existing["prefix"] if existing else predicted_prefix(conda, env_name)
    orca_module = ""
    if not any(i.get("is_6_1_or_newer") for i in orca.get("installations", [])):
        mods = orca.get("modules_6_1_or_newer") or []
        orca_module = mods[0] if mods else ""
    values = settings_values(partition=part, prefix=prefix, backend=b,
                             orca_module=orca_module)
    lic = licence_notes(b)
    rm = remove_command(conda, env_name)
    plan = {
        "backend": b, "backend_reason": why_b,
        "partition": part, "partition_reason": why_p,
        "cuda_build": tag, "cuda_reason": why_c,
        "env_name": env_name, "env_prefix": prefix,
        "env_exists": bool(existing),
        "licence": lic,
        "commands": [x.to_dict() for x in cmds],
        "remove_command": rm.to_dict() if rm else None,
        "settings": values,
        "download_estimate": {"torch": TORCH_DOWNLOAD_ESTIMATE,
                              "weights": "downloaded on first use into your own cache; "
                                         "size not recorded by DELFIN"},
        "blockers": blockers,
        "note": ("Commands are shown, never run silently: each install step is "
                 "asked before it runs. Run them one per shell call, or with "
                 "`delfin mlp_setup install` in your terminal (asks for each)."),
    }
    if existing:
        plan["note"] += f" An environment with orca-mlips-{b} already exists ({prefix}); install can be skipped."
    return plan


# ------------------------------------------------------------------ install


def run_install(commands: Sequence[dict], *, approve: Callable[[str, str], Optional[bool]],
                runner: Optional[Callable[[list], int]] = None) -> dict:
    """Run each 'run' step only after ``approve(command, purpose)`` is True.

    ``approve`` returning False stops the sequence (later steps depend on
    earlier ones); returning None means nobody can be asked and also stops.
    'manual' steps are never run, only reported. ``runner(argv) -> returncode``
    defaults to a subprocess with the terminal attached.
    """
    if approve is None:
        raise ValueError("run_install needs an approve callback; there is no default yes")
    runner = runner or (lambda argv: subprocess.run(list(argv)).returncode)
    done: list = []
    for c in commands:
        if c.get("kind") == "manual":
            done.append({"command": c["command"], "status": "manual",
                         "note": "for you to run yourself; DELFIN does not run it"})
            continue
        answer = approve(c["command"], c.get("purpose", ""))
        if answer is None:
            done.append({"command": c["command"], "status": "not asked",
                         "note": "nobody can be asked here; nothing ran"})
            return {"ok": False, "steps": done, "stopped": "no approval channel"}
        if not answer:
            done.append({"command": c["command"], "status": "declined"})
            return {"ok": False, "steps": done, "stopped": "declined"}
        rc = runner(list(c["argv"]))
        done.append({"command": c["command"], "status": "ran", "returncode": rc})
        if rc != 0:
            return {"ok": False, "steps": done, "stopped": f"exit code {rc}"}
    return {"ok": True, "steps": done, "stopped": ""}


# ---------------------------------------------------------------- configure


def _validate_values(values: dict) -> None:
    v = values
    if v.get("runtime.slurm.gpu_partitions") and not _PARTITION_RE.match(
            v["runtime.slurm.gpu_partitions"]):
        raise ValueError("gpu_partitions must be comma-separated partition names")
    if v.get("runtime.mlp.gres") and not _GRES_RE.match(v["runtime.mlp.gres"]):
        raise ValueError("gres must look like gpu:1 or gpu:a100:1")
    if v.get("runtime.mlp.time_limit") and not _TIME_RE.match(v["runtime.mlp.time_limit"]):
        raise ValueError("time_limit must be [D-]HH[:MM[:SS]]")
    for k in ("runtime.mlp.cpus_per_task", "runtime.mlp.mem_mb"):
        if k in v and not (0 < int(v[k]) <= 10_000_000):
            raise ValueError(f"{k} must be a positive integer")
    if "\n" in str(v.get("runtime.mlp.env_activation", "")):
        raise ValueError("env_activation must be one line (steps separated by ';')")


def _get(d: dict, dotted: str):
    cur = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _set(d: dict, dotted: str, value) -> None:
    parts = dotted.split(".")
    cur = d
    for part in parts[:-1]:
        nxt = cur.get(part)
        if not isinstance(nxt, dict):
            nxt = {}
            cur[part] = nxt
        cur = nxt
    cur[parts[-1]] = value


def configure(values: dict, *, write: bool = False, settings_path=None) -> dict:
    """The settings diff for ``values``; written only with ``write=True``.

    The write goes through :func:`delfin.user_settings.save_settings` (never
    a direct file write) and is read back with ``load_settings``; a key that
    did not survive the round trip is reported in ``not_kept``.
    """
    from delfin.user_settings import get_settings_path, load_settings, save_settings

    _validate_values(values)
    current = load_settings(settings_path)
    new = copy.deepcopy(current)
    diff = []
    for key, value in values.items():
        old = _get(current, key)
        if old != value:
            diff.append({"key": key, "old": old, "new": value})
            _set(new, key, value)
    out = {"diff": diff, "written": False,
           "settings_path": str(get_settings_path(settings_path)), "not_kept": []}
    if write and diff:
        save_settings(new, settings_path)
        back = load_settings(settings_path)
        out["not_kept"] = [d["key"] for d in diff if _get(back, d["key"]) != d["new"]]
        out["written"] = True
    return out


def format_diff(diff: Sequence[dict]) -> str:
    if not diff:
        return "settings: no change"
    return "\n".join(f"  {d['key']}: {d['old']!r} -> {d['new']!r}" for d in diff)


# ------------------------------------------------------------------- doctor


def mlip_readiness(settings: Optional[dict] = None,
                   exists: Callable[[str], bool] = os.path.exists) -> dict:
    """What the doctor reports: env, partition and backend configured and
    the environment's ProgExt command present. File checks only."""
    if settings is None:
        settings = _default_settings()
    runtime = (settings or {}).get("runtime") or {}
    part = str((runtime.get("slurm") or {}).get("gpu_partitions", "") or "")
    mlp = runtime.get("mlp") or {}
    backend = str(mlp.get("backend", "") or "")
    prefix = str(mlp.get("env_prefix", "") or "") or _env_prefix_from_activation(
        str(mlp.get("env_activation", "") or ""))
    if not (part or mlp):
        return {"configured": False, "ready": False, "missing": [],
                "detail": "MLIP: not set up (optional) -- delfin mlp_setup plan"}
    missing = []
    if not prefix:
        missing.append("env")
    elif backend and not exists(str(Path(prefix) / "bin" / f"orca-mlips-{backend}")):
        missing.append(f"env (no orca-mlips-{backend} in {prefix}/bin)")
    elif not backend and not exists(prefix):
        missing.append(f"env ({prefix} not found)")
    if not part:
        missing.append("partition")
    if not backend:
        missing.append("backend")
    ready = not missing
    detail = ("MLIP: env/partition/backend ready" if ready
              else "MLIP: missing " + ", ".join(missing))
    return {"configured": True, "ready": ready, "missing": missing, "detail": detail}


# ------------------------------------------------------------------- verify


VERIFY_SCRIPT = r'''"""orca-mlips verify job (written by delfin.mlp_tools.guided_setup).

Standard library only; runs in the MLIP environment on the GPU node.
1. a fresh interpreter imports torch (start + import cost, CUDA visible?);
2. the orca-mlips ProgExt command is called as ORCA calls it, once cold and
   then warm: cold = process start + model load + one gradient, warm = one
   gradient with whatever the plugin keeps resident;
3. a short ORCA ExtOpt Opt of water;
then writes verify_result.json.
"""
import argparse, json, os, shlex, statistics, subprocess, sys, time

XYZ = "3\nwater\nO 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\nH 0.0 -0.7572 -0.4692\n"
TORCH_PROBE = (
    "import json,time;t=time.monotonic();import torch;"
    "d={'torch':torch.__version__,'cuda_build':torch.version.cuda,"
    "'cuda_available':torch.cuda.is_available()};"
    "d['device']=torch.cuda.get_device_name(0) if d['cuda_available'] else '';"
    "d['import_s']=time.monotonic()-t;print(json.dumps(d))")


def tail(text, n=3000):
    return (text or "")[-n:]


def parse_engrad(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            s = line.split("#", 1)[0].strip()
            if s:
                vals.extend(s.split())
    n = int(vals[0])
    energy = float(vals[1])
    grad = [float(x) for x in vals[2:2 + 3 * n]]
    return n, energy, grad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--progext", required=True)
    ap.add_argument("--ext-params", default="")
    ap.add_argument("--orca", default="orca")
    ap.add_argument("--opt-input", default="water_opt.inp")
    ap.add_argument("--warm-calls", type=int, default=2)
    a = ap.parse_args()
    res = {"progext": a.progext, "ext_params": a.ext_params}
    t = time.monotonic()
    r = subprocess.run([sys.executable, "-c", TORCH_PROBE], capture_output=True, text=True)
    res["torch"] = {"wall_s": time.monotonic() - t, "rc": r.returncode,
                    "stderr_tail": tail(r.stderr)}
    try:
        res["torch"].update(json.loads(r.stdout.strip().splitlines()[-1]))
    except Exception:
        pass
    with open("probe_EXT.xyz", "w") as fh:
        fh.write(XYZ)
    with open("probe_EXT.extinp.tmp", "w") as fh:
        fh.write("probe_EXT.xyz\n0\n1\n1\n1\n")
    calls = []
    for i in range(1 + max(0, a.warm_calls)):
        if os.path.exists("probe_EXT.engrad"):
            os.remove("probe_EXT.engrad")
        t = time.monotonic()
        try:
            r = subprocess.run([a.progext, "probe_EXT.extinp.tmp"] + shlex.split(a.ext_params),
                               capture_output=True, text=True, timeout=1800)
            rc, err = r.returncode, r.stderr
        except Exception as exc:
            rc, err = -1, "%s: %s" % (type(exc).__name__, exc)
        call = {"wall_s": time.monotonic() - t, "rc": rc, "stderr_tail": tail(err)}
        try:
            n, e, g = parse_engrad("probe_EXT.engrad")
            call.update(n_atoms=n, energy_eh=e, n_gradient=len(g),
                        grad_norm=sum(x * x for x in g) ** 0.5)
        except Exception as exc:
            call["engrad_error"] = "%s: %s" % (type(exc).__name__, exc)
        calls.append(call)
        if rc != 0:
            break
    res["gradient_calls"] = calls
    warm = [c["wall_s"] for c in calls[1:] if c.get("rc") == 0]
    if calls and calls[0].get("rc") == 0 and warm:
        res["timing"] = {"cold_call_s": calls[0]["wall_s"],
                         "warm_call_s": statistics.median(warm),
                         "start_and_load_s": calls[0]["wall_s"] - statistics.median(warm)}
    t = time.monotonic()
    try:
        with open("water_opt.out", "w") as out:
            r = subprocess.run([a.orca, a.opt_input], stdout=out, stderr=subprocess.PIPE,
                               text=True, timeout=1500)
        rc, err = r.returncode, r.stderr
    except Exception as exc:
        rc, err = -1, "%s: %s" % (type(exc).__name__, exc)
    opt = {"wall_s": time.monotonic() - t, "rc": rc, "stderr_tail": tail(err)}
    try:
        with open("water_opt.out") as fh:
            text = fh.read()
    except Exception:
        text = ""
    opt["cycles"] = text.count("GEOMETRY OPTIMIZATION CYCLE")
    opt["converged"] = "THE OPTIMIZATION HAS CONVERGED" in text
    opt["terminated_normally"] = "ORCA TERMINATED NORMALLY" in text
    energies = [ln.split()[-1] for ln in text.splitlines() if "FINAL SINGLE POINT ENERGY" in ln]
    opt["final_energy_eh"] = float(energies[-1]) if energies else None
    opt["out_tail"] = tail(text, 2000)
    res["extopt_opt"] = opt
    try:
        subprocess.run([a.progext, "--stop-server"], capture_output=True, timeout=60)
    except Exception:
        pass
    with open("verify_result.json", "w") as fh:
        json.dump(res, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''


def _water_xyz_lists() -> tuple:
    return [a for a, _ in _WATER], [list(c) for _, c in _WATER]


def render_verify_sbatch(cfg, *, job_dir: Path, progext: str, ext_params: str,
                         orca: str = "orca") -> str:
    """The verify job's sbatch text (partition, gres, node and CPU from cfg)."""
    if not cfg.has_gpu_partition():
        raise ValueError("no GPU partition configured: set runtime.slurm.gpu_partitions "
                         "(delfin mlp_setup configure) before the verify job")
    lines = [
        "#!/bin/bash",
        "# delfin orca-mlips verify job -- generated by delfin.mlp_tools.guided_setup",
        f"#SBATCH --job-name={VERIFY_JOB_NAME}",
        f"#SBATCH --partition={cfg.gpu_partitions}",
        f"#SBATCH --gres={cfg.gres}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(cfg.cpus_per_task)}",
        f"#SBATCH --mem={int(cfg.mem_mb)}M",
        f"#SBATCH --time={VERIFY_TIME_LIMIT}",
        f"#SBATCH --chdir={job_dir}",
        f"#SBATCH --output={VERIFY_JOB_NAME}_%j.out",
        f"#SBATCH --error={VERIFY_JOB_NAME}_%j.err",
        "",
        "set -uo pipefail",
        "export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2",
        "export TORCH_NUM_THREADS=2",
        "export PYTHONNOUSERSITE=1",
    ]
    for chunk in cfg.env_activation.split(";"):
        if chunk.strip():
            lines.append(chunk.strip())
    lines += [
        "",
        f"cd {shlex.quote(str(job_dir))}",
        "if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi > nvidia_smi.txt 2>&1; "
        "else echo 'nvidia-smi not found on this node' > nvidia_smi.txt; fi",
        "python mlip_verify.py --progext " + shlex.quote(progext)
        + " --ext-params " + shlex.quote(ext_params)
        + " --orca " + shlex.quote(orca) + " --opt-input water_opt.inp",
        "",
    ]
    return "\n".join(lines)


def prepare_verify_job(cfg, *, backend: str, job_dir, orca: str = "orca") -> dict:
    """Dry run: write the verify job directory. Submits nothing.

    Writes ``mlip_verify.py``, the ExtOpt Opt input (through
    :func:`delfin.mlp_tools.orca_mlips.build_orca_mlips_input`, device
    cuda) and ``mlp_verify.sbatch``. Returns the paths and the submit
    command as text.
    """
    from delfin.mlp_tools.orca_mlips import build_orca_mlips_input, ext_params_for

    b = backend.strip().lower()
    if b not in ORCA_MLIPS_BACKENDS:
        raise ValueError(f"backend {backend!r} is not an orca-mlips plugin")
    if orca != "orca" and not _SAFE_PATH_RE.match(orca):
        raise ValueError(f"ORCA path {orca!r} holds shell characters")
    job_dir = Path(job_dir).expanduser().resolve()
    progext = f"orca-mlips-{b}"
    ext_params = ext_params_for(b, device="cuda")
    script_text = render_verify_sbatch(cfg, job_dir=job_dir, progext=progext,
                                       ext_params=ext_params, orca=orca)
    job_dir.mkdir(parents=True, exist_ok=True)
    elements, coords = _water_xyz_lists()
    inp = build_orca_mlips_input("opt", elements=elements, coordinates=coords,
                                 backend=b, device="cuda", work_dir=job_dir,
                                 base_name="water_opt")
    (job_dir / "mlip_verify.py").write_text(VERIFY_SCRIPT, encoding="utf-8")
    script = job_dir / f"{VERIFY_JOB_NAME}.sbatch"
    script.write_text(script_text, encoding="utf-8")
    script.chmod(0o755)
    return {"job_dir": str(job_dir), "script": str(script), "input": str(inp),
            "backend": b, "progext": progext, "ext_params": ext_params,
            "submit_command": shlex.join([_SBATCH_PROG, str(script)]),
            "submitted": False, "dry_run": True}


def submit_verify_job(prepared: dict, *, active_jobs: Optional[int] = None,
                      runner: Optional[Callable] = None) -> dict:
    """Submit a prepared verify job through scheduler_client (job cap held).

    The caller has the user's approval; this function does not ask.
    """
    from delfin.mlp_tools.cluster_dispatch import assert_can_submit, count_active_mlp_jobs

    active = active_jobs if active_jobs is not None else count_active_mlp_jobs()
    assert_can_submit(active=active)
    argv = [_SBATCH_PROG, str(prepared["script"])]
    seen: dict = {}

    def _run(cmd):
        rc, out, err = _default_run(cmd, timeout=60)
        seen.update(rc=rc, out=out, err=err)
        return out if rc == 0 else None

    answer = scheduler_client.query(scheduler_client.kind_of_argv(argv), argv,
                                    runner=runner or _run)
    out = dict(prepared, dry_run=False)
    if answer is scheduler_client.THROTTLED:
        out.update(submitted=False, error="sbatch throttled by the scheduler client; retry in 25 s")
        return out
    text = str(answer or "")
    m = re.search(r"Submitted batch job (\d+)", text)
    out.update(submitted=bool(m), job_id=m.group(1) if m else "",
               sbatch_output=text.strip() or str(seen.get("err", "")).strip())
    if not m:
        out["error"] = "sbatch did not submit: " + (out["sbatch_output"] or "no output")
        out["likely_causes"] = diagnose(out["sbatch_output"])
    return out


_CAUSES = (
    ("cuda", re.compile(r"CUDA driver version is insufficient|no kernel image is available|"
                        r"NVIDIA driver on your system is too old|Found no NVIDIA driver|"
                        r"CUDA error|\"cuda_available\":\s*false", re.I),
     "CUDA mismatch: the torch CUDA build does not run on this node's driver, or the job "
     "has no GPU; re-run the plan with the driver's CUDA version (--cuda) and reinstall torch"),
    ("login", re.compile(r"GatedRepoError|Cannot access gated repo|401 Client Error|"
                         r"is restricted|LocalEntryNotFoundError|huggingface_hub\.errors", re.I),
     "gated checkpoint not reachable: accept the model terms on Hugging Face and run "
     "`hf auth login` yourself, then re-probe"),
    ("partition", re.compile(r"invalid partition|Invalid generic resource|"
                             r"Requested node configuration is not available|"
                             r"Invalid account or account/partition", re.I),
     "partition or gres refused by Slurm: choose another GPU partition or gres"),
    ("module", re.compile(r"orca: command not found|module: command not found|"
                          r"Unable to locate a modulefile|No such file or directory: 'orca'", re.I),
     "ORCA not found in the job: add the ORCA module or path to runtime.mlp.env_activation"),
    ("env", re.compile(r"orca-mlips-\w+: command not found|No such file or directory: "
                       r"'orca-mlips|ModuleNotFoundError|No module named", re.I),
     "MLIP environment not active in the job: check runtime.mlp.env_activation and the env"),
)


def diagnose(text: str) -> list:
    """Likely causes named by patterns in a job's output."""
    return [{"cause": key, "fix": fix} for key, rx, fix in _CAUSES if rx.search(text or "")]


def _finite(x) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def read_verify_result(job_dir) -> dict:
    """Pass/fail of a verify job with its numbers and, on failure, the
    likely causes. Reads files only."""
    job_dir = Path(job_dir)
    logs = ""
    for pat in (f"{VERIFY_JOB_NAME}_*.err", f"{VERIFY_JOB_NAME}_*.out", "nvidia_smi.txt"):
        for f in sorted(job_dir.glob(pat)):
            try:
                logs += f.read_text(encoding="utf-8", errors="replace")[-4000:] + "\n"
            except OSError:
                pass
    path = job_dir / VERIFY_RESULT
    if not path.is_file():
        state = "failed before writing a result" if logs.strip() else "no result yet (queued or running)"
        return {"status": "pending" if not logs.strip() else "fail", "detail": state,
                "likely_causes": diagnose(logs), "numbers": {}}
    try:
        res = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return {"status": "fail", "detail": f"unreadable {VERIFY_RESULT}: {exc}",
                "likely_causes": diagnose(logs), "numbers": {}}
    torch = res.get("torch") or {}
    calls = res.get("gradient_calls") or []
    opt = res.get("extopt_opt") or {}
    problems = []
    if not torch.get("cuda_available"):
        problems.append("torch sees no CUDA device")
    grads_ok = bool(calls) and all(
        c.get("rc") == 0 and _finite(c.get("energy_eh")) and c.get("n_gradient") == 9
        for c in calls)
    if not grads_ok:
        problems.append("energy+gradient call failed or returned no 3x3 gradient")
    if not (opt.get("terminated_normally") and opt.get("converged")):
        problems.append("ExtOpt Opt did not converge / terminate normally")
    timing = res.get("timing") or {}
    numbers = {
        "torch_import_s": torch.get("import_s"),
        "interpreter_start_and_torch_import_s": torch.get("wall_s"),
        "device": torch.get("device", ""),
        "torch": torch.get("torch", ""), "torch_cuda_build": torch.get("cuda_build", ""),
        "cold_call_s": timing.get("cold_call_s"),
        "warm_call_s": timing.get("warm_call_s"),
        "start_and_model_load_s": timing.get("start_and_load_s"),
        "energy_eh": calls[0].get("energy_eh") if calls else None,
        "opt_wall_s": opt.get("wall_s"), "opt_cycles": opt.get("cycles"),
        "opt_s_per_cycle": (opt["wall_s"] / opt["cycles"]
                            if opt.get("cycles") and _finite(opt.get("wall_s")) else None),
        "opt_final_energy_eh": opt.get("final_energy_eh"),
    }
    text = logs + json.dumps(res)
    return {"status": "pass" if not problems else "fail",
            "detail": "all checks passed" if not problems else "; ".join(problems),
            "numbers": numbers,
            "likely_causes": [] if not problems else diagnose(text)}


def verify_summary(result: dict) -> str:
    n = result.get("numbers") or {}

    def f(x, unit="s"):
        return f"{x:.2f} {unit}" if _finite(x) else "n/a"

    lines = [f"verify: {result.get('status', '?').upper()} -- {result.get('detail', '')}"]
    if n:
        lines.append(f"  device {n.get('device') or 'none'}, torch {n.get('torch') or '?'} "
                     f"(CUDA {n.get('torch_cuda_build') or '?'})")
        lines.append(f"  interpreter start + torch import {f(n.get('interpreter_start_and_torch_import_s'))}")
        lines.append(f"  first ProgExt call {f(n.get('cold_call_s'))}, warm call "
                     f"{f(n.get('warm_call_s'))} -> start + model load ~{f(n.get('start_and_model_load_s'))}")
        lines.append(f"  ExtOpt Opt {f(n.get('opt_wall_s'))} over {n.get('opt_cycles') or 0} cycles "
                     f"({f(n.get('opt_s_per_cycle'))}/cycle)")
    for c in result.get("likely_causes") or []:
        lines.append(f"  likely cause [{c['cause']}]: {c['fix']}")
    return "\n".join(lines)


def load_mlp_backend(settings: Optional[dict] = None) -> str:
    """The configured orca-mlips backend (``runtime.mlp.backend``) or ''."""
    settings = settings if settings is not None else _default_settings()
    b = str((((settings or {}).get("runtime") or {}).get("mlp") or {}).get("backend", "") or "")
    return b if b in ORCA_MLIPS_BACKENDS else ""


def cluster_config(settings: Optional[dict] = None, env: Optional[dict] = None):
    """MLPClusterConfig from the user's settings and the environment."""
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig

    settings = settings if settings is not None else _default_settings()
    return MLPClusterConfig.from_sources(settings=(settings or {}).get("runtime") or {},
                                         env=env)


__all__ = [
    "Host", "probe", "summarize_probe", "propose", "cuda_build_for",
    "licence_notes", "default_backend", "install_commands", "remove_command",
    "run_install", "configure", "format_diff", "settings_values",
    "mlip_readiness", "prepare_verify_job", "submit_verify_job",
    "read_verify_result", "verify_summary", "diagnose", "inspect_env",
    "cluster_config",
]
