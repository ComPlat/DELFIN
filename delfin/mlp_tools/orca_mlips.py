"""orca-mlips as a second ExtOpt provider next to ``delfin-mlp-extopt``.

orca-mlips (github.com/t-0hmura/orca-mlips, MIT) ships four MLIP plugins
for ORCA's ExtTool/ProgExt interface with console entry points named
``uma``, ``orb``, ``mace``, ``aimnet2`` (prefixed aliases
``orca-mlips-uma`` ...).  This module

* detects the installed package and lists its real ProgExt entry points
  (never invents them, never installs anything — install hint only);
* builds the ``! ExtOpt`` ORCA input (opt, optts, freq, GOAT, NEB,
  NEB-TS, IRC, relaxed scan) through the same OPI machinery as
  :mod:`delfin.mlp_tools.orca_extopt` (parity-tested in
  tests/test_opi_mlp_extopt_builder.py) with the ProgExt wrapper being the
  orca-mlips console command;
* supports the xTB implicit-solvent correction (``--solvent``,
  ``--solvent-model alpb|cpcmx``, OPTIONS.md) and the two-step analytical
  Hessian for OptTS/Freq: the MLIP writes ``.hess`` via ``--dump-hessian``
  and ORCA reads it back with ``%geom inhess Read / InHessName``
  (ORCA manual optimizations_TS.txt lines 394-400).

``delfin-mlp-extopt`` (orca_extopt.py) stays the default provider; this
module is used when the user asks for an orca-mlips backend explicitly.
"""

from __future__ import annotations

import importlib.metadata
import re
import shutil
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.mlp_tools.orca_extopt import (
    EXTOPT_RUN_TYPES,
    _build_with_opi,
)

#: The ProgExt name DELFIN passes for an orca-mlips backend: the prefixed
#: console-script alias.  orca-mlips README troubleshooting: 'ProgExt "uma"
#: runs the wrong plugin -- Use ProgExt "orca-mlips-uma" to avoid alias
#: conflicts.'
_MLIPS_BACKENDS = ("uma", "orb", "mace", "aimnet2")
_SOLVENT_MODELS = ("alpb", "cpcmx")
_BASE_NAME = "extopt"
#: ORCA reads the dumped Hessian back from this file name (relative to the
#: job directory, which is the ORCA working directory at run time).
_HESSIAN_FILE = f"{_BASE_NAME}.hess"


#: orca-mlips plugin -> the registry backend whose weight licence it loads
#: (delfin.mlp_tools.registry, LICENSE_GUARD.md family 4): UMA weights are
#: gated (FAIR Chemistry terms), orca-mlips' MACE plugin loads MACE
#: foundation weights (Academic Software License, non-commercial), ORB
#: Apache-2.0, AIMNet2 MIT.
MLIPS_REGISTRY_BACKEND = {
    "uma": "fairchem_uma",
    "mace": "mace_omol",
    "orb": "orb",
    "aimnet2": "aimnet2",
}


def mlips_weight_licence(backend: str) -> dict:
    """The weight-licence record of an orca-mlips plugin (registry record)."""
    from delfin.mlp_tools import registry

    name = backend.strip().lower()
    if name not in MLIPS_REGISTRY_BACKEND:
        raise ValueError(
            f"backend {backend!r} is not an orca-mlips plugin: use one of "
            f"{', '.join(_MLIPS_BACKENDS)}"
        )
    return dict(registry.get_backend(MLIPS_REGISTRY_BACKEND[name]).weight_licence)


def mlips_progext_names(backends: Sequence[str]) -> tuple[str, ...]:
    """Prefixed ProgExt names for the given orca-mlips backends.

    ``uma`` -> ``orca-mlips-uma`` (README alias-conflict troubleshooting).
    """
    out = []
    for b in backends:
        name = b.strip().lower()
        if name not in _MLIPS_BACKENDS:
            raise ValueError(
                f"backend {b!r} is not an orca-mlips plugin: use one of "
                f"{', '.join(_MLIPS_BACKENDS)}"
            )
        out.append(f"orca-mlips-{name}")
    return tuple(out)


#: The distribution and the package it installs: ``orca-mlips`` ships its
#: plugins as ``orca_mlips_plugins`` (pyproject ``[tool.setuptools]
#: packages``); there is no ``orca_mlips`` module.
ORCA_MLIPS_DIST = "orca-mlips"
ORCA_MLIPS_MODULE = "orca_mlips_plugins"
#: What installs it: DELFIN's installer, into an environment of its own
#: (delfin.common.external_builders.PROGRAMS["orca_mlips"]).
ORCA_MLIPS_INSTALL = "python -m delfin.installer --install orca-mlips"


def _env_dists(python: str):
    """The ``orca-mlips`` distributions in the site-packages beside *python*.

    Metadata only, read from the files; that interpreter is never started.
    """
    prefix = Path(python).resolve().parent.parent
    sites = [str(site) for site in sorted(prefix.glob("lib/python*/site-packages"))]
    if not sites:
        return []
    return list(importlib.metadata.distributions(name=ORCA_MLIPS_DIST, path=sites))


def _own_env_python(environ: Optional[Mapping[str, str]] = None) -> Optional[str]:
    """The interpreter of orca-mlips' own environment, or None when there is none.

    ``DELFIN_ORCA_MLIPS_PYTHON``, else the environment DELFIN's installer
    builds (``<AI tools root>/.mamba_env/orca_mlips``) when it exists.
    """
    import os
    import sys

    from delfin.common import external_builders

    python = external_builders.tool_python("orca_mlips", environ)
    if python == sys.executable or not os.access(python, os.X_OK):
        return None
    return python


def detect_orca_mlips(environ: Optional[Mapping[str, str]] = None) -> dict:
    """Detect the installed orca-mlips package and its ProgExt entry points.

    Looks first in orca-mlips' own environment (the one DELFIN's installer
    builds, or ``DELFIN_ORCA_MLIPS_PYTHON``), then in this interpreter.
    Reads distribution metadata only; never imports the MLIP stack (torch
    etc. must not load into this process), never installs anything.

    ``progext_paths`` maps each plugin to the absolute path of its
    ``orca-mlips-<backend>`` command where that file exists: the name ORCA's
    ``ProgExt`` needs on a compute node whose ``PATH`` does not hold the
    environment.
    """
    import sys
    import sysconfig

    dist = None
    bin_dir: Optional[Path] = None
    python = _own_env_python(environ)
    if python:
        found = _env_dists(python)
        if found:
            dist, bin_dir = found[0], Path(python).parent
    if dist is None:
        try:
            dist = importlib.metadata.distribution(ORCA_MLIPS_DIST)
            scripts = sysconfig.get_path("scripts")
            bin_dir = Path(scripts) if scripts else None
            python = sys.executable
        except importlib.metadata.PackageNotFoundError:
            dist = None
    names: list[str] = []
    paths: dict[str, str] = {}
    version = ""
    if dist is not None:
        version = dist.version or ""
        names = sorted(
            ep.name for ep in dist.entry_points
            if ep.group == "console_scripts" and ep.name.startswith("orca-mlips-")
        )
        for name in names:
            exe = bin_dir / name if bin_dir is not None else None
            if exe is not None and exe.is_file():
                paths[name[len("orca-mlips-"):]] = str(exe)
    installed = dist is not None
    return {
        "available": installed,
        "progext_names": names or (
            [f"orca-mlips-{b}" for b in _MLIPS_BACKENDS] if installed else []
        ),
        "progext_paths": paths,
        "python": python if installed else "",
        "version": version,
        "install_hint": (
            f"{ORCA_MLIPS_INSTALL}   (its own environment with a CUDA PyTorch "
            "and the aimnet2 + orb plugins; UMA and MACE need their own "
            "environment and licence terms)"
        ),
        "xtb_required_for_solvent": bool(shutil.which("xtb")) is False,
    }


def resolve_progext(backend: str, environ: Optional[Mapping[str, str]] = None) -> str:
    """What ExtOpt's ``ProgExt`` names for an orca-mlips plugin.

    The absolute ``orca-mlips-<backend>`` of the installed environment when
    there is one, else the command found on this ``PATH`` (absolute), else
    the bare name -- which ORCA then has to find on the ``PATH`` of the node
    it runs on.  ORCA starts ``ProgExt`` from its own process, and a bare
    name was "not found" on a GPU node whose ``PATH`` did not hold the
    environment (docs/MLIP_GOAT_GPU.md, D3).
    """
    (name,) = mlips_progext_names([backend])
    found = detect_orca_mlips(environ).get("progext_paths", {})
    path = found.get(name[len("orca-mlips-"):])
    if path:
        return path
    on_path = shutil.which(name)
    return str(Path(on_path).absolute()) if on_path else name


def ext_params_for(
    backend: str,
    *,
    model: Optional[str] = None,
    device: str = "cpu",
    solvent: Optional[str] = None,
    solvent_model: str = "alpb",
    precision: Optional[str] = None,
    extra: Sequence[str] = (),
) -> str:
    """The ``Ext_Params`` argument string forwarded to the orca-mlips command.

    Options per OPTIONS.md (github.com/t-0hmura/orca-mlips): --model,
    --device auto|cpu|cuda, --solvent <name|none>, --solvent-model
    alpb|cpcmx, --dump-hessian <path>, --precision (ORB) / --dtype (MACE).
    """
    name = backend.strip().lower()
    if name not in _MLIPS_BACKENDS:
        raise ValueError(
            f"backend {backend!r} is not an orca-mlips plugin: use one of "
            f"{', '.join(_MLIPS_BACKENDS)}"
        )
    if solvent_model not in _SOLVENT_MODELS:
        raise ValueError(
            f"solvent-model {solvent_model!r} not one of alpb|cpcmx "
            "(orca-mlips OPTIONS.md)"
        )
    if precision and name != "orb":
        raise ValueError(
            f"--precision is an option of the orb plugin only, not of {name!r} "
            "(orca-mlips cli_orca.py)"
        )
    # No ``--backend``: the command itself is the backend (orca-mlips-<name>),
    # and no plugin accepts that option -- argparse stops the plugin with
    # "unrecognized arguments" and ORCA reports "external program failed"
    # (docs/MLIP_GOAT_GPU.md, D4).
    parts: list[str] = []
    if model:
        parts.append(f"--model {model}")
    parts.append(f"--device {device}")
    if solvent is not None:
        parts.append(f"--solvent {solvent}")
        parts.append(f"--solvent-model {solvent_model}")
    if precision:
        parts.append(f"--precision {precision}")
    parts.extend(extra)
    return " ".join(parts)


#: Run types the orca-mlips builder writes, with their simple keywords.
#: Manual 4.1.10 ("Optimizing with External Methods"): ``!ExtOpt`` "must
#: be combined with other keywords like !Opt or !GOAT to select the
#: specific run type"; it names optimisations, NEB-TS and GOAT and adds
#: "Other calculations which exclusively rely on energies and gradients
#: are also supported" -- NEB, IRC (initial Hessian computed numerically
#: from gradients, or read from the MLIP's dumped Hessian) and the
#: relaxed surface scan.  ``freq`` reads the Hessian the MLIP dumps.
MLIPS_RUN_TYPES: dict[str, tuple[str, ...]] = {
    "opt": tuple(EXTOPT_RUN_TYPES["opt"]),
    "optts": tuple(EXTOPT_RUN_TYPES["optts"]),
    "freq": ("Freq",),
    "goat": ("GOAT",),
    "neb": ("NEB",),
    "neb-ts": ("NEB-TS",),
    "irc": ("IRC",),
    "scan": ("Opt",),
}

#: GOAT simple keywords, manual section 4.10.11 Table 4.8.  GOAT-COARSE is
#: left out: it needs user-defined fragments (section 4.10.9), which this
#: builder does not write.
GOAT_VARIANTS = ("GOAT", "GOAT-ENTROPY", "GOAT-EXPLORE", "GOAT-REACT", "GOAT-DIVERSITY")

#: Readable names -> ``%goat`` option, and the manual's spelling -> the
#: field name of OPI's ``BlockGoat``.  ORCA 6.1.1's parser accepts both
#: KEEPWORKERDATA and KEEPDATA, SKIPINITIALOPT and SKIPINITOPT (measured
#: with delfin.orca_knowledge.probe.probe_input on ``! XTB GOAT``).
GOAT_OPTION_ALIASES: dict[str, str] = {
    "temperature": "conftemp",       # Table 4.9 CONFTEMP, K, default 298.15
    "energy_window": "maxen",        # Table 4.9 MAXEN, kcal/mol
    "workers": "nworkers",           # Table 4.9 NWORKERS, default auto
    "keepworkerdata": "keepdata",    # Table 4.9 KEEPWORKERDATA
    "skipinitialopt": "skipinitopt",  # Table 4.9 SKIPINITIALOPT
}

#: GOAT's default number of workers: "the number of necessary GOAT
#: iterations will be computed and divided among NWorkers (8 by default)"
#: (manual 4.10.1).
GOAT_DEFAULT_WORKERS = 8

_SAFE_VALUE = re.compile(r"^[A-Za-z0-9_.,+\-]{1,64}$")


def normalize_goat_options(options: Optional[Mapping[str, Any]]) -> dict[str, Any]:
    """Validated ``%goat`` options, keyed by OPI ``BlockGoat`` field name.

    Input: a mapping of option -> value; keys are case-blind and may be a
    readable alias (:data:`GOAT_OPTION_ALIASES`).  Output: a new dict with
    canonical lower-case keys.  Only what the caller gives is written;
    every unset option keeps ORCA's default (manual Table 4.9).

    Refused (``ValueError``, before any file is written): a key ORCA does
    not know (``delfin.orca_knowledge.known_block_option``: OPI's block
    model plus measured parser acceptance), a value that is not a bool,
    number or short token (no quotes, spaces or newlines reach the input),
    and an ``nworkers`` that is not a positive multiple of 4 ("give any
    number multiple of 4 (number of temperatures)", Table 4.9) unless a
    ``templist`` defines the temperatures.
    """
    from delfin.orca_knowledge import known_block_option

    out: dict[str, Any] = {}
    for raw_key, value in dict(options or {}).items():
        key = str(raw_key).strip().lower()
        key = GOAT_OPTION_ALIASES.get(key, key)
        if not known_block_option("goat", key):
            raise ValueError(
                f"%goat option {raw_key!r} is not one ORCA knows "
                "(OPI block model / ORCA parser); refused"
            )
        if isinstance(value, bool) or isinstance(value, (int, float)):
            pass
        elif isinstance(value, str) and _SAFE_VALUE.match(value.strip()):
            value = value.strip()
        else:
            raise ValueError(
                f"%goat option {raw_key!r}: value {value!r} is not a bool, "
                "number or short token; refused"
            )
        out[key] = value
    if "nworkers" in out and "templist" not in out:
        n = out["nworkers"]
        if isinstance(n, bool) or not isinstance(n, int) or n <= 0 or n % 4:
            raise ValueError(
                f"%goat nworkers {n!r} must be a positive multiple of 4 "
                "(manual Table 4.9: number of temperatures)"
            )
    return out


def goat_output_files(base_name: str) -> dict[str, str]:
    """The files a GOAT run named ``base_name`` writes (manual section 4.10).

    ``finalensemble``: "BaseName.finalensemble.xyz" (4.10.3); ``global_minimum``:
    "Writing structure to HIS.globalminimum.xyz" (4.10.2); ``output``: the
    ORCA output; ``worker_outputs``: "BaseName.goat.globaliteration.workernumber.out",
    deleted after the run unless KEEPWORKERDATA is TRUE (4.10.2).
    """
    return {
        "output": f"{base_name}.out",
        "finalensemble": f"{base_name}.finalensemble.xyz",
        "global_minimum": f"{base_name}.globalminimum.xyz",
        "worker_outputs": f"{base_name}.goat.<global_iter>.<worker>.out",
    }


def read_mlips_goat(work_dir: Path | str, base_name: str = _BASE_NAME) -> dict:
    """Read an orca-mlips GOAT run: final ensemble and global minimum.

    Reuses :func:`delfin.tools.adapters.orca_goat_opi.read_goat_ensemble`
    (no second ensemble parser).  Output: that reader's dict plus
    ``energy_source`` and ``final`` -- MLIP energies are a screen, never a
    final number (mlp-method skill; a DFT re-ranking of the lowest
    conformers decides).
    """
    from delfin.tools.adapters.orca_goat_opi import read_goat_ensemble

    out = read_goat_ensemble(Path(work_dir), base_name)
    out["energy_source"] = "mlip"
    out["final"] = False
    return out


def build_orca_mlips_input(
    run_type: str,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    backend: str = "uma",
    model: Optional[str] = None,
    device: str = "cpu",
    solvent: Optional[str] = None,
    solvent_model: str = "alpb",
    precision: Optional[str] = None,
    hessian: Optional[str] = None,
    wrapper_path: Optional[str] = None,
    base_name: str = _BASE_NAME,
    work_dir: Path | str = ".",
    extra_keywords: Sequence[str] = (),
    goat_options: Optional[Mapping[str, Any]] = None,
    goat_variant: str = "GOAT",
    ncores: int = 1,
    neb_end_coordinates: Optional[Sequence[Sequence[float]]] = None,
    scan: Optional[dict] = None,
    xyz_name: Optional[str] = None,
    ext_params_text: Optional[str] = None,
) -> Path:
    """Write a ``! ExtOpt`` ORCA input driven by the orca-mlips command.

    ``wrapper_path`` + ``ext_params_text`` write another provider's ProgExt and
    Ext_Params instead (orca-external-tools, delfin/mlp_tools/oet.py); the
    rest of the input is the same.

    ``run_type``: one of :data:`MLIPS_RUN_TYPES`.
    ``hessian="two_step"`` (optts | freq | irc): the MLIP dumps its
    analytical Hessian (``--dump-hessian``) and the input reads it back
    (``%geom inhess Read / InHessName``, ORCA manual optimizations_TS.txt
    lines 394-400; for IRC ``%irc InitHess read / Hess_Filename``, irc.txt).
    ``goat``: ``goat_variant`` (Table 4.8) and ``goat_options`` (``%goat``,
    :func:`normalize_goat_options`); ``ncores`` becomes ``%pal nprocs``,
    which GOAT shares among its workers (manual 4.10.2).
    ``neb`` / ``neb-ts``: ``neb_end_coordinates`` (same atom order) is
    written to ``<base_name>_end.xyz`` and named in ``%neb
    NEB_End_XYZFile`` (neb.txt).  ``irc`` without a dumped Hessian:
    ``%irc InitHess calc_numfreq`` (gradients only).  ``scan``: the
    ``%geom Scan`` dict of :func:`delfin.mlp_tools.orca_extopt.build_extopt_input`.
    ``xyz_name``: the start-geometry file named in ``*xyzfile`` (default
    ``<base_name>.xyz``; a run that must not share the name with ORCA's
    own ``<base_name>.xyz`` passes another).

    Returns the written ``<base_name>.inp`` path.  opt/optts/freq with the
    new arguments at their defaults write the same bytes as before
    (tests/test_opi_orca_mlips_parity.py).
    """
    run_type = str(run_type).lower().strip()
    if run_type not in MLIPS_RUN_TYPES:
        raise ValueError(
            f"run_type {run_type!r} not one of {sorted(MLIPS_RUN_TYPES)}"
        )
    name = backend.strip().lower()
    if name not in _MLIPS_BACKENDS:
        raise ValueError(
            f"backend {backend!r} is not an orca-mlips plugin: use one of "
            f"{', '.join(_MLIPS_BACKENDS)}"
        )
    two_step = hessian == "two_step"
    if hessian is not None and not two_step:
        raise ValueError(f"hessian {hessian!r} not one of None | 'two_step'")
    if two_step and run_type not in ("optts", "freq", "irc"):
        raise ValueError(
            f"hessian='two_step' applies to optts | freq | irc, not {run_type!r}"
        )
    if goat_options and run_type != "goat":
        raise ValueError("goat_options apply to run_type='goat' only")
    if run_type in ("neb", "neb-ts"):
        if neb_end_coordinates is None or len(neb_end_coordinates) != len(elements):
            raise ValueError(
                f"run_type {run_type!r} needs neb_end_coordinates with one row "
                f"per atom ({len(elements)}): the product geometry (manual "
                "neb.txt, NEB_End_XYZFile)"
            )
    if run_type == "scan" and scan is None:
        raise ValueError("run_type='scan' requires the scan= dict")
    if int(ncores) < 1:
        raise ValueError(f"ncores {ncores!r} must be >= 1")

    progext = wrapper_path or resolve_progext(name)

    extra_ext: list[str] = []
    if run_type == "goat":
        variant = str(goat_variant).strip().upper()
        if variant not in GOAT_VARIANTS:
            raise ValueError(
                f"goat_variant {goat_variant!r} not one of {', '.join(GOAT_VARIANTS)}"
            )
        keywords = [variant]
    else:
        keywords = list(MLIPS_RUN_TYPES[run_type])
    goat_block = normalize_goat_options(goat_options) if run_type == "goat" else {}
    if two_step:
        extra_ext.append(f"--dump-hessian {_HESSIAN_FILE}")
        if run_type == "optts":
            keywords.append("Freq")  # Hessian-driven TS search: Freq + OptTS
    if run_type == "freq" and not two_step:
        # A bare MLIP Freq request still needs the dumped Hessian to read.
        extra_ext.append(f"--dump-hessian {_HESSIAN_FILE}")
    keywords.extend(str(k) for k in extra_keywords)

    ext_params = ext_params_text if ext_params_text is not None else ext_params_for(
        name, model=model, device=device, solvent=solvent,
        solvent_model=solvent_model, precision=precision, extra=extra_ext,
    )

    geom_block: Optional[list[str]] = None
    if two_step and run_type != "irc":
        geom_block = [
            "%geom",
            "  inhess Read",
            f'  InHessName "{_HESSIAN_FILE}"',
            "end",
        ]

    blocks: list[tuple[str, dict[str, Any]]] = []
    if goat_block:
        blocks.append(("goat", goat_block))
    if run_type in ("neb", "neb-ts"):
        end_name = f"{base_name}_end.xyz"
        _write_xyz(Path(work_dir) / end_name, elements, neb_end_coordinates)
        blocks.append(("neb", {"neb_end_xyzfile": end_name}))
    if run_type == "irc":
        if two_step:
            blocks.append(("irc", {"inithess": "read", "hess_filename": _HESSIAN_FILE}))
        else:
            blocks.append(("irc", {"inithess": "calc_numfreq"}))

    return _write_mlips_input(
        keywords, elements, coordinates, charge, mult, progext,
        ext_params, geom_block, base_name, work_dir,
        ncores=int(ncores), blocks=blocks, scan=scan, xyz_name=xyz_name,
    )


#: Option values written as quoted file names in the plain fallback (OPI
#: models them as ``InputFilePath`` and quotes them itself).
_FILE_OPTIONS = frozenset({"neb_end_xyzfile", "hess_filename", "ensemblefile"})


def _opi_blocks(blocks: Sequence[tuple[str, dict[str, Any]]]) -> list:
    """OPI block objects for ``(name, options)`` pairs; ImportError without OPI."""
    from opi.input.blocks.block_goat import BlockGoat
    from opi.input.blocks.block_irc import BlockIrc
    from opi.input.blocks.block_neb import BlockNeb

    classes = {"goat": BlockGoat, "neb": BlockNeb, "irc": BlockIrc}
    return [classes[name](**opts) for name, opts in blocks]


def _plain_block(name: str, opts: Mapping[str, Any]) -> list[str]:
    lines = [f"%{name}"]
    for key, value in opts.items():
        text = f'"{value}"' if key in _FILE_OPTIONS else str(value)
        lines.append(f"    {key} {text}")
    lines.append("end")
    return lines


def _write_xyz(path: Path, elements, coordinates) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    xyz_lines = [str(len(elements)), ""]
    for el, (x, y, z) in zip(elements, coordinates):
        xyz_lines.append(f"{el} {float(x):.10f} {float(y):.10f} {float(z):.10f}")
    path.write_text("\n".join(xyz_lines) + "\n")


def _write_mlips_input(keywords, elements, coordinates, charge, mult,
                       progext, ext_params, geom_block, base_name,
                       work_dir, *, ncores: int = 1,
                       blocks: Sequence[tuple[str, dict[str, Any]]] = (),
                       scan: Optional[dict] = None,
                       xyz_name: Optional[str] = None) -> Path:
    """OPI path (py3.12 venv) or plain-text fallback with identical content.

    The OPI path reuses orca_extopt's writer for the %method block, the
    %pal/%goat/%neb/%irc blocks and the *xyzfile rewrite; the %geom inhess
    block is appended after it.  With ``ncores=1`` and no blocks the
    fallback writes no ``%pal`` and the bytes of the earlier builder.
    """
    try:
        out_path = _build_with_opi(
            keywords, elements, coordinates, charge, mult, progext,
            ext_params, scan, base_name, work_dir,
            ncores=ncores, extra_blocks=_opi_blocks(blocks), xyz_name=xyz_name,
        )
        if geom_block:
            text = out_path.read_text()
            if "%geom" in text.lower():
                raise RuntimeError(
                    "OPI already wrote a %geom block; the inhess append "
                    "would create a duplicate block"
                )
            out_path.write_text(text.rstrip("\n") + "\n\n" + "\n".join(geom_block) + "\n")
        return out_path
    except ImportError:
        xyz_name = xyz_name or f"{base_name}.xyz"
        if scan is not None:
            from delfin.mlp_tools.orca_extopt import _write_plain_input

            out_path = Path(work_dir) / f"{base_name}.inp"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            _write_plain_input(
                keywords, elements, coordinates, charge, mult, progext,
                ext_params, scan, xyz_name, out_path,
            )
            _write_xyz(Path(work_dir) / xyz_name, elements, coordinates)
            if ncores > 1:
                head, rest = out_path.read_text().split("\n", 1)
                out_path.write_text(
                    f"{head}\n\n%pal\n    nprocs {ncores}\nend\n" + rest)
            return out_path
        lines = [f"! {' '.join(keywords)} ExtOpt", ""]
        if ncores > 1:
            lines.extend(["%pal", f"    nprocs {ncores}", "end", ""])
        lines.append("%method")
        lines.append(f'  ProgExt "{progext}"')
        lines.append(f'  Ext_Params "{ext_params}"')
        lines.append("end")
        lines.append("")
        for block_name, opts in blocks:
            lines.extend(_plain_block(block_name, opts))
            lines.append("")
        if geom_block:
            lines.extend(geom_block)
            lines.append("")
        lines.append(f"*xyzfile {charge} {mult} {xyz_name}")
        out_path = Path(work_dir) / f"{base_name}.inp"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(lines) + "\n")
        _write_xyz(Path(work_dir) / xyz_name, elements, coordinates)
        return out_path
