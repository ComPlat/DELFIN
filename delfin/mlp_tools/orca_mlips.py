"""orca-mlips as a second ExtOpt provider next to ``delfin-mlp-extopt``.

orca-mlips (github.com/t-0hmura/orca-mlips, MIT) ships four MLIP plugins
for ORCA's ExtTool/ProgExt interface with console entry points named
``uma``, ``orb``, ``mace``, ``aimnet2`` (prefixed aliases
``orca-mlips-uma`` ...).  This module

* detects the installed package and lists its real ProgExt entry points
  (never invents them, never installs anything — install hint only);
* builds the ``! ExtOpt`` ORCA input through the same OPI machinery as
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
import shutil
from pathlib import Path
from typing import Optional, Sequence

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


def _find_spec(name: str):
    import importlib.util

    return importlib.util.find_spec(name)


def detect_orca_mlips() -> dict:
    """Detect the installed orca-mlips package and its ProgExt entry points.

    Reads the distribution metadata only; never imports the MLIP stack
    (torch etc. must not load into this process), never installs anything.
    """
    installed = _find_spec("orca_mlips") is not None
    names: list[str] = []
    version = ""
    if installed:
        try:
            dist = importlib.metadata.distribution("orca-mlips")
            version = dist.version or ""
            eps = importlib.metadata.entry_points(group="console_scripts")
            names = sorted(
                ep.name for ep in eps if ep.name.startswith("orca-mlips-")
            )
        except importlib.metadata.PackageNotFoundError:
            installed = False
    return {
        "available": installed,
        "progext_names": names or (
            [f"orca-mlips-{b}" for b in _MLIPS_BACKENDS] if installed else []
        ),
        "version": version,
        "install_hint": (
            "pip install 'orca-mlips[uma]'   (or [orb], [mace], [aimnet2]; "
            "UMA needs its own environment -- conflicts with MACE/ORB)"
        ),
        "xtb_required_for_solvent": bool(shutil.which("xtb")) is False,
    }


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
    parts = [f"--backend {name}"]
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
) -> Path:
    """Write a ``! ExtOpt`` ORCA input driven by the orca-mlips command.

    ``run_type``: opt | optts | freq (neb/irc/goat/scan go through
    :func:`delfin.mlp_tools.orca_extopt.build_extopt_input` unchanged).
    ``hessian="two_step"``: the MLIP dumps its analytical Hessian
    (``--dump-hessian``) and the input reads it back (``%geom inhess Read /
    InHessName``, ORCA manual optimizations_TS.txt lines 394-400) for the
    OptTS search / Freq run.

    Returns the written ``<base_name>.inp`` path.
    """
    run_type = str(run_type).lower().strip()
    supported = ("opt", "optts", "freq")
    if run_type not in supported:
        raise ValueError(
            f"run_type {run_type!r} not one of {sorted(supported)} "
            "(for neb/irc/goat/scan use delfin-mlp-extopt)"
        )
    name = backend.strip().lower()
    if name not in _MLIPS_BACKENDS:
        raise ValueError(
            f"backend {backend!r} is not an orca-mlips plugin: use one of "
            f"{', '.join(_MLIPS_BACKENDS)}"
        )

    progext = wrapper_path or f"orca-mlips-{name}"

    extra_ext: list[str] = []
    # freq is not an EXTOPT_RUN_TYPES entry (orca_extopt has no freq run
    # type): a plain ! Freq ExtOpt run reads the MLIP's dumped Hessian.
    if run_type == "freq":
        keywords = ["Freq"]
    else:
        keywords = list(EXTOPT_RUN_TYPES[run_type])
    two_step = hessian == "two_step"
    if two_step:
        extra_ext.append(f"--dump-hessian {_HESSIAN_FILE}")
        if run_type == "optts":
            keywords.append("Freq")  # Hessian-driven TS search: Freq + OptTS
    if run_type == "freq" and not two_step:
        # A bare MLIP Freq request still needs the dumped Hessian to read.
        extra_ext.append(f"--dump-hessian {_HESSIAN_FILE}")

    ext_params = ext_params_for(
        name, model=model, device=device, solvent=solvent,
        solvent_model=solvent_model, precision=precision, extra=extra_ext,
    )

    geom_block: Optional[list[str]] = None
    if two_step:
        geom_block = [
            "%geom",
            "  inhess Read",
            f'  InHessName "{_HESSIAN_FILE}"',
            "end",
        ]

    return _write_mlips_input(
        keywords, elements, coordinates, charge, mult, progext,
        ext_params, geom_block, base_name, work_dir,
    )


def _write_mlips_input(keywords, elements, coordinates, charge, mult,
                       progext, ext_params, geom_block, base_name,
                       work_dir) -> Path:
    """OPI path (py3.12 venv) or plain-text fallback with identical content.

    The OPI path reuses orca_extopt's writer for the %method block and the
    *xyzfile rewrite; the %geom inhess block is appended after it.
    """
    try:
        out_path = _build_with_opi(
            keywords, elements, coordinates, charge, mult, progext,
            ext_params, None, base_name, work_dir,
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
        lines = [f"! {' '.join(keywords)} ExtOpt", ""]
        lines.append("%method")
        lines.append(f'  ProgExt "{progext}"')
        lines.append(f'  Ext_Params "{ext_params}"')
        lines.append("end")
        lines.append("")
        if geom_block:
            lines.extend(geom_block)
            lines.append("")
        xyz_name = f"{base_name}.xyz"
        lines.append(f"*xyzfile {charge} {mult} {xyz_name}")
        out_path = Path(work_dir) / f"{base_name}.inp"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(lines) + "\n")
        xyz_lines = [str(len(elements)), ""]
        for el, (x, y, z) in zip(elements, coordinates):
            xyz_lines.append(f"{el} {float(x):.10f} {float(y):.10f} {float(z):.10f}")
        (Path(work_dir) / xyz_name).write_text("\n".join(xyz_lines) + "\n")
        return out_path
