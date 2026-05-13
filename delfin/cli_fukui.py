"""CLI tool: compute atomic Fukui indices via three ORCA single-points.

Implements the OPI workflow described at
https://www.faccts.de/docs/opi/2.0/docs/contents/notebooks/atomic_fukui_indices.html

Input: SMILES string OR path to an .xyz file.

Workflow:
    1. Build / load geometry (optional pre-OPT for raw input).
    2. Run three ORCA single-points at identical geometry:
         neutral.inp (charge=0,  mult=1)
         anion.inp   (charge=-1, mult=2, UKS)
         cation.inp  (charge=+1, mult=2, UKS)
    3. Parse Mulliken or Loewdin atomic charges per state.
    4. Compute atomic f+/f-/f0 Fukui indices.
    5. Generate 3 density cubes via orca_plot.
    6. Subtract → 3 Fukui difference cubes.
    7. Write fukui_result.json + .fukui_job marker.

Pipeline mode is always *classic* regardless of any CONTROL.txt setting:
no OCCUPIER branching, no smart-recalc skipping, no special workflow
handling. The CLI is intentionally self-contained so that it can be
invoked locally OR on a SLURM compute node via submit_delfin.sh.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin import fukui as fukui_lib
from delfin.common.logging import get_logger
from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks, resolve_maxiter

logger = get_logger(__name__)


# Default ORCA settings mirror the OPI tutorial (B3LYP / def2-SVP / D3).
# These get used when neither CLI args nor CONTROL.txt supply a value.
DEFAULT_FUNCTIONAL = "B3LYP"
DEFAULT_BASIS = "def2-SVP"
DEFAULT_DISPERSION = "D3BJ"
DEFAULT_SOLVATION = ""        # "" = gas phase
DEFAULT_SOLVENT = ""
DEFAULT_PAL = 4
DEFAULT_MAXCORE = 3000

# ORCA's plot_type IDs for the total electron density vary slightly
# between versions; we try this sequence in order until one succeeds.
DENSITY_PLOT_TYPE_IDS: Tuple[int, ...] = (41, 40, 42, 44)


# ---------------------------------------------------------------------------
# Input detection
# ---------------------------------------------------------------------------

def _looks_like_xyz_path(value: str) -> bool:
    if not value:
        return False
    p = Path(value)
    return p.exists() and p.suffix.lower() == ".xyz"


def _looks_like_smiles(value: str) -> bool:
    """Heuristic: not a path AND contains at least one element-y character."""
    if not value:
        return False
    if Path(value).exists():
        return False
    if "/" in value or "\\" in value:
        return False
    return any(c.isalpha() for c in value)


def _smiles_to_xyz_text(smiles: str) -> str:
    """Convert SMILES → DELFIN-format xyz text via the quick converter."""
    from delfin.smiles_converter import smiles_to_xyz_quick

    xyz, err = smiles_to_xyz_quick(smiles)
    if xyz is None or err:
        raise RuntimeError(f"smiles_to_xyz_quick failed: {err!r}")
    return xyz


def _ensure_xyz_header(xyz_text: str) -> str:
    """Add the standard XYZ header (count + comment) if missing.

    The DELFIN quick converter returns header-less text; ORCA's
    ``* xyzfile`` directive expects a standard file with the leading
    atom-count and comment lines.
    """
    lines = xyz_text.strip().splitlines()
    if not lines:
        raise ValueError("empty xyz text")
    first = lines[0].strip()
    if first.isdigit():
        return xyz_text if xyz_text.endswith("\n") else xyz_text + "\n"
    n = sum(1 for line in lines if line.strip() and len(line.split()) >= 4)
    return f"{n}\nfukui geometry\n" + "\n".join(lines) + "\n"


def _strip_xyz_header(xyz_text: str) -> str:
    """Return coordinate-only lines (no atom-count, no comment)."""
    lines = xyz_text.strip().splitlines()
    if not lines:
        return ""
    if lines[0].strip().isdigit():
        return "\n".join(lines[2:])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CONTROL.txt loading (lenient: missing file = use defaults)
# ---------------------------------------------------------------------------

def _maybe_load_control(workdir: Path) -> Dict[str, object]:
    """Try to find and parse CONTROL.txt and return its dict (or minimal fallback).

    Looked up in: ``workdir``, ``workdir.parent``, current directory.
    Failures are logged but never raised. DELFIN placeholder values
    (``"[SOLVENT]"``, ``"[METAL]"`` etc.) are stripped here so the
    canonical bang-line builder never emits e.g. ``CPCM([SOLVENT])``.
    """
    candidates = [workdir / "CONTROL.txt", workdir.parent / "CONTROL.txt", Path.cwd() / "CONTROL.txt"]
    for path in candidates:
        try:
            if path.exists():
                from delfin.config import parse_control_text, _is_placeholder_value

                cfg = parse_control_text(
                    path.read_text(encoding="utf-8"), keep_steps_literal=True
                )
                filtered = {
                    k: v for k, v in cfg.items()
                    if not _is_placeholder_value(v)
                }
                logger.info("loaded ORCA defaults from %s", path)
                return filtered
        except Exception as exc:  # noqa: BLE001
            logger.warning("failed to parse %s: %s", path, exc)
    logger.info("no CONTROL.txt found; using built-in OPI-style defaults")
    return {}


def _resolve_setting(
    cli_value: Optional[str],
    control: Dict[str, object],
    control_key: str,
    default: str,
) -> str:
    if cli_value not in (None, ""):
        return str(cli_value)
    raw = control.get(control_key)
    if raw not in (None, ""):
        return str(raw)
    return default


# ---------------------------------------------------------------------------
# ORCA input builders
# ---------------------------------------------------------------------------

def _resolve_fukui_config(
    workdir: Path,
    settings: Dict[str, str],
    geom_xyz_path: Optional[Path] = None,
) -> Dict[str, object]:
    """Merge CLI-resolved settings with CONTROL.txt-style keys for the canonical builder.

    The returned dict is shaped to feed :func:`delfin.xyz_io._build_bang_line`
    and :class:`delfin.common.orca_blocks.OrcaInputBuilder` directly.
    """
    raw = _maybe_load_control(workdir)
    merged: Dict[str, object] = {}
    for k, v in raw.items():
        merged[k] = v

    # CLI-resolved values override CONTROL.txt
    merged["functional"] = settings["functional"]
    merged["main_basisset"] = settings["basis"]
    merged["PAL"] = int(settings["pal"])
    merged["maxcore"] = int(settings["maxcore"])

    # Map our short flags onto canonical CONTROL.txt keys
    if settings.get("dispersion"):
        merged["disp_corr"] = settings["dispersion"]
    if settings.get("solvation"):
        merged["implicit_solvation_model"] = settings["solvation"]
    if settings.get("solvent"):
        merged["solvent"] = settings["solvent"]

    # Required-but-typical defaults for SP jobs
    merged.setdefault("ri_jkx", "")
    merged.setdefault("initial_guess", "")
    merged.setdefault("metal_basisset", merged.get("metal_basisset"))

    if geom_xyz_path is not None:
        merged["_fukui_geom_path"] = str(geom_xyz_path)
    return merged


def _detect_metals(geom_xyz_path: Optional[Path]) -> List[str]:
    """Return transition-metal symbols present in ``geom_xyz_path`` (or []).

    Uses :func:`delfin.utils.search_transition_metals` — the same scanner
    that the classic pipeline runs before building ``initial.inp``.
    """
    if geom_xyz_path is None or not Path(geom_xyz_path).exists():
        return []
    from delfin.utils import search_transition_metals

    return list(search_transition_metals(str(geom_xyz_path)) or [])


def _build_orca_input(
    *,
    job_type: str,
    charge: int,
    mult: int,
    coords: str,
    config: Dict[str, object],
    found_metals: List[str],
    cube_basename: Optional[str] = None,
    cube_grid: int = 80,
) -> str:
    """Build a Fukui ORCA input via DELFIN's canonical bang-line + OrcaInputBuilder.

    For ``job_type='SP'`` we drop FREQ and the geom-opt token (Fukui needs a
    pure single-point). For ``job_type='OPT'`` we re-enable both — that
    matches the behaviour of ``initial.inp`` so a Fukui pre-OPT inherits
    every CONTROL.txt knob the classic pipeline respects (RI, dispersion,
    implicit solvation, relativity, metal basis assignment, ...).

    If ``cube_basename`` is set, append a ``%plots`` block that asks ORCA
    to write the total electron density to ``<cube_basename>.cube`` as a
    side product of the SCF — avoids the brittle post-hoc ``orca_plot``
    dance whose menu IDs differ between ORCA versions.
    """
    from delfin.utils import resolve_level_of_theory
    from delfin.xyz_io import _build_bang_line, _implicit_token

    solvent = str(config.get("solvent", "")).strip()
    main, metal, rel_token, aux_jk = resolve_level_of_theory(
        found_metals,
        config,
        config.get("main_basisset"),
        config.get("metal_basisset"),
    )
    implicit = _implicit_token(config, solvent)

    is_opt = job_type.upper() == "OPT"
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=is_opt,
        geom_key="geom_opt" if is_opt else "",
    )

    # Open-shell SPs (anion / cation of a closed-shell neutral) need UKS.
    if mult != 1 and "UKS" not in bang.upper():
        bang = bang.replace("!", "! UKS", 1)

    output_blocks = collect_output_blocks(config, allow=True)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(int(config["maxcore"]), int(config["PAL"]), resolve_maxiter(config))
    builder.add_blocks(output_blocks)

    if cube_basename and not is_opt:
        plots_block = (
            "%plots\n"
            "  Format Gaussian_Cube\n"
            f"  dim1 {int(cube_grid)}\n"
            f"  dim2 {int(cube_grid)}\n"
            f"  dim3 {int(cube_grid)}\n"
            f'  ElDens("{cube_basename}.cube");\n'
            "end\n"
        )
        builder.add_block(plots_block)

    lines = builder.lines
    lines.append(f"* xyz {charge} {mult}\n")
    lines.append(coords.rstrip() + "\n")
    lines.append("*\n")
    return "".join(lines)


# ---------------------------------------------------------------------------
# ORCA execution wrappers
# ---------------------------------------------------------------------------

def _run_orca_job(inp_path: Path, out_path: Path) -> bool:
    """Execute a single ORCA job, returning True on success."""
    from delfin.orca import run_orca

    logger.info("running ORCA: %s -> %s", inp_path.name, out_path.name)
    t0 = time.time()
    ok = run_orca(str(inp_path), str(out_path), working_dir=inp_path.parent)
    elapsed = time.time() - t0
    if ok:
        logger.info("ORCA %s finished in %.1fs", inp_path.name, elapsed)
    else:
        logger.error("ORCA %s FAILED after %.1fs", inp_path.name, elapsed)
    return ok


def _extract_optimized_xyz(out_path: Path) -> Optional[str]:
    """Locate the ``<base>.xyz`` file ORCA writes after a successful OPT.

    ORCA writes the final optimized geometry to ``<base>.xyz`` in the
    same directory as the input file.
    """
    base = out_path.with_suffix("")
    candidate = base.with_suffix(".xyz")
    if candidate.exists():
        return candidate.read_text(encoding="utf-8")
    return None


def _run_orca_plot_density(gbw_path: Path, out_cube: Path) -> bool:
    """Generate an electron-density cube next to ``gbw_path``.

    Uses :func:`delfin.reporting.esp_report._orca_plot_run` and tries the
    plot_type IDs that ORCA uses across versions. After success, copies
    the produced cube to ``out_cube``.
    """
    from delfin.reporting.esp_report import _orca_plot_run, _pick_density_cube

    if not gbw_path.exists():
        logger.error("gbw not found: %s", gbw_path)
        return False

    before = time.time()
    for plot_type in DENSITY_PLOT_TYPE_IDS:
        ok = _orca_plot_run(gbw_path, f"1\n{plot_type}\n0\n11\n12\n")
        if not ok:
            continue
        produced = _pick_density_cube(gbw_path.parent, modified_after=before)
        if produced and produced != out_cube:
            shutil.move(str(produced), str(out_cube))
        if out_cube.exists():
            return True
    return False


# ---------------------------------------------------------------------------
# Workflow driver
# ---------------------------------------------------------------------------

def _prepare_geometry_via_pipeline(
    workdir: Path,
    *,
    pre_opt: bool,
    config: Dict[str, object],
    settings: Dict[str, str],
) -> Tuple[Path, str]:
    """Materialize the Fukui geometry by reusing DELFIN's canonical pipeline:

    1. ``normalize_input_file`` reads ``input.txt`` (SMILES or XYZ block)
       and produces ``start.txt`` with the configured smiles_converter
       (QUICK / ARCHITECTOR / GUPPY).
    2. If ``pre_opt`` is on, build ``initial.inp`` via
       :func:`delfin.xyz_io.read_and_modify_file_1` — the same builder
       the classic pipeline uses for ``initial.inp`` — then run ORCA and
       use ORCA's final ``initial.xyz`` as the Fukui geometry.
    3. Otherwise, build ``initial.xyz`` from start.txt directly.
    """
    control_path = workdir / "CONTROL.txt"
    if not control_path.exists():
        raise SystemExit(
            f"input.txt mode requires CONTROL.txt in workdir; not found at {control_path}"
        )

    config.setdefault("input_file", "input.txt")
    config.setdefault("smiles_converter", config.get("smiles_converter", "QUICK"))

    from delfin.workflows.pipeline import normalize_input_file

    normalized_input = normalize_input_file(config, control_path)
    normalized_path = Path(normalized_input)
    if not normalized_path.is_absolute():
        normalized_path = workdir / normalized_path

    # `normalize_input_file` puts the XYZ geometry into start.txt for
    # SMILES inputs; for plain XYZ inputs the geometry stays in input.txt
    # (or its renamed form). Use start.txt if present, else fall back.
    start_txt = workdir / "start.txt"
    geom_source = start_txt if start_txt.exists() else normalized_path
    if not geom_source.exists():
        raise SystemExit(
            f"could not locate normalized geometry (looked at start.txt and {normalized_path})"
        )

    if not pre_opt:
        # Convert the DELFIN-format start.txt to a standard XYZ file.
        initial_xyz = workdir / "initial.xyz"
        raw_lines = [
            ln for ln in geom_source.read_text(encoding="utf-8").splitlines()
            if ln.strip() and ln.strip() != "*"
        ]
        initial_xyz.write_text(
            f"{len(raw_lines)}\nfukui geometry (no pre-OPT)\n" + "\n".join(raw_lines) + "\n",
            encoding="utf-8",
        )
        return initial_xyz, "user_xyz"

    # ---- Pre-OPT via the canonical initial.inp builder ----
    from delfin.utils import search_transition_metals, resolve_level_of_theory
    from delfin.xyz_io import read_and_modify_file_1

    metals = list(search_transition_metals(str(geom_source)) or [])
    main_basis, metal_basis, _rel, _aux = resolve_level_of_theory(
        metals, config, config.get("main_basisset"), config.get("metal_basisset"),
    )
    solvent = str(config.get("solvent", "")).strip()
    try:
        base_charge = int(str(config.get("charge", 0)).strip())
    except ValueError:
        base_charge = 0
    try:
        base_multiplicity = int(str(config.get("multiplicity_global_opt", 1)).strip())
    except ValueError:
        base_multiplicity = 1

    initial_inp = workdir / "initial.inp"
    read_and_modify_file_1(
        str(geom_source),
        str(initial_inp),
        base_charge,
        base_multiplicity,
        solvent,
        metals,
        metal_basis,
        main_basis,
        config,
        "",   # broken_sym
    )

    initial_out = workdir / "initial.out"
    if not _run_orca_job(initial_inp, initial_out):
        raise SystemExit("pre-OPT ORCA job failed (see initial.out)")

    # ORCA writes the final optimised geometry to initial.xyz next to initial.out
    final_xyz = workdir / "initial.xyz"
    if not final_xyz.exists():
        raise SystemExit(f"pre-OPT finished but optimized initial.xyz not found at {final_xyz}")
    return final_xyz, "opt"


def _run_three_singlepoints(
    workdir: Path,
    geom_xyz: Path,
    *,
    settings: Dict[str, str],
    request_cubes: bool = True,
) -> Dict[str, Path]:
    """Build, submit, and verify the neutral / anion / cation SPs.

    When ``request_cubes`` is True (the default), each SP's input file
    asks ORCA to emit the total electron density as a Gaussian-Cube via
    a ``%plots`` block. That sidesteps ``orca_plot``'s version-specific
    interactive menu and lets the density-cube files appear as a
    deterministic side-product of each SCF.
    """
    coords = _strip_xyz_header(geom_xyz.read_text(encoding="utf-8"))
    config = _resolve_fukui_config(workdir, settings, geom_xyz_path=geom_xyz)
    found_metals = _detect_metals(geom_xyz)
    states = [
        ("neutral", 0, 1),
        ("anion", -1, 2),
        ("cation", +1, 2),
    ]
    outputs: Dict[str, Path] = {}
    for name, charge, mult in states:
        sub = workdir / name
        sub.mkdir(parents=True, exist_ok=True)
        inp = sub / f"{name}.inp"
        out = sub / f"{name}.out"
        cube_basename = f"density_{name}" if request_cubes else None
        inp.write_text(
            _build_orca_input(
                job_type="SP",
                charge=charge, mult=mult,
                coords=coords,
                config=config,
                found_metals=found_metals,
                cube_basename=cube_basename,
            ),
            encoding="utf-8",
        )
        if not _run_orca_job(inp, out):
            raise SystemExit(f"ORCA SP for {name} failed")
        outputs[name] = out
    return outputs


def _generate_density_cubes(workdir: Path) -> Dict[str, Path]:
    """Collect the density cubes ORCA emitted via the SP ``%plots`` block.

    ORCA writes ``<base>.cube`` next to the ``.inp`` (so inside the
    per-state subdirectory). We move/rename it next to fukui_result.json
    so the viewer can pick it up without traversing subfolders.
    """
    cubes: Dict[str, Path] = {}
    for name in ("neutral", "anion", "cation"):
        sub = workdir / name
        target = workdir / f"density_{name}.cube"
        # ORCA may put the cube in sub/, sub/with-basename suffix, or at
        # the top level depending on path quoting in %plots.
        candidates = [
            sub / f"density_{name}.cube",
            sub / f"{name}.cube",
            workdir / f"density_{name}.cube",
            workdir / f"{name}.cube",
        ]
        located = next((c for c in candidates if c.exists()), None)
        if located is None:
            # As a last resort, search the subdir for any *.cube ORCA wrote.
            extras = list(sub.glob("*.cube"))
            located = extras[0] if extras else None
        if located is None:
            logger.warning("density cube not found for %s in %s", name, sub)
            continue
        if located != target:
            shutil.move(str(located), str(target))
        cubes[name] = target
    return cubes


def _generate_fukui_cubes(workdir: Path, density: Dict[str, Path]) -> Dict[str, Path]:
    """Compute Fukui-function difference cubes from the three density cubes."""
    out: Dict[str, Path] = {}
    if "anion" in density and "neutral" in density:
        target = workdir / "fukui_plus.cube"
        fukui_lib.subtract_cubes(
            density["anion"], density["neutral"], target,
            title="Fukui f+ = rho(anion) - rho(neutral)",
        )
        out["f_plus"] = target
    if "neutral" in density and "cation" in density:
        target = workdir / "fukui_minus.cube"
        fukui_lib.subtract_cubes(
            density["neutral"], density["cation"], target,
            title="Fukui f- = rho(neutral) - rho(cation)",
        )
        out["f_minus"] = target
    if "anion" in density and "cation" in density:
        target = workdir / "fukui_zero.cube"
        fukui_lib.subtract_cubes(
            density["anion"], density["cation"], target,
            scale=0.5,
            title="Fukui f0 = (rho(anion) - rho(cation)) / 2",
        )
        out["f_zero"] = target
    return out


# ---------------------------------------------------------------------------
# CLI entry
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute atomic Fukui indices via three ORCA single-points "
            "(N, N+1, N-1) at identical geometry, following the OPI recipe."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # SMILES input (always pre-OPT, B3LYP/def2-SVP defaults):
  delfin-fukui --input "C=O" --workdir ./fukui_HCHO

  # XYZ input, skip pre-OPT (geometry is already optimized):
  delfin-fukui --input mol.xyz --workdir ./fukui_run --no-pre-opt

  # Loewdin charges, custom method:
  delfin-fukui --input mol.xyz --scheme loewdin \\
      --functional PBE0 --basis def2-TZVP --workdir ./fukui_run
""",
    )
    parser.add_argument(
        "--input", default=None,
        help=(
            "Input geometry (SMILES or .xyz). Optional: if --workdir already "
            "contains input.txt (dashboard convention), the CLI uses that."
        ),
    )
    parser.add_argument(
        "--workdir", default=None,
        help="Working directory for ORCA jobs and results (default: ./fukui_<auto>).",
    )
    parser.add_argument(
        "--pre-opt", dest="pre_opt", action="store_true",
        help="Run an OPT+FREQ before the Fukui SPs (forced for SMILES input).",
    )
    parser.add_argument(
        "--no-pre-opt", dest="pre_opt", action="store_false",
        help="Skip OPT; assume input geometry is already optimized.",
    )
    parser.set_defaults(pre_opt=None)
    parser.add_argument(
        "--scheme", choices=list(fukui_lib.CHARGE_SCHEMES), default="mulliken",
        help="Atomic-charge scheme (default: mulliken).",
    )
    parser.add_argument("--functional", default=None, help="ORCA functional (default: from CONTROL.txt or B3LYP).")
    parser.add_argument("--basis", default=None, help="ORCA basis set (default: def2-SVP).")
    parser.add_argument("--dispersion", default=None, help="ORCA dispersion correction (default: D3BJ).")
    parser.add_argument("--solvation", default=None, help="Solvation model (CPCM/SMD) or empty for gas phase.")
    parser.add_argument("--solvent", default=None, help="Solvent name for the solvation model.")
    parser.add_argument("--pal", type=int, default=None, help="Number of ORCA worker threads (default: 4).")
    parser.add_argument("--maxcore", type=int, default=None, help="ORCA per-core RAM in MB (default: 3000).")
    parser.add_argument(
        "--skip-cubes", action="store_true",
        help="Skip orca_plot density-cube generation (faster, no isosurface viz).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    return parser


def _resolve_workdir(input_value: Optional[str], cli_workdir: Optional[str]) -> Path:
    if cli_workdir:
        return Path(cli_workdir).expanduser().resolve()
    if not input_value:
        # No --input and no --workdir → assume current directory is the
        # dashboard-style workdir (must contain input.txt + CONTROL.txt).
        return Path.cwd().resolve()
    if _looks_like_xyz_path(input_value):
        stem = Path(input_value).stem
    else:
        safe = "".join(c if c.isalnum() else "_" for c in input_value)[:32] or "fukui"
        stem = safe
    return (Path.cwd() / f"fukui_{stem}").resolve()


def _resolve_pre_opt(input_value: str, cli_pre_opt: Optional[bool]) -> bool:
    if _looks_like_smiles(input_value):
        if cli_pre_opt is False:
            logger.warning("SMILES input requires pre-OPT; --no-pre-opt ignored.")
        return True
    if cli_pre_opt is None:
        return False
    return cli_pre_opt


def _seed_input_txt(workdir: Path, input_value: Optional[str]) -> None:
    """Ensure ``workdir/input.txt`` exists in DELFIN-canonical form.

    DELFIN's convention for ``input.txt`` (matching every other workflow):
        * SMILES: single line (no leading int, no comment)
        * XYZ:    coord lines only (header stripped — no atom count, no comment)

    The dashboard already pre-seeds it correctly via ``clean_input_data``;
    here we mirror that for direct CLI users so the file ends up identical.
    """
    target = workdir / "input.txt"
    if target.exists() and target.stat().st_size > 0:
        return
    if not input_value:
        raise SystemExit(
            f"no input.txt found in {workdir} and --input not given; "
            "either pre-seed input.txt (dashboard convention) or pass --input."
        )
    if _looks_like_xyz_path(input_value):
        raw = Path(input_value).read_text(encoding="utf-8")
    else:
        raw = input_value
    from delfin.dashboard.input_processing import clean_input_data

    cleaned, kind = clean_input_data(raw)
    if kind == "empty" or not cleaned:
        raise SystemExit(f"could not classify --input {input_value!r}: empty after cleaning")
    target.write_text(cleaned + "\n", encoding="utf-8")
    logger.info("seeded input.txt (%s) at %s", kind, target)


def _seed_control_txt(
    workdir: Path,
    args_settings: Dict[str, Optional[str]],
) -> None:
    """Synthesize a minimal CONTROL.txt from CLI flags when none is present.

    Required so :func:`delfin.config.read_control_file` validation passes
    when the user only supplies CLI arguments (e.g. test or scripted use).
    """
    control_path = workdir / "CONTROL.txt"
    if control_path.exists():
        return
    lines = [
        f"functional={args_settings.get('functional') or DEFAULT_FUNCTIONAL}",
        f"main_basisset={args_settings.get('basis') or DEFAULT_BASIS}",
        f"disp_corr={args_settings.get('dispersion') or DEFAULT_DISPERSION}",
        f"PAL={args_settings.get('pal') or DEFAULT_PAL}",
        f"maxcore={args_settings.get('maxcore') or DEFAULT_MAXCORE}",
        "charge=0",
        "multiplicity_global_opt=1",
        "method=classic",
        "smiles_converter=QUICK",
        "calc_initial=yes",
        "geom_opt=Opt",
        "freq_type=FREQ",
    ]
    solvation = args_settings.get("solvation")
    solvent = args_settings.get("solvent")
    if solvation:
        lines.append(f"implicit_solvation_model={solvation}")
    if solvent:
        lines.append(f"solvent={solvent}")
    control_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("seeded minimal CONTROL.txt at %s", control_path)


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point for ``delfin-fukui``."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    workdir = _resolve_workdir(args.input, args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    logger.info("workdir: %s", workdir)

    # Materialize input.txt + CONTROL.txt in the workdir so the rest of
    # the routine can rely on the standard DELFIN file layout. The
    # dashboard already pre-seeds both; CLI users get sensible defaults.
    _seed_input_txt(workdir, args.input)

    control = _maybe_load_control(workdir)
    settings = {
        "functional":  _resolve_setting(args.functional, control, "functional", DEFAULT_FUNCTIONAL),
        "basis":       _resolve_setting(args.basis, control, "main_basisset", DEFAULT_BASIS),
        "dispersion":  _resolve_setting(args.dispersion, control, "disp_corr", DEFAULT_DISPERSION),
        "solvation":   _resolve_setting(args.solvation, control, "implicit_solvation_model", DEFAULT_SOLVATION),
        "solvent":     _resolve_setting(args.solvent, control, "solvent", DEFAULT_SOLVENT),
        "pal":         _resolve_setting(str(args.pal) if args.pal else None, control, "PAL", str(DEFAULT_PAL)),
        "maxcore":     _resolve_setting(str(args.maxcore) if args.maxcore else None, control, "maxcore", str(DEFAULT_MAXCORE)),
    }
    logger.info("ORCA settings: %s", settings)
    _seed_control_txt(workdir, settings)

    # Reload full + validated config now that CONTROL.txt is guaranteed to exist.
    try:
        from delfin.config import read_control_file
        config = read_control_file(str(workdir / "CONTROL.txt"))
    except Exception as exc:  # noqa: BLE001
        logger.warning("read_control_file failed (%s); falling back to lenient parse", exc)
        config = dict(_maybe_load_control(workdir))

    pre_opt = _resolve_pre_opt(args.input or str(workdir / "input.txt"), args.pre_opt)
    logger.info("pre-OPT: %s", pre_opt)

    geom_xyz, geometry_origin = _prepare_geometry_via_pipeline(
        workdir, pre_opt=pre_opt, config=config, settings=settings,
    )
    fukui_geom = workdir / "fukui_geom.xyz"
    if geom_xyz != fukui_geom:
        fukui_geom.write_text(geom_xyz.read_text(encoding="utf-8"), encoding="utf-8")

    outputs = _run_three_singlepoints(workdir, fukui_geom, settings=settings)

    symbols, q_neutral = fukui_lib.read_atoms_and_charges(outputs["neutral"], scheme=args.scheme)
    q_anion = fukui_lib.read_charges(outputs["anion"], scheme=args.scheme)
    q_cation = fukui_lib.read_charges(outputs["cation"], scheme=args.scheme)
    indices = fukui_lib.compute_fukui_from_charges(q_neutral, q_anion, q_cation)

    cubes_summary: Dict[str, str] = {}
    if not args.skip_cubes:
        density = _generate_density_cubes(workdir)
        fukui_cubes = _generate_fukui_cubes(workdir, density)
        cubes_summary = {k: str(v.name) for k, v in {**density, **fukui_cubes}.items()}
    else:
        logger.info("skipping cube generation (--skip-cubes)")

    fukui_lib.write_fukui_result_json(
        workdir,
        atoms=symbols,
        scheme=args.scheme,
        q_neutral=q_neutral, q_anion=q_anion, q_cation=q_cation,
        fukui=indices,
        orca_settings=settings,
        geometry_origin=geometry_origin,
        extra={"cubes": cubes_summary, "input": args.input},
    )
    fukui_lib.write_marker(workdir)

    print(f"\n✓ Fukui analysis written to {workdir}/fukui_result.json")
    print(f"  marker: {workdir}/{fukui_lib.FUKUI_MARKER}")
    if cubes_summary:
        print(f"  cubes:  {len(cubes_summary)} files in {workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
