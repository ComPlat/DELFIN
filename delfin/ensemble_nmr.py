"""Helpers for CENSO/ANMR ensemble NMR workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence
import re


CENSO_NMR_SOLVENT_CHOICES: tuple[tuple[str, str], ...] = (
    ("chloroform", "chcl3"),
    ("water", "h2o"),
    ("acetonitrile", "acetonitrile"),
    ("dmso", "dmso"),
    ("methanol", "methanol"),
    ("thf", "thf"),
    ("toluene", "toluene"),
    ("dichloromethane", "ch2cl2"),
    ("acetone", "acetone"),
    ("gas phase", "gas"),
)

_CPCM_SOLVENT_ALIASES: dict[str, str] = {
    "chcl3": "chloroform",
    "chloroform": "chloroform",
    "h2o": "water",
    "water": "water",
    "ch2cl2": "dichloromethane",
    "dichloromethane": "dichloromethane",
    "acetonitrile": "acetonitrile",
    "dmso": "dmso",
    "methanol": "methanol",
    "thf": "thf",
    "toluene": "toluene",
    "acetone": "acetone",
    "gas": "gas",
}

_ANMR_SOLVENT_ALIASES: dict[str, str] = {
    "chloroform": "chcl3",
    "water": "h2o",
    "dichloromethane": "ch2cl2",
}


def normalize_cpcm_solvent_name(solvent: str) -> str:
    """Normalize UI solvent aliases to names accepted by CENSO/ORCA CPCM."""
    solvent_name = str(solvent or "chcl3").strip().lower() or "chcl3"
    return _CPCM_SOLVENT_ALIASES.get(solvent_name, solvent_name)


def normalize_anmr_solvent_name(solvent: str) -> str:
    """Normalize solvent names to conventional ANMR/TMS short labels."""
    cpcm_name = normalize_cpcm_solvent_name(solvent)
    return _ANMR_SOLVENT_ALIASES.get(cpcm_name, cpcm_name)


def xyz_body_to_coord_text(xyz_body: str | Iterable[str]) -> str:
    """Convert XYZ atom lines to TURBOMOLE ``coord`` format."""
    lines = xyz_body.splitlines() if isinstance(xyz_body, str) else list(xyz_body)
    bohr = 1.8897259886
    coord_lines: list[str] = []
    for line in lines:
        parts = str(line).split()
        if len(parts) < 4:
            continue
        try:
            elem = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except (TypeError, ValueError):
            continue
        coord_lines.append(
            f"  {x * bohr:14.8f}  {y * bohr:14.8f}  {z * bohr:14.8f}  {elem.lower()}"
        )
    if not coord_lines:
        raise ValueError("XYZ input contains no coordinate lines")
    return "$coord\n" + "\n".join(coord_lines) + "\n$end\n"


def build_censo_anmr_rc(
    *,
    solvent: str,
    resonance_frequency: float,
    active_nuclei: Sequence[str] = ("h",),
    orca_path: str = "",
    xtb_path: str = "",
    refinement_threshold: float | None = None,
) -> str:
    """Build a CENSO 3 rcfile for an ORCA-based 1H ensemble NMR run."""
    solvent_name = normalize_cpcm_solvent_name(solvent)
    gas_phase = solvent_name == "gas"
    active = ",".join(
        str(nucleus).strip().lower()
        for nucleus in active_nuclei
        if str(nucleus).strip()
    ) or "h"
    if refinement_threshold is None:
        refinement_threshold = 0.0
    orca_version = ""
    match = re.search(r"orca[_-]?(\d+)_(\d+)_(\d+)", str(orca_path or ""), re.IGNORECASE)
    if match:
        orca_version = ".".join(match.groups())

    path_lines = [
        f"orcapath = {orca_path}".rstrip(),
        f"orcaversion = {orca_version}".rstrip(),
        f"xtbpath = {xtb_path}".rstrip(),
        "tm =",
        "cosmotherm =",
        "cosmorssetup =",
    ]
    return (
        "[general]\n"
        "temperature = 298.15\n"
        "evaluate_rrho = True\n"
        "sm_rrho = gbsa\n"
        "imagthr = -100.0\n"
        "sthr = 50.0\n"
        f"solvent = {solvent_name}\n"
        f"gas_phase = {'True' if gas_phase else 'False'}\n"
        "copy_mo = True\n"
        "balance = True\n"
        "ignore_failed = False\n"
        "\n"
        "[prescreening]\n"
        "prog = orca\n"
        "func = pbe-d3\n"
        "basis = def2-sv(p)\n"
        "gfnv = gfn2\n"
        "threshold = 4.0\n"
        "template = False\n"
        "\n"
        "[screening]\n"
        "prog = orca\n"
        "func = r2scan-3c\n"
        "basis = def2-mtzvpp\n"
        "sm = cpcm\n"
        "gfnv = gfn2\n"
        "threshold = 3.5\n"
        "gsolv_included = True\n"
        "template = False\n"
        "\n"
        "[optimization]\n"
        "prog = orca\n"
        "func = r2scan-3c\n"
        "basis = def2-mtzvpp\n"
        "sm = cpcm\n"
        "gfnv = gfn2\n"
        "optcycles = 8\n"
        "maxcyc = 200\n"
        "optlevel = normal\n"
        "threshold = 3.0\n"
        "gradthr = 0.01\n"
        "hlow = 0.01\n"
        "macrocycles = True\n"
        "constrain = False\n"
        "xtb_opt = True\n"
        "template = False\n"
        "\n"
        "[refinement]\n"
        "prog = orca\n"
        "func = wb97x-d3\n"
        "basis = def2-tzvp\n"
        "sm = cpcm\n"
        "gfnv = gfn2\n"
        f"threshold = {float(refinement_threshold):.3f}\n"
        "gsolv_included = True\n"
        "template = False\n"
        "\n"
        "[nmr]\n"
        "run = True\n"
        "prog = orca\n"
        "func_s = pbe0-d4\n"
        "basis_s = def2-tzvp\n"
        "sm_s = cpcm\n"
        "func_j = pbe0-d4\n"
        "basis_j = def2-tzvp\n"
        "sm_j = cpcm\n"
        f"resonance_frequency = {float(resonance_frequency):.1f}\n"
        "ss_cutoff = 8.0\n"
        "fc_only = True\n"
        "shieldings = True\n"
        "couplings = True\n"
        f"h_active = {'True' if 'h' in active.split(',') else 'False'}\n"
        f"c_active = {'True' if 'c' in active.split(',') else 'False'}\n"
        f"f_active = {'True' if 'f' in active.split(',') else 'False'}\n"
        f"si_active = {'True' if 'si' in active.split(',') else 'False'}\n"
        f"p_active = {'True' if 'p' in active.split(',') else 'False'}\n"
        "template = False\n"
        "\n"
        "[rot]\n"
        "prog = tm\n"
        "func = pbe-d4\n"
        "basis = def2-svpd\n"
        "freq = [589.0, 633.0]\n"
        "template = False\n"
        "\n"
        "[uvvis]\n"
        "prog = orca\n"
        "func = wb97x-d4\n"
        "basis = def2-tzvp\n"
        "sm = smd\n"
        "nroots = 20\n"
        "template = False\n"
        "\n"
        "[paths]\n"
        + "\n".join(path_lines)
        + "\n"
    )


def build_anmrrc_text(
    *,
    solvent: str,
    resonance_frequency: float,
    shielding_ref_h: float,
    shielding_ref_c: float,
    linewidth: float = 1.0,
    active_h: bool = True,
    active_c: bool = False,
) -> str:
    """Build a minimal ``.anmrrc`` matching the CENSO ORCA NMR setup."""
    solvent_name = normalize_anmr_solvent_name(solvent)
    return (
        "7 8 XH acid atoms\n"
        f"ENSO qm= ORCA mf= {float(resonance_frequency):.1f} lw= {float(linewidth):.2f}  J= on S= on\n"
        f"TMS[{solvent_name}] pbe0-d4[CPCM]/def2-TZVP//r2scan-3c[CPCM]/def2-mTZVPP\n"
        f"1  {float(shielding_ref_h):.6f}    0.0    {1 if active_h else 0}\n"
        f"6  {float(shielding_ref_c):.6f}    0.0    {1 if active_c else 0}\n"
        "9  0.000000    0.0    0\n"
        "15 0.000000    0.0    0\n"
    )


def build_orca_reference_input(
    xyz_body: str,
    *,
    solvent: str,
    pal: int,
    maxcore: int,
    charge: int = 0,
    multiplicity: int = 1,
) -> str:
    """Build an ORCA NMR input for the TMS reference calculation."""
    solvent_name = normalize_cpcm_solvent_name(solvent)
    cpcm = "" if solvent_name == "gas" else f" CPCM({solvent_name})"
    coords_text = "\n".join(
        f"  {line.strip()}"
        for line in xyz_body.splitlines()
        if str(line).strip()
    )
    return (
        f"! PBE0 D4 DEF2-TZVP NMR{cpcm}\n"
        "\n"
        "%pal\n"
        f"  nprocs {max(1, int(pal))}\n"
        "end\n"
        "\n"
        f"%maxcore {max(100, int(maxcore))}\n"
        "\n"
        f"* xyz {int(charge)} {max(1, int(multiplicity))}\n"
        f"{coords_text}\n"
        "*\n"
    )


def write_text_file(path: str | Path, text: str) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text.rstrip() + "\n", encoding="utf-8")
    return target
