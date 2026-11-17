"""Chain-setup for the proton coordinator.

Prepares per-donor coordinator folders (``H3O+_coordinator/``,
``TFE_coordinator/``, ...) inside a finished DELFIN job directory, mirroring
``delfin/co2/chain_setup.py``:

1. Read the DELFIN ``CONTROL.txt``.
2. Take over species xyz + charge/multiplicity/BrokenSym for the requested
   redox species (state JSON -> OCCUPIER.txt -> CONTROL MANUALLY fallbacks,
   same precedence as the CO2 chain setup).
3. Charge/spin bookkeeping for complex + donor -> transfer fragments.
4. Emit the sub-CONTROL for the coordinator run.
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Species naming (identical convention to delfin/co2/chain_setup.py)
# ---------------------------------------------------------------------------

def species_delta_to_name(delta: int) -> str:
    """Map a redox delta to the DELFIN species folder/file name."""
    if delta == 0:
        return "initial"
    if delta < 0:
        return f"red_step_{abs(delta)}"
    return f"ox_step_{delta}"


# ---------------------------------------------------------------------------
# Spin handover (same precedence chain as CO2 chain_setup)
# ---------------------------------------------------------------------------

_MANUAL_KEY_MAP = {
    0: ("multiplicity_0", "BrokenSym_0"),
    -1: ("multiplicity_red1", "BrokenSym_red1"),
    -2: ("multiplicity_red2", "BrokenSym_red2"),
    +1: ("multiplicity_ox1", "BrokenSym_ox1"),
    +2: ("multiplicity_ox2", "BrokenSym_ox2"),
}


def _read_delfin_control(path: Path) -> Dict[str, str]:
    ctrl: Dict[str, str] = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("----"):
                continue
            if "=" in line:
                k, _, v = line.partition("=")
                ctrl[k.strip()] = v.strip()
    return ctrl


def _spin_from_state_json(job_dir: Path, species_delta: int
                          ) -> Tuple[Optional[int], Optional[str]]:
    """Primary source: DELFIN state JSON preferred spin per species."""
    state_path = job_dir / "DELFIN_Data.json"
    if not state_path.is_file():
        return None, None
    try:
        state = json.loads(state_path.read_text())
    except (OSError, ValueError):
        return None, None
    species_name = species_delta_to_name(species_delta)
    try:
        entry = state["species"][species_name]
        mult = entry.get("multiplicity")
        bs = entry.get("broken_sym") or None
        if mult is None:
            return None, None
        return int(mult), bs
    except (KeyError, TypeError, ValueError):
        return None, None


def _spin_from_occupier_txt(job_dir: Path, species_name: str
                            ) -> Tuple[Optional[int], Optional[str]]:
    """Fallback 1: extract spin from OCCUPIER.txt."""
    from delfin.copy_helpers import extract_preferred_spin

    occ_folder = job_dir / f"{species_name}_OCCUPIER"
    if not occ_folder.is_dir():
        return None, None
    return extract_preferred_spin(occ_folder)


def _spin_from_control_manual(ctrl: Dict[str, str], delta: int
                              ) -> Tuple[Optional[int], Optional[str]]:
    """Fallback 2: CONTROL.txt MANUALLY multiplicity/BrokenSym keys."""
    keys = _MANUAL_KEY_MAP.get(delta)
    if not keys:
        return None, None
    mult_key, broken_sym_key = keys
    mult_raw = ctrl.get(mult_key, "").strip()
    if not mult_raw:
        return None, None
    try:
        m = int(mult_raw)
    except ValueError:
        return None, None
    broken_sym = ctrl.get(broken_sym_key, "").strip() or None
    return m, broken_sym


def resolve_species_spin(job_dir: Path, species_delta: int
                         ) -> Tuple[Optional[int], Optional[str], str]:
    """Resolve (multiplicity, BrokenSym, source) with CO2-chain precedence."""
    job_dir = Path(job_dir)
    species_name = species_delta_to_name(species_delta)

    mult, bs = _spin_from_state_json(job_dir, species_delta)
    if mult is not None:
        return mult, bs, "state JSON"

    mult, bs = _spin_from_occupier_txt(job_dir, species_name)
    if mult is not None:
        return mult, bs, "OCCUPIER.txt"

    ctrl_path = job_dir / "CONTROL.txt"
    if ctrl_path.is_file():
        mult, bs = _spin_from_control_manual(_read_delfin_control(ctrl_path),
                                             species_delta)
        if mult is not None:
            return mult, bs, "CONTROL MANUALLY"

    return None, None, "none"


# ---------------------------------------------------------------------------
# Charge / spin bookkeeping: complex + donor -> transfer fragments
# ---------------------------------------------------------------------------

def bookkeeping(donor_key: str, complex_charge: int = 0,
                complex_multiplicity: int = 1,
                transfer_model: Optional[str] = None) -> Dict[str, Any]:
    """Charge/spin bookkeeping for the combined system and the fragments.

    Rules (chemically validated in the feature discussion):
    - protic donors: combined = complex + NEUTRAL donor; the H+ leaves, so
      the donor product is the conjugate base (charge -1) and the complex
      carries one MORE proton: [M-H] with complex_charge unchanged for the
      combined frame (total charge conserved: [M-H]+base charge sum equals
      complex+donor charge).
    - HAT donors (BIH, phenol default): combined = complex + neutral donor;
      the H atom leaves with one electron: donor product = neutral radical
      (multiplicity 2), complex keeps its charge.

    Returns combined charge/multiplicity plus per-fragment specs so the
    post-processing can compute reaction energies against relaxed fragments.
    """
    from delfin.proton.proton_coordinator import DONOR_TEMPLATES

    tpl = DONOR_TEMPLATES.get(donor_key)
    if tpl is None:
        raise KeyError(f"Unbekannter Donor '{donor_key}'.")

    transfer = transfer_model or tpl["transfer"]
    if transfer not in {"protic", "hat"}:
        raise ValueError(f"Unbekanntes Transfer-Modell '{transfer}'.")

    # Combined frame: complex + donor molecule. The donor is neutral for
    # the acid/HAT donors and +1 for H3O+ (the template says which).
    donor_neutral_charge = int(tpl.get("donor_charge", 0))
    combined_charge = complex_charge + donor_neutral_charge
    combined_multiplicity = _combine_multiplicity(complex_multiplicity, 1)

    if transfer == "protic":
        # H+ transferred: the complex GAINS the proton (charge +1 relative
        # to the complex), the donor becomes its conjugate base (-1).
        # Conservation: (complex+1) + (-1) = complex charge == combined.
        complex_product_charge = complex_charge + 1
        # protic_product overrides product when the template's default
        # transfer is HAT but the user forces the proton-transfer model
        # (e.g. phenol -> phenolate instead of phenoxyl radical).
        donor_product = dict(tpl.get("protic_product")
                             if tpl.get("protic_product") is not None
                             else tpl["product"])
    else:
        # H atom (H+ + e-) transferred: donor radical neutral (mult 2),
        # complex product keeps complex charge, multiplicity from pairing
        # handled by _combine_multiplicity with the radical not involved.
        complex_product_charge = complex_charge
        donor_product = dict(tpl["product"])   # charge 0, mult 2

    fragment_charge_sum = (complex_product_charge
                           + donor_product["charge"])
    if fragment_charge_sum != combined_charge:
        raise AssertionError(
            f"Charge bookkeeping broken for donor '{donor_key}': "
            f"combined {combined_charge} vs fragments {fragment_charge_sum}.")

    return {
        "scheme": transfer,
        "combined_charge": combined_charge,
        "combined_multiplicity": combined_multiplicity,
        "complex_product": {
            "name": "metal_hydride",
            "charge": complex_product_charge,
            "multiplicity": complex_multiplicity,
        },
        "donor_product": donor_product,
        "fragment_reaction": (
            f"complex({complex_charge}, {complex_multiplicity}) + donor "
            f"-> [M-H]({complex_product_charge}) + "
            f"{donor_product['name']}({donor_product['charge']})"),
    }


def _combine_multiplicity(mult_a: int, mult_b: int) -> int:
    """Couple two multiplicities (2S+1) to the spin-aligned range bounds.

    Returns the LOWEST possible multiplicity of the coupled pair (the
    antiferromagnetic bound), which is the conservative default for a scan
    start; BrokenSym in the CONTROL is the tool to reach the aligned state.
    """
    s_a = (mult_a - 1) / 2.0
    s_b = (mult_b - 1) / 2.0
    s_min = abs(s_a - s_b)
    return int(2 * s_min + 1)


# ---------------------------------------------------------------------------
# Sub-CONTROL generation
# ---------------------------------------------------------------------------

_TRANSFER_KEYS = (
    "xyz", "charge", "multiplicity", "broken_sym",
    "functional", "disp_corr", "ri_jkx", "aux_jk",
    "main_basisset", "metal_basisset",
    "first_coordination_sphere_metal_basisset",
    "first_coordination_sphere_scale",
    "second_coordination_sphere_metal_basisset",
    "second_coordination_sphere_scale",
    "implicit_solvation_model", "solvent",
    "PAL", "maxcore",
    "orientation_distance", "rot_step_deg", "rot_range_deg",
    "scan_end", "scan_steps",
    "parallel_orientation_scan", "max_workers",
)


def build_control(parent_args: Dict[str, str], donor_key: str,
                  control_path: Optional[Path] = None,
                  outdir: Optional[str] = None,
                  species_delta: int = 0) -> Path:
    """Write the per-donor coordinator CONTROL with ``proton_``-prefixed
    forwarded keys (mirrors the CO2 chain setup forwarding)."""
    from delfin.proton.proton_coordinator import CANONICAL_DONOR_DIRS

    if control_path is None:
        control_path = Path(outdir or ".") / "CONTROL.txt"
    control_path = Path(control_path)
    control_path.parent.mkdir(parents=True, exist_ok=True)

    lines = ["# Proton coordinator CONTROL (auto-generated)",
             "# " + CANONICAL_DONOR_DIRS.get(donor_key, donor_key),
             "------------------------------------"]
    # Plain keys, NOT proton_-prefixed: the coordinator's minimal CONTROL
    # parser reads them un-prefixed (CO2 convention, _build_co2_control) —
    # a prefixed key would silently fall back to defaults (charge=0 etc.).
    # Only proton_donor / proton_species_delta are namespaced.
    for key in _TRANSFER_KEYS:
        val = parent_args.get(key)
        if val not in (None, ""):
            lines.append(f"{key}={val}")
    lines.append(f"proton_donor={donor_key}")
    lines.append(f"proton_species_delta={species_delta}")
    control_path.write_text("\n".join(lines) + "\n")
    return control_path


def write_default_control(control_path: str = "CONTROL.txt",
                          overwrite: bool = False) -> str:
    """Standalone CONTROL template (``--define``)."""
    template = """# Proton coordinator - standalone CONTROL template
# Donor selection: comma-separated; each gets its own <Donor>_coordinator/
# folder. Valid: tfe | bih | phenol | hoac | h3o  (aliases accepted)
------------------------------------
proton_donor=tfe,phenol,bih,h3o
xyz=complex.xyz
charge=0
multiplicity=1

# Solvation (never a bare proton: donor + implicit solvent)
------------------------------------
implicit_solvation_model=CPCM
solvent=MeCN

# Method
------------------------------------
functional=PBE0
disp_corr=D4
main_basisset=def2-SVP
metal_basisset=def2-TZVP

# Geometry / scan
------------------------------------
orientation_distance=3.0
rot_step_deg=10
rot_range_deg=180
scan_end=1.5
scan_steps=15

# Resources
------------------------------------
PAL=8
maxcore=3800
"""
    if os_exists(control_path) and not overwrite:
        print(f"[INFO] {control_path} existiert bereits - nichts geschrieben.")
        return control_path
    with open(control_path, "w", newline="\n") as f:
        f.write(template)
    print(f"[OK] geschrieben: {control_path}")
    return control_path


def os_exists(path: str) -> bool:
    import os
    return os.path.exists(path)


# ---------------------------------------------------------------------------
# prepare_proton_coordination: build the folders in a finished job dir
# ---------------------------------------------------------------------------

def prepare_proton_coordination(job_dir, species_delta: int = 0,
                                donor: str = "tfe,phenol,bih,h3o") -> list:
    """Create per-donor coordinator folders inside *job_dir* after a DELFIN
    run (mirror of delfin.co2.chain_setup.prepare_co2_coordination)."""
    from delfin.proton.proton_coordinator import CANONICAL_DONOR_DIRS

    job_dir = Path(job_dir)
    ctrl = _read_delfin_control(job_dir / "CONTROL.txt")
    mult, bs, source = resolve_species_spin(job_dir, species_delta)
    if mult is None:
        raise RuntimeError(
            f"Keine Multiplizitaet fuer Species-Delta {species_delta} "
            f"gefunden (state JSON / OCCUPIER.txt / MANUALLY leer).")

    species_name = species_delta_to_name(species_delta)
    species_xyz = job_dir / f"{species_name}.xyz"
    if not species_xyz.is_file():
        raise FileNotFoundError(f"Species-XYZ nicht gefunden: {species_xyz}")

    parent_args: Dict[str, str] = dict(ctrl)
    parent_args["charge"] = ctrl.get("charge", "0")
    parent_args["multiplicity"] = str(mult)
    if bs:
        parent_args["broken_sym"] = bs

    created = []
    for raw in str(donor).split(","):
        raw = raw.strip()
        if not raw:
            continue
        from delfin.proton.proton_coordinator import resolve_donor_name
        canon = resolve_donor_name(raw)
        if canon is None:
            print(f"[WARN] Unbekannter Donor '{raw}' - uebersprungen.")
            continue
        folder = job_dir / CANONICAL_DONOR_DIRS[canon]
        folder.mkdir(exist_ok=True)
        shutil.copy(species_xyz, folder / "complex.xyz")
        bk = bookkeeping(canon, complex_charge=int(parent_args["charge"]),
                         complex_multiplicity=mult)
        (folder / "bookkeeping.json").write_text(json.dumps(bk, indent=2))
        build_control(parent_args, canon, control_path=folder / "CONTROL.txt",
                      species_delta=species_delta)
        created.append(folder)
    return created
