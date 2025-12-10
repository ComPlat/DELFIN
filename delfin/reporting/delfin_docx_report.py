"""Combined DELFIN Word report assembled from DELFIN_Data.json and generated plots."""

from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Optional
import json

from delfin.common.logging import get_logger
from delfin.ir_spectrum import parse_ir_spectrum

logger = get_logger(__name__)

try:
    from docx import Document
    from docx.shared import Inches, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    DOCX_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    DOCX_AVAILABLE = False

HARTREE_TO_EV = 27.211386245988


def _translate_state(orca_state: str) -> str:
    """Translate ORCA notation (0-1A/1-3A) to S0/S1/T1 for readability."""
    match = re.match(r'(\d+)-([13])A', str(orca_state))
    if not match:
        return str(orca_state)

    root_number = int(match.group(1))
    multiplicity = match.group(2)
    if multiplicity == "1":
        return f"S{root_number}"
    if multiplicity == "3":
        return f"T{root_number}"
    return str(orca_state)


def _wavelength_to_rgb(wavelength: float) -> tuple[int, int, int]:
    """Map wavelength (nm) to approximate RGB (copied from AFP plotting logic)."""
    if wavelength < 380:
        r, g, b = 0.5, 0.0, 0.5
    elif wavelength < 440:
        r = -(wavelength - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif wavelength < 490:
        r = 0.0
        g = (wavelength - 440) / (490 - 440)
        b = 1.0
    elif wavelength < 510:
        r = 0.0
        g = 1.0
        b = -(wavelength - 510) / (510 - 490)
    elif wavelength < 580:
        r = (wavelength - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif wavelength < 645:
        r = 1.0
        g = -(wavelength - 645) / (645 - 580)
        b = 0.0
    elif wavelength < 781:
        r, g, b = 1.0, 0.0, 0.0
    else:
        r, g, b = 0.5, 0.0, 0.0

    if wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif wavelength > 700:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
    else:
        factor = 1.0

    r = int(max(0, min(1, r * factor)) * 255)
    g = int(max(0, min(1, g * factor)) * 255)
    b = int(max(0, min(1, b * factor)) * 255)
    return r, g, b


@dataclass
class ReportAssets:
    """Container for optional plot assets referenced in the final report."""

    afp_png: Optional[Path] = None
    smiles_png: Optional[Path] = None
    uv_vis_pngs: Dict[str, Path] | None = None  # keyed by state name, e.g., "S0", "S1", "T1"
    ir_png: Optional[Path] = None
    energy_level_png: Optional[Path] = None


def _add_key_value_table(doc: Document, title: str, rows: Iterable[tuple[str, str]]) -> None:
    """Render a simple two-column table."""
    doc.add_heading(title, level=2)
    table = doc.add_table(rows=0, cols=2)
    table.style = "Light Grid Accent 1"
    for key, value in rows:
        row_cells = table.add_row().cells
        row_cells[0].text = str(key)
        row_cells[1].text = str(value)


def _add_plot_if_exists(doc: Document, title: str, image_path: Optional[Path]) -> None:
    if not image_path:
        return
    if not image_path.exists():
        logger.warning("Plot '%s' not found at %s", title, image_path)
        return
    doc.add_heading(title, level=2)
    doc.add_picture(str(image_path), width=Inches(6.5))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER


def _xyz_body_to_mol(xyz_lines: list[str]):
    """Convert XYZ body (no header lines) to RDKit Mol via bond perception."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Geometry import Point3D

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return None

    atoms = []
    coords = []
    for line in xyz_lines:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 4:
            return None
        symbol = parts[0]
        try:
            x, y, z = map(float, parts[1:4])
        except Exception:
            return None
        atoms.append(symbol)
        coords.append(Point3D(x, y, z))

    if not atoms:
        return None

    try:
        mol = Chem.RWMol()
        for sym in atoms:
            mol.AddAtom(Chem.Atom(sym))
        conf = Chem.Conformer(len(atoms))
        for idx, pt in enumerate(coords):
            conf.SetAtomPosition(idx, pt)
        mol.AddConformer(conf)
        rdDetermineBonds.DetermineBonds(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def _mol_to_smiles(mol):
    try:
        from rdkit import Chem
        smi = Chem.MolToSmiles(mol)
        return smi
    except Exception:
        return None


def _generate_smiles_image(project_dir: Path) -> Optional[Path]:
    """Render SMILES from input.txt (SMILES or XYZ body) to a PNG using RDKit."""
    smiles_file = project_dir / "input.txt"
    if not smiles_file.exists():
        logger.warning("SMILES source 'input.txt' not found for RDKit rendering")
        return None

    try:
        lines = [line.strip() for line in smiles_file.read_text(encoding="utf-8", errors="ignore").splitlines()]
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read input.txt: %s", exc)
        return None

    smiles = ""
    mol_from_xyz = None

    # Detect if content looks like XYZ body (no header lines)
    xyz_candidate = [ln for ln in lines if ln and not ln.startswith("#")]
    if xyz_candidate and all(len(ln.split()) >= 4 and ln.split()[0][0].isalpha() for ln in xyz_candidate):
        # Try XYZ -> Mol
        mol_from_xyz = _xyz_body_to_mol(xyz_candidate)

    if mol_from_xyz:
        smiles = _mol_to_smiles(mol_from_xyz) or ""
    else:
        # Fallback: treat first non-empty line as SMILES
        smiles = next((ln for ln in lines if ln), "")

    if not smiles:
        logger.warning("No SMILES or usable XYZ detected in input.txt")
        return None

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, Draw

        RDLogger.DisableLog("rdApp.*")
    except Exception as exc:  # noqa: BLE001
        logger.warning("RDKit not available for SMILES rendering: %s", exc)
        return None

    mol = Chem.MolFromSmiles(smiles) if not mol_from_xyz else mol_from_xyz
    if mol is None:
        logger.warning("Invalid SMILES in input.txt: %s", smiles)
        return None

    try:
        # Hide hydrogens in the depiction for a cleaner schematic
        mol_no_h = Chem.RemoveHs(mol, updateExplicitCount=True)
        AllChem.Compute2DCoords(mol_no_h)
        img = Draw.MolToImage(mol_no_h, size=(500, 400))
        out_path = project_dir / "SMILES.png"
        img.save(out_path)
        return out_path
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to render SMILES to image: %s", exc)
        return None


def _render_color_chip(rgb: tuple[int, int, int], project_dir: Path, label: str) -> Optional[Path]:
    """Create a small PNG color chip for Word compatibility."""
    try:
        from PIL import Image
    except Exception:
        return None

    try:
        chip = Image.new("RGB", (60, 40), rgb)
        safe_label = (
            label.replace("→", "-")
            .replace("₀", "0")
            .replace("₁", "1")
            .replace("₂", "2")
            .replace("₃", "3")
            .replace("₄", "4")
            .replace("₅", "5")
            .replace("₆", "6")
        )
        out_path = project_dir / f"color_chip_{safe_label}.png"
        chip.save(out_path)
        return out_path
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to create color chip: %s", exc)
        return None


def _format_energy_block(optimization: Dict[str, Any]) -> str:
    def _fmt(value, fmt: str) -> Optional[str]:
        try:
            return fmt.format(float(value))
        except Exception:
            return None

    parts: list[str] = []
    hartree = _fmt(optimization.get("hartree"), "{:.6f}")
    if hartree:
        parts.append(f"{hartree} Eh")
    ev = _fmt(optimization.get("eV"), "{:.3f}")
    if ev:
        parts.append(f"{ev} eV")
    kj = _fmt(optimization.get("kJ_mol"), "{:.2f}")
    if kj:
        parts.append(f"{kj} kJ/mol")
    return ", ".join(parts) if parts else "n/a"


def _energy_ev(optimization: Dict[str, Any]) -> Optional[float]:
    if not optimization:
        return None
    if optimization.get("eV") is not None:
        try:
            return float(optimization["eV"])
        except Exception:
            return None
    if optimization.get("hartree") is not None:
        try:
            return float(optimization["hartree"]) * HARTREE_TO_EV
        except Exception:
            return None
    return None


def _build_summary_text(data: Dict[str, Any], project_dir: Path) -> tuple[Optional[str], list[tuple[str, float, tuple[int, int, int]]]]:
    meta = data.get("metadata", {}) or {}
    name = meta.get("NAME") or meta.get("name") or project_dir.name
    functional = meta.get("functional") or "unknown functional"
    basis = meta.get("basis_set") or "unknown basis"
    aux_basis = meta.get("auxiliary_basis")
    ri_method = meta.get("ri_method")
    solvent_model = meta.get("implicit_solvation")
    solvent = meta.get("solvent")
    dispersion = meta.get("dispersion_correction")

    gs = data.get("ground_state_S0", {}) or {}
    gs_opt = gs.get("optimization", {}) or {}
    gs_orb = gs.get("orbitals", {}) or {}
    gs_dip = gs.get("dipole_moment", {}) or {}
    s0_e = _energy_ev(gs_opt)

    homo = gs_orb.get("homo_eV")
    lumo = gs_orb.get("lumo_eV")
    gap = gs_orb.get("gap_eV")

    # Vibrational highlights from IR spectrum if available
    vib_frequencies: list[float] = []
    negative_freqs = 0
    ir_file = None
    for candidate in ["ESD/S0.out", "S0.out", "initial.out"]:
        cand_path = project_dir / candidate
        if cand_path.exists():
            ir_file = cand_path
            break
    if ir_file:
        try:
            modes = parse_ir_spectrum(ir_file)
            negative_freqs = sum(1 for m in modes if m.frequency_cm1 < 0)
            vib_frequencies = sorted(
                [m.frequency_cm1 for m in modes],
                key=lambda x: x,
            )
            # Keep top 5 by intensity for summary
            modes_sorted = sorted(modes, key=lambda m: m.intensity_km_mol, reverse=True)
            vib_frequencies = [m.frequency_cm1 for m in modes_sorted[:5]]
        except Exception as exc:  # noqa: BLE001
            logger.warning("IR parsing failed for summary: %s", exc)

    excited = data.get("excited_states", {}) or {}
    excited_count = len(excited)
    s1_e = _energy_ev((excited.get("S1") or {}).get("optimization", {}) or {})
    t1_e = _energy_ev((excited.get("T1") or {}).get("optimization", {}) or {})

    def rel_energy(ref, val):
        if ref is None or val is None:
            return None
        return val - ref

    s1_gap = rel_energy(s0_e, s1_e)
    t1_gap = rel_energy(s0_e, t1_e)
    est = None
    if s1_gap is not None and t1_gap is not None:
        est = abs(s1_gap - t1_gap)

    # Absorption highlights from S0 transitions
    s0_abs = (gs.get("tddft_absorption") or {}).get("transitions", []) or []
    if s0_abs:
        top_abs = sorted(
            s0_abs,
            key=lambda t: t.get("oscillator_strength", 0.0) or 0.0,
            reverse=True,
        )[:2]
        abs_peaks = [t.get("wavelength_nm") for t in top_abs if t.get("wavelength_nm")]
    else:
        abs_peaks = []

    def _pick_transition_nm(transitions: list[Dict[str, Any]], state_a: str, state_b: str) -> Optional[float]:
        """Pick strongest transition between two states (order-agnostic)."""
        best: tuple[float, float] | None = None
        target_set = {state_a, state_b}
        for t in transitions or []:
            from_state = _translate_state(t.get("from_state", ""))
            to_state = _translate_state(t.get("to_state", ""))
            if {from_state, to_state} != target_set:
                continue
            try:
                fosc = float(t.get("oscillator_strength", 0.0) or 0.0)
            except Exception:
                fosc = 0.0
            wl = t.get("wavelength_nm")
            if wl is None:
                continue
            if best is None or fosc > best[0]:
                best = (fosc, wl)
        return best[1] if best else None

    # Key transitions for AFP plot: S0↔S1, S1↔S0, T1↔S0 (order-agnostic)
    s0_to_s1_nm = _pick_transition_nm(s0_abs, "S0", "S1")

    s1_transitions = (excited.get("S1", {}) or {}).get("tddft_from_geometry", {}).get("transitions", [])
    s1_to_s0_nm = _pick_transition_nm(s1_transitions, "S0", "S1")

    t1_transitions = (excited.get("T1", {}) or {}).get("tddft_from_geometry", {}).get("transitions", [])
    t1_to_s0_nm = _pick_transition_nm(t1_transitions, "S0", "T1")

    # Redox potentials from DELFIN.txt and OCCUPIER data
    delfin_summary = data.get("delfin_summary", {}) or {}
    redox_potentials_fc = delfin_summary.get("redox_potentials_vs_fc", {})
    oxidized_states = data.get("oxidized_states", {}) or {}
    reduced_states = data.get("reduced_states", {}) or {}
    occupier_data = data.get("occupier", {}) or {}

    def fmt(val, suffix=""):
        if val is None:
            return "n/a"
        try:
            return f"{float(val):.2f}{suffix}"
        except Exception:
            return str(val) + suffix

    # Build summary
    parts = [
        "The SMILES of the target system were converted into 3D coordinates using RDKit. "
        "A global geometry optimization with the GOAT algorithm and GFN-xTB in ORCA provided the starting structure for the high-level DFT calculations automated by DELFIN.",
        f"The calculation of optimized structure, vibrational frequencies and excited states for the system '{name}' is presented, with automated analysis provided by the DELFIN software package."
    ]

    # Theory level with all details
    theory_parts = [f"{functional}/{basis}"]
    if ri_method:
        theory_parts.append(ri_method)
    if aux_basis:
        theory_parts.append(aux_basis)
    if dispersion:
        theory_parts.append(dispersion)
    if solvent_model and solvent:
        theory_parts.append(f"{solvent_model}({solvent})")
    elif solvent_model:
        theory_parts.append(solvent_model)

    parts.append(f"The calculations were performed at the {' '.join(theory_parts)} level of theory.")

    # Ground state properties (NO SCF energy)
    parts.append(f"HOMO/LUMO energies are {fmt(homo, ' eV')} and {fmt(lumo, ' eV')}, yielding a gap of {fmt(gap, ' eV')}.")
    parts.append(f"The permanent dipole moment is {fmt(gs_dip.get('magnitude_debye'), ' D')}.")

    if vib_frequencies:
        vib_list = ", ".join(f"{freq:.0f}" for freq in vib_frequencies)
        parts.append(f"The most intense vibrational frequencies are approximately {vib_list} cm⁻¹.")
    if ir_file:
        parts.append(f"Negative frequencies detected: {negative_freqs}.")

    # Excited states + TDDFT settings from S1
    tddft_settings = (excited.get("S1", {}) or {}).get("tddft_settings", {})
    if tddft_settings:
        tokens = []
        for key in ["nroots", "maxdim", "followiroot"]:
            if key in tddft_settings:
                tokens.append(f"{key.capitalize()}={tddft_settings[key]}")
        settings_str = f" using TDDFT ({', '.join(tokens)})" if tokens else " using TDDFT"
    else:
        settings_str = ""
    if excited:
        state_list = ", ".join(sorted(excited.keys()))
        parts.append(f"In total, {excited_count} excited states were calculated{settings_str}: {state_list}.")
    else:
        parts.append(f"In total, {excited_count} excited states were calculated{settings_str}.")
    if est is not None:
        parts.append(f"The singlet–triplet energy gap ΔEₛₜ(S₁–T₁) is {fmt(est, ' eV')} (absolute difference).")

    # Absorption peaks
    if abs_peaks:
        peak_str = " and ".join(f"{w:.0f} nm" for w in abs_peaks if isinstance(w, (int, float)))
        if peak_str:
            parts.append(f"The most intense absorption peaks (S₀→Sₙ) are at {peak_str}.")

    # AFP spectrum wavelengths (absorption/fluorescence/phosphorescence)
    color_boxes: list[tuple[str, float, tuple[int, int, int]]] = []
    if any(v is not None for v in (s0_to_s1_nm, s1_to_s0_nm, t1_to_s0_nm)):
        def fmt_nm(val):
            try:
                return f"{float(val):.0f} nm"
            except Exception:
                return "n/a"

        afp_parts = [
            f"S₀→S₁ (absorption): {fmt_nm(s0_to_s1_nm)}",
            f"S₁→S₀ (fluorescence): {fmt_nm(s1_to_s0_nm)}",
            f"T₁→S₀ (phosphorescence): {fmt_nm(t1_to_s0_nm)}",
        ]
        parts.append(f"Transitions: {', '.join(afp_parts)}.")

        # Collect color chips for these wavelengths
        for label, wl in [
            ("S₀→S₁", s0_to_s1_nm),
            ("S₁→S₀", s1_to_s0_nm),
            ("T₁→S₀", t1_to_s0_nm),
        ]:
            if wl:
                color_boxes.append((label, wl, _wavelength_to_rgb(float(wl))))

    # Redox potentials
    if oxidized_states or reduced_states:
        redox_parts = []
        def _step_label(name: str) -> str:
            match = re.search(r"_step_(\d+)", name)
            return f" step {match.group(1)}" if match else ""

        for state_name, state_data in oxidized_states.items():
            opt = state_data.get("optimization", {}) or {}
            charge = opt.get("charge")
            mult = opt.get("multiplicity")
            energy = opt.get("hartree")
            if charge is not None and energy is not None:
                step_info = _step_label(state_name)
                # Try to get OCCUPIER info for this state
                occ_key = state_name.replace("_step_", "")
                occ_info = occupier_data.get(occ_key, {})
                preferred_idx = occ_info.get("preferred_index")
                if preferred_idx is not None and occ_info.get("entries"):
                    preferred_entry = next((e for e in occ_info["entries"] if e.get("index") == preferred_idx), None)
                    if preferred_entry:
                        mult = preferred_entry.get("multiplicity", mult)
                        unpaired_alpha = preferred_entry.get("unpaired_alpha")
                        unpaired_beta = preferred_entry.get("unpaired_beta")
                        spin_cont = preferred_entry.get("spin_contamination")
                        info = f"oxidation{step_info} (charge {charge:+d}, multiplicity {mult}"
                        if unpaired_alpha is not None and unpaired_beta is not None:
                            info += f", α/β electrons: {unpaired_alpha}/{unpaired_beta}"
                        if spin_cont is not None:
                            info += f", spin contamination: {spin_cont:.3f}"
                        info += ")"
                        redox_parts.append(info)
                    else:
                        redox_parts.append(f"oxidation{step_info} (charge {charge:+d}, multiplicity {mult})")
                else:
                    redox_parts.append(f"oxidation{step_info} (charge {charge:+d}, multiplicity {mult})")

        for state_name, state_data in reduced_states.items():
            opt = state_data.get("optimization", {}) or {}
            charge = opt.get("charge")
            mult = opt.get("multiplicity")
            energy = opt.get("hartree")
            if charge is not None and energy is not None:
                step_info = _step_label(state_name)
                occ_key = state_name.replace("_step_", "")
                occ_info = occupier_data.get(occ_key, {})
                preferred_idx = occ_info.get("preferred_index")
                if preferred_idx is not None and occ_info.get("entries"):
                    preferred_entry = next((e for e in occ_info["entries"] if e.get("index") == preferred_idx), None)
                    if preferred_entry:
                        mult = preferred_entry.get("multiplicity", mult)
                        unpaired_alpha = preferred_entry.get("unpaired_alpha")
                        unpaired_beta = preferred_entry.get("unpaired_beta")
                        spin_cont = preferred_entry.get("spin_contamination")
                        info = f"reduction{step_info} (charge {charge:+d}, multiplicity {mult}"
                        if unpaired_alpha is not None and unpaired_beta is not None:
                            info += f", α/β electrons: {unpaired_alpha}/{unpaired_beta}"
                        if spin_cont is not None:
                            info += f", spin contamination: {spin_cont:.3f}"
                        info += ")"
                        redox_parts.append(info)
                    else:
                        redox_parts.append(f"reduction{step_info} (charge {charge:+d}, multiplicity {mult})")
                else:
                    redox_parts.append(f"reduction{step_info} (charge {charge:+d}, multiplicity {mult})")

        if redox_parts:
            parts.append(f"Redox calculations were performed for {', '.join(redox_parts)}.")

    # Add redox potentials vs Fc+/Fc from DELFIN.txt
    if redox_potentials_fc:
        redox_strs = []
        label_map = {
            "E_red": "E_red1",
            "E_red_2": "E_red2",
            "E_red_3": "E_red3",
            "E_ox": "E_ox1",
            "E_ox_2": "E_ox2",
            "E_ox_3": "E_ox3",
        }
        # Reductions
        for key in ['E_red', 'E_red_2', 'E_red_3']:
            if key in redox_potentials_fc:
                val = redox_potentials_fc[key]
                label = label_map.get(key, key)
                redox_strs.append(f"{label}: {val:.3f} V")
        # Oxidations
        for key in ['E_ox', 'E_ox_2', 'E_ox_3']:
            if key in redox_potentials_fc:
                val = redox_potentials_fc[key]
                label = label_map.get(key, key)
                redox_strs.append(f"{label}: {val:.3f} V")
        # Reference
        if 'E_ref' in redox_potentials_fc:
            redox_strs.append(f"E_ref: {redox_potentials_fc['E_ref']:.3f} V")

        if redox_strs:
            parts.append(f"Redox potentials (vs. Fc⁺/Fc): {', '.join(redox_strs)}.")

    return " ".join(parts), color_boxes


def _style_header_row(row) -> None:
    for cell in row.cells:
        for run in cell.paragraphs[0].runs:
            run.font.bold = True


def _format_state_with_subscript(paragraph, state_text: str) -> None:
    """Format state text with subscript numbers (e.g., S₀, S₁, T₁)."""
    import re
    # Match patterns like S0, S1, T1, T2, etc.
    match = re.match(r'([ST])(\d+)', state_text)
    if match:
        letter = match.group(1)
        number = match.group(2)
        # Add letter without bold
        run1 = paragraph.add_run(letter)
        run1.font.bold = False
        # Add number as subscript
        run2 = paragraph.add_run(number)
        run2.font.subscript = True
        run2.font.bold = False
    else:
        # Fallback for non-standard formats
        run = paragraph.add_run(state_text)
        run.font.bold = False


def _add_heading_with_subscript(doc: Document, title: str, level: int = 2) -> None:
    """Add a heading with subscript formatting for state labels (e.g., S₀, T₁)."""
    import re

    heading = doc.add_heading(level=level)

    # Split title by state patterns (S0, S1, T1, etc.)
    parts = re.split(r'([ST]\d+)', title)

    for part in parts:
        if re.match(r'[ST]\d+', part):
            # This is a state label - format with subscript
            match = re.match(r'([ST])(\d+)', part)
            if match:
                letter = match.group(1)
                number = match.group(2)
                run1 = heading.add_run(letter)
                run2 = heading.add_run(number)
                run2.font.subscript = True
        else:
            # Regular text
            heading.add_run(part)


def _format_significant_figures(value, sig_figs: int = 3) -> str:
    """Format a number to specified significant figures."""
    try:
        val = float(value)
        if val == 0:
            return "0.00"
        from math import log10, floor
        return f"{val:.{sig_figs - 1 - int(floor(log10(abs(val))))}f}"
    except (ValueError, TypeError):
        return str(value)


def _add_state_table(doc: Document, title: str, states: Dict[str, Any]) -> None:
    if not states:
        return
    doc.add_heading(title, level=2)
    table = doc.add_table(rows=1, cols=6)
    table.style = "Light Grid Accent 1"
    headers = ["State", "Type", "Charge", "Multiplicity", "Energy (Eh)", "ZPE (Eh)"]
    for idx, text in enumerate(headers):
        cell = table.rows[0].cells[idx]
        cell.text = text
    _style_header_row(table.rows[0])

    for state_name, entry in sorted(states.items()):
        opt = entry.get("optimization", {}) or {}
        thermo = entry.get("thermochemistry", {}) or {}
        row = table.add_row().cells

        # Format State with subscript
        row[0].text = ""
        _format_state_with_subscript(row[0].paragraphs[0], state_name)

        row[1].text = str(entry.get("_type", ""))
        row[2].text = str(opt.get("charge", ""))
        row[3].text = str(opt.get("multiplicity", ""))

        # Only show Hartree energy with 3 significant figures
        hartree = opt.get("hartree")
        row[4].text = _format_significant_figures(hartree, sig_figs=3) if hartree is not None else ""

        row[5].text = (
            _format_significant_figures(thermo.get('zero_point_energy_hartree'), sig_figs=3)
            if "zero_point_energy_hartree" in thermo
            else ""
        )


def _add_transition_table(
    doc: Document,
    title: str,
    transitions: list[Dict[str, Any]],
    limit: int | None = None,
) -> None:
    if not transitions:
        return
    _add_heading_with_subscript(doc, title, level=3)
    table = doc.add_table(rows=1, cols=5)
    table.style = "Light Grid Accent 1"
    headers = ["From", "To", "Energy (eV)", "Wavelength (nm)", "fosc"]
    for idx, text in enumerate(headers):
        cell = table.rows[0].cells[idx]
        cell.text = text
    _style_header_row(table.rows[0])

    def _wl_nm(entry: Dict[str, Any]) -> float:
        try:
            return float(entry.get("wavelength_nm", float("inf")))
        except Exception:
            return float("inf")

    # Sort by wavelength (nm) ascending
    sorted_transitions = sorted(transitions, key=_wl_nm)
    max_items = len(sorted_transitions) if limit is None else limit
    for trans in sorted_transitions[:max_items]:
        row = table.add_row().cells

        # Format "From" state with subscript
        from_state = _translate_state(trans.get("from_state", ""))
        row[0].text = ""  # Clear default text
        _format_state_with_subscript(row[0].paragraphs[0], from_state)

        # Format "To" state with subscript
        to_state = _translate_state(trans.get("to_state", ""))
        row[1].text = ""  # Clear default text
        _format_state_with_subscript(row[1].paragraphs[0], to_state)

        # Format numeric values with 3 significant figures
        energy = trans.get('energy_eV', '')
        row[2].text = _format_significant_figures(energy) if energy != '' else ''

        wavelength = trans.get('wavelength_nm', '')
        row[3].text = _format_significant_figures(wavelength) if wavelength != '' else ''

        fosc = trans.get('oscillator_strength', '')
        row[4].text = _format_significant_figures(fosc) if fosc != '' else ''


def _add_rate_table(doc: Document, title: str, entries: Dict[str, Any]) -> None:
    if not entries:
        return
    doc.add_heading(title, level=2)
    table = doc.add_table(rows=1, cols=7)
    table.style = "Light Grid Accent 1"
    headers = ["Transition", "Rate (s⁻¹)", "Temperature (K)", "Δ0-0 (cm⁻¹)", "SOC (cm⁻¹)", "FC (%)", "HT (%)"]
    for idx, text in enumerate(headers):
        cell = table.rows[0].cells[idx]
        cell.text = text
    _style_header_row(table.rows[0])

    def _fmt_sci(val) -> str:
        try:
            return f"{float(val):.3e}"
        except Exception:
            return str(val) if val is not None else ""

    def _format_soc(rec: Dict[str, Any]) -> str:
        re_part = rec.get("soc_re_cm1")
        im_part = rec.get("soc_im_cm1")
        if re_part is None and im_part is None:
            return ""
        if re_part is not None and im_part is not None:
            return f"{_fmt_sci(re_part)}+i{_fmt_sci(im_part)}"
        if re_part is not None:
            return _fmt_sci(re_part)
        return f"i{_fmt_sci(im_part)}"

    def _add_row(label: str, rec: Dict[str, Any]) -> None:
        row = table.add_row().cells
        row[0].text = label
        rate_val = rec.get("rate_s1") or rec.get("total_rate_s1") or rec.get("rate")
        row[1].text = _fmt_sci(rate_val) if rate_val is not None else ""
        row[2].text = str(rec.get("temperature_K", ""))
        row[3].text = str(rec.get("delta_E_cm1", ""))
        row[4].text = _format_soc(rec)
        row[5].text = str(rec.get("fc_percent", ""))
        row[6].text = str(rec.get("ht_percent", ""))

    for name, record in sorted(entries.items()):
        ms_comps = record.get("ms_components", {})
        if ms_comps:
            ms_rec_list = []
            for ms_key, ms_rec in sorted(ms_comps.items()):
                pretty = {
                    "ms_0": "Ms=0",
                    "ms_p1": "Ms=+1",
                    "ms_m1": "Ms=-1",
                }.get(ms_key, ms_key.replace("_", "="))
                _add_row(f"{name} ({pretty})", ms_rec)
                ms_rec_list.append(ms_rec)
            # Add total if available
            total_rate = record.get("total_rate_s1") or record.get("rate_s1")
            if total_rate is not None:
                total_rec = {"total_rate_s1": total_rate}
                if ms_rec_list:
                    sample = ms_rec_list[0]
                    for key in ["temperature_K", "delta_E_cm1"]:
                        if sample.get(key) is not None:
                            total_rec[key] = sample.get(key)
                _add_row(f"{name} (total)", total_rec)
        else:
            _add_row(name, record)


def _extract_frontier_orbitals(orbital_data: Optional[Dict[str, Any]]) -> list[tuple[str, str, str]]:
    """Extract LUMO+3 to HOMO-3 orbital energies from orbital data.

    Returns:
        List of tuples (Index, Energy, Orbital) where Orbital is empty string.
    """
    if not orbital_data or "orbital_list" not in orbital_data:
        return []

    orbital_list = orbital_data["orbital_list"]

    # Find HOMO index (last occupied orbital)
    homo_idx = None
    for orbital in orbital_list:
        if orbital.get("occupancy", 0) > 1e-3:
            homo_idx = orbital.get("index")

    if homo_idx is None:
        return []

    # Define target orbitals relative to HOMO (LUMO is HOMO+1)
    target_orbitals = [
        ("LUMO+3", homo_idx + 4),
        ("LUMO+2", homo_idx + 3),
        ("LUMO+1", homo_idx + 2),
        ("LUMO", homo_idx + 1),
        ("HOMO", homo_idx),
        ("HOMO-1", homo_idx - 1),
        ("HOMO-2", homo_idx - 2),
        ("HOMO-3", homo_idx - 3),
    ]

    # Build index to orbital mapping
    idx_to_orbital = {orb["index"]: orb for orb in orbital_list}

    # Extract energies
    rows = []
    for label, idx in target_orbitals:
        if idx in idx_to_orbital:
            energy_ev = idx_to_orbital[idx].get("energy_eV")
            if energy_ev is not None:
                energy_str = f"{energy_ev:.4f}"
            else:
                energy_str = "n/a"
        else:
            energy_str = "n/a"
        rows.append((label, energy_str, ""))

    return rows


def _add_frontier_orbital_table(doc: Document, orbital_data: Optional[Dict[str, Any]]) -> None:
    """Add table showing frontier orbital energies (LUMO+3 to HOMO-3)."""
    rows = _extract_frontier_orbitals(orbital_data)
    if not rows:
        return

    doc.add_heading("Frontier Orbital Energies", level=2)
    table = doc.add_table(rows=1, cols=3)
    table.style = "Light Grid Accent 1"

    # Header row
    headers = ["Index", "Energy (eV)", "Orbital"]
    for idx, text in enumerate(headers):
        cell = table.rows[0].cells[idx]
        cell.text = text
    _style_header_row(table.rows[0])

    # Data rows
    for label, energy, orbital in rows:
        row_cells = table.add_row().cells
        row_cells[0].text = label
        row_cells[1].text = energy
        row_cells[2].text = orbital


def _create_energy_level_plot(data: Dict[str, Any], output_path: Path) -> Optional[Path]:
    """Create energy level diagram for S and T states.

    Args:
        data: DELFIN data dictionary
        output_path: Path to save the PNG plot

    Returns:
        Path to saved plot or None if creation failed
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
    except ImportError:
        logger.warning("matplotlib not available; skipping energy level plot")
        return None

    # Extract S and T states with their energies
    s_states = []
    t_states = []

    # Ground state S0
    gs = data.get("ground_state_S0", {}) or {}
    if gs:
        opt = gs.get("optimization", {}) or {}
        hartree = opt.get("hartree")
        if hartree is not None:
            s_states.append(("S0", float(hartree)))

    # Excited states
    excited = data.get("excited_states", {}) or {}
    for state_name, entry in excited.items():
        if not (state_name.startswith("S") or state_name.startswith("T")):
            continue
        opt = entry.get("optimization", {}) or {}
        hartree = opt.get("hartree")
        if hartree is not None:
            if state_name.startswith("S"):
                s_states.append((state_name, float(hartree)))
            elif state_name.startswith("T"):
                t_states.append((state_name, float(hartree)))

    if not s_states and not t_states:
        logger.warning("No S or T states found for energy level plot")
        return None

    # Sort by energy
    s_states.sort(key=lambda x: x[1])
    t_states.sort(key=lambda x: x[1])

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Define x positions for the two lanes
    s_x = 1.0
    t_x = 2.0
    line_width = 0.3

    # Plot S states
    for state_name, energy in s_states:
        ax.plot([s_x - line_width, s_x + line_width], [energy, energy], 'b-', linewidth=2)
        # Format label with subscript
        label = state_name.replace("S", "$S_").replace("0", "{0}") \
                          .replace("1", "{1}").replace("2", "{2}") \
                          .replace("3", "{3}").replace("4", "{4}") \
                          .replace("5", "{5}").replace("6", "{6}") + "$"
        ax.text(s_x + line_width + 0.05, energy, label, va='center', fontsize=10)

    # Plot T states
    for state_name, energy in t_states:
        ax.plot([t_x - line_width, t_x + line_width], [energy, energy], 'r-', linewidth=2)
        # Format label with subscript
        label = state_name.replace("T", "$T_").replace("1", "{1}") \
                          .replace("2", "{2}").replace("3", "{3}") \
                          .replace("4", "{4}").replace("5", "{5}") \
                          .replace("6", "{6}") + "$"
        ax.text(t_x + line_width + 0.05, energy, label, va='center', fontsize=10)

    # Styling
    ax.set_xlim(0.3, 3.2)
    ax.set_ylabel("Energy (Eh)", fontsize=12)
    ax.set_xticks([s_x, t_x])
    ax.set_xticklabels(["Singlet States", "Triplet States"], fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.title("Energy Level Diagram", fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save plot
    try:
        plt.savefig(str(output_path), dpi=300, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Energy level plot saved to {output_path}")
        return output_path
    except Exception as exc:
        logger.error(f"Failed to save energy level plot: {exc}")
        plt.close(fig)
        return None


def generate_combined_docx_report(
    project_dir: Path,
    json_path: Path,
    output_docx: Path,
    assets: ReportAssets | None = None,
) -> Optional[Path]:
    """
    Build DELFIN.docx using the collected JSON data and available plots.

    Args:
        project_dir: Workspace directory (used for relative paths)
        json_path: Path to DELFIN_Data.json
        output_docx: Destination DOCX path
        assets: Optional plot paths to embed

    Returns:
        Path to generated DOCX or None if generation failed.
    """
    if not DOCX_AVAILABLE:
        logger.error("python-docx not installed; cannot build DELFIN.docx")
        return None

    if assets is None:
        assets = ReportAssets()

    if assets.smiles_png is None:
        assets.smiles_png = _generate_smiles_image(project_dir)

    if not json_path.exists():
        logger.error("JSON data file not found: %s", json_path)
        return None

    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to read %s: %s", json_path, exc)
        return None

    doc = Document()

    # Title and metadata
    meta = data.get("metadata", {}) or {}
    molecule_name = meta.get("name") or meta.get("NAME") or project_dir.name
    heading = doc.add_heading(f"DELFIN Report – {molecule_name}", level=1)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER

    summary_text, color_boxes = _build_summary_text(data, project_dir)
    if summary_text:
        summary_para = doc.add_paragraph(summary_text)
        summary_para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    # Color chips for key transitions (render before structure)
    if color_boxes:
        chip_para = doc.add_paragraph()
        chip_para.add_run("Transition colors: ")
        for label, wl, rgb in color_boxes:
            chip_path = _render_color_chip(rgb, project_dir, label)
            if chip_path and chip_path.exists():
                run = chip_para.add_run()
                run.add_picture(str(chip_path), width=Inches(0.25))
            else:
                run_box = chip_para.add_run("[]")
                run_box.font.color.rgb = RGBColor(*rgb)
                run_box.font.bold = True
            chip_para.add_run(f" {label} ({wl:.0f} nm)  ")

    # SMILES picture (if available) near the top
    if assets.smiles_png and assets.smiles_png.exists():
        doc.add_heading("Structure from SMILES (input.txt)", level=2)
        doc.add_picture(str(assets.smiles_png), width=Inches(2.8))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

    meta_rows = [
        ("Functional", meta.get("functional", "")),
        ("Basis set", meta.get("basis_set", "")),
        ("Auxiliary basis", meta.get("auxiliary_basis", "")),
        ("RI method", meta.get("ri_method", "")),
        ("Dispersion correction", meta.get("dispersion_correction", "")),
        ("Solvent model", meta.get("implicit_solvation", "")),
        ("Solvent", meta.get("solvent", "")),
        ("Charge", meta.get("charge", "")),
        ("Multiplicity", meta.get("multiplicity", "")),
    ]
    _add_key_value_table(doc, "Calculation Setup", meta_rows)

    # Add frontier orbital energies table for S0
    gs = data.get("ground_state_S0", {}) or {}
    gs_orbitals = gs.get("orbitals") if gs else None
    _add_frontier_orbital_table(doc, gs_orbitals)

    # Collect state summaries for one consolidated table
    state_rows: Dict[str, Dict[str, Any]] = {}

    if gs:
        gs_entry = dict(gs)
        gs_entry["_type"] = "Ground"
        state_rows["S0"] = gs_entry
        # Absorption transitions from S0_TDDFT
        s0_abs = (gs.get("tddft_absorption") or {}).get("transitions", [])
        _add_transition_table(doc, "Absorption transitions (S0)", s0_abs)

    excited = data.get("excited_states", {}) or {}
    for name, entry in excited.items():
        entry = dict(entry)
        entry["_type"] = "Excited"
        state_rows[name] = entry
        transitions = (entry.get("tddft_from_geometry") or {}).get("transitions", [])
        if transitions:
            _add_transition_table(doc, f"Vertical transitions from {name}", transitions)

    for name, entry in (data.get("oxidized_states", {}) or {}).items():
        row = dict(entry)
        row["_type"] = "Oxidized"
        state_rows[name] = row

    for name, entry in (data.get("reduced_states", {}) or {}).items():
        row = dict(entry)
        row["_type"] = "Reduced"
        state_rows[name] = row

    _add_state_table(doc, "Energetics overview", state_rows)

    # Add energy level diagram
    _add_plot_if_exists(doc, "Energy Level Diagram", assets.energy_level_png)

    # Rates
    _add_rate_table(doc, "Intersystem crossing", data.get("intersystem_crossing", {}) or {})
    _add_rate_table(doc, "Internal conversion", data.get("internal_conversion", {}) or {})

    # Plots
    _add_plot_if_exists(doc, "AFP spectrum", assets.afp_png)
    if assets.uv_vis_pngs:
        for state, png_path in sorted(assets.uv_vis_pngs.items()):
            label = {
                "S0": "Absorption spectrum (S0)",
                "S1": "Fluorescence spectrum (S1)",
                "T1": "Phosphorescence spectrum (T1)",
            }.get(state, f"Spectrum ({state})")
            # Use heading with subscript for spectrum labels
            if label:
                _add_heading_with_subscript(doc, label, level=2)
                doc.add_picture(str(png_path), width=Inches(6.5))
                doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    # IR spectrum
    if assets.ir_png:
        _add_heading_with_subscript(doc, "IR spectrum (S0)", level=2)
        doc.add_picture(str(assets.ir_png), width=Inches(6.5))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

    output_docx.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(output_docx))
    logger.info("DELFIN.docx written to %s", output_docx)
    return output_docx
