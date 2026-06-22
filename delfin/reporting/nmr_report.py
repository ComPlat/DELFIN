"""
Generate NMR spectrum plots with molecule visualization and atom assignment.

Creates a two-panel figure:
  - Left:  1H NMR spectrum (Lorentzian-broadened) with labeled peaks
  - Right: 2D molecule structure with numbered H atoms, color-coded to match peaks
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import to_hex
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from delfin.common.logging import get_logger
from delfin.nmr_spectrum import NMRResult, NMRShielding

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Spectrum helpers
# ---------------------------------------------------------------------------

def lorentzian_broadening(
    shifts: List[float],
    intensities: Optional[List[float]] = None,
    ppm_range: Tuple[float, float] = (-1.0, 14.0),
    num_points: int = 4000,
    fwhm: float = 0.02,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply Lorentzian broadening to a list of chemical shifts."""
    ppm = np.linspace(ppm_range[0], ppm_range[1], num_points)
    spectrum = np.zeros(num_points)
    gamma = fwhm / 2.0

    if intensities is None:
        intensities = [1.0] * len(shifts)

    for delta, intensity in zip(shifts, intensities):
        spectrum += intensity * (gamma**2) / ((ppm - delta)**2 + gamma**2)

    return ppm, spectrum


def group_equivalent_hydrogens(
    shieldings: List[NMRShielding],
    tolerance_ppm: float = 0.15,
) -> List[List[NMRShielding]]:
    """Group H atoms with nearly identical chemical shifts (likely equivalent)."""
    sorted_s = sorted(shieldings, key=lambda s: s.chemical_shift)
    groups: List[List[NMRShielding]] = []
    current: List[NMRShielding] = []

    for s in sorted_s:
        if not current or abs(s.chemical_shift - current[-1].chemical_shift) <= tolerance_ppm:
            current.append(s)
        else:
            groups.append(current)
            current = [s]
    if current:
        groups.append(current)

    return groups


# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------

def _get_colors(n: int) -> List[str]:
    """Return *n* distinct colors as hex strings."""
    if not MATPLOTLIB_AVAILABLE:
        return ["#000000"] * n
    cmap = plt.cm.get_cmap("tab20", max(n, 1))
    return [to_hex(cmap(i)) for i in range(n)]


# ---------------------------------------------------------------------------
# Molecule 2D image (RDKit)
# ---------------------------------------------------------------------------

def _mol_from_coords(result: NMRResult) -> Optional["Chem.Mol"]:
    """Build an RDKit mol from parsed Cartesian coordinates."""
    if not RDKIT_AVAILABLE:
        return None

    from rdkit.Chem import rdDetermineBonds

    rwmol = Chem.RWMol()
    conf = Chem.Conformer(len(result.atoms))
    for atom in result.atoms:
        idx = rwmol.AddAtom(Chem.Atom(atom.element))
        conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
    rwmol.AddConformer(conf, assignId=True)

    # Let RDKit determine bonds from 3D coordinates
    mol = rwmol.GetMol()
    try:
        rdDetermineBonds.DetermineBonds(mol)
    except Exception as exc:
        logger.warning("rdDetermineBonds failed: %s — falling back to xyz2mol", exc)
        return None

    return mol


def _draw_mol_with_highlights(
    mol: "Chem.Mol",
    h_groups: List[List[NMRShielding]],
    colors: List[str],
    img_size: Tuple[int, int] = (900, 900),
) -> Optional[np.ndarray]:
    """Render 2D structure with highlighted & numbered H atoms."""
    if not RDKIT_AVAILABLE or not MATPLOTLIB_AVAILABLE:
        return None

    from rdkit.Chem import Draw
    from PIL import Image
    import io

    # Compute 2D coords
    AllChem.Compute2DCoords(mol)

    # Build highlight maps
    atom_colors: Dict[int, tuple] = {}
    atom_radii: Dict[int, float] = {}
    highlight_atoms: List[int] = []

    for gi, group in enumerate(h_groups):
        color_rgb = tuple(int(colors[gi].lstrip("#")[j:j+2], 16) / 255.0 for j in (0, 2, 4))
        for s in group:
            aidx = s.orca_index
            if aidx < mol.GetNumAtoms():
                highlight_atoms.append(aidx)
                atom_colors[aidx] = color_rgb
                atom_radii[aidx] = 0.4

    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = True
    opts.annotationFontScale = 0.7
    opts.addStereoAnnotation = False

    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=atom_colors,
        highlightAtomRadii=atom_radii,
    )
    drawer.FinishDrawing()

    png_bytes = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_bytes))
    return np.array(img)


# ---------------------------------------------------------------------------
# Main plot function
# ---------------------------------------------------------------------------

def create_nmr_report(
    result: NMRResult,
    output_png: Path,
    ppm_range: Tuple[float, float] = (-1.0, 14.0),
    fwhm: float = 0.02,
    equiv_tol: float = 0.15,
    dpi: int = 300,
    title: Optional[str] = None,
) -> None:
    """
    Create a two-panel NMR report figure.

    Left panel:  broadened 1H NMR spectrum with labeled peaks.
    Right panel: 2D molecular structure with color-coded H atoms.
    """
    if not MATPLOTLIB_AVAILABLE:
        logger.error("matplotlib is required for NMR plotting")
        return

    h_data = result.h_shieldings
    if not h_data:
        logger.error("No hydrogen shieldings found — nothing to plot")
        return

    # Group equivalent H
    groups = group_equivalent_hydrogens(h_data, tolerance_ppm=equiv_tol)
    colors = _get_colors(len(groups))

    # Build broadened spectrum
    shifts_all = [s.chemical_shift for s in h_data]
    ppm, spectrum = lorentzian_broadening(shifts_all, ppm_range=ppm_range, fwhm=fwhm)

    # Normalize
    if spectrum.max() > 0:
        spectrum = spectrum / spectrum.max()

    # --- Try to build molecule image ---
    mol = _mol_from_coords(result)
    mol_img = None
    if mol is not None:
        mol_img = _draw_mol_with_highlights(mol, groups, colors)

    has_mol = mol_img is not None

    if has_mol:
        fig, (ax_spec, ax_mol) = plt.subplots(
            1, 2, figsize=(18, 7), gridspec_kw={"width_ratios": [2, 1]}
        )
    else:
        fig, ax_spec = plt.subplots(figsize=(12, 6))

    # ---- Spectrum panel ----
    ax_spec.plot(ppm, spectrum, "k-", linewidth=0.8)
    ax_spec.fill_between(ppm, spectrum, alpha=0.05, color="black")

    # Color-coded sticks + labels for each group
    for gi, group in enumerate(groups):
        avg_shift = np.mean([s.chemical_shift for s in group])
        n_h = len(group)
        label_indices = ",".join(str(s.orca_index) for s in group)
        label = f"H{label_indices} ({n_h}H)"

        ax_spec.axvline(avg_shift, color=colors[gi], alpha=0.6, linewidth=1.2)

    ax_spec.set_xlabel("Chemical Shift (ppm)", fontsize=11)
    ax_spec.set_ylabel("Relative Intensity", fontsize=11)
    ax_spec.set_xlim(ppm_range[1], ppm_range[0])  # NMR convention: high ppm left
    ax_spec.set_ylim(-0.02, 1.15)
    ax_spec.set_title(title or "Calculated 1H NMR Spectrum", fontsize=13)

    # ---- Molecule panel ----
    if has_mol:
        ax_mol.imshow(mol_img)
        ax_mol.set_axis_off()
        ax_mol.set_title("Molecular Structure (atom indices)", fontsize=11)

    # ---- Legend ----
    legend_patches = []
    for gi, group in enumerate(groups):
        avg_shift = np.mean([s.chemical_shift for s in group])
        indices = ",".join(str(s.orca_index) for s in group)
        heavy = result.h_bonded_to(group[0].orca_index)
        heavy_label = f"-{heavy.element}{heavy.index}" if heavy else ""
        legend_patches.append(mpatches.Patch(
            color=colors[gi],
            label=f"H{indices}{heavy_label}: {avg_shift:.2f} ppm ({len(group)}H)",
        ))

    ax_spec.legend(
        handles=legend_patches,
        fontsize=5.5,
        loc="upper left",
        ncol=2,
        framealpha=0.7,
    )

    fig.tight_layout()
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("NMR report saved to %s", output_png)


# ---------------------------------------------------------------------------
# Assignment table (text)
# ---------------------------------------------------------------------------

def print_assignment_table(result: NMRResult, equiv_tol: float = 0.15) -> str:
    """Return a formatted text table of H assignments."""
    h_data = result.h_shieldings
    groups = group_equivalent_hydrogens(h_data, tolerance_ppm=equiv_tol)

    lines = [
        f"{'Group':>6}  {'H indices':<25}  {'delta (ppm)':>12}  {'nH':>4}  {'Bonded to':<12}",
        "-" * 75,
    ]
    for gi, group in enumerate(groups, 1):
        avg_shift = np.mean([s.chemical_shift for s in group])
        indices = ", ".join(str(s.orca_index) for s in group)
        heavy = result.h_bonded_to(group[0].orca_index)
        heavy_label = f"{heavy.element}{heavy.index}" if heavy else "?"
        lines.append(
            f"{gi:>6}  {indices:<25}  {avg_shift:>12.2f}  {len(group):>4}  {heavy_label:<12}"
        )

    return "\n".join(lines)
