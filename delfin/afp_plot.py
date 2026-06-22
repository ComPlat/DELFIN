"""
AFP (Absorption-Fluorescence-Phosphorescence) combined spectrum plot.

Creates a single plot showing the main absorption, fluorescence, and phosphorescence transitions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, List
import re
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from delfin.common.logging import get_logger
from delfin.uv_vis_spectrum import parse_absorption_spectrum, Transition

logger = get_logger(__name__)


def _translate_state(orca_state: str) -> str:
    """Translate ORCA notation (0-1A/1-3A) to S0/S1/T1 for matching/labels."""
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


def find_first_allowed_singlet(transitions: List[Transition], min_fosc: float = 0.001) -> Optional[Transition]:
    """
    Find the first allowed singlet transition with significant oscillator strength.

    Args:
        transitions: List of transitions
        min_fosc: Minimum oscillator strength threshold

    Returns:
        First singlet transition with fosc >= min_fosc, or None
    """
    for trans in transitions:
        # Check if it's a singlet transition (ends with -1A)
        if trans.to_state.endswith('-1A') and trans.fosc >= min_fosc:
            return trans

    # If no allowed singlet found, try with lower threshold
    logger.warning(f"No singlet transition found with fosc >= {min_fosc}, trying lower threshold")
    for trans in transitions:
        if trans.to_state.endswith('-1A') and trans.fosc >= min_fosc / 10:
            return trans

    # Fallback: return first singlet even if fosc is very small
    for trans in transitions:
        if trans.to_state.endswith('-1A'):
            logger.warning(f"Using singlet transition with very low fosc: {trans.fosc:.6f}")
            return trans

    return None


def wavelength_to_rgb(wavelength: float) -> tuple:
    """
    Convert wavelength in nm to RGB color.

    Args:
        wavelength: Wavelength in nanometers

    Returns:
        RGB tuple (r, g, b) with values 0-1
    """
    if wavelength < 380:
        # UV - show as violet
        return (0.5, 0.0, 0.5)
    elif wavelength < 440:
        # Violet to Blue
        r = -(wavelength - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif wavelength < 490:
        # Blue to Cyan
        r = 0.0
        g = (wavelength - 440) / (490 - 440)
        b = 1.0
    elif wavelength < 510:
        # Cyan to Green
        r = 0.0
        g = 1.0
        b = -(wavelength - 510) / (510 - 490)
    elif wavelength < 580:
        # Green to Yellow
        r = (wavelength - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif wavelength < 645:
        # Yellow to Red
        r = 1.0
        g = -(wavelength - 645) / (645 - 580)
        b = 0.0
    elif wavelength < 781:
        # Red
        r = 1.0
        g = 0.0
        b = 0.0
    else:
        # IR - show as dark red
        return (0.5, 0.0, 0.0)

    # Intensity factor for edges of visible spectrum
    if wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif wavelength > 700:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
    else:
        factor = 1.0

    return (r * factor, g * factor, b * factor)


def _pick_by_states(transitions: List[Transition], from_state: str, to_state: str) -> Optional[Transition]:
    """Pick the first transition matching translated from/to states (keeps file order)."""
    for trans in transitions:
        if _translate_state(trans.from_state) == from_state and _translate_state(trans.to_state) == to_state:
            return trans
    return None


def gaussian(x: np.ndarray, center: float, amplitude: float, fwhm: float = 20.0) -> np.ndarray:
    """
    Generate Gaussian curve.

    Args:
        x: X values (wavelengths)
        center: Peak center (nm)
        amplitude: Peak amplitude
        fwhm: Full width at half maximum (nm)

    Returns:
        Y values (intensities)
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))


def create_afp_plot(
    workspace_dir: Path,
    output_png: Optional[Path] = None,
    fwhm: float = 50.0,
    dpi: int = 300
) -> Optional[Path]:
    """
    Create combined Absorption-Fluorescence-Phosphorescence plot.

    Args:
        workspace_dir: Working directory containing ESD/ folder or S0.out, S1.out, T1.out
        output_png: Output PNG path (default: workspace_dir/AFP_spectrum.png)
        fwhm: Full width at half maximum for Gaussian broadening (nm)
        dpi: Resolution of output image

    Returns:
        Path to saved PNG or None if failed
    """
    if not MATPLOTLIB_AVAILABLE:
        logger.error("matplotlib not available, cannot create AFP plot")
        return None

    # Look for output files in ESD/ or current directory
    esd_dir = workspace_dir / "ESD"
    search_dirs = [esd_dir, workspace_dir] if esd_dir.exists() else [workspace_dir]

    # Find S0, S1, T1 output files
    s0_file = None
    s1_file = None
    t1_file = None

    for search_dir in search_dirs:
        if s0_file is None:
            s0_candidates = [search_dir / "S0.out", search_dir / "S0_TDDFT.out"]
            for candidate in s0_candidates:
                if candidate.exists():
                    s0_file = candidate
                    break

        if s1_file is None:
            s1_candidates = [
                search_dir / "S1_first_TDDFT.out",  # hybrid1 mode
                search_dir / "S1_TDDFT.out",        # deltaSCF mode
                search_dir / "S1.out"               # TDDFT mode
            ]
            for candidate in s1_candidates:
                if candidate.exists():
                    s1_file = candidate
                    break

        if t1_file is None and (search_dir / "T1.out").exists():
            t1_file = search_dir / "T1.out"

    # Parse first transition from each file
    transitions_data = []

    if s0_file:
        logger.info(f"Parsing S0 absorption from {s0_file.name}")
        s0_transitions = parse_absorption_spectrum(s0_file)
        if s0_transitions:
            # Pick the S0 → S1 transition (use file order)
            s0_singlet = _pick_by_states(s0_transitions, "S0", "S1")
            if s0_singlet:
                transitions_data.append({
                    'transition': s0_singlet,
                    'label': f"S0 → S1 (Absorption)",
                    'color': 'blue',
                    'state': 'S0'
                })
                logger.info(f"  S0 → S1 (absorption): {s0_singlet.wavelength_nm:.1f} nm (fosc = {s0_singlet.fosc:.4f})")
            else:
                logger.warning("No S0 → S1 transition found in S0.out")

    if s1_file:
        logger.info(f"Parsing S1 fluorescence from {s1_file.name}")
        s1_transitions = parse_absorption_spectrum(s1_file)
        if s1_transitions:
            # Pick the S0 → S1 transition at S1 geometry (mirror image approximation for emission)
            s1_singlet = _pick_by_states(s1_transitions, "S0", "S1")
            if s1_singlet:
                transitions_data.append({
                    'transition': s1_singlet,
                    'label': f"S1 → S0 (Fluorescence)",
                    'color': 'green',
                    'state': 'S1'
                })
                logger.info(f"  S1 → S0 (emission): {s1_singlet.wavelength_nm:.1f} nm (fosc = {s1_singlet.fosc:.4f})")
            else:
                logger.warning("No S0 → S1 transition found in S1.out")

    if t1_file:
        logger.info(f"Parsing T1 phosphorescence from {t1_file.name}")
        t1_transitions = parse_absorption_spectrum(t1_file)
        if t1_transitions:
            # Pick the S0 → T1 transition at T1 geometry (mirror image approximation for emission)
            t1_to_singlet = _pick_by_states(t1_transitions, "S0", "T1")
            if t1_to_singlet:
                transitions_data.append({
                    'transition': t1_to_singlet,
                    'label': f"T1 → S0 (Phosphorescence)",
                    'color': 'red',
                    'state': 'T1'
                })
                logger.info(f"  T1 → S0 (emission): {t1_to_singlet.wavelength_nm:.1f} nm (fosc = {t1_to_singlet.fosc:.6f})")
            else:
                logger.warning("No S0 → T1 transition found in T1.out")

    if not transitions_data:
        logger.error("No transitions found in S0.out, S1.out, or T1.out")
        return None

    # Determine wavelength range
    wavelengths = [t['transition'].wavelength_nm for t in transitions_data]
    min_wl = min(wavelengths)
    max_wl = max(wavelengths)
    padding = max(50, (max_wl - min_wl) * 0.2)
    wl_range = (max(200, min_wl - padding), max_wl + padding)

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 7))

    # Generate wavelength array
    x = np.linspace(wl_range[0], wl_range[1], 3000)

    # Plot each transition with Gaussian curve
    max_intensity = 0
    # Different vertical offsets for labels to avoid overlap
    label_offsets = [15, 45, 75]  # Different heights for S0, S1, T1

    # Find maximum fosc for normalization
    max_fosc = max(data['transition'].fosc for data in transitions_data)
    if max_fosc <= 0:
        logger.warning("All fosc values are zero or negative, using uniform intensity")
        max_fosc = 1.0

    for idx, data in enumerate(transitions_data):
        trans = data['transition']
        color = data['color']
        label = data['label']

        # Normalize intensity based on fosc (highest fosc = 1.0)
        intensity = trans.fosc / max_fosc if max_fosc > 0 else 1.0

        # Generate Gaussian curve
        y = gaussian(x, trans.wavelength_nm, intensity, fwhm)
        max_intensity = max(max_intensity, y.max())

        # Plot Gaussian curve
        ax.plot(x, y, color=color, linewidth=2.5, label=label, alpha=0.8)

        # Add stick line
        ax.axvline(trans.wavelength_nm, color=color, linestyle='--',
                   linewidth=1.5, alpha=0.5)

        # Add peak label with wavelength (use different offset for each)
        peak_y = gaussian(np.array([trans.wavelength_nm]), trans.wavelength_nm, intensity, fwhm)[0]
        y_offset = label_offsets[idx] if idx < len(label_offsets) else 15

        # Get the actual color of this wavelength
        wl_color = wavelength_to_rgb(trans.wavelength_nm)

        # Add wavelength label with white background, including fosc
        ax.annotate(f'{trans.wavelength_nm:.1f} nm\n(f={trans.fosc:.4f})',
                    xy=(trans.wavelength_nm, peak_y),
                    xytext=(-15, y_offset),
                    textcoords='offset points',
                    ha='center',
                    fontsize=9,
                    fontweight='bold',
                    color=color,
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                             edgecolor=color, alpha=0.8))

        # Add small colored box showing the actual wavelength color next to the label
        ax.annotate('  ',  # Empty space with colored background
                    xy=(trans.wavelength_nm, peak_y),
                    xytext=(25, y_offset),
                    textcoords='offset points',
                    ha='center',
                    fontsize=10,
                    color='black',
                    bbox=dict(boxstyle='square,pad=0.35', facecolor=wl_color,
                             edgecolor='black', linewidth=0.8))

    # Formatting
    ax.set_xlabel('Wavelength (nm)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Intensity (normalized)', fontsize=14, fontweight='bold')
    ax.set_title('Absorption, Fluorescence, and Phosphorescence Spectrum',
                 fontsize=16, fontweight='bold', pad=20)
    ax.set_xlim(wl_range)
    ax.set_ylim(0, 1.6)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='upper right', fontsize=12, framealpha=0.9)

    # Add info text
    info_text = f"FWHM = {fwhm} nm"
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Determine output path
    if output_png is None:
        output_png = workspace_dir / "AFP_spectrum.png"

    fig.savefig(output_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    logger.info(f"AFP spectrum plot saved to {output_png}")
    return output_png


def generate_afp_report(
    workspace_dir: Path,
    output_png: Optional[Path] = None,
    fwhm: float = 50.0
) -> Optional[Path]:
    """
    Generate AFP (Absorption-Fluorescence-Phosphorescence) spectrum plot.

    Args:
        workspace_dir: Working directory
        output_png: Optional output PNG path
        fwhm: Full width at half maximum for Gaussian broadening (nm)

    Returns:
        Path to generated PNG or None
    """
    logger.info("Generating AFP spectrum plot")
    return create_afp_plot(workspace_dir, output_png, fwhm=fwhm)
