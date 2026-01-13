from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from delfin.common.logging import get_logger
from delfin.reporting.delfin_collector import save_esd_data_json
from delfin.reporting.delfin_docx_report import ReportAssets, generate_combined_docx_report

logger = get_logger(__name__)


def _find_output_file(workspace: Path, candidates: list[str]) -> Optional[Path]:
    for name in candidates:
        for folder in (workspace / "ESD", workspace):
            candidate = folder / name
            if candidate.exists():
                return candidate
    return None


def _expected_uvvis_png(output_file: Path, spectrum_type: str, state_name: str) -> Path:
    return output_file.parent / f"{spectrum_type.replace(' ', '_')}_{state_name}.png"


def _generate_uvvis_outputs(workspace: Path, fwhm: float) -> Dict[str, Path]:
    from delfin.reporting.uv_vis_report import generate_uv_vis_word_report

    state_to_file = {
        "S0": ["S0_TDDFT.out", "S0.out"],
        "S1": ["S1_TDDFT.out", "S1.out"],
        "T1": ["T1.out"],
    }
    pngs: Dict[str, Path] = {}

    for state, names in state_to_file.items():
        output_file = _find_output_file(workspace, names)
        if not output_file:
            logger.warning("No output file found for %s spectrum", state)
            continue

        logger.info("Generating UV/Vis report for %s using %s", state, output_file)
        try:
            generate_uv_vis_word_report(output_file, output_docx=None, fwhm=fwhm)
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to generate UV/Vis report for %s: %s", state, exc, exc_info=True)
            continue

        state_name = output_file.stem.replace("_TDDFT", "")
        if state_name == "S0":
            spectrum_type = "Absorption Spectrum"
        elif state_name.startswith("S"):
            spectrum_type = "Fluorescence Spectrum"
        elif state_name.startswith("T"):
            spectrum_type = "Phosphorescence Spectrum"
        else:
            spectrum_type = "Absorption Spectrum"

        png_path = _expected_uvvis_png(output_file, spectrum_type, state_name)
        if png_path.exists():
            pngs[state_name] = png_path
        else:
            logger.warning("Expected spectrum PNG not found for %s at %s", state_name, png_path)

    return pngs


def _generate_ir_outputs(workspace: Path) -> Optional[Path]:
    from delfin.ir_spectrum import parse_ir_spectrum
    from delfin.reporting.ir_report import create_ir_spectrum_plot, generate_ir_report

    # Try S0.out first; if it exists but has no IR section, fall back to initial.out
    candidates = ["S0.out", "initial.out"]
    output_file: Optional[Path] = None
    modes = None

    for name in candidates:
        candidate = _find_output_file(workspace, [name])
        if not candidate:
            continue
        logger.info("Parsing IR spectrum from %s", candidate)
        modes = parse_ir_spectrum(candidate)
        if modes:
            output_file = candidate
            break
        logger.warning("No IR modes found in %s", candidate)

    if output_file is None or not modes:
        logger.warning("No IR-capable output (S0.out or initial.out) yielded IR modes")
        return None

    output_docx = output_file.parent / f"IR_{output_file.stem}.docx"
    output_png = output_docx.parent / f"{output_docx.stem}.png"

    create_ir_spectrum_plot(modes, output_png=output_png)
    generate_ir_report(
        output_file=output_docx,
        modes=modes,
        spectrum_png=output_png,
        source_file=output_file,
        molecule_name=output_file.stem,
    )
    if output_png.exists():
        return output_png
    return None


def run_docx_report_mode(
    workspace_root: Path,
    config,
    afp_fwhm: float,
    json_output_path: Optional[Path] = None,
) -> int:
    """Build DELFIN.docx by collecting JSON data and embedding generated plots."""
    try:
        json_path = save_esd_data_json(workspace_root, json_output_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to build DELFIN_Data.json: %s", exc, exc_info=True)
        return 1

    assets = ReportAssets()

    # AFP plot
    try:
        from delfin.afp_plot import generate_afp_report

        assets.afp_png = generate_afp_report(workspace_root, fwhm=afp_fwhm)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate AFP spectrum: %s", exc, exc_info=True)

    # UV/Vis reports for S0, S1, T1
    try:
        assets.uv_vis_pngs = _generate_uvvis_outputs(workspace_root, fwhm=afp_fwhm)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate UV/Vis spectra: %s", exc, exc_info=True)

    # IR report for S0
    try:
        assets.ir_png = _generate_ir_outputs(workspace_root)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate IR spectrum: %s", exc, exc_info=True)

    # Energy level diagram
    try:
        import json
        from delfin.reporting.delfin_docx_report import _create_energy_level_plot

        data = json.loads(json_path.read_text(encoding="utf-8"))
        energy_plot_path = workspace_root / "Energy_Level_Diagram.png"
        assets.energy_level_png = _create_energy_level_plot(data, energy_plot_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate energy level diagram: %s", exc, exc_info=True)

    # Vertical excitation diagram
    try:
        import json
        from delfin.reporting.delfin_docx_report import _create_vertical_excitation_plot

        data = json.loads(json_path.read_text(encoding="utf-8"))
        vertical_plot_path = workspace_root / "Vertical_Excitation_Energies.png"
        assets.vertical_excitation_png = _create_vertical_excitation_plot(data, vertical_plot_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate vertical excitation diagram: %s", exc, exc_info=True)

    # Correlation diagram
    try:
        import json
        from delfin.reporting.delfin_docx_report import _create_correlation_plot

        data = json.loads(json_path.read_text(encoding="utf-8"))
        correlation_plot_path = workspace_root / "Correlation_Diagram.png"
        assets.correlation_png = _create_correlation_plot(data, correlation_plot_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate correlation diagram: %s", exc, exc_info=True)

    # Dipole moment visualization
    try:
        import json
        from delfin.reporting.delfin_docx_report import _create_dipole_moment_plot

        data = json.loads(json_path.read_text(encoding="utf-8"))
        dipole_plot_path = workspace_root / "Dipole_Moment_Visualization.png"
        assets.dipole_moment_png = _create_dipole_moment_plot(workspace_root, data, dipole_plot_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate dipole moment visualization: %s", exc, exc_info=True)

    # Molecular orbital visualizations
    try:
        from delfin.reporting.delfin_docx_report import _create_mo_visualizations

        # Generate MOs for HOMO-3 to LUMO+3
        orbital_indices = [-3, -2, -1, 0, 1, 2, 3, 4]
        assets.mo_pngs = _create_mo_visualizations(workspace_root, orbital_indices)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate MO visualizations: %s", exc, exc_info=True)

    # Electrostatic potential (ESP) plots
    try:
        from delfin.reporting.esp_report import generate_esp_png, generate_esp_png_for_state

        esp_pngs: Dict[str, Path] = {}
        for state in ("S0", "S1", "T1"):
            png = generate_esp_png_for_state(workspace_root, state)
            if png:
                esp_pngs[state] = png
        assets.esp_pngs = esp_pngs or None
        if not assets.esp_pngs:
            assets.esp_png = generate_esp_png(workspace_root)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to generate ESP plot: %s", exc, exc_info=True)

    # Assemble DELFIN.docx
    output_docx = workspace_root / "DELFIN.docx"
    result = generate_combined_docx_report(workspace_root, json_path, output_docx, assets=assets)
    if result is None:
        return 1

    logger.info("DELFIN.docx generated at %s", result)
    print(f"DELFIN.docx written to: {result}")
    return 0
