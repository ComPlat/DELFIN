# delfin_reports.py
# Main DELFIN report generation functions

from collections import defaultdict
from typing import Optional
from pathlib import Path

from ..common.banners import build_standard_banner
from ..utils import (
    search_transition_metals,
    select_rel_and_aux,
    get_git_commit_info,
)
from ..esd_results import ESDSummary, ISCResult


def generate_summary_report_DELFIN(charge, multiplicity, solvent, E_ox, E_ox_2, E_ox_3,
                                   E_red, E_red_2, E_red_3, E_00_t1, E_00_s1,
                                   metals, metal_basisset, NAME, main_basisset,
                                   config, duration, E_ref, esd_summary: Optional[ESDSummary] = None,
                                   output_dir: Optional[Path] = None):
    import logging

    # If output_dir is provided, use it; otherwise use current directory
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)

    xyz = ""
    try:
        xyz_path = output_dir / 'initial.xyz'
        with open(xyz_path, 'r', encoding='utf-8') as xyz_file:
            xyz_lines = xyz_file.readlines()
            xyz = "".join(xyz_lines[2:])
    except FileNotFoundError:
        logging.warning(
            "File 'initial.xyz' not found; DELFIN.txt will omit coordinates."
        )
        xyz = "[initial.xyz missing]\n"

    # --- helpers for pretty printing of numbers --------------------------------
    def fmt_num(x):
        if x is None:
            return None
        try:
            v = float(x)
            return f" {v:.3f}" if v >= 0 else f"{v:.3f}"
        except Exception:
            s = str(x).strip()
            return f" {s}" if s and not s.startswith("-") else s

    def fmt_ev(x):
        if x is None:
            return None
        try:
            v = float(x)
            return f" {v:.3f} eV" if v >= 0 else f"{v:.3f} eV"
        except Exception:
            s = str(x).strip()
            return f" {s} eV" if s and not s.startswith("-") else f"{s} eV"

    def fmt_hartree(x):
        if x is None:
            return "n/a"
        try:
            return f"{float(x):.6f} Eh"
        except Exception:
            s = str(x).strip()
            return s or "n/a"

    def fmt_rate(rate):
        if rate is None:
            return "n/a"
        try:
            return f"{float(rate):.6e} s^-1"
        except Exception:
            return f"{rate} s^-1"

    # --- timing ----------------------------------------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    # --- derived properties ----------------------------------------------------
    E_red_star_s1 = (E_red + E_00_s1) if (E_red is not None and E_00_s1 is not None) else None
    E_red_star_t1 = (E_red + E_00_t1) if (E_red is not None and E_00_t1 is not None) else None
    E_ox_star_s1  = (E_ox  - E_00_s1) if (E_ox  is not None and E_00_s1 is not None) else None
    E_ox_star_t1  = (E_ox  - E_00_t1) if (E_ox  is not None and E_00_t1 is not None) else None
    ZFS = None  # placeholder

    # ---- method tokens (relativity + aux via utils policy) --------------------
    # If 'metals' argument is empty, re-detect from input for robustness
    if not metals:
        input_file = str(config.get("input_file", "input.txt")).strip() or "input.txt"
        try:
            metals = search_transition_metals(input_file)
        except Exception:
            metals = []
    rel_token, aux_jk_token, _use_rel = select_rel_and_aux(metals, config)

    # solvent helper (argument OR config fallback)
    def _solv_name():
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return s

    def implicit_token():
        model = str(config.get('implicit_solvation_model', '')).strip()
        if not model:
            return ""
        s = _solv_name()
        return f"{model}({s})" if s else model

    # Frequency method line
    method_freq_line = (
        f"Method freq: {config['functional']} {rel_token} {main_basisset} "
        f"{config.get('disp_corr','')} {config.get('ri_jkx','')} {aux_jk_token} {implicit_token()} "
        f"{config.get('geom_opt','OPT')} FREQ PAL{config.get('PAL','')} MAXCORE({config.get('maxcore','')})"
    ).replace("  ", " ").strip()

    # ESD method block (shown if ESD module is enabled)
    from delfin.esd_module import parse_esd_config
    esd_enabled, _, _, _ = parse_esd_config(config)
    if esd_enabled:
        esd_modus = str(config.get('ESD_modus', 'TDDFT')).strip().lower()

        if esd_modus == 'deltascf':
            # deltaSCF mode: show deltaSCF-specific keywords (no relativity, no RI-SOMF, no RI-J)
            # Note: deltaSCF uses exact same keywords as initial calculation but adds deltaSCF
            method_esd_block = (
                f"Method ESD (deltaSCF): {config['functional']} {main_basisset} {config.get('disp_corr','')} "
                f"{config.get('ri_jkx','')} {aux_jk_token} {implicit_token()} "
                f"OPT FREQ deltaSCF PAL{config['PAL']} MAXCORE({config['maxcore']})"
            ).replace("  ", " ").strip()
        else:
            # TDDFT mode: show TDDFT-specific keywords
            method_esd_block = (
                f"Method ESD (TDDFT): {config['functional']} {rel_token} {main_basisset} "
                f"{config.get('disp_corr','')} {config.get('ri_jkx','')} {aux_jk_token} {implicit_token()} "
                f"{config.get('ri_soc','')} PAL{config['PAL']} "
                f"NROOTS {config.get('ESD_nroots', config.get('NROOTS', 15))} "
                f"DOSOC {config.get('DOSOC', 'FALSE')} "
                f"TDA {config.get('ESD_TDA', config.get('TDA', 'FALSE'))} "
                f"MAXCORE({config['maxcore']})"
            ).replace("  ", " ").strip()

        # Add metals line if present
        if metals and metal_basisset:
            method_esd_block += f"\n                       {', '.join(metals)} {metal_basisset}"
    else:
        method_esd_block = ""

    # ---- blocks formatting -----------------------------------------------------
    def format_block(d):
        items = [(k, d[k]) for k in d if d[k] is not None and d[k] != ""]
        if not items:
            return ""
        width = max(len(k) for k, _ in items)
        return "\n".join(f"{k:<{width}} = {val}" for k, val in items)

    calculated_properties = {
        "E_00 (S1)": fmt_ev(E_00_s1),
        "E_00 (T1)": fmt_ev(E_00_t1),
        "E_red": fmt_num(E_red),
        "E_red_2": fmt_num(E_red_2),
        "E_red_3": fmt_num(E_red_3),
        "E_ox": fmt_num(E_ox),
        "E_ox_2": fmt_num(E_ox_2),
        "E_ox_3": fmt_num(E_ox_3),
        "*E_red (S1)": fmt_num(E_red_star_s1),
        "*E_red (T1)": fmt_num(E_red_star_t1),
        "*E_ox (S1)": fmt_num(E_ox_star_s1),
        "*E_ox (T1)": fmt_num(E_ox_star_t1),
        "ZFS": fmt_num(ZFS),
        "E_ref": (fmt_num(E_ref) if E_ref is not None else "need to be referenced!"),
    }
    calculated_block = format_block(calculated_properties)

    exp_pairs = [
        ("E_00", str(config.get("E_00_exp", "")).strip()),
        ("E_red", str(config.get("E_red_exp", "")).strip()),
        ("E_red_2", str(config.get("E_red_2_exp", "")).strip()),
        ("E_red_3", str(config.get("E_red_3_exp", "")).strip()),
        ("E_ox", str(config.get("E_ox_exp", "")).strip()),
        ("E_ox_2", str(config.get("E_ox_2_exp", "")).strip()),
        ("E_ox_3", str(config.get("E_ox_3_exp", "")).strip()),
        ("*E_red", str(config.get("*E_red_exp", "")).strip()),
        ("*E_ox", str(config.get("*E_ox_exp", "")).strip()),
    ]
    experimental_properties = {k: v for k, v in exp_pairs if v}
    experimental_block = format_block(experimental_properties)

    literature_reference = str(config.get('Literature_reference', '')).strip()
    smiles_info = str(config.get('SMILES', '')).strip()

    cfg_eref_str = str(config.get('E_ref', '')).strip()
    if cfg_eref_str == "":
        header_scale = "V vs. Fc+/Fc"
    else:
        try:
            cfg_eref_val = float(cfg_eref_str.replace(",", "."))
            if abs(cfg_eref_val - 4.345) <= 0.01:
                header_scale = "V vs. SCE"
            else:
                header_scale = "User defined"
        except ValueError:
            header_scale = "User defined"
    calc_properties_header = f"Calculated properties ({header_scale}):"

    # names
    if isinstance(NAME, (list, tuple, set)):
        name_str = ", ".join(map(str, NAME))
    else:
        name_str = str(NAME) if NAME is not None else ""

    # ---- assemble middle sections -------------------------------------------
    sections = []
    sections.append(f"{calc_properties_header}\n{calculated_block}")
    if experimental_block:
        sections.append(f"Experimental properties ({config.get('reference_CV', 'Unknown reference electrode')}):\n{experimental_block}")
    if literature_reference:
        sections.append(f"Literature References:\n(1): {literature_reference}")
    if smiles_info:
        sections.append(f"Informations:\nSMILES: {smiles_info}")
    if esd_summary and esd_summary.has_data:
        esd_lines: list[str] = []
        if esd_summary.states:
            esd_lines.append("Final single point energies (Hartree):")
            for state, record in sorted(esd_summary.states.items()):
                esd_lines.append(f"  {state} = {fmt_hartree(record.fspe)}")
            # Add Zero Point Energies if available
            zpe_data = [(state, record.zpe) for state, record in sorted(esd_summary.states.items()) if record.zpe is not None]
            if zpe_data:
                esd_lines.append("")  # Empty line
                esd_lines.append("Zero Point Energies (Hartree):")
                for state, zpe in zpe_data:
                    esd_lines.append(f"  {state} = {fmt_hartree(zpe)}")
            # Add Gibbs Free Energies if available
            gibbs_data = [(state, record.gibbs) for state, record in sorted(esd_summary.states.items()) if record.gibbs is not None]
            if gibbs_data:
                esd_lines.append("")  # Empty line
                esd_lines.append("Gibbs Free Energies (Hartree):")
                for state, gibbs in gibbs_data:
                    esd_lines.append(f"  {state} = {fmt_hartree(gibbs)}")
        if esd_summary.isc:
            if esd_lines:  # Add empty line before ISC section if there's already content
                esd_lines.append("")
            esd_lines.append("ISC rate constants (s^-1):")
            isc_grouped: dict[str, list[tuple[str, ISCResult]]] = defaultdict(list)

            # Group transitions by base (without Ms suffix) to allow summed output
            for transition, record in esd_summary.isc.items():
                base_transition = transition.split("(Ms=", 1)[0] if "(Ms=" in transition else transition
                isc_grouped[base_transition].append((transition, record))

            def _format_isc_line(label: str, record) -> str:
                extras: list[str] = []
                if record.temperature is not None:
                    extras.append(f"T={record.temperature:.2f} K")
                if record.delta_cm1 is not None:
                    extras.append(f"Δ0-0={record.delta_cm1:.2f} cm^-1")
                if record.soc is not None:
                    re_part, im_part = record.soc
                    if re_part is not None or im_part is not None:
                        soc_str = f"SOC={re_part:.6e}" if re_part is not None else "SOC=n/a"
                        if im_part is not None:
                            soc_str += f"+i{im_part:.6e}"
                        extras.append(soc_str)
                if record.fc_percent is not None:
                    extras.append(f"FC={record.fc_percent:.2f}%")
                if record.ht_percent is not None:
                    extras.append(f"HT={record.ht_percent:.2f}%")
                detail = f" ({', '.join(extras)})" if extras else ""
                return f"  {label} = {fmt_rate(record.rate)}{detail}"

            for base_transition in sorted(isc_grouped):
                records = sorted(isc_grouped[base_transition])
                for transition, record in records:
                    esd_lines.append(_format_isc_line(transition, record))

                # If multiple Ms components exist, append summed rate
                if len(records) > 1:
                    rates = [rec.rate for _, rec in records if rec.rate is not None]
                    total_rate = sum(rates) if rates else None
                    esd_lines.append(
                        f"  {base_transition} (total) = {fmt_rate(total_rate)}"
                    )
        if esd_summary.ic:
            if esd_lines:  # Add empty line before IC section if there's already content
                esd_lines.append("")
            esd_lines.append("IC rate constants (s^-1):")
            for transition, record in sorted(esd_summary.ic.items()):
                extras: list[str] = []
                if record.temperature is not None:
                    extras.append(f"T={record.temperature:.2f} K")
                if record.delta_cm1 is not None:
                    extras.append(f"Δ0-0={record.delta_cm1:.2f} cm^-1")
                detail = f" ({', '.join(extras)})" if extras else ""
                esd_lines.append(f"  {transition} = {fmt_rate(record.rate)}{detail}")
        if getattr(esd_summary, "fluor", None):
            fluor_items = [(k, v) for k, v in sorted(esd_summary.fluor.items()) if v and v.rate is not None]
            if fluor_items:
                if esd_lines:
                    esd_lines.append("")
                esd_lines.append("Fluorescence rate constants (s^-1):")
                for transition, record in fluor_items:
                    extras: list[str] = []
                    if record.temperature is not None:
                        extras.append(f"T={record.temperature:.2f} K")
                    if record.delta_cm1 is not None:
                        extras.append(f"Δ0-0={record.delta_cm1:.2f} cm^-1")
                    detail = f" ({', '.join(extras)})" if extras else ""
                    esd_lines.append(f"  {transition} = {fmt_rate(record.rate)}{detail}")
        if getattr(esd_summary, "phosp", None):
            phosp_items = [(k, v) for k, v in sorted(esd_summary.phosp.items()) if v and v.rate_mean is not None]
            if phosp_items:
                if esd_lines:
                    esd_lines.append("")
                esd_lines.append("Phosphorescence rate constants (s^-1):")
                for transition, record in phosp_items:
                    extras: list[str] = []
                    if record.temperature is not None:
                        extras.append(f"T={record.temperature:.2f} K")
                    if record.delta_cm1 is not None:
                        extras.append(f"Δ0-0={record.delta_cm1:.2f} cm^-1")
                    detail = f" ({', '.join(extras)})" if extras else ""
                    # Print k1/k2/k3 lines if present, then arithmetic mean
                    sub = [r for r in record.sublevel_rates if r is not None]
                    for i, r in enumerate(sub, start=1):
                        esd_lines.append(f"  {transition} (k{i}) = {fmt_rate(r)}{detail}")
                    esd_lines.append(f"  {transition} (mean) = {fmt_rate(record.rate_mean)}{detail}")
        if esd_lines:
            sections.append("ESD:\n" + "\n".join(esd_lines))
    middle = "\n\n".join(sections)

    # ---- write file ----------------------------------------------------------
    banner = build_standard_banner(header_indent=4, info_indent=4)

    # Get git commit info for reproducibility tracking
    git_commit = get_git_commit_info()
    if git_commit:
        # Remove -dirty suffix for GitHub link
        clean_commit = git_commit.replace("-dirty", "")
        git_info_line = f"Git commit: {git_commit} (https://github.com/ComPlat/DELFIN/commit/{clean_commit})"
    else:
        git_info_line = "Git commit: unknown (not in git repository)"

    delfin_output_path = output_dir / 'DELFIN.txt'
    with open(delfin_output_path, 'w', encoding='utf-8') as file:
        multiplicity_display = config.get('_multiplicity_display')
        if not multiplicity_display:
            multiplicity_display = str(multiplicity)
        # Format metals line for frequency method
        metals_line = f"                       {', '.join(metals)} {metal_basisset}" if metals and metal_basisset else ""

        # Build output with proper formatting
        output_parts = [
            banner,
            "",
            git_info_line,
            "",
            f"Compound name (NAME): {name_str}",
            "",
            method_freq_line,
        ]

        if metals_line:
            output_parts.append(metals_line)

        if method_esd_block:
            output_parts.append("")
            output_parts.append(method_esd_block)

        output_parts.extend([
            "",
            f"used Method: {config['method']}",
            f"Charge:        {charge}",
            f"Multiplicity:  {multiplicity_display}",
            "",
            f"Coordinates:\n{xyz}",
            "",
            middle,
            "",
            f"TOTAL RUN TIME: {duration_format}",
        ])

        file.write("\n".join(output_parts))
