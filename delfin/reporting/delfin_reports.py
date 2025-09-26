# delfin_reports.py
# Main DELFIN report generation functions

from decimal import Decimal, ROUND_DOWN
from typing import Optional
import os, re

from ..utils import (
    search_transition_metals,
    select_rel_and_aux,
)
from ..parser import extract_last_uhf_deviation, extract_last_J3


def generate_summary_report_DELFIN(charge, multiplicity, solvent, E_ox, E_ox_2, E_ox_3,
                                   E_red, E_red_2, E_red_3, E_00_t1, E_00_s1,
                                   metals, metal_basisset, NAME, main_basisset,
                                   config, duration, E_ref):
    import logging
    xyz = ""
    try:
        with open('initial.xyz', 'r', encoding='utf-8') as xyz_file:
            xyz_lines = xyz_file.readlines()
            xyz = "".join(xyz_lines[2:])
    except FileNotFoundError:
        logging.error("File 'initial.xyz' not found.")
        return

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

    # TDDFT method block (only if E_00 requested)
    if str(config.get("E_00", "")).lower() == "yes":
        method_tddft_block = (
            f"Method TDDFT: {config['functional']} {rel_token} {main_basisset} "
            f"{implicit_token()} {config.get('ri_soc','')} "
            f"PAL{config['PAL']} NROOTS {config['NROOTS']} DOSOC {config['DOSOC']} TDA {config['TDA']} MAXCORE({config['maxcore']})\n"
            f"        {', '.join(metals)} {metal_basisset if metal_basisset else ''}"
        ).replace("  ", " ").strip()
    else:
        method_tddft_block = ""

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
    middle = "\n\n".join(sections)

    # ---- write file ----------------------------------------------------------
    with open('DELFIN.txt', 'w', encoding='utf-8') as file:
        file.write(f"""
                          ******************
                          *     DELFIN     *
                          ******************

    #############################################################
    #                           -***-                           #
    #                          ComPlat                          #
    #                        M. Hartmann                        #
    #           Karlsruhe Institute of Technology (KIT)         #
    #   Automates ORCA 6.1.0, xTB 6.7.1 and CREST 3.0.2 runs    #
    #                       Version 1.0.2                       #
    #                           -***-                           #
    #############################################################

Compound name (NAME): {name_str}

{method_freq_line}
        {', '.join(metals)} {metal_basisset if metal_basisset else ''}
{method_tddft_block}

used Method: {config['method']}
Charge:        {charge}
Multiplicity:  {multiplicity}

Coordinates:
{xyz}

{middle}

TOTAL RUN TIME: {duration_format}
""")