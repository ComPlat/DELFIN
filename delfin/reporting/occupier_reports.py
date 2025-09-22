# occupier_reports.py
# OCCUPIER-specific report generation functions

from decimal import Decimal, ROUND_DOWN
from pathlib import Path
from typing import Optional
import os, re

from ..utils import (
    search_transition_metals,
    select_rel_and_aux,
)
from ..parser import extract_last_uhf_deviation, extract_last_J3


# This file will contain the OCCUPIER report generation functions
# Functions will be moved here from report.py

def generate_summary_report_OCCUPIER(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence):
    """Generate the OCCUPIER summary report using the revised selection logic.

    Selection logic (revised):
      1) Energy is the primary criterion; in tolerance mode include all candidates within epsilon.
      2) If there is only a single energy winner, apply the *clean override* and add significantly
         "cleaner" candidates (better spin quality) within a small Hartree window around the minimum.
      3) Tie breaks:
         • Prefer lower spin contamination ("quality"): non-BS uses the effective deviation,
           BS uses |Dev − expected BS pair count|.
         • No explicit bias on multiplicity/BS beyond quality.
      4) Apply the *clean bias* before the energy bias:
         • If energies are comparable, choose the cleaner solution.
         • If quality is comparable, choose the lower-energy solution.

    Additional CONTROL switches (optional, meaningful defaults):
      clean_override_window_h   (float, default 0.003): energy window for pulling in cleaner candidates
                                  when there is only one energy winner.
      clean_quality_improvement (float, default 0.05): minimum quality improvement required for the override.
      clean_quality_good        (float, default 0.05): absolute threshold regarded as "good" quality.
      clean_bias_window_h       (float, default 0.003): energy proximity used by the clean bias.
      quality_bias_window       (float, default 0.05): minimum quality difference for the clean bias to trigger.

    Further configurable knobs:
      dev_similarity            (float, default 0.15): similarity threshold when comparing deviations.
      bs_override_window_h      (float, default 0.002): AF override energy window.

    Existing switches that continue to apply:
      occupier_selection = tolerance | truncation | rounding
      occupier_precision
      occupier_epsilon
      energy_bias_window_h
      mismatch_bias_window
      dev_max
      approximate_spin_projection_APMethod
    """
    # ----------------------- imports & tiny helpers ----------------------------
    import os, re
    from typing import Optional
    from decimal import Decimal, ROUND_DOWN

    def truncate(x: float, d: int) -> float:
        q = Decimal(10) ** -d
        return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))

    # pretty-print with the active selection precision & mode
    def fmt_truncate(x: float, d: int) -> str:
        return f"{truncate(x, d):.{d}f}"

    # ----------------------- locate output files -------------------------------
    def file_for_index(idx: int) -> str:
        out_files = config.get('out_files')
        if isinstance(out_files, dict):
            return out_files.get(idx) or ("output.out" if idx == 1 else f"output{idx}.out")
        if isinstance(out_files, (list, tuple)) and 1 <= idx <= len(out_files):
            return out_files[idx-1] or ("output.out" if idx == 1 else f"output{idx}.out")
        return "output.out" if idx == 1 else f"output{idx}.out"

    # ----------------------- input method line (optional) ----------------------
    def _inp_candidates_for_index(idx: int) -> list[str]:
        cands: list[str] = []
        inp_files = config.get('inp_files')
        if isinstance(inp_files, dict):
            v = inp_files.get(idx)
            if v: cands.append(str(v))
        elif isinstance(inp_files, (list, tuple)) and 1 <= idx <= len(inp_files):
            v = inp_files[idx-1]
            if v: cands.append(str(v))

        ofile = file_for_index(idx)
        base = Path(ofile).stem
        cands += [base + ".inp"]
        if base.lower().startswith("output"):
            cands.append("input" + base[6:] + ".inp")
        cands.append(ofile.replace("output", "input").rsplit(".", 1)[0] + ".inp")
        if idx == 1: cands.append("input.inp")
        cands += [f"input{idx}.inp", f"input_{idx}.inp"]

        seen, uniq = set(), []
        for p in cands:
            p2 = p.strip()
            if p2 and p2 not in seen:
                uniq.append(p2); seen.add(p2)
        return uniq

    def _inp_for_index(idx: int) -> Optional[str]:
        for p in _inp_candidates_for_index(idx):
            if Path(p).exists():
                return p
        return None

    def _parse_method_from_inp(inp_path: str, pal: str, include_freq: bool) -> Optional[str]:
        try:
            with open(inp_path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    ls = line.lstrip()
                    if ls.startswith("!"):
                        tokens = ls[1:].strip().split()
                        drop = {"PMODEL", "MOREAD"}
                        tokens = [t for t in tokens if t.upper() not in drop]
                        if include_freq and not any(t.upper() == "FREQ" for t in tokens):
                            tokens.append("FREQ")
                        method = "Method: " + " ".join(tokens).strip()
                        if str(pal).strip():
                            method = (method + f" PAL{pal}").replace("  ", " ").strip()
                        return method
        except Exception:
            pass
        return None

    def _parse_metal_basis_from_inp(inp_path: str) -> Optional[str]:
        try:
            with open(inp_path, "r", encoding="utf-8", errors="replace") as fh:
                in_xyz = False
                for line in fh:
                    s = line.strip()
                    if s.startswith("*"):
                        in_xyz = not in_xyz
                        continue
                    if in_xyz:
                        m = re.search(r'NewGTO\s+"([^"]+)"', s)
                        if m:
                            return m.group(1).strip()
        except Exception:
            pass
        return None

    # ----------------------- CONTROL fallbacks for method line -----------------
    rel = str(config.get("relativity", "none")).strip().lower()
    rel_token_ctl = "ZORA" if rel == "zora" else ""
    aux_jk_token_ctl = str(config.get('aux_jk_rel' if rel == 'zora' else 'aux_jk', '')).strip()

    if rel == "zora":
        metal_basisset_rep_ctl = str(config.get("metal_basisset_rel", config.get("metal_orbital_basis_zora", ""))).strip()
    else:
        metal_basisset_rep_ctl = str(config.get("metal_basisset",     config.get("metal_orbital_basis",      ""))).strip()

    try:
        metals = search_transition_metals(str(config.get("input_file", "input.txt")).strip() or "input.txt")
    except Exception:
        metals = []

    def implicit_token():
        model = str(config.get('implicit_solvation_model', '')).strip()
        if not model: return ""
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return f"{model}({s})" if s else model

    def method_line_from_control(include_freq: bool) -> str:
        tokens = [
            "Method:", str(config.get('functional','')).strip(), rel_token_ctl,
            str(main_basisset).strip(),
            str(config.get('disp_corr','')).strip(),
            str(config.get('ri_jkx','')).strip(),
            aux_jk_token_ctl, implicit_token(),
            str(config.get('geom_opt_OCCUPIER','')).strip()
        ]
        if include_freq:
            tokens.append("FREQ")
        tokens.append(f"PAL{config.get('PAL','')}")
        return " ".join(t for t in tokens if t).replace("  ", " ").strip()

    # ----------------------- timing & labels -----------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""
    parity = "is_even" if is_even else "is_odd"
    use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    energy_label = "Final Gibbs free energy" if use_gibbs else "FINAL SINGLE POINT ENERGY"
    lowest_label = f"LOWEST {energy_label}:"

    # selection mode
    raw_sel = str(config.get('occupier_selection', 'tolerance')).lower()
    sel = raw_sel.split('|')[0].strip()
    ap = config.get('approximate_spin_projection_APMethod')
    ap_str = str(ap) if ap not in (None, 0, '0') else 'none'

    if sel in {'rounding','round','gerundet','runden'}:
        method = 'rounding'
    elif sel in {'tolerance','toleranz','toleranzband','epsilon'}:
        method = 'tolerance'
    else:
        method = 'truncation'

    # precision / epsilon
    try:
        prec = int(config.get('occupier_precision', config.get('occupier_rounded_value', 6)))
    except (TypeError, ValueError):
        prec = 6
    prec = max(0, min(10, prec))

    try:
        epsilon = float(config.get('occupier_epsilon', 10.0**(-prec)))
        if not (epsilon > 0): raise ValueError
    except (TypeError, ValueError):
        epsilon = 10.0**(-prec)

    # optional dev_max (file-local, not CONTROL-driven) – filters heavily contaminated entries early
    try:
        dev_max = float(config.get('dev_max')) if config.get('dev_max') is not None else None
    except (TypeError, ValueError):
        dev_max = None

    # ----------------------- Schwellen/Parameter -------------------------------
    DEV_TINY       = 1e-3
    DEV_SIMILARITY = float(config.get('dev_similarity', 0.15))
    DEV_GOOD_MARGIN= 0.30
    DEV_HIGH       = 0.50
    DEV_MATCH_WINDOW = 0.30

    # AF override (energy window to pull BS candidates back in when contamination is high)
    EPS_AF = float(config.get('bs_override_window_h', 0.002))

    # Energy-bias (wie bisher)
    E_BIAS_H = float(config.get('energy_bias_window_h', 0.002))
    MIS_BIAS = float(config.get('mismatch_bias_window', 0.05))

    # Clean-preference knobs
    CLEAN_OVERRIDE_H = float(config.get('clean_override_window_h', 0.003))
    CLEAN_Q_IMPROVE  = float(config.get('clean_quality_improvement', 0.05))
    CLEAN_Q_GOOD     = float(config.get('clean_quality_good', 0.05))
    CLEAN_BIAS_H     = float(config.get('clean_bias_window_h', 0.003))
    QUAL_BIAS_WIN    = float(config.get('quality_bias_window', 0.05))

    # ----------------------- caches & parsers ----------------------------------
    dev_cache: dict[int, Optional[float]] = {}
    j3_cache: dict[int, Optional[float]] = {}

    def get_dev(idx: int) -> Optional[float]:
        if idx not in dev_cache:
            try:
                dev_cache[idx] = extract_last_uhf_deviation(file_for_index(idx))
            except Exception:
                dev_cache[idx] = None
        return dev_cache[idx]

    def get_j3(idx: int, is_bs_flag: bool) -> Optional[float]:
        """Parse J(3) for printing only."""
        if not is_bs_flag:
            return None
        if idx not in j3_cache:
            try:
                j3_cache[idx] = extract_last_J3(file_for_index(idx))
            except Exception:
                j3_cache[idx] = None
        return j3_cache[idx]

    # sequence accessors
    def entry_record(idx: int) -> dict:
        return next((e for e in sequence if e["index"] == idx), {})

    def entry_is_bs(idx: int) -> bool:
        return bool(entry_record(idx).get("BS"))

    def entry_mult(idx: int) -> int:
        try:    return int(entry_record(idx).get("m", 0))
        except: return 0

    def bs_pair_count(idx: int) -> Optional[int]:
        lab = str(entry_record(idx).get("BS") or "").strip()
        if not lab or "," not in lab:
            return None
        try:
            _, n = lab.split(",", 1)
            return int(n.strip())
        except Exception:
            return None

    # deviation utilities
    def effective_dev(idx: int) -> float:
        dev = get_dev(idx)
        if dev is None:
            if (not entry_is_bs(idx)) and (entry_mult(idx) == 1):
                return 0.0  # closed-shell RKS
            return float("inf")
        return dev

    def bs_mismatch(idx: int) -> float:
        n = bs_pair_count(idx)
        if n is None:
            return float("inf")
        return abs(effective_dev(idx) - n)

    def is_pseudo_closed_shell(idx: int) -> bool:
        return entry_is_bs(idx) and (entry_mult(idx) == 1) and (effective_dev(idx) <= DEV_TINY)

    def within_dev_limit(idx: int) -> bool:
        if dev_max is None:
            return True
        d = get_dev(idx)
        return (d is None) or (d <= dev_max)

    # ----------------------- energy band ---------------------------------------
    valid_all = [(e["index"], f) for e, f in zip(sequence, fspe_values) if f is not None]
    valid = [pair for pair in valid_all if within_dev_limit(pair[0])] or valid_all
    energies_by_idx = {i: f for i, f in valid_all}

    def energy_key(val: float) -> float:
        if method == "rounding":   return round(val, prec)
        if method == "truncation": return truncate(val, prec)
        return val  # tolerance -> raw

    if method in ("rounding", "truncation"):
        processed = [(i, f, energy_key(f)) for i, f in valid]
        best_key = min(v for _, _, v in processed)
        cands = [(i, f) for i, f, v in processed if v == best_key]
        min_raw = min(f for _, f in valid)
    else:
        min_raw = min(f for _, f in valid)
        eps_eff = epsilon + 1e-12
        cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

    # --- AF override (unchanged logic, but parameters configurable) --------------
    if len(cands) == 1:
        best_idx, _best_val = cands[0]
        if effective_dev(best_idx) >= DEV_HIGH:
            extra = []
            for i, f in valid:
                if i == best_idx or not entry_is_bs(i):
                    continue
                d = effective_dev(i)
                if d >= DEV_HIGH and abs(d - effective_dev(best_idx)) <= DEV_SIMILARITY and (f - min_raw) <= EPS_AF + 1e-12:
                    extra.append((i, f))
            if extra:
                cands = cands + extra  # unify flow below

    # --- Clean override: pull in cleaner candidates that are close in energy ----
    if len(cands) == 1:
        best_idx, _ = cands[0]
        q_best = (abs(effective_dev(best_idx) - bs_pair_count(best_idx))
                  if entry_is_bs(best_idx) else effective_dev(best_idx))
        extra = []
        for i, f in valid:
            if i == best_idx:
                continue
            if (f - min_raw) <= CLEAN_OVERRIDE_H + 1e-12:
                q_i = (abs(effective_dev(i) - bs_pair_count(i))
                       if entry_is_bs(i) else effective_dev(i))
                if ((q_best - q_i) >= CLEAN_Q_IMPROVE) or (q_i <= CLEAN_Q_GOOD):
                    if (i, f) not in cands:
                        extra.append((i, f))
        if extra:
            cands = cands + extra

    # ----------------------- tie-breaking (single flow) ------------------------
    min_fspe_index: Optional[int] = None
    min_fspe_value: Optional[float] = None

    if len(cands) == 1:
        min_fspe_index, min_fspe_value = cands[0]
    elif cands:
        # reclassify pseudo-CS BS to non-BS
        raw_bs = [(i, f) for i, f in cands if entry_is_bs(i)]
        raw_nb = [(i, f) for i, f in cands if not entry_is_bs(i)]
        nb_cands = raw_nb + [(i, f) for i, f in raw_bs if is_pseudo_closed_shell(i)]
        bs_cands = [(i, f) for i, f in raw_bs if not is_pseudo_closed_shell(i)]

        # scores
        def score_nb(i: int):
            # prefer smaller deviation; prefer genuine non-BS over pseudo (flag 0 < 1)
            pseudo_flag = 1 if is_pseudo_closed_shell(i) else 0
            return (effective_dev(i), pseudo_flag, i)

        def score_bs(i: int):
            # deviation-prior then mismatch then index
            mis = bs_mismatch(i)
            return (mis > DEV_MATCH_WINDOW, mis, effective_dev(i), i)

        # regime check
        def _min_eff_dev(pairs):
            vals = [effective_dev(i) for i, _ in pairs]
            return min(vals) if vals else float("inf")

        min_dev_bs = _min_eff_dev(bs_cands)
        min_dev_nb = _min_eff_dev(nb_cands)

        if bs_cands and nb_cands:
            # A) non-BS clearly cleaner and not high-contamination → pick best non-BS candidate
            if (min_dev_nb + DEV_GOOD_MARGIN < min_dev_bs) and (min_dev_nb < DEV_HIGH):
                pick = min(nb_cands, key=lambda p: score_nb(p[0]))
            # B) similar in quality and both high → pick the BS candidate with the better mismatch score
            elif (abs(min_dev_nb - min_dev_bs) <= DEV_SIMILARITY) and (min_dev_nb > DEV_HIGH) and (min_dev_bs > DEV_HIGH):
                pick = min(bs_cands, key=lambda p: score_bs(p[0]))
            # C) Fallback → kombiniere, bevorzuge sauberer (non-BS first), dann Dev, dann Index
            else:
                def fallback_score(i: int):
                    return (1 if entry_is_bs(i) else 0, effective_dev(i), i)
                pick = min(cands, key=lambda p: fallback_score(p[0]))
        elif bs_cands:
            pick = min(bs_cands, key=lambda p: score_bs(p[0]))
        else:
            pick = min(nb_cands, key=lambda p: score_nb(p[0]))

        min_fspe_index, min_fspe_value = pick

    # ----------------------- clean-bias vor energy-bias ------------------------
    def _qual_metric(i: int) -> float:
        # Spinsauberkeit: kleiner ist besser
        return bs_mismatch(i) if entry_is_bs(i) else effective_dev(i)

    if (min_fspe_index is not None) and cands:
        pick_i = min_fspe_index
        pick_E = energies_by_idx.get(pick_i, float("inf"))
        pick_Q = _qual_metric(pick_i)

        # 1) Clean-Bias: wenn E nah, nimm deutlich saubereren
        for j, _ in cands:
            if j == pick_i:
                continue
            Ej = energies_by_idx.get(j, float("inf"))
            Qj = _qual_metric(j)
            energy_close = abs(Ej - pick_E) <= CLEAN_BIAS_H
            if energy_close and (pick_Q - Qj) >= QUAL_BIAS_WIN:
                pick_i, pick_E, pick_Q = j, Ej, Qj

        # 2) Energy bias (only when quality is practically identical)
        for j, _ in cands:
            if j == pick_i:
                continue
            Ej = energies_by_idx.get(j, float("inf"))
            Qj = _qual_metric(j)
            close_in_quality = abs(Qj - pick_Q) <= MIS_BIAS
            close_in_energy  = abs(Ej - pick_E) <= E_BIAS_H
            if close_in_quality and close_in_energy and Ej < pick_E:
                pick_i, pick_E, pick_Q = j, Ej, Qj

        min_fspe_index = pick_i
        min_fspe_value = pick_E

    # ----------------------- printing helpers ----------------------------------
    def fmt_display(x: Optional[float]) -> str:
        if x is None: return "N/A"
        return f"{x:.{prec}f}" if method == 'rounding' else fmt_truncate(x, prec)

    fspe_lines = ""
    for entry, fspe in zip(sequence, fspe_values):
        idx = entry["index"]
        multiplicity = entry["m"]
        bs = entry.get("BS", "")
        is_bs_flag = bool(bs)

        fspe_value = (fmt_display(fspe) + " (H)") if fspe is not None else "Not a valid value"
        dev_txt = fmt_display(get_dev(idx))
        mark = " <-- PREFERRED VALUE" if (min_fspe_index is not None and idx == min_fspe_index) else ""
        bs_line = f", BrokenSym {bs}" if bs else ""

        fspe_lines += f"{energy_label} ({idx})   = {fspe_value}{mark}\n"
        fspe_lines += f"multiplicity {multiplicity}{bs_line}\n"
        fspe_lines += f"Spin Contamination (⟨S²⟩ - S(S+1))   : {dev_txt}\n"
        if is_bs_flag:
            j3_val = get_j3(idx, True)
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm**-1"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(<S**2>HS-<S**2>BS) : {j3_str}\n"
        fspe_lines += "----------------------------------------------------------------\n"

    # ----------------------- method line / metals block ------------------------
    include_freq = (freq_flag == 'FREQ')
    pal = config.get('PAL', '')

    method_str = None
    metal_basis_from_inp = None
    prefer_idx = (min_fspe_index if min_fspe_index is not None else 1)
    inp_pref = _inp_for_index(prefer_idx) or _inp_for_index(1)

    if inp_pref:
        parsed = _parse_method_from_inp(inp_pref, pal, include_freq)
        if parsed: method_str = parsed
        metal_basis_from_inp = _parse_metal_basis_from_inp(inp_pref)

    if not method_str:
        method_str = method_line_from_control(include_freq)

    metal_basis_print = metal_basis_from_inp or metal_basisset_rep_ctl or ""

    lowest_str = f"{fmt_truncate(min_fspe_value, prec)} (H)" if min_fspe_value is not None else "No valid FSPE values found"

    # ----------------------- write report --------------------------------------
    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(f"""
                     *******************
                     *     OCCUPIER    *
                     *******************

      #####################################################
      #                     -***-                         #
      #                    ComPlat                        #
      #     Karlsruhe Institute of Technology (KIT)       #
      #         Automates ORCA 6.1.0 calculations         #
      #                 Version 1.0.0                     #
      #                     -***-                         #
      #####################################################

{method_str}
        {', '.join(metals)} {metal_basis_print}

Charge: {charge}
-------------
{fspe_lines}

TOTAL RUN TIME: {duration_format}

{lowest_label} {lowest_str}

(Selection: {method}, APmethod {ap_str}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})
(Energy-bias: energy_bias_window_h={E_BIAS_H}, mismatch_bias_window={MIS_BIAS} -> prefer lower E only when qualities are similar)
(Clean-bias: clean_bias_window_h={CLEAN_BIAS_H}, quality_bias_window={QUAL_BIAS_WIN} -> prefer cleaner when energies are close)
(Clean-override: clean_override_window_h={CLEAN_OVERRIDE_H}, clean_quality_improvement={CLEAN_Q_IMPROVE}, clean_quality_good={CLEAN_Q_GOOD})
(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})
(Electron number: {parity})
""")


def generate_summary_report_OCCUPIER_safe(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence):
    """
    Selection is driven by energy first; if tied, use spin-contamination only.
    J(3) is parsed/printed for BS entries but NEVER used in selection.

    Tie-break summary:
      - Prefer clearly smaller contamination.
      - Prior toward the expected BS pair count: minimize |Deviation - n_BS|.
      - BS singlet with ~0 contamination is treated as pseudo closed-shell (non-BS).
      - Final tie-breakers: smaller deviation, then lower index.
      - AF-override (no J): if a single best-energy candidate is highly contaminated,
        temporarily add nearby BS candidates within ΔE <= bs_override_window_h and
        with similarly high contamination; then re-run the same tie-break logic.
      - NEW: Energy-bias (optional): if two candidates are very close in quality (mismatch/Dev)
        and energy (windows configurable), prefer the lower energy.
        Keys: energy_bias_window_h (default 0.002), mismatch_bias_window (default 0.05).
    """
    # ----------------------- imports & tiny helpers ----------------------------
    import os, re
    from typing import Optional
    from decimal import Decimal, ROUND_DOWN

    def truncate(x: float, d: int) -> float:
        q = Decimal(10) ** -d
        return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))

    # pretty-print with the active selection precision & mode
    def fmt_truncate(x: float, d: int) -> str:
        return f"{truncate(x, d):.{d}f}"

    # ----------------------- locate output files -------------------------------
    def file_for_index(idx: int) -> str:
        out_files = config.get('out_files')
        if isinstance(out_files, dict):
            return out_files.get(idx) or ("output.out" if idx == 1 else f"output{idx}.out")
        if isinstance(out_files, (list, tuple)) and 1 <= idx <= len(out_files):
            return out_files[idx-1] or ("output.out" if idx == 1 else f"output{idx}.out")
        return "output.out" if idx == 1 else f"output{idx}.out"

    # ----------------------- input method line (optional) ----------------------
    def _inp_candidates_for_index(idx: int) -> list[str]:
        cands: list[str] = []
        inp_files = config.get('inp_files')
        if isinstance(inp_files, dict):
            v = inp_files.get(idx)
            if v: cands.append(str(v))
        elif isinstance(inp_files, (list, tuple)) and 1 <= idx <= len(inp_files):
            v = inp_files[idx-1]
            if v: cands.append(str(v))

        ofile = file_for_index(idx)
        base = Path(ofile).stem
        cands += [base + ".inp"]
        if base.lower().startswith("output"):
            cands.append("input" + base[6:] + ".inp")
        cands.append(ofile.replace("output", "input").rsplit(".", 1)[0] + ".inp")
        if idx == 1: cands.append("input.inp")
        cands += [f"input{idx}.inp", f"input_{idx}.inp"]

        seen, uniq = set(), []
        for p in cands:
            p2 = p.strip()
            if p2 and p2 not in seen:
                uniq.append(p2); seen.add(p2)
        return uniq

    def _inp_for_index(idx: int) -> Optional[str]:
        for p in _inp_candidates_for_index(idx):
            if Path(p).exists():
                return p
        return None

    def _parse_method_from_inp(inp_path: str, pal: str, include_freq: bool) -> Optional[str]:
        try:
            with open(inp_path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    ls = line.lstrip()
                    if ls.startswith("!"):
                        tokens = ls[1:].strip().split()
                        drop = {"PMODEL", "MOREAD"}
                        tokens = [t for t in tokens if t.upper() not in drop]
                        if include_freq and not any(t.upper() == "FREQ" for t in tokens):
                            tokens.append("FREQ")
                        method = "Method: " + " ".join(tokens).strip()
                        if str(pal).strip():
                            method = (method + f" PAL{pal}").replace("  ", " ").strip()
                        return method
        except Exception:
            pass
        return None

    def _parse_metal_basis_from_inp(inp_path: str) -> Optional[str]:
        try:
            with open(inp_path, "r", encoding="utf-8", errors="replace") as fh:
                in_xyz = False
                for line in fh:
                    s = line.strip()
                    if s.startswith("*"):
                        in_xyz = not in_xyz
                        continue
                    if in_xyz:
                        m = re.search(r'NewGTO\s+"([^"]+)"', s)
                        if m:
                            return m.group(1).strip()
        except Exception:
            pass
        return None

    # ----------------------- CONTROL fallbacks for method line -----------------
    rel = str(config.get("relativity", "none")).strip().lower()
    rel_token_ctl = "ZORA" if rel == "zora" else ""
    aux_jk_token_ctl = str(config.get('aux_jk_rel' if rel == 'zora' else 'aux_jk', '')).strip()

    if rel == "zora":
        metal_basisset_rep_ctl = str(config.get("metal_basisset_rel", config.get("metal_orbital_basis_zora", ""))).strip()
    else:
        metal_basisset_rep_ctl = str(config.get("metal_basisset",     config.get("metal_orbital_basis",      ""))).strip()

    try:
        metals = search_transition_metals(str(config.get("input_file", "input.txt")).strip() or "input.txt")
    except Exception:
        metals = []

    def implicit_token():
        model = str(config.get('implicit_solvation_model', '')).strip()
        if not model: return ""
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return f"{model}({s})" if s else model

    def method_line_from_control(include_freq: bool) -> str:
        tokens = [
            "Method:", str(config.get('functional','')).strip(), rel_token_ctl,
            str(main_basisset).strip(),
            str(config.get('disp_corr','')).strip(),
            str(config.get('ri_jkx','')).strip(),
            aux_jk_token_ctl, implicit_token(),
            str(config.get('geom_opt_OCCUPIER','')).strip()
        ]
        if include_freq:
            tokens.append("FREQ")
        tokens.append(f"PAL{config.get('PAL','')}")
        return " ".join(t for t in tokens if t).replace("  ", " ").strip()

    # ----------------------- timing & labels -----------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""
    parity = "is_even" if is_even else "is_odd"
    use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    energy_label = "Final Gibbs free energy" if use_gibbs else "FINAL SINGLE POINT ENERGY"
    lowest_label = f"LOWEST {energy_label}:"

    # selection mode
    raw_sel = str(config.get('occupier_selection', 'tolerance')).lower()
    sel = raw_sel.split('|')[0].strip()
    ap = config.get('approximate_spin_projection_APMethod')
    ap_str = str(ap) if ap not in (None, 0, '0') else 'none'

    if sel in {'rounding','round','gerundet','runden'}:
        method = 'rounding'
    elif sel in {'tolerance','toleranz','toleranzband','epsilon'}:
        method = 'tolerance'
    else:
        method = 'truncation'

    # precision / epsilon
    try:
        prec = int(config.get('occupier_precision', config.get('occupier_rounded_value', 6)))
    except (TypeError, ValueError):
        prec = 6
    prec = max(0, min(10, prec))

    try:
        epsilon = float(config.get('occupier_epsilon', 10.0**(-prec)))
        if not (epsilon > 0): raise ValueError
    except (TypeError, ValueError):
        epsilon = 10.0**(-prec)

    # optional dev_max (file-local, not CONTROL)
    try:
        dev_max = float(config.get('dev_max')) if config.get('dev_max') is not None else None
    except (TypeError, ValueError):
        dev_max = None
    #----------------------------------------------------------------------------
    # thresholds
    DEV_TINY         = 1e-3
    DEV_SIMILARITY   = 0.15
    DEV_GOOD_MARGIN  = 0.30
    DEV_HIGH         = 0.50
    DEV_MATCH_WINDOW = 0.30
    EPS_AF = float(config.get('bs_override_window_h', 0.002))

    # Energy-bias knobs (prefer lower E if "close")
    E_BIAS_H = float(config.get('energy_bias_window_h', 0.002))
    MIS_BIAS = float(config.get('mismatch_bias_window', 0.05))

    # ----------------------- caches & parsers ----------------------------------
    dev_cache: dict[int, Optional[float]] = {}
    j3_cache: dict[int, Optional[float]] = {}

    def get_dev(idx: int) -> Optional[float]:
        if idx not in dev_cache:
            try:
                dev_cache[idx] = extract_last_uhf_deviation(file_for_index(idx))
            except Exception:
                dev_cache[idx] = None
        return dev_cache[idx]

    def get_j3(idx: int, is_bs_flag: bool) -> Optional[float]:
        """Parse J(3) for printing only."""
        if not is_bs_flag:
            return None
        if idx not in j3_cache:
            try:
                j3_cache[idx] = extract_last_J3(file_for_index(idx))
            except Exception:
                j3_cache[idx] = None
        return j3_cache[idx]

    # sequence accessors
    def entry_record(idx: int) -> dict:
        return next((e for e in sequence if e["index"] == idx), {})

    def entry_is_bs(idx: int) -> bool:
        return bool(entry_record(idx).get("BS"))

    def entry_mult(idx: int) -> int:
        try:    return int(entry_record(idx).get("m", 0))
        except: return 0

    def bs_pair_count(idx: int) -> Optional[int]:
        lab = str(entry_record(idx).get("BS") or "").strip()
        if not lab or "," not in lab:
            return None
        try:
            _, n = lab.split(",", 1)
            return int(n.strip())
        except Exception:
            return None

    # deviation utilities
    def effective_dev(idx: int) -> float:
        dev = get_dev(idx)
        if dev is None:
            if (not entry_is_bs(idx)) and (entry_mult(idx) == 1):
                return 0.0  # closed-shell RKS
            return float("inf")
        return dev

    def bs_mismatch(idx: int) -> float:
        n = bs_pair_count(idx)
        if n is None:
            return float("inf")
        return abs(effective_dev(idx) - n)

    def is_pseudo_closed_shell(idx: int) -> bool:
        return entry_is_bs(idx) and (entry_mult(idx) == 1) and (effective_dev(idx) <= DEV_TINY)

    def within_dev_limit(idx: int) -> bool:
        if dev_max is None:
            return True
        d = get_dev(idx)
        return (d is None) or (d <= dev_max)

    # ----------------------- energy band ---------------------------------------
    valid_all = [(e["index"], f) for e, f in zip(sequence, fspe_values) if f is not None]
    valid = [pair for pair in valid_all if within_dev_limit(pair[0])] or valid_all
    energies_by_idx = {i: f for i, f in valid_all}

    def energy_key(val: float) -> float:
        if method == "rounding":   return round(val, prec)
        if method == "truncation": return truncate(val, prec)
        return val  # tolerance -> raw

    if method in ("rounding", "truncation"):
        processed = [(i, f, energy_key(f)) for i, f in valid]
        best_key = min(v for _, _, v in processed)
        cands = [(i, f) for i, f, v in processed if v == best_key]
        min_raw = min(f for _, f in valid)
    else:
        min_raw = min(f for _, f in valid)
        eps_eff = epsilon + 1e-12
        cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

    # --- AF-override simplified: extend candidate set, then reuse same logic ---
    if len(cands) == 1:
        best_idx, _best_val = cands[0]
        if effective_dev(best_idx) >= DEV_HIGH:
            extra = []
            for i, f in valid:
                if i == best_idx or not entry_is_bs(i):
                    continue
                d = effective_dev(i)
                if d >= DEV_HIGH and abs(d - effective_dev(best_idx)) <= DEV_SIMILARITY and (f - min_raw) <= EPS_AF + 1e-12:
                    extra.append((i, f))
            if extra:
                cands = cands + extra  # unify flow below

    # ----------------------- tie-breaking (single flow) ------------------------
    min_fspe_index: Optional[int] = None
    min_fspe_value: Optional[float] = None

    if len(cands) == 1:
        min_fspe_index, min_fspe_value = cands[0]
    elif cands:
        # reclassify pseudo-CS BS to non-BS
        raw_bs = [(i, f) for i, f in cands if entry_is_bs(i)]
        raw_nb = [(i, f) for i, f in cands if not entry_is_bs(i)]
        nb_cands = raw_nb + [(i, f) for i, f in raw_bs if is_pseudo_closed_shell(i)]
        bs_cands = [(i, f) for i, f in raw_bs if not is_pseudo_closed_shell(i)]

        # scores
        def score_nb(i: int):
            # prefer smaller deviation; prefer genuine non-BS over pseudo (flag 0 < 1)
            pseudo_flag = 1 if is_pseudo_closed_shell(i) else 0
            return (effective_dev(i), pseudo_flag, i)

        def score_bs(i: int):
            # deviation-prior then mismatch then index
            mis = bs_mismatch(i)
            return (mis > DEV_MATCH_WINDOW, mis, effective_dev(i), i)

        # regime check
        def _min_eff_dev(pairs):
            vals = [effective_dev(i) for i, _ in pairs]
            return min(vals) if vals else float("inf")

        min_dev_bs = _min_eff_dev(bs_cands)
        min_dev_nb = _min_eff_dev(nb_cands)

        if bs_cands and nb_cands:
            # A) non-BS clearly cleaner and not high → choose best non-BS
            if (min_dev_nb + DEV_GOOD_MARGIN < min_dev_bs) and (min_dev_nb < DEV_HIGH):
                pick = min(nb_cands, key=lambda p: score_nb(p[0]))
            # B) similar & both high → choose best BS by mismatch prior
            elif (abs(min_dev_nb - min_dev_bs) <= DEV_SIMILARITY) and (min_dev_nb > DEV_HIGH) and (min_dev_bs > DEV_HIGH):
                pick = min(bs_cands, key=lambda p: score_bs(p[0]))
            # C) fallback → combine and pick by class preference (non-BS first), then deviation, then index
            else:
                def fallback_score(i: int):
                    return (1 if entry_is_bs(i) else 0, effective_dev(i), i)
                pick = min(cands, key=lambda p: fallback_score(p[0]))
        elif bs_cands:
            pick = min(bs_cands, key=lambda p: score_bs(p[0]))
        else:
            pick = min(nb_cands, key=lambda p: score_nb(p[0]))

        min_fspe_index, min_fspe_value = pick

    # ----------------------- energy-bias: prefer lower E when "close" ----------
    def _qual_metric(i: int) -> float:
        # Use same notion as tie-break: BS -> mismatch to expected n; non-BS -> effective deviation
        return bs_mismatch(i) if entry_is_bs(i) else effective_dev(i)

    if (min_fspe_index is not None) and cands:
        pick_i = min_fspe_index
        pick_E = energies_by_idx.get(pick_i, float("inf"))
        pick_Q = _qual_metric(pick_i)
        for j, _ in cands:
            if j == pick_i:
                continue
            Ej = energies_by_idx.get(j, float("inf"))
            Qj = _qual_metric(j)
            close_in_quality = abs(Qj - pick_Q) <= MIS_BIAS
            close_in_energy  = abs(Ej - pick_E) <= E_BIAS_H
            if close_in_quality and close_in_energy and Ej < pick_E:
                pick_i, pick_E = j, Ej
        # apply potential energy-biased pick
        min_fspe_index = pick_i
        min_fspe_value = pick_E

    # ----------------------- printing helpers ----------------------------------
    def fmt_display(x: Optional[float]) -> str:
        if x is None: return "N/A"
        return f"{x:.{prec}f}" if method == 'rounding' else fmt_truncate(x, prec)

    fspe_lines = ""
    for entry, fspe in zip(sequence, fspe_values):
        idx = entry["index"]
        multiplicity = entry["m"]
        bs = entry.get("BS", "")
        is_bs_flag = bool(bs)

        fspe_value = (fmt_display(fspe) + " (H)") if fspe is not None else "Not a valid value"
        dev_txt = fmt_display(get_dev(idx))
        mark = " <-- PREFERRED VALUE" if (min_fspe_index is not None and idx == min_fspe_index) else ""
        bs_line = f", BrokenSym {bs}" if bs else ""

        fspe_lines += f"{energy_label} ({idx})   = {fspe_value}{mark}\n"
        fspe_lines += f"multiplicity {multiplicity}{bs_line}\n"
        fspe_lines += f"Spin Contamination (⟨S²⟩ - S(S+1))   : {dev_txt}\n"
        if is_bs_flag:
            j3_val = get_j3(idx, True)
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm**-1"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(<S**2>HS-<S**2>BS) : {j3_str}\n"
        fspe_lines += "----------------------------------------------------------------\n"

    # ----------------------- method line / metals block ------------------------
    include_freq = (freq_flag == 'FREQ')
    pal = config.get('PAL', '')

    method_str = None
    metal_basis_from_inp = None
    prefer_idx = (min_fspe_index if min_fspe_index is not None else 1)
    inp_pref = _inp_for_index(prefer_idx) or _inp_for_index(1)

    if inp_pref:
        parsed = _parse_method_from_inp(inp_pref, pal, include_freq)
        if parsed: method_str = parsed
        metal_basis_from_inp = _parse_metal_basis_from_inp(inp_pref)

    if not method_str:
        method_str = method_line_from_control(include_freq)

    metal_basis_print = metal_basis_from_inp or metal_basisset_rep_ctl or ""

    lowest_str = f"{fmt_truncate(min_fspe_value, prec)} (H)" if min_fspe_value is not None else "No valid FSPE values found"

    # ----------------------- write report --------------------------------------
    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(f"""
                     *******************
                     *     OCCUPIER    *
                     *******************

      #####################################################
      #                     -***-                         #
      #                    ComPlat                        #
      #     Karlsruhe Institute of Technology (KIT)       #
      #         Automates ORCA 6.1.0 calculations         #
      #                 Version 1.0.0                     #
      #                     -***-                         #
      #####################################################

{method_str}
        {', '.join(metals)} {metal_basis_print}

Charge: {charge}
-------------
{fspe_lines}

TOTAL RUN TIME: {duration_format}

{lowest_label} {lowest_str}

(Selection: {method}, APmethod {ap_str}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})
(Energy-bias: energy_bias_window_h={E_BIAS_H}, mismatch_bias_window={MIS_BIAS} -> prefer lower E when close)
(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})
(Electron number: {parity})
""")
