from delfin.common.banners import build_occupier_banner

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
    """
    # ----------------------- imports & tiny helpers ----------------------------
    from decimal import Decimal, ROUND_DOWN
    from pathlib import Path
    from typing import Optional
    import re

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
            # deviation-prior then deviation then index
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
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm⁻¹"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(⟨S²⟩HS-⟨S²⟩BS) : {j3_str}\n"
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
    banner = build_occupier_banner(header_indent=6, info_indent=6)

    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(
            f"{banner}\n\n"
            f"{method_str}\n"
            f"        {', '.join(metals)} {metal_basis_print}\n\n"
            f"Charge: {charge}     \n"
            f"-------------  \n"
            f"{fspe_lines}\n\n"
            f"TOTAL RUN TIME: {duration_format}\n\n"
            f"{lowest_label} {lowest_str}\n\n"
            f"(Selection: {method}, APmethod {ap_str}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})\n"
            f"(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})\n"
            f"(Electron number: {parity})\n"
        )

def generate_summary_report_OCCUPIER_safe3(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence):
    """
    Decide the preferred occupation (entry) among energetically equivalent candidates
    using only spin contamination. J(3) is parsed and printed but NOT used for selection.

    Selection rules (summary):
      1) Energy decides first (with user-chosen rounding/truncation/tolerance).
      2) If energies tie:
         a) Prefer clearly smaller spin contamination.
         b) Use a deviation "prior" toward the expected BS pair count:
            - Deviation ~ 0  -> prefer non-BS (restricted-like).
            - Deviation ~ 1  -> prefer BS(*,1).
            - Deviation ~ 2  -> prefer BS(*,2).
            - etc.
            Concretely: minimize |Deviation - n_BS| for BS entries, with a window to treat "matches" equally.
         c) BS singlet with tiny deviation (<= DEV_TINY) is treated as pseudo closed-shell (non-BS).
         d) Final tie-breakers: smaller deviation, then lower index.
      3) AF-override (no J): if the only candidate in the energy band is highly contaminated (>= DEV_HIGH),
         consider BS candidates within ΔE <= bs_override_window_h (default 0.002 Ha) whose deviations are also
         high and similar; among them prefer smallest |Deviation - n_BS|, then smaller deviation, then index.

    Note: J(3) (Yamaguchi) is printed for BS entries but not used for selection.
    """
    # ----------------------- numeric formatting helpers ------------------------
    from decimal import Decimal, ROUND_DOWN
    from typing import Optional
    import re

    def truncate(x: float, d: int) -> float:
        q = Decimal(10) ** -d
        return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))

    def fmt_truncate(x: float, d: int) -> str:
        return f"{truncate(x, d):.{d}f}"

    # ----------------------- extraction helpers (outputs) ----------------------
    def file_for_index(idx: int) -> str:
        """
        Map sequence index -> out file. Supports dict/list in config['out_files'],
        else falls back to 'output.out' / f'output{idx}.out'.
        """
        out_files = config.get('out_files')
        if isinstance(out_files, dict):
            return out_files.get(idx) or ("output.out" if idx == 1 else f"output{idx}.out")
        if isinstance(out_files, (list, tuple)) and 1 <= idx <= len(out_files):
            return out_files[idx-1] or ("output.out" if idx == 1 else f"output{idx}.out")
        return "output.out" if idx == 1 else f"output{idx}.out"

    # ----------------------- extraction helpers (inputs) -----------------------
    def _inp_candidates_for_index(idx: int) -> list[str]:
        """
        Build a list of plausible .inp names for a given index.
        1) Respect config['inp_files'] if present (dict/list like out_files).
        2) Derive from out filename by swapping extension / 'output'->'input'.
        3) Common fallbacks: input.inp, input{idx}.inp, input_{idx}.inp
        """
        cands: list[str] = []

        inp_files = config.get('inp_files')
        if isinstance(inp_files, dict):
            v = inp_files.get(idx)
            if v:
                cands.append(str(v))
        elif isinstance(inp_files, (list, tuple)) and 1 <= idx <= len(inp_files):
            v = inp_files[idx-1]
            if v:
                cands.append(str(v))

        ofile = file_for_index(idx)
        base = Path(ofile).stem
        cands.append(base + ".inp")
        if base.lower().startswith("output"):
            cands.append("input" + base[6:] + ".inp")
        cands.append(ofile.replace("output", "input").rsplit(".", 1)[0] + ".inp")

        if idx == 1:
            cands.append("input.inp")
        cands.append(f"input{idx}.inp")
        cands.append(f"input_{idx}.inp")

        seen = set()
        unique = []
        for p in cands:
            p2 = p.strip()
            if p2 and p2 not in seen:
                unique.append(p2); seen.add(p2)
        return unique

    def _inp_for_index(idx: int) -> Optional[str]:
        for p in _inp_candidates_for_index(idx):
            if Path(p).exists():
                return p
        return None

    def _parse_method_from_inp(inp_path: str, pal: str, include_freq: bool) -> Optional[str]:
        """
        Parse the first '!' line from the .inp file.
        Drop tokens like PMODEL/MOREAD; ensure FREQ if include_freq=True; append PAL.
        """
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
                        pal_str = str(pal).strip()
                        if pal_str:
                            method = (method + f" PAL{pal_str}").replace("  ", " ").strip()
                        return method
        except Exception:
            pass
        return None

    def _parse_metal_basis_from_inp(inp_path: str) -> Optional[str]:
        """
        Return the first NewGTO "BASISNAME" found inside the first * xyz ... * block.
        """
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

    # ----------------------- CONTROL fallbacks ---------------------------------
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
        if not model:
            return ""
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return f"{model}({s})" if s else model

    def method_line_from_control(include_freq: bool) -> str:
        tokens = [
            "Method:",
            str(config.get('functional','')).strip(),
            rel_token_ctl,
            str(main_basisset).strip(),
            str(config.get('disp_corr','')).strip(),
            str(config.get('ri_jkx','')).strip(),
            aux_jk_token_ctl,
            implicit_token(),
            str(config.get('geom_opt_OCCUPIER','')).strip()
        ]
        if include_freq:
            tokens.append("FREQ")
        tokens.append(f"PAL{config.get('PAL','')}")
        return " ".join(t for t in tokens if t).replace("  ", " ").strip()

    # ----------------------- time formatting -----------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    # ----------------------- energy / selection setup --------------------------
    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""
    parity = "is_even" if is_even else "is_odd"
    use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    energy_label = "Final Gibbs free energy" if use_gibbs else "FINAL SINGLE POINT ENERGY"
    lowest_label = f"LOWEST {energy_label}:"

    raw_sel = str(config.get('occupier_selection', 'tolerance')).lower()
    sel = raw_sel.split('|')[0].strip()
    ap = config.get('approximate_spin_projection_APMethod')
    ap_str = str(ap) if ap not in (None, 0, '0') else 'none'

    if sel in {'rounding', 'round', 'gerundet', 'runden'}:
        method = 'rounding'
    elif sel in {'tolerance', 'toleranz', 'toleranzband', 'epsilon'}:
        method = 'tolerance'
    else:
        method = 'truncation'

    try:
        prec = int(config.get('occupier_precision', config.get('occupier_rounded_value', 6)))
    except (TypeError, ValueError):
        prec = 6
    prec = max(0, min(10, prec))

    try:
        epsilon = float(config.get('occupier_epsilon', 10.0**(-prec)))
        if not (epsilon > 0):
            raise ValueError
    except (TypeError, ValueError):
        epsilon = 10.0**(-prec)

    dev_max_raw = config.get('dev_max', None)
    try:
        dev_max = float(dev_max_raw) if dev_max_raw is not None else None
    except (TypeError, ValueError):
        dev_max = None

    # ----------------------- decision thresholds -------------------------------
    DEV_TINY         = 1e-3   # treat BS singlet with ~zero deviation as closed-shell
    DEV_SIMILARITY   = 0.15   # "similar deviation" band
    DEV_GOOD_MARGIN  = 0.30   # "clearly smaller" margin
    DEV_HIGH         = 0.50   # high contamination regime
    DEV_MATCH_WINDOW = 0.30   # how close deviation should be to n_BS for the prior

    # AF-override extra energy window in Hartree (default 0.002 Ha)
    EPS_AF = float(config.get('bs_override_window_h', 0.002))

    # ----------------------- caches for deviation and J(3) ---------------------
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
        """Parse J(3) for printing only; NEVER used in selection."""
        if not is_bs_flag:
            return None
        if idx not in j3_cache:
            try:
                j3_cache[idx] = extract_last_J3(file_for_index(idx))
            except Exception:
                j3_cache[idx] = None
        return j3_cache[idx]

    # ----------------------- sequence accessors --------------------------------
    def entry_record(idx: int) -> dict:
        return next((e for e in sequence if e["index"] == idx), {})

    def entry_is_bs(idx: int) -> bool:
        return bool(entry_record(idx).get("BS"))

    def entry_mult(idx: int) -> int:
        try:
            return int(entry_record(idx).get("m", 0))
        except Exception:
            return 0

    def bs_pair_count(idx: int) -> Optional[int]:
        """Return n from 'a,n' in BS label, e.g. '2,1' -> 1; None if not BS or unparsable."""
        lab = str(entry_record(idx).get("BS") or "").strip()
        if not lab or "," not in lab:
            return None
        try:
            _, n = lab.split(",", 1)
            return int(n.strip())
        except Exception:
            return None

    # ----------------------- deviation utilities --------------------------------
    def effective_dev(idx: int) -> float:
        """
        Effective deviation used for ranking:
          - For non-BS singlets (m == 1) with no UHF block → 0.0 (closed-shell RKS).
          - Otherwise use parsed deviation; if missing → +inf.
        """
        dev = get_dev(idx)
        if dev is None:
            if (not entry_is_bs(idx)) and (entry_mult(idx) == 1):
                return 0.0
            return float("inf")
        return dev

    def bs_mismatch(idx: int) -> float:
        """|Deviation - n_BS|; inf if n_BS unknown."""
        n = bs_pair_count(idx)
        if n is None:
            return float("inf")
        return abs(effective_dev(idx) - n)

    def is_pseudo_closed_shell(idx: int) -> bool:
        """BS singlet with tiny deviation behaves like closed-shell → treat like non-BS."""
        return entry_is_bs(idx) and (entry_mult(idx) == 1) and (effective_dev(idx) <= DEV_TINY)

    def within_dev_limit(idx: int) -> bool:
        if dev_max is None:
            return True
        d = get_dev(idx)
        return (d is None) or (d <= dev_max)

    # ----------------------- candidate set by energy policy --------------------
    valid_all = [(e["index"], f) for e, f in zip(sequence, fspe_values) if f is not None]
    valid = [pair for pair in valid_all if within_dev_limit(pair[0])]
    if not valid and dev_max is not None:
        valid = valid_all

    def energy_key(val: float) -> float:
        if method == "rounding":
            return round(val, prec)
        if method == "truncation":
            return truncate(val, prec)
        return val  # tolerance => raw energies

    min_fspe_index: Optional[int] = None
    min_fspe_value: Optional[float] = None

    if valid:
        # Build the energy band according to policy
        if method in ("rounding", "truncation"):
            processed = [(i, f, energy_key(f)) for i, f in valid]
            best_key = min(v for _, _, v in processed)
            cands = [(i, f) for i, f, v in processed if v == best_key]
            min_raw = min(f for _, f in valid)  # for AF-override window
        else:
            min_raw = min(f for _, f in valid)
            eps_eff = epsilon + 1e-12
            cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

        # --- AF-override: single best but highly contaminated -> look for BS nearby (no J)
        if len(cands) == 1:
            best_idx, best_val = cands[0]
            best_dev = effective_dev(best_idx)
            if best_dev >= DEV_HIGH:
                bs_pool = []
                for i, f in valid:
                    if i == best_idx or not entry_is_bs(i):
                        continue
                    d = effective_dev(i)
                    if d >= DEV_HIGH and abs(d - best_dev) <= DEV_SIMILARITY and (f - min_raw) <= EPS_AF + 1e-12:
                        bs_pool.append((i, f))
                if bs_pool:
                    # Prefer smallest mismatch to n_BS, then smaller deviation, then index
                    def bs_pick_score(i):
                        return (bs_mismatch(i) > DEV_MATCH_WINDOW, bs_mismatch(i), effective_dev(i), i)
                    picked = min(bs_pool, key=lambda p: bs_pick_score(p[0]))
                    cands = [picked]

        if len(cands) == 1:
            min_fspe_index, min_fspe_value = cands[0]
        else:
            # Split candidates; reclassify pseudo-closed-shell BS to non-BS
            raw_bs = [(i, f) for i, f in cands if entry_is_bs(i)]
            raw_nb = [(i, f) for i, f in cands if not entry_is_bs(i)]
            nb_cands = raw_nb + [(i, f) for i, f in raw_bs if is_pseudo_closed_shell(i)]
            bs_cands = [(i, f) for i, f in raw_bs if not is_pseudo_closed_shell(i)]

            def _min_eff_dev(pairs):
                vals = [effective_dev(i) for i, _ in pairs]
                return min(vals) if vals else float("inf")

            min_dev_bs = _min_eff_dev(bs_cands)
            min_dev_nb = _min_eff_dev(nb_cands)

            chosen: Optional[tuple[int, float]] = None

            if bs_cands and nb_cands:
                # A) Non-BS clearly cleaner and not in high regime → prefer non-BS
                if (min_dev_nb + DEV_GOOD_MARGIN < min_dev_bs) and (min_dev_nb < DEV_HIGH):
                    def nb_score(i):
                        return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                    chosen = min(nb_cands, key=lambda p: nb_score(p[0]))
                # B) Similar and both high → prefer BS; among BS use mismatch prior (no J)
                elif (abs(min_dev_nb - min_dev_bs) <= DEV_SIMILARITY) and (min_dev_nb > DEV_HIGH) and (min_dev_bs > DEV_HIGH):
                    def bs_score(i):
                        return (bs_mismatch(i) > DEV_MATCH_WINDOW, bs_mismatch(i), effective_dev(i), i)
                    chosen = min(bs_cands, key=lambda p: bs_score(p[0]))
                else:
                    # C) Fallback: smallest deviation; tie → prefer non-BS; then lowest index
                    def fb_score(i):
                        return (effective_dev(i), 1 if entry_is_bs(i) else 0, i)
                    chosen = min(cands, key=lambda p: fb_score(p[0]))
            elif bs_cands:
                # Only BS → minimize mismatch, then deviation, then index
                def bs_only_score(i):
                    return (bs_mismatch(i) > DEV_MATCH_WINDOW, bs_mismatch(i), effective_dev(i), i)
                chosen = min(bs_cands, key=lambda p: bs_only_score(p[0]))
            else:
                # Only non-BS → smallest deviation; prefer genuine non-BS over pseudo
                def nb_only_score(i):
                    return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                chosen = min(nb_cands, key=lambda p: nb_only_score(p[0]))

            if chosen is not None:
                min_fspe_index, min_fspe_value = chosen

    # ----------------------- build FSPE table text -----------------------------
    def fmt_display(x: Optional[float]) -> str:
        if x is None:
            return "N/A"
        return f"{x:.{prec}f}" if method == 'rounding' else fmt_truncate(x, prec)

    fspe_lines = ""
    for entry, fspe in zip(sequence, fspe_values):
        idx = entry["index"]
        multiplicity = entry["m"]
        bs = entry.get("BS", "")
        is_bs_flag = bool(bs)

        fspe_value = (fmt_display(fspe) + " (H)") if fspe is not None else "Not a valid value"
        dev_val = get_dev(idx)
        dev_txt = fmt_display(dev_val)
        mark = " <-- PREFERRED VALUE" if (min_fspe_index is not None and idx == min_fspe_index) else ""
        bs_line = f", BrokenSym {bs}" if bs else ""

        fspe_lines += f"{energy_label} ({idx})   = {fspe_value}{mark}\n"
        fspe_lines += f"multiplicity {multiplicity}{bs_line}\n"
        fspe_lines += f"Spin Contamination (⟨S²⟩ - S(S+1))   : {dev_txt}\n"

        # Print J(3) only for BS entries (for information only)
        if is_bs_flag:
            j3_val = get_j3(idx, True)
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm⁻¹"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(⟨S²⟩HS-⟨S²⟩BS) : {j3_str}\n"

        fspe_lines += "----------------------------------------------------------------\n"

    print("\n" + fspe_lines)

    # ----------------------- compose method & metal-basis lines ----------------
    include_freq = (freq_flag == 'FREQ')
    pal = config.get('PAL', '')

    method_str = None
    metal_basis_from_inp = None

    prefer_idx = (min_fspe_index if min_fspe_index is not None else 1)
    inp_pref = _inp_for_index(prefer_idx) or _inp_for_index(1)

    if inp_pref:
        parsed = _parse_method_from_inp(inp_pref, pal, include_freq)
        if parsed:
            method_str = parsed
        metal_basis_from_inp = _parse_metal_basis_from_inp(inp_pref)

    if not method_str:
        method_str = method_line_from_control(include_freq)

    metal_basis_print = metal_basis_from_inp or metal_basisset_rep_ctl or ""

    lowest_str = (
        f"{fmt_truncate(min_fspe_value, prec)} (H)" if min_fspe_value is not None
        else "No valid FSPE values found"
    )

    # ---- write report ----------------------------------------------------------
    banner = build_occupier_banner(header_indent=6, info_indent=6)

    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(
            f"{banner}\n\n"
            f"{method_str}\n"
            f"        {', '.join(metals)} {metal_basis_print}\n\n"
            f"Charge: {charge}     \n"
            f"-------------  \n"
            f"{fspe_lines}\n\n"
            f"TOTAL RUN TIME: {duration_format}\n\n"
            f"{lowest_label} {lowest_str}\n\n"
            f"(Selection: {method}, APmethod {ap_str}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})\n"
            f"(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})\n"
            f"(Electron number: {parity})\n"
        )


#-------------------------------------------------------OCCUPIER Report---------------------------------------------------------------------------
def generate_summary_report_OCCUPIER_safe2(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence):
    """
    Decide the preferred occupation (entry) among energetically equivalent candidates,
    using spin contamination and, for BS jobs, the Yamaguchi J(3) coupling as tie-breakers.

    Selection rules (summary):
      1) Energy decides first (with user-chosen rounding/truncation/tolerance).
      2) If energies tie:
         a) Prefer the one with clearly smaller spin contamination.
         b) If both deviations are high (> DEV_HIGH) and similar, prefer BS (AF-like open-shell).
         c) Otherwise prefer non-BS; among equals, prefer lower deviation.
      3) BS with m == 1 and tiny deviation (<= DEV_TINY) is treated like non-BS (pseudo closed-shell).
      4) Among BS ties, prefer (i) smaller |Deviation - n_BS|, (ii) larger |J(3)|, (iii) smaller deviation.
      5) AF-override: if the sole best-energy candidate is highly contaminated (>= DEV_HIGH),
         consider BS candidates within ΔE <= bs_override_window_h (default 0.001 Ha)
         whose deviations are also high and similar to the best candidate.
    """
    # ----------------------- numeric formatting helpers ------------------------
    from decimal import Decimal, ROUND_DOWN
    from typing import Optional
    import re

    def truncate(x: float, d: int) -> float:
        q = Decimal(10) ** -d
        return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))

    def fmt_truncate(x: float, d: int) -> str:
        return f"{truncate(x, d):.{d}f}"

    # ----------------------- extraction helpers (outputs) ----------------------
    def file_for_index(idx: int) -> str:
        """
        Map sequence index -> out file. Supports dict/list in config['out_files'],
        else falls back to 'output.out' / f'output{idx}.out'.
        """
        out_files = config.get('out_files')
        if isinstance(out_files, dict):
            return out_files.get(idx) or ("output.out" if idx == 1 else f"output{idx}.out")
        if isinstance(out_files, (list, tuple)) and 1 <= idx <= len(out_files):
            return out_files[idx-1] or ("output.out" if idx == 1 else f"output{idx}.out")
        return "output.out" if idx == 1 else f"output{idx}.out"

    # ----------------------- extraction helpers (inputs) -----------------------
    def _inp_candidates_for_index(idx: int) -> list[str]:
        """
        Build a list of plausible .inp names for a given index.
        1) Respect config['inp_files'] if present (dict/list like out_files).
        2) Derive from out filename by swapping extension / 'output'->'input'.
        3) Common fallbacks: input.inp, input{idx}.inp, input_{idx}.inp
        """
        cands: list[str] = []

        inp_files = config.get('inp_files')
        if isinstance(inp_files, dict):
            v = inp_files.get(idx)
            if v:
                cands.append(str(v))
        elif isinstance(inp_files, (list, tuple)) and 1 <= idx <= len(inp_files):
            v = inp_files[idx-1]
            if v:
                cands.append(str(v))

        ofile = file_for_index(idx)
        base = Path(ofile).stem
        cands.append(base + ".inp")
        if base.lower().startswith("output"):
            cands.append("input" + base[6:] + ".inp")
        cands.append(ofile.replace("output", "input").rsplit(".", 1)[0] + ".inp")

        if idx == 1:
            cands.append("input.inp")
        cands.append(f"input{idx}.inp")
        cands.append(f"input_{idx}.inp")

        seen = set()
        unique = []
        for p in cands:
            p2 = p.strip()
            if p2 and p2 not in seen:
                unique.append(p2); seen.add(p2)
        return unique

    def _inp_for_index(idx: int) -> Optional[str]:
        for p in _inp_candidates_for_index(idx):
            if Path(p).exists():
                return p
        return None

    def _parse_method_from_inp(inp_path: str, pal: str, include_freq: bool) -> Optional[str]:
        """
        Parse the first '!' line from the .inp file.
        Drop tokens like PMODEL/MOREAD; ensure FREQ if include_freq=True; append PAL.
        """
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
                        pal_str = str(pal).strip()
                        if pal_str:
                            method = (method + f" PAL{pal_str}").replace("  ", " ").strip()
                        return method
        except Exception:
            pass
        return None

    def _parse_metal_basis_from_inp(inp_path: str) -> Optional[str]:
        """
        Return the first NewGTO "BASISNAME" found inside the first * xyz ... * block.
        """
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

    # ----------------------- CONTROL fallbacks ---------------------------------
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
        if not model:
            return ""
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return f"{model}({s})" if s else model

    def method_line_from_control(include_freq: bool) -> str:
        tokens = [
            "Method:",
            str(config.get('functional','')).strip(),
            rel_token_ctl,
            str(main_basisset).strip(),
            str(config.get('disp_corr','')).strip(),
            str(config.get('ri_jkx','')).strip(),
            aux_jk_token_ctl,
            implicit_token(),
            str(config.get('geom_opt_OCCUPIER','')).strip()
        ]
        if include_freq:
            tokens.append("FREQ")
        tokens.append(f"PAL{config.get('PAL','')}")
        return " ".join(t for t in tokens if t).replace("  ", " ").strip()

    # ----------------------- time formatting -----------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    # ----------------------- energy / selection setup --------------------------
    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""
    parity = "is_even" if is_even else "is_odd"
    use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    energy_label = "Final Gibbs free energy" if use_gibbs else "FINAL SINGLE POINT ENERGY"
    lowest_label = f"LOWEST {energy_label}:"

    raw_sel = str(config.get('occupier_selection', 'tolerance')).lower()
    sel = raw_sel.split('|')[0].strip()
    ap = config.get('approximate_spin_projection_APMethod')
    ap_str = str(ap) if ap not in (None, 0, '0') else 'none'

    if sel in {'rounding', 'round', 'gerundet', 'runden'}:
        method = 'rounding'
    elif sel in {'tolerance', 'toleranz', 'toleranzband', 'epsilon'}:
        method = 'tolerance'
    else:
        method = 'truncation'

    try:
        prec = int(config.get('occupier_precision', config.get('occupier_rounded_value', 6)))
    except (TypeError, ValueError):
        prec = 6
    prec = max(0, min(10, prec))

    try:
        epsilon = float(config.get('occupier_epsilon', 10.0**(-prec)))
        if not (epsilon > 0):
            raise ValueError
    except (TypeError, ValueError):
        epsilon = 10.0**(-prec)

    dev_max_raw = config.get('dev_max', None)
    try:
        dev_max = float(dev_max_raw) if dev_max_raw is not None else None
    except (TypeError, ValueError):
        dev_max = None

    # ----------------------- decision thresholds -------------------------------
    DEV_TINY         = 1e-3   # treat BS singlet with ~zero deviation as closed-shell
    DEV_SIMILARITY   = 0.15   # "similar deviation" band
    DEV_GOOD_MARGIN  = 0.30   # "clearly smaller" margin
    DEV_HIGH         = 0.50   # high contamination regime
    DEV_MATCH_WINDOW = 0.30   # how close deviation should be to n_BS for the prior

    # AF-override extra energy window in Hartree (default 0.002 Ha)
    EPS_AF = float(config.get('bs_override_window_h', 0.002))

    # ----------------------- caches for deviation and J(3) ---------------------
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
        if not is_bs_flag:
            return None
        if idx not in j3_cache:
            try:
                j3_cache[idx] = extract_last_J3(file_for_index(idx))
            except Exception:
                j3_cache[idx] = None
        return j3_cache[idx]

    # ----------------------- sequence accessors --------------------------------
    def entry_record(idx: int) -> dict:
        return next((e for e in sequence if e["index"] == idx), {})

    def entry_is_bs(idx: int) -> bool:
        return bool(entry_record(idx).get("BS"))

    def entry_mult(idx: int) -> int:
        try:
            return int(entry_record(idx).get("m", 0))
        except Exception:
            return 0

    def bs_pair_count(idx: int) -> Optional[int]:
        """Return n from 'a,n' in BS label, e.g. '2,1' -> 1; None if not BS or unparsable."""
        lab = str(entry_record(idx).get("BS") or "").strip()
        if not lab or "," not in lab:
            return None
        try:
            _, n = lab.split(",", 1)
            return int(n.strip())
        except Exception:
            return None

    # ----------------------- deviation utilities --------------------------------
    def effective_dev(idx: int) -> float:
        """
        Effective deviation used for ranking:
          - For non-BS singlets (m == 1) with no UHF block → 0.0 (closed-shell RKS).
          - Otherwise use parsed deviation; if missing → +inf.
        """
        dev = get_dev(idx)
        if dev is None:
            if (not entry_is_bs(idx)) and (entry_mult(idx) == 1):
                return 0.0
            return float("inf")
        return dev

    def bs_mismatch(idx: int) -> float:
        """|Deviation - n_BS|; inf if n_BS unknown."""
        n = bs_pair_count(idx)
        if n is None:
            return float("inf")
        return abs(effective_dev(idx) - n)

    def is_pseudo_closed_shell(idx: int) -> bool:
        """BS singlet with tiny deviation behaves like closed-shell → treat like non-BS."""
        return entry_is_bs(idx) and (entry_mult(idx) == 1) and (effective_dev(idx) <= DEV_TINY)

    def within_dev_limit(idx: int) -> bool:
        if dev_max is None:
            return True
        d = get_dev(idx)
        return (d is None) or (d <= dev_max)

    # ----------------------- candidate set by energy policy --------------------
    valid_all = [(e["index"], f) for e, f in zip(sequence, fspe_values) if f is not None]
    valid = [pair for pair in valid_all if within_dev_limit(pair[0])]
    if not valid and dev_max is not None:
        valid = valid_all

    def energy_key(val: float) -> float:
        if method == "rounding":
            return round(val, prec)
        if method == "truncation":
            return truncate(val, prec)
        return val  # tolerance => raw energies

    min_fspe_index: Optional[int] = None
    min_fspe_value: Optional[float] = None

    if valid:
        # Build the energy band according to policy
        if method in ("rounding", "truncation"):
            processed = [(i, f, energy_key(f)) for i, f in valid]
            best_key = min(v for _, _, v in processed)
            cands = [(i, f) for i, f, v in processed if v == best_key]
            min_raw = min(f for _, f in valid)  # for AF-override window
        else:
            min_raw = min(f for _, f in valid)
            eps_eff = epsilon + 1e-12
            cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

        # --- AF-override: single best but highly contaminated -> look for BS near by
        if len(cands) == 1:
            best_idx, best_val = cands[0]
            best_dev = effective_dev(best_idx)
            if best_dev >= DEV_HIGH:
                bs_pool = []
                for i, f in valid:
                    if i == best_idx or not entry_is_bs(i):
                        continue
                    d = effective_dev(i)
                    if d >= DEV_HIGH and abs(d - best_dev) <= DEV_SIMILARITY and (f - min_raw) <= EPS_AF + 1e-12:
                        bs_pool.append((i, f))
                if bs_pool:
                    # prioritize BS with small |dev - n_BS| first, then larger |J|, then smaller dev
                    def bs_pick_score(i):
                        mis = bs_mismatch(i)
                        mis_in = 0 if mis <= DEV_MATCH_WINDOW else 1
                        j = get_j3(i, True)
                        mag = abs(j) if (j is not None) else -float("inf")
                        return (mis_in, mis, -mag, effective_dev(i), i)
                    picked = min(bs_pool, key=lambda p: bs_pick_score(p[0]))
                    cands = [picked]

        if len(cands) == 1:
            min_fspe_index, min_fspe_value = cands[0]
        else:
            # Split candidates; reclassify pseudo-closed-shell BS to non-BS
            raw_bs = [(i, f) for i, f in cands if entry_is_bs(i)]
            raw_nb = [(i, f) for i, f in cands if not entry_is_bs(i)]
            nb_cands = raw_nb + [(i, f) for i, f in raw_bs if is_pseudo_closed_shell(i)]
            bs_cands = [(i, f) for i, f in raw_bs if not is_pseudo_closed_shell(i)]

            def _min_eff_dev(pairs):
                vals = [effective_dev(i) for i, _ in pairs]
                return min(vals) if vals else float("inf")

            min_dev_bs = _min_eff_dev(bs_cands)
            min_dev_nb = _min_eff_dev(nb_cands)

            chosen: Optional[tuple[int, float]] = None

            if bs_cands and nb_cands:
                # A) Non-BS clearly cleaner and not in high regime → prefer non-BS
                if (min_dev_nb + DEV_GOOD_MARGIN < min_dev_bs) and (min_dev_nb < DEV_HIGH):
                    def nb_score(i):
                        return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                    chosen = min(nb_cands, key=lambda p: nb_score(p[0]))
                # B) Similar and both high → prefer BS; among BS add Deviation-Prior and |J(3)|
                elif (abs(min_dev_nb - min_dev_bs) <= DEV_SIMILARITY) and (min_dev_nb > DEV_HIGH) and (min_dev_bs > DEV_HIGH):
                    def bs_score(i):
                        mis = bs_mismatch(i)
                        mis_in = 0 if mis <= DEV_MATCH_WINDOW else 1
                        j = get_j3(i, True)
                        mag = abs(j) if (j is not None) else -float("inf")
                        return (mis_in, mis, -mag, effective_dev(i), i)
                    chosen = min(bs_cands, key=lambda p: bs_score(p[0]))
                else:
                    # C) Fallback: smallest deviation; tie → prefer non-BS; then lowest index
                    def fb_score(i):
                        return (effective_dev(i), 1 if entry_is_bs(i) else 0, i)
                    chosen = min(cands, key=lambda p: fb_score(p[0]))
            elif bs_cands:
                # Only BS → Deviation-Prior then |J(3)| then deviation
                def bs_only_score(i):
                    mis = bs_mismatch(i)
                    mis_in = 0 if mis <= DEV_MATCH_WINDOW else 1
                    j = get_j3(i, True)
                    mag = abs(j) if (j is not None) else -float("inf")
                    return (mis_in, mis, -mag, effective_dev(i), i)
                chosen = min(bs_cands, key=lambda p: bs_only_score(p[0]))
            else:
                # Only non-BS → smallest deviation; prefer genuine non-BS over pseudo
                def nb_only_score(i):
                    return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                chosen = min(nb_cands, key=lambda p: nb_only_score(p[0]))

            if chosen is not None:
                min_fspe_index, min_fspe_value = chosen

    # ----------------------- build FSPE table text -----------------------------
    def fmt_display(x: Optional[float]) -> str:
        if x is None:
            return "N/A"
        return f"{x:.{prec}f}" if method == 'rounding' else fmt_truncate(x, prec)

    fspe_lines = ""
    for entry, fspe in zip(sequence, fspe_values):
        idx = entry["index"]
        multiplicity = entry["m"]
        bs = entry.get("BS", "")
        is_bs_flag = bool(bs)

        fspe_value = (fmt_display(fspe) + " (H)") if fspe is not None else "Not a valid value"
        dev_val = get_dev(idx)
        dev_txt = fmt_display(dev_val)
        mark = " <-- PREFERRED VALUE" if (min_fspe_index is not None and idx == min_fspe_index) else ""
        bs_line = f", BrokenSym {bs}" if bs else ""

        fspe_lines += f"{energy_label} ({idx})   = {fspe_value}{mark}\n"
        fspe_lines += f"multiplicity {multiplicity}{bs_line}\n"
        fspe_lines += f"Spin Contamination (⟨S²⟩ - S(S+1))   : {dev_txt}\n"

        if is_bs_flag:
            j3_val = get_j3(idx, True)
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm⁻¹"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(⟨S²⟩HS-⟨S²⟩BS) : {j3_str}\n"

        fspe_lines += "----------------------------------------------------------------\n"

    print("\n" + fspe_lines)

    # ----------------------- compose method & metal-basis lines ----------------
    include_freq = (freq_flag == 'FREQ')
    pal = config.get('PAL', '')

    method_str = None
    metal_basis_from_inp = None

    prefer_idx = (min_fspe_index if min_fspe_index is not None else 1)
    inp_pref = _inp_for_index(prefer_idx) or _inp_for_index(1)

    if inp_pref:
        parsed = _parse_method_from_inp(inp_pref, pal, include_freq)
        if parsed:
            method_str = parsed
        metal_basis_from_inp = _parse_metal_basis_from_inp(inp_pref)

    if not method_str:
        method_str = method_line_from_control(include_freq)

    metal_basis_print = metal_basis_from_inp or metal_basisset_rep_ctl or ""

    lowest_str = (
        f"{fmt_truncate(min_fspe_value, prec)} (H)" if min_fspe_value is not None
        else "No valid FSPE values found"
    )

    # ---- write report ----------------------------------------------------------
    banner = build_occupier_banner(header_indent=6, info_indent=6)

    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(
            f"{banner}\n\n"
            f"{method_str}\n"
            f"        {', '.join(metals)} {metal_basis_print}\n\n"
            f"Charge: {charge}     \n"
            f"-------------  \n"
            f"{fspe_lines}\n\n"
            f"TOTAL RUN TIME: {duration_format}\n\n"
            f"{lowest_label} {lowest_str}\n\n"
            f"(Selection: {method}, APmethod {ap_str}, dev_max {dev_max if dev_max is not None else 'none'}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})\n"
            f"(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})\n"
            f"(Electron number: {parity})\n"
        )


#-------------------------------------------------------OCCUPIER Report---------------------------------------------------------------------------
def generate_summary_report_OCCUPIER_safe(duration, fspe_values, is_even, charge, solvent, config, main_basisset, sequence):
    """
    Decide the preferred occupation (entry) among energetically equivalent candidates,
    using spin contamination and, for BS jobs, the Yamaguchi J(3) coupling as tie-breakers.

    Selection rules (summary):
      1) Energy decides first (with user-chosen rounding/truncation/tolerance).
      2) If energies tie:
         a) Prefer the one with clearly smaller spin contamination.
         b) If both deviations are high (> DEV_HIGH) and similar, prefer BS (AF-like open-shell).
         c) Otherwise prefer non-BS; among equals, prefer lower deviation.
      3) BS with m == 1 and tiny deviation (<= DEV_TINY) is treated like non-BS (pseudo closed-shell).
      4) Among BS ties, prefer larger |J(3)|; missing J(3) is considered worst.
      5) AF-override: if the single best-energy candidate shows high contamination (>= DEV_HIGH),
         consider BS candidates within an extra window ΔE <= bs_override_window_h (default 0.001 Ha)
         provided their deviations are similarly high and similar to the best candidate.
    """
    # ----------------------- numeric formatting helpers ------------------------
    from decimal import Decimal, ROUND_DOWN
    from typing import Optional
    import re

    def truncate(x: float, d: int) -> float:
        q = Decimal(10) ** -d
        return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))

    def fmt_truncate(x: float, d: int) -> str:
        return f"{truncate(x, d):.{d}f}"

    # ----------------------- extraction helpers (outputs) ----------------------
    def file_for_index(idx: int) -> str:
        """
        Map sequence index -> out file. Supports dict/list in config['out_files'],
        else falls back to 'output.out' / f'output{idx}.out'.
        """
        out_files = config.get('out_files')
        if isinstance(out_files, dict):
            return out_files.get(idx) or ("output.out" if idx == 1 else f"output{idx}.out")
        if isinstance(out_files, (list, tuple)) and 1 <= idx <= len(out_files):
            return out_files[idx-1] or ("output.out" if idx == 1 else f"output{idx}.out")
        return "output.out" if idx == 1 else f"output{idx}.out"

    # ----------------------- extraction helpers (inputs) -----------------------
    def _inp_candidates_for_index(idx: int) -> list[str]:
        """
        Build a list of plausible .inp names for a given index.
        1) Respect config['inp_files'] if present (dict/list like out_files).
        2) Derive from out filename by swapping extension / 'output'->'input'.
        3) Common fallbacks: input.inp, input{idx}.inp, input_{idx}.inp
        """
        cands: list[str] = []

        inp_files = config.get('inp_files')
        if isinstance(inp_files, dict):
            v = inp_files.get(idx)
            if v:
                cands.append(str(v))
        elif isinstance(inp_files, (list, tuple)) and 1 <= idx <= len(inp_files):
            v = inp_files[idx-1]
            if v:
                cands.append(str(v))

        ofile = file_for_index(idx)
        base = Path(ofile).stem
        cands.append(base + ".inp")
        if base.lower().startswith("output"):
            cands.append("input" + base[6:] + ".inp")
        cands.append(ofile.replace("output", "input").rsplit(".", 1)[0] + ".inp")

        if idx == 1:
            cands.append("input.inp")
        cands.append(f"input{idx}.inp")
        cands.append(f"input_{idx}.inp")

        seen = set()
        unique = []
        for p in cands:
            p2 = p.strip()
            if p2 and p2 not in seen:
                unique.append(p2); seen.add(p2)
        return unique

    def _inp_for_index(idx: int) -> Optional[str]:
        for p in _inp_candidates_for_index(idx):
            if Path(p).exists():
                return p
        return None

    def _parse_method_from_inp(inp_path: str, pal: str, include_freq: bool) -> Optional[str]:
        """
        Parse the first '!' line from the .inp file.
        Drop tokens like PMODEL/MOREAD; ensure FREQ if include_freq=True; append PAL.
        """
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
                        pal_str = str(pal).strip()
                        if pal_str:
                            method = (method + f" PAL{pal_str}").replace("  ", " ").strip()
                        return method
        except Exception:
            pass
        return None

    def _parse_metal_basis_from_inp(inp_path: str) -> Optional[str]:
        """
        Return the first NewGTO "BASISNAME" found inside the first * xyz ... * block.
        """
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

    # ----------------------- CONTROL fallbacks ---------------------------------
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
        if not model:
            return ""
        s = (solvent or str(config.get('solvent', '')).strip() or '').strip()
        return f"{model}({s})" if s else model

    def method_line_from_control(include_freq: bool) -> str:
        tokens = [
            "Method:",
            str(config.get('functional','')).strip(),
            rel_token_ctl,
            str(main_basisset).strip(),
            str(config.get('disp_corr','')).strip(),
            str(config.get('ri_jkx','')).strip(),
            aux_jk_token_ctl,
            implicit_token(),
            str(config.get('geom_opt_OCCUPIER','')).strip()
        ]
        if include_freq:
            tokens.append("FREQ")
        tokens.append(f"PAL{config.get('PAL','')}")
        return " ".join(t for t in tokens if t).replace("  ", " ").strip()

    # ----------------------- time formatting -----------------------------------
    hours, remainder = divmod(duration, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_format = f"{int(hours):02d} hours {int(minutes):02d} minutes {seconds:05.2f} seconds"

    # ----------------------- energy / selection setup --------------------------
    freq_flag = "FREQ" if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes' else ""
    parity = "is_even" if is_even else "is_odd"
    use_gibbs = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == 'yes'
    energy_label = "Final Gibbs free energy" if use_gibbs else "FINAL SINGLE POINT ENERGY"
    lowest_label = f"LOWEST {energy_label}:"

    raw_sel = str(config.get('occupier_selection', 'tolerance')).lower()
    sel = raw_sel.split('|')[0].strip()
    ap = config.get('approximate_spin_projection_APMethod')
    ap_str = str(ap) if ap not in (None, 0, '0') else 'none'

    if sel in {'rounding', 'round', 'gerundet', 'runden'}:
        method = 'rounding'
    elif sel in {'tolerance', 'toleranz', 'toleranzband', 'epsilon'}:
        method = 'tolerance'
    else:
        method = 'truncation'

    try:
        prec = int(config.get('occupier_precision', config.get('occupier_rounded_value', 6)))
    except (TypeError, ValueError):
        prec = 6
    prec = max(0, min(10, prec))

    try:
        epsilon = float(config.get('occupier_epsilon', 10.0**(-prec)))
        if not (epsilon > 0):
            raise ValueError
    except (TypeError, ValueError):
        epsilon = 10.0**(-prec)

    dev_max_raw = config.get('dev_max', None)
    try:
        dev_max = float(dev_max_raw) if dev_max_raw is not None else None
    except (TypeError, ValueError):
        dev_max = None

    #-----------------------------------------------------------------------------
    # Decision thresholds
    DEV_TINY        = 1e-3   # treat BS singlet with ~zero deviation as closed-shell
    DEV_SIMILARITY  = 0.15   # "similar deviation" band
    DEV_GOOD_MARGIN = 0.30   # "clearly smaller" margin
    DEV_HIGH        = 0.50   # high contamination regime

    # AF-override extra energy window in Hartree (user-configurable)
    EPS_AF = float(config.get('bs_override_window_h', 0.002))

    # ----------------------- caches for deviation and J(3) ---------------------
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
        if not is_bs_flag:
            return None
        if idx not in j3_cache:
            try:
                j3_cache[idx] = extract_last_J3(file_for_index(idx))
            except Exception:
                j3_cache[idx] = None
        return j3_cache[idx]

    def entry_record(idx: int) -> dict:
        return next((e for e in sequence if e["index"] == idx), {})

    def entry_is_bs(idx: int) -> bool:
        return bool(entry_record(idx).get("BS"))

    def entry_mult(idx: int) -> int:
        try:
            return int(entry_record(idx).get("m", 0))
        except Exception:
            return 0

    def effective_dev(idx: int) -> float:
        """
        Effective deviation used for ranking:
          - For non-BS singlets (m == 1) with no UHF block → 0.0 (closed-shell RKS).
          - Otherwise use parsed deviation; if missing → +inf.
        """
        dev = get_dev(idx)
        if dev is None:
            if (not entry_is_bs(idx)) and (entry_mult(idx) == 1):
                return 0.0
            return float("inf")
        return dev

    def is_pseudo_closed_shell(idx: int) -> bool:
        return entry_is_bs(idx) and (entry_mult(idx) == 1) and (effective_dev(idx) <= DEV_TINY)

    def within_dev_limit(idx: int) -> bool:
        if dev_max is None:
            return True
        d = get_dev(idx)
        return (d is None) or (d <= dev_max)

    # ----------------------- candidate set by energy policy --------------------
    valid_all = [(e["index"], f) for e, f in zip(sequence, fspe_values) if f is not None]
    valid = [pair for pair in valid_all if within_dev_limit(pair[0])]
    if not valid and dev_max is not None:
        valid = valid_all

    def energy_key(val: float) -> float:
        if method == "rounding":
            return round(val, prec)
        if method == "truncation":
            return truncate(val, prec)
        return val  # tolerance => raw energies

    min_fspe_index: Optional[int] = None
    min_fspe_value: Optional[float] = None

    if valid:
        # Build the energy band according to policy
        if method in ("rounding", "truncation"):
            processed = [(i, f, energy_key(f)) for i, f in valid]
            best_key = min(v for _, _, v in processed)
            cands = [(i, f) for i, f, v in processed if v == best_key]
            min_raw = min(f for _, f in valid)  # for AF-override window
        else:
            min_raw = min(f for _, f in valid)
            eps_eff = epsilon + 1e-12
            cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

        # --- AF-override when only a single candidate remains but is highly contaminated ---
        if len(cands) == 1:
            best_idx, best_val = cands[0]
            best_dev = effective_dev(best_idx)
            if best_dev >= DEV_HIGH:
                # search BS candidates within extra AF window and with similar high deviation
                bs_pool = []
                for i, f in valid:
                    if i == best_idx:
                        continue
                    if not entry_is_bs(i):
                        continue
                    d = effective_dev(i)
                    if d >= DEV_HIGH and abs(d - best_dev) <= DEV_SIMILARITY and (f - min_raw) <= EPS_AF + 1e-12:
                        bs_pool.append((i, f))
                if bs_pool:
                    # choose among BS by larger |J(3)|, then smaller deviation, then lowest index
                    def bs_pick_score(i):
                        j = get_j3(i, True)
                        mag = abs(j) if (j is not None) else -float("inf")
                        return (-mag, effective_dev(i), i)
                    picked = min(bs_pool, key=lambda p: bs_pick_score(p[0]))
                    cands = [picked]  # override the single contaminated candidate

        if len(cands) == 1:
            min_fspe_index, min_fspe_value = cands[0]
        else:
            # Split candidates; reclassify pseudo-closed-shell BS to non-BS
            raw_bs = [(i, f) for i, f in cands if entry_is_bs(i)]
            raw_nb = [(i, f) for i, f in cands if not entry_is_bs(i)]
            nb_cands = raw_nb + [(i, f) for i, f in raw_bs if is_pseudo_closed_shell(i)]
            bs_cands = [(i, f) for i, f in raw_bs if not is_pseudo_closed_shell(i)]

            def _min_eff_dev(pairs):
                vals = [effective_dev(i) for i, _ in pairs]
                return min(vals) if vals else float("inf")

            min_dev_bs = _min_eff_dev(bs_cands)
            min_dev_nb = _min_eff_dev(nb_cands)

            chosen: Optional[tuple[int, float]] = None

            if bs_cands and nb_cands:
                # A) Non-BS clearly cleaner and not in high regime → prefer non-BS
                if (min_dev_nb + DEV_GOOD_MARGIN < min_dev_bs) and (min_dev_nb < DEV_HIGH):
                    def nb_score(i):
                        return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                    chosen = min(nb_cands, key=lambda p: nb_score(p[0]))
                # B) Similar and both high → prefer BS; among BS prefer larger |J(3)|
                elif (abs(min_dev_nb - min_dev_bs) <= DEV_SIMILARITY) and (min_dev_nb > DEV_HIGH) and (min_dev_bs > DEV_HIGH):
                    def bs_score(i):
                        j = get_j3(i, True)
                        mag = abs(j) if (j is not None) else -float("inf")
                        return (-mag, effective_dev(i), i)
                    chosen = min(bs_cands, key=lambda p: bs_score(p[0]))
                else:
                    # C) Fallback: smallest deviation; tie → prefer non-BS; then lowest index
                    def fb_score(i):
                        return (effective_dev(i), 1 if entry_is_bs(i) else 0, i)
                    chosen = min(cands, key=lambda p: fb_score(p[0]))
            elif bs_cands:
                # Only BS → |J(3)| first, then deviation
                def bs_only_score(i):
                    j = get_j3(i, True)
                    mag = abs(j) if (j is not None) else -float("inf")
                    return (-mag, effective_dev(i), i)
                chosen = min(bs_cands, key=lambda p: bs_only_score(p[0]))
            else:
                # Only non-BS → smallest deviation; prefer genuine non-BS over pseudo
                def nb_only_score(i):
                    return (effective_dev(i), 1 if is_pseudo_closed_shell(i) else 0, i)
                chosen = min(nb_cands, key=lambda p: nb_only_score(p[0]))

            if chosen is not None:
                min_fspe_index, min_fspe_value = chosen

    # ----------------------- build FSPE table text -----------------------------
    def fmt_display(x: Optional[float]) -> str:
        if x is None:
            return "N/A"
        return f"{x:.{prec}f}" if method == 'rounding' else fmt_truncate(x, prec)

    fspe_lines = ""
    for entry, fspe in zip(sequence, fspe_values):
        idx = entry["index"]
        multiplicity = entry["m"]
        bs = entry.get("BS", "")
        is_bs_flag = bool(bs)

        fspe_value = (fmt_display(fspe) + " (H)") if fspe is not None else "Not a valid value"
        dev_val = get_dev(idx)
        dev_txt = fmt_display(dev_val)
        mark = " <-- PREFERRED VALUE" if (min_fspe_index is not None and idx == min_fspe_index) else ""
        bs_line = f", BrokenSym {bs}" if bs else ""

        fspe_lines += f"{energy_label} ({idx})   = {fspe_value}{mark}\n"
        fspe_lines += f"multiplicity {multiplicity}{bs_line}\n"
        fspe_lines += f"Spin Contamination (⟨S²⟩ - S(S+1))   : {dev_txt}\n"

        if is_bs_flag:
            j3_val = get_j3(idx, True)
            j3_str = "N/A" if j3_val is None else f"{j3_val:.2f} cm⁻¹"
            fspe_lines += f"J(Yamaguchi)=-(E[HS]-E[BS])/(⟨S²⟩HS-⟨S²⟩BS) : {j3_str}\n"

        fspe_lines += "----------------------------------------------------------------\n"

    print("\n" + fspe_lines)

    # ----------------------- compose method & metal-basis lines ----------------
    include_freq = (freq_flag == 'FREQ')
    pal = config.get('PAL', '')

    method_str = None
    metal_basis_from_inp = None

    prefer_idx = (min_fspe_index if min_fspe_index is not None else 1)
    inp_pref = _inp_for_index(prefer_idx) or _inp_for_index(1)

    if inp_pref:
        parsed = _parse_method_from_inp(inp_pref, pal, include_freq)
        if parsed:
            method_str = parsed
        metal_basis_from_inp = _parse_metal_basis_from_inp(inp_pref)

    if not method_str:
        method_str = method_line_from_control(include_freq)

    metal_basis_print = metal_basis_from_inp or metal_basisset_rep_ctl or ""

    lowest_str = (
        f"{fmt_truncate(min_fspe_value, prec)} (H)" if min_fspe_value is not None
        else "No valid FSPE values found"
    )

    # ---- write report ----------------------------------------------------------
    banner = build_occupier_banner(header_indent=6, info_indent=6)

    with open('OCCUPIER.txt', 'w', encoding='utf-8') as file:
        file.write(
            f"{banner}\n\n"
            f"{method_str}\n"
            f"        {', '.join(metals)} {metal_basis_print}\n\n"
            f"Charge: {charge}     \n"
            f"-------------  \n"
            f"{fspe_lines}\n\n"
            f"TOTAL RUN TIME: {duration_format}\n\n"
            f"{lowest_label} {lowest_str}\n\n"
            f"(Selection: {method}, APmethod {ap_str}, dev_max {dev_max if dev_max is not None else 'none'}, precision={prec}, epsilon={epsilon if method=='tolerance' else 'n/a'})\n"
            f"(Preferred Index: {min_fspe_index if min_fspe_index is not None else 'N/A'})\n"
            f"(Electron number: {parity})\n"
        )
