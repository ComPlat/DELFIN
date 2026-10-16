"""Index DELFIN calculation directories for structured search.

Scans ``calc_dir``, ``archive_dir``, and ``remote_archive_dir`` for
DELFIN calculations.  Each calculation is a directory that contains at
least one of: ``DELFIN_Data.json``, ``CONTROL.txt``, or ORCA ``.inp``
files.

The resulting index is a flat list of records with normalised metadata
(method, basis set, solvent, charge, energies, modules enabled, …)
that can be searched by keyword or structured query.

**Read-only**: this module never writes to calc/archive directories.
"""

from __future__ import annotations

import re

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Per-calc extraction
# ---------------------------------------------------------------------------

def _extract_from_delfin_data(data: dict) -> dict[str, Any]:
    """Extract searchable fields from a DELFIN_Data.json dict."""
    rec: dict[str, Any] = {}

    # -- metadata block (archive calcs) ------------------------------------
    # DELFIN writes the method under "metadata"; the benchmark fixture and
    # older files write it flat. Both are read, or a flat file files its
    # run under no method at all.
    meta = data.get("metadata") if isinstance(data.get("metadata"), dict) else {}
    for key in ("functional", "basis_set", "solvent", "charge"):
        if key not in meta and key in data:
            meta = {**meta, key: data.get(key)}
    if isinstance(meta, dict):
        rec["functional"] = meta.get("functional", "")
        rec["basis_set"] = meta.get("basis_set", "")
        rec["aux_basis"] = meta.get("auxiliary_basis", "")
        rec["ri_method"] = meta.get("ri_method", "")
        rec["dispersion"] = meta.get("dispersion_correction", "")
        rec["solvation"] = meta.get("implicit_solvation", "")
        rec["solvent"] = meta.get("solvent", "")
        rec["charge"] = meta.get("charge", "")

    # -- input block -------------------------------------------------------
    inp = data.get("input", {})
    if isinstance(inp, dict):
        rec["smiles"] = inp.get("smiles", "")

    # -- control.parsed (always present) -----------------------------------
    ctrl = data.get("control", {})
    parsed = ctrl.get("parsed", {}) if isinstance(ctrl, dict) else {}
    if isinstance(parsed, dict):
        if not rec.get("smiles"):
            rec["smiles"] = parsed.get("SMILES", "")
        rec["name"] = parsed.get("NAME", "")
        rec["charge"] = rec.get("charge") or parsed.get("charge", "")
        rec["solvent"] = rec.get("solvent") or parsed.get("solvent", "")
        rec["solvation"] = rec.get("solvation") or parsed.get(
            "implicit_solvation_model", ""
        )
        rec["xtb_method"] = parsed.get("xTB_method", "")

        # Modules enabled
        modules = []
        if _is_yes(parsed.get("IMAG")):
            modules.append("IMAG")
        if _is_yes(parsed.get("ESD_modul")):
            modules.append("ESD")
            rec["esd_modus"] = parsed.get("ESD_modus", "")
        if _is_yes(parsed.get("calc_prop_of_interest")):
            modules.append("properties")
            rec["properties_of_interest"] = parsed.get(
                "properties_of_interest", ""
            )
        if parsed.get("oxidation_steps"):
            modules.append("oxidation")
            rec["oxidation_steps"] = parsed.get("oxidation_steps", "")
        if parsed.get("reduction_steps"):
            modules.append("reduction")
            rec["reduction_steps"] = parsed.get("reduction_steps", "")
        parsed_optimizer = str(parsed.get("global_optimizer", "")).strip().upper()
        if _is_yes(parsed.get("XTB_GOAT")) or parsed_optimizer == "GOAT":
            modules.append("GOAT")
        if _is_yes(parsed.get("CREST")) or parsed_optimizer == "CREST":
            modules.append("CREST")
        # The archive holds 126 calculations tagged GUPPY; that tag stays
        # readable so they do not fall out of the index.  A run configured the
        # new way is tagged MANTA, and a legacy GUPPY=yes run is both -- it is
        # the same builder either way.
        if _is_yes(parsed.get("GUPPY")):
            modules.append("GUPPY")
        if str(parsed.get("smiles_converter", "")).strip().upper() in ("MANTA", "GUPPY") \
                or _is_yes(parsed.get("GUPPY")):
            modules.append("MANTA")
        if _is_yes(_flag(parsed, "hyperpol_xTB", "hyperpol_xtb_module")):
            modules.append("hyperpol_xtb")
        if _is_yes(_flag(parsed, "tadf_xTB", "TADF_xTB_module")):
            modules.append("TADF_xTB")
        if _is_yes(_flag(parsed, "NMR_module")):
            modules.append("NMR")
        if _is_yes(_flag(parsed, "OCCUPIER")) or _flag(parsed, "method").strip().upper() == "OCCUPIER":
            modules.append("OCCUPIER")
        if _is_yes(_flag(parsed, "co2_coordination", "CO2_coordinator")):
            modules.append("CO2")
        rec["modules"] = modules

        # ORCA-level settings from CONTROL
        for key in (
            "method", "method_rel", "basis_set", "basis_set_rel",
            "auxiliary_basis_set", "auxiliary_basis_set_rel",
            "RI_method", "dispersion_correction",
        ):
            val = parsed.get(key, "")
            if val and isinstance(val, str) and val.strip():
                ctrl_key = f"ctrl_{key}"
                rec[ctrl_key] = val.strip()
        # CONTROL.txt "method" = DELFIN workflow (classic/OCCUPIER/…), not DFT
        rec["workflow_method"] = parsed.get("method", "")
        # Fill basis from CONTROL if metadata was empty
        if not rec.get("basis_set"):
            rec["basis_set"] = parsed.get("basis_set", "")

    # -- energies from ground_state_S0 -------------------------------------
    gs = data.get("ground_state_S0", {})
    if isinstance(gs, dict):
        thermo = gs.get("thermochemistry", {})
        if isinstance(thermo, dict):
            rec["gibbs_energy"] = thermo.get("gibbs_free_energy")
            rec["electronic_energy"] = thermo.get("electronic_energy")
            rec["zpe"] = thermo.get("zero_point_energy")

    # -- summary -----------------------------------------------------------
    summary = data.get("delfin_summary", {})
    if isinstance(summary, dict):
        rt = summary.get("total_run_time", {})
        if isinstance(rt, dict):
            rec["run_time_seconds"] = rt.get("total_seconds")

    # -- which result sections have data -----------------------------------
    has_data = []
    for key in (
        "ground_state_S0", "excited_states", "oxidized_state",
        "oxidized_states", "reduced_state", "reduced_states",
        "vibrational_frequencies", "emission", "intersystem_crossing",
        "internal_conversion", "fluorescence_rates",
        "phosphorescence_rates", "properties_of_interest",
        "reorganisation_energy", "occupier", "photophysical_rates",
        "esd_results", "hyperpol_xtb", "tadf_xtb",
    ):
        val = data.get(key, {})
        if isinstance(val, dict) and val:
            has_data.append(key)
    rec["has_data"] = has_data

    return rec


def latest_exit_code_marker(d):
    """The newest ``.exit_code_*`` file in a folder, or None.

    A resubmitted calculation leaves one marker per run; the last run
    decides how the folder ended.
    """
    from pathlib import Path as _P
    markers = list(_P(d).glob(".exit_code_*"))
    if not markers:
        return None
    try:
        return max(markers, key=lambda p: p.stat().st_mtime)
    except OSError:
        return markers[0]


def exit_code_of_marker(marker):
    """The exit code a ``.exit_code_<job id>`` marker records.

    The runner writes the code as the file's CONTENT and the scheduler's
    job id into the name (dashboard/local_runner.py), so the name's number
    is a job id, not a code. Only an empty marker — older folders and the
    benchmark fixtures — falls back to the number in the name.
    """
    try:
        text = marker.read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return "unknown"
    if text:
        return int(text) if text.lstrip("-").isdigit() else "unknown"
    code = marker.name.rsplit("_", 1)[-1]
    return int(code) if code.isdigit() else "unknown"


def completion_of(d) -> tuple:
    """(completed, exit_code) for a calculation folder, from what it left.

    An exit-code file means the run ENDED (the code says how); a run log
    without one means running or crashed; neither means unknown. One
    reading, shared by the index and the energy tools, so the two never
    call the same folder two different things.
    """
    from pathlib import Path as _P
    d = _P(d)
    marker = latest_exit_code_marker(d)
    if marker is not None:
        return True, exit_code_of_marker(marker)
    if (d / "delfin_run.log").is_file():
        return False, None
    return None, None


def outcome_of_folder(d) -> str:
    """The outcome phrase for a folder, in one call.

    When neither an exit-code file nor a run log exists, the pipeline
    state file still may: a status of running / finished there is
    evidence, and "unknown" beside a file that says "finished" was the
    answer a model called misleading.
    """
    import json as _json
    from pathlib import Path as _P
    completed, exit_code = completion_of(d)
    if completed is not None:
        return _outcome_of(completed, exit_code)
    d = _P(d)
    for name in ("DELFIN_Data.json", "DELFIN_data.json"):
        f = d / name
        if not f.is_file():
            continue
        try:
            status = str((_json.loads(f.read_text(encoding="utf-8")) or {}).get("status") or "").strip().lower()
        except Exception:
            break
        if status in ("running", "active", "started", "in_progress"):
            return f"running per {name} (no exit code yet)"
        if status in ("finished", "completed", "done", "ok", "success"):
            return f"finished per {name} (no exit code file)"
        if status in ("failed", "error", "crashed"):
            return f"failed per {name} (no exit code file)"
        break
    # An input with no output and no log: not started, or running on a
    # machine that writes no log here. Either way it is not "unknown"
    # in the sense of "nothing to go on" -- there is an input waiting.
    if any(d.glob("*.inp")) and not any(d.glob("*.out")):
        return "no output yet (input present; not started or still running)"
    # An output with no bookkeeping around it -- an archived run copied
    # without its state file, a calculation run by hand -- still says how
    # it ended: ORCA writes its own termination line. Read from the tail,
    # where it is; an output of a million lines is not read to find it.
    # Driven 2026-09-10 over such a folder: "unknown" beside an output
    # that says TERMINATED NORMALLY.
    out = _largest_output(d)
    if out is not None:
        tail = _tail_text(out)
        if "ORCA TERMINATED NORMALLY" in tail:
            return "finished per ORCA output (no exit code file)"
        if "ORCA finished by error termination" in tail or "ABORTING THE RUN" in tail:
            return "failed per ORCA output (error termination)"
        return "running or crashed (output has no termination line)"
    return _outcome_of(completed, exit_code)


def _largest_output(d):
    """The .out most likely to be the run's own: the largest."""
    from pathlib import Path as _P
    best = None
    try:
        for f in _P(d).glob("*.out"):
            try:
                size = f.stat().st_size
            except OSError:
                continue
            if best is None or size > best[0]:
                best = (size, f)
    except OSError:
        return None
    return best[1] if best else None


def _tail_text(path, size: int = 16384) -> str:
    """The last ``size`` bytes of a file as text, never raising."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, 2)
            n = fh.tell()
            fh.seek(max(0, n - size))
            return fh.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def method_of_folder(d) -> tuple:
    """(functional, basis) from what the folder holds besides the output.

    An ORCA output need not state its method -- a DELFIN run's does not
    -- so a parser that reads only the .out returns "", and every tool
    downstream then files the run under method None. The folder knows:
    DELFIN_Data.json first, then CONTROL.txt, then the .inp header, the
    order the index itself uses. See :func:`method_parts_of_folder` for
    the solvent and dispersion that complete the method.
    """
    parts = method_parts_of_folder(d)
    return parts["functional"], parts["basis"]


def _outcome_of(completed, exit_code) -> str:
    """One phrase a reader cannot mistake: did the run succeed?"""
    if completed is None:
        return "unknown (no run log and no exit code)"
    if completed is False:
        return "running or crashed (run log present, no exit code)"
    if exit_code == 0:
        return "succeeded (exit code 0)"
    if isinstance(exit_code, int):
        return f"failed (exit code {exit_code})"
    return "ended with an unreadable exit code"



# The CONTROL.txt keys read below, as delfin.define.TEMPLATE writes them
# (pinned against TEMPLATE in
# tests/test_the_calc_index_reads_the_control_template_keys.py).
CONTROL_TEMPLATE_KEYS = (
    "input_file", "NAME", "SMILES", "charge", "multiplicity_0",
    "implicit_solvation_model", "solvent", "functional", "disp_corr",
    "relativity", "main_basisset", "main_basisset_rel",
    "metal_basisset", "metal_basisset_rel",
)


# A geometry file larger than this is not read: 16 MiB is ~300k atoms of
# XYZ, far above any DFT input, and bounds the scan on a stray dump.
_GEOMETRY_MAX_BYTES = 16 * 1024 * 1024

# First line of a geometry: an XYZ atom count, or an atom line (element
# token, then a number). DELFIN's input.txt has no atom count
# (define.convert_xyz_to_input_txt drops it), so the count alone misses it.
_GEOMETRY_FIRST_LINE = re.compile(r"^(\d+$|[A-Za-z]{1,2}\S*\s+[-+]?\.?\d)")


def _metals_in(folder: Path | None, kv: dict) -> list[str] | None:
    """The transition metals of the structure, as delfin.utils.
    search_transition_metals reads them from a geometry (input.txt
    coordinates or XYZ), or from the bracket atoms of a SMILES; None when
    no structure is readable (no SMILES in CONTROL, no geometry file
    inside *folder*, or one above _GEOMETRY_MAX_BYTES)."""
    from delfin.utils import _TM_LIST, search_transition_metals
    smiles = [kv.get("SMILES", "")]
    if folder is not None:
        geom = folder / (kv.get("input_file") or "input.txt")
        try:
            if geom.is_file() and geom.resolve().is_relative_to(folder.resolve()):
                if geom.stat().st_size > _GEOMETRY_MAX_BYTES:
                    return None
                lines = [ln.strip() for ln in geom.read_text(
                    encoding="utf-8", errors="replace").splitlines()
                    if ln.strip() and not ln.strip().startswith("*")]
                if lines and _GEOMETRY_FIRST_LINE.match(lines[0]):
                    return search_transition_metals(str(geom))
                if lines:
                    smiles.append(lines[0])
        except OSError:
            pass
    smiles = [s for s in smiles if s.strip()]
    if not smiles:
        return None
    return sorted({m for s in smiles
                   for m in re.findall(r"\[\d*([A-Z][a-z]?)", s) if m in _TM_LIST})


def _extract_from_control_txt(text: str, folder: Path | None = None) -> dict[str, Any]:
    """Extract searchable fields from a CONTROL.txt file.

    Level of theory from the define.TEMPLATE keys: ``functional``,
    ``disp_corr``, ``main_basisset``, ``metal_basisset``; for a structure
    with a 4d/5d/4f/5f metal the ``*_rel`` basis values and the
    ``relativity`` Hamiltonian, chosen as delfin.utils.set_main_basisset
    and select_rel_and_aux choose them (``relativity`` alone does not
    switch the basis). The atoms come from ``SMILES`` or the geometry file
    ``input_file`` inside *folder*; with neither, the template basis stays
    empty and the .inp header supplies it. The non-template ``basis_set``
    key remains the fallback (the benchmark fixtures write it).
    """
    rec: dict[str, Any] = {}
    kv: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if "=" in line and not line.startswith("#") and not line.startswith("-"):
            key, _, val = line.partition("=")
            kv[key.strip()] = val.strip()

    rec["name"] = kv.get("NAME", "")
    rec["smiles"] = kv.get("SMILES", "")
    rec["charge"] = kv.get("charge", "")
    rec["solvent"] = kv.get("solvent", "")
    rec["solvation"] = kv.get("implicit_solvation_model", "")
    rec["xtb_method"] = kv.get("xTB_method", "")
    # CONTROL.txt "method" is the DELFIN workflow method (classic/OCCUPIER/...),
    # NOT the DFT functional.  We store it separately.
    rec["workflow_method"] = kv.get("method", "")
    rec["functional"] = kv.get("functional", "")
    rec["dispersion"] = _clean_dispersion(kv.get("disp_corr")) or ""
    rec["multiplicity"] = kv.get("multiplicity_0", "")
    metals = _metals_in(folder, kv)
    main, metal, relativity = "", "", ""
    if metals is not None and ("main_basisset" in kv or "main_basisset_rel" in kv):
        # DELFIN's own choice, defaults for an empty key included
        # (def2-SVP / ZORA-def2-SVP); a CONTROL without the template keys
        # keeps the basis_set fallback below.
        from delfin.utils import select_rel_and_aux, set_main_basisset
        main, metal = set_main_basisset(metals, kv)
        metal = metal or ""
        relativity = select_rel_and_aux(metals, kv)[0]
    rec["basis_set"] = main or kv.get("basis_set", "")
    rec["metal_basis"] = metal
    rec["relativity"] = relativity

    modules = []
    if _is_yes(kv.get("IMAG")):
        modules.append("IMAG")
    if _is_yes(kv.get("ESD_modul")):
        modules.append("ESD")
        rec["esd_modus"] = kv.get("ESD_modus", "")
    if _is_yes(kv.get("GUPPY")):
        modules.append("GUPPY")
    if str(kv.get("smiles_converter", "")).strip().upper() in ("MANTA", "GUPPY") \
            or _is_yes(kv.get("GUPPY")):
        modules.append("MANTA")
    global_optimizer = str(kv.get("global_optimizer", "")).strip().upper()
    if _is_yes(kv.get("XTB_GOAT")) or global_optimizer == "GOAT":
        modules.append("GOAT")
    if _is_yes(kv.get("CREST")) or global_optimizer == "CREST":
        modules.append("CREST")
    if kv.get("oxidation_steps"):
        modules.append("oxidation")
    if kv.get("reduction_steps"):
        modules.append("reduction")
    if _is_yes(_flag(kv, "hyperpol_xTB", "hyperpol_xtb_module")):
        modules.append("hyperpol_xtb")
    if _is_yes(_flag(kv, "tadf_xTB", "TADF_xTB_module")):
        modules.append("TADF_xTB")
    if _is_yes(_flag(kv, "NMR_module")):
        modules.append("NMR")
    if _is_yes(_flag(kv, "OCCUPIER")) or _flag(kv, "method").strip().upper() == "OCCUPIER":
        modules.append("OCCUPIER")
    if _is_yes(_flag(kv, "co2_coordination", "CO2_coordinator")):
        modules.append("CO2")
    rec["modules"] = modules

    return rec


_KNOWN_FUNCTIONALS = {
    "hf", "rhf", "uhf",
    "b3lyp", "pbe", "pbe0", "bp86", "tpss", "tpssh", "m06", "m06-2x", "m06-l",
    "cam-b3lyp", "wb97x", "wb97x-d3", "wb97x-d3bj", "wb97x-v",
    "wb97m-v", "wb97m-d3bj", "wb97x-3c",
    "b2plyp", "ri-b2plyp", "dlpno-ccsd(t)", "ccsd(t)",
    "r2scan", "r2scan-3c", "b97-3c",
    "revpbe", "blyp", "pw6b95",
    "b97m-v", "b97m-d3bj", "b97m-d4", "wb97x-d", "wb97x-d4", "wb97m-d4",
    "pbeh-3c", "hf-3c",
}
# Dispersion written into the functional token (B3LYP-D3BJ, PBE0-D4).
_KNOWN_FUNCTIONALS |= {
    f"{base}-{disp}"
    for base in ("b3lyp", "pbe", "pbe0", "bp86", "blyp", "tpss", "tpssh",
                 "revpbe", "r2scan", "b2plyp", "cam-b3lyp", "m06-2x", "pw6b95")
    for disp in ("d3", "d3bj", "d3zero", "d4")
}
# Semiempirical tight-binding methods; ORCA's keywords and the xtb names.
# A header naming one of these is DELFIN's pre-optimisation (XTB.inp,
# XTB_GOAT.inp), not the run's DFT method: it fills a functional nothing
# else states and replaces none.
_XTB_METHODS = frozenset({"xtb", "xtb0", "xtb1", "xtb2", "gfn0-xtb",
                          "gfn1-xtb", "gfn2-xtb", "gfn-ff", "gfnff"})
_KNOWN_FUNCTIONALS |= _XTB_METHODS


def _is_xtb(header: dict) -> bool:
    """Whether an .inp header's method is a tight-binding one (_XTB_METHODS)."""
    return str(header.get("inp_functional", "")).lower() in _XTB_METHODS

_KNOWN_BASIS = {
    "def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvp", "def2-qzvpp",
    "ma-def2-svp", "ma-def2-tzvp", "ma-def2-tzvpp",
    "cc-pvdz", "cc-pvtz", "cc-pvqz", "aug-cc-pvdz", "aug-cc-pvtz",
    "6-31g*", "6-31g**", "6-311g**", "6-31+g*", "6-311+g(2d,p)",
    "sto-3g", "def2-svp(d)",
    "sarc-dkh-tzvp", "sarc-zora-tzvp", "x2c-tzvp", "x2c-svp",
    "zora-def2-svp", "zora-def2-tzvp", "zora-def2-tzvpp",
    "vdzp",
}


def _extract_from_inp_header(text: str) -> dict[str, Any]:
    """Extract method/basis from the first ``!`` line of an ORCA .inp file."""
    rec: dict[str, Any] = {}
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("!"):
            keywords = line[1:].split()
            rec["orca_keywords"] = keywords
            # Try to identify functional and basis set
            for kw in keywords:
                kw_lower = kw.lower()
                if kw_lower in _KNOWN_FUNCTIONALS and "inp_functional" not in rec:
                    rec["inp_functional"] = kw
                if kw_lower in _KNOWN_BASIS and "inp_basis" not in rec:
                    rec["inp_basis"] = kw
                # Solvation and dispersion ride on the same line -- CPCM(DMF),
                # SMD(water), D3BJ -- and they decide comparability as much
                # as the functional does.
                m = re.match(r"(cpcm|smd|cpcmc|alpb|ddcosmo)(?:\((.+)\))?$", kw_lower)
                if m and "inp_solvation" not in rec:
                    rec["inp_solvation"] = kw.split("(", 1)[0].upper()
                    if m.group(2):
                        rec["inp_solvent"] = kw[kw.index("(") + 1:-1]
                if kw_lower in _KNOWN_DISPERSION and "inp_dispersion" not in rec:
                    rec["inp_dispersion"] = kw.upper()
            break
    return rec


# Record field -> .inp header field that, when present, overrides the
# CONTROL.txt value (not a DELFIN_Data.json one).
_INP_STATED = {"functional": "inp_functional", "basis_set": "inp_basis",
               "dispersion": "inp_dispersion"}

_KNOWN_DISPERSION = frozenset({"d2", "d3", "d3bj", "d3zero", "d4", "d3(bj)", "d3(zero)"})
_NO_SOLVENT = frozenset({"", "none", "no", "gas", "gasphase", "gas phase", "vacuum", "vac", "false", "off"})


def _clean_solvent(value) -> str | None:
    v = str(value or "").strip()
    return None if v.lower() in _NO_SOLVENT else v


def _clean_dispersion(value) -> str | None:
    v = str(value or "").strip()
    if v.lower() in ("", "none", "no", "false", "off"):
        return None
    return v.upper()


def method_parts_of_folder(d) -> dict:
    """functional, basis, dispersion and solvent from what the folder holds.

    The four a total energy is only comparable within. Read in the order
    the index uses: DELFIN_Data.json first; functional, basis and
    dispersion then from the .inp header (what ran) before CONTROL.txt;
    the solvent from CONTROL.txt before the .inp. A solvent of
    "none" is the gas phase and comes back as None, so two gas-phase
    runs compare and a CPCM run compares with neither.
    """
    import json as _json
    from pathlib import Path as _P
    d = _P(d)
    parts: dict = {"functional": None, "basis": None, "solvent": None,
                   "dispersion": None, "solvation": None}
    seen_solvent = False
    ctrl: dict = {}
    xtb = None  # a tight-binding header's method: fills only, after CONTROL
    try:
        if (d / "DELFIN_Data.json").is_file():
            rec = _extract_from_delfin_data(_json.loads((d / "DELFIN_Data.json").read_text(encoding="utf-8")))
            parts["functional"] = rec.get("functional") or None
            parts["basis"] = rec.get("basis_set") or None
            parts["dispersion"] = _clean_dispersion(rec.get("dispersion"))
            if "solvent" in rec:
                parts["solvent"] = _clean_solvent(rec.get("solvent"))
                seen_solvent = True
            parts["solvation"] = str(rec.get("solvation") or "").strip() or None
    except Exception:
        pass
    try:
        if (d / "CONTROL.txt").is_file():
            text = (d / "CONTROL.txt").read_text(encoding="utf-8")
            rec = _extract_from_control_txt(text, d)
            kv: dict = {}
            for line in text.splitlines():
                line = line.strip()
                if "=" in line and not line.startswith("#") and not line.startswith("-"):
                    k, _, v = line.partition("=")
                    kv[k.strip()] = v.strip()
            ctrl = {"basis": rec.get("basis_set") or None,
                    "functional": rec.get("functional") or None,
                    "dispersion": (rec.get("dispersion")
                                   or _clean_dispersion(kv.get("dispersion_correction")))}
            if not seen_solvent and "solvent" in kv:
                parts["solvent"] = _clean_solvent(kv.get("solvent"))
                seen_solvent = True
            parts["solvation"] = parts["solvation"] or (rec.get("solvation") or None)
    except Exception:
        pass
    try:
        if not parts["functional"] or not parts["basis"] or not seen_solvent or not parts["dispersion"]:
            for inp in sorted(d.glob("*.inp")):
                rec = _extract_from_inp_header(inp.read_text(encoding="utf-8", errors="replace"))
                if _is_xtb(rec):
                    xtb = xtb or rec["inp_functional"]
                else:
                    parts["functional"] = parts["functional"] or rec.get("inp_functional") or None
                parts["basis"] = parts["basis"] or rec.get("inp_basis") or None
                if not seen_solvent and ("inp_solvent" in rec or "inp_solvation" in rec):
                    parts["solvent"] = _clean_solvent(rec.get("inp_solvent"))
                    parts["solvation"] = parts["solvation"] or rec.get("inp_solvation")
                    seen_solvent = True
                parts["dispersion"] = parts["dispersion"] or _clean_dispersion(rec.get("inp_dispersion"))
                if parts["functional"] and parts["basis"] and seen_solvent:
                    break
    except Exception:
        pass
    for key, value in ctrl.items():
        parts[key] = parts[key] or value
    parts["functional"] = parts["functional"] or xtb
    return parts


def _is_yes(val: Any) -> bool:
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        # "on" is what define.py documents for co2_coordination.
        return val.strip().lower() in ("yes", "true", "1", "on")
    return False


def _norm_key(raw: str) -> str:
    """A CONTROL key the way config.py compares it: lowercase, letters and
    digits only. ``co2_coordination``, ``Co2-Coordination`` and
    ``CO2 coordination`` are one key there, and must be one here."""
    return re.sub(r"[^a-z0-9]+", "", str(raw or "").strip().lower())


def _flag(kv: dict, *names: str) -> str:
    """The value of the first of *names* present in ``kv``, compared by
    normalised key, or "".

    The module flags used to be read by exact, case-sensitive names that
    CONTROL never writes: ``CO2_coordinator`` for the key that is
    ``co2_coordination``, ``TADF_xTB_module`` for ``tadf_xTB``,
    ``hyperpol_xtb_module`` for ``hyperpol_xTB``. No indexed run ever
    carried those modules (found by an operator interview, 2026-09-11).
    The legacy spellings stay accepted; the canonical ones come first.
    """
    normed = {_norm_key(k): v for k, v in (kv or {}).items()}
    for name in names:
        v = normed.get(_norm_key(name))
        if v is not None and str(v).strip():
            return str(v)
    return ""



# ---------------------------------------------------------------------------
# Directory scanner
# ---------------------------------------------------------------------------

def _scan_calc_dir(
    calc_dir: Path,
    source: str,
    quiet: bool = False,
    truncated: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Scan a directory for DELFIN calculations.

    Each immediate subdirectory that contains CONTROL.txt, DELFIN_Data.json,
    or .inp files is treated as a calculation.  Also recurses one level for
    archive structures like ``archive/<user>/calc_name/``.

    Directories still holding subdirectories at the depth limit are appended
    to *truncated* — a calculation deeper than the limit is silently outside
    every count otherwise, and ``calc`` is scanned one level deep while the
    archives get three.
    """
    records: list[dict[str, Any]] = []
    if not calc_dir.is_dir():
        return records

    visited: set[str] = set()

    def _process_dir(d: Path, rel_prefix: str = "") -> None:
        key = str(d.resolve())
        if key in visited:
            return
        visited.add(key)

        has_control = (d / "CONTROL.txt").is_file()
        has_data = (d / "DELFIN_Data.json").is_file()
        has_inp = any(d.glob("*.inp"))

        if not (has_control or has_data or has_inp):
            return

        rec: dict[str, Any] = {
            "calc_id": d.name,
            "path": str(d),
            "rel_path": f"{rel_prefix}{d.name}" if rel_prefix else d.name,
            "source": source,
        }

        # Priority: DELFIN_Data.json > CONTROL.txt > .inp files
        if has_data:
            try:
                data = json.loads(
                    (d / "DELFIN_Data.json").read_text(encoding="utf-8")
                )
                rec.update(_extract_from_delfin_data(data))
                rec["index_source"] = "DELFIN_Data.json"
            except Exception as exc:
                if not quiet:
                    print(f"  [warn] {d.name}/DELFIN_Data.json: {exc}",
                          file=sys.stderr)

        from_json = {k for k in _INP_STATED if rec.get(k)}

        if has_control and "index_source" not in rec:
            try:
                text = (d / "CONTROL.txt").read_text(encoding="utf-8")
                rec.update(_extract_from_control_txt(text, d))
                rec["index_source"] = "CONTROL.txt"
            except Exception:
                pass
        elif has_control and "functional" not in rec:
            # Supplement from CONTROL if DELFIN_Data didn't have metadata
            try:
                text = (d / "CONTROL.txt").read_text(encoding="utf-8")
                ctrl_rec = _extract_from_control_txt(text, d)
                for k, v in ctrl_rec.items():
                    if k not in rec or not rec[k]:
                        rec[k] = v
            except Exception:
                pass

        if has_inp and not rec.get("orca_keywords"):
            # ORCA keywords from the first .inp whose header names a DFT
            # functional (XTB.inp, XTB_GOAT.inp sort before initial.inp),
            # else from the first .inp, as method_parts_of_folder reads them.
            first = None
            for inp_file in sorted(d.glob("*.inp")):
                try:
                    header = _extract_from_inp_header(
                        inp_file.read_text(encoding="utf-8", errors="replace")[:500])
                except Exception:
                    continue
                first = first or header
                if header.get("inp_functional") and not _is_xtb(header):
                    first = header
                    break
            rec.update(first or {})

        # What the .inp header states wins over CONTROL.txt: DELFIN
        # rewrites every .inp from CONTROL at run time, so they differ only
        # when CONTROL was edited after the run, and the .inp is then what
        # produced the outputs here. DELFIN_Data.json still wins.
        for field, inp_key in _INP_STATED.items():
            if rec.get(inp_key) and field not in from_json and not (
                    _is_xtb(rec) and rec.get(field)):
                rec[field] = rec[inp_key]

        # List .out files present
        out_files = sorted(p.name for p in d.glob("*.out")
                           if not p.name.startswith("delfin_"))
        rec["out_files"] = out_files

        # Check completion status.
        #
        # `completed` means the run ENDED -- an exit-code file exists --
        # not that it succeeded; a run that died with exit code 1025 is
        # "completed: true" here. Driven through get_calc_info, that
        # record read as a success. `outcome` says which it was, in
        # words, so the two are never read as one.
        completed, exit_code = completion_of(d)
        rec["completed"] = completed
        if completed:
            rec["exit_code"] = exit_code
        rec["outcome"] = _outcome_of(completed, exit_code)

        records.append(rec)

    def _recurse(d: Path, prefix: str, depth: int, max_depth: int) -> None:
        """Recursively scan for calc dirs up to *max_depth* levels."""
        if depth > max_depth:
            return
        try:
            children = sorted(d.iterdir())
        except OSError:
            return
        for child in children:
            if not child.is_dir() or child.name.startswith("."):
                continue
            _process_dir(child, rel_prefix=prefix)
            # Go deeper if this wasn't a calc dir (i.e. it's a grouping folder)
            child_key = str(child.resolve())
            if child_key not in visited or depth < max_depth:
                _recurse(
                    child,
                    prefix=f"{prefix}{child.name}/",
                    depth=depth + 1,
                    max_depth=max_depth,
                )
            elif truncated is not None and _has_subdirs(child):
                # The scan stops here and there is still tree below it. A
                # calculation deeper than the limit is outside every count
                # otherwise, silently, and calc/ has a limit of one level
                # against the archives' three.
                truncated.append(f"{prefix}{child.name}")

    _recurse(calc_dir, "", 0, scan_depth_for(source))

    return records


def scan_depth_for(source: str) -> int:
    """How many levels below a root are scanned. calc/ is flat, archives nest."""
    return 3 if source in ("archive", "remote_archive") else 1


def _has_subdirs(d: Path) -> bool:
    try:
        return any(c.is_dir() and not c.name.startswith(".")
                   for c in d.iterdir())
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Full index builder
# ---------------------------------------------------------------------------

def build_calc_index(
    calc_dir: Path | None = None,
    archive_dir: Path | None = None,
    remote_archive_dir: Path | None = None,
    quiet: bool = False,
) -> dict[str, Any]:
    """Build a searchable index of all DELFIN calculations.

    Parameters
    ----------
    calc_dir : Path, optional
        Active calculations directory (default ``~/calc``).
    archive_dir : Path, optional
        Local archive (default ``~/archive``).
    remote_archive_dir : Path, optional
        Remote/shared archive (if configured).
    quiet : bool
        Suppress progress output.

    Returns
    -------
    dict
        Index with ``records`` list, metadata, and ``provenance`` — which
        roots were scanned, which were missing, and where the scan was cut
        off by its depth limit. Without it a home with no ``calc``
        directory produced a SUCCESS with zero records, so "you have no
        calculations" and "I scanned nothing" were the same answer: the
        explicit "index could not be built" error was unreachable,
        ``search_calcs`` returned ``[]`` and ``calc_summary`` returned all
        zeros, exactly as for a populated home whose query simply missed.
    """
    all_records: list[dict[str, Any]] = []
    roots: list[dict[str, Any]] = []

    def _scan(directory: Path | None, source: str) -> None:
        if directory is None:
            roots.append({"source": source, "path": "", "status": "not_configured",
                          "records": 0, "max_depth": scan_depth_for(source)})
            return
        if not directory.is_dir():
            roots.append({"source": source, "path": str(directory),
                          "status": "missing", "records": 0,
                          "max_depth": scan_depth_for(source)})
            return
        if not quiet:
            print(f"  scanning {source}: {directory}", file=sys.stderr)
        truncated: list[str] = []
        found = _scan_calc_dir(directory, source, quiet=quiet,
                               truncated=truncated)
        all_records.extend(found)
        roots.append({
            "source": source, "path": str(directory), "status": "scanned",
            "records": len(found), "max_depth": scan_depth_for(source),
            "not_scanned_below_depth": sorted(set(truncated))[:20],
        })

    _scan(calc_dir, "calc")
    _scan(archive_dir, "archive")
    _scan(remote_archive_dir, "remote_archive")

    # Build search text for each record
    for rec in all_records:
        rec["_search_text"] = _build_search_text(rec)

    built_at = datetime.now(timezone.utc).isoformat()
    scanned = [r for r in roots if r["status"] == "scanned"]
    return {
        "version": 1,
        "built_at": built_at,
        "record_count": len(all_records),
        "records": all_records,
        "provenance": {
            "built_at": built_at,
            "roots": roots,
            "roots_scanned": len(scanned),
            "any_root_scanned": bool(scanned),
            "truncated": sorted({
                t for r in roots for t in r.get("not_scanned_below_depth", [])
            })[:20],
        },
    }


def _build_search_text(rec: dict[str, Any]) -> str:
    """Build a flat text string for keyword matching."""
    parts = [
        rec.get("calc_id", ""),
        rec.get("name", ""),
        rec.get("smiles", ""),
        rec.get("functional", ""),
        rec.get("basis_set", ""),
        rec.get("aux_basis", ""),
        rec.get("ri_method", ""),
        rec.get("dispersion", ""),
        rec.get("metal_basis", ""),
        rec.get("relativity", ""),
        rec.get("solvation", ""),
        rec.get("solvent", ""),
        rec.get("xtb_method", ""),
        rec.get("esd_modus", ""),
        rec.get("workflow_method", ""),
        rec.get("source", ""),
    ]
    # ORCA keywords
    kw = rec.get("orca_keywords", [])
    if isinstance(kw, list):
        parts.extend(kw)
    # Modules
    mods = rec.get("modules", [])
    if isinstance(mods, list):
        parts.extend(mods)
    # Has-data sections
    hd = rec.get("has_data", [])
    if isinstance(hd, list):
        parts.extend(hd)
    # Out files
    of = rec.get("out_files", [])
    if isinstance(of, list):
        parts.extend(of)

    return " ".join(str(p) for p in parts if p).lower()


# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------

def get_default_calc_index_path() -> Path:
    return Path.home() / ".delfin" / "calc_index.json"
