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
        if _is_yes(parsed.get("hyperpol_xtb_module")):
            modules.append("hyperpol_xtb")
        if _is_yes(parsed.get("TADF_xTB_module")):
            modules.append("TADF_xTB")
        if _is_yes(parsed.get("NMR_module")):
            modules.append("NMR")
        if _is_yes(parsed.get("OCCUPIER")):
            modules.append("OCCUPIER")
        if _is_yes(parsed.get("CO2_coordinator")):
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


def completion_of(d) -> tuple:
    """(completed, exit_code) for a calculation folder, from what it left.

    An exit-code file means the run ENDED (the code says how); a run log
    without one means running or crashed; neither means unknown. One
    reading, shared by the index and the energy tools, so the two never
    call the same folder two different things.
    """
    from pathlib import Path as _P
    d = _P(d)
    exit_codes = list(d.glob(".exit_code_*"))
    if exit_codes:
        try:
            code = exit_codes[0].name.split("_")[-1]
            return True, (int(code) if code.isdigit() else code)
        except Exception:
            return True, "unknown"
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
    return _outcome_of(completed, exit_code)


def method_of_folder(d) -> tuple:
    """(functional, basis) from what the folder holds besides the output.

    An ORCA output need not state its method -- a DELFIN run's does not
    -- so a parser that reads only the .out returns "", and every tool
    downstream then files the run under method None. The folder knows:
    DELFIN_Data.json first, then CONTROL.txt, then the .inp header, the
    order the index itself uses.
    """
    import json as _json
    from pathlib import Path as _P
    d = _P(d)
    functional = basis = None
    try:
        if (d / "DELFIN_Data.json").is_file():
            rec = _extract_from_delfin_data(_json.loads((d / "DELFIN_Data.json").read_text(encoding="utf-8")))
            functional = rec.get("functional") or None
            basis = rec.get("basis_set") or None
    except Exception:
        pass
    try:
        if (not functional or not basis) and (d / "CONTROL.txt").is_file():
            rec = _extract_from_control_txt((d / "CONTROL.txt").read_text(encoding="utf-8"))
            basis = basis or rec.get("basis_set") or None
            functional = functional or rec.get("functional") or None
    except Exception:
        pass
    try:
        if not functional or not basis:
            for inp in sorted(d.glob("*.inp")):
                rec = _extract_from_inp_header(inp.read_text(encoding="utf-8", errors="replace"))
                functional = functional or rec.get("inp_functional") or None
                basis = basis or rec.get("inp_basis") or None
                if functional and basis:
                    break
    except Exception:
        pass
    return functional, basis


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



def _extract_from_control_txt(text: str) -> dict[str, Any]:
    """Extract searchable fields from a CONTROL.txt file."""
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
    rec["basis_set"] = kv.get("basis_set", "")

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
    if _is_yes(kv.get("hyperpol_xtb_module")):
        modules.append("hyperpol_xtb")
    if _is_yes(kv.get("TADF_xTB_module")):
        modules.append("TADF_xTB")
    if _is_yes(kv.get("NMR_module")):
        modules.append("NMR")
    if _is_yes(kv.get("OCCUPIER")):
        modules.append("OCCUPIER")
    if _is_yes(kv.get("CO2_coordinator")):
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
}

_KNOWN_BASIS = {
    "def2-svp", "def2-tzvp", "def2-tzvpp", "def2-qzvp", "def2-qzvpp",
    "ma-def2-svp", "ma-def2-tzvp", "ma-def2-tzvpp",
    "cc-pvdz", "cc-pvtz", "cc-pvqz", "aug-cc-pvdz", "aug-cc-pvtz",
    "6-31g*", "6-31g**", "6-311g**", "6-31+g*", "6-311+g(2d,p)",
    "sto-3g", "def2-svp(d)",
    "sarc-dkh-tzvp", "sarc-zora-tzvp", "x2c-tzvp", "x2c-svp",
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
            break
    return rec


def _is_yes(val: Any) -> bool:
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.strip().lower() in ("yes", "true", "1")
    return False


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
    archive structures like ``archive/Fritz/calc_name/``.

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

        if has_control and "index_source" not in rec:
            try:
                text = (d / "CONTROL.txt").read_text(encoding="utf-8")
                rec.update(_extract_from_control_txt(text))
                rec["index_source"] = "CONTROL.txt"
            except Exception:
                pass
        elif has_control and "functional" not in rec:
            # Supplement from CONTROL if DELFIN_Data didn't have metadata
            try:
                text = (d / "CONTROL.txt").read_text(encoding="utf-8")
                ctrl_rec = _extract_from_control_txt(text)
                for k, v in ctrl_rec.items():
                    if k not in rec or not rec[k]:
                        rec[k] = v
            except Exception:
                pass

        if has_inp and not rec.get("orca_keywords"):
            # Read first .inp for ORCA keywords
            for inp_file in sorted(d.glob("*.inp"))[:1]:
                try:
                    text = inp_file.read_text(encoding="utf-8")[:500]
                    rec.update(_extract_from_inp_header(text))
                except Exception:
                    pass

        # Fill functional/basis from .inp if still missing
        if not rec.get("functional") and rec.get("inp_functional"):
            rec["functional"] = rec["inp_functional"]
        if not rec.get("basis_set") and rec.get("inp_basis"):
            rec["basis_set"] = rec["inp_basis"]

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
