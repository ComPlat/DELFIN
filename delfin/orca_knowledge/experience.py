"""ORCA experience extraction: one record per run of what was computed, how
it ended and what trouble occurred.

Input: an ORCA input file (``.inp`` or ``input.txt``) and its output file
(``.out``). Output: a record with the settings used (simple-input keywords,
block names and options, charge, multiplicity), the system class (metals,
atom counts, charge, multiplicity), the job type(s) and method family, how
the run ended (normal termination, SCF iterations and convergence,
optimisation cycles and convergence, imaginary frequencies, wall time), and
the error/warning signature lines ORCA printed.

Privacy: records never contain coordinates, element-pair lists beyond
counts, molecule names or user names taken from paths, or file contents
beyond the fields below. A run is identified by a hash of its directory
path (:func:`dir_hash`), never by the path itself.

Used by DELFIN to later learn which settings work for which kind of system
and which fix healed which failure (:func:`link`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

__all__ = [
    "extract", "link", "dir_hash", "system_class", "main",
    "MAX_FILE_BYTES",
]

MAX_FILE_BYTES = 200 * 1024 * 1024  # skip files over 200 MB

# Metals: all elements DELFIN users would call metals (alkali/alkaline-earth,
# transition, lanthanides, post-transition). Reported via system_class.
_METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi", "Po",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
}
# Transition metals and heavy metals (Z >= 21), reported separately.
_TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
}
_HEAVY_ELEMENTS_Z21 = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
}

_LIGHT_ELEMENTS = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"}

# Every element symbol (1-118), for validating geometry tokens: only real
# element symbols may be counted in element_counts -- never input tokens
# like "$new_job", "!", "end" or "nroots" that an unterminated geometry
# section would otherwise swallow.
_ELEMENT_SYMBOLS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}

# Method families, from the simple-input keywords.
_METHOD_FAMILIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("xtb", ("XTB", "XTB2", "GFN1", "GFN2", "XTB1")),
    ("semiempirical", ("AM1", "PM3", "MNDO", "PM6", "OM2", "OM3", "PPP", "INDO")),
    ("casscf", ("CASSCF", "NEVPT2", "CASPT2", "DMRG", "ICE-CI", "ICECI")),
    ("coupled-cluster", ("CCSD", "CCSD(T)", "CCSDT", "DLPNO-CCSD", "DLPNO-CCSD(T)",
                         "DLPNO-MP2", "CCSD(T)-F12")),
    ("double-hybrid", ("B2PLYP", "DSD-BLYP", "PWPB95", "DOD-PBE", "PT2B95")),
    ("mp2", ("MP2", "RI-MP2", "DLPNO-MP2", "SCS-MP2", "MP2-GRAD")),
    ("compound", ("COMPOUND",)),
    ("goat", ("GOAT",)),
    ("tddft", ("TDDFT", "TDA", "TDDFRT")),
    ("esd", ("ESD",)),
    ("hf", ("HF",)),
    ("dft", ()),
)
# First keyword-only families; "dft" is the fallback for anything with a
# DFT functional keyword and no higher family.

_FUNCTIONALS = re.compile(
    r"^(PBE|BP86|B3LYP|PBE0|TPSS|TPSSh|M06|M06L|M06-2X|BHandHLYP|PW6B95|"
    r"B97-3|B3LYP5|CAM-B3LYP|ωB97X|WB97X|WB97X-D|WB97X-V|BMK|B1LYP|"
    r"B2PLYP|DSD-BLYP|PWPB95|DOD-PBE|PT2B95|revDSD|PBE0-2|SOGGA11-X|"
    r"TPSS0|B97M-V|PBEh-3c|HF-3c|B97-3c|r2SCAN|SCAN|R2SCAN|OLYP|OPBE|"
    r"KT1|KT2|BHandH|BHandHLC|PW91|RPBE|REVPBE|PBEx|OPTX|OLYP|KYT|"
    r"B1B95|MPW1K|MPWPW91|N12|GAM|BMK|MN15|MN15-L|MN12-L|MN12-SX|"
    r"LC-ωPBE|LRC-ωPBE|CAM-QTP-01|QTP)$", re.IGNORECASE)

_JOB_KEYWORDS: dict[str, tuple[str, ...]] = {
    "opt": ("OPT", "OPTTS", "BESSER", "LOOSEOPT", "TIGHTOPT", "NROPT", "NRS",
            "NR", "Z-OPT"),
    "ts": ("OPTTS", "TS", "TSRC", "NEB", "NEB-TS", "QSHING"),
    "freq": ("FREQ", "NUMFREQ", "ANFREQ", "VPT2", "AHAS", "RCA"),
    "sp": (),  # default: no other job keyword
    "goat": ("GOAT",),
    "esd": ("ESD",),
    "tddft": ("TDDFT", "TDA", "TD-CIS", "TDDFRT", "CIS", "ROCIS", "EOM-CCSD"),
    "compound": ("COMPOUND",),
    "scan": ("SCAN",),
    "irc": ("IRC",),
    "md": ("MD", "AIMD", "MDRR"),
    "npcrelax": ("NPCR",),
    "xtb-only": ("XTB", "XTB2", "GFN1", "GFN2"),
}
# A run is a "sp" run when none of opt/ts/freq/goat/esd/tddft/... keywords
# are present.

_ERROR_SIGNATURES: tuple[tuple[str, str], ...] = (
    # (regex over the output text, normalised signature)
    (r"SCF NOT CONVERGED AFTER\s+\d+\s+CYCLES", "SCF NOT CONVERGED"),
    (r"unfortunately, the SCF has not converged", "SCF NOT CONVERGED"),
    (r"ORCA finished by error termination in\s+(\w+)",
     "ORCA finished by error termination in {1}"),
    (r"ORCA finished by error termination", "ORCA finished by error termination"),
    (r"SCF convergence could not be reached", "SCF NOT CONVERGED"),
    (r"GEOMETRY OPTIMIZATION DID NOT (?:CONVERGE|FINISH)", "OPTIMIZATION DID NOT CONVERGE"),
    (r"The optimization did not converge", "OPTIMIZATION DID NOT CONVERGE"),
    (r"OPTIMIZATION HAS CONVERGED", "OPTIMIZATION CONVERGED"),
    (r"Sorry, no answer from the server", "SORRY, NO ANSWER FROM THE SERVER"),
    (r"Convergence failure - run restarted", "CONVERGENCE FAILURE - RUN RESTARTED"),
    (r"SCF failed", "SCF FAILED"),
    (r"MP2 NOT CONVERGED", "MP2 NOT CONVERGED"),
    (r"Gradient convergence criteria cannot be met", "GRADIENT CRITERIA CANNOT BE MET"),
    (r"Geometry optimization did not converge", "OPTIMIZATION DID NOT CONVERGE"),
    (r"ERROR: the \w+ basis was not found", "BASIS SET NOT FOUND"),
    (r"Error in the input file", "INPUT ERROR"),
    (r"ORCA has detected (.*)failed", "ORCA DETECTED FAILURE: {1}"),
    (r"membarray", "ORCA INTERNAL ERROR (MEMBARRAY)"),
    (r"Segmentation fault", "SEGMENTATION FAULT"),
    (r"Fatal Python error", "FATAL PYTHON ERROR"),
    (r"NEB.*did not converge", "NEB DID NOT CONVERGE"),
    (r"TDDFT.*not converged", "TDDFT NOT CONVERGED"),
    (r"converge, but only slowly", "SCF CONVERGES SLOWLY"),
    (r"Slow convergence", "SCF CONVERGES SLOWLY"),
    (r"the optimum geometry could not be found", "GEOMETRY NOT FOUND"),
    (r"Stopped - The (?:(?:\S+ ?){1,4})force constant matrix", "HESSIAN STOPPED"),
    (r"Number of imaginary modes\s*\|?\s*[:=]?\s*(\d+)", "IMAGINARY MODES {1}"),
    (r"Possible imaginary frequencies", "POSSIBLE IMAGINARY FREQUENCIES"),
    (r"mpirun noticed that process rank", "MPI PROCESS EXITED ABNORMALLY"),
    (r"Primary job\s+terminated normally", "MPI JOB TERMINATED ABNORMALLY"),
    (r"Too many bad steps", "TOO MANY BAD STEPS"),
    (r"Failed to converge in the (?:(?:\S+ ?){1,3})microiterations", "MICROITERATION FAILURE"),
    (r"Keyword cannot be used with keyword", "INCOMPATIBLE KEYWORDS"),
    (r"The EXCHANGESCALING", "EXCHANGESCALING ERROR"),
    (r"was not found in the library", "KEYWORD NOT IN LIBRARY"),
    (r"Not implemented for this (?:(?:\S+ ?){1,3})type", "NOT IMPLEMENTED"),
    (r"Time limit reached", "WALLTIME LIMIT REACHED"),
    (r"JOB DONE", "JOB DONE"),
    (r"Attempt to allocate memory", "MEMORY ALLOCATION FAILURE"),
    (r"Out of memory", "OUT OF MEMORY"),
    (r"orca_\w+_mpi was terminated by signal", "ORCA KILLED BY SIGNAL"),
    (r"Please submit a bug report", "ORCA BUG REPORT REQUESTED"),
)

# Signatures that are informational, not trouble: not reported as errors.
_INFO_SIGNATURES = {
    "JOB DONE", "OPTIMIZATION CONVERGED", "DEGREES OF FREEDOM",
    "POSSIBLE IMAGINARY FREQUENCIES", "IMAGINARY MODES 0",
}

_WARNING_SIGNATURES: tuple[tuple[str, str], ...] = (
    (r"Warning.*: (.*)", "WARNING: {1}"),
)


def dir_hash(path: str | os.PathLike) -> str:
    """A short hash of the directory path -- the path itself never leaves."""
    h = hashlib.sha256(os.fspath(path).encode("utf-8", "surrogateescape"))
    return h.hexdigest()[:16]


def _read_text(path: Path) -> str | None:
    try:
        if path.stat().st_size > MAX_FILE_BYTES:
            return None
        with open(path, "rb") as fh:
            data = fh.read()
        return data.decode("utf-8", "replace")
    except OSError:
        return None


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------

_XYZ_HEAD = re.compile(r"^\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE)
_BLOCK = re.compile(r"^\s*%\s*(\w+)\b(.*?)$", re.IGNORECASE)


def parse_input(inp_text: str) -> dict:
    """Keywords, blocks, charge/multiplicity and element counts from an input.

    Public alias of :func:`_parse_input_scanner`; parses every job of a
    multi-job input and never returns coordinates.
    """
    return _parse_input_scanner(inp_text)


def _parse_input_scanner(inp_text: str) -> dict:
    """Scan an ORCA input, parsing EVERY job of a multi-job input.

    A state machine over the lines, split at ``$new_job``:

    - a block (``%name ...``) collects option lines until its ``end``;
      ``$new_job`` also closes an open block (defensive: real inputs
      sometimes omit the block's ``end`` before a new job);
    - the geometry section (``* xyz ... *``) is NEVER read into blocks or
      keywords -- only element symbols from its atom lines are counted, and
      only tokens that are real element symbols;
    - ``$``-directive lines, ``!`` keyword lines and ``%`` block lines are
      never counted as elements.

    No coordinates are returned -- only counts per element.
    """
    keywords: list[str] = []
    blocks: dict[str, list[str]] = defaultdict(list)
    charge = None
    multiplicity = None
    elements: Counter = Counter()

    def finish_job() -> None:
        nonlocal in_xyz, current
        in_xyz = False
        current = None

    in_xyz = False
    current = None
    for raw in inp_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        # $new_job starts a new job: close whatever is open.
        if line.startswith("$"):
            finish_job()
            continue
        if in_xyz:
            if line.startswith("*"):
                in_xyz = False
            else:
                parts = line.split()
                token = parts[0]
                if token in _ELEMENT_SYMBOLS:
                    elements[token] += 1
            continue
        m = _XYZ_HEAD.match(line)
        if m:
            charge = int(m.group(1))
            multiplicity = int(m.group(2))
            in_xyz = True
            continue
        if line.lower() == "end":
            current = None
            continue
        if line.startswith("!"):
            for kw in line[1:].split():
                keywords.append(kw)
            continue
        m = _BLOCK.match(line)
        if m:
            name = m.group(1).upper()
            opts = m.group(2).strip()
            if opts.lower() == "end":
                current = None
                continue
            current = name
            if opts:
                blocks[name].append(opts)
            continue
        if current and not line.startswith("!"):
            blocks[current].append(line)

    _scrub_blocks(blocks)
    return {
        "simple_input_keywords": sorted({k.upper() for k in keywords}),
        "blocks": blocks,
        "charge": charge,
        "multiplicity": multiplicity,
        "element_counts": dict(elements),
    }


def _scrub_blocks(blocks: dict[str, list[str]]) -> None:
    """Privacy scrub of collected block options, in place.

    - drops any option line that looks like a coordinate (element followed
      by numbers) or a geometry directive -- a mis-terminated block would
      otherwise leak the geometry into the record;
    - replaces path tokens with ``<path>`` so no directory, user name or
      molecule name from a path survives in the record.
    """
    coord_like = re.compile(
        r"^[A-Za-z]{1,3}\s+-?\d+\.?\d*\s+-?\d+\.?\d*(\s+-?\d+\.?\d*)?$")
    path_like = re.compile(r"\S*/\S*")
    for name in list(blocks):
        kept = []
        for o in blocks[name]:
            if coord_like.match(o.strip()) or o.strip().startswith("*"):
                continue
            kept.append(path_like.sub("<path>", o).strip())
        if kept:
            blocks[name] = kept
        else:
            del blocks[name]


# --------------------------------------------------------------------------
# Output parsing
# --------------------------------------------------------------------------

_SCF_DONE = re.compile(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES")
_SCF_FAIL = re.compile(r"SCF NOT CONVERGED AFTER\s+(\d+)\s+CYCLES")
_OPT_CYCLE = re.compile(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)")
_OPT_DONE = re.compile(r"THE OPTIMIZATION HAS CONVERGED")
_OPT_FAIL = re.compile(r"THE OPTIMIZATION DID NOT CONVERGE")
_NATOMS = re.compile(r"Number of atoms\s+\.+\s+(\d+)")
_CHARGE = re.compile(r"Total Charge\s+Charge\s+\.+\s+(-?\d+)")
_MULT = re.compile(r"Multiplicity\s+Mult\s+\.+\s+(\d+)")
_VERSION = re.compile(r"Program Version\s+([\d.]+)")
_NIMAG = re.compile(
    r"(?:imaginary modes?|imaginary frequencies?)\s*[:=|]?\s*(\d+)",
    re.IGNORECASE)
_FREQ_PRE = re.compile(r"(?:are|is)\s+(\d+)\s+imaginary frequencies?",
                       re.IGNORECASE)
_FREQ_WARNING = re.compile(r"\*\*\*\s+(\d+)\s+IMAGINARY FREQUENC", re.IGNORECASE)
_WALL = re.compile(r"Total Wall time\s+:?\s*([\d.]+)s.*", re.IGNORECASE)
_WALL_ALT = re.compile(r"Total wall time\s+(\d+):(\d+):(\d+)", re.IGNORECASE)
_NORMAL = re.compile(r"\\*\s*V I R A L   T H E O R E M", re.IGNORECASE)
_NORMAL2 = re.compile(r"ORCA TERMINATED NORMALLY", re.IGNORECASE)


def _parse_out(text: str) -> dict:
    out: dict = {
        "orca_version": None,
        "normal_termination": False,
        "scf_iterations": None,
        "scf_converged": None,
        "optimization_cycles": None,
        "optimization_converged": None,
        "imaginary_frequencies": None,
        "wall_time_seconds": None,
        "error_signatures": [],
        "warning_signatures": [],
        "atom_count": None,
        "charge": None,
        "multiplicity": None,
    }
    m = _VERSION.search(text)
    if m:
        out["orca_version"] = m.group(1)
    m = _NATOMS.search(text)
    if m:
        out["atom_count"] = int(m.group(1))
    m = _CHARGE.search(text)
    if m:
        out["charge"] = int(m.group(1))
    m = _MULT.search(text)
    if m:
        out["multiplicity"] = int(m.group(1))
    for m in _SCF_DONE.finditer(text):
        out["scf_converged"] = True
        out["scf_iterations"] = int(m.group(1))
    for m in _SCF_FAIL.finditer(text):
        out["scf_converged"] = False
        out["scf_iterations"] = int(m.group(1))
    for m in _OPT_CYCLE.finditer(text):
        out["optimization_cycles"] = int(m.group(1))
    if _OPT_DONE.search(text):
        out["optimization_converged"] = True
    elif _OPT_FAIL.search(text):
        out["optimization_converged"] = False
    for pattern in (_FREQ_WARNING, _FREQ_PRE, _NIMAG):
        m = pattern.search(text)
        if m:
            out["imaginary_frequencies"] = int(m.group(1))
            break
    m = _WALL.search(text)
    if m:
        out["wall_time_seconds"] = float(m.group(1))
    else:
        m = _WALL_ALT.search(text)
        if m:
            out["wall_time_seconds"] = (
                int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
            )
    out["normal_termination"] = bool(
        _NORMAL2.search(text) or _NORMAL.search(text)
    )
    sigs: list[str] = []
    warns: list[str] = []
    for pattern, signature in _ERROR_SIGNATURES:
        m = re.search(pattern, text)
        if not m:
            continue
        sig = signature
        for i in range(1, 4):
            sig = sig.replace("{%d}" % i, (m.group(i) or "") if m.lastindex and m.lastindex >= i else "")
        sig = sig.strip()
        if sig in _INFO_SIGNATURES or sig == "IMAGINARY MODES 0":
            continue
        sigs.append(sig)
    for pattern, signature in _WARNING_SIGNATURES:
        for m in re.finditer(pattern, text):
            w = signature.replace("{1}", m.group(1)).strip()[:200]
            # Privacy: ORCA warnings can quote file paths; scrub them.
            warns.append(re.sub(r"\S*/\S*", "<path>", w))
    # dedup, keep order
    out["error_signatures"] = sorted(set(sigs))
    out["warning_signatures"] = sorted(set(warns))
    return out


# --------------------------------------------------------------------------
# Method family, job types, system class
# --------------------------------------------------------------------------

def _method_family(keywords: list[str]) -> str:
    ks = {k.upper() for k in keywords}
    for family, kws in _METHOD_FAMILIES:
        if family == "dft":
            continue
        if any(k in ks for k in kws):
            return family
    if any(_FUNCTIONALS.match(k) for k in ks):
        return "dft"
    if "HF" in ks:
        return "hf"
    return "unknown"


def _job_types(keywords: list[str]) -> list[str]:
    ks = {k.upper() for k in keywords}
    jobs: list[str] = []
    for job, kws in _JOB_KEYWORDS.items():
        if job == "sp":
            continue
        if any(k in ks for k in kws):
            jobs.append(job)
    if not jobs:
        jobs = ["sp"]
    return jobs


def _metals_of(elements: dict[str, int]) -> list[dict]:
    """Metals present, as [{element, count}]. Order: most abundant first."""
    found = [
        {"element": el, "count": cnt}
        for el, cnt in sorted(elements.items(), key=lambda kv: (-kv[1], kv[0]))
        if el in _METALS
    ]
    return found


def system_class(record: dict) -> dict:
    """A privacy-safe summary of the system: metals, charge, mult, size."""
    elements = record.get("element_counts") or {}
    natoms = sum(elements.values())
    metals = _metals_of(elements)
    transition = {m["element"]: m["count"] for m in metals
                  if m["element"] in _TRANSITION_METALS}
    light = sorted(el for el in elements if el in _LIGHT_ELEMENTS)
    heavy = sorted(el for el in elements if el in _HEAVY_ELEMENTS_Z21)
    return {
        "total_atoms": natoms,
        "total_charge": record.get("charge"),
        "multiplicity": record.get("multiplicity"),
        "metals": metals,
        "transition_metals": transition,
        "light_elements": light,
        "heavy_elements": [e for e in heavy if e not in transition],
        "element_counts": dict(sorted(elements.items())),
    }


# --------------------------------------------------------------------------
# extract()
# --------------------------------------------------------------------------

def extract(inp_path: str | os.PathLike, out_path: str | os.PathLike) -> dict:
    """One record for one ORCA run: settings, system, outcome, trouble.

    Neither path is recorded -- only a hash of the containing directory.
    """
    inp_path = Path(inp_path)
    out_path = Path(out_path)
    inp_text = _read_text(inp_path)
    out_text = _read_text(out_path)
    if inp_text is None or out_text is None:
        return {}
    record: dict = {
        "dir_hash": dir_hash(inp_path.parent),
        # A hash, never the name: file names can carry molecule names.
        "run_hash": hashlib.sha256(str(inp_path.resolve()).encode()).hexdigest()[:16],
    }
    parsed = parse_input(inp_text)
    record.update({
        "orca_version": None,
        "simple_input_keywords": parsed["simple_input_keywords"],
        "blocks": parsed["blocks"],
        "charge": parsed["charge"],
        "multiplicity": parsed["multiplicity"],
        "element_counts": parsed["element_counts"],
        "job_types": _job_types(parsed["simple_input_keywords"]),
        "method_family": _method_family(parsed["simple_input_keywords"]),
    })
    from_out = _parse_out(out_text) if out_text else {}
    for key, value in from_out.items():
        if value is not None:
            record[key] = value
    record["normal_termination"] = from_out.get("normal_termination", False)
    record["system_class"] = system_class(record)
    return record


# --------------------------------------------------------------------------
# link(): pair a failed run with the later run in the same directory that
# succeeded, recording what changed (the "healing")
# --------------------------------------------------------------------------

def _diff_keys(before: dict, after: dict, key: str) -> list[str]:
    """Added / removed entries of a keywords list or block options, as
    human-readable change strings."""
    b = set(before or [])
    a = set(after or [])
    return sorted([f"+{x}" for x in a - b] + [f"-{x}" for x in b - a])


def link(records: list[dict]) -> list[dict]:
    """Healings: for every failed record, the next later successful record in
    the same directory (same ``dir_hash``) with what changed in keywords and
    block options."""
    by_dir: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_dir[r["dir_hash"]].append(r)
    heals: list[dict] = []
    for dhash, rs in by_dir.items():
        # Keep input order; pair each failed record with the next success.
        for i, failed in enumerate(rs):
            if failed.get("normal_termination") or not failed.get("error_signatures"):
                continue
            for later in rs[i + 1:]:
                if later.get("normal_termination") and not later.get("error_signatures"):
                    heals.append({
                        "dir_hash": dhash,
                        "failed_run": failed.get("run_hash"),
                        "healed_run": later.get("run_hash"),
                        "error_signatures": failed.get("error_signatures", []),
                        "keyword_changes": _diff_keys(
                            failed.get("simple_input_keywords"),
                            later.get("simple_input_keywords"), "keywords"),
                        "block_changes": {
                            name: _diff_keys(
                                (failed.get("blocks") or {}).get(name, []),
                                (later.get("blocks") or {}).get(name, []), name)
                            for name in sorted(set(failed.get("blocks") or {})
                                               | set(later.get("blocks") or {}))
                        },
                        "changed_blocks": sorted(
                            name for name in
                            set(failed.get("blocks") or {}) | set(later.get("blocks") or {})
                            if (failed.get("blocks") or {}).get(name)
                            != (later.get("blocks") or {}).get(name)),
                    })
                    break
    return heals


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

_INP_NAMES = re.compile(r".*\.(inp|txt|in)$", re.IGNORECASE)
_OUT_NAMES = re.compile(r".*\.(out)$", re.IGNORECASE)


def _candidates(dirpath: str, filenames: list[str]) -> list[tuple[Path, Path]]:
    """Match .inp (or input.txt) files with .out files of the same stem.
    DELFIN inputs are input.txt; ORCA outputs are <stem>.out."""
    outs = {Path(f).stem: f for f in filenames if _OUT_NAMES.match(f)}
    pairs = []
    for f in filenames:
        if not _INP_NAMES.match(f):
            continue
        stem = Path(f).stem
        if stem in outs:
            pairs.append((Path(dirpath) / f, Path(dirpath) / outs[stem]))
    return pairs


def _walk_records(root: str, limit: int | None = None) -> list[dict]:
    """Walk ROOT read-only, never following symlinks out of ROOT, skipping
    files over 200 MB."""
    root_real = os.path.realpath(root)
    records: list[dict] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        # Never follow a symlinked directory out of ROOT.
        pruned = []
        for name in dirnames:
            full = os.path.join(dirpath, name)
            if os.path.islink(full) and not os.path.realpath(full).startswith(root_real + os.sep):
                continue
            pruned.append(name)
        dirnames[:] = pruned
        for inp, out in _candidates(dirpath, filenames):
            if limit is not None and len(records) >= limit:
                return records
            rec = extract(inp, out)
            if rec:
                records.append(rec)
    return records


def _summary(records: list[dict], heals: list[dict]) -> dict:
    err_counts: Counter = Counter()
    for r in records:
        for s in r.get("error_signatures", []):
            err_counts[s] += 1
    heal_counts: Counter = Counter()
    for h in heals:
        key = ",".join(h["error_signatures"]) + " -> " + ",".join(
            h["keyword_changes"] + [f"%{b}" for b in h["changed_blocks"]])
        heal_counts[key] += 1
    return {
        "record_count": len(records),
        "normal_terminations": sum(1 for r in records if r.get("normal_termination")),
        "error_signature_counts": dict(err_counts.most_common()),
        "healing_count": len(heals),
        "top_healings": dict(heal_counts.most_common(10)),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="delfin.orca_knowledge.experience",
        description="Extract ORCA run experience records (read-only walk).",
    )
    ap.add_argument("root", help="root directory to walk (read-only)")
    ap.add_argument("--out", required=True, help="output .jsonl file")
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after N records")
    args = ap.parse_args(argv)
    records = _walk_records(args.root, args.limit)
    heals = link(records)
    summary = _summary(records, heals)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
        fh.write(json.dumps({"summary": summary}, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
