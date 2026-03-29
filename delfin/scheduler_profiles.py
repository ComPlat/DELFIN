"""Scheduler stage profiles, stage inference, and molecule feature helpers."""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

_DEFAULT_PROFILE_PATH = Path(__file__).with_name("data") / "scheduler_stage_profiles.json"
_PROFILE_CACHE: Dict[str, Dict[str, Any]] | None = None
_PROFILE_MTIME_NS: int | None = None
_PROFILE_PATH_STR: str | None = None
_FEATURE_CACHE: Dict[str, Dict[str, Any]] = {}

_XYZ_COORD_RE = re.compile(r"^\s*([A-Z][a-z]?)\s+[-+]?\d")
_SMILES_BRACKET_RE = re.compile(r"^\[([A-Z][a-z]?)")
_AROMATIC_SMILES_MAP = {
    "b": "B",
    "c": "C",
    "n": "N",
    "o": "O",
    "p": "P",
    "s": "S",
}
_TWO_CHAR_ELEMENTS = {
    "Ag", "Al", "As", "Au", "Br", "Ca", "Cl", "Co", "Cr", "Cu", "Fe", "Hg",
    "Ir", "Li", "Mg", "Mn", "Mo", "Na", "Ni", "Os", "Pd", "Pt", "Rh", "Ru",
    "Sc", "Se", "Si", "Sn", "Ti", "Zn", "Zr",
}
_TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}


def get_stage_profile_path() -> Path:
    override = str(os.environ.get("DELFIN_SCHEDULER_STAGE_PROFILE_FILE", "")).strip()
    if override:
        return Path(override).expanduser()
    return _DEFAULT_PROFILE_PATH


def get_stage_profiles(force_reload: bool = False) -> Dict[str, Dict[str, Any]]:
    global _PROFILE_CACHE, _PROFILE_MTIME_NS, _PROFILE_PATH_STR

    profile_path = get_stage_profile_path()

    try:
        stat = profile_path.stat()
    except OSError:
        _PROFILE_CACHE = {}
        _PROFILE_MTIME_NS = None
        _PROFILE_PATH_STR = str(profile_path)
        return {}

    if (
        not force_reload
        and _PROFILE_CACHE is not None
        and _PROFILE_MTIME_NS == stat.st_mtime_ns
        and _PROFILE_PATH_STR == str(profile_path)
    ):
        return _PROFILE_CACHE

    raw = json.loads(profile_path.read_text(encoding="utf-8"))
    profiles: Dict[str, Dict[str, Any]] = {}
    for stage, payload in raw.items():
        if str(stage).startswith("_") or not isinstance(payload, dict):
            continue
        profiles[str(stage)] = dict(payload)

    _PROFILE_CACHE = profiles
    _PROFILE_MTIME_NS = stat.st_mtime_ns
    _PROFILE_PATH_STR = str(profile_path)
    return profiles


def infer_stage_key(job_id: str = "", description: str = "") -> str:
    text = f"{job_id} {description}".lower()

    if "fluorescence" in text:
        return "fluorescence_S1_S0"
    if "phosphorescence" in text:
        return "phosphorescence_T1_S0"
    if "tddft check" in text:
        if "s1" in text:
            return "tddft_check_S1"
        if "s2" in text:
            return "tddft_check_S2"
        if "t1" in text:
            return "tddft_check_T1"
        if "t2" in text:
            return "tddft_check_T2"
        return "tddft_check"
    if "esd s0" in text or "esd_s0" in text:
        return "esd_S0"
    if "esd s1" in text or "esd_s1" in text:
        return "esd_S1"
    if "esd s2" in text or "esd_s2" in text:
        return "esd_S2"
    if "esd t1" in text or "esd_t1" in text:
        return "esd_T1"
    if "esd t2" in text or "esd_t2" in text:
        return "esd_T2"
    if "esd t3" in text or "esd_t3" in text:
        return "esd_T3"
    if "esd t4" in text or "esd_t4" in text:
        return "esd_T4"
    if text.startswith("ic ") or " ic " in text:
        return "ic"
    if text.startswith("isc ") or " isc " in text:
        return "isc"
    if "single_point" in text or "sp " in text or text.endswith(" sp"):
        return "single_point"
    if "red_step_3" in text or "red_3" in text:
        return "red_step_3"
    if "red_step_2" in text or "red_2" in text:
        return "red_step_2"
    if "red_step_1" in text or "red_1" in text:
        return "red_step_1"
    if "ox_step_1" in text or "ox_1" in text:
        return "ox_step_1"
    if "initial" in text or "s0 optimization" in text:
        return "initial"
    return "other"


def _extract_elements_from_xyz_lines(lines: Iterable[str]) -> list[str]:
    elements: list[str] = []
    for line in lines:
        match = _XYZ_COORD_RE.match(line)
        if match:
            elements.append(match.group(1))
    return elements


def _parse_xyz_or_coord_text(text: str) -> Optional[Dict[str, Any]]:
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    if not lines:
        return None

    if len(lines) >= 3:
        try:
            atom_count = int(lines[0].strip())
        except ValueError:
            atom_count = None
        if atom_count and atom_count > 0 and len(lines) >= atom_count + 2:
            elements = _extract_elements_from_xyz_lines(lines[2: 2 + atom_count])
            if len(elements) == atom_count:
                return {"atom_count": atom_count, "elements": elements, "source": "xyz"}

    elements = _extract_elements_from_xyz_lines(lines)
    if len(elements) >= 2:
        return {"atom_count": len(elements), "elements": elements, "source": "coords"}

    return None


def _parse_smiles_text(text: str) -> Optional[Dict[str, Any]]:
    line = text.strip()
    if not line or "\n" in line:
        return None

    try:
        from rdkit import Chem  # type: ignore

        mol = Chem.MolFromSmiles(line)
        if mol is not None:
            elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
            return {"atom_count": len(elements), "elements": elements, "source": "smiles"}
    except Exception:
        pass

    elements: list[str] = []
    i = 0
    while i < len(line):
        ch = line[i]
        if ch == "[":
            j = line.find("]", i + 1)
            if j == -1:
                break
            content = line[i + 1: j]
            match = _SMILES_BRACKET_RE.match(f"[{content}")
            if match:
                elements.append(match.group(1))
            i = j + 1
            continue
        if ch in _AROMATIC_SMILES_MAP:
            elements.append(_AROMATIC_SMILES_MAP[ch])
            i += 1
            continue
        if ch.isupper():
            if i + 1 < len(line) and line[i:i + 2] in _TWO_CHAR_ELEMENTS:
                elements.append(line[i:i + 2])
                i += 2
                continue
            if ch != "H":
                elements.append(ch)
        i += 1

    if not elements:
        return None
    return {"atom_count": len(elements), "elements": elements, "source": "smiles_fallback"}


def parse_molecule_text(text: str) -> Optional[Dict[str, Any]]:
    parsed = _parse_xyz_or_coord_text(text)
    if parsed is None:
        parsed = _parse_smiles_text(text)
    if parsed is None:
        return None

    elements = list(parsed.get("elements", []))
    parsed["has_transition_metal"] = any(element in _TRANSITION_METALS for element in elements)
    parsed["transition_metals"] = sorted({element for element in elements if element in _TRANSITION_METALS})
    return parsed


def _candidate_input_files(base_dir: Path) -> list[Path]:
    roots = [base_dir]
    if base_dir.parent != base_dir:
        roots.append(base_dir.parent)

    names = [
        "input.txt",
        "start.txt",
        "initial.xyz",
        "S0.xyz",
        "S1.xyz",
        "T1.xyz",
        "red_step_1.xyz",
        "red_step_2.xyz",
        "red_step_3.xyz",
        "ox_step_1.xyz",
    ]
    candidates: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        for name in names:
            path = (root / name).resolve()
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            candidates.append(path)
    return candidates


def get_molecule_features(base_dir: Optional[Path]) -> Dict[str, Any]:
    root = Path(base_dir) if base_dir is not None else Path.cwd()
    cache_key = str(root.resolve())
    cached = _FEATURE_CACHE.get(cache_key)
    if cached is not None:
        return dict(cached)

    for candidate in _candidate_input_files(root):
        if not candidate.exists():
            continue
        parsed = parse_molecule_text(candidate.read_text(encoding="utf-8", errors="ignore"))
        if parsed is None:
            continue
        parsed["path"] = str(candidate)
        _FEATURE_CACHE[cache_key] = dict(parsed)
        return parsed

    fallback = {
        "atom_count": None,
        "elements": [],
        "has_transition_metal": False,
        "transition_metals": [],
        "source": "missing",
        "path": None,
    }
    _FEATURE_CACHE[cache_key] = dict(fallback)
    return fallback


def _apply_profile_modifiers(
    profile: Dict[str, Any],
    features: Dict[str, Any],
    *,
    duration_s: float,
    cores: float,
) -> tuple[float, float]:
    modifiers = profile.get("modifiers")
    if not isinstance(modifiers, dict):
        return duration_s, cores

    for key, payload in modifiers.items():
        if not isinstance(payload, dict):
            continue
        active = False
        if key == "has_transition_metal":
            active = bool(features.get("has_transition_metal"))
        elif key.startswith("transition_metal:"):
            metal = key.split(":", 1)[1]
            active = metal in set(features.get("transition_metals", []))
        if not active:
            continue
        duration_s *= float(payload.get("duration_mult", 1.0) or 1.0)
        cores *= float(payload.get("cores_mult", 1.0) or 1.0)

    return duration_s, cores


def predict_stage_profile(
    stage: str,
    *,
    atom_count: Optional[int] = None,
    features: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    profiles = get_stage_profiles()
    profile = dict(profiles.get(stage) or profiles.get("other") or {})
    if not profile:
        return {
            "duration_s": None,
            "avg_cores": None,
            "observed_avg_cores": None,
            "recommended_cores": None,
            "profile": {},
            "stage": stage,
        }

    duration_s = float(profile.get("duration_s", 0.0) or 0.0)
    observed_avg_cores = float(
        profile.get("observed_avg_cores", profile.get("avg_cores", 0.0)) or 0.0
    )
    recommended_cores = float(profile.get("recommended_cores", 0.0) or 0.0)
    reference_atoms = float(profile.get("reference_atoms", 0.0) or 0.0)
    atom_scaling_exp = float(profile.get("atom_scaling_exp", 1.0) or 1.0)
    core_scaling_exp = float(profile.get("core_scaling_exp", min(1.0, max(0.0, atom_scaling_exp * 0.35))) or 0.0)

    if atom_count is not None and atom_count > 0 and reference_atoms > 0:
        atom_ratio = max(0.25, float(atom_count) / reference_atoms)
        duration_s *= math.pow(atom_ratio, atom_scaling_exp)
        observed_avg_cores *= math.pow(atom_ratio, core_scaling_exp)
        if recommended_cores > 0:
            recommended_cores *= math.pow(atom_ratio, core_scaling_exp)

    if features:
        duration_s, observed_avg_cores = _apply_profile_modifiers(
            profile,
            features,
            duration_s=duration_s,
            cores=observed_avg_cores,
        )
        if recommended_cores > 0:
            _duration_unused, recommended_cores = _apply_profile_modifiers(
                profile,
                features,
                duration_s=duration_s,
                cores=recommended_cores,
            )

    return {
        "stage": stage,
        "profile": profile,
        "duration_s": duration_s if duration_s > 0 else None,
        "avg_cores": observed_avg_cores if observed_avg_cores > 0 else None,
        "observed_avg_cores": observed_avg_cores if observed_avg_cores > 0 else None,
        "recommended_cores": recommended_cores if recommended_cores > 0 else None,
    }
