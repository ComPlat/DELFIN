"""Ensemble utilities for CREST / CENSO conformer ensembles.

Complements the existing ``crest_conformers`` / ``censo_sort`` adapters
(delfin/tools/adapters/crest.py, censo.py) with the reading side:

- ``parse_ensemble`` reads a multi-frame XYZ with the energy in the comment
  line, as CREST writes ``crest_conformers.xyz`` (energies in Eh on the
  xtb/level-of-the-search surface).
- ``relative_kcal`` derives relative energies in kcal/mol from the lowest
  frame (627.5094740631 kcal/mol per Eh, the conversion factor of the
  hartree atomic unit of energy; and the CENSO convention for its table).
- ``boltzmann_weights`` weights an ensemble at a stated temperature. Every
  weight uses ONE level of theory; a frame tagged with a different method
  is refused with ValueError instead of being silently ranked (the
  two-levels trap). Degeneracy (g_i) enters as a multiplier:
  p_i ~ g_i * exp(-dE_i / kT).
- ``collapse_degeneracies`` groups frames that are degenerate within a
  tolerance instead of collapsing them silently.

No binary is required: these are pure reading/weighting helpers.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence

#: 1 hartree (Eh) in kcal/mol — the exact atomic-unit conversion value.
KCAL_PER_EH = 627.5094740631

#: R in kcal/(mol*K), for dE in kcal/mol in the Boltzmann factor.
R_KCAL_MOL_K = 1.9872036e-3

#: Tolerance under which two frames count as degenerate (kcal/mol).
DEGENERACY_TOLERANCE_KCAL = 0.01


def parse_ensemble(path: Path) -> List[Dict[str, object]]:
    """Read a multi-frame XYZ; energy taken from the comment line.

    Accepts CREST comment shapes like ``energy: -12.34 Eh`` and bare
    numbers (CREST also writes plain floats). Raises ValueError when no
    energy is readable in a comment line, or OSError if the file is
    missing.
    """
    frames: List[Dict[str, object]] = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        try:
            natom = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(
                f"line {i + 1}: expected the atom count of a frame, "
                f"got {lines[i].strip()!r}"
            ) from exc
        if natom <= 0:
            raise ValueError(f"line {i + 1}: atom count must be positive, got {natom}")
        comment = lines[i + 1]
        atoms = lines[i + 2 : i + 2 + natom]
        if len(atoms) < natom:
            raise ValueError(
                f"frame at line {i + 1}: declares {natom} atoms but the "
                "file ends after " + str(len(atoms))
            )
        frames.append(
            {
                "n_atoms": natom,
                "energy_Eh": _energy_from_comment(comment, i + 2),
                "symbols": [a.split()[0] for a in atoms],
                "positions": [
                    tuple(float(x) for x in a.split()[1:4]) for a in atoms
                ],
            }
        )
        i += 2 + natom
    if not frames:
        raise ValueError("no conformer frames found in " + str(path))
    return frames


def _energy_from_comment(comment: str, line_no: int) -> float:
    tokens = comment.split()
    method = None
    for j, tok in enumerate(tokens):
        if tok.lower() == "energy:" and j + 1 < len(tokens):
            value = float(tokens[j + 1])
            # "energy: -12.34 Eh" or "energy: -12.34" (Eh is CREST's unit)
            return value
        if tok.lower().startswith("gfn") or tok.lower().startswith("xtb"):
            method = tok.lower()
        if tok.lower() in ("r2scan-3c", "b3lyp", "wb97x", "pbe0"):
            method = tok.lower()
    # fall back: bare number comment
    for tok in tokens:
        try:
            value = float(tok)
        except ValueError:
            continue
        if method is not None:
            return value
    raise ValueError(
        f"line {line_no}: no readable energy in comment {comment.strip()!r}"
    )


def relative_kcal(frames: Sequence[Dict[str, object]]) -> List[float]:
    """Relative energies in kcal/mol against the lowest frame."""
    energies = [float(f["energy_Eh"]) for f in frames]
    lowest = min(energies)
    return [(e - lowest) * KCAL_PER_EH for e in energies]


def boltzmann_weights(
    frames: Sequence[Dict[str, object]], temperature_K: float
) -> List[float]:
    """Boltzmann weights at ``temperature_K``; degeneracy enters as g_i.

    Raises ValueError when the frames carry energies from more than one
    level of theory (tagged via the ``method`` key) or when the
    temperature is not positive.
    """
    if temperature_K <= 0:
        raise ValueError(f"temperature must be positive, got {temperature_K}")
    methods = {f.get("method") for f in frames if f.get("method") is not None}
    if len(methods) > 1:
        raise ValueError(
            "the ensemble mixes energies from more than one level of "
            f"theory ({', '.join(sorted(str(m) for m in methods))}); "
            "refusing to weight across methods"
        )
    energies = [float(f["energy_Eh"]) for f in frames]
    degeneracies = [int(f.get("degeneracy", 1)) for f in frames]
    lowest = min(energies)
    rt = R_KCAL_MOL_K * temperature_K
    raw = [
        g * math.exp(-((e - lowest) * KCAL_PER_EH) / rt)
        for e, g in zip(energies, degeneracies)
    ]
    total = sum(raw)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Boltzmann normalisation failed (total " + str(total) + ")")
    return [w / total for w in raw]


def collapse_degeneracies(
    frames: Sequence[Dict[str, object]],
    tolerance_kcal: float = DEGENERACY_TOLERANCE_KCAL,
) -> List[Dict[str, object]]:
    """Group frames degenerate within ``tolerance_kcal``; groups stay
    ordered by energy and carry ``degeneracy`` = group size.
    """
    ordered = sorted(frames, key=lambda f: float(f["energy_Eh"]))
    groups: List[Dict[str, object]] = []
    for frame in ordered:
        e = float(frame["energy_Eh"])
        if groups and (
            (e - float(groups[-1]["energy_Eh"])) * KCAL_PER_EH <= tolerance_kcal
        ):
            groups[-1]["degeneracy"] = int(groups[-1].get("degeneracy", 1)) + 1
        else:
            new = dict(frame)
            new["degeneracy"] = 1
            groups.append(new)
    return groups
