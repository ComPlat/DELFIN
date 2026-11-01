"""TS -> GOAT: reactive-center detection and constrained GOAT input.

After an OptTS calculation DELFIN can propose which atoms belong to the
reactive centre -- the atoms directly involved in the chemical
transformation -- and fix them with ORCA cartesian constraints for a
GOAT transition-state conformer search (ORCA 6 manual, "Transition State
Conformers": keep the reacting atoms fixed so GOAT explores the
conformational space *around* the TS instead of relaxing it to a
minimum).

Everything here is a pure library: no widgets, no ORCA invocation.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

logger = __import__("logging").getLogger(__name__)

# ORCA's own bond criterion (GOAT manual): a bond exists when the distance
# is below 1.3 * (sum of covalent radii).
BOND_SCALE = 1.3

# Covalent radii (Angstrom), Cordero et al. 2008 -- only the elements that
# appear in the organic/organometallic test corpus; unknown elements fall
# back to a generous default so no pair is silently treated as unbonded.
COVALENT_RADII = {
    "H": 0.31, "D": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96,
    "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76,
    "Ti": 1.60, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Br": 1.20, "I": 1.39,
    "Pd": 1.39, "Ag": 1.45, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
}
_DEFAULT_RADIUS = 1.5


# --------------------------------------------------------------------- data

@dataclass
class TSFiles:
    """The files one OptTS calculation offers, matched by atom content."""

    out_path: Path | None = None
    inp_path: Path | None = None
    xyz_path: Path | None = None
    hess_path: Path | None = None
    reactant_path: Path | None = None
    product_path: Path | None = None
    neb_dir: Path | None = None
    problems: list[str] = field(default_factory=list)

    @property
    def complete(self) -> bool:
        return all((self.out_path, self.inp_path, self.xyz_path, self.hess_path))


@dataclass
class TSStatus:
    """Whether the selected .out is a usable converged TS."""

    is_ts: bool = False
    converged: bool = False
    n_imaginary: int = 0
    imaginary_freq: float | None = None
    #: (mode index, frequency) of every imaginary mode in the FINAL
    #: frequency block, most negative first -- drives the mode picker.
    imaginary_modes: tuple[tuple[int, float], ...] = ()
    charge: int | None = None
    mult: int | None = None
    elements: tuple[str, ...] = ()
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def usable(self) -> bool:
        """OptTS + convergence + at least one imaginary mode.

        More than one imaginary mode no longer blocks the workflow: a
        very small extra negative frequency (numerical noise, a flat
        torsion) is common in practice and the user can pick the reactive
        mode explicitly.  0 imaginary modes still blocks, with a clear
        message; >1 adds a warning instead of a rejection.
        """
        return self.is_ts and self.converged and self.n_imaginary >= 1

    @property
    def auto_mode_index(self) -> int | None:
        """The reactive mode auto-selected when several are imaginary.

        If exactly one mode is CLEARLY negative and the others are only
        very small negatives, that mode is pre-selected and the small
        ones become warnings.  Threshold: a secondary mode counts as
        noise below |freq| < 50 cm-1 (about the magnitude of numerical
        Hessian noise and clearly below any chemical TS mode, which sits
        in the several-hundred cm-1 range); deliberate and documented so
        it can be discussed in review.
        """
        if self.n_imaginary <= 1:
            return self.imaginary_modes[0][0] if self.imaginary_modes else None
        noise_floor_cm1 = 50.0
        clear = [
            (i, f) for i, f in self.imaginary_modes if abs(f) >= noise_floor_cm1
        ]
        if len(clear) == 1:
            self.warnings.append(
                "multiple imaginary modes in the final frequency block: "
                + ", ".join(f"mode {i} ({f:.1f} cm-1)"
                            for i, f in self.imaginary_modes)
                + " -- modes below |50| cm-1 treated as near-flat noise; "
                "the clearly reactive mode is pre-selected, please verify."
            )
            return clear[0][0]
        # several clearly negative modes: ambiguous, user must choose
        self.warnings.append(
            "several clearly imaginary modes in the final frequency block: "
            + ", ".join(f"mode {i} ({f:.1f} cm-1)"
                        for i, f in self.imaginary_modes)
            + " -- no automatic choice is safe; please select the reactive "
            "mode explicitly below."
        )
        return None


@dataclass
class ReactiveCenter:
    """The proposed constraint atoms, with the evidence for the choice."""

    indices: tuple[int, ...] = ()
    elements: tuple[str, ...] = ()
    reasons: list[str] = field(default_factory=list)
    reliability: str = "low"  # "high" | "medium" | "low"
    ambiguous: bool = False
    # per-atom score contributions, 0-based index -> labelled evidence
    scores: dict[int, float] = field(default_factory=dict)
    formed_bonds: tuple[tuple[int, int], ...] = ()
    broken_bonds: tuple[tuple[int, int], ...] = ()


# --------------------------------------------------------------- output scan

_CONVERGED_RE = re.compile(r"THE OPTIMIZATION HAS CONVERGED", re.IGNORECASE)
_HURRAY_RE = re.compile(r"HURRAY\s+.*\s+THE OPTIMIZATION HAS CONVERGED", re.IGNORECASE)
_CHARGE_RE = re.compile(r"^\s*%?charge\s+(\-?\d+)", re.IGNORECASE | re.MULTILINE)
_MULT_RE = re.compile(r"^\s*%?mult(?:iplicity)?\s+(\d+)", re.IGNORECASE | re.MULTILINE)
_XYZ_TAG_RE = re.compile(r"^\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE)
_IMAG_RE = re.compile(r"^\s*(\d+):\s*(-?\d+(?:\.\d+)?)\s+cm\*\*-1")


def read_inp_charge_mult(inp_path: Path) -> tuple[int | None, int | None]:
    """Charge/multiplicity from an ORCA input's ``* xyz N M`` tag."""
    try:
        text = Path(inp_path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None, None
    for line in text.splitlines():
        m = _XYZ_TAG_RE.search(line)
        if m:
            return int(m.group(1)), int(m.group(2))
    charge = _CHARGE_RE.search(text)
    mult = _MULT_RE.search(text)
    return (
        int(charge.group(1)) if charge else None,
        int(mult.group(1)) if mult else None,
    )


def read_ts_status(out_path: Path) -> TSStatus:
    """Is this .out an OptTS job, did it converge, how many imaginary modes?"""
    status = TSStatus()
    out_path = Path(out_path)
    try:
        text = out_path.read_text(encoding="utf-8", errors="ignore")
    except OSError as exc:
        status.reasons.append(f"output not readable: {exc}")
        return status

    status.is_ts = bool(re.search(r"\bOptTS\b", text))
    if not status.is_ts:
        status.reasons.append("no OptTS keyword in the output")

    status.converged = bool(_CONVERGED_RE.search(text) or _HURRAY_RE.search(text))
    if status.is_ts and not status.converged:
        status.reasons.append("TS optimization did not report convergence")

    # Imaginary modes: ONLY the last complete frequency block of the final
    # structure counts.  ORCA prints one such block per frequency job, and
    # earlier (historical) blocks along the OptTS trajectory also contain
    # imaginary modes -- a global scan would count those stale listings and
    # wrongly reject an otherwise fine TS.  The last block is the one whose
    # mode index can be looked up in the .hess, so prefer the .hess data
    # when present: it is the structured source for the final geometry.
    last_imag: list[tuple[int, float]] = []
    block_imag: list[tuple[int, float]] = []
    block_indices: list[int] = []
    for line in text.splitlines():
        m = _IMAG_RE.match(line)
        if m:
            block_indices.append(int(m.group(1)))
            if "***imaginary mode***" in line:
                block_imag.append((int(m.group(1)), float(m.group(2))))
            continue
        # a line that is not a frequency entry closes any open block
        if block_indices:
            last_imag, block_imag, block_indices = block_imag, [], []
    if block_indices:
        last_imag = block_imag
    status.n_imaginary = len(last_imag)
    if last_imag:
        status.imaginary_freq = min(freq for _, freq in last_imag)
        status.imaginary_modes = tuple(last_imag)
    if status.is_ts and status.converged and status.n_imaginary == 0:
        status.reasons.append(
            "converged TS reported no imaginary mode in the output "
            "(frequencies may not have been recomputed at the final geometry)"
        )
    if status.n_imaginary > 1:
        status.reasons.append(
            f"{status.n_imaginary} imaginary modes -- not a plain TS"
        )
    return status


def read_xyz(path: Path) -> list[tuple[str, float, float, float]] | None:
    """First frame of an .xyz as (element, x, y, z); None if unparseable.

    ORCA writes a comment line that may itself contain 4+ tokens
    ("Coordinates from ORCA-job ... E -401.250208837948" -- the energy is
    part of the comment), so a non-atom line must be SKIPPED, not treated
    as a parse error.  A line counts as an atom line only when the second
    through fourth tokens parse as floats; anything else is a comment.
    """
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None
    if not lines:
        return None
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return None
    atoms: list[tuple[str, float, float, float]] = []
    for row in lines[1:]:
        if len(atoms) >= n:
            break
        parts = row.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue  # comment / non-atom line -- skip, do not abort
        atoms.append((parts[0], x, y, z))
    return atoms if len(atoms) == n else None


def _same_composition(a: Sequence, b: Sequence) -> bool:
    """Same atom count AND same element order (ORCA is order-sensitive)."""
    if len(a) != len(b):
        return False
    return all(
        str(x[0]).capitalize() == str(y[0]).capitalize() for x, y in zip(a, b)
    )


def _composition_mismatch(a: Sequence, b: Sequence) -> str:
    """Human-readable first difference between two atom orders.

    Coordinates are deliberately NOT compared: a Hessian stores positions
    in Bohr, an .xyz in Angstrom, and a small numeric mismatch must not
    reject a genuinely matching geometry.  Only count and element order
    (order-sensitive, as ORCA requires) are compared.
    """
    if len(a) != len(b):
        return (f"atom count {len(a)} != Hessian {len(b)}")
    for i, (x, y) in enumerate(zip(a, b)):
        if str(x[0]).capitalize() != str(y[0]).capitalize():
            return (f"element order differs at index {i}: "
                    f"xyz has {x[0]}, Hessian has {y[0]}")
    return "identical"


# -------------------------------------------------------------- file finding

_HESS_SUFFIXES = (".hess",)
_EDUKT_RE = re.compile(r"(?:edukt|reactant|reac|start|endpoint[_-]?0)", re.IGNORECASE)
_PRODUKT_RE = re.compile(r"(?:produkt|product|prod|endpoint[_-]?1)", re.IGNORECASE)


def collect_ts_files(calc_dir: Path, out_path: Path | None = None) -> TSFiles:
    """Locate .out/.inp/.xyz/.hess (+ reactant/product) for one OptTS run.

    Matching is by atom content, not by name alone: a candidate geometry
    is accepted only when its element sequence matches the Hessian's
    atoms (or, without a Hessian yet, the output's own).  Ambiguity is
    reported in ``problems`` instead of silently picking one file.
    """
    calc_dir = Path(calc_dir)
    found = TSFiles(out_path=out_path)

    # -- the output's own base name first; different ORCA runs name things
    #    differently (<base>.out/<base>.inp/<base>.hess, orca_ts_*_TS_*.out...)
    base_stems: list[str] = []
    if out_path is not None:
        base_stems.append(out_path.stem)

    def _by_stem(suffix: str) -> Path | None:
        for stem in base_stems:
            cand = calc_dir / f"{stem}{suffix}"
            if cand.is_file():
                return cand
        return None

    found.inp_path = _by_stem(".inp")
    found.hess_path = _by_stem(".hess")

    # -- geometry: prefer an optimised-structure file that matches the
    #    Hessian's element order
    reference = None
    xyz_rejections: list[str] = []
    if found.hess_path is not None:
        from delfin.imag import _hessian_atoms  # noqa: PLC0415 -- reuse of the parser
        try:
            hess_text = found.hess_path.read_text(encoding="utf-8", errors="ignore")
            reference = _hessian_atoms(hess_text.splitlines())
        except (ValueError, OSError):
            reference = None

    xyz_candidates: list[Path] = []
    if found.inp_path is not None:
        # ORCA writes the optimised structure as <base>.xyz next to the input.
        cand = found.inp_path.with_suffix(".xyz")
        if cand.is_file():
            xyz_candidates.append(cand)
    xyz_candidates.extend(sorted(calc_dir.glob("*.xyz")))
    # unique, order-preserving
    seen: set = set()
    xyz_candidates = [p for p in xyz_candidates
                      if not (str(p) in seen or seen.add(str(p)))]

    for cand in xyz_candidates:
        atoms = read_xyz(cand)
        if atoms is None:
            xyz_rejections.append(f"{cand.name}: unparseable (no atom lines)")
            continue
        if reference is not None and not _same_composition(atoms, reference):
            mismatch = _composition_mismatch(atoms, reference)
            xyz_rejections.append(f"{cand.name}: {mismatch}")
            continue
        # Do not mistake reactant/product endpoint files for the TS itself.
        if _EDUKT_RE.search(cand.stem) or _PRODUKT_RE.search(cand.stem):
            if _EDUKT_RE.search(cand.stem) and found.reactant_path is None:
                found.reactant_path = cand
            if _PRODUKT_RE.search(cand.stem) and found.product_path is None:
                found.product_path = cand
            continue
        found.xyz_path = cand
        break

    # -- reactant / product endpoints (NEB naming and Ed/Prt naming)
    if found.reactant_path is None or found.product_path is None:
        for cand in xyz_candidates:
            atoms = read_xyz(cand)
            if atoms is None:
                continue
            if reference is not None and not _same_composition(atoms, reference):
                continue
            if _EDUKT_RE.search(cand.stem) and found.reactant_path is None:
                found.reactant_path = cand
            elif _PRODUKT_RE.search(cand.stem) and found.product_path is None:
                found.product_path = cand

    # -- NEB leftovers (e.g. *.final.interp, TS_*_NEB* folders are parents)
    neb_hits = sorted(calc_dir.glob("*.final.interp"))
    if neb_hits:
        found.neb_dir = calc_dir

    # -- report what is missing
    if found.hess_path is None:
        found.problems.append("no Hessian (.hess) found -- imaginary-mode "
                              "analysis is not possible without it")
    if found.xyz_path is None:
        msg = ("no matching TS geometry (.xyz) found "
               "(atom order must match the Hessian)")
        if xyz_rejections:
            msg += " -- candidates checked: " + "; ".join(xyz_rejections)
        found.problems.append(msg)
    if found.inp_path is None:
        found.problems.append("no ORCA input (.inp) found -- charge and "
                              "multiplicity must be set manually")
    if found.reactant_path is None and found.product_path is None:
        found.problems.append("no reactant/product endpoint geometry found "
                              "-- connectivity-change evidence unavailable")
    return found


# --------------------------------------------------- reactive-centre scoring

def _bond_list(atoms: Sequence[tuple[str, float, float, float]],
               scale: float = BOND_SCALE) -> set:
    """Pairs (i, j), i < j, bonded by ORCA's distance criterion."""
    bonds = set()
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            cut = scale * (
                COVALENT_RADII.get(str(atoms[i][0]).capitalize(), _DEFAULT_RADIUS)
                + COVALENT_RADII.get(str(atoms[j][0]).capitalize(), _DEFAULT_RADIUS)
            )
            dx = atoms[i][1] - atoms[j][1]
            dy = atoms[i][2] - atoms[j][2]
            dz = atoms[i][3] - atoms[j][3]
            if math.sqrt(dx * dx + dy * dy + dz * dz) < cut:
                bonds.add((i, j))
    return bonds


def connectivity_changes(ts_atoms, reactant_atoms=None, product_atoms=None):
    """(formed, broken) bond pairs (0-based) from endpoint -> TS -> endpoint."""
    formed: set = set()
    broken: set = set()
    ts_bonds = _bond_list(ts_atoms)
    if reactant_atoms is not None and len(reactant_atoms) == len(ts_atoms):
        reac_bonds = _bond_list(reactant_atoms)
        broken |= reac_bonds - ts_bonds
    if product_atoms is not None and len(product_atoms) == len(ts_atoms):
        prod_bonds = _bond_list(product_atoms)
        formed |= prod_bonds - ts_bonds
        # a bond present in the TS but gone in the product is also a
        # reacting contact (seen from the other side)
        broken |= ts_bonds - prod_bonds
    return formed, broken


def analyze_reactive_center(
    hess_path: Path,
    ts_atoms: Sequence[tuple[str, float, float, float]],
    reactant_atoms=None,
    product_atoms=None,
    *,
    mode_index: int | None = None,
    max_atoms: int = 8,
) -> ReactiveCenter:
    """Propose the reactive-centre atoms from the imaginary mode + endpoints.

    ``mode_index`` selects WHICH imaginary mode drives the analysis;
    None means the most negative one (the classical reactive mode).
    Scoring combines three independent pieces of evidence so that the
    choice is NOT driven by the largest Cartesian displacement alone
    (a torsional or umbrella-like motion moves many innocent atoms):

    1. **Mode displacement** -- |eigenvector| per atom of the most
       imaginary mode (scaled 0..1).
    2. **Bond changes** -- atoms of every bond formed or broken between
       the endpoints and the TS (strongest evidence; a proton transfer,
       bond break or ring closure shows up here even when its atoms move
       modestly).
    3. **Distance response along the mode** -- for each atom pair whose
       bond status differs between the endpoints and the TS, how much
       that distance changes when the structure is displaced along the
       mode in both directions (displacement *increases* a forming or
       *decreases* a breaking contact -> the mode really acts on it).
    """
    ts_atoms = list(ts_atoms)
    n = len(ts_atoms)
    formed, broken = connectivity_changes(ts_atoms, reactant_atoms, product_atoms)
    reacting_pairs = sorted(formed | broken)

    scores: dict[int, float] = {}
    reasons: list[str] = []

    # -- mode eigenvector (selected or most imaginary column) -------------
    mode_displacement: dict[int, float] = {}
    try:
        from delfin.imag import _normal_mode, imaginary_modes  # noqa: PLC0415
        if mode_index is None:
            modes = imaginary_modes(Path(hess_path), {"allow_imaginary_freq": 0.0})
            mode_index = modes[0][0] if modes else None
        if mode_index is not None:
            text = Path(hess_path).read_text(encoding="utf-8", errors="ignore")
            vector = _normal_mode(text.splitlines(), mode_index)
        norms = [
            math.sqrt(sum(v * v for v in vector[3 * a:3 * a + 3]))
            for a in range(n)
        ]
        biggest = max(norms) if norms else 0.0
        if biggest > 0:
            mode_displacement = {
                a: norms[a] / biggest for a in range(n) if norms[a] > 1e-6
            }
    except (ValueError, OSError) as exc:
        reasons.append(f"imaginary-mode eigenvector unavailable ({exc})")

    for a, rel in mode_displacement.items():
        scores[a] = scores.get(a, 0.0) + 1.0 * rel

    # -- bond changes: the reacting contacts and their atoms -------------
    reacting_atoms: set = set()
    for i, j in reacting_pairs:
        reacting_atoms.update((i, j))
        for a in (i, j):
            scores[a] = scores.get(a, 0.0) + 2.0
    if reacting_pairs:
        formed_txt = ", ".join(
            f"{ts_atoms[i][0]}{i}-{ts_atoms[j][0]}{j}" for i, j in sorted(formed)
        ) or "none"
        broken_txt = ", ".join(
            f"{ts_atoms[i][0]}{i}-{ts_atoms[j][0]}{j}" for i, j in sorted(broken)
        ) or "none"
        reasons.append(
            "connectivity change vs. endpoints: formed " + formed_txt
            + ", broken " + broken_txt
        )
    else:
        reasons.append(
            "no endpoint geometries matched the TS atom order -- "
            "connectivity-change evidence unavailable"
        )

    # -- distance response along the mode --------------------------------
    if mode_displacement and reacting_pairs:
        from delfin.imag import displaced_geometries  # noqa: PLC0415
        try:
            displaced = displaced_geometries(Path(hess_path), mode_index, 0.25)
        except (ValueError, OSError) as exc:
            displaced = None
            reasons.append(f"displacement along the mode unavailable ({exc})")
        if displaced:
            for i, j in reacting_pairs:
                d0 = _distance(ts_atoms[i], ts_atoms[j])
                response = max(
                    abs(_distance(displaced["pos"][i], displaced["pos"][j]) - d0),
                    abs(_distance(displaced["neg"][i], displaced["neg"][j]) - d0),
                )
                if response > 0.05:  # the mode actually opens/closes this contact
                    for a in (i, j):
                        scores[a] = scores.get(a, 0.0) + 1.5
            reasons.append(
                "distance response along the imaginary mode checked for every "
                "reacting contact"
            )

    # -- pick the atoms ---------------------------------------------------
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    # Everything that participates in a reacting bond is in; mode-only
    # atoms follow while their score stays clearly above the tail.
    selected = {a for a in reacting_atoms}
    threshold = 1.0  # one full-strength mode contribution
    for a, score in ranked:
        if len(selected) >= max_atoms:
            break
        if score >= threshold:
            selected.add(a)
    indices = tuple(sorted(selected))
    elements = tuple(str(ts_atoms[a][0]) for a in indices)

    if not indices:
        # last resort: the atoms that move most along the mode -- flagged
        indices = tuple(sorted(mode_displacement, key=lambda a: -mode_displacement[a])[:4])
        elements = tuple(str(ts_atoms[a][0]) for a in indices)
        reasons.append(
            "fallback: the four atoms with the largest mode displacement "
            "(no connectivity evidence was available)"
        )

    # -- reliability ------------------------------------------------------
    reliability = "low"
    if reacting_pairs and mode_displacement:
        reliability = "high" if len(reacting_pairs) >= 2 else "medium"
    elif reacting_pairs:
        reliability = "medium"
    ambiguous = (
        not reacting_pairs
        or (mode_displacement and len(selected) > max_atoms // 2)
        or len(selected) > max_atoms
    )
    if ambiguous:
        reasons.append(
            "the proposal is ambiguous -- please check and adjust the "
            "selection in the viewer"
        )
    return ReactiveCenter(
        indices=indices,
        elements=elements,
        reasons=reasons,
        reliability=reliability,
        ambiguous=ambiguous,
        scores={a: round(s, 3) for a, s in sorted(scores.items())},
        formed_bonds=tuple(sorted(formed)),
        broken_bonds=tuple(sorted(broken)),
    )


def _distance(a, b) -> float:
    return math.sqrt(
        (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2
    )


# ------------------------------------------------------------- GOAT input

def render_constraints_block(indices: Sequence[int],
                             elements: Sequence[str]) -> str:
    """The ORCA ``%geom Constraints`` block, 0-based, cartesian 'C' fixed."""
    lines = ["%geom", "  Constraints"]
    for idx, el in zip(indices, elements):
        lines.append(f"    {{ {el} {int(idx)} C }}")
    lines.append("  end")
    lines.append("end")
    return "\n".join(lines)


def build_ts_goat_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    constraint_indices: Sequence[int],
    charge: int = 0,
    mult: int = 1,
    method: str = "XTB",
    functional: str = "",
    basis: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "goat_ts",
) -> Path:
    """Write the constrained GOAT input for a TS conformer search.

    Delegates to the existing ``build_goat_input`` (OPI-based) and then
    inserts the ``%geom Constraints`` block before the coordinate
    section, because the OPI goat block does not model ORCA geometry
    constraints.  Returns the written input path.
    """
    from delfin.tools.adapters.orca_goat_opi import build_goat_input  # noqa: PLC0415

    inp_path = build_goat_input(
        Path(work_dir),
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        method=method,
        goat_options=None,
        extra_keywords=None,
        functional=functional,
        basis=basis,
        aux_basis=aux_basis,
        dispersion=dispersion,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    block = render_constraints_block(
        constraint_indices,
        [str(elements[i]) for i in constraint_indices],
    )
    text = inp_path.read_text(encoding="utf-8")
    # The coordinate section starts at the '*' card; ORCA accepts %geom
    # blocks anywhere before it.
    star = re.search(r"^\s*\*", text, re.MULTILINE)
    if star is None:
        raise ValueError("generated GOAT input has no coordinate section")
    patched = text[:star.start()] + block + "\n\n" + text[star.start():]
    inp_path.write_text(patched, encoding="utf-8")
    return inp_path
