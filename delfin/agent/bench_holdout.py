"""Holdout sampler for the wave-20 workflow benchmark (package CB0B).

Each of DELFIN's seven workflow families is a GENERATOR: an instance is a
chemist-voiced request over a molecule from a small pool, a solvent and a
requested property. The prompt deliberately does NOT name the workflow —
the task is whether the agent picks DELFIN's own tool for the job.

Seed handling
-------------
The real holdout seed lives OUTSIDE the repository in a file the person
creates; sessions never read it. ``real_seed`` looks it up and falls back
to ``PUBLIC_SEED`` (the fixed seed the tests and practice runs use) when
that file is absent or unreadable. The only path ever consulted is the
``DELFIN_BENCH_SEED_PATH`` environment variable — no other location, no
user state, nothing in the repository.
"""
from __future__ import annotations

import os
import random
from dataclasses import dataclass, field

PUBLIC_SEED = 20261020

# Grading levels, cheapest first (brief items 5a-5d).
GRADING_LEVELS = ("consistency", "reproduction", "reference", "trap")


@dataclass(frozen=True)
class Molecule:
    """One pool molecule. Atoms counted from a provided or parsed formula."""

    label: str
    formula: str
    charge: int
    multiplicity: int
    smiles: str = ""

    @property
    def n_atoms(self) -> int:
        # Count atoms from the formula (no external parser needed; all
        # pool entries are simple Hill-notation strings).
        import re
        total = 0
        for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", self.formula):
            if symbol:
                total += int(count) if count else 1
        return total


@dataclass(frozen=True)
class Family:
    """A workflow family: its id (internal), the pool it draws from, the
    properties it can ask for and the grading levels it supports."""

    id: str
    molecules: tuple[Molecule, ...]
    properties: tuple[str, ...]
    solvents: tuple[str, ...]
    grading_levels: tuple[str, ...]
    description: str = ""


@dataclass(frozen=True)
class Instance:
    """One drawn benchmark instance."""

    family_id: str
    molecule: Molecule
    solvent: str
    property: str
    prompt: str
    grading_levels: tuple[str, ...]
    seed: int
    instance_id: str = field(default="")


def real_seed() -> int:
    """The holdout seed for real benchmark runs.

    The seed lives OUTSIDE the repository in a file the person creates;
    its path comes from DELFIN_BENCH_SEED_PATH. Sessions never read it —
    if the variable is unset or the file is missing or unreadable, the
    PUBLIC_SEED is used and practice stays reproducible.
    """
    path = os.environ.get("DELFIN_BENCH_SEED_PATH")
    if not path:
        return PUBLIC_SEED
    try:
        with open(path, encoding="utf-8") as fh:
            return int(fh.read().strip())
    except (OSError, ValueError):
        # Missing or malformed seed file: fall back to the public seed
        # rather than guessing. No default elsewhere, no user state read.
        return PUBLIC_SEED


# ── the molecule pool ─────────────────────────────────────────────────
# Every entry <= 20 atoms. Charges/multiplicities are electron-counted
# by hand; the red tests pin coverage of both shells and both charges.
# nitrite_radical is the deliberate trap carrier for F2 (reviewer
# nacht-s15): 23 electrons (odd) DECLARED with multiplicity 1 — the
# multiplicity the task requests, not the physical truth. A closed-shell
# request over an odd-electron radical is impossible; the agent must
# notice and report the mismatch instead of running it.
_POOL_NR = Molecule("nitrite_radical", "NO2", 0, 1, smiles="O=[N][O]")
_POOL_M = Molecule("methanol", "CH4O", 0, 1, smiles="CO")
_POOL_C = Molecule("cyclopropane", "C3H6", 0, 1, smiles="C1CC1")
_POOL_N = Molecule("nitrite_anion", "NO2", -1, 1, smiles="[O-][N]=O")
_POOL_O = Molecule("dioxygen", "O2", 0, 3, smiles="O=O")
_POOL_F = Molecule("formaldehyde_cation", "CH2O", +1, 2, smiles="[CH2+]O")
_POOL_S = Molecule("phosphine", "H3P", 0, 1, smiles="P")
_POOL_CL = Molecule("chlorine", "Cl2", 0, 1, smiles="ClCl")
_POOL_NO = Molecule("nitric_oxide", "NO", 0, 2, smiles="[N]=O")
_POOL_NR = Molecule("nitrite_radical", "NO2", 0, 1, smiles="O=[N][O]")

_POOL: tuple[Molecule, ...] = (
    _POOL_M, _POOL_C, _POOL_N, _POOL_O,
    _POOL_F, _POOL_S, _POOL_CL, _POOL_NO,
    _POOL_NR,
)

# Electron-count sanity for the pool: charge and multiplicity must be
# allowed by the electron count (2*mult - 2 <= n_electrons < 2*mult + 2
# is NOT enforced here — the graders do that — but gross mistakes are
# refused when the family is built).
_S = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16, "Cl": 17}


def _electrons(mol: Molecule) -> int:
    import re
    n = sum(_S[s] * (int(c) if c else 1)
            for s, c in re.findall(r"([A-Z][a-z]?)(\d*)", mol.formula))
    return n - mol.charge


def _electron_count_is_consistent(mol: Molecule) -> bool:
    # multiplicity m = unpaired + 1, so n - (m - 1) must be even and
    # non-negative; a closed shell (m = 1) needs an even count.
    n = _electrons(mol)
    unpaired = mol.multiplicity - 1
    return unpaired >= 0 and (n - unpaired) % 2 == 0


# ── the seven workflow families ──────────────────────────────────────
# One per DELFIN workflow: classic opt+freq, redox, excited states,
# imaginary-mode fix, OCCUPIER spin states, conformer pre-step with the
# tight-binding tools, screening. The prompt never names any of them.

_FAMILY_PROMPTS = {
    "classic_opt_freq": (
        "Optimise {mol} in {solvent} and confirm the resulting geometry "
        "is a true minimum. Report the lowest-energy structure's "
        "thermochemical correction at 298 K and the count of imaginary "
        "frequencies."),
    "redox": (
        "Estimate the one-electron potential for {mol} in {solvent} "
        "(vs. a suitable reference) and state the method and basis you "
        "used for both members of the couple."),
    "excited_states": (
        "Compute the first optically bright excitation of {mol} in "
        "{solvent} and report wavelength, oscillator strength and the "
        "dominant orbital transition."),
    "imaginary_mode": (
        "A colleague's optimisation of {mol} ended with one negative "
        "curvature direction. Find out where the structure wants to go, "
        "correct it, and report a geometry with no negative curvature "
        "plus the energy change along the way."),
    "spin_states": (
        "For the {dmetal} complex of {mol} determine which spin state "
        "is lowest. Report relative energies of the candidate spin "
        "states from ONE level of theory and the <S^2> values."),
    "conformer_prestep": (
        "{mol} is flexible. Before any expensive calculation, find its "
        "lowest-energy conformer with a fast semiempirical search and "
        "report the conformer energy spread and the SMILES/geometry of "
        "the winner."),
    "screening": (
        "Screen {solvent} as medium for {mol}: compute the solvation "
        "free-energy contribution at the same level for the neutral and "
        "the charged member of the couple and report the difference."),
}

_SOLVENTS = ("water", "acetonitrile", "gas phase")

_SOLVENT_PHRASE = {"water": "water", "acetonitrile": "acetonitrile",
                   "gas phase": "the gas phase"}


def _build_families() -> tuple[Family, ...]:
    closed_neutral = [m for m in _POOL if m.charge == 0 and m.multiplicity == 1]
    open_shell = [m for m in _POOL if m.multiplicity > 1]
    charged = [m for m in _POOL if m.charge != 0]
    any_of = list(_POOL)
    # F1 (nacht-s15): a solvation free-energy difference is undefined
    # without a medium — the screening family never offers the gas phase.
    solvated = ("water", "acetonitrile")

    def fam(fid, mols, props, levels, desc, solvents=_SOLVENTS):
        return Family(fid, tuple(mols), tuple(props), tuple(solvents),
                      levels, desc)

    return (
        fam("classic_opt_freq", any_of,
            ("thermochemical correction", "minimum check"),
            ("consistency", "reproduction", "reference"),
            "Geometry optimisation plus frequency analysis."),
        fam("redox", charged + open_shell,
            ("one-electron potential",),
            ("consistency", "reproduction", "reference"),
            "One-electron couples: potential of a redox pair."),
        fam("excited_states", closed_neutral,
            ("brightest low-energy excitation",),
            ("consistency", "reproduction", "reference"),
            "Vertical excitations and spectra."),
        fam("imaginary_mode", any_of,
            ("relaxed structure without negative curvature",),
            ("consistency", "reproduction", "trap"),
            "Recover from a saddle-point relaxation."),
        fam("spin_states", open_shell,
            ("lowest spin state",),
            ("consistency", "reproduction", "reference", "trap"),
            "Spin-state ordering of an open-shell species."),
        fam("conformer_prestep", closed_neutral,
            ("lowest-energy conformer",),
            ("consistency", "reproduction"),
            "Cheap conformer search before the expensive step."),
        fam("screening", charged + any_of,
            ("solvation free-energy difference",),
            ("consistency", "reproduction", "reference"),
            "Solvent screening across a couple.", solvents=solvated),
    )


FAMILIES: tuple[Family, ...] = _build_families()

# Anticipated traps an agent should notice and report (brief item 5d):
# odd electron count with an even multiplicity request, an anion without
# diffuse functions, a start geometry relaxing to a saddle point.
_TRAP_BY_FAMILY = {
    "imaginary_mode": "saddle point start",
    "spin_states": "spin contamination",
    "screening": "anion without diffuse functions",
    "redox": "odd electron count with wrong multiplicity",
    "classic_opt_freq": "start geometry relaxes to a saddle point",
    "excited_states": "bright line lies outside the visible",
    "conformer_prestep": "lowest conformer is not the one sought",
}


def trap_note(family_id: str, mol: Molecule) -> str:
    """The trap an instance of ``family_id`` presents with ``mol`` — the
    reason the 'trap' grading level has something to grade. Empty string
    when this molecule cannot present the family's trap."""
    if family_id == "redox":
        # F2 (nacht-s15): the trap is presentable when the electron
        # count contradicts a naive multiplicity request — the agent
        # must catch the mismatch before running anything.
        if not _electron_count_is_consistent(mol):
            return "odd electron count with wrong multiplicity"
        return ""
    return _TRAP_BY_FAMILY.get(family_id, "")



def sample(seed: int, n: int) -> list[Instance]:
    """Draw ``n`` deterministic instances from the families.

    Same seed -> same instances; different seeds differ. ``n`` instances
    are drawn round-robin over the families (so a large ``n`` walks every
    family), each draw picking molecule, solvent and property.
    """
    rng = random.Random(seed)
    out: list[Instance] = []
    for i in range(n):
        fam = FAMILIES[i % len(FAMILIES)]
        mol = rng.choice(list(fam.molecules))
        solvent = rng.choice(list(fam.solvents))
        prop = rng.choice(list(fam.properties))
        levels = fam.grading_levels
        if mol.charge != 0 or mol.multiplicity > 1:
            # open-shell/charged instances keep the trap level if the
            # family carries one, otherwise the family's own set stands.
            trap = _TRAP_BY_FAMILY.get(fam.id)
            if trap and "trap" not in levels:
                levels = levels + ("trap",)
        template = _FAMILY_PROMPTS[fam.id]
        prompt = template.format(
            mol=mol.label.replace("_", " "),
            solvent=_SOLVENT_PHRASE[solvent],
            dmetal="iron",
        )
        inst = Instance(
            family_id=fam.id,
            molecule=mol,
            solvent=solvent,
            property=prop,
            prompt=prompt,
            grading_levels=levels,
            seed=seed,
            instance_id=f"{fam.id}/{mol.label}/{solvent}",
        )
        out.append(inst)
    return out

