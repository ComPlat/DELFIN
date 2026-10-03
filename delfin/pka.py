"""pKa estimation from isodesmic proton-transfer cycles.

Computes the pKa of a mono-protic acid in solution relative to a
reference acid of known experimental pKa:

    pKa(target) = pKa(ref) + (dG_deprot(target) - dG_deprot(ref))
                          / (ln(10) * R * T)

with dG_deprot = G(conjugate base) - G(acid).  All four species are
treated identically (same functional, basis set, solvation model and
concentration standard state), so systematic errors cancel inside the
isodesmic difference to first order -- this is why no absolute
proton free energy is needed.

Standard state: both phases at 1 mol/L in solution.  When geometries and
vibrations are computed with a continuum solvation model (ORCA SMD or
CPCM, Manual 6.1.1 section 2.13.3), the 1 atm -> 1 M concentration
correction of the gas phase (+1.89 kcal/mol) cancels between the target
and the reference acid, so no explicit correction is applied here.

References:
    Ho, J.; Coote, M. L. A universal approach for calculating the pKa of
    organic molecules. ChemBioChem 2009, 10, 1837-1843.
    Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. Universal solvation
    model (SMD). J. Phys. Chem. B 2009, 113, 6378-6396.

ORCA execution reuses :func:`delfin.orca.run_orca`; output parsing
reuses :func:`delfin.energies.find_gibbs_energy` and
:func:`delfin.energies.find_electronic_energy`.  Structure generation
reuses the DELFIN SMILES converter via
:func:`delfin.tadf_xtb._write_smiles_inputs`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from delfin.common.logging import get_logger
from delfin.energies import find_electronic_energy, find_gibbs_energy


def _write_smiles_inputs(smiles: str, label: str, workdir: Path):
    """Delegate to DELFIN's SMILES converter through tadf_xtb's helper."""
    from delfin.tadf_xtb import _write_smiles_inputs as _delfin_write

    return _delfin_write(smiles, label, workdir)


def run_orca(input_file_path, output_log, **kwargs):
    """Delegate to DELFIN's ORCA runner (resolved, CONTROL-overridden).

    Patch point for tests: monkeypatching ``delfin.pka.run_orca`` swaps
    out the ORCA execution for the whole pKa cycle.
    """
    from delfin.orca import run_orca as _delfin_run_orca

    return _delfin_run_orca(input_file_path, output_log, **kwargs)

logger = get_logger(__name__)

# ln(10) * R in kcal/(mol*K), CODATA 2018 molar gas constant.
R_KCAL_PER_MOL_K = 1.98720425864083e-3
PKA_KELVIN_DEFAULT = 298.15

#: Experimental pKa values in water at 25 C (I = 0 dilute), monoprotic.
#: HA = neutral acid, A = conjugate base (deprotonated, charge -1).
KNOWN_ACIDS: Dict[str, Dict[str, Any]] = {
    "acetic": {
        "pka": 4.756,
        "ha_smiles": "CC(=O)O",
        "a_smiles": "CC(=O)[O-]",
    },
    "formic": {
        "pka": 3.745,
        "ha_smiles": "OC=O",
        "a_smiles": "[O-]C=O",
    },
    "phenol": {
        "pka": 9.99,
        "ha_smiles": "Oc1ccccc1",
        "a_smiles": "[O-]c1ccccc1",
    },
    "cyanoacetic": {
        "pka": 2.45,
        "ha_smiles": "N#CCC(=O)O",
        "a_smiles": "N#CCC(=O)[O-]",
    },
    "benzoic": {
        "pka": 4.202,
        "ha_smiles": "O=C(O)c1ccccc1",
        "a_smiles": "O=C([O-])c1ccccc1",
    },
}

#: The reference acid every target is anchored to.
REFERENCE_ACID = "acetic"

# CODATA 2018: 1 Hartree = 627.5094740631 kcal/mol.
_HARTREE_TO_KCAL = 627.5094740631

#: Execution defaults shared by every species of one cycle so the
#: isodesmic difference stays internally consistent.
_DEFAULT_FUNCTIONAL = "B3LYP"
_DEFAULT_BASIS = "def2-SVP"
_DEFAULT_SP_BASIS = "def2-TZVP"
_DEFAULT_SOLVENT = "water"

_HA = "HA"
_A = "A"


@dataclass(frozen=True)
class Species:
    """One species of a proton-transfer cycle.

    label: unique identifier used as file stem and output key
    charge: molecular charge (neutral acid 0, conjugate base -1)
    multiplicity: spin multiplicity (closed shell 1 by default)
    smiles: optional source structure; when given, DELFIN's SMILES
        converter generates the 3D coordinates
    """

    label: str
    charge: int
    multiplicity: int = 1
    smiles: Optional[str] = None


def pka_kcal_per_pk_unit(temperature_k: float = PKA_KELVIN_DEFAULT) -> float:
    """Free energy in kcal/mol that shifts the pKa by one unit: ln(10)*R*T."""
    return math.log(10.0) * R_KCAL_PER_MOL_K * temperature_k


def deprotonation_free_energy(gibbs_acid: float, gibbs_base: float) -> float:
    """dG_deprot = G(conjugate base) - G(acid) in Hartree."""
    return gibbs_base - gibbs_acid


def pka_from_cycle(
    gibbs: Mapping[str, Optional[float]],
    reference_pka: float,
    temperature_k: float = PKA_KELVIN_DEFAULT,
) -> float:
    """Anchor the target acid to the reference acid's experimental pKa.

    Args:
        gibbs: mapping with keys ``target_HA``, ``target_A``,
            ``reference_HA``, ``reference_A`` holding Gibbs free
            energies in Hartree (solution phase, same protocol for all
            four species).
        reference_pka: experimental pKa of the reference acid under the
            same conditions.
        temperature_k: temperature of the thermochemistry, default 298.15 K.

    Returns:
        The estimated pKa of the target acid.

    Raises:
        ValueError: when any required Gibbs energy is missing, so the
            cycle can never silently guess from partial data.
    """
    required = ("target_HA", "target_A", "reference_HA", "reference_A")
    missing = [key for key in required if gibbs.get(key) is None]
    if missing:
        raise ValueError(
            "pKa cycle incomplete: missing Gibbs energy for "
            + ", ".join(missing)
            + " ( Hartree ); refusing to guess"
        )
    d_deprot_target = gibbs["target_A"] - gibbs["target_HA"]
    d_deprot_reference = gibbs["reference_A"] - gibbs["reference_HA"]
    delta_kcal = (
        d_deprot_target - d_deprot_reference
    ) * _HARTREE_TO_KCAL
    logger.info(
        "Isodesmic pKa cycle: dG_deprot(target)=%.6f Eh, "
        "dG_deprot(reference)=%.6f Eh, delta=%.3f kcal/mol "
        "(%.3f kcal per pK unit at %.2f K)",
        d_deprot_target,
        d_deprot_reference,
        delta_kcal,
        pka_kcal_per_pk_unit(temperature_k),
        temperature_k,
    )
    return reference_pka + delta_kcal / pka_kcal_per_pk_unit(temperature_k)


def plan_species(acid: str, reference_acid: str = REFERENCE_ACID) -> List[Species]:
    """Return the four Species of the isodesmic cycle for one target acid.

    Order matters for callers: the two target species first, then the
    two reference species, all closed-shell singlets with charge 0 for
    the neutral acid and -1 for the conjugate base.
    """
    target = KNOWN_ACIDS.get(acid)
    if target is None:
        raise ValueError(
            f"Unknown acid '{acid}'; known acids: {sorted(KNOWN_ACIDS)}"
        )
    reference = KNOWN_ACIDS[reference_acid]
    return [
        Species(label=f"{acid}_{_HA}", charge=0, multiplicity=1,
                smiles=target["ha_smiles"]),
        Species(label=f"{acid}_{_A}", charge=-1, multiplicity=1,
                smiles=target["a_smiles"]),
        Species(label=f"{reference_acid}_{_HA}", charge=0, multiplicity=1,
                smiles=reference["ha_smiles"]),
        Species(label=f"{reference_acid}_{_A}", charge=-1, multiplicity=1,
                smiles=reference["a_smiles"]),
    ]


def build_opt_freq_input(
    xyz_text: str,
    charge: int,
    multiplicity: int,
    solvent: str,
    *,
    functional: str = _DEFAULT_FUNCTIONAL,
    basis: str = _DEFAULT_BASIS,
    pal: int = 8,
    maxcore: int = 4000,
    solvation_model: str = "SMD",
    extra_keywords: Optional[Iterable[str]] = None,
) -> str:
    """Write the ORCA Opt+Freq input that yields Gibbs energies in solution.

    Solvation via SMD by default (ORCA Manual 6.1.1 section 2.13.3);
    ``solvation_model="CPCM"`` switches the continuum model while the
    solvent name stays the same.  The thermochemistry block gives the
    'Final Gibbs free energy' line that
    :func:`delfin.energies.find_gibbs_energy` reads.
    """
    keywords = [
        functional,
        basis,
        "Opt",
        "Freq",
        f"{solvation_model}({solvent})",
    ]
    if extra_keywords:
        keywords.extend(extra_keywords)
    return _render_input(
        xyz_text=xyz_text,
        charge=charge,
        multiplicity=multiplicity,
        keywords=keywords,
        pal=pal,
        maxcore=maxcore,
    )


def build_single_point_input(
    xyz_text: str,
    charge: int,
    multiplicity: int,
    solvent: str,
    *,
    functional: str = _DEFAULT_FUNCTIONAL,
    basis: str = _DEFAULT_SP_BASIS,
    pal: int = 8,
    maxcore: int = 4000,
    solvation_model: str = "SMD",
    extra_keywords: Optional[Iterable[str]] = None,
) -> str:
    """Write the ORCA single-point input refining the energy at a larger basis.

    The geometry comes from the Opt+Freq step; this step refines the
    electronic energy only.  The printed Gibbs correction is computed at
    ``basis`` and must NOT be mixed into cycles whose Freq ran on a
    different basis -- combine refined energies with the small-basis
    thermochemistry correction explicitly instead.
    """
    keywords = [functional, basis, f"{solvation_model}({solvent})", "TightSCF"]
    if extra_keywords:
        keywords.extend(extra_keywords)
    return _render_input(
        xyz_text=xyz_text,
        charge=charge,
        multiplicity=multiplicity,
        keywords=keywords,
        pal=pal,
        maxcore=maxcore,
    )


def _render_input(
    *,
    xyz_text: str,
    charge: int,
    multiplicity: int,
    keywords: Iterable[str],
    pal: int,
    maxcore: int,
) -> str:
    lines = [
        "! " + " ".join(kw for kw in keywords if kw),
        f"%maxcore {maxcore}",
        "%pal nprocs " + str(int(pal)) + " end",
    ]
    lines.append(f"* xyz {charge} {multiplicity}")
    # XYZ body: skip the atom-count header and the comment line, keep
    # every coordinate line verbatim.
    body = [line.strip() for line in str(xyz_text).splitlines() if line.strip()]
    if body:
        try:
            atom_count = int(body[0])
        except ValueError:
            atom_count = None
        if atom_count is not None and len(body) >= atom_count + 2:
            body = body[2:]
    lines.extend(body)
    lines.append("*")
    return "\n".join(lines) + "\n"


def read_cycle_gibbs_energies(
    outputs: Mapping[str, Path],
) -> Dict[str, Optional[float]]:
    """Read the Gibbs energy of each cycle species from its ORCA output.

    Reuses :func:`delfin.energies.find_gibbs_energy`, the same reader the
    classic pipeline and the reporting collector use, so a pKa cycle and
    a DELFIN run never disagree about what an output file says.

    Args:
        outputs: mapping from species label to ORCA output file.

    Returns:
        Mapping from label to Gibbs energy in Hartree; a missing output
        file or an output without thermochemistry yields None rather
        than a guess.
    """
    values: Dict[str, Optional[float]] = {}
    for label, path in outputs.items():
        if path is None or not Path(path).is_file():
            logger.info(
                "pKa cycle: no output for %s yet ( %s ); Gibbs energy None",
                label,
                path,
            )
            values[label] = None
            continue
        values[label] = find_gibbs_energy(str(path))
        if values[label] is None:
            logger.info(
                "pKa cycle: output %s has no Gibbs thermochemistry", path
            )
    return values


def _write_xyz_input_from_file(xyz_path: Path, workdir: Path) -> str:
    """Read the coordinate block of a generated XYZ as input text."""
    lines = xyz_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"XYZ file is incomplete: {xyz_path}")
    atom_count = int(lines[0].strip())
    return "\n".join(lines[2 : atom_count + 2]) + "\n"


def prepare_cycle(
    species_list: List[Species],
    cycle_dir: Path,
    *,
    write_smiles_inputs=None,
) -> List[Species]:
    """Generate the 3D structure of every species of the cycle.

    Structure generation goes through DELFIN's SMILES converter (via
    :func:`delfin.tadf_xtb._write_smiles_inputs`); a species without a
    SMILES string is passed through unchanged so pre-made XYZ inputs
    remain usable.

    Args:
        species_list: the four Species of one cycle, in plan order.
        cycle_dir: working directory receiving one XYZ per species.
        write_smiles_inputs: test seam; defaults to the DELFIN helper.

    Returns:
        The species list in the same order (a caller writes the inputs).
    """
    cycle_dir = Path(cycle_dir)
    cycle_dir.mkdir(parents=True, exist_ok=True)
    if write_smiles_inputs is None:
        write_smiles_inputs = _write_smiles_inputs
    for species in species_list:
        if species.smiles:
            _start, xyz_path = write_smiles_inputs(
                species.smiles, species.label, cycle_dir
            )
            if not xyz_path.is_file():
                raise RuntimeError(
                    f"SMILES conversion produced no XYZ for {species.label}"
                )
    return species_list


def run_cycle(
    acid: str,
    cycle_dir: Path,
    *,
    solvent: str = _DEFAULT_SOLVENT,
    solvation_model: str = "SMD",
    pal: int = 8,
    maxcore: int = 4000,
    reference_acid: str = REFERENCE_ACID,
    temperature_k: float = PKA_KELVIN_DEFAULT,
) -> Dict[str, Any]:
    """Run the full isodesmic cycle of one acid through ORCA.

    Per species: XYZ from DELFIN's SMILES converter, one Opt+Freq ORCA
    run (continuum solvation, SMD by default, CPCM selectable via
    ``solvation_model``), Gibbs energy read with
    energies.find_gibbs_energy through read_cycle_gibbs_energies.  The
    pKa is anchored to the reference acid's experimental value with
    pka_from_cycle.

    Returns:
        Dict with keys:
        - ``pka``: the estimated pKa, or None when any species failed
          (the cycle refuses to guess from partial data);
        - ``gibbs``: label -> Gibbs energy in Hartree (None where absent);
        - ``failures``: labels whose ORCA run failed or produced no
          thermochemistry.
    """
    cycle_dir = Path(cycle_dir)
    species_list = plan_species(acid, reference_acid=reference_acid)
    prepare_cycle(species_list, cycle_dir)

    gibbs: Dict[str, Optional[float]] = {}
    failures: List[str] = []
    for species in species_list:
        xyz_path = cycle_dir / f"{species.label}.xyz"
        inp_path = cycle_dir / f"{species.label}.inp"
        out_path = cycle_dir / f"{species.label}.out"
        xyz_text = _write_xyz_input_from_file(xyz_path, cycle_dir)
        inp_path.write_text(
            build_opt_freq_input(
                xyz_text,
                charge=species.charge,
                multiplicity=species.multiplicity,
                solvent=solvent,
                solvation_model=solvation_model,
                pal=pal,
                maxcore=maxcore,
            ),
            encoding="utf-8",
        )
        # Call the module attribute, not a saved alias: tests patch
        # delfin.pka.run_orca, and an alias would silently bypass the
        # patch and hit the real ORCA resolver (module probe, process
        # budget).
        ok = run_orca(str(inp_path), str(out_path))
        if not ok:
            logger.error(
                "pKa cycle: ORCA run for %s failed ( %s )", species.label, out_path
            )
            failures.append(species.label)
            gibbs[species.label] = None
            continue
        gibbs[species.label] = find_gibbs_energy(str(out_path))
        if gibbs[species.label] is None:
            failures.append(species.label)

    # Map per-species Gibbs energies onto the cycle keys.
    target = acid
    reference = reference_acid
    cycle_gibbs = {
        "target_HA": gibbs.get(f"{target}_{_HA}"),
        "target_A": gibbs.get(f"{target}_{_A}"),
        "reference_HA": gibbs.get(f"{reference}_{_HA}"),
        "reference_A": gibbs.get(f"{reference}_{_A}"),
    }
    if failures:
        logger.warning(
            "pKa cycle for %s incomplete: failed species: %s",
            acid,
            ", ".join(failures),
        )
        return {"pka": None, "gibbs": cycle_gibbs, "failures": failures}

    try:
        pka = pka_from_cycle(
            cycle_gibbs,
            reference_pka=KNOWN_ACIDS[reference]["pka"],
            temperature_k=temperature_k,
        )
    except ValueError as exc:
        logger.warning("pKa cycle for %s incomplete: %s", acid, exc)
        return {"pka": None, "gibbs": cycle_gibbs, "failures": failures}

    return {"pka": pka, "gibbs": cycle_gibbs, "failures": []}


