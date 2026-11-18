"""Fragment reference energies for the proton coordinator (HER competition).

After the relaxed M...H scan the reaction energy needs a fragment
reference: the scan endpoint (combined complex+donor frame at the
stretched distance) minus the total energy of the separated products
[M-H] and the donor product (conjugate base / radical), each relaxed at
the same level of theory. Only the electronic + free-energy bookkeeping
is implemented here; the actual ORCA runs are driven by the caller
(CLI `delfin proton` pipeline or manual ORCA invocation).
"""
from __future__ import annotations

import json
import os
from typing import Dict, Optional

from ase import Atoms
from ase.io import read as ase_read, write as ase_write


def run_fragment_reference(
    *,
    combined_path: str,
    complex_n_atoms: int,
    donor_h_global: int,
    bookkeeping: Dict,
    keywords: str,
    PAL: int,
    maxcore: int,
    outdir: str = "fragment_reference",
) -> Optional[Dict]:
    """End-to-end fragment reference: split, write inputs, parse energies,
    build DeltaG. Returns the reactions dict (also written to
    reactions.json) or None if the ORCA outputs are not there yet.

    The caller (proton_coordinator) knows donor_h_global from
    combine_complex_donor; this module keeps all energy parsing and
    reaction-energy arithmetic.
    """
    combined = ase_read(combined_path)
    bk = bookkeeping
    os.makedirs(outdir, exist_ok=True)

    # --- split ------------------------------------------------------------
    hydride_atoms = Atoms()   # [M-H]: complex + transferred H
    donor_atoms = Atoms()     # donor product: donor minus transferred H
    for i, atom in enumerate(combined):
        if i < complex_n_atoms or i == donor_h_global:
            hydride_atoms += atom
        else:
            donor_atoms += atom

    cspec = bk["complex_product"]
    dspec = bk["donor_product"]

    hydride_inp = os.path.join(outdir, f"{cspec['name']}.inp")
    donor_inp = os.path.join(outdir, f"{dspec['name']}.inp")
    _write_opt_input(hydride_atoms, hydride_inp, cspec, keywords, PAL,
                     maxcore)
    _write_opt_input(donor_atoms, donor_inp, dspec, keywords, PAL, maxcore)
    ase_write(os.path.join(outdir, f"{cspec['name']}.xyz"), hydride_atoms)
    ase_write(os.path.join(outdir, f"{dspec['name']}.xyz"), donor_atoms)

    # --- energies (if the OPTs have run) ------------------------------------
    from delfin.proton.proton_coordinator import parse_orca_energy
    e_hydride = parse_orca_energy(os.path.join(outdir, f"{cspec['name']}.out"))
    e_donor = parse_orca_energy(os.path.join(outdir, f"{dspec['name']}.out"))
    e_scan_end = parse_orca_energy(
        os.path.join("relaxed_surface_scan", "scan.relaxscanact.out"))

    reaction = {
        "reaction": bk["fragment_reaction"],
        "scheme": bk["scheme"],
        "combined_charge": bk["combined_charge"],
        "combined_multiplicity": bk["combined_multiplicity"],
        "complex_product": cspec,
        "donor_product": dspec,
        "energies_Eh": {
            "scan_end": e_scan_end,
            f"{cspec['name']}": e_hydride,
            f"{dspec['name']}": e_donor,
        },
    }

    if None in (e_hydride, e_donor, e_scan_end):
        reaction["status"] = "incomplete"
        reaction["note"] = (
            "Run the two fragment OPTs (and finish the scan), then rerun "
            "to evaluate the reaction energy. Missing: "
            + ", ".join(name for name, e in
                        [("scan_end", e_scan_end),
                         (cspec["name"], e_hydride),
                         (dspec["name"], e_donor)] if e is None))
    else:
        # DeltaE = E(products) - E(scan end). The scan end is the combined
        # frame at the stretched M...H distance; the products are the
        # relaxed fragments. A negative value means the separated,
        # relaxed products are LOWER in energy than the scan-end point
        # (proton transfer completed and favourable at this distance).
        dE = (e_hydride + e_donor) - e_scan_end
        reaction["status"] = "complete"
        reaction["deltaE_products_minus_scanend_Eh"] = dE
        reaction["deltaE_products_minus_scanend_kcalmol"] = dE * 627.509474

    with open(os.path.join(outdir, "reactions.json"), "w") as f:
        json.dump(reaction, f, indent=2)
    print(f"[INFO] fragment reference: status={reaction['status']}"
          + (f" dE={reaction.get('deltaE_products_minus_scanend_kcalmol'):.2f}"
             " kcal/mol" if reaction["status"] == "complete" else ""))
    return reaction


def _write_opt_input(atoms: Atoms, path: str, spec: Dict, keywords: str,
                     PAL: int, maxcore: int) -> None:
    lines = [f"! {keywords}", "", "%pal nprocs {0} end".format(PAL),
             "%maxcore {0}".format(maxcore), ""]
    if len(atoms) == 0:
        raise ValueError(f"Empty fragment for {spec['name']} - "
                         "donor_h_global outside the donor range?")
    symbols = atoms.get_chemical_symbols()
    coords = atoms.positions
    lines.append(f"* xyz {spec['charge']} {spec['multiplicity']}")
    for s, c in zip(symbols, coords):
        lines.append(f"  {s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}")
    lines.append("*")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
