"""ESD-style excited state pipeline using reactive step population.

Recreates the ESD (Excited State Dynamics) pattern where the set of
calculations depends on which states are requested and which mode
(TDDFT vs deltaSCF) is used.

Usage::

    from delfin.tools.templates import esd_states

    # TDDFT mode: all excited states independent (parallel)
    pipe = esd_states(
        charge=0, mult=1,
        method="B3LYP", basis="def2-SVP",
        states=["S0", "S1", "S2", "T1"],
        mode="tddft",
    )
    result = pipe.run(geometry="optimized.xyz", cores=16)

    # deltaSCF mode: sequential dependency chain
    pipe = esd_states(
        charge=0, mult=1,
        method="B3LYP", basis="def2-TZVP",
        states=["S0", "S1", "T1"],
        mode="deltascf",
    )
    result = pipe.run(geometry="optimized.xyz", cores=8)
"""

from __future__ import annotations

from typing import List, Optional

from delfin.tools.pipeline import Pipeline


def esd_states(
    *,
    charge: int,
    mult: int = 1,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    states: Optional[List[str]] = None,
    mode: str = "tddft",
    include_freq: bool = False,
    include_isc: bool = False,
    **orca_kwargs,
) -> Pipeline:
    """Build an ESD-style excited state pipeline.

    Parameters
    ----------
    charge, mult : molecular charge and spin multiplicity
    method, basis : DFT method and basis set
    states : list of state labels to compute (default: ["S0", "S1", "T1"])
    mode : "tddft" (all states as branches, parallel) or
           "deltascf" (sequential dependency chain)
    include_freq : also run frequency calculations for each state
    include_isc : add inter-system crossing rate calculations
    **orca_kwargs : additional ORCA parameters

    Returns
    -------
    Pipeline ready to run with ``pipe.run(geometry=..., cores=...)``
    """
    if states is None:
        states = ["S0", "S1", "T1"]

    defaults = {
        "charge": charge, "mult": mult,
        "method": method, "basis": basis,
        **orca_kwargs,
    }

    pipe = Pipeline("esd", defaults=defaults)

    if mode == "tddft":
        # TDDFT mode: S0 as trunk, excited states as parallel branches
        pipe.add("orca_opt", label="S0_opt")
        if include_freq:
            pipe.add("orca_freq", label="S0_freq")

        for state in states:
            if state == "S0":
                continue
            branch = pipe.branch(state)
            # TDDFT excited state optimization
            root_idx = int(state[1:]) if state[1:].isdigit() else 1
            is_triplet = state.startswith("T")
            branch.add(
                "orca_tddft",
                nroots=max(root_idx, 1),
                iroot=root_idx,
                triplet=is_triplet,
                label=f"{state}_tddft",
            )
            if include_freq:
                branch.add(
                    "orca_freq",
                    label=f"{state}_freq",
                )

    elif mode == "deltascf":
        # deltaSCF mode: sequential chain S0 → S1 → S2 → ...
        # Each state depends on the previous one for orbital guess
        pipe.add("orca_opt", label="S0_opt")
        if include_freq:
            pipe.add("orca_freq", label="S0_freq")

        for state in states:
            if state == "S0":
                continue
            root_idx = int(state[1:]) if state[1:].isdigit() else 1
            is_triplet = state.startswith("T")
            pipe.add(
                "orca_opt",
                label=f"{state}_opt",
            )
            if include_freq:
                pipe.add("orca_freq", label=f"{state}_freq")

    else:
        raise ValueError(f"Unknown ESD mode: {mode!r} (use 'tddft' or 'deltascf')")

    # Optional: ISC rate calculations between singlet/triplet pairs
    if include_isc:
        singlets = [s for s in states if s.startswith("S") and s != "S0"]
        triplets = [s for s in states if s.startswith("T")]
        if singlets and triplets:
            def _add_isc(results, last, inserter):
                for s in singlets:
                    for t in triplets:
                        inserter.add(
                            "orca_sp",
                            label=f"isc_{s}_{t}",
                        )
            pipe.add_reactive(_add_isc, label="isc_rates")

    return pipe
