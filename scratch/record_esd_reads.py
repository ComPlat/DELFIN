#!/usr/bin/env python3
"""Baseline / post-switch recorder for the ESD reading chain.

Records, for every run directory given on the command line, the exact
outputs of every function the ESD chain uses to read ORCA output:

  - delfin.esd_results.collect_esd_results  (states fspe/zpe/gibbs,
    ISC/IC/fluor/phosp rate records) -- the workflow-level result
  - delfin.esd_input_generator.calculate_dele_cm1
  - the underlying readers find_electronic_energy / find_ZPE /
    find_gibbs_energy per state output (direct record)
  - delfin.uv_vis_spectrum.parse_absorption_spectrum (not switching --
    recorded to prove the workflow output is unchanged)
  - delfin.esd_input_generator._parse_tddft_and_derive_deltascf

States/ISCs/ICs are derived from the files actually present in the run
directory, so the recorder measures the real content of each run.

Floats are recorded with repr() (exact double).  Writes JSON to the path
given as last argument.  Runs from the worktree root (PYTHONPATH=.) so
``delfin`` resolves to THIS tree.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def record_run(esd_dir: Path, run_id: str) -> dict:
    from delfin.esd_results import collect_esd_results
    from delfin.esd_input_generator import (
        calculate_dele_cm1, _parse_tddft_and_derive_deltascf,
        _resolve_state_filename,
    )
    from delfin.energies import (
        find_electronic_energy, find_ZPE, find_gibbs_energy,
    )
    from delfin.uv_vis_spectrum import parse_absorption_spectrum

    rec: dict = {"run_id": run_id, "dir": str(esd_dir.name)}

    # discover state outputs (S0.out, S1.out, T1.out, ... and hybrid1
    # "<state>.out" resolved via _resolve_state_filename for tddft mode)
    states: list[str] = []
    for f in sorted(esd_dir.glob("*.out")):
        m = re.match(r"([ST]\d+)(?:\+.*)?\.out$", f.name)
        if m:
            s = m.group(1)
            if s not in states:
                states.append(s)
    # hybrid1 names resolve to "<state>.out" in tddft mode; try '+' variants
    for f in sorted(esd_dir.glob("*.out")):
        m = re.match(r"([ST]\d+)\+\.out$", f.name)
        if m and m.group(1) not in states:
            states.append(m.group(1))

    # discover ISC / IC transition files: <A>_<B>_ISC_<ms>.out, <A>_<B>_IC.out
    iscs: list[str] = []
    ics: list[str] = []
    seen_isc: set[tuple] = set()
    seen_ic: set[tuple] = set()
    for f in sorted(esd_dir.glob("*_ISC_*.out")):
        m = re.match(r"([ST]\d+)_([ST]\d+)_ISC_.+\.out$", f.name)
        if m and (m.group(1), m.group(2)) not in seen_isc:
            seen_isc.add((m.group(1), m.group(2)))
            iscs.append(f"{m.group(1)}>{m.group(2)}")
    for f in sorted(esd_dir.glob("*_IC.out")):
        m = re.match(r"([ST]\d+)_([ST]\d+)_IC\.out$", f.name)
        if m and (m.group(1), m.group(2)) not in seen_ic:
            seen_ic.add((m.group(1), m.group(2)))
            ics.append(f"{m.group(1)}>{m.group(2)}")

    rec["states_discovered"] = states
    rec["iscs_discovered"] = iscs
    rec["ics_discovered"] = ics

    # workflow-level collect
    try:
        summary = collect_esd_results(esd_dir, states, iscs, ics, config=None)
        wf: dict = {"states": {}, "isc": {}, "ic": {}, "fluor": {}, "phosp": {}}
        for state, sr in summary.states.items():
            wf["states"][state] = {
                "fspe": repr(sr.fspe), "zpe": repr(sr.zpe),
                "gibbs": repr(sr.gibbs),
            }
        for key, ir in summary.isc.items():
            wf["isc"][key] = {
                "rate": repr(ir.rate),
                "temperature": repr(ir.temperature),
                "delta_cm1": repr(ir.delta_cm1),
                "soc": (None if ir.soc is None
                        else [repr(ir.soc[0]), repr(ir.soc[1])]),
                "fc_percent": repr(ir.fc_percent),
                "ht_percent": repr(ir.ht_percent),
                "problem": ir.problem,
            }
        for key, cr in summary.ic.items():
            wf["ic"][key] = {"rate": repr(cr.rate), "problem": cr.problem}
        for key, fr in summary.fluor.items():
            wf["fluor"][key] = {"rate": repr(fr.rate), "problem": fr.problem}
        for key, pr in summary.phosp.items():
            wf["phosp"][key] = {
                "rate_mean": repr(pr.rate_mean),
                "sublevel_rates": [repr(v) for v in (pr.sublevel_rates or ())],
                "problem": pr.problem,
            }
        rec["collect_esd_results"] = wf
    except Exception as exc:
        rec["collect_esd_results"] = {"error": f"{type(exc).__name__}: {exc}"}

    # direct readers per state output (the exact functions called inside)
    direct: dict = {}
    for state in states:
        fname = _resolve_state_filename(state, "out", "tddft")
        p = esd_dir / fname
        if p.exists():
            direct[state] = {
                "fspe": repr(find_electronic_energy(str(p))),
                "zpe": repr(find_ZPE(str(p))),
                "gibbs": repr(find_gibbs_energy(str(p))),
            }
    rec["direct_readers"] = direct

    # DELE where both S0 and T1 exist
    if (esd_dir / "S0.out").exists() and (esd_dir / "T1.out").exists():
        try:
            dele = calculate_dele_cm1(str(esd_dir / "S0.out"),
                                      str(esd_dir / "T1.out"))
            rec["dele_cm1"] = repr(dele)
        except Exception as exc:
            rec["dele_cm1"] = f"error {type(exc).__name__}: {exc}"

    # deltaSCF derivation from S0.out (transition assignment reading)
    s0_out = esd_dir / "S0.out"
    if s0_out.exists():
        for state in ("S1", "T1"):
            try:
                got = _parse_tddft_and_derive_deltascf(s0_out, state)
                rec[f"deltascf_{state}"] = (None if got is None
                                            else [got[0], got[1]])
            except Exception as exc:
                rec[f"deltascf_{state}"] = (f"error {type(exc).__name__}: {exc}")

    # absorption spectrum (not switching; workflow-level proof)
    spectrum_files = {f for f in sorted(esd_dir.glob("*.out"))
                      if "TDDFT" in f.name}
    rec.setdefault("uv_vis", {})
    for spectrum_file in sorted(spectrum_files):
        try:
            trans = parse_absorption_spectrum(spectrum_file)
            rec["uv_vis"][spectrum_file.name] = [
                [t.from_state, t.to_state, repr(t.energy_ev),
                 repr(t.energy_cm1), repr(t.wavelength_nm), repr(t.fosc),
                 repr(t.d2), repr(t.dx), repr(t.dy), repr(t.dz)]
                for t in trans
            ]
        except Exception as exc:
            rec["uv_vis"][spectrum_file.name] = (
                f"error {type(exc).__name__}: {exc}")
    return rec


def main() -> None:
    out_path = Path(sys.argv[-1])
    run_dirs = [Path(a) for a in sys.argv[1:-1]]
    records = []
    for d in run_dirs:
        run_id = str(d)
        records.append(record_run(d, run_id))
    out_path.write_text(json.dumps(records, indent=1, sort_keys=True))
    print(f"recorded {len(records)} runs -> {out_path}")


if __name__ == "__main__":
    main()
