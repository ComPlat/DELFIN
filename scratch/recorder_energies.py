"""Probe: current readers vs orca_read on the four fixtures + archive sample.

Records, for each run and each candidate quantity, the value from the
current caller path and from delfin.orca_read.read(...), plus the source
so we can see when OPI is involved. This is the RECORDER for step 1 (baseline)
reused after the switch (step 3). Write-only into ../scratch/; never touches
the fixtures or the archive.

Usage: run from the DELFIN worktree root.
"""
import json
import sys
from pathlib import Path

from delfin import energies as E
from delfin.orca_read import Quantity, read as orca_read

WORKTREE = Path(__file__).resolve().parents[1]
FIXTURES = WORKTREE / "tests" / "fixtures" / "orca611"


def read_via_energies(out_path: Path):
    """The current caller reading path (all 5 quantities where defined)."""
    s = str(out_path)
    rec = {
        "final": E.find_electronic_energy(s),
        "gibbs": E.find_gibbs_energy(s),
        "zpe": E.find_ZPE(s),
    }
    # enthalpy / inner energy are read by api.extract_thermochem / the
    # collector; record them via the same text functions orca_read uses.
    try:
        from delfin.api import extract_thermochem
        thermo = extract_thermochem(str(out_path.parent))
        rec["enthalpy"] = thermo.enthalpy_corr
    except Exception as exc:
        rec["enthalpy"] = {"err": f"{type(exc).__name__}: {exc}"}
    try:
        from delfin.orca_read import _text_inner_energy
        rec["inner"] = _text_inner_energy(out_path)
    except Exception as exc:
        rec["inner"] = {"err": f"{type(exc).__name__}: {exc}"}
    return rec


Q_BY_NAME = {
    "final": Quantity.FINAL_SINGLE_POINT_ENERGY,
    "gibbs": Quantity.GIBBS,
    "zpe": Quantity.ZPE,
    "enthalpy": Quantity.ENTHALPY,
    "inner": Quantity.INNER_ENERGY,
}


def read_via_orca(out_path: Path):
    rec = {}
    for name, q in Q_BY_NAME.items():
        r = orca_read(out_path.parent, out_path.stem, q)
        rec[name] = r.value
        rec[f"{name}_src"] = r.source
    return rec


def collect(dir_: Path):
    results = {}
    if not dir_.is_dir():
        return results
    run_dirs = [d for d in dir_.iterdir()]
    # every run is DIRECTLY one level down: <basename>/<basename>.out
    for d in sorted(run_dirs):
        if not d.is_dir():
            continue
        out = d / f"{d.name}.out"
        if not out.is_file():
            continue
        cur = read_via_energies(out)
        new = read_via_orca(out)
        results[d.name] = {
            "status": "ok",
            "inp_exists": (d / f"{d.name}.inp").is_file(),
            "out_size": out.stat().st_size,
            "current": cur,
            "orca": new,
        }
    return results


def main(out_json: str, roots):
    all_results = {"fixtures": collect(FIXTURES)}
    for label, root in roots:
        all_results[label] = collect(Path(root))
    out_path = Path(out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(all_results, indent=1, sort_keys=True))
    print(f"wrote {out_path} with {sum(len(v) for v in all_results.values())} runs")


if __name__ == "__main__":
    out_json = sys.argv[1]
    roots = []
    for i in range(2, len(sys.argv), 2):
        roots.append((sys.argv[i], sys.argv[i + 1]))
    main(out_json, roots)
