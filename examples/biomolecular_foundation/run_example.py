from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from delfin.tools import Pipeline


def run_example(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("The example requires an empty output directory")
    fixtures = Path(__file__).resolve().parent
    inputs = {name: fixtures / name for name in ("receptor.pdb", "reference.sdf", "poses.sdf")}
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()}
    pipeline = Pipeline("biomolecular_foundation_example")
    pipeline.add("biostructure_summary", structure_file=str(inputs["receptor.pdb"]))
    pipeline.add("ligand_prepare", ligand_file=str(inputs["reference.sdf"]))
    pipeline.add("pose_rmsd", pose_file=str(inputs["poses.sdf"]), score_tag="DEMO_SCORE")
    result = pipeline.run(work_dir=output_dir, cores=1)
    if not result.ok:
        errors = [step.error for step in result.all_results if not step.ok]
        raise RuntimeError("; ".join(errors) or "The example pipeline did not complete")
    steps = result.all_results
    report = {
        "kind": "synthetic_offline_io_example",
        "live_docking_performed": False,
        "score_meaning": "arbitrary demonstration units, not an affinity or energy",
        "input_sha256": hashes,
        "steps": [
            {"name": step.step_name, "data": step.data,
             "artifacts": {key: str(Path(path).relative_to(output_dir))
                           for key, path in step.artifacts.items()}}
            for step in steps
        ],
        "expected_rmsd_angstrom": [0.0, 1.0, 3.0],
        "observed_n_under_cutoff": steps[-1].data["n_under_cutoff"],
    }
    (output_dir / "example_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the offline biomolecular I/O example")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run_example(args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
