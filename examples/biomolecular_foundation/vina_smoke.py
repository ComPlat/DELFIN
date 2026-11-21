"""Explicit installed-engine I/O probe for a compute job, not a benchmark.

Not part of pytest: ordinary test collection must never discover and launch
whatever scientific engine happens to be installed on the host.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from delfin.tools import run_step
from delfin.tools.bio import BindingBox


def run_smoke(output: Path):
    binary = shutil.which('vina')
    if binary is None:
        raise RuntimeError('Vina must already be installed; this probe does not install it')
    version = subprocess.run([binary, '--version'], capture_output=True, text=True, timeout=10, check=True)
    banner = (version.stdout + version.stderr).strip()
    output.mkdir(parents=True, exist_ok=False)
    receptor = output / 'receptor.pdbqt'
    receptor.write_text('ATOM      1  C   ALA A   1       4.000   0.000   0.000  1.00  0.00     0.000 C\n')
    ligand = output / 'ligand.pdbqt'
    ligand.write_text('ROOT\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n')
    site = output / 'box.json'
    site.write_text(json.dumps(BindingBox(center=(0, 0, 0), size=(10, 10, 10)).to_dict()))
    result = run_step('vina_dock', receptor_file=receptor, ligand_file=ligand, site_file=site,
                      expected_version=banner, seed=7, exhaustiveness=1, num_modes=1, timeout=30,
                      work_dir=output / 'dock', cores=1)
    if not result.ok:
        raise RuntimeError(result.error)
    if len(result.data['scores']) != 1 or result.data['scientific_validation'] != 'not_performed':
        raise RuntimeError('unexpected Vina smoke result')
    return {'status': 'success', 'output': str(output), 'scientific_validation': 'not_performed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='new directory for synthetic inputs and engine output')
    args = parser.parse_args()
    print(json.dumps(run_smoke(args.output.resolve()), indent=2))


if __name__ == '__main__':
    main()
