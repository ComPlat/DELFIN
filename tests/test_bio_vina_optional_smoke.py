"""Tiny synthetic I/O smoke test; never a docking benchmark."""
import json
import shutil
import subprocess

import pytest

from delfin.tools import run_step
from delfin.tools.bio import BindingBox


@pytest.mark.skipif(shutil.which('vina') is None, reason='optional Vina executable not installed')
def test_installed_vina_can_produce_a_checked_artifact(tmp_path):
    version = subprocess.run([shutil.which('vina'), '--version'], capture_output=True, text=True, timeout=10)
    banner = (version.stdout + version.stderr).strip()
    receptor = tmp_path / 'receptor.pdbqt'
    receptor.write_text('ATOM      1  C   ALA A   1       4.000   0.000   0.000  1.00  0.00     0.000 C\n')
    ligand = tmp_path / 'ligand.pdbqt'
    ligand.write_text('ROOT\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n')
    site = tmp_path / 'box.json'
    site.write_text(json.dumps(BindingBox(center=(0, 0, 0), size=(10, 10, 10)).to_dict()))
    result = run_step('vina_dock', receptor_file=receptor, ligand_file=ligand, site_file=site,
                      expected_version=banner, seed=7, exhaustiveness=1, num_modes=1, timeout=30,
                      work_dir=tmp_path / 'dock', cores=1)
    assert result.ok, result.error
    assert len(result.data['scores']) == 1
    assert result.data['scientific_validation'] == 'not_performed'
