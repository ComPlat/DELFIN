import os
import tempfile
from delfin.co2.CO2_Coordinator6 import write_default_files

def test_multiplicity_overrides_control_template(tmp_path=None):
    if tmp_path is None:
        tmp_path = tempfile.mkdtemp()
    control_path = os.path.join(tmp_path, 'CONTROL.txt')
    co2_path = os.path.join(tmp_path, 'co2.xyz')
    write_default_files(control_path=control_path, co2_path=co2_path, charge=0, multiplicity=3, solvent='DMF', overwrite=True)
    with open(control_path, 'r') as f:
        content = f.read()
    assert 'multiplicity=3' in content
    assert os.path.exists(co2_path)
