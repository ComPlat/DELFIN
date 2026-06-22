import time

import pytest

pytest.importorskip("rdkit")

from delfin.config import _load_template_defaults
from delfin.pipeline import PipelineContext
from delfin.stability_constant import THERMODYNAMICS_DIRNAME, build_stability_reaction_plan


def test_reaction_plan_uses_thermodynamics_output_dir(tmp_path):
    config = _load_template_defaults()
    config["SMILES"] = "O"
    config["solvent"] = "water"
    config["stability_constant"] = "yes"
    config["stability_constant_mode"] = "reaction"
    config["stability_reaction"] = "1*{N}>>>1*{N}"
    config["thdy_smiles_converter"] = "NORMAL"
    config["thdy_preopt"] = "xtb"
    config["PAL"] = "1"

    ctx = PipelineContext(
        config=config,
        control_file_path=tmp_path / "CONTROL.txt",
        input_file="input.txt",
        charge=0,
        PAL=1,
        multiplicity=1,
        solvent="water",
        metals=[],
        main_basisset="def2-SVP",
        metal_basisset="def2-SVP",
        number_explicit_solv_molecules=0,
        total_electrons_txt=10,
        start_time=time.time(),
        name="test",
    )

    plan = build_stability_reaction_plan(ctx)

    assert plan.sc_dir == tmp_path / THERMODYNAMICS_DIRNAME
    assert plan.sc_dir.is_dir()
    assert all(THERMODYNAMICS_DIRNAME in str(job.working_dir) for job in plan.jobs if job.working_dir is not None)
