"""The ExtOpt job script is a shell script: no ORCA input block in it.

render_extopt_job_script wrote ``%pal / nprocs N / end`` into the bash
script for the delfin-mlp-extopt provider.  bash reads ``%pal`` as a
job specification; measured with bash under ``set -euo pipefail``: "fg:
no job control", exit status 1, the job stops before ORCA starts.  %pal belongs in the ORCA
input, never in the script.
"""

from __future__ import annotations

import re

from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.mlp_tools.extopt_submit import ExtOptJobRequest, plan_extopt_job

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


def test_delfin_provider_script_has_no_pal_block(tmp_path):
    req = ExtOptJobRequest(run_type="opt", elements=WATER[0], coordinates=WATER[1],
                           provider="delfin", backend="ani2x", work_dir=tmp_path)
    plan = plan_extopt_job(req, cfg=MLPClusterConfig(gpu_partitions="gpu", cpus_per_task=8))
    script = plan.script.read_text()
    assert not re.search(r"(?m)^\s*(%pal|nprocs\b|end\s*$)", script)
