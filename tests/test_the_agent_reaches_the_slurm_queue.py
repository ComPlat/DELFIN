"""On a SLURM host the agent's queue tools get a SLURM backend.

api._resolve_backend imported ``SLURMBackend``, a name backend_slurm never
defined -- the class is ``SlurmJobBackend``. On every cluster node, each agent
call that needed the queue (list_active_calculations, submit, cancel, recalc)
answered with "cannot import name 'SLURMBackend'". The local branch had been
fixed for the same kind of rename; this one was only reachable where sbatch
exists, which the test machines of that fix were not.
"""

from __future__ import annotations

import shutil

from delfin import api
from delfin.dashboard import backend_slurm


def test_a_slurm_host_gets_the_backend_the_dashboard_builds(monkeypatch):
    built = {}

    class _Backend:
        def __init__(self, **kwargs):
            built.update(kwargs)

    monkeypatch.setattr(shutil, "which",
                        lambda name: f"/usr/bin/{name}"
                        if name in ("sbatch", "squeue") else None)
    monkeypatch.setattr(backend_slurm, "SlurmJobBackend", _Backend)
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"runtime": {"slurm": {"profile": "bwunicluster3"}}})

    backend = api._resolve_backend()

    assert isinstance(backend, _Backend)
    assert built["slurm_profile"] == "bwunicluster3"
    assert built["submit_templates_dir"].name == "submit_templates"
