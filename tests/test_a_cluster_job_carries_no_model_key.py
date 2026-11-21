# SPDX-License-Identifier: MIT
"""A submitted job's environment carries no model-provider credential.

``sbatch --export=ALL`` propagates the whole submission environment to the
compute node, and that flag is required here: passing explicit variables
with ``--export`` makes slurmd fetch the user's login environment at launch,
which times out on this site and holds the job. So the environment handed to
``sbatch`` IS what the job gets, and everything the job runs can read it.

The bash cage already removes these variables from every command it runs
(``api_client._scrubbed_bash_env``: "Nothing the agent legitimately runs
needs the model provider's own key"). This submission path built its
environment from ``os.environ`` directly and never had that removal, so a
key the operator exported in the shell that started DELFIN travelled to the
node with every job -- and supplying a key through the environment is one of
the two documented ways to configure one.

Both directions are asserted, because a removal with no measured limit is
how a security change costs function: the model and repository credentials
must be gone, and the credentials a job genuinely uses must still be there.
Measured 2026-11-21 over ``delfin/submit_templates/``,
``delfin/mlp_tools/`` and ``delfin/cluster_bench/``: those files reference
``HF_TOKEN``, ``HUGGING_FACE_HUB_TOKEN``, ``HF_TOKEN_PATH`` and
``DELFIN_RUN_TOKEN``, and none of the removed names.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.dashboard.backend_slurm import _CREDENTIALS_NO_JOB_NEEDS

_ROOT = Path(__file__).resolve().parents[1]

#: What a job legitimately reads from its environment. Removing one of these
#: would break a GPU MLIP job's weight download or the run's own reporting.
_A_JOB_NEEDS = ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", "HF_TOKEN_PATH",
                "DELFIN_RUN_TOKEN")

#: A DISTINCT stand-in per name. One shared value cannot tell "removed"
#: from "present under another name": the kept credentials carry it
#: legitimately, so an assertion that the value appears nowhere is false by
#: construction. Measured the hard way.
def _fake(name: str) -> str:
    return f"zzz-not-a-real-credential-{name.lower()}-0123456789"


def _submit_env(monkeypatch, tmp_path) -> dict:
    """The environment the backend would hand to sbatch, captured.

    Driven through ``_sbatch`` itself rather than rebuilt here: the removal
    has to be in the path that submits, not in a helper beside it.
    """
    from delfin.dashboard import backend_slurm as B

    seen: dict = {}

    class _Done:
        returncode = 0
        stdout = "Submitted batch job 1\n"
        stderr = ""

    calls: list = []

    def _fake_run(cmd, **kw):
        calls.append((list(cmd), dict(kw.get("env") or {})))
        return _Done()

    monkeypatch.setattr(B.subprocess, "run", _fake_run)
    for name in (*_CREDENTIALS_NO_JOB_NEEDS, *_A_JOB_NEEDS):
        monkeypatch.setenv(name, _fake(name))

    script = tmp_path / "run.sh"
    script.write_text("#!/bin/bash\n#SBATCH --partition=dev\necho hi\n")
    backend = B.SlurmJobBackend(str(tmp_path))
    backend._sbatch(tmp_path, [], "01:00:00", 1, 1024, "job", script,
                    partition="dev")
    # The sbatch call, not merely the last subprocess the backend ran:
    # it also asks the queue afterwards, and scoring the wrong call is how
    # a probe measures something nobody asked about.
    subs = [(c, e) for c, e in calls
            if c and c[0] == "sbatch" and "--test-only" not in c]
    assert len(subs) == 1, [c for c, _ in calls]
    seen["cmd"], seen["env"] = subs[0]
    return seen


def test_the_flag_that_makes_this_matter_is_still_there(monkeypatch, tmp_path):
    """--export=ALL is why the submission environment IS the job's.

    If this flag ever goes, the reasoning behind the removal changes and
    this file should be re-read rather than silently kept.
    """
    seen = _submit_env(monkeypatch, tmp_path)
    assert "--export=ALL" in seen["cmd"], seen["cmd"]


@pytest.mark.parametrize("name", sorted(_CREDENTIALS_NO_JOB_NEEDS))
def test_a_model_or_repository_credential_does_not_reach_the_node(
        name, monkeypatch, tmp_path):
    seen = _submit_env(monkeypatch, tmp_path)
    assert name not in seen["env"], (
        f"{name} would travel to the compute node with every job")
    assert _fake(name) not in "\n".join(seen["env"].values()), (
        f"{name}'s value is still in the environment under another name")


@pytest.mark.parametrize("name", _A_JOB_NEEDS)
def test_a_credential_a_job_uses_is_still_passed(name, monkeypatch, tmp_path):
    """The other direction: the removal must not be a blanket one."""
    seen = _submit_env(monkeypatch, tmp_path)
    assert seen["env"].get(name) == _fake(name), (
        f"{name} was removed, which breaks the job that reads it")


def test_the_kept_names_are_the_ones_the_job_files_actually_reference():
    """Why the list is names and not the cage's suffix markers.

    Re-derived from the job-side sources, so the exemption cannot drift
    away from what a job really uses: a name kept here must appear there,
    and a name removed must not.
    """
    text = ""
    for folder in ("delfin/submit_templates", "delfin/mlp_tools",
                   "delfin/cluster_bench"):
        for path in sorted((_ROOT / folder).rglob("*")):
            if path.is_file() and path.suffix in (".py", ".sh", ""):
                text += path.read_text(encoding="utf-8", errors="replace")
    referenced = set(re.findall(
        r"[A-Z][A-Z0-9_]*(?:API_KEY|TOKEN|SECRET|PASSWORD)[A-Z0-9_]*", text))
    for name in _A_JOB_NEEDS:
        assert name in referenced, (
            f"{name} is kept for the job but no job file references it")
    leaked = sorted(set(_CREDENTIALS_NO_JOB_NEEDS) & referenced)
    assert not leaked, (
        "removed from the job environment although a job file reads it: "
        + ", ".join(leaked))
