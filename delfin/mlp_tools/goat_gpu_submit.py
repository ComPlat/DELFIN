"""Where a DELFIN job whose CONTROL.txt asks for an MLIP GOAT is sent.

``global_opt_method=mlip:<backend>`` with ``global_opt_mlip_device=cuda``
(the default) runs GOAT inside the DELFIN job, on the node the job lands on,
with the MLIP on that node's GPU (delfin/xtb_crest.py refuses to start it
without one).  So the DELFIN job itself needs a GPU partition and
``--gres=gpu:1``; a CPU partition is no fallback -- the job would stop at
the GOAT step.  The job submission (``SlurmJobBackend.submit_delfin``) asks
:func:`mlip_goat_method` whether a folder's CONTROL.txt is such a run and
:func:`gpu_partitions_for` which configured GPU partitions can hold its time
limit; it never changes PAL, maxcore or the time limit to fit one.
"""

from __future__ import annotations

import re
from typing import Callable, Optional, Sequence, Tuple


def _control_values(control_text: str) -> dict:
    """``key=value`` lines of a CONTROL text, first occurrence per key."""
    values: dict = {}
    for raw in (control_text or "").splitlines():
        line = raw.strip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key and key not in values:
            values[key] = value.strip()
    return values


def mlip_goat_method(control_text: str) -> str:
    """``mlip:<backend>`` when the CONTROL runs GOAT with an MLIP on a GPU, else "".

    GOAT is on when ``global_optimizer=GOAT``, or, without that choice, the
    older ``XTB_GOAT=yes`` (the precedence of delfin/config.py).  An MLIP
    method with ``global_opt_mlip_device=cpu`` needs no GPU, and a value the
    validator refuses is left to the validator (the run reports it).
    """
    from delfin.common.control_validator import (
        _as_global_opt_method,
        _as_global_opt_mlip_device,
    )
    from delfin.config import _is_yes_token, _normalize_global_optimizer

    values = _control_values(control_text)
    try:
        method = _as_global_opt_method(values.get("global_opt_method"))
        device = _as_global_opt_mlip_device(values.get("global_opt_mlip_device"))
    except ValueError:
        return ""
    if not method.startswith("mlip:") or device != "cuda":
        return ""
    chosen = _normalize_global_optimizer(values.get("global_optimizer"))
    goat = chosen == "GOAT" if chosen else _is_yes_token(values.get("XTB_GOAT"))
    return method if goat else ""


#: A partition named for AMD GPUs (MI250, MI300A, ...).  orca-mlips' stacks
#: are CUDA builds of PyTorch: on such a node the job sees no /dev/nvidia*
#: and GOAT stops before ORCA starts (measured: an ethanol MLIP GOAT routed
#: to an MI300 partition failed in 30 s).  The partition's TRES carries no
#: GPU vendor and a user may not read the nodes' Gres type, so the name is
#: what there is to go by.
_NON_CUDA_PARTITION = re.compile(r"(?i)(?:^|[_\-])(?:mi\d{2,3}[a-z]?|amd|rocm|instinct)(?:[_\-]|$)")


def cuda_partitions(candidates: Sequence[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """``(kept, dropped)``: the candidates split into CUDA and AMD-named partitions."""
    kept, dropped = [], []
    for name in candidates:
        name = str(name).strip()
        if name:
            (dropped if _NON_CUDA_PARTITION.search(name) else kept).append(name)
    return tuple(kept), tuple(dropped)


def gpu_partitions_for(candidates: Sequence[str], time_limit: str, *,
                       ask: Optional[Callable] = None) -> Tuple[str, str]:
    """``(partitions, reason)``: the CUDA candidates whose MaxTime holds *time_limit*.

    ``partitions`` is a comma-separated ``--partition`` value, "" with the
    reason when there is no candidate at all or none can hold the time
    limit.  Partitions named for AMD GPUs are left out
    (:func:`cuda_partitions`).  One scontrol query for all of them
    (:func:`delfin.mlp_tools.extopt_submit.partitions_for_time`).
    """
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import partitions_for_time

    names, dropped = cuda_partitions(candidates)
    if not names and dropped:
        return "", (
            f"the only GPU partitions ({', '.join(dropped)}) have AMD GPUs; the "
            "orca-mlips stacks are CUDA builds and need an NVIDIA GPU"
        )
    if not names:
        return "", (
            "no GPU partition is configured (runtime.slurm.gpu_partitions, "
            "DELFIN_SLURM_GPU_PARTITIONS or DELFIN_MLP_CLUSTER_GPU_PARTITIONS) "
            "and none was found; GOAT with an MLIP on cuda has no CPU fallback"
        )
    cfg = MLPClusterConfig(gpu_partitions=",".join(names), time_limit=str(time_limit))
    return partitions_for_time(cfg, ask)
