"""A cluster's Slurm, asked over its tunnel.

``delfin.scheduler_client`` drives the Slurm of the machine DELFIN runs on;
this module asks a cluster's, through a ``delfin.remote.connection``
Connection, with argv lists. Its callers throttle it -- the window's
details (``clusters_feed``) reuse a job's facts for a few seconds -- which
is why tests/test_r3_direct_scheduler_call_guard.py lists it among the
sanctioned scheduler callers.
"""

from __future__ import annotations

import re

#: sacct's columns when scontrol no longer knows a job (it forgets a job a
#: few minutes after it ended), named as scontrol names them.
_SACCT = (("JobName", "JobName"), ("State", "JobState"), ("Reason", "Reason"),
          ("Partition", "Partition"), ("NNodes", "NumNodes"), ("NCPUS", "NumCPUs"),
          ("AllocTRES", "AllocTRES"), ("Submit", "SubmitTime"), ("Start", "StartTime"),
          ("Elapsed", "RunTime"), ("Timelimit", "TimeLimit"), ("WorkDir", "WorkDir"),
          ("Comment", "Comment"))


def job_fields(conn, job_id: str) -> dict[str, str]:
    """What Slurm says of ``job_id``, keyed as ``scontrol show job`` names it
    (sacct's answer, renamed, once scontrol has forgotten the job); empty
    when neither knows it."""
    shown = conn.run(["scontrol", "show", "job", "-o", job_id], timeout=30)
    fields = parse_scontrol(shown.stdout or "") if shown.ok else {}
    if fields:
        return fields
    columns = [column for column, _ in _SACCT]
    got = conn.run(["sacct", "-j", job_id, "-X", "-n", "-P", "-o", ",".join(columns)],
                   timeout=30)
    row = (got.stdout or "").strip().splitlines()[:1] if got.ok else []
    if not row:
        return {}
    values = dict(zip(columns, row[0].split("|")))
    return {name: values.get(column, "") for column, name in _SACCT}


def parse_scontrol(text: str) -> dict[str, str]:
    """``scontrol show job -o``'s one line as key -> value (a value may hold
    spaces: it runs to the next ``Key=``)."""
    line = text.strip().splitlines()[0] if text.strip() else ""
    out = {}
    for part in re.split(r" (?=[A-Za-z][A-Za-z/:]*=)", line):
        key, sep, value = part.partition("=")
        if sep:
            out[key] = value
    return out
