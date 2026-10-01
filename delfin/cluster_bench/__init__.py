"""Batch construction runs on a Slurm cluster: MANTA and, for comparison, other builders.

``delfin cluster`` (also ``delfin-cluster`` and ``python -m delfin.cluster_bench``) cuts a list
of ``ID;SMILES`` lines into shards, writes a Slurm array script, builds one shard per array task,
reports progress and merges the chunks into one archive in DELFIN's multi-frame xyz format:

    prepare    -> run directory with shards, per-shard specs and a manifest (sha256 of every file,
                  code commit, environment versions, settings)
    slurm      -> the sbatch array script (48 cores, 72 h, array throttle), optionally submitted
    run-shard  -> the worker of one array task; resumable
    status     -> progress per shard
    collect    -> one archive + per-system class (ok / empty / timeout / not_expressible / fail)
    repeat-stats -> byte identity of a subset built twice

Tools: ``manta`` (the shipped MANTA construction, one child process per system), and the
external builders ``architector``, ``molsimplify`` and ``mace`` (epic-MACE).  An external builder
runs in its own Python environment: its worker (``tool_workers/``) is started with that
environment's interpreter and imports nothing from DELFIN.

No input data ships with this package; lists, selections and results live in the run directory.
See docs/CLUSTER_BATCH_JUSTUS.md.
"""
