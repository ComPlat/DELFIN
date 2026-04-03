# Job Prioritization in DELFIN

DELFIN gives higher priority to bottleneck jobs, meaning jobs that unblock many
downstream jobs. By default, a job with at least 3 downstream dependents is
promoted to `HIGH` priority.

DELFIN also detects exclusive bottlenecks. If all other pending jobs depend on
one ready job, the scheduler can wait for currently running work to finish and
then launch that bottleneck with the maximum practical core allocation.

This improves throughput without changing scientific results. It only affects
job ordering and resource assignment.

Relevant implementation points:

- `delfin.job_priority.count_downstream_jobs(...)`
- `delfin.job_priority.is_exclusive_bottleneck(...)`
- `delfin.job_priority.adjust_job_priorities(...)`
- `parallel_classic_manually._adjust_priorities_for_bottlenecks(...)`
- `parallel_classic_manually._check_exclusive_bottleneck_boost(...)`

Typical log messages:

```text
[priority] Job occ_proc_initial boosted to HIGH priority (5 downstream jobs depend on it)
[occupier] Waiting for running jobs to finish before starting exclusive bottleneck occ_proc_initial
[occupier] Job occ_proc_initial is exclusive bottleneck -> allocating max cores (64)
```

This is a heuristic optimization. The threshold and waiting behavior may need
future tuning for different workflow shapes.
