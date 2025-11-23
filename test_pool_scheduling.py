#!/usr/bin/env python3
"""Test script to debug pool scheduling behavior."""

import time
import logging
from delfin.dynamic_pool import DynamicCorePool, PoolJob, JobPriority

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)

def make_job(job_id: str, cores_min: int, cores_opt: int, cores_max: int):
    """Create a test job that just sleeps."""
    def work(cores):
        print(f"[{job_id}] Starting with {cores} cores")
        time.sleep(5)
        print(f"[{job_id}] Completed")

    return PoolJob(
        job_id=job_id,
        cores_min=cores_min,
        cores_optimal=cores_opt,
        cores_max=cores_max,
        memory_mb=cores_opt * 1000,
        priority=JobPriority.NORMAL,
        execute_func=work,
        args=(),
        kwargs={},
        estimated_duration=5.0,
    )

# Create pool with 64 cores, max 4 jobs
pool = DynamicCorePool(total_cores=64, total_memory_mb=64000, max_jobs=4)
print(f"\nPool created: total_cores=64, max_concurrent_jobs={pool.max_concurrent_jobs}")

# Submit 4 jobs similar to OCCUPIER FoBs
jobs = [
    make_job("fob_1", cores_min=10, cores_opt=10, cores_max=16),
    make_job("fob_4", cores_min=10, cores_opt=10, cores_max=16),
    make_job("fob_6", cores_min=10, cores_opt=10, cores_max=16),
    make_job("fob_7", cores_min=10, cores_opt=10, cores_max=16),
]

print("\nSubmitting 4 jobs...")
for job in jobs:
    pool.submit_job(job)
    print(f"Submitted {job.job_id}")

# Wait a bit for scheduler to process
print("\nWaiting for jobs to start...")
time.sleep(2)

# Check status
status = pool.get_status()
print(f"\nPool status after 2s:")
print(f"  Running jobs: {status['running_jobs']}")
print(f"  Queued jobs: {status['queued_jobs']}")
print(f"  Allocated cores: {status['allocated_cores']}/{status['total_cores']}")
print(f"  Job details:")
for job_id, details in status['job_details'].items():
    print(f"    {job_id}: {details['cores']} cores, {details['duration']:.1f}s")

# Wait for completion
print("\nWaiting for all jobs to complete...")
pool.wait_for_completion()

# Final status
status = pool.get_status()
print(f"\nFinal status:")
print(f"  Completed jobs: {status['completed_jobs']}")
print(f"  Failed jobs: {status['failed_jobs']}")

pool.shutdown()
print("\nTest complete!")
