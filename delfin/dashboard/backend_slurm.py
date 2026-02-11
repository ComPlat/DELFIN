"""SLURM job backend: sbatch / squeue / scancel."""

import os
import subprocess
from pathlib import Path
from typing import List, Tuple

from .backend_base import JobBackend, JobInfo, SubmitResult
from .input_processing import parse_resource_settings, parse_inp_resources


class SlurmJobBackend(JobBackend):
    """SLURM cluster backend (sbatch/squeue/scancel)."""

    def __init__(self, submit_templates_dir, orca_base=None):
        self.submit_templates_dir = Path(submit_templates_dir)
        self.orca_base = orca_base

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _resolve_resources(self, job_dir, inp_file=None, pal=40, maxcore=6000):
        """Determine PAL and memory from CONTROL / .inp / defaults."""
        pal_used = None
        maxcore_used = None

        # 1) Try CONTROL.txt
        try:
            control_path = Path(job_dir) / 'CONTROL.txt'
            if control_path.exists():
                pal_used, maxcore_used = parse_resource_settings(control_path.read_text())
        except Exception:
            pal_used, maxcore_used = None, None

        # 2) Fallback to .inp file
        if pal_used is None or maxcore_used is None:
            try:
                inp_candidate = None
                if inp_file:
                    cand = Path(job_dir) / inp_file
                    if cand.exists():
                        inp_candidate = cand
                if inp_candidate is None:
                    inp_files = sorted(Path(job_dir).glob('*.inp'))
                    inp_candidate = inp_files[0] if inp_files else None
                if inp_candidate is not None and inp_candidate.exists():
                    pal_inp, maxcore_inp = parse_inp_resources(inp_candidate.read_text())
                    if pal_used is None:
                        pal_used = pal_inp
                    if maxcore_used is None:
                        maxcore_used = maxcore_inp
            except Exception:
                pass

        if pal_used is None:
            pal_used = int(pal)
        if maxcore_used is None:
            maxcore_used = int(maxcore)

        mem_used = int(pal_used) * int(maxcore_used)
        return pal_used, mem_used

    def _sbatch(self, job_dir, env_vars, time_limit, pal, mem_mb,
                job_name, submit_script):
        """Run sbatch with the given parameters."""
        return subprocess.run(
            [
                'sbatch',
                f'--export=ALL,{env_vars}',
                f'--time={time_limit}',
                '--ntasks=1',
                f'--cpus-per-task={pal}',
                f'--mem={mem_mb}M',
                f'--job-name={job_name}',
                str(submit_script),
            ],
            cwd=str(job_dir),
            capture_output=True,
            text=True,
        )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def submit_delfin(self, job_dir, job_name, mode='delfin',
                      time_limit='48:00:00', pal=40, maxcore=6000,
                      override=None, build_mult=None) -> SubmitResult:
        env_vars = f'DELFIN_MODE={mode},DELFIN_JOB_NAME={job_name}'
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        if override:
            env_vars += f',DELFIN_OVERRIDE={override}'
        if build_mult is not None:
            env_vars += f',BUILD_MULTIPLICITY={build_mult}'

        pal_used, mem_used = self._resolve_resources(job_dir, pal=pal, maxcore=maxcore)
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=40, maxcore=6000) -> SubmitResult:
        env_vars = f'DELFIN_MODE=orca,DELFIN_JOB_NAME={job_name}'
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        if inp_file:
            env_vars += f',DELFIN_INP_FILE={inp_file}'

        pal_used, mem_used = self._resolve_resources(
            job_dir, inp_file=inp_file, pal=pal, maxcore=maxcore,
        )
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def submit_turbomole(self, job_dir, job_name, module='ridft',
                         time_limit='48:00:00', nprocs=40, mem_per_cpu=6000,
                         para_arch='SMP') -> SubmitResult:
        mem_mb = nprocs * mem_per_cpu
        env_vars = (
            f'TM_JOB_NAME={job_name},'
            f'TM_MODULE={module},'
            f'TM_NPROCS={nprocs},'
            f'TM_PARA_ARCH={para_arch}'
        )
        result = self._sbatch(
            job_dir, env_vars, time_limit, nprocs, mem_mb,
            job_name, self.submit_templates_dir / 'submit_turbomole.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def list_jobs(self) -> List[JobInfo]:
        try:
            result = subprocess.run(
                ['squeue', '-u', os.environ.get('USER', ''),
                 '-o', '%.12i %.12P %.35j %.10u %.3t %.12M %.6D %R'],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            return []

        if result.returncode != 0:
            return []

        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return []

        jobs = []
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.split()
            if not parts:
                continue
            job_id = parts[0].strip()
            partition = parts[1].strip() if len(parts) > 1 else ''
            name = parts[2].strip() if len(parts) > 2 else 'unknown'
            user = parts[3].strip() if len(parts) > 3 else ''
            status = parts[4].strip() if len(parts) > 4 else ''
            time_used = parts[5].strip() if len(parts) > 5 else ''
            nodes = parts[6].strip() if len(parts) > 6 else ''
            reason = parts[7].strip() if len(parts) > 7 else ''

            jobs.append(JobInfo(
                job_id=job_id,
                name=name,
                status=status,
                extra={
                    'partition': partition,
                    'user': user,
                    'time_used': time_used,
                    'nodes': nodes,
                    'reason': reason,
                    'raw_line': line,
                },
            ))
        return jobs

    def cancel_job(self, job_id) -> Tuple[bool, str]:
        try:
            result = subprocess.run(
                ['scancel', str(job_id)],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                return True, f'Job {job_id} cancelled successfully.'
            else:
                return False, f'Error cancelling job {job_id}: {result.stderr or result.stdout}'
        except Exception as e:
            return False, f'Error cancelling job {job_id}: {e}'

    def get_pending_start_times(self):
        try:
            result = subprocess.run(
                ['squeue', '-u', os.environ.get('USER', ''), '--start',
                 '-o', '%.12i %.35j %.20S'],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            return []

        if result.returncode != 0:
            return []

        lines = result.stdout.strip().split('\n')
        if len(lines) <= 1:
            return []

        pending = []
        for line in lines[1:]:
            if not line.strip() or 'N/A' in line:
                continue
            parts = line.split(None, 2)
            if len(parts) >= 3:
                pending.append({
                    'id': parts[0].strip(),
                    'name': parts[1].strip(),
                    'start': parts[2].strip(),
                })
        return pending

    @property
    def supports_turbomole(self):
        return True

    @property
    def backend_name(self):
        return 'SLURM (BwUniCluster)'
