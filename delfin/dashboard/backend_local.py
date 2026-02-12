"""Local job backend: JSON-based queue with a background worker thread."""

import json
import os
import signal
import subprocess
import threading
import time as _time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

from .backend_base import JobBackend, JobInfo, SubmitResult
from .helpers import parse_time_to_seconds


class LocalJobBackend(JobBackend):
    """Local server backend using a JSON job queue and subprocess execution."""

    def __init__(self, run_script, orca_base=None,
                 jobs_file=None,
                 max_cores=384, max_ram_mb=1_400_000):
        self.run_script = Path(run_script)
        self.orca_base = orca_base
        self.jobs_file = Path(jobs_file) if jobs_file else Path.home() / '.delfin_jobs.json'
        self.max_cores = max_cores
        self.max_ram_mb = max_ram_mb
        self._next_job_id_key = '_next_job_id'
        self._lock = threading.Lock()
        self._worker_running = True
        self._worker_thread = threading.Thread(
            target=self._queue_worker, daemon=True, name='delfin-queue'
        )
        self._worker_thread.start()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------
    def _load_jobs(self):
        if self.jobs_file.exists():
            try:
                return json.loads(self.jobs_file.read_text())
            except (json.JSONDecodeError, IOError):
                return {self._next_job_id_key: 1001, 'jobs': []}
        return {self._next_job_id_key: 1001, 'jobs': []}

    def _save_jobs(self, data):
        self.jobs_file.write_text(json.dumps(data, indent=2, default=str))

    # ------------------------------------------------------------------
    # Status helpers
    # ------------------------------------------------------------------
    def _update_job_status(self, job):
        if job['status'] not in ('RUNNING',):
            return job
        pid = job.get('pid')
        if pid is None:
            job['status'] = 'FAILED'
            return job
        try:
            wpid, _ = os.waitpid(pid, os.WNOHANG)
            if wpid == 0:
                return job
        except ChildProcessError:
            try:
                os.kill(pid, 0)
                return job
            except ProcessLookupError:
                pass
            except PermissionError:
                return job

        job_dir = job.get('job_dir', '')
        job_id = job.get('job_id', 0)
        exit_code_file = Path(job_dir) / f'.exit_code_{job_id}'
        if exit_code_file.exists():
            try:
                code = int(exit_code_file.read_text().strip())
                if code == 0:
                    job['status'] = 'COMPLETED'
                elif code == 124:
                    job['status'] = 'TIMEOUT'
                else:
                    job['status'] = 'FAILED'
            except (ValueError, IOError):
                job['status'] = 'FAILED'
        else:
            job['status'] = 'FAILED'
        return job

    def _get_running_resources(self):
        data = self._load_jobs()
        total_cores = 0
        total_ram_mb = 0
        for job in data.get('jobs', []):
            if job['status'] == 'RUNNING':
                total_cores += job.get('pal', 0)
                total_ram_mb += job.get('pal', 0) * job.get('maxcore', 0)
        return total_cores, total_ram_mb

    # ------------------------------------------------------------------
    # Job start / queue management
    # ------------------------------------------------------------------
    def _start_job(self, job, data):
        timeout_secs = parse_time_to_seconds(job.get('time_limit', '48:00:00'))
        env = os.environ.copy()
        env['DELFIN_MODE'] = job['mode']
        env['DELFIN_JOB_NAME'] = job['name']
        env['DELFIN_JOB_ID'] = str(job['job_id'])
        if self.orca_base:
            env['DELFIN_ORCA_BASE'] = self.orca_base
        if job.get('inp_file'):
            env['DELFIN_INP_FILE'] = job['inp_file']
        if job.get('override'):
            env['DELFIN_OVERRIDE'] = job['override']
        if job.get('build_mult') is not None:
            env['BUILD_MULTIPLICITY'] = str(job['build_mult'])
        if job.get('co2_species_delta') is not None:
            env['DELFIN_CO2_SPECIES_DELTA'] = str(job['co2_species_delta'])

        try:
            log_file = Path(job['job_dir']) / f'delfin_{job["job_id"]}.out'
            log_fd = open(log_file, 'w')
            proc = subprocess.Popen(
                ['timeout', str(timeout_secs), 'bash', str(self.run_script)],
                cwd=str(job['job_dir']),
                env=env,
                stdout=log_fd,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            job['pid'] = proc.pid
            job['pgid'] = os.getpgid(proc.pid)
            job['status'] = 'RUNNING'
            job['start_time'] = datetime.now().isoformat()
            job['log_file'] = str(log_file)
            return True
        except Exception as e:
            job['status'] = 'FAILED'
            job['error'] = str(e)
            return False

    def _try_start_pending_jobs(self):
        with self._lock:
            data = self._load_jobs()
            jobs = data.get('jobs', [])

            for job in jobs:
                if job['status'] == 'RUNNING':
                    self._update_job_status(job)

            used_cores = sum(j.get('pal', 0) for j in jobs if j['status'] == 'RUNNING')
            used_ram = sum(j.get('pal', 0) * j.get('maxcore', 0)
                          for j in jobs if j['status'] == 'RUNNING')

            started_any = False
            for job in jobs:
                if job['status'] != 'PENDING':
                    continue
                needed_cores = job.get('pal', 40)
                needed_ram = job.get('pal', 40) * job.get('maxcore', 6000)
                if (used_cores + needed_cores <= self.max_cores and
                        used_ram + needed_ram <= self.max_ram_mb):
                    if self._start_job(job, data):
                        used_cores += needed_cores
                        used_ram += needed_ram
                        started_any = True
            self._save_jobs(data)
            return started_any

    def _queue_worker(self):
        while self._worker_running:
            try:
                self._try_start_pending_jobs()
            except Exception:
                pass
            _time.sleep(5)

    # ------------------------------------------------------------------
    # Internal submit (enqueue)
    # ------------------------------------------------------------------
    def _enqueue(self, job_dir, mode, job_name, inp_file=None,
                 time_limit='48:00:00', override=None, build_mult=None,
                 pal=40, maxcore=6000, co2_species_delta=None):
        if not self.run_script.exists():
            return SubmitResult(1, stderr=f'Run script not found: {self.run_script}')

        with self._lock:
            data = self._load_jobs()
            job_id = data.get(self._next_job_id_key, 1001)
            data[self._next_job_id_key] = job_id + 1

            job_entry = {
                'job_id': job_id,
                'name': job_name,
                'mode': mode,
                'pid': None,
                'pgid': None,
                'status': 'PENDING',
                'submit_time': datetime.now().isoformat(),
                'start_time': None,
                'job_dir': str(job_dir),
                'time_limit': time_limit,
                'pal': pal,
                'maxcore': maxcore,
                'inp_file': inp_file,
                'override': override,
                'build_mult': build_mult,
                'co2_species_delta': co2_species_delta,
                'log_file': None,
            }
            data['jobs'].append(job_entry)
            self._save_jobs(data)

        self._try_start_pending_jobs()
        return SubmitResult(0, stdout=f'Submitted local job {job_id}')

    # ------------------------------------------------------------------
    # Public JobBackend interface
    # ------------------------------------------------------------------
    def submit_delfin(self, job_dir, job_name, mode='delfin',
                      time_limit='48:00:00', pal=40, maxcore=6000,
                      override=None, build_mult=None,
                      co2_species_delta=None) -> SubmitResult:
        return self._enqueue(
            job_dir, mode, job_name,
            time_limit=time_limit, override=override,
            build_mult=build_mult, pal=pal, maxcore=maxcore,
            co2_species_delta=co2_species_delta,
        )

    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=40, maxcore=6000) -> SubmitResult:
        return self._enqueue(
            job_dir, 'orca', job_name, inp_file=inp_file,
            time_limit=time_limit, pal=pal, maxcore=maxcore,
        )

    def list_jobs(self) -> List[JobInfo]:
        with self._lock:
            data = self._load_jobs()
            jobs = data.get('jobs', [])
            for job in jobs:
                if job['status'] == 'RUNNING':
                    self._update_job_status(job)
            self._save_jobs(data)

        result = []
        for job in jobs:
            if job['status'] in ('RUNNING', 'PENDING'):
                result.append(JobInfo(
                    job_id=job.get('job_id'),
                    name=job.get('name', ''),
                    mode=job.get('mode', ''),
                    status=job.get('status', ''),
                    submit_time=job.get('submit_time', ''),
                    start_time=job.get('start_time', ''),
                    pal=job.get('pal', 0),
                    maxcore=job.get('maxcore', 0),
                    time_limit=job.get('time_limit', ''),
                    job_dir=job.get('job_dir', ''),
                ))
        return result

    def cancel_job(self, job_id) -> Tuple[bool, str]:
        with self._lock:
            data = self._load_jobs()
            for job in data['jobs']:
                if job['job_id'] == job_id:
                    if job['status'] == 'PENDING':
                        job['status'] = 'CANCELLED'
                        self._save_jobs(data)
                        return True, f'Pending job {job_id} cancelled'
                    elif job['status'] != 'RUNNING':
                        return False, f'Job {job_id} is not running (status: {job["status"]})'
                    pgid = job.get('pgid')
                    pid = job.get('pid')
                    try:
                        if pgid:
                            os.killpg(pgid, signal.SIGTERM)
                        elif pid:
                            os.kill(pid, signal.SIGTERM)
                        job['status'] = 'CANCELLED'
                        self._save_jobs(data)
                        return True, f'Job {job_id} cancelled'
                    except ProcessLookupError:
                        job['status'] = 'CANCELLED'
                        self._save_jobs(data)
                        return True, f'Job {job_id} already finished, marked cancelled'
                    except Exception as e:
                        return False, f'Error cancelling job {job_id}: {e}'
            return False, f'Job {job_id} not found'

    def get_pending_start_times(self):
        """Return pending jobs with their queue position."""
        jobs = self.list_jobs()
        pending = [j for j in jobs if j.status == 'PENDING']
        return [{'id': j.job_id, 'name': j.name, 'position': i + 1}
                for i, j in enumerate(pending)]
