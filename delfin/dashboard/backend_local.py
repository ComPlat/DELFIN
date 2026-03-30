"""Local job backend: JSON-based queue with a background worker thread."""

import importlib.util
import json
import os
import signal
import subprocess
import sys
import threading
import time as _time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

from delfin.qm_runtime import binary_env_var_name, canonical_tool_name

from .backend_base import JobBackend, JobInfo, SubmitResult
from .helpers import parse_time_to_seconds


class LocalJobBackend(JobBackend):
    """Local server backend using a JSON job queue and subprocess execution."""

    def __init__(self, run_script, orca_base=None,
                 jobs_file=None, tool_binaries=None,
                 max_cores=384, max_ram_mb=1_400_000,
                 allow_oversubscribe=False, oversubscribe_factor=1.0):
        self.run_script = Path(run_script)
        self.orca_base = orca_base
        self.jobs_file = Path(jobs_file) if jobs_file else Path.home() / '.delfin_jobs.json'
        self.tool_binaries = {
            canonical_tool_name(name): str(value).strip()
            for name, value in (tool_binaries or {}).items()
            if str(value or '').strip()
        }
        self.max_cores = max_cores
        self.max_ram_mb = max_ram_mb
        self.allow_oversubscribe = bool(allow_oversubscribe)
        try:
            factor = float(oversubscribe_factor)
        except (TypeError, ValueError):
            factor = 1.0
        self.oversubscribe_factor = max(1.0, factor)
        self._next_job_id_key = '_next_job_id'
        self._lock = threading.Lock()
        self._worker_running = True
        self._worker_thread = threading.Thread(
            target=self._queue_worker, daemon=True, name='delfin-queue'
        )
        self._worker_thread.start()

    def _core_budget(self) -> int:
        budget = int(self.max_cores)
        if not self.allow_oversubscribe:
            return budget
        boosted = int(round(float(self.max_cores) * float(self.oversubscribe_factor)))
        return max(budget, boosted)

    def _ram_budget(self) -> int:
        budget = int(self.max_ram_mb)
        if not self.allow_oversubscribe:
            return budget
        boosted = int(round(float(self.max_ram_mb) * float(self.oversubscribe_factor)))
        return max(budget, boosted)

    def _has_launch_target(self):
        return self.run_script.exists() or importlib.util.find_spec(
            'delfin.dashboard.local_runner'
        ) is not None

    def _build_launch_command(self, timeout_secs):
        if self.run_script.exists():
            return ['timeout', str(timeout_secs), 'bash', str(self.run_script)]
        return [
            'timeout',
            str(timeout_secs),
            sys.executable,
            '-m',
            'delfin.dashboard.local_runner',
        ]

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
    def _list_processes(self):
        """Return a best-effort snapshot of current processes."""
        try:
            result = subprocess.run(
                ['ps', '-eo', 'pid=,pgid=,args='],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        except Exception:
            return []
        if result.returncode != 0:
            return []

        processes = []
        for line in result.stdout.splitlines():
            parts = line.strip().split(None, 2)
            if len(parts) < 3:
                continue
            try:
                proc_pid = int(parts[0])
                proc_pgid = int(parts[1])
            except ValueError:
                continue
            try:
                cwd = os.readlink(Path('/proc') / str(proc_pid) / 'cwd')
            except Exception:
                cwd = None
            processes.append((proc_pid, proc_pgid, parts[2], cwd))
        return processes

    def _job_has_active_processes(self, job, process_table=None):
        """Detect live descendants or subprocesses that still belong to a job."""
        job_dir = str(job.get('job_dir') or '').strip()
        if not job_dir:
            return False
        job_dir_prefix = job_dir.rstrip(os.sep) + os.sep

        try:
            wrapper_pid = int(job.get('pid')) if job.get('pid') is not None else None
        except (TypeError, ValueError):
            wrapper_pid = None
        try:
            wrapper_pgid = int(job.get('pgid')) if job.get('pgid') is not None else None
        except (TypeError, ValueError):
            wrapper_pgid = None

        current_pid = os.getpid()
        for proc_pid, proc_pgid, args, cwd in (process_table or self._list_processes()):
            if proc_pid == current_pid:
                continue
            if wrapper_pid is not None and proc_pid == wrapper_pid:
                return True
            if wrapper_pgid is not None and proc_pgid == wrapper_pgid:
                return True
            if job_dir in args:
                return True
            if cwd == job_dir or cwd.startswith(job_dir_prefix):
                return True
        return False

    def _update_job_status(self, job, process_table=None):
        if job['status'] not in ('RUNNING',):
            return job
        pid = job.get('pid')
        if pid is not None:
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

        if self._job_has_active_processes(job, process_table=process_table):
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
        if job.get('pal') is not None:
            env['DELFIN_PAL'] = str(job['pal'])
        if job.get('maxcore') is not None:
            env['DELFIN_MAXCORE'] = str(job['maxcore'])
        if job.get('xyz_file'):
            env['DELFIN_XYZ_FILE'] = str(job['xyz_file'])
        if job.get('workflow_label'):
            env['DELFIN_WORKFLOW_LABEL'] = str(job['workflow_label'])
        for key, value in (job.get('extra_env') or {}).items():
            if value is None:
                continue
            env[str(key)] = str(value)
        job_tool_binaries = job.get('tool_binaries', {}) or self.tool_binaries
        for tool_name, binary_path in job_tool_binaries.items():
            env[binary_env_var_name(tool_name)] = str(binary_path)

        try:
            log_file = Path(job['job_dir']) / f'delfin_{job["job_id"]}.out'
            log_fd = open(log_file, 'w')
            proc = subprocess.Popen(
                self._build_launch_command(timeout_secs),
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
            process_table = self._list_processes()

            changed = False
            for job in jobs:
                if job['status'] == 'RUNNING':
                    old_status = job['status']
                    self._update_job_status(job, process_table=process_table)
                    if job['status'] != old_status:
                        changed = True

            used_cores = sum(j.get('pal', 0) for j in jobs if j['status'] == 'RUNNING')
            used_ram = sum(j.get('pal', 0) * j.get('maxcore', 0)
                          for j in jobs if j['status'] == 'RUNNING')
            core_budget = self._core_budget()
            ram_budget = self._ram_budget()

            for job in jobs:
                if job['status'] != 'PENDING':
                    continue
                needed_cores = job.get('pal', 40)
                needed_ram = job.get('pal', 40) * job.get('maxcore', 6000)
                if (used_cores + needed_cores <= core_budget and
                        used_ram + needed_ram <= ram_budget):
                    if self._start_job(job, data):
                        used_cores += needed_cores
                        used_ram += needed_ram
                        changed = True
            if changed:
                self._save_jobs(data)
            return changed

    def _queue_worker(self):
        while self._worker_running:
            try:
                self._try_start_pending_jobs()
            except Exception:
                pass
            _time.sleep(30)

    # ------------------------------------------------------------------
    # Internal submit (enqueue)
    # ------------------------------------------------------------------
    def _enqueue(self, job_dir, mode, job_name, inp_file=None,
                 time_limit='48:00:00', override=None, build_mult=None,
                 pal=40, maxcore=6000, co2_species_delta=None,
                 xyz_file=None, workflow_label=None, extra_env=None):
        if not self._has_launch_target():
            return SubmitResult(
                1,
                stderr=(
                    'Neither the local run script nor the bundled Python '
                    'runner is available.'
                ),
            )

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
                'xyz_file': xyz_file,
                'workflow_label': workflow_label,
                'extra_env': dict(extra_env or {}),
                'tool_binaries': dict(self.tool_binaries),
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
                      co2_species_delta=None, extra_env=None) -> SubmitResult:
        return self._enqueue(
            job_dir, mode, job_name,
            time_limit=time_limit, override=override,
            build_mult=build_mult, pal=pal, maxcore=maxcore,
            co2_species_delta=co2_species_delta, extra_env=extra_env,
        )

    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=40, maxcore=6000) -> SubmitResult:
        return self._enqueue(
            job_dir, 'orca', job_name, inp_file=inp_file,
            time_limit=time_limit, pal=pal, maxcore=maxcore,
        )

    def submit_hyperpol_xtb(self, job_dir, job_name, xyz_file, label,
                            time_limit='48:00:00', pal=4, maxcore=1000) -> SubmitResult:
        return self._enqueue(
            job_dir, 'hyperpol_xtb', job_name,
            time_limit=time_limit, pal=pal, maxcore=maxcore,
            xyz_file=xyz_file, workflow_label=label,
        )

    def submit_tadf_xtb(self, job_dir, job_name, xyz_file, label,
                        time_limit='48:00:00', pal=4, maxcore=1000) -> SubmitResult:
        return self._enqueue(
            job_dir, 'tadf_xtb', job_name,
            time_limit=time_limit, pal=pal, maxcore=maxcore,
            xyz_file=xyz_file, workflow_label=label,
            extra_env={
                'DELFIN_TADF_XTB_PREOPT': 'xtb',
                'DELFIN_TADF_XTB_T1_OPT': 'yes',
            },
        )

    def list_jobs(self) -> List[JobInfo]:
        with self._lock:
            data = self._load_jobs()
            jobs = data.get('jobs', [])
            process_table = self._list_processes()
            for job in jobs:
                if job['status'] == 'RUNNING':
                    self._update_job_status(job, process_table=process_table)
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
