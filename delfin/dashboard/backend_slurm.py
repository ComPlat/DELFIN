"""SLURM job backend: sbatch / squeue / scancel."""

import os
import re
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple

from delfin.qm_runtime import (
    binary_env_var_name,
    canonical_tool_name,
    resolve_tool,
    settings_selectable_tools,
)
from delfin.runtime_setup import _tool_environment_overrides
from delfin.slurm_submit import (  # noqa: F401 -- normalize_time_limit is re-exported
    PROFILE_PARTITIONS,
    detect_site_profile,
    normalize_time_limit,
    partition_accepts,
    partition_list,
    time_limit_seconds,
)

from .backend_base import JobBackend, JobInfo, SubmitResult
from .input_processing import parse_resource_settings, parse_inp_resources


#: What a pending job's reason code means, in words a user can act on.
_PENDING_REASONS: dict[str, str] = {
    'Priority': 'Waiting for jobs with higher priority to start first.',
    'Resources': 'Next in line: waiting for enough free cores and memory on one node.',
    'None': 'Being scheduled; it should start shortly.',
    'Dependency': 'Waiting for a job it depends on to finish.',
    'DependencyNeverSatisfied': 'Will never start: the job it depends on failed. Cancel it.',
    'BeginTime': 'Not to start before its requested begin time.',
    'JobHeldUser': 'Held by you. Release it with: scontrol release <jobid>',
    'JobHeldAdmin': 'Held by the cluster administrators.',
    'ReqNodeNotAvail': 'The nodes it needs are unavailable, often reserved for maintenance.',
    'Reservation': 'Waiting for a reservation to become active.',
    'PartitionTimeLimit': 'Will never start: its time limit is above the partition maximum.',
    'PartitionNodeLimit': 'Will never start: it asks for more nodes than the partition allows.',
    'BadConstraints': 'Will never start: no node matches what it asks for.',
    'InvalidAccount': 'Will never start: the account is not valid.',
    'InvalidQOS': 'Will never start: the QOS is not valid.',
    'NodeDown': 'A node it needs is down.',
    'Prolog': 'Starting: its node is being prepared.',
    'Cleaning': 'Waiting for a node to finish cleaning up after a previous job.',
}


def describe_pending_reason(reason: str) -> str:
    """The reason SLURM gives for a job's state, explained.

    Takes the bare code (``Priority``), squeue's %R form (``(Priority)``) or a
    code with detail (``ReqNodeNotAvail, Reserved for maintenance``). A code
    nobody explained comes back as it is -- still better than no reason.
    """
    code = str(reason or '').strip()
    if code.startswith('(') and code.endswith(')'):
        code = code[1:-1].strip()
    if not code:
        return ''
    if SlurmJobBackend._is_env_hold(code):
        return ('Its launch failed and SLURM held it. DELFIN releases such a '
                'job on the next refresh of the job list.')
    head = code.split(',', 1)[0].strip()
    if head in _PENDING_REASONS:
        return _PENDING_REASONS[head]
    if head.startswith('QOS') and 'Limit' in head:
        return ('A limit of the queue (QOS) is reached, for example cores per '
                'user; it starts when your other jobs finish.')
    if head.startswith('Assoc') and 'Limit' in head:
        return ("Your account's limit is reached, for example total cores of "
                'the group; it starts when running jobs of the account finish.')
    return code


class SlurmJobBackend(JobBackend):
    """SLURM cluster backend (sbatch/squeue/scancel)."""

    _auto_detected_cache: dict[str, str] | None = None

    # Minimum seconds between two ``squeue`` calls. Several widgets poll job
    # state on their own timer -- the agent activity panel ticks every 5 s --
    # and each uncached call forks a query at the cluster controller, which is
    # a shared, single-instance service. Queue state does not change on that
    # scale, so all callers share one result per interval; submit and cancel
    # drop the cache so the UI still reacts to the user immediately. Sites can
    # widen or (at their own risk) narrow this via DELFIN_SQUEUE_MIN_INTERVAL.
    _JOBS_CACHE_SECONDS = 25.0

    # Tools with expensive system-wide scans (module spider, deep /opt
    # traversal) that are not needed for auto-export to SLURM jobs.
    # Users configure these explicitly in the Settings tab instead.
    _SKIP_AUTO_DETECT = frozenset({"turbomole", "gaussian"})

    @classmethod
    def _auto_detected_tools(cls) -> dict[str, str]:
        """Resolve tool binaries once per process and cache the result."""
        if cls._auto_detected_cache is None:
            detected: dict[str, str] = {}
            for tool_name in settings_selectable_tools():
                if tool_name in cls._SKIP_AUTO_DETECT:
                    continue
                resolved = resolve_tool(tool_name)
                if resolved is not None:
                    detected[tool_name] = resolved.path
            cls._auto_detected_cache = detected
        return cls._auto_detected_cache

    # What every SLURM site gets. Python started from a venv on a network HOME
    # is an I/O problem on any cluster, not on one: the job runs the venv from
    # node-local disk when that disk exists and has room, and from where it is
    # otherwise (Section 7 of submit_delfin.sh). A variable the user exported
    # themselves is left as they set it.
    _GENERIC_ENV: dict[str, str] = {
        'DELFIN_STAGE_VENV': '1',
        'DELFIN_RUNTIME_CACHE': '1',
    }

    # Site-specific environment variables injected into every sbatch call
    # based on the SLURM profile name from settings.
    _PROFILE_ENV: dict[str, dict[str, str]] = {
        'bwunicluster3': {
            'DELFIN_MODULES': 'devel/python/3.11.7-gnu-14.2',
            'DELFIN_STAGE_ORCA': '1',
            'DELFIN_STAGE_VENV': '1',
            'DELFIN_RUNTIME_CACHE': '1',
            'DELFIN_NODE_CORES': '96',
            'DELFIN_NODE_MEM_MB': str(384 * 1024),
            'DELFIN_HIGHMEM_MB': str(2304 * 1024),
        },
    }

    # SLURM holds a job (state PENDING, reason "…requeued held") when the
    # controller's user-environment retrieval fails: it runs a login shell as
    # the user, bounded by GetEnvTimeout (default 2 s), and a slow ~/.bashrc on
    # a busy login node pushes it over. The hold is transient — releasing the
    # job re-attempts retrieval, which almost always succeeds. DELFIN releases
    # such holds automatically so submits are not silently stuck.
    _ENV_HOLD_MARKERS = (
        'user_env_retrieval_failed',
        'launch_failed_requeued_held',
    )

    # Partitions a site's CPU jobs may start in; see delfin.slurm_submit.
    _PROFILE_PARTITIONS: dict[str, tuple[str, ...]] = PROFILE_PARTITIONS

    @staticmethod
    def _detect_profile() -> str:
        """Auto-detect the site profile from hostname / FQDN."""
        return detect_site_profile()

    def __init__(self, submit_templates_dir, orca_base=None, tool_binaries=None,
                 slurm_profile=None, partitions=None):
        self.submit_templates_dir = Path(submit_templates_dir)
        explicit_profile = str(slurm_profile or '').strip()
        self.slurm_profile = explicit_profile or self._detect_profile()
        # Settings first, then the environment, then the site profile. Empty
        # everywhere means the template's own #SBATCH --partition decides.
        self.partitions = (
            self._partition_list(partitions)
            or self._partition_list(os.environ.get('DELFIN_SLURM_PARTITIONS', ''))
            or self._PROFILE_PARTITIONS.get(self.slurm_profile, ())
        )
        self._partition_verdicts: dict[tuple, bool] = {}
        self.orca_base = orca_base
        explicit = {
            canonical_tool_name(name): str(value).strip()
            for name, value in (tool_binaries or {}).items()
            if str(value or '').strip()
        }
        # Only auto-detect if no tools were configured in settings.
        # Users configure tools via Settings tab "Scan" buttons which
        # persist paths to ~/.delfin_settings.json.
        if not explicit:
            for tool_name, path in self._auto_detected_tools().items():
                explicit[tool_name] = path
        self.tool_binaries = explicit
        self._jobs_cache: Optional[List[JobInfo]] = None
        self._jobs_cache_at: Optional[float] = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _resolve_resources(self, job_dir, inp_file=None, pal=40, maxcore=6000):
        """PAL and memory for a job, from the file that job actually reads.

        An ORCA job reads the .inp it is submitted with, so its %pal and
        %maxcore come first. A CONTROL.txt beside it belongs to the DELFIN
        run the folder came from: taken first, a recalc whose input was edited
        down to 12 processes still reserved CONTROL's 40 cores and 240 GB and
        waited for them -- and one edited up ran more processes than it had
        cores. A DELFIN job (no inp_file) reads CONTROL.txt, so there it comes
        first. Whatever is still missing falls back to the other source, then
        to the caller's values.
        """
        pal_used = None
        maxcore_used = None
        job_path = Path(job_dir)

        def fill(pal_found, maxcore_found):
            nonlocal pal_used, maxcore_used
            if pal_used is None:
                pal_used = pal_found
            if maxcore_used is None:
                maxcore_used = maxcore_found

        def from_submitted_inp():
            try:
                candidate = job_path / inp_file if inp_file else None
                if candidate is not None and candidate.exists():
                    fill(*parse_inp_resources(candidate.read_text()))
            except Exception:
                pass

        def from_control():
            try:
                control_path = job_path / 'CONTROL.txt'
                if control_path.exists():
                    fill(*parse_resource_settings(control_path.read_text()))
            except Exception:
                pass

        if inp_file:
            from_submitted_inp()
            from_control()
        else:
            from_control()

        # Still missing: the first .inp of the folder.
        if pal_used is None or maxcore_used is None:
            try:
                inp_files = sorted(job_path.glob('*.inp'))
                if inp_files:
                    fill(*parse_inp_resources(inp_files[0].read_text()))
            except Exception:
                pass

        if pal_used is None:
            pal_used = int(pal)
        if maxcore_used is None:
            maxcore_used = int(maxcore)

        mem_used = int(pal_used) * int(maxcore_used)
        return pal_used, mem_used

    def _append_tool_exports(self, env_vars: str) -> str:
        # Compute all tool env vars: DELFIN_*_BINARY + XTB4STDAHOME, STD2HOME, TURBODIR
        tool_env, _path_entries = _tool_environment_overrides(self.tool_binaries)
        exports = [f'{key}={value}' for key, value in tool_env.items()]
        if not exports:
            return env_vars
        return f'{env_vars},{",".join(exports)}'

    @staticmethod
    def _append_extra_env(env_vars: str, extra_env: dict | None) -> str:
        exports = []
        for key, value in (extra_env or {}).items():
            if value is None:
                continue
            exports.append(f'{key}={value}')
        if not exports:
            return env_vars
        return f'{env_vars},{",".join(exports)}'

    def _resolve_turbomole_command(self, module: str) -> str | None:
        selected = str(self.tool_binaries.get('turbomole') or '').strip()
        if not selected:
            return None
        path = Path(selected).expanduser()
        try:
            path = path.resolve()
        except Exception:
            pass
        roots: list[Path] = []
        if path.is_dir():
            roots.append(path)
        roots.extend(path.parents)
        for root in roots:
            bin_dir = root / 'bin'
            if not bin_dir.is_dir():
                continue
            direct = bin_dir / module
            if direct.is_file() and direct.stat().st_mode & 0o111:
                return str(direct)
            try:
                subdirs = [item for item in bin_dir.iterdir() if item.is_dir()]
            except Exception:
                subdirs = []
            for subdir in subdirs:
                candidate = subdir / module
                if candidate.is_file() and candidate.stat().st_mode & 0o111:
                    return str(candidate)
        return None

    @staticmethod
    def _time_limit_seconds(time_limit: str) -> int:
        """SLURM's time formats in seconds; see delfin.slurm_submit."""
        return time_limit_seconds(time_limit)

    def _append_profile_env(self, env_vars: str) -> str:
        """Inject the generic job env, and the site profile's on top of it."""
        merged = {
            key: value for key, value in self._GENERIC_ENV.items()
            if key not in os.environ
        }
        merged.update(self._PROFILE_ENV.get(self.slurm_profile, {}))
        if not merged:
            return env_vars
        exports = [f'{key}={value}' for key, value in merged.items()]
        return f'{env_vars},{",".join(exports)}'

    @staticmethod
    def _runtime_location_env() -> dict[str, str]:
        """Where this DELFIN runs from, so the job runs the same one.

        The job used to find it by walking up from the submit directory to a
        ``software/delfin`` -- one layout. A pip install, a conda environment
        or a checkout anywhere else was never found, and the job took whatever
        ``python`` its PATH offered.
        """
        env: dict[str, str] = {}
        prefix = Path(sys.prefix)
        if sys.prefix != getattr(sys, 'base_prefix', sys.prefix) or (prefix / 'conda-meta').is_dir():
            env['DELFIN_VENV'] = str(prefix)
        repo = Path(__file__).resolve().parents[2]
        if (repo / 'pyproject.toml').is_file() and (repo / '.git').exists():
            env['DELFIN_REPO'] = str(repo)
        mpirun = shutil.which('mpirun')
        if mpirun:
            ompi_home = Path(mpirun).resolve().parent.parent
            # Only an OpenMPI that stands on its own: the job copies this
            # directory to node-local disk, and /usr is not something to copy.
            if (ompi_home / 'bin' / 'ompi_info').is_file() and str(ompi_home) not in ('/', '/usr', '/usr/local'):
                env['DELFIN_OMPI_HOME'] = str(ompi_home)
        return {key: value for key, value in env.items() if key not in os.environ}

    def _append_runtime_location_env(self, env_vars: str) -> str:
        location = self._runtime_location_env()
        if not location:
            return env_vars
        exports = [f'{key}={value}' for key, value in location.items()]
        return f'{env_vars},{",".join(exports)}'

    @staticmethod
    def _partition_list(value) -> tuple:
        return partition_list(value)

    def _partitions_for(self, job_dir, submit_script, time_limit, pal, mem_mb,
                        submit_env) -> str:
        """The candidate partitions this request fits, as one --partition value.

        A partition whose nodes cannot hold the request -- cpu_il has 64 cores
        and 256 GB per node -- makes SLURM reject the whole job when it is
        listed, so each candidate is asked first with ``sbatch --test-only``,
        which submits nothing. One question per partition per request shape
        for the life of the dashboard. '' keeps the template's own partition.
        """
        candidates = self.partitions
        if not candidates:
            return ''
        if len(candidates) == 1:
            return candidates[0]
        fitting = []
        for name in candidates:
            key = (name, str(time_limit), int(pal), int(mem_mb), str(submit_script))
            verdict = self._partition_verdicts.get(key)
            if verdict is None:
                verdict = self._partition_accepts(
                    name, job_dir, submit_script, time_limit, pal, mem_mb, submit_env)
                if verdict is not None:
                    self._partition_verdicts[key] = verdict
            if verdict:
                fitting.append(name)
        return ','.join(fitting)

    @staticmethod
    def _partition_accepts(name, job_dir, submit_script, time_limit, pal, mem_mb,
                           submit_env) -> Optional[bool]:
        """True or False as SLURM answered, None when it could not be asked."""
        return partition_accepts(
            name, job_dir, submit_script,
            resource_args=(f'--time={time_limit}', '--ntasks=1',
                           f'--cpus-per-task={pal}', f'--mem={mem_mb}M'),
            submit_env=submit_env,
        )

    def _sbatch(self, job_dir, env_vars, time_limit, pal, mem_mb,
                job_name, submit_script, *, gpu=None, partition=None,
                array=None):
        """Run sbatch with the given parameters.

        Parameters
        ----------
        gpu : str or None
            GRES GPU spec, e.g. ``"gpu:1"`` or ``"gpu:a100:1"``.
            Passed as ``--gres=<gpu>`` when set.
        partition : str or None
            SLURM partition name, e.g. ``"gpu"`` or ``"gpu_4"``.
            Passed as ``--partition=<partition>`` when set.
        array : str or None
            SLURM job-array spec, e.g. ``"0-9"`` or ``"0-9%4"`` (max 4
            concurrent). Passed as ``--array=<array>`` when set.
        """
        try:
            time_limit = normalize_time_limit(time_limit)
        except ValueError as problem:
            # Said to the user the way sbatch's own refusals are: callers show
            # stderr when the return code is not 0.
            return subprocess.CompletedProcess(
                args=['sbatch'], returncode=1, stdout='',
                stderr=(f'Invalid time limit {str(time_limit)!r}: {problem}. Use HH:MM:SS, '
                        'D-HH:MM:SS, or a duration such as 48h, 2d or 90min.'),
            )
        env_vars = self._append_profile_env(env_vars)
        env_vars = self._append_runtime_location_env(env_vars)
        # USR1 for preemptive sync: half the walltime, capped 30s–300s
        total_secs = (self._time_limit_seconds(time_limit)
                      if time_limit[:1].isdigit() else 3 * 86400)
        usr1_offset = max(30, min(300, total_secs // 2))
        # Pass DELFIN_* variables via the sbatch process environment, NOT via
        # `--export=ALL,VAR=…`. Passing explicit variables with --export makes
        # slurmd fetch the user's login environment on the compute node at
        # launch (bounded by GetEnvTimeout, default 2 s) — on bwUniCluster3
        # that times out intermittently and the job is requeued/held with
        # "user env retrieval failed". A plain `--export=ALL` propagates the
        # full submission environment (including these variables) without
        # triggering that retrieval. Verified by an A/B test pinned to the
        # same node (uc3n003): plain/ALL-only completed, ALL,VARS failed.
        cmd = [
            'sbatch',
            '--export=ALL',
            f'--time={time_limit}',
            f'--signal=B:USR1@{usr1_offset}',
            '--ntasks=1',
            f'--cpus-per-task={pal}',
            f'--mem={mem_mb}M',
            f'--job-name={job_name}',
        ]
        submit_env = {**os.environ, **self._env_vars_to_dict(env_vars)}
        # Inherited SBATCH_* input variables act like sbatch CLI options.
        # SBATCH_GET_USER_ENV has NO command-line negation, so if a user's
        # shell exports it, env retrieval would come back despite
        # --export=ALL. SBATCH_EXPORT could likewise fight our flag.
        submit_env.pop('SBATCH_GET_USER_ENV', None)
        submit_env.pop('SBATCH_EXPORT', None)
        if gpu:
            cmd.append(f'--gres={gpu}')
        if not partition and not gpu:
            partition = self._partitions_for(
                job_dir, submit_script, time_limit, pal, mem_mb, submit_env)
        if partition:
            cmd.append(f'--partition={partition}')
        if array:
            cmd.append(f'--array={array}')
        cmd.append(str(submit_script))
        result = subprocess.run(
            cmd,
            cwd=str(job_dir),
            capture_output=True,
            text=True,
            env=submit_env,
        )
        if result.returncode == 0:
            # Best-effort recovery for the transient env-retrieval hold.
            # Silent side effect: never raises, never mutates result.stdout
            # (downstream extracts the job id from its last token).
            self._release_env_hold(result.stdout)
            # The queue just changed by our own doing -- show it at once
            # instead of making the user wait out the squeue rate limit.
            self._invalidate_jobs_cache()
        return result

    @staticmethod
    def _env_vars_to_dict(env_vars: str) -> dict[str, str]:
        """Parse the ``KEY=val,KEY2=val2`` export string into a dict.

        Same comma/equals grammar sbatch used for ``--export=ALL,<str>``, so
        the constraint (no commas inside values) is unchanged.
        """
        parsed: dict[str, str] = {}
        for chunk in (env_vars or '').split(','):
            chunk = chunk.strip()
            if not chunk or '=' not in chunk:
                continue
            key, value = chunk.split('=', 1)
            key = key.strip()
            if key:
                parsed[key] = value
        return parsed

    @staticmethod
    def _job_id_from_submit(stdout: str) -> str | None:
        """Extract the numeric job id from ``Submitted batch job N``."""
        match = re.search(r'Submitted batch job (\d+)', stdout or '')
        return match.group(1) if match else None

    @classmethod
    def _is_env_hold(cls, text: str) -> bool:
        """True if ``text`` (an scontrol Reason token or a squeue reason
        field) describes the transient user-env-retrieval hold. ``scontrol``
        writes it with underscores, ``squeue`` with spaces — normalise both."""
        normalized = (text or '').lower().replace('_', ' ')
        return any(marker.replace('_', ' ') in normalized
                   for marker in cls._ENV_HOLD_MARKERS)

    @staticmethod
    def _release_job_quietly(job_id: str) -> None:
        """Best-effort ``scontrol release`` — never raises. Idempotent: a
        release on an already-pending job is a harmless no-op."""
        try:
            subprocess.run(
                ['scontrol', 'release', str(job_id)],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            pass

    def _scontrol_job_state(self, job_id: str) -> str | None:
        """Return ``scontrol show job`` output (lowercased) or None if the job
        cannot be queried (finished, gone, or scontrol unavailable)."""
        try:
            result = subprocess.run(
                ['scontrol', 'show', 'job', str(job_id), '-o'],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            return None
        if result.returncode != 0:
            return None
        text = (result.stdout or '').lower()
        return text if 'jobstate=' in text else None

    def _release_env_hold(self, stdout: str, *, attempts: int = 6,
                          delay: float = 0.5, max_releases: int = 2) -> None:
        """Release a job SLURM held because user-env retrieval failed.

        Polls briefly after submit for the hold to appear, then runs
        ``scontrol release``. Best-effort: never raises, never touches the
        SubmitResult so job-id parsing downstream stays intact.
        """
        job_id = self._job_id_from_submit(stdout)
        if not job_id:
            return
        releases = 0
        try:
            for _ in range(max(1, attempts)):
                info = self._scontrol_job_state(job_id)
                if info is None:
                    return  # not queryable -> nothing to recover
                reason_match = re.search(r'reason=(\S+)', info)
                reason = reason_match.group(1) if reason_match else ''
                if self._is_env_hold(reason):
                    subprocess.run(
                        ['scontrol', 'release', str(job_id)],
                        capture_output=True, text=True, timeout=10,
                    )
                    releases += 1
                    print(
                        f'[DELFIN] Auto-released SLURM env-retrieval hold on '
                        f'job {job_id} (attempt {releases}).',
                        file=sys.stderr, flush=True,
                    )
                    if releases >= max_releases:
                        return
                    time.sleep(delay)
                    continue
                # A concrete non-hold scheduler reason (Priority, Resources, …)
                # or a running/completing state means the controller accepted
                # the job past env retrieval — nothing to do.
                if reason and reason not in ('none', 'null', '(null)'):
                    return
                if ('jobstate=running' in info or 'jobstate=completing' in info
                        or 'jobstate=completed' in info):
                    return
                time.sleep(delay)  # reason not set yet; re-check shortly
        except Exception:
            return

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def submit_delfin(self, job_dir, job_name, mode='delfin',
                      time_limit='48:00:00', pal=40, maxcore=6000,
                      override=None, build_mult=None,
                      co2_species_delta=None, extra_env=None) -> SubmitResult:
        env_vars = f'DELFIN_MODE={mode},DELFIN_JOB_NAME={job_name}'
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        if override:
            env_vars += f',DELFIN_OVERRIDE={override}'
        if build_mult is not None:
            env_vars += f',BUILD_MULTIPLICITY={build_mult}'
        if co2_species_delta is not None:
            env_vars += f',DELFIN_CO2_SPECIES_DELTA={co2_species_delta}'

        pal_used, mem_used = self._resolve_resources(job_dir, pal=pal, maxcore=maxcore)
        maxcore_used = max(1, int(mem_used) // max(1, int(pal_used)))
        env_vars += f',DELFIN_PAL={pal_used},DELFIN_MAXCORE={maxcore_used}'
        env_vars = self._append_extra_env(env_vars, extra_env)
        env_vars = self._append_tool_exports(env_vars)
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def submit_guppy_batch(
        self,
        job_dir,
        job_name,
        smiles_csv,
        *,
        array_size: int,
        array_concurrency: int | None = None,
        time_limit: str = '24:00:00',
        pal: int = 16,
        maxcore: int = 6000,
        start_strategy: str = 'isomers',
        max_isomers: int = 100,
        goat_topk: int = 0,
        rmsd_cutoff: float = 0.3,
        energy_window_kcal: float = 25.0,
        extra_env: dict | None = None,
    ) -> SubmitResult:
        """Submit a SLURM job array: one task per SMILES row.

        Each array task reads the CSV, picks its row via ``SLURM_ARRAY_TASK_ID``
        (1-based), and runs GUPPY in a per-entry subdirectory. Charge is
        always derived from that row's SMILES.
        """
        if array_size < 1:
            raise ValueError('array_size must be >= 1')

        pal_used = max(1, int(pal))
        maxcore_used = max(1, int(maxcore))
        mem_used = pal_used * maxcore_used

        env_vars = (
            f'DELFIN_MODE=guppy_batch,DELFIN_JOB_NAME={job_name},'
            f'GUPPY_BATCH_CSV={smiles_csv},'
            f'GUPPY_START_STRATEGY={start_strategy},'
            f'GUPPY_MAX_ISOMERS={int(max_isomers)},'
            f'GUPPY_GOAT_TOPK={int(goat_topk)},'
            f'GUPPY_RMSD_CUTOFF={float(rmsd_cutoff)},'
            f'GUPPY_ENERGY_WINDOW_KCAL={float(energy_window_kcal)},'
            f'DELFIN_PAL={pal_used},DELFIN_MAXCORE={maxcore_used}'
        )
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        env_vars = self._append_extra_env(env_vars, extra_env)
        env_vars = self._append_tool_exports(env_vars)

        # Array runs as 1-based tasks so SLURM_ARRAY_TASK_ID == CSV row.
        hi = int(array_size)
        if array_concurrency and array_concurrency > 0:
            array_spec = f'1-{hi}%{int(array_concurrency)}'
        else:
            array_spec = f'1-{hi}'

        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
            array=array_spec,
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
        env_vars = self._append_tool_exports(env_vars)
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def submit_hyperpol_xtb(self, job_dir, job_name, xyz_file, label,
                            time_limit='48:00:00', pal=4, maxcore=1000,
                            use_bfw: bool = False) -> SubmitResult:
        pal_used = max(1, int(pal))
        maxcore_used = max(1, int(maxcore))
        mem_used = pal_used * maxcore_used
        env_vars = (
            f'DELFIN_MODE=hyperpol_xtb,DELFIN_JOB_NAME={job_name},'
            f'DELFIN_XYZ_FILE={xyz_file},DELFIN_WORKFLOW_LABEL={label},'
            f'DELFIN_PAL={pal_used},DELFIN_MAXCORE={maxcore_used}'
        )
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        env_vars = self._append_extra_env(env_vars, {
            'DELFIN_HYPERPOL_XTB_BFW': '1' if use_bfw else '0',
        })
        env_vars = self._append_tool_exports(env_vars)
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def submit_tadf_xtb(self, job_dir, job_name, xyz_file, label,
                        time_limit='48:00:00', pal=4, maxcore=1000,
                        use_bfw: bool = False) -> SubmitResult:
        pal_used = max(1, int(pal))
        maxcore_used = max(1, int(maxcore))
        mem_used = pal_used * maxcore_used
        env_vars = (
            f'DELFIN_MODE=tadf_xtb,DELFIN_JOB_NAME={job_name},'
            f'DELFIN_XYZ_FILE={xyz_file},DELFIN_WORKFLOW_LABEL={label},'
            f'DELFIN_PAL={pal_used},DELFIN_MAXCORE={maxcore_used}'
        )
        if self.orca_base:
            env_vars += f',DELFIN_ORCA_BASE={self.orca_base}'
        env_vars = self._append_extra_env(env_vars, {
            'DELFIN_TADF_XTB_PREOPT': 'xtb',
            'DELFIN_TADF_XTB_T1_OPT': 'yes',
            'DELFIN_TADF_XTB_BFW': '1' if use_bfw else '0',
        })
        env_vars = self._append_tool_exports(env_vars)
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
        env_vars = self._append_tool_exports(env_vars)
        tm_command = self._resolve_turbomole_command(module)
        if tm_command:
            env_vars += f',TM_COMMAND={tm_command}'
        result = self._sbatch(
            job_dir, env_vars, time_limit, nprocs, mem_mb,
            job_name, self.submit_templates_dir / 'submit_turbomole.sh',
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    @staticmethod
    def _detect_gpu_partition() -> str:
        """Auto-detect the first available GPU partition on this cluster."""
        try:
            result = subprocess.run(
                ['sinfo', '-h', '-o', '%P', '--partition=gpu,gpu_4,gpu_8,gpu_a100'],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                # sinfo returns partition names, possibly with '*' for default
                return result.stdout.strip().split('\n')[0].rstrip('*')
        except Exception:
            pass
        # Fallback: try generic 'gpu'
        try:
            result = subprocess.run(
                ['sinfo', '-h', '-o', '%P'],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    name = line.strip().rstrip('*')
                    if 'gpu' in name.lower():
                        return name
        except Exception:
            pass
        return ''

    def submit_mlp(self, job_dir, job_name, xyz_file, backend='ani2x',
                   time_limit='24:00:00', pal=4, maxcore=4000,
                   charge=0, mult=1) -> SubmitResult:
        """Submit an MLP job — automatically requests GPU if available."""
        pal_used = max(1, int(pal))
        maxcore_used = max(1, int(maxcore))
        mem_used = pal_used * maxcore_used
        env_vars = (
            f'DELFIN_MODE=mlp,DELFIN_JOB_NAME={job_name},'
            f'DELFIN_XYZ_FILE={xyz_file},'
            f'DELFIN_MLP_BACKEND={backend},'
            f'DELFIN_CHARGE={charge},DELFIN_MULT={mult},'
            f'DELFIN_PAL={pal_used},DELFIN_MAXCORE={maxcore_used}'
        )
        # Auto-detect GPU partition — if found, request 1 GPU
        gpu_partition = self._detect_gpu_partition()
        env_vars = self._append_tool_exports(env_vars)
        result = self._sbatch(
            job_dir, env_vars, time_limit, pal_used, mem_used,
            job_name, self.submit_templates_dir / 'submit_delfin.sh',
            gpu='gpu:1' if gpu_partition else None,
            partition=gpu_partition or None,
        )
        return SubmitResult(result.returncode, result.stdout, result.stderr)

    def _jobs_cache_seconds(self) -> float:
        raw = os.environ.get('DELFIN_SQUEUE_MIN_INTERVAL', '').strip()
        if raw:
            try:
                return max(0.0, float(raw))
            except ValueError:
                pass
        return self._JOBS_CACHE_SECONDS

    def _invalidate_jobs_cache(self) -> None:
        """Force the next list_jobs() to ask SLURM again."""
        self._jobs_cache = None
        self._jobs_cache_at = None

    def list_jobs(self, force: bool = False) -> List[JobInfo]:
        """Active jobs of this user, at most one squeue per cache interval.

        Pass ``force=True`` for an explicit user-triggered refresh; periodic
        widgets must not, or the rate limit has no effect.
        """
        # The gate is the time of the last *attempt*, not of the last success:
        # a failing squeue must be rate-limited too, otherwise a busy
        # controller gets hit at full tick rate exactly when it can least
        # afford it.
        if not force and self._jobs_cache_at is not None:
            age = time.monotonic() - self._jobs_cache_at
            if 0.0 <= age < self._jobs_cache_seconds():
                return list(self._jobs_cache or [])

        jobs = self._query_jobs()
        self._jobs_cache_at = time.monotonic()
        if jobs is None:
            # Hold the last good answer rather than flashing an empty queue
            # at the user over one hiccup.
            return list(self._jobs_cache or [])
        self._jobs_cache = jobs
        return list(jobs)

    #: Everything the job table shows, from one squeue. Why a job waits and
    #: SLURM's estimate of its start (%r, %S) needed a second, unthrottled
    #: ``squeue --start`` on every refresh before. '|' separates the fields:
    #: a node list or a reason can contain spaces.
    _SQUEUE_FORMAT = '%i|%P|%j|%u|%t|%M|%D|%l|%C|%m|%S|%r|%R'

    def _query_jobs(self) -> Optional[List[JobInfo]]:
        """Ask SLURM. Returns None when the query failed, [] when idle."""
        try:
            result = subprocess.run(
                ['squeue', '-u', os.environ.get('USER', ''),
                 '-o', self._SQUEUE_FORMAT],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            return None

        if result.returncode != 0:
            return None

        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return []

        jobs = []
        for line in lines[1:]:
            if not line.strip():
                continue
            fields = self._parse_squeue_line(line)
            if fields is None:
                continue
            if self._is_env_hold(line):
                # Dispatch-time env-retrieval holds recur after submit; clear
                # them on every refresh so DELFIN jobs don't stay stuck.
                self._release_job_quietly(fields['job_id'])

            jobs.append(JobInfo(
                job_id=fields['job_id'],
                name=fields['name'],
                status=fields['status'],
                start_time=fields['start'],
                time_limit=fields['time_limit'],
                extra={
                    'partition': fields['partition'],
                    'user': fields['user'],
                    'time_used': fields['time_used'],
                    'nodes': fields['nodes'],
                    'reason': fields['reason'],
                    'reason_code': fields['reason_code'],
                    'time_limit': fields['time_limit'],
                    'cpus': fields['cpus'],
                    'memory': fields['memory'],
                    'start_estimate': fields['start'],
                    'raw_line': line,
                },
            ))
        return jobs

    @staticmethod
    def _parse_squeue_line(line: str) -> Optional[dict]:
        """One squeue line, in the '|' format or the plain one of older output."""
        if line.count('|') >= 12:
            parts = line.split('|')
            job_id, partition = (s.strip() for s in parts[:2])
            name = '|'.join(parts[2:-10]).strip()
            (user, status, time_used, nodes, time_limit, cpus, memory, start,
             reason_code, reason) = (s.strip() for s in parts[-10:])
        else:
            parts = line.split(None, 7)
            if not parts:
                return None
            parts += [''] * (8 - len(parts))
            job_id, partition, name, user, status, time_used, nodes, reason = (
                s.strip() for s in parts)
            time_limit = cpus = memory = start = ''
            reason_code = reason[1:-1] if reason.startswith('(') and reason.endswith(')') else ''
        if start in ('N/A', 'Unknown', 'None'):
            start = ''
        return {
            'job_id': job_id, 'partition': partition, 'name': name or 'unknown',
            'user': user, 'status': status, 'time_used': time_used, 'nodes': nodes,
            'time_limit': time_limit, 'cpus': cpus, 'memory': memory, 'start': start,
            'reason_code': reason_code, 'reason': reason,
        }

    def cancel_job(self, job_id) -> Tuple[bool, str]:
        try:
            result = subprocess.run(
                ['scancel', str(job_id)],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                self._invalidate_jobs_cache()
                return True, f'Job {job_id} cancelled successfully.'
            else:
                return False, f'Error cancelling job {job_id}: {result.stderr or result.stdout}'
        except Exception as e:
            return False, f'Error cancelling job {job_id}: {e}'

    def get_pending_start_times(self):
        """Pending jobs with SLURM's start estimate and why they wait.

        Read from the job listing, which carries both, so asking costs no
        squeue of its own.
        """
        pending = []
        for job in self.list_jobs():
            if (job.status or '').upper() not in ('PD', 'PENDING'):
                continue
            extra = job.extra or {}
            reason = extra.get('reason_code') or extra.get('reason', '')
            pending.append({
                'id': str(job.job_id),
                'name': job.name,
                'start': extra.get('start_estimate', ''),
                'reason': reason,
                'explanation': describe_pending_reason(reason),
                'partition': extra.get('partition', ''),
            })
        return pending

    @property
    def supports_turbomole(self):
        return True

    @property
    def backend_name(self):
        return 'SLURM (BwUniCluster)'
