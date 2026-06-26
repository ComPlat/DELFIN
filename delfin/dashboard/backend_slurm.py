"""SLURM job backend: sbatch / squeue / scancel."""

import os
import socket
import subprocess
from pathlib import Path
from typing import List, Tuple

from delfin.qm_runtime import (
    binary_env_var_name,
    canonical_tool_name,
    resolve_tool,
    settings_selectable_tools,
)
from delfin.runtime_setup import _tool_environment_overrides

from .backend_base import JobBackend, JobInfo, SubmitResult
from .input_processing import parse_resource_settings, parse_inp_resources


class SlurmJobBackend(JobBackend):
    """SLURM cluster backend (sbatch/squeue/scancel)."""

    _auto_detected_cache: dict[str, str] | None = None

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

    @staticmethod
    def _detect_profile() -> str:
        """Auto-detect the site profile from hostname / FQDN."""
        try:
            fqdn = socket.getfqdn().lower()
        except Exception:
            fqdn = ""
        hostname = os.environ.get("HOSTNAME", "").lower() or fqdn.split(".")[0]
        if "scc.kit.edu" in fqdn or hostname.startswith("uc3"):
            return "bwunicluster3"
        return ""

    def __init__(self, submit_templates_dir, orca_base=None, tool_binaries=None,
                 slurm_profile=None):
        self.submit_templates_dir = Path(submit_templates_dir)
        explicit_profile = str(slurm_profile or '').strip()
        self.slurm_profile = explicit_profile or self._detect_profile()
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
        """Convert HH:MM:SS or MM:SS or bare minutes to total seconds."""
        parts = str(time_limit).split(':')
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
        elif len(parts) == 2:
            return int(parts[0]) * 60 + int(parts[1])
        return int(parts[0]) * 60  # SLURM treats bare number as minutes

    def _append_profile_env(self, env_vars: str) -> str:
        """Inject site-specific env vars based on the SLURM profile."""
        profile_env = self._PROFILE_ENV.get(self.slurm_profile, {})
        if not profile_env:
            return env_vars
        exports = [f'{key}={value}' for key, value in profile_env.items()]
        return f'{env_vars},{",".join(exports)}'

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
        env_vars = self._append_profile_env(env_vars)
        # USR1 for preemptive sync: half the walltime, capped 30s–300s
        total_secs = self._time_limit_seconds(time_limit)
        usr1_offset = max(30, min(300, total_secs // 2))
        cmd = [
            'sbatch',
            f'--export=ALL,{env_vars}',
            f'--time={time_limit}',
            f'--signal=B:USR1@{usr1_offset}',
            '--ntasks=1',
            f'--cpus-per-task={pal}',
            f'--mem={mem_mb}M',
            f'--job-name={job_name}',
        ]
        if gpu:
            cmd.append(f'--gres={gpu}')
        if partition:
            cmd.append(f'--partition={partition}')
        if array:
            cmd.append(f'--array={array}')
        cmd.append(str(submit_script))
        return subprocess.run(
            cmd,
            cwd=str(job_dir),
            capture_output=True,
            text=True,
        )

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
