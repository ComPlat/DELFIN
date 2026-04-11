"""Central local job runner for all DELFIN execution modes."""

from __future__ import annotations

import os
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path


_ACTIVE_CHILD = None
_ACTIVE_CHILD_LOCK = threading.RLock()


def _job_id() -> str:
    return str(os.environ.get('DELFIN_JOB_ID', '0') or '0')


def _write_exit_code(code: int) -> None:
    (Path.cwd() / f'.exit_code_{_job_id()}').write_text(f'{int(code)}\n', encoding='utf-8')


def _set_active_child(proc) -> None:
    global _ACTIVE_CHILD
    with _ACTIVE_CHILD_LOCK:
        _ACTIVE_CHILD = proc


def _clear_active_child(proc=None) -> None:
    global _ACTIVE_CHILD
    with _ACTIVE_CHILD_LOCK:
        if proc is None or _ACTIVE_CHILD is proc:
            _ACTIVE_CHILD = None


def _signal_active_child(signum: int, *, kill_after_seconds: float = 10.0) -> None:
    with _ACTIVE_CHILD_LOCK:
        proc = _ACTIVE_CHILD
    if proc is None:
        return
    try:
        if proc.poll() is not None:
            _clear_active_child(proc)
            return
    except Exception:
        return

    pgid = None
    try:
        pgid = os.getpgid(proc.pid)
    except Exception:
        pgid = None

    delivered = False
    if pgid and pgid > 0:
        try:
            os.killpg(pgid, signum)
            delivered = True
        except ProcessLookupError:
            _clear_active_child(proc)
            return
        except Exception:
            pass
    if not delivered:
        try:
            proc.send_signal(signum)
            delivered = True
        except ProcessLookupError:
            _clear_active_child(proc)
            return
        except Exception:
            pass

    deadline = time.time() + max(0.0, float(kill_after_seconds))
    while time.time() < deadline:
        try:
            if proc.poll() is not None:
                _clear_active_child(proc)
                return
        except Exception:
            break
        time.sleep(0.1)

    try:
        if proc.poll() is None:
            if pgid and pgid > 0:
                os.killpg(pgid, signal.SIGKILL)
            else:
                proc.kill()
            proc.wait(timeout=5)
    except Exception:
        pass
    finally:
        _clear_active_child(proc)


def _handle_termination(signum, _frame):
    print(f'Received signal {signum}; stopping local DELFIN job.')
    _signal_active_child(signum)
    _write_exit_code(124)
    raise SystemExit(124)


def _configure_environment() -> None:
    # Configure logging to stdout so messages appear in delfin_*.out on SLURM.
    # Without this, logger.error() calls in orca.py are silently dropped.
    try:
        from delfin.common.logging import configure_logging
        configure_logging(stream=sys.stdout, force=True)
    except Exception:
        import logging
        logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                            format='%(levelname)s: %(message)s', force=True)

    os.environ['OMPI_MCA_pml'] = 'ob1'
    os.environ['OMPI_MCA_btl'] = 'self,tcp,vader'
    os.environ['OMPI_MCA_mpi_show_mca_params_file'] = '0'
    os.environ['OMPI_MCA_mpi_yield_when_idle'] = '1'
    os.environ['OMPI_MCA_coll_hcoll_enable'] = '0'
    os.environ.pop('OMPI_MCA_btl_tcp_if_include', None)
    os.environ.pop('OMPI_MCA_btl_tcp_if_exclude', None)
    os.environ['OMPI_MCA_hwloc_base_binding_policy'] = 'none'
    os.environ['OMPI_MCA_rmaps_base_mapping_policy'] = 'core'
    os.environ['OMPI_MCA_rmaps_base_oversubscribe'] = 'true'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['DELFIN_ORCA_PROGRESS'] = '0'
    os.environ['MPLBACKEND'] = 'Agg'


def _resolve_orca_bin() -> str | None:
    orca_base = str(os.environ.get('DELFIN_ORCA_BASE') or '').strip()
    if orca_base:
        candidate = Path(orca_base).expanduser() / 'orca'
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    candidate = shutil.which('orca')
    return candidate if candidate else None


def _delfin_version_string() -> str:
    """Return DELFIN version including runtime key (git hash) when available."""
    try:
        from delfin import __version__
        ver = __version__
    except Exception:
        ver = 'unknown'
    runtime_key = os.environ.get('DELFIN_RUNTIME_KEY', '')
    return f'{ver} ({runtime_key})' if runtime_key else ver


def _print_job_banner(mode: str, job_name: str, inp_file: str) -> None:
    print('========================================')
    print(f'DELFIN Local Job - {job_name}')
    print(f'Mode: {mode}')
    print('========================================')
    print(f'Job ID:      {_job_id()}')
    print(f'Host:        {os.uname().nodename}')
    print(f'CPUs avail:  {os.cpu_count() or 1}')
    print(f'Working Dir: {Path.cwd()}')
    print(f'Started:     {datetime.now().isoformat(timespec="seconds")}')
    if inp_file:
        print(f'Input File:  {inp_file}')
    print(f'DELFIN:      {_delfin_version_string()}')
    print('========================================')
    print('')
    sys.stdout.flush()


def _run_command(cmd: list[str], *, cwd: Path | None = None) -> int:
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        start_new_session=True,
    )
    _set_active_child(proc)
    try:
        return proc.wait()
    finally:
        _clear_active_child(proc)


def _run_and_tee(cmd: list[str], output_path: Path, *, cwd: Path | None = None) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as handle:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        _set_active_child(proc)
        assert proc.stdout is not None
        try:
            for line in proc.stdout:
                sys.stdout.write(line)
                handle.write(line)
            return proc.wait()
        finally:
            _clear_active_child(proc)


def _detect_mode(mode: str) -> str:
    if mode != 'auto':
        return mode
    if Path('CONTROL.txt').exists() and Path('input.txt').exists():
        print('Auto-detected mode: DELFIN (CONTROL.txt + input.txt found)')
        return 'delfin'
    inp_files = sorted(Path.cwd().glob('*.inp'))
    if inp_files:
        print('Auto-detected mode: ORCA (*.inp files found)')
        return 'orca'
    raise RuntimeError(
        'No valid input files found. Expected CONTROL.txt + input.txt (DELFIN) or *.inp (ORCA).'
    )


def _print_orca_summary(out_file: str, ok: bool, elapsed: float) -> None:
    """Print ORCA result summary to stdout for the DELFIN job log."""
    hours, rem = divmod(int(elapsed), 3600)
    mins, secs = divmod(rem, 60)
    status = 'SUCCESS' if ok else 'FAILED'
    print('')
    print(f'ORCA Result:   {status}')
    print(f'ORCA Elapsed:  {hours:02d}:{mins:02d}:{secs:02d}')
    try:
        out_path = Path(out_file)
        if out_path.is_file():
            size = out_path.stat().st_size
            print(f'Output Size:   {size:,} bytes')
            with open(out_path, 'r', errors='replace') as f:
                if size > 8192:
                    f.seek(size - 8192)
                tail = f.read()
                if 'ORCA TERMINATED NORMALLY' in tail:
                    print('ORCA Marker:   TERMINATED NORMALLY')
                else:
                    print('ORCA Marker:   NOT FOUND (output may be incomplete)')
        else:
            print(f'Output File:   NOT FOUND')

        # On failure: classify error and show relevant lines
        if not ok and out_path.is_file():
            _print_orca_error_diagnosis(out_path)
    except Exception:
        pass


# Keywords that flag an ORCA error line (case-insensitive check).
_ORCA_ERROR_KEYWORDS = (
    'error', 'aborting', 'not converged', 'segmentation fault',
    'signal:', 'mpirun noticed', 'cannot allocate', 'insufficient memory',
    'out of memory', 'no space left', 'disk quota',
)


def _print_orca_error_diagnosis(out_path: Path) -> None:
    """Classify the ORCA error and print relevant lines from the output."""
    # 1) Use the existing error detector for classification
    try:
        from delfin.orca_recovery import OrcaErrorDetector
        error_type = OrcaErrorDetector.analyze_output(out_path)
        if error_type is not None:
            print(f'Error Type:    {error_type.value}')
    except Exception:
        pass

    # 2) Detect the last ORCA computation phase from the output
    try:
        _print_last_orca_phase(out_path)
    except Exception:
        pass

    # 3) Extract the last ~50 error-relevant lines from the output tail
    try:
        with open(out_path, 'r', errors='replace') as f:
            size = out_path.stat().st_size
            # Read last 32 KB — enough to capture error context
            if size > 32768:
                f.seek(size - 32768)
                f.readline()  # skip partial first line
            lines = f.readlines()

        error_lines: list[str] = []
        for line in lines:
            low = line.lower()
            if any(kw in low for kw in _ORCA_ERROR_KEYWORDS):
                stripped = line.rstrip()
                if stripped and stripped not in error_lines:
                    error_lines.append(stripped)

        if error_lines:
            # Show up to 15 most relevant lines
            print('Error Lines:')
            for line in error_lines[-15:]:
                print(f'  | {line}')
        elif not any('ORCA TERMINATED NORMALLY' in l for l in lines[-20:]):
            # No error keywords but also no success marker → abrupt termination
            print('Error Lines:   (none found — output ends abruptly, likely killed by signal/timeout)')
            # Show last 5 non-empty lines for context
            tail_lines = [l.rstrip() for l in lines[-10:] if l.strip()]
            if tail_lines:
                print('Last output:')
                for line in tail_lines[-5:]:
                    print(f'  | {line}')
    except Exception:
        pass

    # 4) Disk space — helps spot "no space left" before reading the error
    try:
        stat = os.statvfs(out_path.parent)
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
        print(f'Disk Free:     {free_gb:.1f} GB ({out_path.parent})')
    except Exception:
        pass


# ORCA computation phases, ordered by typical execution sequence.
# The LAST match in the output tells us where ORCA was when it died.
_ORCA_PHASE_MARKERS = (
    ('SCF ITERATIONS', 'SCF iterations'),
    ('GEOMETRY OPTIMIZATION CYCLE', 'geometry optimization'),
    ('THE OPTIMIZATION HAS CONVERGED', 'optimization converged'),
    ('OPTIMIZATION RUN DONE', 'optimization done'),
    ('SCF HESSIAN', 'analytical Hessian'),
    ('COSX Hessian', 'COSX Hessian evaluation'),
    ('CP-SCF DRIVER', 'CPSCF / frequency perturbations'),
    ('VIBRATIONAL FREQUENCIES', 'frequency analysis'),
    ('NORMAL MODES', 'normal mode analysis'),
    ('IR SPECTRUM', 'IR spectrum'),
    ('THERMOCHEMISTRY', 'thermochemistry'),
)


def _print_last_orca_phase(out_path: Path) -> None:
    """Detect and print the last ORCA computation phase from output."""
    size = out_path.stat().st_size
    # Read last 64 KB to find the latest phase marker
    with open(out_path, 'r', errors='replace') as f:
        if size > 65536:
            f.seek(size - 65536)
            f.readline()
        content = f.read()

    last_phase = None
    for marker, label in _ORCA_PHASE_MARKERS:
        if marker in content:
            last_phase = label
    if last_phase:
        print(f'Last Phase:    {last_phase}')


def _run_mode(mode: str) -> int:
    inp_file = str(os.environ.get('DELFIN_INP_FILE') or '').strip()
    override = str(os.environ.get('DELFIN_OVERRIDE') or '').strip()
    xyz_file = str(os.environ.get('DELFIN_XYZ_FILE') or '').strip()
    workflow_label = str(os.environ.get('DELFIN_WORKFLOW_LABEL') or '').strip()
    hyperpol_use_bfw = str(os.environ.get('DELFIN_HYPERPOL_XTB_BFW') or '').strip().lower() in (
        '1', 'true', 'yes', 'on'
    )
    tadf_use_bfw = str(os.environ.get('DELFIN_TADF_XTB_BFW') or '').strip().lower() in (
        '1', 'true', 'yes', 'on'
    )
    build_mult = str(os.environ.get('BUILD_MULTIPLICITY') or '1')
    delfin_pal = str(os.environ.get('DELFIN_PAL') or '4')
    delfin_maxcore = str(os.environ.get('DELFIN_MAXCORE') or '1000')
    delfin_cmd = [sys.executable, '-m', 'delfin']
    resolved_mode = _detect_mode(mode)
    orca_bin = _resolve_orca_bin()

    print(f'ORCA: {orca_bin or "not found"}')
    print(f'Python: {sys.executable}')
    print(f'DELFIN: {" ".join(shlex.quote(part) for part in delfin_cmd)}')
    print('')

    if resolved_mode == 'delfin':
        print('Starting DELFIN...')
        return _run_command(delfin_cmd)
    if resolved_mode == 'delfin-recalc':
        print('Starting DELFIN --recalc...')
        os.environ.setdefault('DELFIN_SMART_RECALC', '1')
        return _run_command(delfin_cmd + ['--recalc'])
    if resolved_mode == 'delfin-recalc-classic':
        print('Starting DELFIN --recalc (classic marker mode)...')
        os.environ['DELFIN_SMART_RECALC'] = '0'
        return _run_command(delfin_cmd + ['--recalc'])
    if resolved_mode == 'delfin-recalc-override':
        if not override:
            print('ERROR: DELFIN_OVERRIDE not set for delfin-recalc-override mode')
            return 1
        print(f'Starting DELFIN --recalc --occupier-override {override}...')
        return _run_command(delfin_cmd + ['--recalc', '--occupier-override', override])
    if resolved_mode == 'orca':
        from delfin.orca import run_orca

        if not inp_file:
            inp_candidates = sorted(Path.cwd().glob('*.inp'))
            inp_file = inp_candidates[0].name if inp_candidates else ''
        if not orca_bin:
            print('ERROR: ORCA executable not found.')
            print('       Set DELFIN_ORCA_BASE to a valid ORCA installation path or add ORCA to PATH.')
            return 1
        if not inp_file:
            print('ERROR: No .inp file found for ORCA mode')
            return 1
        if not Path(inp_file).is_file():
            print(f'ERROR: Input file not found: {inp_file}')
            return 1
        out_file = f'{Path(inp_file).stem}.out'
        print(f'Starting ORCA: {inp_file} -> {out_file}')
        t0 = time.monotonic()
        ok = run_orca(inp_file, out_file)
        elapsed = time.monotonic() - t0
        _print_orca_summary(out_file, ok, elapsed)
        return 0 if ok else 1
    if resolved_mode == 'build':
        print('Starting delfin-build (complex build-up)...')
        print(f'  Multiplicity: {build_mult}')
        return _run_command([
            sys.executable,
            '-m',
            'delfin.build_up_complex',
            'input.txt',
            '--goat',
            '--directory',
            'builder',
            '--multiplicity',
            build_mult,
            '--verbose',
        ])
    if resolved_mode == 'guppy':
        guppy_runs = str(os.environ.get('GUPPY_RUNS') or '20')
        guppy_pal = str(os.environ.get('GUPPY_PAL') or os.environ.get('DELFIN_PAL') or '12')
        guppy_maxcore = str(os.environ.get('GUPPY_MAXCORE') or os.environ.get('DELFIN_MAXCORE') or '500')
        guppy_parallel_jobs = str(os.environ.get('GUPPY_PARALLEL_JOBS') or '4')
        guppy_goat_topk = str(os.environ.get('GUPPY_GOAT_TOPK') or '0')
        guppy_goat_parallel = str(os.environ.get('GUPPY_GOAT_PARALLEL_JOBS') or guppy_parallel_jobs)
        print('Starting GUPPY SMILES sampling...')
        return _run_command([
            sys.executable,
            '-m',
            'delfin.guppy_sampling',
            'input.txt',
            '--runs',
            guppy_runs,
            '--pal',
            guppy_pal,
            '--maxcore',
            guppy_maxcore,
            '--parallel-jobs',
            guppy_parallel_jobs,
            '--goat-topk',
            guppy_goat_topk,
            '--goat-parallel-jobs',
            guppy_goat_parallel,
            '--output',
            'GUPPY_try.xyz',
        ])
    if resolved_mode in {'hyperpol_xtb', 'tadf_xtb', 'censo_anmr'}:
        if not xyz_file:
            print(f'ERROR: DELFIN_XYZ_FILE not set for {resolved_mode} mode')
            return 1
        if not Path(xyz_file).is_file():
            print(f'ERROR: XYZ file not found for {resolved_mode} mode: {xyz_file}')
            return 1
        label = workflow_label or resolved_mode
        summary_name = f'{resolved_mode}_summary.json'
        output_name = f'{resolved_mode}.output'
        print(f'Starting browser {resolved_mode} workflow...')
        cmd = [
            sys.executable,
            '-m',
            'delfin.dashboard.browser_workflows',
            resolved_mode,
            '--xyz-file',
            xyz_file,
            '--label',
            label,
            '--workdir',
            str(Path.cwd()),
            '--pal',
            delfin_pal,
            '--maxcore',
            delfin_maxcore,
            '--json-out',
            str(Path.cwd() / summary_name),
        ]
        if resolved_mode == 'hyperpol_xtb':
            cmd.extend(['--engine', 'std2', '--preopt', 'none', '--static-only', '--energy-window', '15'])
            if hyperpol_use_bfw:
                cmd.append('--bfw')
        elif resolved_mode == 'tadf_xtb':
            tadf_preopt = str(os.environ.get('DELFIN_TADF_XTB_PREOPT') or 'xtb').strip().lower() or 'xtb'
            tadf_excited_method = (
                str(os.environ.get('DELFIN_TADF_XTB_EXCITED_METHOD') or 'stda').strip().lower() or 'stda'
            )
            tadf_energy_window = str(os.environ.get('DELFIN_TADF_XTB_ENERGY_WINDOW') or '10').strip() or '10'
            tadf_t1_opt = str(os.environ.get('DELFIN_TADF_XTB_T1_OPT') or 'yes').strip().lower()
            cmd.extend([
                '--preopt',
                tadf_preopt,
                '--excited-method',
                tadf_excited_method,
                '--energy-window',
                tadf_energy_window,
            ])
            if tadf_t1_opt in {'1', 'true', 'yes', 'on'}:
                cmd.append('--t1-opt')
            if tadf_use_bfw:
                cmd.append('--bfw')
        elif resolved_mode == 'censo_anmr':
            censo_solvent = str(os.environ.get('DELFIN_CENSO_NMR_SOLVENT') or 'chcl3').strip().lower() or 'chcl3'
            censo_charge = str(os.environ.get('DELFIN_CENSO_NMR_CHARGE') or '0').strip() or '0'
            censo_multiplicity = str(os.environ.get('DELFIN_CENSO_NMR_MULTIPLICITY') or '1').strip() or '1'
            censo_mhz = str(os.environ.get('DELFIN_CENSO_NMR_MHZ') or '400').strip() or '400'
            censo_resume = str(os.environ.get('DELFIN_CENSO_NMR_RESUME') or '').strip().lower()
            cmd.extend([
                '--solvent',
                censo_solvent,
                '--charge',
                censo_charge,
                '--multiplicity',
                censo_multiplicity,
                '--mhz',
                censo_mhz,
            ])
            if censo_resume in {'1', 'true', 'yes', 'on'}:
                cmd.append('--resume')
        return _run_and_tee(cmd, Path.cwd() / output_name)
    if resolved_mode == 'delfin-co2-chain':
        species_delta = str(os.environ.get('DELFIN_CO2_SPECIES_DELTA') or '0')
        print('Starting DELFIN + CO2 Coordinator chain...')
        exit_code = _run_command(delfin_cmd)
        if exit_code != 0:
            print(f'ERROR: DELFIN failed (exit {exit_code}). CO2 Coordinator skipped.')
            return exit_code
        setup_code = _run_command([sys.executable, '-m', 'delfin.co2.chain_setup', species_delta])
        if setup_code != 0:
            print(f'ERROR: CO2 setup failed (exit {setup_code})')
            return setup_code
        return _run_command(delfin_cmd + ['co2'], cwd=Path.cwd() / 'CO2_coordination')

    print(f'ERROR: Unknown mode: {resolved_mode}')
    return 1


def main(argv: list[str] | None = None) -> int:
    del argv
    signal.signal(signal.SIGTERM, _handle_termination)
    signal.signal(signal.SIGINT, _handle_termination)
    signal.signal(signal.SIGHUP, _handle_termination)
    _configure_environment()

    mode = str(os.environ.get('DELFIN_MODE') or 'auto').strip() or 'auto'
    job_name = str(os.environ.get('DELFIN_JOB_NAME') or 'local_job').strip() or 'local_job'
    inp_file = str(os.environ.get('DELFIN_INP_FILE') or '').strip()
    _print_job_banner(mode, job_name, inp_file)

    exit_code = 1
    try:
        exit_code = _run_mode(mode)
        return int(exit_code)
    except SystemExit as exc:
        exit_code = exc.code if isinstance(exc.code, int) else 1
        raise
    finally:
        _write_exit_code(int(exit_code))


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        # Catch-all so import errors or unexpected crashes are ALWAYS visible
        # in the SLURM delfin_*.out / delfin_*.err files.
        import traceback
        print(f'[local_runner] FATAL: {exc}', file=sys.stderr, flush=True)
        traceback.print_exc()
        raise SystemExit(1)
