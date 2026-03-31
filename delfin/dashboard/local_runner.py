"""Built-in local job runner used when ``run_local.sh`` is unavailable."""

from __future__ import annotations

import os
import shlex
import shutil
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def _job_id() -> str:
    return str(os.environ.get('DELFIN_JOB_ID', '0') or '0')


def _write_exit_code(code: int) -> None:
    (Path.cwd() / f'.exit_code_{_job_id()}').write_text(f'{int(code)}\n', encoding='utf-8')


def _handle_termination(signum, _frame):
    print(f'Received signal {signum}; stopping local DELFIN job.')
    _write_exit_code(124)
    raise SystemExit(124)


def _configure_environment() -> None:
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
    print('========================================')
    print('')


def _run_command(cmd: list[str], *, cwd: Path | None = None) -> int:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=False).returncode


def _run_and_tee(cmd: list[str], output_path: Path, *, cwd: Path | None = None) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as handle:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            handle.write(line)
        return proc.wait()


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
        with Path(out_file).open('w', encoding='utf-8') as handle:
            proc = subprocess.run([orca_bin, inp_file], stdout=handle, stderr=subprocess.STDOUT, check=False)
        return proc.returncode
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
    if resolved_mode in {'hyperpol_xtb', 'tadf_xtb'}:
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
    finally:
        _write_exit_code(int(exit_code))


if __name__ == '__main__':
    raise SystemExit(main())
