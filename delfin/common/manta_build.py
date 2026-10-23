"""One way to build a MANTA manifold, shared by every entry point.

``delfin-manta`` (cli_manta), the dashboard's MANTA button (structure_editor via
input_processing) and the CONTROL pipeline (guppy_sampling) all build through
:func:`build_isomers_isolated`, so they cannot drift apart again in how the build
is run:

* the construction environment (``cli_manta.construction_env``) is handed to the
  build's own subprocess only -- the calling process (a Voila kernel, a pipeline)
  is never mutated, so a concurrent build elsewhere in it cannot pick it up;
* the subprocess always runs with ``PYTHONHASHSEED=0``, whatever the shell or the
  kernel has, so set / dict iteration order is the same for every entry point;
* the hapto fail-fast retry is decided once (:func:`hapto_retry_wanted`);
* the opt-in post-processing -- single-point energy ranking and top-N geometry
  optimisation -- is the same two functions for the CLI and the dashboard.

Deliberately free of heavy imports: the CLI imports it before anything that
reads the construction switches at import time.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

#: PYTHONHASHSEED of every MANTA build subprocess.  Fixed, and not taken from the
#: caller's environment: the dashboard pinned 0 while a shell running delfin-manta
#: had whatever the user had, which made the entry points differ by construction.
HASH_SEED = "0"

#: Opt-in geometry optimisation: default top-N and the parallel worker count
#: (laptop-bounded; tune here).
OPT_TOPN = 10
OPT_WORKERS = 4

HAPTO_FAILFAST_TOKEN = "Hapto (eta) coordination detected"
_FALSE = {"0", "false", "no", "off"}


def kill_process_group(proc):
    """Kill a timed-out build and everything it spawned.

    The builder's own children are the point: it opens a ProcessPoolExecutor
    for batch UFF, so a bare ``proc.kill()`` can leave dozens of workers
    reparented to init.  SIGTERM to the group first so anything with a handler
    can tidy up, SIGKILL after a short grace period, then reap.
    """
    import signal
    import time

    try:
        pgid = os.getpgid(proc.pid)
    except (OSError, ProcessLookupError):
        pgid = None

    for sig, wait in ((signal.SIGTERM, 2.0), (signal.SIGKILL, 1.0)):
        try:
            if pgid is not None:
                os.killpg(pgid, sig)
            else:
                proc.send_signal(sig)
        except (OSError, ProcessLookupError):
            break
        deadline = time.monotonic() + wait
        while time.monotonic() < deadline:
            if proc.poll() is not None:
                break
            time.sleep(0.05)
        if proc.poll() is not None:
            break
    try:
        proc.communicate(timeout=1.0)
    except Exception:                              # noqa: BLE001
        pass


def hapto_retry_wanted(error, hapto_approx, environ=None) -> bool:
    """Whether a hapto fail-fast answer is retried with ``hapto_approx=True``.

    Only when the caller left hapto on auto (``None``) AND the user did not ask
    for fail-fast explicitly with ``DELFIN_HAPTO_APPROX=0`` (or false/no/off) --
    retrying then would overrule exactly the setting that produced the error.
    """
    if hapto_approx is not None or not error or HAPTO_FAILFAST_TOKEN not in str(error):
        return False
    environ = os.environ if environ is None else environ
    return str(environ.get("DELFIN_HAPTO_APPROX", "")).strip().lower() not in _FALSE


def run_isomers_isolated(smiles, kwargs, *, timeout=None, env=None):
    """Run ``smiles_to_xyz_isomers`` in a fresh subprocess.

    Returns ``(results, error)`` where ``results`` is a list of
    ``(xyz_string, label)`` tuples and ``error`` is an optional string.
    Subprocess exits after returning -> all RDKit / OB memory is
    released back to the OS, which keeps the Voila kernel's RSS
    bounded across many conversions.

    Input ``kwargs`` are JSON-serialised, so only primitives / lists
    are accepted.  Output XYZ strings + labels are also JSON-safe.
    """
    # timeout: seconds, or None for no limit (the caller resolves its own budget).
    payload = json.dumps({"smiles": smiles, "kwargs": kwargs})
    script = (
        "import json, sys\n"
        "from delfin.smiles_converter import smiles_to_xyz_isomers\n"
        "req = json.loads(sys.stdin.read())\n"
        "res, err = smiles_to_xyz_isomers(req['smiles'], **req['kwargs'])\n"
        "out = {'r': [list(t) for t in (res or [])], 'e': err}\n"
        "sys.stdout.write('__DELFIN_RESULT__' + json.dumps(out))\n"
    )
    try:
        # The caller's environment, plus this build's construction switches (env),
        # plus the fixed hash seed.  Only the child sees them.  PYTHONHASHSEED is
        # forced, not defaulted: without a fixed value the enumerator's
        # candidate-retention order can follow hash randomisation (once measured:
        # 1 vs 9 vs 30 isomers between runs), and the dashboard and the CLI would
        # differ by whatever the kernel or the shell happened to carry.
        child_env = dict(os.environ)
        child_env.update(env or {})
        child_env["PYTHONHASHSEED"] = HASH_SEED
        # ``start_new_session`` puts the child in its own process group so the
        # timeout can kill the group.  ``subprocess.run(timeout=)`` kills only
        # the direct child, and this child is not a leaf: the builder opens a
        # ProcessPoolExecutor for batch UFF with up to DELFIN_MAX_PROCESS_WORKERS
        # (64) workers.  Killing the parent alone leaves those reparented to
        # init, holding RAM until somebody notices -- measured elsewhere in this
        # tree as 128 orphans alive for 3-5 hours after one such kill.
        proc = subprocess.Popen(
            [sys.executable, "-c", script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=child_env,
            start_new_session=True,
        )
        try:
            out, err = proc.communicate(input=payload, timeout=timeout)
        except subprocess.TimeoutExpired:
            kill_process_group(proc)
            return [], f"Subprocess timed out after {timeout}s"
        proc = subprocess.CompletedProcess(
            proc.args, proc.returncode, stdout=out, stderr=err)
    except subprocess.TimeoutExpired:
        return [], f"Subprocess timed out after {timeout}s"
    if proc.returncode != 0:
        tail = (proc.stderr or "").splitlines()[-10:]
        return [], f"Subprocess failed (exit {proc.returncode}): {' | '.join(tail)}"
    marker = "__DELFIN_RESULT__"
    for line in proc.stdout.splitlines():
        idx = line.find(marker)
        if idx >= 0:
            j = json.loads(line[idx + len(marker):])
            results = [tuple(x) for x in j.get("r", [])]
            return results, j.get("e")
    return [], "Subprocess returned no result marker"


def build_isomers_isolated(smiles, kwargs, *, timeout=None, env=None):
    """Build the manifold in a fresh subprocess, with the shared hapto retry.

    ``kwargs`` go to ``smiles_to_xyz_isomers`` (JSON primitives only); ``env`` is
    the construction environment for this build (``cli_manta.construction_env``),
    applied to the subprocess only; ``timeout`` in seconds, ``None`` = no limit.
    Returns ``([(xyz, label), ...], error)``.
    """
    results, error = run_isomers_isolated(smiles, kwargs, timeout=timeout, env=env)
    merged = dict(os.environ)
    merged.update(env or {})
    if hapto_retry_wanted(error, kwargs.get("hapto_approx"), merged):
        results, error = run_isomers_isolated(
            smiles, dict(kwargs, hapto_approx=True), timeout=timeout, env=env)
    return results, error


def rank_by_single_point(isomers, charge, method="gfn2", spin="auto"):
    """RANK the manifold by xtb SINGLE-POINT energy: reorder best (lowest-energy) first WITHOUT
    changing any geometry.  Each item is ``(xyz_string, num_atoms, label)``; the emitted structures
    stay byte-identical to construction — only their ORDER changes.  spin='auto' -> parity-correct
    uhf per structure (even electrons=singlet, odd=doublet); a fixed multiplicity sets uhf=mult-1.
    Best-effort: any structure whose energy eval fails sinks to the end keeping its geometry.
    Returns the list unchanged if xtb is unavailable or there is nothing to reorder."""
    if not isomers or len(isomers) < 2:
        return isomers
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    import concurrent.futures as _cf

    def _uhf_for(xyz):
        if str(spin) != "auto":
            return max(0, int(spin) - 1)
        try:
            return _gff._n_electrons(xyz, int(charge)) % 2   # parity-correct ground-state multiplicity
        except Exception:
            return 0

    def _energy_one(item):
        xyz = item[0]
        try:
            return _gff.gfnff_energy(xyz, charge=int(charge), uhf=_uhf_for(xyz), method=method)
        except Exception:
            return None
    _max_workers = max(1, min(len(isomers), (os.cpu_count() or 4)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=_max_workers) as ex:
            energies = list(ex.map(_energy_one, isomers))
    except Exception:
        return isomers
    # Ascending by energy; failed evals (None) sink to the end preserving their relative order.
    order = sorted(range(len(isomers)),
                   key=lambda i: (energies[i] is None, energies[i] if energies[i] is not None else 0.0, i))
    return [isomers[i] for i in order]

def optimise_top(isomers, charge, topn=None, method="gfn2", spin="auto"):
    """Geometry-optimize the top-N ranked isomers in parallel (laptop-bounded),
    replace their geometry + label, re-sort the optimized head by opt energy. The
    opt ``method`` FOLLOWS the Rank selection (gfn2/gfnff/gfn1/gfn0) so one switch
    controls both. Each item is ``(xyz_string, num_atoms, label)``. Best-effort:
    any structure whose optimization fails keeps its unrelaxed geometry.
    ``topn`` (user-settable): None -> OPT_TOPN; 0 -> ALL structures (optimise
    the complete ranked manifold, slowest/best); N>0 -> top-N; N<0 -> none."""
    if not isomers:
        return isomers
    if topn is None:
        _n = OPT_TOPN
    elif int(topn) == 0:
        _n = len(isomers)                  # 0 = ALL (optimise everything)
    elif int(topn) < 0:
        return isomers                     # negative = none
    else:
        _n = int(topn)
    import concurrent.futures as _cf
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    head = list(isomers[:_n])
    tail = list(isomers[_n:])

    def _opt_one(item):
        xyz, _na, label = item
        try:
            if str(spin) == "auto":
                # auto-spin: scan multiplicity -> GFN2 ground state (parity-correct)
                r = _gff.gfnff_optimize_autospin(xyz, charge=int(charge), method=method)
            else:
                # fixed multiplicity chosen by the user: uhf = multiplicity - 1
                _uhf = max(0, int(spin) - 1)
                r = _gff.gfnff_optimize(xyz, charge=int(charge), uhf=_uhf, method=method)
        except Exception:
            r = None
        if r and r[0]:
            opt_xyz, e = r
            na = len([ln for ln in opt_xyz.splitlines() if ln.strip()])
            _m = method.upper()
            tag = (" [%s-opt %.1f kcal]" % (_m, e)) if e is not None else " [%s-opt]" % _m
            return ((opt_xyz, na, (label or "isomer") + tag), e)
        return (item, None)

    workers = max(1, min(OPT_WORKERS, len(head)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=workers) as ex:
            opted = list(ex.map(_opt_one, head))
    except Exception:
        return isomers
    # optimized structures first, sorted by GFN2-opt energy (failed/None last)
    opted.sort(key=lambda t: (t[1] is None, t[1] if t[1] is not None else 0.0))
    return [it for (it, _e) in opted] + tail

