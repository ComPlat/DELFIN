"""Input processing helpers: SMILES wrappers, resource parsing, sanitisation."""

import json
import os
import re
import subprocess
import sys

from delfin.smiles_converter import (
    smiles_to_xyz as _delfin_smiles_to_xyz,
    smiles_to_xyz_isomers as _delfin_smiles_to_xyz_isomers,
    smiles_to_xyz_quick as _delfin_smiles_to_xyz_quick,
    is_smiles_string as _delfin_is_smiles_string,
    contains_metal,
)


# Subprocess-isolation of isomers calls: each conversion runs in a
# fresh Python process so the Voila / Jupyter kernel never accumulates
# RDKit / OpenBabel memory across calls.  Enabled by default; set
# ``DELFIN_UI_INLINE=1`` to bypass (e.g. when profiling or for
# cooperative debugging).
#: Fallback only.  The live answer comes from :func:`_isolation_wanted`, which
#: reads the environment at call time -- see the note on ``_UI_ISOLATE_TIMEOUT``
#: below for why a module-level read is not good enough.
_UI_ISOLATE_DEFAULT = os.environ.get("DELFIN_UI_INLINE", "0") != "1"
# Per-conversion subprocess timeout (s).  Set DELFIN_UI_ISOLATE_TIMEOUT=0 (or negative) for NO
# timeout, so a large complete manifold (huge, heavily-substituted macrocycles) can run to completion
# instead of being killed — the construction is deterministic, not a hang.
#
# Raised from 1800 to 3600.  A cut build returns *nothing*, not a smaller
# answer, and 1800 s was cutting the band DELFIN users actually submit.
# Measured over 2000 systems built at extreme with max_isomers=0, by build time
# against total atom count:
#
#     atoms     n     mean    >1800s   >3600s   >7200s
#     <=30    442      35 s     0.0%     0.0%     0.0%
#     31-50   379     333 s     2.9%     1.1%     0.0%
#     51-80   757     824 s     9.6%     2.5%     0.0%
#     >80     422    1802 s    39.6%    15.6%     0.0%
#
# A cyclam complex is ~40 atoms, a bis-terpyridine ~60-75, a substituted
# porphyrin 50-90.  At 1800 s two of those bands lost 9.6 % and 39.6 % of their
# builds to the clock and returned an empty answer.  Nothing in that sample
# exceeded 7200 s, so the tail is bounded, not infinite.
#: The default, not the answer.  Read at call time by
#: :func:`_resolve_isolate_timeout`, because a module-level read freezes the
#: value at first import: CONTROL sets DELFIN_UI_ISOLATE_TIMEOUT from
#: MANTA_TIME_BUDGET, and if anything imports this module before that runs, a
#: user who asked for 7200 s is killed at the default and nothing says so.
#: That is the same failure as a provenance record written from a second
#: environment lookup -- it does not raise, it just quietly does something else.
_UI_ISOLATE_TIMEOUT = int(os.environ.get("DELFIN_UI_ISOLATE_TIMEOUT", "3600"))


def _resolve_isolate_timeout(explicit=None):
    """The build timeout in seconds, or ``None`` for no limit.

    ``0`` means no limit and has to survive as such: a complete manifold on a
    heavy macrocycle is a long deterministic construction, not a hang.  An
    explicit ``0`` from a caller used to fall through to the default, because
    ``timeout or _UI_ISOLATE_TIMEOUT`` cannot tell ``0`` from unset.
    """
    if explicit is None:
        raw = os.environ.get("DELFIN_UI_ISOLATE_TIMEOUT", _UI_ISOLATE_TIMEOUT)
    else:
        raw = explicit
    try:
        seconds = int(float(raw))
    except (TypeError, ValueError):
        seconds = _UI_ISOLATE_TIMEOUT
    return None if seconds <= 0 else seconds


def _isolation_wanted():
    """Whether the build runs in its own subprocess, read at call time.

    ``DELFIN_UI_INLINE=1`` turns the isolation off, and with it the only place
    a build timeout can be enforced -- the budget then has no effect at all,
    not even a wrong one.  Inheriting that variable from a parent process is
    enough to do it, so a run that asked for a budget and cannot get one is
    told rather than left to assume.
    """
    inline = os.environ.get("DELFIN_UI_INLINE", "0") == "1"
    if inline and os.environ.get("DELFIN_UI_ISOLATE_TIMEOUT"):
        print("WARNING: DELFIN_UI_INLINE=1 runs the structure build in this "
              "process, so MANTA_TIME_BUDGET / DELFIN_UI_ISOLATE_TIMEOUT "
              "cannot be enforced and a long build will not be stopped.",
              file=sys.stderr)
    return not inline


def _kill_process_group(proc):
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


def _run_isomers_subprocess(smiles, kwargs, timeout=None):
    """Run ``smiles_to_xyz_isomers`` in a fresh subprocess.

    Returns ``(results, error)`` where ``results`` is a list of
    ``(xyz_string, label)`` tuples and ``error`` is an optional string.
    Subprocess exits after returning -> all RDKit / OB memory is
    released back to the OS, which keeps the Voila kernel's RSS
    bounded across many conversions.

    Input ``kwargs`` are JSON-serialised, so only primitives / lists
    are accepted.  Output XYZ strings + labels are also JSON-safe.
    """
    timeout = _resolve_isolate_timeout(timeout)
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
        # Fixed PYTHONHASHSEED in the child so set / dict iteration order
        # is stable across invocations.  Without this, the same SMILES can
        # produce wildly different isomer counts (e.g. 1 vs 9 vs 30) between
        # runs because the enumerator's candidate-retention order depends
        # on Python hash randomisation.
        env = dict(os.environ)
        env.setdefault("PYTHONHASHSEED", "0")
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
            env=env,
            start_new_session=True,
        )
        try:
            out, err = proc.communicate(input=payload, timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_process_group(proc)
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
try:
    from delfin.smiles_converter import (
        smiles_to_xyz_quick_hapto_previews as _delfin_smiles_to_xyz_quick_hapto_previews,
    )
except ImportError:
    _delfin_smiles_to_xyz_quick_hapto_previews = None

_HAPTO_FAILFAST_TOKEN = "Hapto (eta) coordination detected"


def _is_hapto_failfast(error: str) -> bool:
    return bool(error) and _HAPTO_FAILFAST_TOKEN in error


def smiles_to_xyz(smiles, apply_uff=True, hapto_approx=None):
    """Convert a SMILES string to XYZ coordinates.

    Returns ``(xyz_string, num_atoms, method, error)``.
    """
    xyz_string, error = _delfin_smiles_to_xyz(
        smiles, apply_uff=apply_uff, hapto_approx=hapto_approx
    )
    if error and hapto_approx is None and _is_hapto_failfast(error):
        xyz_string, error = _delfin_smiles_to_xyz(
            smiles, apply_uff=apply_uff, hapto_approx=True
        )
    if error:
        return None, 0, None, error
    num_atoms = sum(1 for line in xyz_string.splitlines() if line.strip())
    method = 'delfin.smiles_converter'
    return xyz_string, num_atoms, method, None


def smiles_to_xyz_quick(smiles, hapto_approx=None):
    """Fast single-conformer conversion, no OB, no multi-seed, no UFF.

    Returns ``(xyz_string, num_atoms, method, error)``.
    """
    xyz_string, error = _delfin_smiles_to_xyz_quick(
        smiles, hapto_approx=hapto_approx
    )
    if error and hapto_approx is None and _is_hapto_failfast(error):
        xyz_string, error = _delfin_smiles_to_xyz_quick(
            smiles, hapto_approx=True
        )
    if error:
        return None, 0, None, error
    num_atoms = sum(1 for line in xyz_string.splitlines() if line.strip())
    return xyz_string, num_atoms, 'quick', None


def smiles_to_xyz_quick_with_previews(smiles, hapto_approx=None):
    """Fast single-conformer conversion plus hapto-specific preview structures."""
    xyz_string, num_atoms, method, error = smiles_to_xyz_quick(
        smiles,
        hapto_approx=hapto_approx,
    )
    if error or not xyz_string:
        return xyz_string, num_atoms, method, [], error

    if _delfin_smiles_to_xyz_quick_hapto_previews is None:
        return xyz_string, num_atoms, method, [], None
    previews = _delfin_smiles_to_xyz_quick_hapto_previews(
        smiles,
        hapto_approx=hapto_approx,
    )
    preview_items = []
    seen_keys = {
        "\n".join(line.strip() for line in xyz_string.splitlines() if line.strip())
    }
    for preview_xyz, label in previews:
        key = "\n".join(line.strip() for line in preview_xyz.splitlines() if line.strip())
        if not key or key in seen_keys:
            continue
        seen_keys.add(key)
        preview_num_atoms = sum(1 for line in preview_xyz.splitlines() if line.strip())
        preview_items.append((preview_xyz, preview_num_atoms, label))
    return xyz_string, num_atoms, method, preview_items, None


def append_hapto_previews_to_isomers(
    isomers,
    smiles,
    *,
    include_quick=False,
    hapto_approx=None,
):
    """Append cached hapto preview structures to an isomer list without duplicates."""
    merged = list(isomers)
    seen_keys = {
        "\n".join(line.strip() for line in xyz_string.splitlines() if line.strip())
        for xyz_string, _num_atoms, _label in merged
    }

    xyz_string, num_atoms, _method, preview_items, error = smiles_to_xyz_quick_with_previews(
        smiles,
        hapto_approx=hapto_approx,
    )
    if error or not xyz_string:
        return merged

    extra_items = list(preview_items)
    if include_quick:
        extra_items.insert(0, (xyz_string, num_atoms, 'quick'))

    for preview_xyz, preview_num_atoms, label in extra_items:
        key = "\n".join(line.strip() for line in preview_xyz.splitlines() if line.strip())
        if not key or key in seen_keys:
            continue
        seen_keys.add(key)
        merged.append((preview_xyz, preview_num_atoms, label))
    return merged


def smiles_to_xyz_isomers(
    smiles,
    apply_uff=True,
    collapse_label_variants=True,
    include_binding_mode_isomers=False,
    hapto_approx=None,
    deterministic=True,
    quality_mode="extreme",
    seeds_override=None,
    n_metal_smart=True,
    max_isomers=None,
):
    """Generate distinct coordination isomers for a SMILES string.

    Returns ``([(xyz_string, num_atoms, label), ...], error)``.

    ``quality_mode`` defaults to ``"extreme"`` (60 ETKDG seeds, 5 chelate
    ranks, 5 templates, 12 alt-binding tries).  The pipeline honours
    ``DELFIN_MAX_PROCESS_WORKERS`` / ``_THREAD_WORKERS`` (default 64)
    so CPU / RAM pressure stays bounded on smaller machines.  Pass
    ``"max"`` / ``"normal"`` / ``"fast"`` for progressively cheaper
    candidate pool.

    ``seeds_override`` (int, optional) pins the seed count independently
    of the quality profile — used by the dashboard's custom slider.
    """
    base_kwargs = dict(
        apply_uff=apply_uff,
        collapse_label_variants=collapse_label_variants,
        include_binding_mode_isomers=include_binding_mode_isomers,
        hapto_approx=hapto_approx,
        deterministic=deterministic,
        quality_mode=quality_mode,
        seeds_override=seeds_override,
        n_metal_smart=n_metal_smart,
    )
    # max_isomers: None -> library default (byte-identical); set -> forwarded so the
    # MANTA button can request the COMPLETE manifold (never cut off).
    if max_isomers is not None:
        base_kwargs["max_isomers"] = int(max_isomers)
    if _isolation_wanted():
        results, error = _run_isomers_subprocess(smiles, base_kwargs)
        if error and hapto_approx is None and _is_hapto_failfast(error):
            retry_kwargs = dict(base_kwargs, hapto_approx=True)
            results, error = _run_isomers_subprocess(smiles, retry_kwargs)
    else:
        results, error = _delfin_smiles_to_xyz_isomers(smiles, **base_kwargs)
        if error and hapto_approx is None and _is_hapto_failfast(error):
            retry_kwargs = dict(base_kwargs, hapto_approx=True)
            results, error = _delfin_smiles_to_xyz_isomers(smiles, **retry_kwargs)
    if error:
        return [], error
    out = []
    for xyz_string, label in results:
        num_atoms = sum(1 for line in xyz_string.splitlines() if line.strip())
        out.append((xyz_string, num_atoms, label))
    return out, None


def is_smiles(text):
    """Return *True* if *text* looks like a SMILES string."""
    try:
        return bool(_delfin_is_smiles_string(text))
    except Exception:
        return False


def clean_input_data(input_text):
    """Classify and clean raw input.

    Returns ``(cleaned_text, input_type)`` where *input_type* is one of
    ``'smiles'``, ``'xyz'``, or ``'empty'``.
    """
    text = input_text.strip()
    if not text:
        return '', 'empty'

    if is_smiles(text):
        return text, 'smiles'

    lines = text.split('\n')
    if len(lines) < 2:
        return text, 'xyz'

    first_line = lines[0].strip()
    try:
        int(first_line)
        cleaned_lines = lines[2:]
        return '\n'.join(cleaned_lines).strip(), 'xyz'
    except ValueError:
        return text, 'xyz'


def parse_resource_settings(control_text):
    """Parse PAL and maxcore from CONTROL.txt content.

    Returns ``(pal, maxcore)`` as ints or *None* if not found.
    """
    pal_match = re.search(r'^\s*PAL\s*=\s*(\d+)', control_text, flags=re.MULTILINE)
    maxcore_match = re.search(r'^\s*maxcore\s*=\s*(\d+)', control_text, flags=re.MULTILINE)
    pal = int(pal_match.group(1)) if pal_match else None
    maxcore = int(maxcore_match.group(1)) if maxcore_match else None
    return pal, maxcore


_PAL_NPROCS_RE = re.compile(r'(?i)\bnprocs\s*=?\s*(\d+)')
_PAL_KEYWORD_RE = re.compile(r'(?im)^\s*!.*?\bPAL\s*(\d+)\b')
_MAXCORE_RE = re.compile(r'(?im)^\s*%maxcore\s*=?\s*(\d+)')
_PAL_BLOCK_RE = re.compile(r'(?is)%pal\b.*?\bend\b')


def parse_inp_resources(inp_text):
    """Parse PAL (nprocs) and maxcore from ORCA ``.inp`` text.

    Accepts every form ORCA itself accepts:
      * ``%pal nprocs N end``                    (inline)
      * ``%pal\n  nprocs N\nend``                (multi-line)
      * ``%pal\n  nprocs=N\nend``                (= syntax)
      * ``! PAL<N>``                              (keyword shortcut, e.g. ``! PAL8``)
      * ``%maxcore N`` / ``%MaxCore=N``           (case-insensitive, optional ``=``)

    Returns ``(pal, maxcore)`` as ints or *None* if not found.
    """
    pal = None
    maxcore = None
    if not inp_text:
        return pal, maxcore
    m = _PAL_NPROCS_RE.search(inp_text)
    if m:
        pal = int(m.group(1))
    else:
        m = _PAL_KEYWORD_RE.search(inp_text)
        if m:
            pal = int(m.group(1))
    m = _MAXCORE_RE.search(inp_text)
    if m:
        maxcore = int(m.group(1))
    return pal, maxcore


def sanitize_orca_input(text):
    """Sanitize ORCA input to avoid hidden/invalid characters."""
    if text is None:
        return ''
    text = text.replace('\r\n', '\n').replace('\r', '\n').lstrip('\ufeff')
    text = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f\x7f-\x9f]', '', text)
    text = ''.join(ch for ch in text if ch == '\n' or ch == '\t' or (' ' <= ch <= '~'))
    lines = text.split('\n')
    out_lines = []
    for line in lines:
        if re.search(r'^\s*\*\s*xyzfile\b', line, flags=re.IGNORECASE):
            parts = line.split()
            if len(parts) >= 5:
                filename = parts[4].strip("\"'")
                filename = re.sub(r"[^A-Za-z0-9._/+-]", '', filename)
                m = re.match(r'(.+?\.xyz)', filename, flags=re.IGNORECASE)
                if m:
                    filename = m.group(1)
                parts = parts[:4] + [filename] + parts[5:]
                line = ' '.join(parts)
        out_lines.append(line)
    return '\n'.join(out_lines).strip() + '\n'
