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
_UI_ISOLATE_DEFAULT = os.environ.get("DELFIN_UI_INLINE", "0") != "1"
# Per-conversion subprocess timeout (s).  Set DELFIN_UI_ISOLATE_TIMEOUT=0 (or negative) for NO
# timeout, so a large complete manifold (huge, heavily-substituted macrocycles) can run to completion
# instead of being killed at the default 30 min — the construction is deterministic, not a hang.
_UI_ISOLATE_TIMEOUT = int(os.environ.get("DELFIN_UI_ISOLATE_TIMEOUT", "1800"))


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
    timeout = timeout or _UI_ISOLATE_TIMEOUT
    if not timeout or timeout <= 0:
        timeout = None          # DELFIN_UI_ISOLATE_TIMEOUT<=0 -> run to completion, no kill
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
        proc = subprocess.run(
            [sys.executable, "-c", script],
            input=payload,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
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
    if _UI_ISOLATE_DEFAULT:
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
