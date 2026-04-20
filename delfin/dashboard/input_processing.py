"""Input processing helpers: SMILES wrappers, resource parsing, sanitisation."""

import re

from delfin.smiles_converter import (
    smiles_to_xyz as _delfin_smiles_to_xyz,
    smiles_to_xyz_isomers as _delfin_smiles_to_xyz_isomers,
    smiles_to_xyz_quick as _delfin_smiles_to_xyz_quick,
    is_smiles_string as _delfin_is_smiles_string,
    contains_metal,
)
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
    quality_mode="fast",
):
    """Generate distinct coordination isomers for a SMILES string.

    Returns ``([(xyz_string, num_atoms, label), ...], error)``.

    ``quality_mode`` defaults to ``"fast"`` (12 ETKDG seeds, 1 chelate
    rank, 1 template) for the interactive dashboard so the
    ``CONVERT SMILES + UFF`` button responds within seconds on typical
    σ-complexes.  Pass ``"normal"`` (20 seeds) or ``"max"`` (40 seeds,
    3 ranks, 3 templates) for scripted / batch workflows that can
    afford the extra wall-clock time in exchange for a deeper
    candidate pool.
    """
    results, error = _delfin_smiles_to_xyz_isomers(
        smiles,
        apply_uff=apply_uff,
        collapse_label_variants=collapse_label_variants,
        include_binding_mode_isomers=include_binding_mode_isomers,
        hapto_approx=hapto_approx,
        deterministic=deterministic,
        quality_mode=quality_mode,
    )
    if error and hapto_approx is None and _is_hapto_failfast(error):
        results, error = _delfin_smiles_to_xyz_isomers(
            smiles,
            apply_uff=apply_uff,
            collapse_label_variants=collapse_label_variants,
            include_binding_mode_isomers=include_binding_mode_isomers,
            hapto_approx=True,
            deterministic=deterministic,
            quality_mode=quality_mode,
        )
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


def parse_inp_resources(inp_text):
    """Parse PAL (nprocs) and maxcore from ORCA ``.inp`` text.

    Returns ``(pal, maxcore)`` as ints or *None* if not found.
    """
    pal = None
    maxcore = None
    if not inp_text:
        return pal, maxcore
    pal_match = re.search(r'(?im)^\s*nprocs\s+(\d+)', inp_text)
    if pal_match:
        pal = int(pal_match.group(1))
    maxcore_match = re.search(r'(?im)^\s*%maxcore\s+(\d+)', inp_text)
    if maxcore_match:
        maxcore = int(maxcore_match.group(1))
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
