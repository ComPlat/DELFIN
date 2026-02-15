"""Input processing helpers: SMILES wrappers, resource parsing, sanitisation."""

import re

from delfin.smiles_converter import (
    smiles_to_xyz as _delfin_smiles_to_xyz,
    smiles_to_xyz_isomers as _delfin_smiles_to_xyz_isomers,
    is_smiles_string as _delfin_is_smiles_string,
    contains_metal,
)


def smiles_to_xyz(smiles):
    """Convert a SMILES string to XYZ coordinates.

    Returns ``(xyz_string, num_atoms, method, error)``.
    """
    xyz_string, error = _delfin_smiles_to_xyz(smiles)
    if error:
        return None, 0, None, error
    num_atoms = sum(1 for line in xyz_string.splitlines() if line.strip())
    method = 'delfin.smiles_converter'
    return xyz_string, num_atoms, method, None


def smiles_to_xyz_isomers(smiles):
    """Generate distinct coordination isomers for a SMILES string.

    Returns ``([(xyz_string, num_atoms, label), ...], error)``.
    """
    results, error = _delfin_smiles_to_xyz_isomers(smiles)
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
