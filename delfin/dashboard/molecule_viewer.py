"""3D molecule visualisation helpers using py3Dmol."""

import py3Dmol
from IPython.display import clear_output

DEFAULT_3DMOL_STYLE = {
    'stick': {
        'colorscheme': 'Jmol',
        'radius': 0.11,
        'singleBonds': False,
        'doubleBondScaling': 0.65,
        'tripleBondScaling': 0.65,
    },
    'sphere': {'colorscheme': 'Jmol', 'scale': 0.28},
}
DEFAULT_3DMOL_STYLE_JS = (
    '{'
    'stick:{colorscheme:"Jmol",radius:0.11,singleBonds:false,doubleBondScaling:0.65,tripleBondScaling:0.65},'
    'sphere:{colorscheme:"Jmol",scale:0.28}'
    '}'
)
DEFAULT_3DMOL_BACKGROUND = 'white'
DEFAULT_3DMOL_ZOOM = 0.90


def apply_molecule_view_style(view, zoom=DEFAULT_3DMOL_ZOOM):
    """Apply a shared ChemDarwin/MSILES-like style to a py3Dmol viewer."""
    view.setStyle({}, DEFAULT_3DMOL_STYLE)
    view.setBackgroundColor(DEFAULT_3DMOL_BACKGROUND)
    view.zoomTo()
    view.center()
    if zoom is not None:
        view.zoom(zoom)
    view.render()
    return view


def render_xyz_in_output(output_widget, xyz_text, width=560, height=420):
    """Render an XYZ string inside an ipywidgets Output widget."""
    with output_widget:
        clear_output()
        if not xyz_text or not xyz_text.strip():
            print('No coordinates to display.')
            return
        view = py3Dmol.view(width=width, height=height)
        view.addModel(xyz_text, 'xyz')
        apply_molecule_view_style(view)
        view.show()


def coord_to_xyz(coord_text):
    """Convert TURBOMOLE ``coord`` format to XYZ.

    Handles both ``$coord`` blocks (Bohr) and plain Cartesian lines.
    Returns an XYZ string with atom count + comment, or *None*.
    """
    lines = coord_text.strip().split('\n')
    xyz_lines = []
    bohr_to_ang = 0.529177249
    in_coord = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('$coord'):
            in_coord = True
            continue
        if stripped.startswith('$'):
            in_coord = False
            continue
        if in_coord:
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    x = float(parts[0]) * bohr_to_ang
                    y = float(parts[1]) * bohr_to_ang
                    z = float(parts[2]) * bohr_to_ang
                    xyz_lines.append(
                        f'{parts[3].capitalize()}  {x:12.6f}  {y:12.6f}  {z:12.6f}'
                    )
                except (ValueError, IndexError):
                    continue
    # Fallback: plain Cartesian lines (element x y z)
    if not xyz_lines:
        for line in lines:
            parts = line.split()
            if len(parts) >= 4 and parts[0][0].isalpha():
                try:
                    float(parts[1])
                    xyz_lines.append(line)
                except (ValueError, IndexError):
                    continue
    if not xyz_lines:
        return None
    return f'{len(xyz_lines)}\nTM Preview\n' + '\n'.join(xyz_lines)


def parse_xyz_frames(content):
    """Parse a multi-frame XYZ file into a list of ``(comment, xyz_block, n_atoms)`` tuples."""
    frames = []
    lines = content.strip().split('\n')
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        if i + 1 >= len(lines):
            break
        comment = lines[i + 1].strip()
        coord_lines = []
        for j in range(n_atoms):
            if i + 2 + j < len(lines):
                coord_lines.append(lines[i + 2 + j])
        if len(coord_lines) == n_atoms:
            xyz_block = '\n'.join(coord_lines)
            frames.append((comment, xyz_block, n_atoms))
        i += 2 + n_atoms
    return frames


def strip_xyz_header(text):
    """Remove the XYZ header (atom count + comment line) if present."""
    text = text.strip()
    if not text:
        return text
    lines = text.split('\n')
    if len(lines) < 2:
        return text
    first_line = lines[0].strip()
    try:
        int(first_line)
        return '\n'.join(lines[2:]).strip()
    except ValueError:
        return text
