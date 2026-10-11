"""Summary of one molecule file for the chat: formula, atom count, frames.

``molecule()`` reads an XYZ text (one or more frames) or a Gaussian cube --
inline, or from a path inside ``workspace_root`` -- and returns the facts the
chat states about it: formula, atom count, frame count, format. The result
carries no markup. The dashboard draws the structure with the chat's existing
3D card (``tab_agent._mol3d_card_html``, bundled 3Dmol, no network); the
model and the terminal receive only the ``text`` line built here.

File content is untrusted. Only tokens that are element symbols (or atomic
numbers) are counted; the XYZ comment line is never read as an atom, so text
placed in a file cannot reach the formula the model is shown.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

# A molecule file larger than this (bytes) is refused. Equal to the chat's
# 3D card limit (tab_agent._MOL3D_MAX_BYTES), so a file accepted here is one
# the card can also draw.
_MAX_FILE_BYTES = 4_000_000

_VALID_KINDS = ("xyz", "multixyz", "cube")

# Element symbols by atomic number, 1..118. Cube atom lines give the atomic
# number; XYZ files written by some programs do too. A truncated table
# reported every transition metal of a complex as "X".
_ELEMENTS = ("", *"""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu
Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs
Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl
Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh
Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split())
_SYMBOLS = {s.lower(): s for s in _ELEMENTS if s}
# "Fe", "FE", "fe", and labelled forms such as "Fe1" or "C12".
_ATOM_TOKEN = re.compile(r"([A-Za-z]{1,2})\d*")


def _atomic_symbol(z: int) -> str:
    return _ELEMENTS[z] if 0 < z < len(_ELEMENTS) else "X"


def _element(token: str) -> str | None:
    """The element a coordinate line's first token names, or None.

    Accepts a symbol in any case, a symbol with a numeric label, or an
    atomic number. Anything else -- a word, a sentence fragment -- is not
    an atom and is not counted.
    """
    if token.isdigit():
        z = int(token)
        return _ELEMENTS[z] if 0 < z < len(_ELEMENTS) else None
    m = _ATOM_TOKEN.fullmatch(token)
    if m is None:
        return None
    return _SYMBOLS.get(m.group(1).lower())


def _formula(counts: dict[str, int]) -> str:
    """Formula from an {element: count} mapping: C first, H second when C is
    present (Hill order), then alphabetically; count '1' omitted."""
    def key(kv):
        el = kv[0]
        if "C" in counts:
            return (0 if el == "C" else 1 if el == "H" else 2, el)
        return (2, el)
    order = sorted(counts.items(), key=key)
    return " ".join(f"{el}{n}" if n > 1 else el for el, n in order)


def _parse_xyz(text: str) -> tuple[list[dict[str, int]], int]:
    """Return (per-frame element counts, number of frames) for XYZ text.

    A frame is a count line, ONE comment line (which may be empty), and
    ``count`` coordinate lines. Blank lines are skipped only where a count
    line is expected, never inside a frame: an empty comment line is a
    comment line, and dropping it shifted every frame by one atom.
    """
    frames: list[dict[str, int]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if not head:
            i += 1
            continue
        if len(head) != 1 or not head[0].isdigit():
            break                    # not a count line: no further frames
        n = int(head[0])
        body = lines[i + 2:i + 2 + n]
        if n <= 0 or len(body) < n:
            break                    # truncated last frame: not counted
        current: dict[str, int] = {}
        for ln in body:
            parts = ln.split()
            el = _element(parts[0]) if len(parts) >= 4 else None
            if el is not None:
                current[el] = current.get(el, 0) + 1
        if not current:
            break
        frames.append(current)
        i += 2 + n
    if not frames:
        raise ValueError("no atom lines found in XYZ content")
    return frames, len(frames)


def _parse_cube(text: str) -> dict:
    """Return {natoms, formula, counts} for Gaussian cube text.

    Cube lines: two comments, ``natoms ox oy oz``, three voxel-count lines,
    then ``|natoms|`` atom lines ``Z charge x y z``. A negative natoms marks
    an orbital cube (an MO-index line follows the atoms); the atoms are the
    same.
    """
    lines = text.splitlines()
    if len(lines) < 6:
        raise ValueError("cube content too short to be a valid cube")
    try:
        natoms = abs(int(lines[2].split()[0]))
    except (ValueError, IndexError):
        raise ValueError("cube file has no atom count") from None
    counts: dict[str, int] = {}
    for ln in lines[6:6 + natoms]:
        parts = ln.split()
        if len(parts) < 5:
            continue
        try:
            z = int(float(parts[0]))
        except ValueError:
            continue
        sym = _atomic_symbol(z)
        counts[sym] = counts.get(sym, 0) + 1
    if not counts:
        raise ValueError("cube atom block is empty")
    return {"natoms": sum(counts.values()), "counts": counts,
            "formula": _formula(counts)}


def _detect_kind(text: str, suffix: str) -> str:
    """Kind from the file suffix, else from the content."""
    sfx = (suffix or "").lower()
    if sfx == ".cube":
        return "cube"
    if sfx != ".xyz":
        # A cube's third line is `<natoms> <ox> <oy> <oz>`.
        lines = text.splitlines()
        if len(lines) >= 6:
            third = lines[2].split()
            if len(third) == 4 and third[0].lstrip("-").isdigit():
                return "cube"
    _, frames = _parse_xyz(text)
    return "multixyz" if frames > 1 else "xyz"


def _load_path(path: str, workspace_root: str | None) -> tuple[str, str]:
    """Read ``path``, refusing it when it resolves outside ``workspace_root``
    (symlinks and ``..`` resolved first) or exceeds ``_MAX_FILE_BYTES``.
    Returns ``(content, suffix)``.
    """
    p = Path(path).expanduser()
    if workspace_root is not None:
        root = Path(workspace_root).resolve()
        try:
            inside = p.resolve().is_relative_to(root)
        except (ValueError, OSError):
            inside = False
        if not inside:
            raise ValueError(f"path outside the workspace is refused: {path}")
    try:
        if not p.is_file():
            raise ValueError(f"no such file: {path}")
        size = p.stat().st_size
    except OSError as exc:
        raise ValueError(f"cannot read {path}: {exc}") from None
    if size > _MAX_FILE_BYTES:
        raise ValueError(
            f"file too large to show ({size} bytes > {_MAX_FILE_BYTES})")
    try:
        content = p.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"cannot read {path}: {exc}") from None
    return content, p.suffix


def molecule(
    path_or_xyz: str | os.PathLike,
    *,
    kind: str | None = None,
    workspace_root: str | None = None,
) -> dict:
    """Summarise one molecule.

    Args:
        path_or_xyz: a path (str or os.PathLike) to an XYZ/cube file, or the
            raw XYZ/cube text.
        kind: ``"xyz"``, ``"multixyz"`` or ``"cube"``; detected from the
            suffix / content when None.
        workspace_root: root that file paths must stay inside.

    Returns ``{text, kind, frames, atom_count, formula, path}``; ``text`` is
    the plain line the model and the terminal are shown.

    Raises:
        ValueError: unknown kind, unparseable content, a path outside
            ``workspace_root``, or a file over the size cap.
    """
    if kind is not None and kind not in _VALID_KINDS:
        raise ValueError(f"unknown kind {kind!r}; expected {_VALID_KINDS}")

    path_or_xyz = os.fspath(path_or_xyz)

    suffix = ""
    path: str | None = None
    if "\n" in path_or_xyz:
        content = path_or_xyz
    else:
        content, suffix = _load_path(path_or_xyz, workspace_root)
        path = path_or_xyz

    eff_kind = kind or _detect_kind(content, suffix)

    frame_count = 1
    if eff_kind in ("xyz", "multixyz"):
        frames, frame_count = _parse_xyz(content)
        atom_count = sum(frames[0].values())
        formula = _formula(frames[0])
        if frame_count > 1:
            eff_kind = "multixyz"
    else:
        info = _parse_cube(content)
        atom_count = info["natoms"]
        formula = info["formula"]

    if eff_kind == "multixyz":
        kind_label = "multi-frame XYZ"
        extra = f"\nFrames: {frame_count}"
    else:
        kind_label = eff_kind.upper()
        extra = ""

    path_line = f"\nFile: {path}" if path else ""
    text = (
        f"Molecule: {formula}\n"
        f"Atoms: {atom_count}\n"
        f"Format: {kind_label}{extra}{path_line}"
    )
    return {
        "text": text,
        "kind": eff_kind,
        "frames": frame_count,
        "atom_count": atom_count,
        "formula": formula,
        "path": path,
    }
