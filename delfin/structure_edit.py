"""Structure-department library (package XE1, SA1 of the wave-27 plan).

Pure helper functions behind the oxyz-editor expert's ``tool_geometry_*``
MCP tools. Every function takes xyz *text* in and returns xyz text / a JSON
object out, matching the existing geometry-tool contract
(``delfin.ops_server.server`` -> ``delfin.api`` -> ``delfin.tools.geometry_edit``).

The library reuses DELFIN's own chemistry primitives -- nothing is
re-implemented:

- xyz <-> Geom and the covalent-radius / connectivity tables live in
  ``delfin.tools.geometry_edit`` (``parse_xyz``, ``to_xyz``, ``Geom``,
  ``_covalent_radii``, ``_bond_graph``).
- renumbering a product onto an educt uses the universal mapper
  ``delfin.atom_mapping.map_atoms`` (a second mapper would be a defect).

Design rules the package target fixes:
- Every mutation returns a NEW structure; no tool writes a file. The
  structure expert writes the new file with the gate-checked ``write_file``
  (which refuses overwriting); an MCP tool writing a file would bypass the
  write-scope gate, so nothing here touches disk.
- A formula mismatch in a renumber is surfaced (ValueError), never hidden.
- ``fix_hydrogens`` only ever ADDS hydrogens, never removes an atom.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from delfin.tools.geometry_edit import (
    Geom,
    GeometryError,
    _bond_graph,
    _covalent_radii,
    parse_xyz,
    to_xyz,
)

try:  # atom_mapping imports numpy at module scope
    from delfin.atom_mapping import map_atoms
except Exception:  # pragma: no cover - import-time guard, exercised in CI
    map_atoms = None  # type: ignore[assignment]


# Distance (A) above which two positions count as "moved" in a diff.
_MOVED_ANGSTROM = 0.1

# Standard maximal valence for common main-group elements. A non-metal atom
# whose heavy-atom degree is below this needs hydrogen(s) placed.
_VALENCE: dict[str, int] = {
    "H": 1,
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "P": 3,
    "S": 2,
    "Cl": 1,
    "Br": 1,
    "I": 1,
}


def _parse(text: str, what: str = "structure") -> Geom:
    try:
        return parse_xyz(text)
    except GeometryError as exc:
        raise ValueError(f"{what} is not a valid xyz block: {exc}") from exc


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


# ---------------------------------------------------------------------------
# renumber (product onto educt)
# ---------------------------------------------------------------------------

def renumber_to_reference(
    xyz_text: str,
    reference_xyz_text: str,
    max_edits: int = 2,
) -> dict[str, Any]:
    """Reorder ``xyz_text`` onto ``reference_xyz_text`` atom order.

    Uses the universal educt<->product mapper ``delfin.atom_mapping.map_atoms``
    (never a second mapper). A formula mismatch raises ``ValueError`` and is
    never hidden. Returns a dict with the reordered structure, the verified
    flag, the renumbering table (product index -> reference index) and checks.
    """
    tgt = _parse(xyz_text, "target")
    ref = _parse(reference_xyz_text, "reference")
    if map_atoms is None:  # pragma: no cover - atom_mapping unavailable
        return {
            "reordered_xyz": "",
            "n_atoms": len(tgt),
            "renumbering": [],
            "verified": False,
            "rmsd": None,
            "error": "atom_mapping is unavailable on this host",
        }
    result = map_atoms(
        list(ref.symbols), list(ref.coords),
        list(tgt.symbols), list(tgt.coords),
        max_edits=max_edits,
    )
    order: list[int] = result["order"]  # target index per reference atom
    reorder: list[int] = [0] * len(tgt)
    for ref_i, tgt_i in enumerate(order):
        reorder[tgt_i] = ref_i
    table = [
        {"product_index": j, "reference_index": ref_i}
        for j, ref_i in enumerate(reorder)
    ]
    reordered = Geom(
        symbols=[tgt.symbols[order[i]] for i in range(len(order))],
        coords=[tgt.coords[order[i]] for i in range(len(order))],
        charge=tgt.charge,
        multiplicity=tgt.multiplicity,
        comment=f"renumbered onto reference ({result['method']})",
    )
    return {
        "reordered_xyz": to_xyz(reordered),
        "n_atoms": len(order),
        "renumbering": table,
        "n_bond_edits": result["n_bond_edits"],
        "rmsd": result["rmsd"],
        "method": result["method"],
        "verified": bool(result["verified"]),
        "checks": ["composition matching via atom_mapping"],
        "error": None,
    }


# ---------------------------------------------------------------------------
# fix hydrogens (add missing coordinating / implicit H)
# ---------------------------------------------------------------------------

def _heavy_neighbour_directions(
    g: Geom, atom: int, nbrs: dict[int, list[int]]
) -> list[tuple[float, float, float]]:
    """Unit vectors from ``atom`` toward each heavy neighbour (away faces)."""
    out: list[tuple[float, float, float]] = []
    for j in nbrs.get(atom, ()):
        if g.symbols[j] == "H":
            continue
        v = tuple(
            (g.coords[j][k] - g.coords[atom][k]) for k in range(3)
        )
        n = sum(x * x for x in v) ** 0.5
        if n > 1e-9:
            out.append(tuple(x / n for x in v))
    return out


def _place_one_h(
    g: Geom, atom: int, away: list[tuple[float, float, float]],
    bond_len: float,
) -> tuple[float, float, float]:
    """A candidate H position on ``atom``, pushed away from heavy neighbours."""
    # Sum of away vectors (may be empty -> fall back to +x hemisphere).
    s = [sum(v[k] for v in away) for k in range(3)]
    n = sum(x * x for x in s) ** 0.5
    if n > 1e-9:
        d = tuple(x / n for x in s)
    else:
        d = (1.0, 0.0, 0.0)
    return tuple(g.coords[atom][k] + d[k] * bond_len for k in range(3))


def fix_hydrogens(xyz_text: str) -> dict[str, Any]:
    """Add missing hydrogens to a structure; never removes an atom.

    A non-metal atom whose heavy-atom degree is below its standard valence
    gets the missing hydrogen(s) placed at DELFIN's covalent-radius bond
    length, pushed away from its heavy neighbours. Metals are never given
    guessed hydrogens (their coordination H is a donor convention, not a
    valence count).
    """
    g = _parse(xyz_text, "structure")
    radii = _covalent_radii(g)
    nbrs = _bond_graph(g)
    added: list[dict[str, Any]] = []
    new_symbols = list(g.symbols)
    new_coords = list(g.coords)
    for i, sym in enumerate(g.symbols):
        if sym == "H":
            continue
        need = _VALENCE.get(sym)
        if not need:
            continue  # metal or unknown -> no guessed H
        heavy_deg = sum(
            1 for j in nbrs.get(i, ()) if g.symbols[j] != "H"
        )
        missing = need - heavy_deg
        if missing <= 0:
            continue
        away = _heavy_neighbour_directions(g, i, nbrs)
        r_h_i = radii[i] + _COSMETIC_H_RADIUS
        for _ in range(missing):
            pos = _place_one_h(g, i, away, r_h_i)
            # push away from atoms already placed this pass too
            for existing_h in new_coords[len(g.symbols):]:
                if _dist(existing_h, pos) < 0.2:
                    d = tuple((pos[k] - existing_h[k]) for k in range(3))
                    dn = sum(x * x for x in d) ** 0.5
                    if dn > 1e-9:
                        pos = tuple(
                            pos[k] + (d[k] / dn) * 0.15 for k in range(3)
                        )
            new_symbols.append("H")
            new_coords.append(pos)
            added.append(
                {"index": len(new_symbols) - 1, "bonded_to": i, "element": "H"}
            )
    out = Geom(
        symbols=new_symbols,
        coords=new_coords,
        charge=g.charge,
        multiplicity=g.multiplicity,
        comment=g.comment or "hydrogens fixed",
    )
    return {
        "xyz": to_xyz(out),
        "n_added": len(added),
        "added_rows": added,
        "checks": [
            f"added {len(added)} hydrogen(s), removed none",
            "H placed at DELFIN covalent-radius bond length",
        ],
        "error": None,
    }


_COSMETIC_H_RADIUS = 0.63  # DELFIN's fallback covalent radius for H


# ---------------------------------------------------------------------------
# structure_diff (per-atom diff + checks)
# ---------------------------------------------------------------------------

def structure_diff(xyz_a_text: str, xyz_b_text: str) -> dict[str, Any]:
    """A per-atom diff of ``xyz_a`` against ``xyz_b`` plus checks.

    Reports changed element and moved atoms (position change above
    ``_MOVED_ANGSTROM``), the composition equality and the charge /
    multiplicity consistency. Atom positions compared positionally (same
    index); a renumbering is a separate operation (``renumber_to_reference``).
    """
    a = _parse(xyz_a_text, "structure A")
    b = _parse(xyz_b_text, "structure B")
    n = min(len(a), len(b))
    rows: list[dict[str, Any]] = []
    for i in range(n):
        sa, sb = a.symbols[i], b.symbols[i]
        da = _dist(a.coords[i], b.coords[i])
        rows.append({
            "index": i,
            "from_symbol": sa,
            "to_symbol": sb,
            "element_changed": sa != sb,
            "moved": da > _MOVED_ANGSTROM,
            "distance": round(da, 3),
        })
    for i in range(n, len(a)):
        rows.append({
            "index": i, "from_symbol": a.symbols[i], "to_symbol": None,
            "element_changed": True, "moved": True, "distance": None,
        })
    for i in range(n, len(b)):
        rows.append({
            "index": i, "from_symbol": None, "to_symbol": b.symbols[i],
            "element_changed": True, "moved": True, "distance": None,
        })
    return {
        "rows": rows,
        "composition_equal": Counter(a.symbols) == Counter(b.symbols),
        "charge_mult_match": (
            a.charge == b.charge and a.multiplicity == b.multiplicity
        ),
        "checks": [
            "composition compared by element count",
            "charge/multiplicity compared (0-based default)",
        ],
        "error": None,
    }
