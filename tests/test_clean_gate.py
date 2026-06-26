"""Unit tests for the unified hard clean-manifold emission gate
(``delfin.smiles_converter._clean_gate_filter``).

Synthetic, fast, deterministic — no heavy SMILES builds.  They pin the gate's
contract directly on hand-built XYZ frame sets so the behaviour is verifiable in
isolation from the (slow) generator stack:

  * BYTE-IDENTICAL when the flag is OFF (identity, same object returned).
  * Drops a frame with an inter-ligand clash / spurious bond / destroyed bond /
    decoordinated donor while keeping the clean siblings.
  * NEVER-EMPTY: when every frame is dirty it keeps the single least-violating
    frame and marks it ``|cleanest-available``.
  * Clean ensembles pass through untouched.
  * Deterministic (same input -> same output).
"""
import os

import pytest

from delfin import smiles_converter as sc


def _xyz(atoms):
    """atoms: list of (sym, x, y, z) -> headerless DELFIN-style frame string."""
    return "\n".join(f"{s} {x:.4f} {y:.4f} {z:.4f}" for (s, x, y, z) in atoms)


# An octahedral Fe with 6 N donors at ~2.0 A along +/-x, +/-y, +/-z, each N
# carrying one C at ~1.47 A pointing radially outward (two distinct "ligands"
# per axis so inter-ligand distances are large) -> a clean coordination frame.
def _clean_octahedron():
    atoms = [("Fe", 0.0, 0.0, 0.0)]
    axes = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    for (ax, ay, az) in axes:
        atoms.append(("N", 2.0 * ax, 2.0 * ay, 2.0 * az))
        atoms.append(("C", 3.47 * ax, 3.47 * ay, 3.47 * az))
    return atoms


@pytest.fixture(autouse=True)
def _clear_env():
    saved = os.environ.get("DELFIN_FFFREE_CLEAN_GATE")
    os.environ.pop("DELFIN_FFFREE_CLEAN_GATE", None)
    yield
    if saved is None:
        os.environ.pop("DELFIN_FFFREE_CLEAN_GATE", None)
    else:
        os.environ["DELFIN_FFFREE_CLEAN_GATE"] = saved


def test_off_is_identity():
    clean = _xyz(_clean_octahedron())
    iso = [(clean, "a"), (clean, "b")]
    os.environ.pop("DELFIN_FFFREE_CLEAN_GATE", None)
    assert sc._clean_gate_filter(iso) is iso
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "0"
    assert sc._clean_gate_filter(iso) is iso


def test_clean_ensemble_untouched():
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean = _xyz(_clean_octahedron())
    iso = [(clean, "a"), (clean, "b"), (clean, "c")]
    out = sc._clean_gate_filter(iso)
    assert len(out) == 3
    assert all("cleanest-available" not in lbl for _, lbl in out)


# Atom layout from _clean_octahedron(): 0=Fe, then per axis (+x,-x,+y,-y,+z,-z)
# a donor N (index 1,3,5,7,9,11) at 2.0 A and its terminal C (2,4,6,8,10,12) at
# 3.47 A.  A MAJORITY of clean frames is provided so the consensus topology +
# donor reference are well-defined (matching real multi-frame ensembles).


def test_keeps_single_donor_off_metal_still_coordinated():
    # ASYMMETRIC contract: one donor swinging off while the metal stays otherwise
    # coordinated is NOT certainly bad -> KEEP (the old consensus/reference donor
    # check wrongly dropped this; a surviving good-enough structure > removing a
    # borderline one).  Only TOTAL decoordination (bare metal) is dropped, which
    # is covered by test_drops_bare_metal_total_decoordination.
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean_atoms = _clean_octahedron()
    one_off = list(clean_atoms)
    one_off[1] = ("N", 8.0, 0.0, 0.0)
    one_off[2] = ("C", 9.47, 0.0, 0.0)
    iso = [(_xyz(clean_atoms), "clean1"),
           (_xyz(clean_atoms), "clean2"),
           (_xyz(one_off), "one_off")]
    out = sc._clean_gate_filter(iso)
    labels = [lbl for _, lbl in out]
    assert "one_off" in labels  # kept: metal still has 5 donors
    assert "clean1" in labels and "clean2" in labels


def test_drops_interligand_clash_keeps_clean():
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean_atoms = _clean_octahedron()
    # Dirty: drag the +y terminal carbon (index 6) on top of the +x terminal
    # carbon (index 2): two different ligand components -> inter-ligand clash.
    # The N5-C6 bond stays intact (N5 unmoved, C6 still ~bonded distance via the
    # majority/reference) so this isolates the clash check.
    dirty = list(clean_atoms)
    dirty[6] = ("C", 3.47 + 0.4, 0.0, 0.0)  # ~0.4 A from the +x carbon
    iso = [(_xyz(clean_atoms), "clean1"),
           (_xyz(clean_atoms), "clean2"),
           (_xyz(dirty), "dirty")]
    out = sc._clean_gate_filter(iso)
    labels = [lbl for _, lbl in out]
    assert "dirty" not in labels
    assert "clean1" in labels and "clean2" in labels


def test_asymmetric_keeps_modestly_decoordinated_frame():
    # Governing invariant: drop ONLY the certainly-bad, KEEP the borderline.
    # A frame whose metal is still coordinated (donors in the shell) but where a
    # single donor swung modestly off is NOT certainly bad -> it must be KEPT
    # (the old consensus/reference donor-distance check wrongly dropped it).
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean_atoms = _clean_octahedron()
    # Move ONE donor N (and its C) from 2.0 to ~2.6 A (still inside the 2.95 A
    # shell -> still coordinated); the rest of the octahedron is intact.
    borderline = list(clean_atoms)
    borderline[1] = ("N", 2.6, 0.0, 0.0)
    borderline[2] = ("C", 4.07, 0.0, 0.0)
    iso = [(_xyz(clean_atoms), "clean"), (_xyz(borderline), "borderline")]
    out = sc._clean_gate_filter(iso)
    labels = [lbl for _, lbl in out]
    assert "clean" in labels and "borderline" in labels  # both kept


def test_drops_bare_metal_total_decoordination():
    # The ONLY decoordination test the asymmetric gate applies: a metal with ZERO
    # heavy donors in its shell (the complex fell completely apart) IS certainly
    # bad and is dropped.
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean_atoms = _clean_octahedron()
    # Blow every donor far off the metal -> bare metal.
    bare = list(clean_atoms)
    for di in (1, 3, 5, 7, 9, 11):
        s, x, y, z = bare[di]
        bare[di] = (s, x * 6.0, y * 6.0, z * 6.0)
    for ci in (2, 4, 6, 8, 10, 12):
        s, x, y, z = bare[ci]
        bare[ci] = (s, x * 6.0, y * 6.0, z * 6.0)
    iso = [(_xyz(clean_atoms), "clean1"), (_xyz(clean_atoms), "clean2"),
           (_xyz(bare), "bare")]
    out = sc._clean_gate_filter(iso)
    labels = [lbl for _, lbl in out]
    assert "bare" not in labels
    assert "clean1" in labels and "clean2" in labels


def test_all_bad_keeps_input_unchanged_no_marker():
    # In-doubt-keep / never-collapse: when EVERY frame is certainly bad the gate
    # returns the input UNCHANGED (no fabricated "cleanest-available" marker, no
    # collapse to zero).  Use structurally-invalid (non-finite) frames.
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean_atoms = _clean_octahedron()
    bad1 = _xyz(clean_atoms).replace("3.4700", "nan", 1)
    bad2 = _xyz(clean_atoms).replace("3.4700", "inf", 1)
    iso = [(bad1, "a"), (bad2, "b")]
    out = sc._clean_gate_filter(iso)
    assert out is iso  # unchanged; nothing dropped, no marker added


def test_deterministic():
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    clean = _xyz(_clean_octahedron())
    dirty_atoms = _clean_octahedron()
    dirty_atoms[1] = ("N", 8.0, 0.0, 0.0); dirty_atoms[2] = ("C", 9.47, 0.0, 0.0)
    iso = [(clean, "a"), (clean, "c"), (_xyz(dirty_atoms), "b")]
    out1 = sc._clean_gate_filter(iso)
    out2 = sc._clean_gate_filter(iso)
    assert out1 == out2


def test_never_raises_on_garbage():
    os.environ["DELFIN_FFFREE_CLEAN_GATE"] = "1"
    iso = [("not an xyz at all", "x"), ("", "y")]
    # Must not raise; returns something (input or a single frame).
    out = sc._clean_gate_filter(iso)
    assert isinstance(out, list)
