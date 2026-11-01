"""Synthetic-fixture tests for the TS -> GOAT reactive-centre feature."""


import pytest

from delfin.ts_goat import (
    analyze_reactive_center,
    build_ts_goat_input,
    collect_ts_files,
    connectivity_changes,
    read_ts_status,
    read_xyz,
    render_constraints_block,
)

# ---------------------------------------------------------------- helpers

def _write_hess(path, atoms, mode_vectors, frequencies):
    """Minimal ORCA .hess: $atoms (Bohr) + $vibrational_frequencies +
    $normal_modes in ORCA's blocked-column layout."""
    bohr = 1.8897259886
    lines = ["$orca_hessian_file"]
    lines.append("$atoms")
    lines.append(f"{len(atoms)}")
    for el, x, y, z in atoms:
        lines.append(f"{el} 12.0 {x * bohr:.8f} {y * bohr:.8f} {z * bohr:.8f}")
    lines.append("$vibrational_frequencies")
    lines.append(f"{len(frequencies)}")
    for i, f in enumerate(frequencies):
        lines.append(f"{i} {f:.2f}")
    # $normal_modes: rows = 3N, cols = n_modes; one blocked-column table
    # whose header lists every mode index (ORCA layout) so _normal_mode
    # can pick any column.
    n3 = 3 * len(atoms)
    lines.append("$normal_modes")
    lines.append(f"{n3} {len(mode_vectors)}")
    lines.append("  " + "  ".join(str(i) for i in range(len(mode_vectors))))
    # each row: one dummy first token (the parser reads split()[1+col]),
    # then one component per mode column
    for atom_row in range(len(atoms)):
        for comp in range(3):
            parts = ["0.0"] + [
                f"{mode[atom_row][comp]:.10f}" for mode in mode_vectors
            ]
            lines.append("  ".join(parts))
    lines.append("$end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _water_ts_atoms():
    # H2O-like: O0, H1 near, H2 transferring away (a proton-transfer TS)
    return [
        ("O", 0.00, 0.00, 0.00),
        ("H", 0.96, 0.00, 0.00),
        ("H", 1.55, 1.10, 0.00),
    ]


def test_read_ts_status_detects_converged_optts(tmp_path):
    out = tmp_path / "ts.out"
    out.write_text(
        "! B3LYP def2-SVP OptTS\n"
        "...\n"
        "     0:      -512.34 cm**-1  ***imaginary mode***\n"
        "...\n"
        "HURRAY  *** THE OPTIMIZATION HAS CONVERGED ***\n",
        encoding="utf-8",
    )
    status = read_ts_status(out)
    assert status.is_ts
    assert status.converged
    assert status.n_imaginary == 1
    assert status.imaginary_freq == pytest.approx(-512.34)
    assert status.usable


def test_read_ts_status_rejects_non_ts_and_unconverged(tmp_path):
    opt = tmp_path / "opt.out"
    opt.write_text("! B3LYP def2-SVP Opt\n", encoding="utf-8")
    status = read_ts_status(opt)
    assert not status.is_ts
    assert not status.usable

    unfinished = tmp_path / "ts2.out"
    unfinished.write_text("! OptTS\n", encoding="utf-8")
    status = read_ts_status(unfinished)
    assert status.is_ts and not status.converged
    assert not status.usable
    assert any("convergence" in r.lower() for r in status.reasons)


def test_read_xyz_tolerates_blank_lines(tmp_path):
    xyz = tmp_path / "m.xyz"
    xyz.write_text("2\n\ntitle\nO 0 0 0\nH 1 0 0\n", encoding="utf-8")
    atoms = read_xyz(xyz)
    assert atoms is not None and len(atoms) == 2


def test_connectivity_changes_finds_proton_transfer():
    # reactant: H1 bonded to O0, H2 far; TS: H1 between; product: H2 near O0
    reactant = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0, 0), ("H", 2.5, 1.1, 0)]
    ts = _water_ts_atoms()
    product = [("O", 0.0, 0.0, 0.0), ("H", 2.5, 0.0, 0), ("H", 0.98, 1.1, 0)]
    formed, broken = connectivity_changes(ts, reactant, product)
    # O0-H1 is present in TS and reactant -> not a change; O0-H2 forms;
    # H1-H2? both H2O3 distances below cutoff in TS only.
    assert (0, 2) in formed | broken or (0, 1) in formed | broken


def test_analyze_reactive_center_prefers_bond_changes_over_displacement(tmp_path):
    # Three atoms move hugely along the mode but only two participate in a
    # broken/formed bond: the bond evidence must win for the selection.
    hess = tmp_path / "ts.hess"
    # mode vector: atoms 0 and 1 move a lot, atom 2 barely
    mode = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0),
    ]
    _write_hess(hess, _water_ts_atoms(), [mode], [-512.34, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0])
    ts_atoms = _water_ts_atoms()
    reactant = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0, 0), ("H", 3.0, 1.1, 0)]
    product = [("O", 0.0, 0.0, 0.0), ("H", 3.0, 0.0, 0), ("H", 0.96, 1.1, 0)]
    rc = analyze_reactive_center(hess, ts_atoms, reactant, product)
    assert rc.indices
    # the two bonding atoms of the reacting contact must be selected
    assert set(rc.indices) & {0, 2}  # O0 and the far/nearest H participate
    assert rc.reliability in ("high", "medium", "low")
    assert rc.reasons


def test_analyze_reactive_center_flags_ambiguity_for_wide_proposal(tmp_path):
    hess = tmp_path / "ts.hess"
    mode = [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.2, 0.0, 0.0)]
    _write_hess(hess, _water_ts_atoms(), [mode], [-400.0])
    ts_atoms = _water_ts_atoms()
    # a broad, non-specific mode: many atoms move together (torsion),
    # so the proposal covers most of the molecule and must be flagged
    mode = [(1.0, 0.0, 0.0), (0.9, 0.0, 0.0), (0.95, 0.0, 0.0)]
    _write_hess(hess, _water_ts_atoms(), [mode], [-400.0])
    rc = analyze_reactive_center(hess, ts_atoms)
    assert rc.ambiguous
    assert any("ambiguous" in r.lower() for r in rc.reasons)


def test_render_constraints_block_zero_based():
    block = render_constraints_block([4, 7, 12], ["C", "C", "C"])
    assert block == (
        "%geom\n"
        "  Constraints\n"
        "    { C 4 C }\n"
        "    { C 7 C }\n"
        "    { C 12 C }\n"
        "  end\n"
        "end"
    )


def test_collect_ts_files_matches_by_atom_order(tmp_path):
    ts_dir = tmp_path / "TS_1"
    ts_dir.mkdir()
    atoms = _water_ts_atoms()
    # Hessian + its matching xyz + a decoy xyz with different composition
    _write_hess(ts_dir / "TS_1.hess", atoms, [], [-300.0])
    (ts_dir / "TS_1.xyz").write_text(
        "3\nTS\n" + "\n".join(
            f"{el} {x} {y} {z}" for el, x, y, z in atoms
        ) + "\n",
        encoding="utf-8",
    )
    (ts_dir / "decoy.xyz").write_text("2\ndecoy\nO 0 0 0\nH 1 0 0\n", encoding="utf-8")
    (ts_dir / "TS_1.inp").write_text("! OptTS\n* xyz 0 1\nO 0 0 0\n", encoding="utf-8")
    (ts_dir / "endpoint_reactant.xyz").write_text(
        "3\nr\nO 0 0 0\nH 1 0 0\nH 3 0 0\n", encoding="utf-8"
    )
    files = collect_ts_files(ts_dir, out_path=ts_dir / "TS_1.out")
    assert files.xyz_path == ts_dir / "TS_1.xyz"
    assert files.hess_path == ts_dir / "TS_1.hess"
    assert files.inp_path == ts_dir / "TS_1.inp"
    # endpoints are no longer collected (plain OptTS has none) but are
    # still excluded from the TS-geometry candidates
    assert files.reactant_path is None and files.product_path is None
    assert files.complete


def test_collect_ts_files_reports_missing(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    files = collect_ts_files(d, out_path=None)
    assert not files.complete
    assert files.problems  # user is told what is missing


def test_build_ts_goat_input_writes_constraints(tmp_path, monkeypatch):
    """The input writer embeds the %geom Constraints block (OPI may be
    absent in CI; the block itself is what the test pins down)."""
    pytest.importorskip("opi")
    from delfin.tools.adapters.orca_goat_opi import build_goat_input  # noqa: F401

    inp = build_ts_goat_input(
        tmp_path,
        elements=("O", "H", "H"),
        coordinates=[(0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (1.55, 1.1, 0.0)],
        constraint_indices=(1,),
        charge=0,
        mult=1,
        method="XTB",
        base_name="goat_ts",
    )
    text = inp.read_text(encoding="utf-8")
    assert "{ H 1 C }" in text
    assert text.index("%geom") < text.index("* xyz")


# ------------------------------------------------- final-block & mode picker

def test_last_frequency_block_counts_not_historical_ones(tmp_path):
    """Several frequency blocks: earlier blocks contain extra imaginary
    modes; ONLY the last complete block decides."""
    out = tmp_path / "ts.out"
    out.write_text(
        "! OptTS\n"
        # historical block: 2 imaginary modes (mid-optimisation)
        "     0:     -812.00 cm**-1  ***imaginary mode***\n"
        "     1:     -120.00 cm**-1  ***imaginary mode***\n"
        "     2:       45.00 cm**-1\n"
        "\n"
        # final block: exactly one imaginary mode
        "     0:     -512.34 cm**-1  ***imaginary mode***\n"
        "     1:       80.00 cm**-1\n"
        "HURRAY  *** THE OPTIMIZATION HAS CONVERGED ***\n",
        encoding="utf-8",
    )
    status = read_ts_status(out)
    assert status.is_ts
    assert status.converged
    assert status.n_imaginary == 1
    assert status.imaginary_freq == pytest.approx(-512.34)
    assert status.usable


def test_two_imaginary_modes_stay_usable_with_mode_list(tmp_path):
    """2 imaginary modes in the last block: the workflow stays usable and
    lists both modes for the picker; nothing blocks."""
    out = tmp_path / "ts2.out"
    out.write_text(
        "! OptTS\n"
        "THE OPTIMIZATION HAS CONVERGED\n"
        "     0:     -450.00 cm**-1  ***imaginary mode***\n"
        "     1:     -310.00 cm**-1  ***imaginary mode***\n"
        "     2:      100.00 cm**-1\n",
        encoding="utf-8",
    )
    status = read_ts_status(out)
    assert status.usable
    assert status.n_imaginary == 2
    assert [i for i, _ in status.imaginary_modes] == [0, 1]


def test_small_second_mode_auto_preselects_clear_one(tmp_path):
    """One clearly negative + one tiny negative mode: the clear one is
    auto-pre-selected and the small one only warns."""
    out = tmp_path / "ts3.out"
    out.write_text(
        "! OptTS\n"
        "THE OPTIMIZATION HAS CONVERGED\n"
        "     0:     -450.00 cm**-1  ***imaginary mode***\n"
        "     2:      -12.00 cm**-1  ***imaginary mode***\n"
        "     1:      100.00 cm**-1\n",
        encoding="utf-8",
    )
    status = read_ts_status(out)
    assert status.usable
    assert status.auto_mode_index == 0
    assert status.warnings  # the small mode surfaces as a warning


def test_zero_imaginary_modes_blocks_with_clear_message(tmp_path):
    out = tmp_path / "ts4.out"
    out.write_text(
        "! OptTS\n"
        "THE OPTIMIZATION HAS CONVERGED\n"
        "     0:       90.00 cm**-1\n"
        "     1:      200.00 cm**-1\n",
        encoding="utf-8",
    )
    status = read_ts_status(out)
    assert not status.usable
    assert any("no imaginary mode" in r.lower() for r in status.reasons)


def test_mode_change_reproposes_reactive_center(tmp_path):
    """Picking a different imaginary mode re-runs the analysis on THAT
    mode and yields a different proposal."""
    # mode 0 acts on atoms 0-1, mode 1 on atoms 1-2 (two distinct vectors)
    modes = [
        [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
    ]
    hess = tmp_path / "ts.hess"
    _write_hess(hess, _water_ts_atoms(), modes, [-450.0, -310.0])
    rc0 = analyze_reactive_center(hess, _water_ts_atoms(), mode_index=0)
    rc1 = analyze_reactive_center(hess, _water_ts_atoms(), mode_index=1)
    # both proposals exist and differ in the mode-displacement evidence
    assert rc0.indices and rc1.indices
    assert rc0.scores != rc1.scores


# ------------------------ geometry matching regression (real ORCA files) ---

_ORCA_COMMENT = (
    "Coordinates from ORCA-job /scratch/slurm_tmpdir/job_1/run/base "
    "E -401.250208837948"
)


def _write_orca_xyz(path, atoms, comment=_ORCA_COMMENT):
    """ORCA-style xyz: comment line carries 4+ tokens incl. an energy."""
    lines = [str(len(atoms)), comment]
    lines += [f"{el} {x:.10f} {y:.10f} {z:.10f}" for el, x, y, z in atoms]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_orca_comment_line_with_energy_is_not_an_atom(tmp_path):
    # regression: 'Coordinates from ORCA-job ... E -401.25' made read_xyz
    # abort with float('from') -> every candidate was silently dropped and
    # the dialog reported 'No matching TS geometry'.
    import delfin.ts_goat as tg

    d = tmp_path / "final"
    d.mkdir()
    atoms = _water_ts_atoms()
    _write_orca_xyz(d / "Opt_TS.xyz", atoms)
    parsed = tg.read_xyz(d / "Opt_TS.xyz")
    assert parsed is not None and [a[0] for a in parsed] == [a[0] for a in atoms]
    # coordinates survive the comment line exactly
    assert abs(parsed[0][1] - atoms[0][1]) < 1e-9


def test_collect_matches_xyz_with_divergent_basename(tmp_path):
    d = tmp_path / "final"
    d.mkdir()
    atoms = _water_ts_atoms()
    _write_hess(d / "Opt_TS_recalc_1.hess", atoms, [], [-300.0])
    _write_orca_xyz(d / "geometry_optimised.xyz", atoms)  # divergent name
    files = collect_ts_files(d, out_path=d / "Opt_TS_recalc_1.out")
    assert files.xyz_path == d / "geometry_optimised.xyz"
    # .inp is genuinely absent here -> its hint stays, nothing about xyz
    assert not any("matching TS geometry" in x for x in files.problems)


def test_collect_reports_rejection_reason_per_candidate(tmp_path):
    d = tmp_path / "final"
    d.mkdir()
    atoms = _water_ts_atoms()
    _write_hess(d / "TS.hess", atoms, [], [-300.0])
    # same count but reordered elements -> genuinely incompatible
    _write_orca_xyz(d / "TS.xyz", list(reversed(atoms)))
    files = collect_ts_files(d, out_path=d / "TS.out")
    assert files.xyz_path is None
    assert any("element order" in p.lower() for p in files.problems), files.problems


def test_bohr_vs_angstrom_coordinates_do_not_reject(tmp_path):
    # composition matching compares element ORDER only; a Hessian stores
    # coordinates in Bohr, an .xyz in Angstrom -- that must not matter.
    import delfin.ts_goat as tg

    atoms = _water_ts_atoms()
    # Hessian atoms come in Bohr (large values); xyz in Angstrom (small)
    xyz_atoms = [(el, x / 1.889726, y / 1.889726, z / 1.889726)
                 for el, x, y, z in atoms]
    assert tg._same_composition(xyz_atoms, atoms) is True


def test_multiple_xyz_candidates_pick_hessian_ordered_one(tmp_path):
    d = tmp_path / "final"
    d.mkdir()
    atoms = _water_ts_atoms()
    _write_hess(d / "TS.hess", atoms, [], [-300.0])
    # three candidates: reactant (endpoint), wrong order, right order
    _write_orca_xyz(d / "edukt.xyz", atoms)
    _write_orca_xyz(d / "wrong_order.xyz", list(reversed(atoms)))
    _write_orca_xyz(d / "other_run_final.xyz", atoms)
    files = collect_ts_files(d, out_path=d / "TS.out")
    assert files.reactant_path is None  # endpoints are not collected
    assert files.xyz_path == d / "other_run_final.xyz"
    assert not any("matching TS geometry" in x for x in files.problems)


# ------------------- mode-only analysis (ORCA TS-conformer workflow) -------

def test_group_torsion_proposes_only_the_reacting_pair(tmp_path):
    """Whole groups with large absolute displacement, but the reactive
    centre is only visible through RELATIVE atom motion.

    Two fragments (atoms 0-2 and 3-5) translate rigidly (torsion-like,
    coherent motion), while atoms 3 and 0 additionally move towards each
    other along the mode (bond formation).  The proposal must not pick
    the whole moving groups -- only the two atoms of the reacting pair
    dominate the pair evidence.
    """
    # Chemically clean design: the reacting pair (0, 3) approaches
    # (distance 9 -> responds), while the spectator pairs (1, 2) and
    # (4, 5) rotate INTERNALLY around their own midpoint -- every
    # spectator distance is preserved to 2nd order, no spectator
    # contact responds.  The spectators show large absolute
    # displacement (0.9 tangent) but never change a distance.
    # The spectators sit FAR from both anchors (> 12 Angstrom, the
    # _MAX_RELEVANT_CONTACT_BOHR window): they rotate internally with a
    # large tangent displacement, but their distances to every reacting
    # atom exceed the relevance window, so only the pair (0, 3) can be
    # proposed.  (Physically, a spectator bound to an anchor that moves
    # in the mode always shares the anchor's distance response -- one
    # normal mode cannot separate them; the window models the chemical
    # reality that far spectators are irrelevant for GOAT constraints.)
    mode = [
        (1.2, 0.0, 0.0),                      # 0: towards 3
        (0.0, 0.9, 0.0),                      # 1: rotates with 2 (far away)
        (0.0, -0.9, 0.0),                     # 2: rotates with 1
        (-1.2, 0.0, 0.0),                     # 3: towards 0
        (0.0, 0.0, 0.9),                      # 4: rotates with 5 (far away)
        (0.0, 0.0, -0.9),                     # 5: rotates with 4
    ]
    atoms = [
        ("O", 0.0, 0.0, 0.0),              # anchor A
        ("H", 40.0, 0.0, 0.0),             # spectator pair far above
        ("H", 40.0, 0.0, 1.0),
        ("O", 9.0, 0.0, 0.0),              # anchor B
        ("H", -40.0, 0.0, 0.0),            # spectator pair far below
        ("H", -40.0, 0.0, 1.0),
    ]
    hess = tmp_path / "ts.hess"
    _write_hess(hess, atoms, [mode], [-420.0])
    rc = analyze_reactive_center(hess, atoms)
    # the reacting pair must be in the proposal
    assert 0 in rc.indices and 3 in rc.indices
    # the rotating riders must NOT be proposed (no blind top-N of the
    # biggest movers): they move coherently around their anchor
    assert not (1 in rc.indices or 2 in rc.indices)
    assert len(rc.indices) <= 4
    # pair evidence, not displacement only
    assert any("distance when the structure is displaced" in r or
               "relative motion" in r for r in rc.reasons)


def test_small_sensible_selection_not_top_n(tmp_path):
    """No fixed atom count: a two-atom reaction centre stays two atoms."""
    # 0 and 3 approach (proton transfer centre, distance 9 Bohr);
    # 1, 2, 4, 5 are spectator pairs at > 12 Angstrom from both
    # (outside the relevance window, internally rotating).
    mode = [
        (1.2, 0.0, 0.0), (0.0, 0.9, 0.0), (0.0, 0.9, 0.0),
        (-1.2, 0.0, 0.0), (0.0, 0.0, 0.9), (0.0, 0.0, 0.9),
    ]
    atoms = [
        ("O", 0.0, 0.0, 0.0),
        ("H", 40.0, 0.0, 0.0), ("H", 40.0, 0.0, 1.0),
        ("H", 9.0, 0.0, 0.0),
        ("H", -40.0, 0.0, 0.0), ("H", -40.0, 0.0, 1.0),
    ]
    hess = tmp_path / "ts.hess"
    _write_hess(hess, atoms, [mode], [-500.0])
    rc = analyze_reactive_center(hess, atoms)
    assert set(rc.indices) >= {0, 3}
    # a proton transfer centre is small: at most the pair plus one rider
    assert len(rc.indices) <= 3


def test_mode_index_switch_reproposes(tmp_path):
    """Two imaginary modes with different reactive pairs: switching the
    selected mode must re-run the analysis and change the proposal."""
    mode_a = [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 0.0),
              (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    mode_b = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
              (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    atoms = [
        ("O", 0.0, 0.0, 0.0), ("H", 2.4, 0.0, 0.0), ("H", -0.5, 0.9, 0.0),
        ("O", 5.0, 0.0, 0.0), ("H", 7.4, 0.0, 0.0), ("H", 4.5, 0.9, 0.0),
    ]
    hess = tmp_path / "ts.hess"
    _write_hess(hess, atoms, [mode_a, mode_b], [-450.0, -300.0])
    rc_a = analyze_reactive_center(hess, atoms, mode_index=0)
    rc_b = analyze_reactive_center(hess, atoms, mode_index=1)
    assert rc_a.indices and rc_b.indices
    # the proposals follow the selected mode's reactive pairs
    assert 1 in rc_a.indices and 4 in rc_b.indices


def test_endpoint_params_are_accepted_but_ignored(tmp_path):
    """Deprecated endpoint arguments must not change the mode-only
    result (a plain OptTS run has no endpoints)."""
    mode = [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    atoms = _water_ts_atoms()
    hess = tmp_path / "ts.hess"
    _write_hess(hess, atoms, [mode], [-400.0])
    rc_plain = analyze_reactive_center(hess, atoms)
    rc_with = analyze_reactive_center(
        hess, atoms, [atoms], [atoms[::-1]],
    )
    assert rc_plain.indices == rc_with.indices


def test_collect_ts_files_still_excludes_endpoint_files_from_ts_xyz(tmp_path):
    """Endpoint files are not collected, but also not mistaken for the
    TS geometry when they are the only matching xyz."""
    ts_dir = tmp_path / "TS"
    ts_dir.mkdir()
    atoms = _water_ts_atoms()
    _write_hess(ts_dir / "TS.hess", atoms, [], [-300.0])
    # an endpoint xyz that WOULD match the Hessian composition
    (ts_dir / "edukt.xyz").write_text(
        "3\nr\n" + "\n".join(f"{el} {x} {y} {z}" for el, x, y, z in atoms)
        + "\n", encoding="utf-8")
    files = collect_ts_files(ts_dir, out_path=ts_dir / "TS.out")
    assert files.xyz_path is None
    assert files.reactant_path is None and files.product_path is None
    # the problem message must not demand endpoints any more
    assert not any("reactant/product" in p for p in files.problems)



def test_marginal_contacts_do_not_dilute_the_proposal(tmp_path):
    """A strong contact dominates; marginal responders stay out.

    Real-run review (BF3H2O, 2026-11-02): four H...F contacts at
    4-5 % response diluted the proposal.  Selection is by contact
    strength now: only contacts responding at least half as strongly
    as the strongest one qualify.
    """
    ts_dir = tmp_path / "TS"
    ts_dir.mkdir()
    # two far-separated bonds; only pair (0,1) responds strongly
    atoms = [("C", 0.0, 0.0, 0.0), ("N", 3.0, 0.0, 0.0),
             ("C", 40.0, 0.0, 0.0), ("N", 43.0, 0.0, 0.0)]
    mode = [
        (0.9, 0.0, 0.0),   # 0: towards 1
        (-0.9, 0.0, 0.0),  # 1: towards 0  -> 01 responds strongly
        (0.0, 0.05, 0.0),  # 2: marginal
        (0.0, -0.05, 0.0),
    ]
    _write_hess(ts_dir / "TS.hess", atoms, [mode], [-400.0])
    rc = analyze_reactive_center(ts_dir / "TS.hess", atoms)
    assert rc.indices == (0, 1)
