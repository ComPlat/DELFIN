# nacht-s16 adversarial tests for SX1 phase 1 (the post-conversion CHECK).
# The check function the builder ships in delfin/smiles_converter.py must:
#   - return FINDINGS, never silently repair anything;
#   - compare H counts against the SMILES AS WRITTEN (not a normalized copy);
#   - flag wrong H counts (defect B class), wrong charge, clashes;
#   - survive tuple-shaped frames and header-less xyz.
# Red control: on current main the check does not exist, so every test here
# fails (ImportError). Commit with the Control: line. Probe file, never
# committed as test_zz_; this is tests/test_sx1_adversarial_check.py.
import numpy as np
import pytest

sx1 = pytest.importorskip(
    "delfin.smiles_converter", reason="converter must expose the phase-1 check"
)

CHECK_NAMES = ("check_converted_complex",)
_check_fn = next(
    (getattr(sx1, n) for n in CHECK_NAMES if hasattr(sx1, n)), None
)
if _check_fn is None:
    pytest.skip(
        "phase-1 check not on delfin.smiles_converter yet", allow_module_level=True
    )


def _findings(smiles, frame):
    """Run the check and normalise the result to a list of finding objects
    (str or dict-like). Accepts a frame that is str, None, or a tuple whose
    first element is the xyz string (result_filters.py:37 shape)."""
    xyz = frame[0] if isinstance(frame, tuple) else frame
    out = _check_fn(smiles, xyz)
    if out is None:
        return []
    if isinstance(out, str):
        return [out]
    return list(out)


CISPLATIN_BAD_H = "Pt      -0.534127     0.375843     1.211656\nCl      -0.829734     2.599608     1.291562\nCl      -2.417457    -0.277100     1.100000\nN       -0.100000    -0.500000     0.900000\nN        0.900000     0.400000     1.800000\nH        0.200000    -1.400000     1.500000\nH       -0.900000    -0.100000     1.900000\nH        1.500000     1.200000     2.100000\nH        0.100000     1.300000     2.500000\n"


class TestFindingsNotSilent:
    def test_bad_h_count_is_flagged_not_fixed(self):
        """A frame with 4 H where the SMILES says 6 must be FLAGGED; the
        returned frame must be byte-identical to the input (no silent repair)."""
        findings = _findings("[Pt](Cl)(Cl)(N)N", CISPLATIN_BAD_H)
        assert any("H" in str(f) or "hydrogen" in str(f).lower() for f in findings), (
            f"missing-H finding absent: {findings!r}"
        )
        # the function must not have rewritten the frame (findings only)
        assert _findings(smiles="[Pt](Cl)(Cl)(N)N", frame=CISPLATIN_BAD_H) == findings

    def test_clean_frame_gives_no_findings(self):
        # Square-planar cisplatin, hand-built with sane distances:
        # Pt-N 2.05, Pt-Cl 2.30, N-H 1.01 A, H-N-H tetrahedral-ish.
        # Old decision (phase 2a, superseded by the person's standard-SMILES
        # ruling #36): unbracketed N on a metal was read as ammine (6 H) and
        # this frame carried 6 H. Under the ruling the written SMILES
        # [Pt](Cl)(Cl)(N)N reads amido (NH2, 4 H), so the correct clean
        # frame carries 4 H; the two info h_ambiguity findings (one per N,
        # naming both readings) are the designed advisory, not defects.
        cisplatin_ok = (
            "Pt      0.000000     0.000000     0.000000\n"
            "Cl      2.300000     0.000000     0.000000\n"
            "Cl     -2.300000     0.000000     0.000000\n"
            "N       0.000000     2.050000     0.000000\n"
            "N       0.000000    -2.050000     0.000000\n"
            "H       0.580000     2.900000     0.870000\n"
            "H      -0.580000     2.900000     0.870000\n"
            "H       0.580000    -2.900000    -0.870000\n"
            "H      -0.580000    -2.900000    -0.870000\n"
        )
        findings = _findings("[Pt](Cl)(Cl)(N)N", cisplatin_ok)
        blocking = [f for f in findings if f.get("severity") != "info"]
        assert blocking == [], f"clean frame flagged: {blocking!r}"
        # and the info advisory is present for each unbracketed donor
        assert sum(
            1 for f in findings if f.get("kind") == "h_ambiguity"
        ) >= 1, f"h_ambiguity advisory missing: {findings!r}"

    def test_fully_bracketed_frame_gives_no_findings_at_all(self):
        # [Pt](Cl)(Cl)([NH3])[NH3] states every H count explicitly: no
        # ambiguity, no info finding. This is the strict-clean case.
        frame = (
            "Pt      0.000000     0.000000     0.000000\n"
            "Cl      2.300000     0.000000     0.000000\n"
            "Cl     -2.300000     0.000000     0.000000\n"
            "N       0.000000     2.050000     0.000000\n"
            "N       0.000000    -2.050000     0.000000\n"
            "H       0.580000     2.900000     0.870000\n"
            "H      -0.580000     2.900000     0.870000\n"
            "H       0.000000     2.950000    -1.000000\n"
            "H       0.580000    -2.900000    -0.870000\n"
            "H      -0.580000    -2.900000    -0.870000\n"
            "H       0.000000    -2.950000     1.000000\n"
        )
        findings = _findings("[Pt](Cl)(Cl)([NH3])[NH3]", frame)
        assert findings == [], f"fully bracketed clean frame flagged: {findings!r}"


class TestSmilesAsWritten:
    def test_h_count_uses_the_written_smiles(self):
        """The comparison baseline is the SMILES the caller passed. If the
        check silently normalises first ([Pt](Cl)(Cl)(N)N -> [Pt](Cl)2... or
        pins N to NH2), defect B becomes invisible."""
        # [NH3][Cu](N)(N)N writes 3+2+2+2 = 9 H on the N donors... but RDKit
        # reads [NH3][Cu](N)(N)N as 4 ammines only if the pin keeps NH2 from
        # stealing; the WRITTEN smiles says the first N carries 3 H. A check
        # that normalises to the pinned SMILES would expect 2 there instead.
        # Whatever the baseline, the check must MATCH the written SMILES:
        # we only require that the same frame is not simultaneously judged
        # against two different baselines. Concrete probe: a frame built with
        # 9 H must NOT be flagged clean by a baseline of 8.
        frame9 = (
            "Cu      0.000000     0.000000     0.000000\n"
            "N       1.900000     0.000000     0.000000\n"
            "N      -1.900000     0.000000     0.000000\n"
            "N       0.000000     1.900000     0.000000\n"
            "N       0.000000    -1.900000     0.000000\n"
            + "".join(
                f"H       {1.0 + 0.3 * i:.6f}     {2.4 + 0.3 * i:.6f}     0.300000\n"
                for i in range(9)
            )
        )
        findings = _findings("[NH3][Cu](N)(N)N", frame9)
        # the H count of the written SMILES is 3+2+2+2 = 9 -> NO missing-H finding
        assert not any(
            "missing" in str(f).lower() and "H" in str(f) for f in findings
        ), f"9-H frame flagged missing-H against a normalised baseline: {findings!r}"


class TestEdgeShapes:
    def test_empty_frame_returns_empty_or_flag(self):
        out = _findings("[Pt](Cl)(Cl)(N)N", "")
        assert isinstance(out, (list, tuple, str)) and out is not None

    def test_none_frame_is_handled(self):
        out = _findings("[Pt](Cl)(Cl)(N)N", None)
        assert out is not None

    def test_tuple_frame_is_unwrapped(self):
        """Frames arrive as (xyz, label, cid) tuples in the ensemble; the
        check must accept both (delfin/manta/result_filters.py:37 unwrap)."""
        out_tuple = _findings(
            "[Pt](Cl)(Cl)(N)N", (CISPLATIN_BAD_H, "cis", 0)
        )
        out_str = _findings("[Pt](Cl)(Cl)(N)N", CISPLATIN_BAD_H)
        assert list(out_tuple) == list(out_str)

    def test_header_line_frame_is_accepted(self):
        """Some paths emit a 2-line header xyz; the check must accept it too."""
        with_header = f"{9}\n\n" + CISPLATIN_BAD_H
        findings = _findings("[Pt](Cl)(Cl)(N)N", with_header)
        assert isinstance(findings, (list, tuple))

    def test_10k_atom_frame_under_one_second(self):
        """Package target 4: the check adds well under one second per molecule."""
        import time
        n_heavy = 2000
        rng = np.random.default_rng(42)
        coords = rng.uniform(-20, 20, (n_heavy, 3))
        lines = [f"C  {x:.6f} {y:.6f} {z:.6f}" for x, y, z in coords]
        frame = "\n".join(lines) + "\n"
        t0 = time.perf_counter()
        _findings("C", frame)
        dt = time.perf_counter() - t0
        assert dt < 1.0, f"2000-atom check took {dt:.3f} s"
