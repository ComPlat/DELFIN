"""F5: the orbital-analysis step reads real UNO occupations from real
ORCA output (fixtures from the CAS1 practice runs), not from typed
test data.

The parser is ``parse_uno_occupations`` in
``delfin/casscf_active_space.py``. It extracts the occupation numbers
from the ``UHF NATURAL ORBITALS`` / occupation-numbers block that
``!UNO`` + ``Print[P_UNO_OccNum]`` makes ORCA 6 write (manual 3.1.2,
5.1.11), and returns them in orbital-index order.

Fixture expectations (real values from the runs):
- ozone: pi pair 1.80143 / 0.19857, everything else near 2/0
- N2 @ 2.5 A: six fractional UNOs (indices 4-9), <S^2>=2.97
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from delfin.casscf_active_space import (
    check_active_occupations_safe,
    parse_casscf_occupations,
    parse_uno_occupations,
)

FIX = Path(__file__).parent / "fixtures" / "orca611"

# Real ORCA outputs are licence-restricted (ORCA EULA 3c-3e) and not
# committed; when the person's environment points this variable at
# the private copies, the real-run checks run instead of skipping.
_PRIVATE = os.environ.get("DELFIN_PRIVATE_FIXTURES")


def _fixture_dir(name: str) -> Path:
    d = FIX / name
    # Provenance is required for the cas1_* fixtures this package
    # adds; pre-existing fixtures are used read-only as-is.
    if name.startswith("cas1_"):
        assert (d / "PROVENANCE.md").is_file(), "fixture provenance missing"
    return d


class TestParseUnoOccupationsOnSyntheticOutput:
    def test_ozone_pi_pair_extracted(self):
        d = _fixture_dir("cas1_ozone_uno")
        occ = parse_uno_occupations(str(d))
        # Synthetic ozone-like pattern: the pi pair carries ~0.2
        # diradical character.
        assert occ[11] == pytest.approx(1.80, abs=1e-5)
        assert occ[12] == pytest.approx(0.20, abs=1e-5)
        # Core orbitals sit at exactly 2.0.
        assert occ[0] == pytest.approx(2.0, abs=1e-5)

    def test_n2_six_fractional_unos(self):
        d = _fixture_dir("cas1_n2_2p5_uno")
        occ = parse_uno_occupations(str(d))
        # Synthetic stretched-N2 pattern (indices 4-9).
        expected = [1.194, 1.047, 1.046, 0.953, 0.953, 0.806]
        got = [occ[i] for i in range(4, 10)]
        assert got == pytest.approx(expected, abs=1e-4)
        # The rest is pinned at 2.0 / 0.0.
        assert occ[0] == pytest.approx(2.0, abs=1e-5)
        assert occ[12] == pytest.approx(0.0, abs=1e-5)

    def test_n2_fractional_window_matches_manual_space(self):
        # The UNOs must be able to JUSTIFY the manual's CAS(6,6):
        # exactly six orbitals between 0.02 and 1.98.
        d = _fixture_dir("cas1_n2_2p5_uno")
        occ = parse_uno_occupations(str(d))
        partial = [o for o in occ if 0.02 < o < 1.98]
        assert len(partial) == 6

    def test_ozone_o2s_nearly_pinned_is_seen(self):
        # CAS(18,12) includes O 2s; in the UNOs those sit near 2.0
        # (synthetic: exactly 2.0 on the core block) -- the analysis
        # must parse them for the resize loop to act.
        d = _fixture_dir("cas1_ozone_uno")
        occ = parse_uno_occupations(str(d))
        assert occ[4] == pytest.approx(2.0, abs=1e-5)

    def test_no_output_returns_none(self):
        # A folder without an .out must return None (not crash), so
        # the workflow can say 'reference not computed yet'.
        assert parse_uno_occupations(str(FIX / "no_such_folder")) is None

    def test_converged_natural_occupations_last_block_wins(self, tmp_path):
        # A parser for the CASSCF "N(occ)=" natural-occupation lines
        # returns the LAST block (the converged iteration) in order.
        src = FIX / "cas1_n2_2p5_cas"
        assert (src / "PROVENANCE.md").is_file()
        synth = (src / "cas1_n2_2p5_cas_synth.out").read_text()
        d = tmp_path / "nocc"
        d.mkdir()
        (d / "cas.out").write_text(synth)
        occ = parse_casscf_occupations(str(d))
        # The LAST N(occ) block (iteration 9), not the first.
        assert occ == [1.27340, 1.07060, 1.07010, 0.92990, 0.92950,
                       0.72660]

    def test_parse_casscf_occupations_none_when_block_missing(self, tmp_path):
        d = tmp_path / "no_nocc"
        d.mkdir()
        (d / "cas.out").write_text("no occupation block here\n")
        assert parse_casscf_occupations(str(d)) is None

    def test_synth_casscf_ozone_occupations_flag_resize(self):
        # Synthetic CAS(18,12) ozone-like pattern: five near-2.0
        # orbitals (O 2s/sigma) must be flagged -> shrink proposal.
        d = _fixture_dir("cas1_ozone_cas")
        occ = parse_casscf_occupations(str(d))
        assert occ[0] == pytest.approx(1.99840, abs=1e-5)
        assert occ[8] == pytest.approx(1.77750, abs=1e-5)
        assert occ[9] == pytest.approx(0.26230, abs=1e-5)
        assert occ[11] == pytest.approx(0.05500, abs=1e-5)
        rep = check_active_occupations_safe(occ)
        assert rep.pinned_indices == [0, 1, 2, 3, 4]
        assert "shrink" in rep.recommendation

    def test_synth_casscf_n2_occupations_healthy(self):
        # Synthetic CAS(6,6) stretched-N2 pattern: healthy.
        d = _fixture_dir("cas1_n2_2p5_cas")
        occ = parse_casscf_occupations(str(d))
        assert occ[0] == pytest.approx(1.27340, abs=1e-5)
        assert occ[5] == pytest.approx(0.72660, abs=1e-5)
        rep = check_active_occupations_safe(occ)
        assert rep.pinned_indices == []
        assert "healthy" in rep.recommendation

    # Optional real-run check: only when the person's environment
    # provides the private ORCA outputs (ORCA EULA restricts
    # redistribution of ORCA-generated data; committed fixtures are
    # synthetic).
    @pytest.mark.skipif(
        _PRIVATE is None,
        reason="DELFIN_PRIVATE_FIXTURES unset: real ORCA outputs are "
               "licence-restricted and not committed")
    def test_real_casscf_ozone_occupations_flag_resize(self):
        d = Path(_PRIVATE) / "cas1_ozone_cas"
        occ = parse_casscf_occupations(str(d))
        assert occ[0] == pytest.approx(1.99836, abs=1e-5)
        rep = check_active_occupations_safe(occ)
        assert rep.pinned_indices == [0, 1, 2, 3, 4]

    @pytest.mark.skipif(
        _PRIVATE is None,
        reason="DELFIN_PRIVATE_FIXTURES unset: real ORCA outputs are "
               "licence-restricted and not committed")
    def test_real_casscf_n2_occupations_healthy(self):
        d = Path(_PRIVATE) / "cas1_n2_2p5_cas"
        occ = parse_casscf_occupations(str(d))
        assert occ[0] == pytest.approx(1.27338, abs=1e-5)
        rep = check_active_occupations_safe(occ)
        assert rep.pinned_indices == []
