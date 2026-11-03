"""The GOAT ensemble reader on the comment-line format ORCA writes.

``delfin.tools.adapters.orca_goat_opi.read_goat_ensemble`` reads
``<base>.finalensemble.xyz``.  Format measured on finished GOAT runs of
ORCA 6.0.1 and 6.1.1 (inspected, not copied): each frame's comment line
is the absolute energy in Eh followed by ``converged=true``, e.g.
``-70.5391993062 converged=true``.  The manual (Table 4.9, READENSEMBLE)
says only "Energy (float)".  The reader took the LAST token as the energy
and so raised ValueError on every real file.

Input: a SYNTHETIC ensemble (invented numbers in that format; never ORCA
output, LICENSE_GUARD.md family 2).  Output: absolute energies in Eh,
relative energies in kcal/mol (lowest = 0), structures, the global
minimum file.  1 Eh = 627.509474 kcal/mol (CODATA 2018).
"""

from __future__ import annotations

import pytest

from delfin.tools.adapters.orca_goat_opi import read_goat_ensemble

# SYNTHETIC: invented coordinates/energies for a water "ensemble" of three.
SYNTHETIC_ENSEMBLE = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
3
-76.4000000000 converged=true
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000
3
-76.3990000000 converged=true
O 0.000000 0.000000 0.100000
H 0.757000 0.586000 0.100000
H -0.757000 0.586000 0.100000
3
-76.3980000000
O 0.000000 0.000000 0.200000
H 0.757000 0.586000 0.200000
H -0.757000 0.586000 0.200000
"""

SYNTHETIC_GLOBAL_MINIMUM = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
3
-76.4000000000
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000
"""

EH_TO_KCAL = 627.509474


@pytest.fixture()
def run_dir(tmp_path):
    (tmp_path / "XTB_GOAT.finalensemble.xyz").write_text(SYNTHETIC_ENSEMBLE)
    (tmp_path / "XTB_GOAT.globalminimum.xyz").write_text(SYNTHETIC_GLOBAL_MINIMUM)
    return tmp_path


def test_energy_with_converged_suffix_is_read(run_dir):
    out = read_goat_ensemble(run_dir, "XTB_GOAT")
    assert out["n_conformers"] == 3
    assert out["energies_eh"] == pytest.approx((-76.4, -76.399, -76.398))


def test_relative_energies_in_kcal_mol(run_dir):
    out = read_goat_ensemble(run_dir, "XTB_GOAT")
    assert out["energies_kcal_mol"] == pytest.approx(
        (0.0, 0.001 * EH_TO_KCAL, 0.002 * EH_TO_KCAL), abs=1e-6)


def test_structures_and_global_minimum(run_dir):
    out = read_goat_ensemble(run_dir, "XTB_GOAT")
    assert len(out["structures"]) == 3
    assert out["structures"][2][0] == ("O", 0.0, 0.0, 0.2)
    assert out["global_minimum_file"] == run_dir / "XTB_GOAT.globalminimum.xyz"


def test_comment_without_a_number_is_refused(tmp_path):
    bad = SYNTHETIC_ENSEMBLE.replace("-76.3990000000 converged=true", "converged=true")
    (tmp_path / "b.finalensemble.xyz").write_text(bad)
    with pytest.raises(ValueError, match="energy"):
        read_goat_ensemble(tmp_path, "b")


def test_keyword_form_keeps_its_relative_kcal_reading(tmp_path):
    """'Energy <x>' (no measured ORCA writes it) stays relative kcal/mol."""
    text = ("# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "1\nEnergy 0.0\nH 0.0 0.0 0.0\n1\nEnergy 1.5\nH 0.0 0.0 0.1\n")
    (tmp_path / "k.finalensemble.xyz").write_text(text)
    out = read_goat_ensemble(tmp_path, "k")
    assert out["energies_kcal_mol"] == (0.0, 1.5)
    assert out["energies_eh"] == ()


def test_mixed_forms_are_refused(tmp_path):
    bad = SYNTHETIC_ENSEMBLE.replace("-76.3980000000\n", "Energy 1.2\n")
    (tmp_path / "m.finalensemble.xyz").write_text(bad)
    with pytest.raises(ValueError, match="mixes"):
        read_goat_ensemble(tmp_path, "m")


def test_missing_ensemble_reports_nothing(tmp_path):
    out = read_goat_ensemble(tmp_path, "none")
    assert out["n_conformers"] == 0 and out["ensemble_file"] is None
