"""Red tests (this review session) for two confirmed real-run defects found in
the ``orca_th_opi`` thermochemistry block (see scratch/REVIEW.md):

1. Finding 2 -- the entropy fallback returns the sign ORCA printed on the
   ``Total entropy correction`` line, which real ORCA prints NEGATIVE (the
   -T·S term applied in G = H - T·S), while the primary ``Final entropy term``
   line prints POSITIVE (+T·S).  The reader's ``ts_Eh``/``ts_kcal`` must be
   sign-stable across both labels.

2. Finding 1+3 -- a Hessian-reusing re-evaluation input must refuse when it
   cannot reproduce the source molecule: a provided-but-missing ``hessian``,
   or a real ``hessian`` with no caller-supplied ``geometry`` (a ``.hess``
   file carries no atomic coordinates, so a nominal single-He ``* xyz`` would
   compute the wrong molecule).

Run with the OPI interpreter (``DELFIN_OPI_REQUIRED=1``); skipped without OPI.
Fixtures are synthetic text -- never real ORCA output.
"""

import pytest

pytest.importorskip("opi")

from delfin.orca_read import read_thermochemistry  # noqa: E402
from delfin.tools import _registry  # noqa: E402

_HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"


def _thermo_body() -> str:
    """A one-temperature thermo block whose GIBBS section carries only the
    ``Total entropy correction`` line (the case a truncated-but-terminated
    output reaches) -- with the NEGATIVE sign real ORCA prints."""
    return _HEADER + "\n" + """------------------
THERMOCHEMISTRY AT 298.15K
------------------
Temperature                         ...     298.15 K
Pressure                            ...      1.00 atm
INNER ENERGY
Electronic energy                   ...   -151.55083691 Eh
Zero point energy                   ...      0.02631509 Eh      16.51 kcal/mol
ENTHALPY
Total Enthalpy                      ...   -151.52034633 Eh
ENTROPY
Thermal entropy correction          ...      0.00072596 Eh
Vibrational entropy computed according to the QRRHO of S. Grimme
Total entropy correction            ...     -0.02583637 Eh     -16.21 kcal/mol
GIBBS FREE ENERGY
Final Gibbs free energy             ...   -151.54618270 Eh
"""


def test_entropy_term_sign_is_stable_across_orca_line_labels(tmp_path):
    """Real ORCA prints ``Final entropy term`` as +T·S and ``Total entropy
    correction`` as -T·S (manual thermochemistry.txt 347 vs 356).  The reader
    must report the entropy TERM (T·S) positive regardless of which label
    survived, so ``ts_Eh``/``ts_kcal`` never flip sign between two otherwise
    identical outputs.  This fixture uses only the correction line (negative,
    real sign) to exercise the fallback."""
    out = tmp_path / "freq.out"
    out.write_text(_thermo_body())
    d = read_thermochemistry(out)[0]
    assert d.ts_Eh == pytest.approx(0.02583637), (
        "fallback must negate the -T·S correction so ts_Eh stays +T·S; "
        f"got {d.ts_Eh}"
    )
    assert d.ts_kcal == pytest.approx(16.21)


def test_reeval_refuses_hessian_reuse_without_source_geometry(tmp_path):
    """A Hessian-reusing re-evaluation input needs the source molecule's
    geometry: a ``.hess`` file carries no atomic coordinates (verified on the
    real h2o_optfreq.hess), so a nominal single-He ``* xyz`` would make ORCA
    compute thermochemistry of helium, not the source molecule.  Refuse rather
    than write a wrong input."""
    out = tmp_path / "freq.out"
    out.write_text(_thermo_body())
    hess = tmp_path / "freq.hess"  # exists on disk
    hess.write_text("$hessian\n3\n$act_atom 0\n")
    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        hessian=str(hess),
    )
    assert not result.ok
    assert result.error_kind.name in ("MISSING_INPUT", "PARSE")
    assert "geometry" in result.error.lower()


def test_reeval_refuses_missing_hessian_file(tmp_path):
    """A provided-but-missing ``hessian`` is a missing prerequisite and must be
    reported, not silently dropped (the block's own comment claims this; the
    code drops it)."""
    out = tmp_path / "freq.out"
    out.write_text(_thermo_body())
    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        hessian=str(tmp_path / "nope.hess"),
    )
    assert not result.ok
    assert result.error_kind.name in ("MISSING_INPUT", "PARSE")
    assert "hess" in result.error.lower()



def _write_source_inp(path, method="B3LYP", basis="def2-SVP",
                      charge=0, mult=1, atoms=None):
    """A SYNTHETIC source ORCA input whose geometry block carries the level the
    re-evaluation input must inherit: method/basis from the ``!`` line and
    charge/multiplicity from the ``* xyz <charge> <mult>`` header."""
    atoms = atoms or [
        ("O", 0.000000, 0.000000, 0.117300),
        ("H", 0.000000, 0.757200, -0.469200),
        ("H", 0.000000, -0.757200, -0.469200),
    ]
    lines = [_HEADER, f"! {method} {basis}", f"* xyz {charge} {mult}"]
    for el, x, y, z in atoms:
        lines.append(f"{el} {x!r} {y!r} {z!r}")
    lines.append("*")
    path.write_text("\n".join(lines) + "\n")
    return path


def test_reeval_input_carries_source_geometry_charge_mult_method(tmp_path):
    """A Hessian-reusing re-evaluation input must inherit the source run's
    geometry, charge, multiplicity AND method/basis (operator-A real-run FAIL:
    a nominal BP86/def2-SVP charge-0 singlet He input returns E_el -2.898 Eh
    for what should be H2O B3LYP/def2-SVP -76.3 Eh).  Everything is read from
    the source .inp.  The inhessname path must be quoted for ORCA to accept
    it (unquoted aborts with 'Invalid assignment in FREQ block')."""
    out = tmp_path / "freq.out"
    out.write_text(_thermo_body())
    hess = tmp_path / "freq.hess"
    hess.write_text("$hessian\n3\n$act_atom 0\n")
    inp = tmp_path / "freq.inp"
    _write_source_inp(inp)
    target = tmp_path / "reval.inp"

    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        re_eval_pressure=2.0,
        hessian=str(hess),
        source_inp=str(inp),
        input_path=target,
    )
    assert result.ok, result.error
    text = target.read_text()
    # the hess path is quoted so ORCA does not abort
    assert f'inhessname "{hess}"' in text, text
    # method/basis of the source run, not the nominal hardcoded level
    assert "!b3lyp" in text and "!def2-svp" in text, text
    # charge + multiplicity from the source run
    assert "* xyz 0 1" in text, text
    # geometry of the source run, not a single He atom
    assert "O" in text and "H" in text and "He" not in text, text


def test_reeval_refuses_without_readable_source_inp(tmp_path):
    """The source level (geometry/charge/mult/method) cannot be invented; if the
    source .inp is missing the block must refuse rather than write the nominal
    placeholder that produced the wrong real-run chemistry."""
    out = tmp_path / "freq.out"
    out.write_text(_thermo_body())
    hess = tmp_path / "freq.hess"
    hess.write_text("$hessian\n3\n$act_atom 0\n")
    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        hessian=str(hess),
        source_inp=str(tmp_path / "missing.inp"),
    )
    assert not result.ok
    assert result.error_kind.name in ("MISSING_INPUT", "PARSE")
    assert "inp" in result.error.lower()
