"""Adversarial tests for the RC2 transition-metal fixtures (reviewer s4).

The builder's tests/test_rc2_tm_fixtures.py pins the happy path. This file
attacks the contract around it: electron-count parity against multiplicity,
the method line in each committed .inp, absence of a stacked basis
override on the -3c composite, output blocks a real run must contain, and
provenance/identity hygiene.

Red control: on the branch as committed before the fixtures land, every
case here fails with "fixture not yet committed" -- it turns green only
when tests/fixtures/orca611/rc2_*/ carry files that satisfy the contract.
"""
import re
from pathlib import Path

import pytest

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir


def _fixtures(name: str):
    """The specific orca611 subfolder: legacy path while its files are
    committed, else the private store (skip without the env var)."""
    return lc1_fixture_dir(f"tests/fixtures/orca611/{name}")

#: name -> (charge, multiplicity, n_atoms, element counts that decide the
#: electron count) -- the physics contract of the package task.
SYSTEMS = {
    "rc2_fe_co5": {"charge": 0, "mult": 1, "n_atoms": 11,
                   "counts": {"Fe": 1, "C": 5, "O": 5}},
    "rc2_cu_nh3_4": {"charge": 2, "mult": 2, "n_atoms": 17,
                     "counts": {"Cu": 1, "N": 4, "H": 12}},
    "rc2_ru_o4": {"charge": 0, "mult": 1, "n_atoms": 5,
                  "counts": {"Ru": 1, "O": 4}},
}

#: valence electron counts for the parity check (neutral atoms)
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "Fe": 26, "Cu": 29, "Ru": 44}


def _inp_text(name: str) -> str:
    p = _fixtures(name) / f"{name}.inp"
    if not p.is_file():
        pytest_fail = f"fixture not yet committed: {p}"
        raise AssertionError(pytest_fail)
    return p.read_text(encoding="utf-8", errors="replace")


def _xyz_block(text: str):
    m = re.search(r"\*\s*xyz\s+(-?\d+)\s+(\d+)\s*\n(.*?)\n\s*\*", text,
                  re.DOTALL)
    assert m, "no '* xyz <charge> <mult>' block in the .inp"
    return int(m.group(1)), int(m.group(2)), m.group(3)


def test_rc2_systems_parities_hold():
    """Electron count vs multiplicity: parity must agree in the .inp."""
    for name, want in SYSTEMS.items():
        text = _inp_text(name)
        charge, mult, xyz = _xyz_block(text)
        assert charge == want["charge"], name
        assert mult == want["mult"], name
        atoms = [ln.split()[0] for ln in xyz.strip().splitlines() if ln.strip()]
        assert len(atoms) == want["n_atoms"], (name, len(atoms))
        from collections import Counter
        got = Counter(atoms)
        assert dict(got) == want["counts"], (name, dict(got))
        n_el = sum(Z[a] for a in atoms) - charge
        if want["mult"] % 2 == 1:  # singlet: even electron count
            assert n_el % 2 == 0, (name, n_el, "odd electrons, singlet")
        else:  # doublet: odd electron count
            assert n_el % 2 == 1, (name, n_el, "even electrons, doublet")


def test_rc2_method_line_is_the_composite_and_no_basis_override():
    """r2SCAN-3c keyword line present; no %basis block, no ZORA stacked.

    The -3c composite fixes its own basis (def2-mTZVPP + def2-ECP; ECPs
    from Rb on, so Ru/Cu/Fe stay all-electron) and its own dispersion.
    A per-atom %basis or a ZORA keyword would silently replace the
    composite's fitted basis.
    """
    for name in SYSTEMS:
        text = _inp_text(name)
        assert "r2SCAN-3c" in text, (name, "no composite keyword line")
        assert not re.search(r"%basis\b", text, re.IGNORECASE), name
        assert "ZORA" not in text, name
        assert re.search(r"%pal\s+nprocs\s+\d+\s+end", text, re.IGNORECASE), \
            (name, "missing %pal block")
        assert re.search(r"%maxcore\s+\d+", text, re.IGNORECASE), name


def test_rc2_out_carries_the_blocks_a_real_run_needs():
    """The committed .out must be a real ORCA output, not a stub.

    A fixture .out that is missing the geometry-optimisation termination,
    the frequency block or the final energy cannot pin the physics. Also
    guards the identity hygiene: the banner host must be the placeholder.
    """
    for name in SYSTEMS:
        p = _fixtures(name) / f"{name}.out"
        if not p.is_file():
            raise AssertionError(f"fixture not yet committed: {p}")
        text = p.read_text(encoding="utf-8", errors="replace")
        assert "Program Version" in text, (name, "not a real ORCA output")
        assert "FINAL SINGLE POINT ENERGY" in text, name
        assert "THE OPTIMIZATION HAS CONVERGED" in text, \
            (name, "geometry optimisation did not report convergence")
        assert "*** OPTIMIZATION RUN DONE ***" in text, \
            (name, "no completed optimisation marker")
        assert "THERMOCHEMISTRY AT" in text, (name, "no thermochemistry block")
        assert "ORCA TERMINATED NORMALLY" in text, (name, "abnormal end")
        for line in text.splitlines():
            if "Host name:" in line:
                assert line.split("Host name:", 1)[1].strip() == "fixture-host", \
                    (name, "real host name leaked into the fixture")
        # Account names are covered repo-wide by
        # tests/test_no_account_name_is_committed.py; this tuple must not
        # name one either (a public repo must not carry it in a check).
        for pattern in ("/home/", "/pfs/"):
            assert pattern not in text, (name, pattern)


def test_rc2_provenance_documents_the_run():
    """PROVENANCE.md per system: job id, ORCA version, core-hours, CPU
    efficiency line, redaction statement; and the relativistic note for
    RuO4 must name the reason (r2SCAN-3c/def2-ECP covers scalar relativity
    for 4d; no ZORA stacked)."""
    for name in SYSTEMS:
        p = _fixtures(name) / "PROVENANCE.md"
        if not p.is_file():
            raise AssertionError(f"fixture not yet committed: {p}")
        text = p.read_text(encoding="utf-8", errors="replace")
        out_text = (_fixtures(name) / f"{name}.out").read_text(
            encoding="utf-8", errors="replace")
        # Field casing varies between the files ("Core-hours:"/"CPU
        # Efficiency" vs "core-hours"/"CPU efficiency"); the substance is
        # the presence of both numbers, so the check is case-insensitive.
        # Each PROVENANCE documents ITS OWN run id (the Cu chain documents
        # v1/v2/v3, Fe/Ru their single job): pin the system's own id.
        own_job = {"rc2_fe_co5": "8172015", "rc2_cu_nh3_4": "8173353",
                   "rc2_ru_o4": "8172021"}[name]
        lower = text.lower()
        for field in (own_job, "core-hours", "cpu efficiency"):
            assert field in lower, (name, field)
        # The run's version must be stated from the .out banner, not the
        # suite directory name: all three rc2 runs print
        # "Program Version 6.1.0 - RELEASE" at :54 (verified in the
        # committed .out files); h2o_optfreq is a true 6.1.1 build, so
        # "orca611" names the suite, not the build.
        assert "ORCA 6.1.0" in text, (name, "ORCA version not stated")
        assert "Program Version 6.1.0" in out_text, \
            (name, "version not pinned to the .out banner")
        if name == "rc2_ru_o4":
            assert re.search(r"ZORA|relativistic|ECP", text, re.IGNORECASE), \
                (name, "relativistic treatment not documented")
