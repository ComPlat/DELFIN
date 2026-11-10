"""ORCA-probe acceptance for the OPI ESD inputs.

These tests hand OPI-built ESD inputs to ORCA's *own* parser through
``delfin.orca_knowledge.probe.probe_input`` so ORCA — the ground truth for
what its grammar accepts — decides whether the input is valid.  They need
OPI (py3.12 job) and the ORCA binary; when ORCA is not installed the probe
reports ``verdict == 'unavailable'`` and the test skips (it does not fail).

The probe guard refuses ``%moinp``/``MOREAD``/``*xyzfile`` and quoted
strings — every ESD state file carries at least one of those (the S1 job
reads S0.gbw, the check job reads S0.xyz, the rate jobs read Hessians) —
so neither variant can be the writer's complete multi-job file.  Both
variants are therefore *job sections the writer itself assembles*:

* ``state`` — the first ``$new_job`` section of the writer's deltaSCF-mode
  S0 file (``write_esd_state_opi_input`` + ``split_jobs``), verbatim: an
  RKS job with inline coordinates and the per-atom ``NewGTO`` path.
* ``rate`` — the writer's IC-rate job shape, assembled through the same
  helpers the writer uses (``_job_in_opi`` + ``_tddft_block_opi`` +
  ``_esd_block_lines``) minus the Hessian file names the guard forbids
  (``_esd_block_lines``'s header is empty).  A live ORCA runs this to the
  ORCA_ESD module and aborts there for the missing Hessian *data* — a
  runtime abort after the input was read, which ``probe_input`` reports as
  ``parsed`` (delfin/orca_knowledge/probe.py:277-280).  The grammar, i.e.
  what this task asks the probe to prove, is accepted.

Both variants: 6 atoms of ideal geometry (never a CSD/crystal structure),
1 nproc, 1000 MB maxcore, GFN2-XTB — a live probe finishes in seconds.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools.esd_opi import (  # noqa: E402
    _esd_block_lines,
    _import_opi,
    _job_in_opi,
    _rate_tokens,
    _tddft_block_opi,
    split_jobs,
    write_esd_state_opi_input,
)

# 6 atoms, ideal geometry (never a CSD/crystal structure).
_ATOMS = [
    ("Fe", (0.0, 0.0, 0.0)),
    ("N", (2.0, 0.0, 0.0)),
    ("N", (-2.0, 0.0, 0.0)),
    ("H", (2.5, 1.0, 0.0)),
    ("H", (-2.5, -1.0, 0.0)),
    ("H", (0.0, 2.5, 1.0)),
]

# Deliberately tiny level of theory (GFN2-XTB keeps the live probe to
# seconds).  No disp_corr/ri_jkx/relativity/aux_jk: resolve_level_of_theory
# would emit tokens ORCA's simple input does not define (verified live:
# "!none" is rejected with UNRECOGNIZED OR DUPLICATED KEYWORD).
_LOT = {
    "functional": "XTB2",
    "geom_opt": "OPT",
    "PAL": 1,
    "maxcore": 1000,
    "ESD_modus": "deltaSCF",
    "ESD_frequency": "no",
    "initial_guess": "",
    "temperature": "298.15",
}


def _write_xyz(path, comment: str = "synthetic probe geometry") -> None:
    lines = ["6", comment] + [
        f"{el} {x} {y} {z}" for el, (x, y, z) in _ATOMS
    ]
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def esd_probe_env(tmp_path, monkeypatch):
    """An ESD dir with the state xyz files the writer reads; cwd = tmp."""
    monkeypatch.chdir(tmp_path)
    esd_dir = tmp_path / "esd"
    esd_dir.mkdir()
    for st in ("S0", "S1", "T1"):
        _write_xyz(esd_dir / f"{st}.xyz")
    _write_xyz(tmp_path / "initial.xyz")
    return esd_dir


def _probed_state_text(esd_probe_env) -> str:
    """First job section of the writer's deltaSCF S0 file, verbatim."""
    path = write_esd_state_opi_input(
        esd_dir=esd_probe_env,
        out_dir=esd_probe_env.parent / "opi_state",
        state="S0",
        charge=0,
        solvent=None,
        metals=[],
        main_basisset="def2-SVP",
        metal_basisset=None,
        config=_LOT,
    )
    text = open(path, encoding="utf-8").read()
    return split_jobs(text)[0]


def _probed_rate_text(esd_probe_env) -> str:
    """The writer's IC-rate job shape, without Hessian file references."""
    opi = _import_opi()
    out_dir = esd_probe_env.parent / "opi_rate"
    out_dir.mkdir(parents=True, exist_ok=True)
    # DOHT=FALSE keeps the probe to the vertical-gradient path: the
    # Herztberg-Teller path requires the Hessian files the guard forbids
    # (ORCA manual esd.txt:105: "set DOHT to TRUE" *with* both Hessians).
    cfg = dict(_LOT, ESD_modus="TDDFT", NROOTS=3, TDDFT_TDA="TRUE",
               TDDFT_SOC="FALSE", DOHT="FALSE")
    toks = _rate_tokens(["XTB2"])
    toks.append("ESD(IC)")
    tddft = _tddft_block_opi(opi, cfg, nroots=3, iroot=1, dosoc=None)
    esd = _esd_block_lines(opi, cfg)  # header-less: no Hessian file names
    atoms = [opi["Atom"](el, coordinates=xyz) for el, xyz in _ATOMS]
    return _job_in_opi(
        opi, out_dir=out_dir, base="probe_rate", tokens=toks,
        maxcore=1000, pal=1, moinp=None,
        blocks=[tddft, esd], atoms=atoms, charge=0, multiplicity=1,
    ).read_text(encoding="utf-8")


@pytest.mark.parametrize("which", ["state", "rate"])
def test_opi_esd_inputs_accepted_by_orca_parser(which, esd_probe_env):
    text = (_probed_state_text(esd_probe_env) if which == "state"
            else _probed_rate_text(esd_probe_env))
    result = probe_input(text, timeout_s=300)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert result["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built ESD '{which}' input; messages: "
        f"{result['messages']}"
    )
