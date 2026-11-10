"""ORCA-probe acceptance for the ``orca_th_opi`` re-evaluation input.

The re-evaluation input the ``orca_th_opi`` block writes is handed to ORCA's
own parser through ``delfin.orca_knowledge.probe.probe_input`` so ORCA — the
ground truth for what its grammar accepts — decides whether it is valid.  The
probe needs OPI (py3.12 job) and the ORCA binary; when ORCA is not installed
the probe reports ``verdict == 'unavailable'`` and the test skips (it does
not fail).  The probed input is deliberately small: a BP86/def2-SVP ``%freq``
re-evaluation at 500 K over a single helium atom at the origin, well within
the probe guard's limits (<= 2 nprocs, <= 2000 maxcore, <= 6 atoms).  The
``inhessname`` variant is NOT probed: the guard refuses quoted file
references (delfin/orca_knowledge/probe.py), and the Hessian-reuse path is
covered by the input-build test; parity is asserted here only for the parts
ORCA's parser can see without a Hessian on disk.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.orca_knowledge.probe import probe_input  # noqa: E402
from delfin.tools import _registry  # noqa: E402

_HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"


def test_opi_orca_th_reeval_input_accepted_by_orca_parser(tmp_path):
    out = tmp_path / "freq.out"
    out.write_text(
        _HEADER + "\n"
        "------------------\nTHERMOCHEMISTRY AT 298.15K\n------------------\n"
        "Temperature                         ...     298.15 K\n"
        "Pressure                            ...      1.00 atm\n"
        "Electronic energy                    ...    -151.55083691 Eh\n"
        "Zero point energy                    ...       0.02631509 Eh      16.51 kcal/mol\n"
        "Total Enthalpy                       ...    -151.52034633 Eh\n"
        "Vibrational entropy computed according to the QRRHO of S. Grimme "
        "(The Journal of Chemical Physics, 2012, 136, 9955)\n"
        "Total entropy correction             ...       0.02583637 Eh      16.21 kcal/mol\n"
        "Final Gibbs free energy              ...    -151.54618270 Eh\n"
    )
    inp = tmp_path / "reval.inp"
    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=500.0,
        re_eval_pressure=2.0,
        input_path=inp,
    )
    assert result.ok
    probe = probe_input(inp.read_text(), timeout_s=300)
    if probe["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert probe["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the orca_th_opi re-evaluation input; messages: "
        f"{probe['messages']}"
    )


def test_opi_orca_th_reeval_input_carries_temp_and_pressure(tmp_path):
    """The re-eval input carries the requested T/p with ORCA's own spelling."""
    out = tmp_path / "freq.out"
    out.write_text(
        _HEADER + "\n"
        "------------------\nTHERMOCHEMISTRY AT 298.15K\n------------------\n"
        "Electronic energy                    ...    -151.55083691 Eh\n"
        "Zero point energy                    ...       0.02631509 Eh      16.51 kcal/mol\n"
        "Total Enthalpy                       ...    -151.52034633 Eh\n"
        "Vibrational entropy computed according to the QRRHO of S. Grimme\n"
        "Total entropy correction             ...       0.02583637 Eh      16.21 kcal/mol\n"
        "Final Gibbs free energy              ...    -151.54618270 Eh\n"
    )
    inp = tmp_path / "reval.inp"
    result = _registry.get("orca_th_opi").execute(
        tmp_path,
        output_file=str(out),
        re_eval_temp=350.0,
        input_path=inp,
    )
    assert result.ok
    text = inp.read_text()
    assert "freq" in text.lower()
    assert "temp" in text.lower()
    assert "350.0" in text
