"""The DLPNO step consumes the chain's ``moread`` wiring.

``delfin.tools.templates._opi_chain`` promises the DLPNO single point is
"started from the chain's orbitals (MOREAD)": the engine's ``gbw -> moread``
wiring rule (``delfin.tools._wiring``, ``inject_kwarg="moread"``) injects the
upstream ``.gbw`` into the step's kwargs, and ``orca_opt_freq_opi`` forwards
the kwarg to ``delfin.tools.opi_chain.run_opi_job`` (manual ground truth:
ORCA 6.1 §2.20.5/§2.20.10 Table 2.68 -- ``MORead`` reads MOs from a previous
calculation via ``%moinp``; ``AutoStart`` only covers same-basename files).

``OrcaDlpnoOpiAdapter.execute`` read the injected ``moread`` kwarg and
forwarded it to ``build_dlpno_input`` / ``run_dlpno``, which wrote ``!
moread`` + ``%moinp`` the same way ``run_opi_job`` does.  These tests pin
that path.  They need no ``opi`` import (the assertions read signatures and
source, never build an input), so they run on every CI job.
"""

from __future__ import annotations

import inspect

from delfin.tools.adapters import orca_dlpno_opi as dlpno
from delfin.tools._wiring import autowire_kwargs


def test_wiring_injects_moread_into_the_dlpno_step():
    # Control (green on both sides): the wiring rule exists and injects
    # ``moread`` for an upstream ``gbw`` -- this is the contract the adapter
    # must satisfy.  It is green before and after the fix; the red/green
    # difference is the next test.
    kw = autowire_kwargs("orca_dlpno_opi", {}, {"gbw": "chain/dft_opt.gbw"})
    assert kw.get("moread") == "chain/dft_opt.gbw"


def test_dlpno_execute_reads_the_moread_kwarg():
    # The fix (red on the baseline): ``execute`` must read ``moread`` --
    # the injected kwarg the chain template promises the DLPNO step
    # consumes.  On the baseline this fails: ``kwargs.get("moread")``
    # appears nowhere in ``execute``'s source, and ``build_dlpno_input``
    # has no ``moread`` parameter, so the injected value is silently
    # dropped and the step runs from a fresh SCF guess.
    src = inspect.getsource(dlpno.OrcaDlpnoOpiAdapter.execute)
    assert "moread" in src, (
        "OrcaDlpnoOpiAdapter.execute ignores the injected 'moread' kwarg; "
        "the chain's gbw->moread wiring is dropped and the DLPNO step "
        "starts from a fresh SCF guess"
    )


def test_build_dlpno_input_takes_and_writes_moread():
    # The fix, continued: ``build_dlpno_input`` accepts ``moread`` and its
    # body writes the manual's MORead form (``! moread`` + ``%moinp``),
    # matching ``run_opi_job``'s application of the same kwarg.
    sig = inspect.signature(dlpno.build_dlpno_input)
    assert "moread" in sig.parameters, (
        "build_dlpno_input has no 'moread' parameter; the chain cannot "
        "start the DLPNO step from the wired orbitals"
    )
    body = inspect.getsource(dlpno.build_dlpno_input)
    for marker in ("moread", "moinp"):
        assert marker in body, f"build_dlpno_input never writes {marker}"


def test_execute_forwards_moread_to_the_runner():
    # ``execute`` passes the kwarg through ``run_dlpno``, whose signature
    # carries ``moread`` and forwards it to ``build_dlpno_input`` -- the
    # same surface ``orca_opt_freq_opi`` uses for its ``moread``.
    runner_sig = inspect.signature(dlpno.run_dlpno)
    assert "moread" in runner_sig.parameters
    src = inspect.getsource(dlpno.OrcaDlpnoOpiAdapter.execute)
    assert "moread" in src
