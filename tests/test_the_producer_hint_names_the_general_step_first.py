"""A missing input was answered with "orca_bs_opi, orca_dlpno_opi, ..."

When a spec step consumes a capability nothing upstream produces, the
validator names the steps that could produce it. It showed the first six
in ALPHABETICAL order and cut the rest with an ellipsis. Measured
2026-11-23 for ``qc_output``: 27 registered producers, and the six shown
were orca_bs_opi, orca_dlpno_opi, orca_eprnmr_opi, orca_freq,
orca_freq_opi, orca_mlip_opt_opi -- while orca_sp, the plain single
point that needs only a charge, was thirteenth and invisible.

It also made the suite say different things on different machines: the
OPI variants register where they can, so on a host with them orca_sp fell
past the cut, and test_platform's test of exactly this hint failed here
and passed in CI.

The order is now GENERAL BEFORE SPECIALISED, read from the registry's own
naming: a variant extends its base name (orca_sp -> orca_sp_opi ->
orca_sp_opi_basis_solv), so fewer name parts is the more general step.
Ties break alphabetically, and the count of the rest is stated instead of
an ellipsis that hides how many there are.
"""
from __future__ import annotations

from delfin.tools import _validate as V


def test_the_general_producers_come_first():
    shown = V._producers_of("qc_output")
    assert shown[:3] == sorted(shown[:3])
    first = shown[:6]
    assert "orca_sp" in first, first
    assert not any(n.endswith("_opi") for n in first), first


def test_the_hint_is_the_same_with_or_without_the_variants(monkeypatch):
    """The host dependence, driven directly: remove every OPI variant from
    the registry and the general producers shown must not change."""
    full = V._list_steps()
    lean = {k: v for k, v in full.items() if "_opi" not in k}
    with_variants = V._producers_of("qc_output")[:6]
    monkeypatch.setattr(V, "_list_steps", lambda: lean)
    without = V._producers_of("qc_output")[:6]
    assert with_variants == without


def test_the_rest_is_counted_not_hidden():
    from delfin.tools import platform
    rep = platform.validate_spec(
        {"name": "m", "steps": [{"step": "cclib_parse"}]}, geometry=False)
    msgs = [m for d in rep.diagnostics for m in d.messages]
    hint = next(m for m in msgs if "producing it" in m)
    assert "orca_sp" in hint
    assert "…" not in hint
    assert "more" in hint
