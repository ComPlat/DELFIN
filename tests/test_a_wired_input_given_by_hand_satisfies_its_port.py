"""A consumed port is satisfied by an upstream producer OR by the parameter
the port wires to, given by hand.

``consumes`` names the ports a step needs; ``wires`` maps a port onto the
keyword that receives it. The static validator accepted a port only from an
upstream step, so a step whose input is a file the person already has could
not stand first: ``cclib_parse`` with ``filepath=run.out`` was refused with
"missing input capability: qc_output" (measured 2026-11-13, also
``censo_sort`` with an ensemble file and ``biostructure_summary`` with a
structure file). Through the MCP server the one-step application failed
static validation and never ran.

The parameter is the same thing the wire would deliver, so giving it is
giving the port. A step's own contract default does not count: a default
is not something anybody supplied.
"""

import pytest

from delfin.tools import platform


def _ok(step, **kwargs):
    rep = platform.validate_spec({"name": "t", "steps": [dict(step=step, **kwargs)]})
    messages = [m for d in rep.diagnostics for m in d.messages]
    return rep.ok, messages


@pytest.mark.parametrize("step, kwargs", [
    ("cclib_parse", {"filepath": "run.out"}),
    ("censo_sort", {"ensemble": "ensemble.xyz"}),
    ("biostructure_summary", {"structure_file": "receptor.pdb"}),
])
def test_the_wired_parameter_given_by_hand_satisfies_the_port(step, kwargs):
    ok, messages = _ok(step, **kwargs)
    assert ok, messages
    assert not any("missing input capability" in m for m in messages)


@pytest.mark.parametrize("step", [
    "cclib_parse", "censo_sort", "biostructure_summary"])
def test_without_the_parameter_the_port_is_still_missing(step):
    """Nothing given and nothing upstream: the refusal stands, with its
    hint about which step could produce the port."""
    ok, messages = _ok(step)
    assert not ok
    assert any("missing input capability" in m for m in messages), messages


def test_an_empty_value_does_not_count_as_given():
    ok, messages = _ok("cclib_parse", filepath="")
    assert not ok
    assert any("missing input capability" in m for m in messages), messages


def test_a_port_without_a_wire_still_needs_an_upstream_step():
    """``geometry`` and the other unwired ports are unchanged: there is no
    parameter that stands for them."""
    ok, messages = _ok("ase_single_point")
    assert not ok
