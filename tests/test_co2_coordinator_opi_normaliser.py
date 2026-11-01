"""Unit tests for the CO2-coordinator OPI input normaliser.

``delfin.tools.co2_coordinator_opi.normalise_inp`` is pure text parsing:
it turns both the classic CO2-Coordinator ``.inp`` and the OPI-built one
into the same canonical structure so parity can be asserted without
string comparison.  These tests run with the py3.11 CI interpreter that
has NO OPI installed, so they must not import OPI anywhere.
"""

from __future__ import annotations

from delfin.tools import co2_coordinator_opi as m


def normalise(text: str) -> dict:
    return m.normalise_inp(text)


# ---------------------------------------------------------------------------
# simple keywords
# ---------------------------------------------------------------------------


def test_tokens_from_classic_one_line():
    inp = "! PBE0 def2-SVP def2/J D4 CPCM(DMF) SP\n"
    d = normalise(inp)
    assert d["tokens"] == frozenset({"pbe0", "def2-svp", "def2/j", "d4",
                                     "cpcm(dmf)", "sp"})


def test_tokens_from_opi_multi_line_case_insensitive():
    # OPI writes each keyword on its own '!' line, lowercased (with
    # 'cpcm(dmf)' and 'ri-socon' surviving as bare strings).
    inp = "!pbe0\n!def2-svp\n!def2/j\n!d4\n!cpcm(dmf)\n!sp\n"
    d = normalise(inp)
    assert d["tokens"] == frozenset({"pbe0", "def2-svp", "def2/j", "d4",
                                     "cpcm(dmf)", "sp"})


def test_tokens_ignore_comment_lines():
    inp = "! pbe0 opt\n# a comment that starts with ! is unusual\n"
    d = normalise(inp)
    assert "pbe0" in d["tokens"] and "opt" in d["tokens"]


# ---------------------------------------------------------------------------
# %maxcore / %pal  (classic inline vs OPI block)
# ---------------------------------------------------------------------------


def test_maxcore_parsed():
    assert normalise("%maxcore 2000\n")["maxcore"] == 2000


def test_pal_inline_classic_and_block_opi_equal():
    classic = "%pal nprocs 4 end\n"
    opi = "%pal\n    nprocs 4\nend\n\n"
    assert normalise(classic)["blocks"]["pal"]["nprocs"] == 4
    assert normalise(opi)["blocks"]["pal"]["nprocs"] == 4


# ---------------------------------------------------------------------------
# %scf (MaxIter / BrokenSym), classic and OPI forms
# ---------------------------------------------------------------------------


def test_scf_maxiter_both_forms():
    classic = "%scf\n  MaxIter 375\nend\n"
    opi = "%scf\n    maxiter 375\nend\n\n"
    assert normalise(classic)["blocks"]["scf"] == {"maxiter": 375}
    assert normalise(opi)["blocks"]["scf"] == {"maxiter": 375}


def test_scf_brokensym_and_maxiter():
    classic = "%scf\n  BrokenSym 1,5\n  MaxIter 375\nend\n"
    opi = "%scf\n    maxiter 375\n    brokensym 1,5\nend\n\n"
    expected = {"maxiter": 375, "brokensym": (1, 5)}
    assert normalise(classic)["blocks"]["scf"] == expected
    assert normalise(opi)["blocks"]["scf"] == expected


# ---------------------------------------------------------------------------
# %geom Scan, classic multi-line vs OPI inline
# ---------------------------------------------------------------------------


def test_geom_scan_classic_and_opi_equal():
    classic = (
        "%geom\n"
        "  MaxIter 400\n"
        "  Scan\n"
        "    B  0  4 = 5.00, 1.7, 15\n"
        "  end\n"
        "end\n"
    )
    opi = (
        "%geom\n"
        "    maxiter 400\n"
        "    scan B 0 4 = 5.0, 1.7, 15 end\n"
        "end\n\n"
    )
    ce = normalise(classic)["blocks"]["geom"]
    oe = normalise(opi)["blocks"]["geom"]
    assert ce == oe
    assert ce["maxiter"] == 400
    assert ce["scan"]["mode"] == "B"
    assert ce["scan"]["atoms"] == (0, 4)
    assert ce["scan"]["start"] == 5.0
    assert ce["scan"]["end"] == 1.7
    assert ce["scan"]["steps"] == 15


# ---------------------------------------------------------------------------
# %qmmm QMATOMS: classic range {0:3} vs OPI enumerated { 0 1 2 3 }
# ---------------------------------------------------------------------------


def test_qmmm_classic_range_and_opi_set_equal():
    classic = "%QMMM\n  QMATOMS {0:3}\nEND\nEND\n"
    opi = "%qmmm\n    qmatoms { 0 1 2 3 } end\nend\n\n"
    assert normalise(classic)["blocks"]["qmmm"]["qmatoms"] == normalise(
        opi
    )["blocks"]["qmmm"]["qmatoms"]
    assert normalise(classic)["blocks"]["qmmm"]["qmatoms"] == frozenset(
        {0, 1, 2, 3}
    )


# ---------------------------------------------------------------------------
# coordinate block
# ---------------------------------------------------------------------------


def test_coords_charge_mult_and_newgto():
    inp = (
        "* xyz -2 3\n"
        "Ni       0.00000000     0.00000000     0.00000000\n"
        'C        0.00000000     0.00000000     3.00000000   NewGTO "def2-TZVP" end\n'
        "*\n"
    )
    d = normalise(inp)
    assert d["charge"] == -2
    assert d["multiplicity"] == 3
    assert d["coords"][0] == ("Ni", (0.0, 0.0, 0.0), None)
    assert d["coords"][1] == ("C", (0.0, 0.0, 3.0), "def2-TZVP")


def test_coords_larger_precision_opi_form():
    # OPI writes 16 decimals; classic 8.  Values must compare equal.
    classic = (
        "* xyz -2 3\n"
        "Ni       0.00000000     0.00000000     0.00000000\n"
        "C        0.00000000     0.00000000     4.16000000\n"
        "*\n"
    )
    opi = (
        "* xyz -2 3\n"
        "Ni             0.0000000000000000             0.0000000000000000             0.0000000000000000\n"
        "C              0.0000000000000000             0.0000000000000000             4.1600000000000001\n"
        "*\n"
    )
    assert m.inputs_equivalent(normalise(classic), normalise(opi)) == []


# ---------------------------------------------------------------------------
# inputs_equivalent
# ---------------------------------------------------------------------------


def test_equivalent_returns_empty_for_same():
    a = normalise("! PBE0 def2-SVP SP\n%maxcore 2000\n%pal nprocs 4 end\n* xyz -2 3\nNi 0 0 0\n*\n")
    b = normalise("!pbe0\n!def2-svp\n!sp\n%maxcore 2000\n%pal\n nprocs 4\nend\n* xyz -2 3\nNi 0.0000000000000000 0.0000000000000000 0.0000000000000000\n*\n")
    assert m.inputs_equivalent(a, b) == []


def test_equivalent_reports_token_diff():
    a = normalise("! PBE0 def2-SVP SP\n* xyz 0 1\nNi 0 0 0\n*\n")
    b = normalise("! PBE0 def2-SVP OPT\n* xyz 0 1\nNi 0 0 0\n*\n")
    diffs = m.inputs_equivalent(a, b)
    assert any("tokens differ" in x for x in diffs)


def test_equivalent_reports_charge_diff():
    a = normalise("* xyz -2 3\nNi 0 0 0\n*\n")
    b = normalise("* xyz 0 1\nNi 0 0 0\n*\n")
    diffs = m.inputs_equivalent(a, b)
    assert any("charge" in x for x in diffs)


def test_equivalent_reports_newgto_diff():
    a = normalise('* xyz 0 1\nC 0 0 0 NewGTO "def2-TZVP" end\n*\n')
    b = normalise('* xyz 0 1\nC 0 0 0 NewGTO "def2-SVP" end\n*\n')
    diffs = m.inputs_equivalent(a, b)
    assert any("NewGTO" in x for x in diffs)
