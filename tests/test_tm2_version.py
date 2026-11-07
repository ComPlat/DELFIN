"""Tests for turbomole_version.py (package TM2, phase 1 -- version safety).

Covers: version parsing from the three never-assume sources (TURBODIR,
module name, binary output), the unknown-version rule (None, no guess), the
MCP-target rule (8.0 for 7.x / unknown, exact for 8.x), the control-directive
compatibility check (unknown directive refused, 8.x-gated directive refused
on 7.x, clean control passes), and the provenance stamp.

No live TURBOMOLE is run: every version source is injected as a string
fixture, as the TM2 package requires.
"""

import importlib


def _mod():
    return importlib.import_module("delfin.turbomole_version")


# --------------------------------------------------------------------------
# Version parsing
# --------------------------------------------------------------------------

def test_turbodir_parses_version_from_path():
    m = _mod()
    v = m.version_from_turbodir({"TURBODIR": "/opt/turbomole/7.9/"})
    assert v is not None
    assert (v.major, v.minor) == (7, 9)
    assert v.source == "TURBODIR"
    assert str(v) == "7.9"


def test_turbodir_variant_separators():
    m = _mod()
    for raw in [
        "/opt/turbomole/8.0",
        "/opt/turbomole-v8.1",
        "/apps/turbomole_8.1",
    ]:
        v = m.version_from_turbodir({"TURBODIR": raw})
        assert v is not None, raw
        assert v.major == 8


def test_turbodir_absent_returns_none_not_guess():
    m = _mod()
    assert m.version_from_turbodir({}) is None
    # A path with no digit-dot-digit version is unparseable -> None (no guess).
    assert m.version_from_turbodir({"TURBODIR": "/opt/turbomole/"}) is None


def test_module_name_parses_version():
    m = _mod()
    v = m.version_from_module_name("chem/turbomole/7.9")
    assert v is not None
    assert str(v) == "7.9"
    assert v.source == "module"
    assert m.version_from_module_name("chem/turbomole") is None
    assert m.version_from_module_name(None) is None


def test_binary_output_parses_version():
    m = _mod()
    out = """\
        TURBOMOLE V7.9   |  copyright (c) 2022
        some build line
        ridft version 7.9
    """
    v = m.version_from_binary_output(out)
    assert v is not None
    assert str(v) == "7.9"
    assert v.source == "binary"
    assert m.version_from_binary_output("") is None
    assert m.version_from_binary_output(None) is None


def test_detect_precedence_turbodir_over_module_over_binary():
    m = _mod()
    # TURBODIR wins even when module/binary disagree.
    v = m.detect_installed_version(
        turbodir="/opt/turbomole/7.9",
        module_name="chem/turbomole/8.0",
        binary_output="TURBOMOLE V8.1",
    )
    assert v is not None and str(v) == "7.9"

    # No TURBODIR -> module wins over binary.
    v = m.detect_installed_version(
        module_name="chem/turbomole/8.0", binary_output="TURBOMOLE V8.1"
    )
    assert v is not None and str(v) == "8.0"

    # Only binary source.
    v = m.detect_installed_version(binary_output="TURBOMOLE V8.1")
    assert v is not None and str(v) == "8.1"


def test_detect_all_sources_absent_returns_none():
    m = _mod()
    assert m.detect_installed_version() is None
    assert (
        m.detect_installed_version(
            turbodir="/opt/turbomole/", module_name="chem/turbomole", binary_output="nothing"
        )
        is None
    )


def test_version_ordering():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    v80 = m.TurbomoleVersion(8, 0, "test", "8.0")
    v81 = m.TurbomoleVersion(8, 1, "test", "8.1")
    assert v79 < v80 < v81
    assert v80 >= v80
    assert v79.is_7x_line
    assert not v80.is_7x_line


# --------------------------------------------------------------------------
# MCP target rule
# --------------------------------------------------------------------------

def test_mcp_target_7x_asks_80_and_flags_differs():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    version, differs = m.mcp_target_version(v79)
    assert version == "8.0"
    assert differs is True


def test_mcp_target_unknown_is_conservative_80_differs():
    m = _mod()
    version, differs = m.mcp_target_version(None)
    assert version == "8.0"
    assert differs is True


def test_mcp_target_8x_matches_installed_no_differs():
    m = _mod()
    v80 = m.TurbomoleVersion(8, 0, "test", "8.0")
    v81 = m.TurbomoleVersion(8, 1, "test", "8.1")
    assert m.mcp_target_version(v80) == ("8.0", False)
    assert m.mcp_target_version(v81) == ("8.1", False)


# --------------------------------------------------------------------------
# Control-directive compatibility
# --------------------------------------------------------------------------

A_CLEAN_CONTROL = """\
$title test molecule
$operating system unix
$symmetry c1
$coord    file=coord
$atoms
    basis=def2-SVP
$charge 0
$dft
    functional b3-lyp
$ricore 2000
$maxcor 8000
$scfconv 7
$ridft
$end
"""


def test_clean_control_passes_on_7x():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    res = m.check_control(A_CLEAN_CONTROL, v79)
    assert res.ok, res.violations


def test_unknown_directive_is_refused():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    res = m.check_control("$title x\n$totally_bogus_directive 1\n$end\n", v79)
    assert not res.ok
    names = [v.directive for v in res.violations]
    assert "totally_bogus_directive" in names
    assert any("not a known TURBOMOLE control directive" in v.message
               for v in res.violations)


def test_unknown_directive_refused_even_when_installed_unknown():
    m = _mod()
    # Unknown installed version must NOT silently pass an unknown directive.
    res = m.check_control("$title x\n$bogus 1\n$end\n", None)
    assert not res.ok
    assert "bogus" in [v.directive for v in res.violations]


def test_8x_gated_directive_refused_on_7x():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    res = m.check_control("$title x\n$cosmo_rs\n$end\n", v79)
    assert not res.ok
    v = [x for x in res.violations if x.directive == "cosmo_rs"]
    assert v and "requires the 8.x line" in v[0].message


def test_8x_gated_directive_allowed_on_8x():
    m = _mod()
    v80 = m.TurbomoleVersion(8, 0, "test", "8.0")
    res = m.check_control("$title x\n$cosmo_rs\n$end\n", v80)
    assert res.ok, res.violations


def test_block_lines_are_not_mistaken_for_directives():
    m = _mod()
    # Lines inside a block carry no leading '$' and must not be read as
    # directives (e.g. the basis= line under $atoms).
    control = "$atoms\nbasis=def2-SVP charge=0\n$end\n"
    res = m.check_control(control, None)
    assert res.ok


def test_known_directives_injectable():
    m = _mod()
    v79 = m.TurbomoleVersion(7, 9, "test", "7.9")
    # The caller may restrict the vocabulary (e.g. to what a specific
    # generator emits); a directive outside it is then refused.
    res = m.check_control("$title x\n$cosmo\n$end\n", v79, known_directives=["title", "end"])
    assert not res.ok
    assert "cosmo" in [v.directive for v in res.violations]


# --------------------------------------------------------------------------
# Provenance stamp
# --------------------------------------------------------------------------

def test_stamp_appends_generator_version():
    m = _mod()
    stamped = m.stamp_generator_version("$title x\n$end\n", "delfin-tm", "1.2.0")
    assert "# generated by delfin-tm 1.2.0" in stamped
    # Idempotent: repeating does not duplicate.
    again = m.stamp_generator_version(stamped, "delfin-tm", "1.2.0")
    assert again.count("generated by delfin-tm 1.2.0") == 1
