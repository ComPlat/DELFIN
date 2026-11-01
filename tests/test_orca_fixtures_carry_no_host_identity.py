"""The committed ORCA output fixtures name no real machine, user or home.

ORCA writes the host name into the banner (``Host name:``) and MPI writes
``[<host>:<pid>]`` before every line of an abort trace. Fixtures copied from
a real run carry the build machine's name unless it is replaced. The
placeholder is ``fixture-host``; parsers do not read either field.
"""
import re
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "orca611"


def _outputs():
    return sorted(FIXTURES.glob("*/*.out"))


def test_there_are_fixtures_to_check():
    assert _outputs()


def test_the_banner_host_is_the_placeholder():
    for out in _outputs():
        for line in out.read_text(errors="replace").splitlines():
            if "Host name:" in line:
                assert line.split("Host name:", 1)[1].strip() == "fixture-host", out.name


def test_mpi_trace_lines_carry_the_placeholder_host():
    tag = re.compile(r"^\[([^:\]]+):\d+\]")
    for out in _outputs():
        for line in out.read_text(errors="replace").splitlines():
            m = tag.match(line)
            if m:
                assert m.group(1) == "fixture-host", out.name


def test_no_home_path_in_any_fixture_file():
    for f in sorted(FIXTURES.rglob("*")):
        if f.is_file() and f.suffix != ".gbw":
            assert "/home/" not in f.read_text(errors="replace"), f.name
