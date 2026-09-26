"""Controls for shared_tmp: one temp directory for the cage and the
file tools.

Red on the previous commit: ``delfin.agent.shared_tmp`` does not exist
(Phase 4 of the LB assignment). Background: bash runs inside bwrap
with ``--tmpfs /tmp`` (mcp_isolation.bwrap_argv), a PRIVATE empty /tmp,
while read_file resolves absolute paths against the real filesystem --
so a file written to /tmp in bash is invisible to read_file, and two
sessions in Welle 4 asked for read access to all of /tmp to bridge
that. The design here: one temp directory per session under DELFIN's
state root, bind-mounted into the cage AS /tmp, and mapped back for
the file tools -- everything else (any other absolute path) is passed
through untouched.

The xfail(strict=True) integration test pins the handler contract: a
path under /tmp handed to read_file must resolve to the session's
shared temp directory instead of the host's /tmp.
"""

import pytest

from pathlib import Path

from delfin.agent import shared_tmp


@pytest.fixture
def base(tmp_path):
    d = tmp_path / "state"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# session_dir: where a session's shared temp directory lives
# ---------------------------------------------------------------------------

def test_session_dir_is_under_the_state_root_and_per_session(base):
    a = shared_tmp.session_dir("sess-a", base=base)
    b = shared_tmp.session_dir("sess-b", base=base)
    assert a.parent == base / "tmp"
    assert b.parent == base / "tmp"
    assert a != b


def test_session_dir_sanitizes_the_session_id(base):
    # A session id with slashes or dots must not escape the tmp root
    # ("../../ssh" is a plausible hostile input).
    d = shared_tmp.session_dir("../../../etc", base=base)
    assert d.parent == base / "tmp"
    assert ".." not in d.parts


def test_session_dir_is_deterministic(base):
    assert (shared_tmp.session_dir("s", base=base)
            == shared_tmp.session_dir("s", base=base))


# ---------------------------------------------------------------------------
# map_to_host: cage paths -> host paths
# ---------------------------------------------------------------------------

def test_tmp_path_maps_into_the_session_dir(base):
    sd = shared_tmp.session_dir("s", base=base)
    host = shared_tmp.map_to_host("/tmp/probe.py", session_dir=sd)
    assert host == sd / "probe.py"


def test_nested_tmp_path_maps_componentwise(base):
    sd = shared_tmp.session_dir("s", base=base)
    host = shared_tmp.map_to_host("/tmp/sub/dir/x.log", session_dir=sd)
    assert host == sd / "sub" / "dir" / "x.log"


def test_non_tmp_absolute_path_passes_through_unchanged(base):
    sd = shared_tmp.session_dir("s", base=base)
    assert shared_tmp.map_to_host("/etc/passwd", session_dir=sd) \
        == Path("/etc/passwd")
    assert shared_tmp.map_to_host("/home/u/ws/f.py", session_dir=sd) \
        == Path("/home/u/ws/f.py")


def test_escape_attempt_through_tmp_is_refused(base):
    # "/tmp/../etc/passwd" must not become <session_dir>/../etc/passwd;
    # normalization happens FIRST, and a normalized path that leaves
    # /tmp is either passed through as its host meaning (it IS the
    # host's /etc/passwd, outside the share) or refused -- never a
    # breakout out of the session dir.
    sd = shared_tmp.session_dir("s", base=base)
    host = shared_tmp.map_to_host("/tmp/../etc/passwd", session_dir=sd)
    assert host == Path("/etc/passwd")        # no <sd>/../ anywhere
    host2 = shared_tmp.map_to_host("/tmp/../../x", session_dir=sd)
    assert host2 == Path("/x") and str(sd) not in str(host2)


def test_relative_path_is_passed_through(base):
    sd = shared_tmp.session_dir("s", base=base)
    assert shared_tmp.map_to_host("plain.txt", session_dir=sd) \
        == Path("plain.txt")


def test_tmp_root_itself_maps_to_the_session_dir(base):
    sd = shared_tmp.session_dir("s", base=base)
    assert shared_tmp.map_to_host("/tmp", session_dir=sd) == sd


# ---------------------------------------------------------------------------
# cage_mount: what bwrap should bind instead of --tmpfs /tmp
# ---------------------------------------------------------------------------

def test_cage_mount_returns_the_dir_to_bind_as_tmp(base):
    sd = shared_tmp.session_dir("s", base=base)
    assert shared_tmp.cage_mount(session_dir=sd) == sd


def test_bwrap_pair_describes_the_substitution():
    old, new = shared_tmp.bwrap_substitution(session_dir="/x/y")
    assert old == ["--tmpfs", "/tmp"]
    assert new == ["--bind", "/x/y", "/tmp"]
