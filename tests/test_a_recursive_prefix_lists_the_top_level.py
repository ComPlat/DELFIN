"""A ``list_files`` pattern with a ``**/`` prefix lists what a glob means.

From a real run of the product on 2026-11-26 (w60 package r2): the model
called ``list_files(path="delfin/agent/pack/skills", pattern="**/*.md")``
on a directory holding 67 ``.md`` files and no subdirectories, and the
tool answered **"No files matching pattern."** It then had to reach for
``ls`` to discover the truth. The model's own summary of that run names
the cause exactly: *"the glob pattern, not the directory, was the
problem"*.

``fnmatch`` is not a glob: ``**`` has no special meaning there, so
``fnmatch("code-mode.md", "**/*.md")`` is ``False`` -- ``*`` does not
cross ``/``, and a top-level file has no ``/`` to cross. Meanwhile
``fnmatch("sub/x.md", "**/*.md")`` is ``True``. So the pattern matched
everywhere EXCEPT the top level of the directory the caller named, which
is the one place they could see.

The fix treats a leading ``**/`` as "this directory and below", which is
what every glob user means by it. A pattern without that prefix goes
through ``fnmatch`` unchanged, so nothing that worked before moves.

Red on the parent commit: the first test gets "No files matching
pattern" instead of the three top-level names. The control cases pass on
both -- they are the shapes that already worked, kept to prove the
prefix rule did not widen into "everything matches".
"""

import pytest

from delfin.agent.api_client import _DocToolExecutor


def _perms(ws):
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "glob-prefix-test"
    return perms


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    (d / "pkg").mkdir(parents=True)
    (d / "pkg" / "one.md").write_text("a", encoding="utf-8")
    (d / "pkg" / "two.md").write_text("b", encoding="utf-8")
    (d / "pkg" / "notes.txt").write_text("c", encoding="utf-8")
    (d / "pkg" / "sub").mkdir()
    (d / "pkg" / "sub" / "deep.md").write_text("d", encoding="utf-8")
    return d


def _list(ws, **arguments):
    return _DocToolExecutor()._execute_list_files(arguments, _perms(ws))


def test_a_recursive_prefix_lists_the_top_level_files_too(ws):
    """``**/*.md`` must find ``pkg/one.md``, not only ``pkg/sub/deep.md``."""
    out = _list(ws, path="pkg", pattern="**/*.md")
    assert "No files matching pattern" not in out, out
    assert "one.md" in out, out
    assert "two.md" in out, out
    assert "deep.md" in out, out
    # A non-matching suffix stays excluded -- the rule widens the DIRECTORY
    # a pattern applies to, it does not stop filtering by suffix.
    assert "notes.txt" not in out, out


def test_the_bare_suffix_still_lists_nested_files(ws):
    """Control: ``*.md`` already reached nested files and must still."""
    out = _list(ws, path="pkg", pattern="*.md")
    assert "one.md" in out, out
    assert "deep.md" in out, out
    assert "notes.txt" not in out, out


def test_a_star_lists_every_file(ws):
    """Control: ``*`` matched everything before and still does."""
    out = _list(ws, path="pkg", pattern="*")
    for name in ("one.md", "two.md", "notes.txt", "deep.md"):
        assert name in out, out


def test_a_missing_pattern_is_still_said_to_be_missing(ws):
    """Control: a genuinely absent pattern must not be papered over."""
    out = _list(ws, path="pkg", pattern="**/*.xyz")
    assert "No files matching pattern" in out, out
