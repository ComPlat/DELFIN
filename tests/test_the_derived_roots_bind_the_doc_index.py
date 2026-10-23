"""The derived containment must bind the doc index, by name.

``delfin_roots`` keeps all of ``~/.delfin`` out because
``credentials.json`` lives there, and binds the two sub-directories the
tools server persists into by name. The doc index
(``~/.delfin/doc_index.json``) is the ops server's working data —
``check_orca_manual_indexed`` reads it and ``index_new_pdf`` rewrites it —
and it was on the wrong side of that line: a contained ops server answered
``indexed: false`` on an account with a healthy index. Measured against
this session's own server (walled vs loose), 2026-09-21.

The fix is in the derivation, not in the server: bind the file the reader
resolves — ``get_default_index_path()``, the same authority the server
uses — writable, and nothing else from ``~/.delfin``.
"""
import shutil

import pytest

from delfin.agent import mcp_isolation


def _make_home(tmp_path, *, with_index=True):
    home = tmp_path / "home"
    dot = home / ".delfin"
    dot.mkdir(parents=True)
    (dot / "credentials.json").write_text("{}")  # must stay invisible
    if with_index:
        (dot / "doc_index.json").write_text("{}")
    return home


def test_the_derived_roots_bind_the_doc_index_by_name(tmp_path):
    home = _make_home(tmp_path)
    iso = mcp_isolation.delfin_roots(workspace=tmp_path, home=home)
    assert iso is not None
    index = str((home / ".delfin" / "doc_index.json").resolve())
    assert index in iso.write_roots


def test_the_doc_index_is_writable_index_new_pdf_rebuilds_it(tmp_path):
    home = _make_home(tmp_path)
    iso = mcp_isolation.delfin_roots(workspace=tmp_path, home=home)
    index = str((home / ".delfin" / "doc_index.json").resolve())
    assert index in iso.roots and index not in iso.read_roots


def test_without_an_index_no_file_root_is_invented(tmp_path):
    home = _make_home(tmp_path, with_index=False)
    iso = mcp_isolation.delfin_roots(workspace=tmp_path, home=home)
    assert iso is not None
    assert not any(p.endswith("doc_index.json") for p in iso.roots)


def test_credentials_stay_outside_the_derived_roots(tmp_path):
    home = _make_home(tmp_path)
    iso = mcp_isolation.delfin_roots(workspace=tmp_path, home=home)
    dot = str((home / ".delfin").resolve())
    assert dot not in iso.roots
    assert str(home.resolve()) not in iso.roots
    assert not any(p.endswith("credentials.json") for p in iso.roots)


@pytest.mark.skipif(shutil.which("bwrap") is None, reason="no bwrap")
def test_the_launch_argv_binds_the_index_file(tmp_path):
    """The bind survives into bwrap's argv: a file target under the
    blanked ``$HOME`` tmpfs needs an explicit ``--bind`` or the server
    sees an empty directory where its index should be."""
    home = _make_home(tmp_path)
    iso = mcp_isolation.delfin_roots(workspace=tmp_path, home=home)
    index = str((home / ".delfin" / "doc_index.json").resolve())
    argv = mcp_isolation.bwrap_argv("/bin/true", [], iso, home=home)
    i = argv.index("--bind") if "--bind" in argv else -1
    assert index in argv
    # ... as a --bind (rw), not an --ro-bind
    for flag in ("--bind",):
        for k in range(len(argv) - 1):
            if argv[k] == flag and argv[k + 1] == index:
                break
        else:
            pytest.fail(f"{index} is not bound rw via {flag}")
