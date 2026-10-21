"""Every file-upload widget of the dashboard writes under the same rules.

The Calculations browser's uploads went through ``upload_limits`` already;
the other three ``FileUpload`` widgets did not. Driven through the real
tabs on tmp directories (the ORCA Builder and the Literature tab; the
agent tab's handler is checked in its source and through the primitive it
calls), with the widget's value set the way the browser front end sets it.

Asserted for each:

* a file over ``DELFIN_UPLOAD_MAX_BYTES`` is not written (the others of the
  same drop are);
* what is written is owner-only (0600);
* the name the browser sends is a file name, not a path: ``../x`` lands
  inside the destination or not at all;
* a link that already sits at the destination is replaced, never written
  through (``upload_limits.write_replacing``).
"""

from __future__ import annotations

import datetime
import pathlib
import stat

import pytest

from delfin.dashboard import upload_limits

WATER = "O 0.000 0.000 0.117\nH 0.000 0.757 -0.467\nH 0.000 -0.757 -0.467\n"


def _walk(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _walk(child)


def _entry(name: str, content: bytes) -> dict:
    return {"name": name, "type": "", "size": len(content),
            "content": memoryview(content),
            "last_modified": datetime.datetime.now(datetime.timezone.utc)}


def _mode(path: pathlib.Path) -> int:
    return stat.S_IMODE(path.lstat().st_mode)


# -- the primitive ----------------------------------------------------------------

def test_write_replacing_does_not_follow_a_planted_link(tmp_path):
    outside = tmp_path / "outside.txt"
    outside.write_text("keep", encoding="utf-8")
    dest = tmp_path / "in" / "file.gbw"
    dest.parent.mkdir()
    dest.symlink_to(outside)
    upload_limits.write_replacing(dest, b"uploaded")
    assert outside.read_text(encoding="utf-8") == "keep"
    assert not dest.is_symlink() and dest.read_bytes() == b"uploaded"
    assert _mode(dest) == 0o600


def test_write_replacing_replaces_a_file_and_refuses_a_directory(tmp_path):
    dest = tmp_path / "file.xyz"
    dest.write_bytes(b"old")
    upload_limits.write_replacing(dest, b"new")
    assert dest.read_bytes() == b"new" and _mode(dest) == 0o600
    (tmp_path / "adir").mkdir()
    with pytest.raises(OSError):
        upload_limits.write_replacing(tmp_path / "adir", b"x")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["adir", "file.xyz"]


def test_write_replacing_refuses_over_the_limit(tmp_path, monkeypatch):
    monkeypatch.setenv(upload_limits.UPLOAD_ENV, "10")
    with pytest.raises(upload_limits.UploadRefused):
        upload_limits.write_replacing(tmp_path / "big", b"x" * 11)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name,expected", [
    ("plain.gbw", "plain.gbw"), ("../up.gbw", "up.gbw"),
    ("/abs/path/abs.xyz", "abs.xyz"), ("dir\\win.hess", "win.hess"),
])
def test_an_upload_name_is_a_file_name(name, expected):
    assert upload_limits.upload_name(name) == expected


@pytest.mark.parametrize("name", ["", ".", "..", "../", "a\0b"])
def test_an_upload_name_that_names_no_file_is_refused(name):
    with pytest.raises(upload_limits.UploadRefused):
        upload_limits.upload_name(name)


# -- the ORCA Builder -------------------------------------------------------------

@pytest.fixture
def orca(tmp_path, monkeypatch):
    widgets = pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_orca_builder
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda *a, **k: None
    tab, refs = tab_orca_builder.create_tab(ctx)
    found = list(_walk(tab))
    upload = next(w for w in found if isinstance(w, widgets.FileUpload))
    job = next(w for w in found if isinstance(w, widgets.Text)
               and getattr(w, "description", "") == "Job Name:")
    refs["orca_coords"].value = WATER
    job.value = "upl"

    def save(*entries):
        upload.value = tuple(entries)
        refs["orca_save_btn"].click()
        return ctx.calc_dir / "upl"

    return {"tmp": tmp_path, "save": save, "calc": ctx.calc_dir}


def test_orca_builder_does_not_write_an_upload_over_the_limit(orca, monkeypatch):
    monkeypatch.setenv(upload_limits.UPLOAD_ENV, "1000")
    job_dir = orca["save"](_entry("big.gbw", b"x" * 2000),
                           _entry("small.xyz", b"3\n\n" + WATER.encode()))
    assert (job_dir / "small.xyz").is_file()
    assert not (job_dir / "big.gbw").exists()


def test_orca_builder_uploads_are_owner_only(orca):
    job_dir = orca["save"](_entry("guess.gbw", b"\0" * 64))
    assert _mode(job_dir / "guess.gbw") == 0o600


def test_orca_builder_upload_name_cannot_leave_the_job_directory(orca):
    job_dir = orca["save"](_entry("../escape.gbw", b"\0" * 8))
    assert not (orca["calc"] / "escape.gbw").exists()
    assert (job_dir / "escape.gbw").is_file()


def test_orca_builder_does_not_write_through_a_planted_link(orca):
    outside = orca["tmp"] / "outside.txt"
    outside.write_text("keep", encoding="utf-8")
    job_dir = orca["calc"] / "upl"
    job_dir.mkdir()
    (job_dir / "guess.gbw").symlink_to(outside)
    orca["save"](_entry("guess.gbw", b"\0" * 8))
    assert outside.read_text(encoding="utf-8") == "keep"
    assert not (job_dir / "guess.gbw").is_symlink()


# -- the Literature tab -----------------------------------------------------------

@pytest.fixture
def literature(tmp_path, monkeypatch):
    widgets = pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_literature
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office", "home", "repo"):
        (tmp_path / name).mkdir()
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.repo_dir = tmp_path / "repo"
    ctx.run_js = lambda *a, **k: None
    tab, _refs = tab_literature.create_tab(ctx)
    upload = next(w for w in _walk(tab) if isinstance(w, widgets.FileUpload))

    def send(*entries):
        upload.value = tuple(entries)
        return tmp_path / "repo" / "literature"

    return send


def test_literature_does_not_write_an_upload_over_the_limit(literature, monkeypatch):
    monkeypatch.setenv(upload_limits.UPLOAD_ENV, "1000")
    lit = literature(_entry("big.pdf", b"%PDF" + b"x" * 2000),
                     _entry("small.pdf", b"%PDF-1.4\n%%EOF\n"))
    assert (lit / "small.pdf").is_file()
    assert not (lit / "big.pdf").exists()


def test_literature_uploads_are_owner_only(literature):
    lit = literature(_entry("paper.pdf", b"%PDF-1.4\n%%EOF\n"))
    assert _mode(lit / "paper.pdf") == 0o600


# -- the agent tab ----------------------------------------------------------------

def _agent_block(start: str, end: str) -> str:
    from delfin.dashboard import tab_agent

    source = pathlib.Path(tab_agent.__file__).read_text(encoding="utf-8")
    return source[source.index(start):source.index(end)]


def test_the_agent_upload_cap_follows_the_upload_limit():
    handler = _agent_block("def _on_image_upload", "def _materialise_uploads")
    assert "max_upload_bytes()" in handler


def test_the_agent_upload_is_written_by_the_shared_primitive():
    writer = _agent_block("def _materialise_uploads",
                          'image_upload.observe(_on_image_upload')
    code = "\n".join(line.split("#", 1)[0] for line in writer.splitlines())
    assert "write_replacing(" in code
    assert ".write_bytes(" not in code


def test_no_upload_handler_writes_bytes_directly():
    """The tabs that hold a FileUpload write uploads only through
    ``upload_limits``."""
    from delfin.dashboard import tab_orca_builder, tab_literature

    for module, start, end in (
            (tab_orca_builder, "def save_uploaded_files", "def reset_orca_builder"),
            (tab_literature, "def on_upload(change)", "def on_delete_click")):
        source = pathlib.Path(module.__file__).read_text(encoding="utf-8")
        block = source[source.index(start):source.index(end)]
        assert ".write_bytes(" not in block, module.__name__
