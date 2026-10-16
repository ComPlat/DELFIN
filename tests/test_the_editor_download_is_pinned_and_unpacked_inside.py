"""The drawing editor is one pinned build, checked by SHA-256, unpacked inside.

Input: ``delfin.dashboard.ketcher.pinned_release`` and ``install`` with the
download replaced by an in-memory archive (no network) and HOME in tmp_path.

Output: failures for
- a release that is asked of the network or not named by version, URL and
  SHA-256: the editor is served on the dashboard's origin, from a route
  without a token, so whoever controls "the newest release" decides what
  script runs beside the dashboard (ASVS V15, SSDF PW.4.4, CIS SSC 3.1);
- a download whose bytes differ from the pin, or that is larger than it,
  being installed;
- an archive member that names a path outside the unpacking directory --
  ``..`` components, an absolute path, a backslash -- or that is a
  symbolic link, being written, instead of the whole install being refused
  with nothing kept (fail closed);
- a clean archive at the pinned digest no longer installing.

The pin is ketcher-standalone-3.18.0.zip, 35 057 618 bytes, SHA-256
484e7f10...3e7e: the asset digest the GitHub release API reports for
epam/ketcher v3.18.0 (read 2026-10-16), recomputed over the downloaded
file with the same result.
"""

from __future__ import annotations

import hashlib
import io
import pathlib
import stat
import zipfile

import pytest

from delfin.dashboard import ketcher as kt

PAGE = '<script defer src="./static/js/main.abc.js"></script>'


@pytest.fixture
def home(tmp_path, monkeypatch):
    where = tmp_path / "home"
    where.mkdir()
    monkeypatch.setenv("HOME", str(where))
    monkeypatch.setattr(pathlib.Path, "home", lambda: where)
    monkeypatch.delenv(kt._URL_ENV, raising=False)
    monkeypatch.delenv("DELFIN_VOILA_ROOT_DIR", raising=False)
    return where


def _archive(*extra) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as bundle:
        bundle.writestr("index.html", PAGE)
        bundle.writestr("static/js/main.abc.js", "// the editor")
        for member, data in extra:
            bundle.writestr(member, data)
    return buf.getvalue()


def _serve(monkeypatch, payload: bytes, *, sha256: str = "", size: int = 0):
    release = {"ok": True, "version": "3.18.0",
               "url": "https://example.invalid/ketcher-standalone-3.18.0.zip",
               "size": size or len(payload),
               "sha256": sha256 or hashlib.sha256(payload).hexdigest(),
               "status": ""}
    monkeypatch.setattr(kt, "latest_release", lambda timeout=20.0: release)
    monkeypatch.setattr(kt, "pinned_release", lambda: release, raising=False)
    monkeypatch.setattr(kt.urllib.request, "urlopen",
                        lambda request, timeout=None: io.BytesIO(payload))


def _nothing_installed() -> bool:
    return not kt.stored_directory().exists()


def test_the_release_is_pinned_and_needs_no_network(monkeypatch):
    def no_network(*args, **kwargs):
        raise AssertionError("the pinned release was asked of the network")

    monkeypatch.setattr(kt.urllib.request, "urlopen", no_network)
    release = kt.pinned_release()
    assert release["ok"] is True
    assert release["version"] == "3.18.0"
    assert release["url"] == ("https://github.com/epam/ketcher/releases/download/"
                              "v3.18.0/ketcher-standalone-3.18.0.zip")
    assert release["sha256"] == (
        "484e7f10a0e74808ae5f43f0e6448db6f4d884bdd2189aa271108474167a3e7e")
    assert release["size"] == 35_057_618
    # What the dashboard asks before offering the download is the pin too.
    assert kt.latest_release()["sha256"] == release["sha256"]


def test_a_download_that_differs_from_the_pin_is_refused(home, monkeypatch):
    _serve(monkeypatch, _archive(), sha256="0" * 64)
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert "SHA-256" in outcome["status"]
    assert _nothing_installed()


def test_a_download_larger_than_the_pin_is_refused(home, monkeypatch):
    payload = _archive()
    _serve(monkeypatch, payload, size=len(payload) - 1)
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert _nothing_installed()


def test_a_member_that_climbs_out_is_refused(home, monkeypatch):
    _serve(monkeypatch, _archive(("static/js/../../../../outside.chunk.js", "// x")))
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert "outside.chunk.js" in outcome["status"]
    assert not (home / ".delfin" / "outside.chunk.js").exists()
    assert _nothing_installed()


def test_an_absolute_member_is_refused(home, tmp_path, monkeypatch):
    target = tmp_path / "absolute.chunk.js"
    _serve(monkeypatch, _archive((str(target), "// x")))
    assert zipfile.ZipFile(io.BytesIO(_archive((str(target), "")))).namelist()[-1] == str(target)
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert not target.exists()
    assert _nothing_installed()


def test_a_backslash_member_is_refused(home, monkeypatch):
    _serve(monkeypatch, _archive(("static\\js\\..\\..\\x.chunk.js", "// x")))
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert _nothing_installed()


def test_a_symlink_member_is_refused(home, monkeypatch):
    link = zipfile.ZipInfo("static/js/link.chunk.js")
    link.external_attr = (stat.S_IFLNK | 0o777) << 16
    _serve(monkeypatch, _archive((link, "index.html")))
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"] is False
    assert "link.chunk.js" in outcome["status"]
    assert _nothing_installed()


def test_a_clean_archive_at_the_pin_installs(home, monkeypatch):
    _serve(monkeypatch, _archive(("static/js/157.chunk.js", "// chunk")))
    outcome = kt.install(folder=kt.stored_directory())
    assert outcome["ok"], outcome["status"]
    assert kt.stored_version() == "3.18.0"
    assert (kt.stored_directory() / "static" / "js" / "157.chunk.js").is_file()
