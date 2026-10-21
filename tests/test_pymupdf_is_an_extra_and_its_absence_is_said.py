"""PyMuPDF (AGPL-3.0) is opt-in, and every place that uses it says why it is missing.

Input: ``pyproject.toml``, ``scripts/check_dependency_licences.py`` with its
allow-list, and the five import sites of PyMuPDF with the module made
unimportable (``sys.modules[name] = None`` makes ``import name`` raise
ImportError). Output: assertions.

PyMuPDF was a core dependency since 0e05f90a; its import was already
optional at every site, so moving it to the ``[pdf]`` extra changes what a
plain install brings in, not what the code requires. ``[dev]`` lists it too,
so the CI jobs that install ``[agent,docs,dev]`` keep running the PDF tests
instead of skipping them.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest
import tomllib

REPO = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "check_dependency_licences", REPO / "scripts" / "check_dependency_licences.py")
CDL = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = CDL
_spec.loader.exec_module(CDL)


def _project() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text())["project"]


def _names(requirements: list[str]) -> set[str]:
    import re
    return {re.split(r"[<>=!~\[; ]", r, maxsplit=1)[0].lower() for r in requirements}


# ---------------------------------------------------------------------------
# Packaging
# ---------------------------------------------------------------------------

def test_pymupdf_is_not_a_core_dependency():
    assert "pymupdf" not in _names(_project()["dependencies"])


def test_the_pdf_extra_installs_it_and_dev_carries_the_same_requirement():
    extras = _project()["optional-dependencies"]
    pdf = [r for r in extras["pdf"] if r.lower().startswith("pymupdf")]
    assert pdf, "the [pdf] extra must install pymupdf"
    assert pdf[0] in extras["dev"], "[dev] must carry the [pdf] requirement verbatim"


def test_no_extra_that_a_plain_agent_install_selects_pulls_it_in():
    extras = _project()["optional-dependencies"]
    for name in ("agent", "office", "docs", "test"):
        assert "pymupdf" not in _names(extras.get(name, [])), name


# ---------------------------------------------------------------------------
# The licence check, driven both ways
# ---------------------------------------------------------------------------

def _allowlist() -> dict:
    return CDL.load_allowlist(CDL.DEFAULT_ALLOWLIST)


@pytest.mark.parametrize("licence", [
    "MIT", "MIT License", "BSD-3-Clause", "Apache-2.0", "ISC License (ISCL)",
    "Python Software Foundation License", "MPL-2.0", "LGPL-2.1-or-later",
    "GNU Library or Lesser General Public License (LGPL)", "MIT-CMU",
    "BSD-2-Clause AND Apache-2.0 WITH LLVM-exception",
])
def test_permissive_and_lgpl_licences_pass(licence):
    assert CDL.classify(licence) == "pass"


@pytest.mark.parametrize("licence", [
    "Dual Licensed - GNU AFFERO GPL 3.0 or Artifex Commercial License",
    "AGPL-3.0-only", "MIT AND AGPL-3.0", "GNU General Public License (GPL)",
    "GPL-2.0-only", "GPLv3", "SSPL-1.0", "Proprietary", "UNKNOWN", "",
    "CC-BY-NC-4.0",
])
def test_copyleft_proprietary_and_unknown_need_review(licence):
    assert CDL.classify(licence) == "review"


def test_pymupdf_in_a_core_environment_fails_the_check():
    packages = [{"Name": "PyMuPDF", "Version": "1.28.2", "License":
                 "Dual Licensed - GNU AFFERO GPL 3.0 or Artifex Commercial License"}]
    failures, _notes = CDL.check(packages, _allowlist())
    assert failures and "PyMuPDF" in failures[0]


def test_open_babel_passes_only_with_the_reviewed_string():
    ok = [{"Name": "openbabel-wheel", "Version": "3.1.1.23",
           "License": "GNU General Public License (GPL)"}]
    assert CDL.check(ok, _allowlist())[0] == []
    changed = [dict(ok[0], License="GNU General Public License v3 (GPLv3)")]
    failures = CDL.check(changed, _allowlist())[0]
    assert failures and "review it again" in failures[0]


def test_an_unknown_licence_needs_an_entry():
    failures = CDL.check([{"Name": "newpkg", "Version": "1", "License": "UNKNOWN"}],
                         _allowlist())[0]
    assert failures and "newpkg" in failures[0]


def test_every_allow_list_entry_is_reviewed_and_explained():
    for name, entry in _allowlist().items():
        assert len(entry["reason"].split()) >= 15, name
        assert entry["spdx"], name


def test_the_command_line_exits_one_on_a_finding(tmp_path, capsys):
    clean = tmp_path / "clean.json"
    clean.write_text(json.dumps([{"Name": "numpy", "Version": "1", "License": "BSD License"}]))
    assert CDL.main([str(clean)]) == 0
    dirty = tmp_path / "dirty.json"
    dirty.write_text(json.dumps([{"Name": "x", "Version": "1", "License": "AGPL-3.0"}]))
    assert CDL.main([str(dirty)]) == 1


def test_the_workflow_runs_the_check_on_a_plain_install():
    text = (REPO / ".github" / "workflows" / "packaging.yml").read_text()
    job = text[text.index("  dependency-licences:"):]
    assert "pip install .\n" in job, "the check must see the core set, no extras"
    # pip-licenses at an exact version, through the CI constraints file.
    assert "pip install -c .github/ci/constraints-ci.txt pip-licenses" in job
    pins = (REPO / ".github" / "ci" / "constraints-ci.txt").read_text()
    assert re.search(r"^pip-licenses==\d", pins, re.M)
    assert "scripts/check_dependency_licences.py" in job


# ---------------------------------------------------------------------------
# Every import site, with PyMuPDF absent
# ---------------------------------------------------------------------------

@pytest.fixture
def no_pymupdf(monkeypatch):
    monkeypatch.setitem(sys.modules, "pymupdf", None)
    monkeypatch.setitem(sys.modules, "fitz", None)


def test_the_pdf_viewer_names_the_extra(no_pymupdf):
    from delfin.dashboard import pdf_view
    with pytest.raises(pdf_view.PdfError) as err:
        pdf_view._fitz()
    message = str(err.value)
    assert "PyMuPDF is not installed" in message
    assert "delfin-complat[pdf]" in message and "AGPL" in message


def test_the_ocr_check_names_the_missing_renderer(no_pymupdf):
    from delfin.agent import office
    status = office.ocr_availability()
    assert status["renderer"] == ""
    assert not status["available"]
    assert any("delfin-complat[pdf]" in d for d in status["detail"])


def test_rendering_a_page_raises_a_reason_not_an_import_error(no_pymupdf, tmp_path):
    from delfin.agent import office
    with pytest.raises(office.OfficeError) as err:
        office._render_page(tmp_path / "a.pdf", 1)
    assert "delfin-complat[pdf]" in str(err.value)


def test_form_labels_degrade_to_none(no_pymupdf, tmp_path):
    from delfin.agent import office
    assert office._field_labels(tmp_path / "a.pdf") == {}


def test_the_chat_preview_degrades_to_no_preview(no_pymupdf, tmp_path):
    from delfin.dashboard import tab_agent
    assert tab_agent._render_pdf_preview(tmp_path / "a.pdf") is None
