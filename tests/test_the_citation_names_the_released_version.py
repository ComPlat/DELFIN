"""CITATION.cff names the version the package declares, and the README cites it.

Input: ``CITATION.cff``, ``pyproject.toml`` (``project.version``),
``delfin/__init__.py`` (``__version__``) and the "Please cite" section of
``README.md``. Output: one assertion per pair that must agree.

Semantics:

* ``version`` in CITATION.cff equals ``project.version`` and
  ``delfin.__version__`` (a leading ``v`` is ignored, as the release tags
  carry one);
* the top-level ``doi`` is one of the ``identifiers``, so the field a
  reader sees first is also listed with its description;
* the README's BibTeX ``@software`` entry carries the same version, and
  every DOI in CITATION.cff (the concept DOI, the version DOI and the
  preferred citation's DOI) appears in the README's citation section.

Observed before this test: CITATION.cff said ``v1.3.2`` (no such tag) while
the package and the latest tag were 1.3.3, and its DOI was the v1.2.0
version record; the README's BibTeX had the key ``hartmann2025delfin`` with
``year = {2026}``.

Rejected: checking the DOIs against Zenodo -- a test must pass offline and
on every machine. Rejected: requiring the version DOI to belong to
``version`` -- Zenodo mints a version DOI only after the tag exists, so a
tagged CITATION.cff cannot carry its own; the identifier's description
names the release it belongs to instead.
Cost: four small file reads.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _cff() -> dict:
    return yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))


def _bare(version: str) -> str:
    return str(version).strip().removeprefix("v")


def _cite_section() -> str:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    start = text.index("## Please cite")
    end = text.index("\n## ", start + 1)
    return text[start:end]


def test_the_citation_version_is_the_package_version():
    # The [project] table's version line, read without tomllib (absent on
    # Python 3.10, which the package still supports).
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    table = re.search(r"^\[project\]\n(.*?)(?=^\[)", pyproject,
                      re.MULTILINE | re.DOTALL)
    assert table, "pyproject.toml has no [project] table"
    pv = re.search(r'^version\s*=\s*"([^"]+)"', table.group(1), re.MULTILINE)
    assert pv, "pyproject.toml [project] has no version"
    init = (ROOT / "delfin" / "__init__.py").read_text(encoding="utf-8")
    m = re.search(r'^__version__\s*=\s*"([^"]+)"', init, re.MULTILINE)
    assert m, "delfin/__init__.py defines no __version__"
    cff_version = _bare(_cff()["version"])
    assert cff_version == _bare(pv.group(1)), (
        f"CITATION.cff version {cff_version} != pyproject {pv.group(1)}")
    assert cff_version == _bare(m.group(1)), (
        f"CITATION.cff version {cff_version} != __version__ {m.group(1)}")


def test_the_top_level_doi_is_a_listed_identifier():
    cff = _cff()
    dois = [i["value"] for i in cff.get("identifiers", [])
            if i.get("type") == "doi"]
    assert cff["doi"] in dois


def test_the_readme_cites_the_same_version_and_dois():
    cff = _cff()
    section = _cite_section()
    software = re.search(r"@software\{[^,]+,(.*?)\n\}", section, re.DOTALL)
    assert software, "README has no @software BibTeX entry"
    v = re.search(r"version\s*=\s*\{([^}]+)\}", software.group(1))
    assert v and _bare(v.group(1)) == _bare(cff["version"]), (
        "README @software version differs from CITATION.cff")
    dois = {cff["doi"], cff["preferred-citation"]["doi"]}
    dois |= {i["value"] for i in cff.get("identifiers", [])
             if i.get("type") == "doi"}
    missing = sorted(d for d in dois if d not in section)
    assert not missing, f"README citation section lacks {missing}"


def test_each_bibtex_key_names_its_year():
    """A key ``name<year>...`` whose year is not the entry's year is how
    ``hartmann2025delfin`` came to stand on a 2026 record."""
    section = _cite_section()
    for m in re.finditer(r"@\w+\{([A-Za-z]+)(\d{4})\w*,(.*?)\n\}", section,
                         re.DOTALL):
        year = re.search(r"year\s*=\s*\{(\d{4})\}", m.group(3))
        assert year and year.group(1) == m.group(2), (
            f"BibTeX key {m.group(1)}{m.group(2)}... has year "
            f"{year.group(1) if year else None}")
