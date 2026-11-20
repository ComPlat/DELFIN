"""The reference constraints file pins the construction's closure, not an interpreter.

``env/reference-environment.lock`` is installed with, as the README says
(``pip install -e ".[agent,docs,dev]" -c env/reference-environment.lock``),
and CI's reference-environment job does the same. pip applies a constraint
to any package it installs anyway, so every line in that file fixes a
version for the person who follows the README.

Until 2026-11-20 the file was ``pip freeze --all`` of the build
interpreter: 143 pins, of which 110 are not in the construction's
dependency closure at all. Measured with pip-audit on 2026-11-20: 59 known
advisories across 14 of the pinned distributions -- urllib3, cryptography,
requests, certifi, anyio, lxml, pygments, setuptools, wheel, pip among them
-- and not one of them a package that computes a coordinate. The
reproducibility of the chemistry never needed them pinned.

Input: the two files in ``env/``. Output: which names each one holds.

Semantics pinned here:

* the constraints file holds only names in the construction's closure, and
  the record holds the whole interpreter. A re-freeze of an interpreter into
  the constraints file brings back the defect, and this test is what
  notices;
* the versions the published structures depend on are the published ones.
  Those are the file's whole reason to exist, and a bump to one of them is a
  statement about the data, not a dependency update;
* the record is never named as something to install with.

Not pinned here: that the pins carry no advisory. That is a question about
the world, not about the repository, and it is asked where the answer can
change without a commit -- .github/workflows/audit.yml, which now resolves
this file and runs pip-audit over it.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ENV = Path(__file__).resolve().parents[1] / "env"
LOCK = ENV / "reference-environment.lock"
RECORD = ENV / "reference-environment-interpreter.txt"

#: What the published MANTA structures were built with. Geometry comes out
#: of rdkit, Open Babel and numpy; stk's optimisers move the atoms.
PUBLISHED = {
    "numpy": "2.2.2",
    "rdkit": "2025.9.3",
    "openbabel-wheel": "3.1.1.22",
    "scipy": "1.16.2",
    "networkx": "3.1",
    "stk": "2026.1.4.0",
}

#: The closure of the construction, resolved on 2026-11-20 with
#: ``pip install --dry-run --report`` for the six packages above on CPython
#: 3.12. Three members are deliberately absent from the pins -- pymongo,
#: dnspython (atomlite's database driver and its SRV resolver) and fonttools
#: (matplotlib's font parsing) -- each with a known advisory and none able to
#: move an atom. They are still installed, at a patched version.
CLOSURE = {
    "atomlite", "contourpy", "cycler", "dill", "dnspython", "fonttools",
    "kiwisolver", "matplotlib", "mchammer", "multiprocess", "networkx",
    "numpy", "openbabel-wheel", "packaging", "pandas", "pathos", "pillow",
    "polars", "polars-runtime-32", "pox", "ppft", "pymongo", "pyparsing",
    "python-dateutil", "pytz", "rdkit", "scipy", "seaborn", "six",
    "spindry", "stk", "tzdata", "vabene",
}

_PIN = re.compile(r"^([A-Za-z0-9._-]+)==([^\s;]+)")


def _pins(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = _PIN.match(line)
        assert m, f"{path.name}: not a pin: {line!r}"
        out[m.group(1).lower().replace("_", "-")] = m.group(2)
    return out


def test_both_files_are_there():
    assert LOCK.is_file(), "the constraints file the README names"
    assert RECORD.is_file(), "the record of the build interpreter"


def test_every_pin_is_in_the_construction_closure():
    """The one that catches a re-freeze. A pin outside the closure fixes a
    version for a package the construction does not use -- which is how ten
    network and crypto libraries came to be held at vulnerable versions."""
    stray = sorted(set(_pins(LOCK)) - CLOSURE)
    assert not stray, (
        "pins outside the construction's dependency closure: "
        f"{stray}. A constraint is applied to any package pip installs, so "
        "each of these fixes a version for a package the reproducibility of "
        "the chemistry does not depend on. Pin it only if it can change a "
        "structure; otherwise leave it out and let it be patched.")


@pytest.mark.parametrize("name,version", sorted(PUBLISHED.items()))
def test_the_published_versions_are_the_pinned_ones(name, version):
    """A bump here is a statement about the published data, not a dependency
    update: the structures were built with these."""
    assert _pins(LOCK).get(name) == version


def test_the_record_still_holds_the_whole_interpreter():
    """It is provenance. Narrowing it would delete the record of the
    environment the published numbers come from."""
    pins = _pins(RECORD)
    assert len(pins) >= 143, f"the record holds {len(pins)} pins"
    # The packages whose pins were the defect are still recorded, because
    # what was installed that day is a fact.
    for name in ("urllib3", "cryptography", "requests", "setuptools"):
        assert name in pins


def test_the_record_says_not_to_install_with_it():
    head = RECORD.read_text()[:2000].upper()
    assert "DO NOT INSTALL" in head
    assert "reference-environment.lock" in RECORD.read_text()[:2000]


def test_the_readme_installs_with_the_constraints_file_and_names_the_record():
    """The section a reader follows has to point at the file whose pins are
    audited, and say where the record went."""
    readme = (ENV.parent / "README.md").read_text()
    assert "-c env/reference-environment.lock" in readme
    assert RECORD.name in readme, "the record must stay findable from the README"
    install_lines = [ln for ln in readme.splitlines()
                     if ln.startswith("pip install")
                     and "reference-environment" in ln]
    assert install_lines, "the README shows how to install against it"
    for line in install_lines:
        assert RECORD.name not in line, (
            f"the README tells a reader to install with the record: {line}")
