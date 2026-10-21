#!/usr/bin/env python3
"""Setup for generalist_large_repo: one defect among ~150 files.

``inventory/`` -- a generated package: core (registry, importer),
validation (sku, names, quantities), 60 plugins, 30 report modules, 20
docs pages and 30 test files. The defect: validation/sku.py compiles
the SKU pattern with re.IGNORECASE and uses ``match`` (a prefix match),
so "abc-1234" and "ABC-12345" are accepted. Both entry points --
``Registry.add`` and ``importer.import_rows`` -- go through
``validate_sku``; a fix at one entry point leaves the other open.
plugins/plugin_017.py has a look-alike ``looks_like_sku`` that nothing
in core uses (a distractor).
Preconditions: >= 150 tracked files; the suite is green; the defect is
reproducible through both entry points.
argv[1] is the workspace.
"""
from __future__ import annotations

import subprocess
import sys

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

SKU = '''"""Stock-keeping unit validation.

A SKU is three upper-case letters, a dash and exactly four digits:
"ABC-1234".
"""

import re

_SKU = re.compile(r"[A-Z]{3}-[0-9]{4}", re.IGNORECASE)


class InvalidSku(ValueError):
    """Raised for a malformed SKU."""


def validate_sku(sku: str) -> str:
    """Return ``sku`` unchanged when it is well formed, else raise."""
    if not isinstance(sku, str) or not _SKU.match(sku):
        raise InvalidSku(f"invalid SKU: {sku!r}")
    return sku
'''

NAMES = '''"""Item name validation."""


def validate_name(name: str) -> str:
    """A stripped, non-empty name of at most 80 characters."""
    name = (name or "").strip()
    if not name or len(name) > 80:
        raise ValueError(f"invalid name: {name!r}")
    return name
'''

QUANTITIES = '''"""Quantity validation."""


def validate_quantity(qty) -> int:
    """A non-negative integer quantity."""
    if not isinstance(qty, int) or isinstance(qty, bool) or qty < 0:
        raise ValueError(f"invalid quantity: {qty!r}")
    return qty
'''

REGISTRY = '''"""The item registry."""

from inventory.validation.names import validate_name
from inventory.validation.quantities import validate_quantity
from inventory.validation.sku import validate_sku


class Registry:
    """Items by SKU."""

    def __init__(self):
        self._items = {}

    def add(self, sku, name, qty=0):
        """Register an item; returns the stored record."""
        sku = validate_sku(sku)
        record = {"sku": sku, "name": validate_name(name),
                  "qty": validate_quantity(qty)}
        self._items[sku] = record
        return record

    def get(self, sku):
        """The record for ``sku`` or None."""
        return self._items.get(sku)

    def __len__(self):
        return len(self._items)
'''

IMPORTER = '''"""Bulk import of CSV rows."""

import csv
import io

from inventory.validation.sku import InvalidSku, validate_sku


def import_rows(registry, text):
    """Add every ``sku,name,qty`` row of ``text`` to ``registry``.

    Returns (added, rejected) counts; a malformed SKU rejects its row.
    """
    added = rejected = 0
    for row in csv.reader(io.StringIO(text)):
        if not row:
            continue
        sku, name, qty = (row + ["", "", "0"])[:3]
        try:
            validate_sku(sku.strip())
        except InvalidSku:
            rejected += 1
            continue
        registry.add(sku.strip(), name, int(qty or 0))
        added += 1
    return added, rejected
'''


def _plugin(i: int) -> str:
    if i == 17:
        return ('"""Plugin 17: label printing."""\n\nimport re\n\n\n'
                'def looks_like_sku(text):\n'
                '    """Loose check used only to pick label templates."""\n'
                '    return bool(re.match(r"[A-Za-z]{3}-\\d+", text or ""))\n\n\n'
                'def label(sku, name):\n'
                '    """One printable label line."""\n'
                '    kind = "sku" if looks_like_sku(sku) else "free"\n'
                '    return f"[{kind}] {sku} {name}"\n')
    return (f'"""Plugin {i}: price rule {i}."""\n\n\n'
            f'def adjust_{i:03d}(cents):\n'
            f'    """Apply rule {i}: add {i % 7} cents per started euro."""\n'
            f'    euros = -(-cents // 100)\n'
            f'    return cents + euros * {i % 7}\n')


def _report(i: int) -> str:
    return (f'"""Report {i}: stock summary variant {i}."""\n\n\n'
            f'def summary_{i:03d}(records):\n'
            f'    """Total quantity of records whose SKU ends in {i % 10}."""\n'
            f'    return sum(r["qty"] for r in records\n'
            f'               if r["sku"].endswith("{i % 10}"))\n')


def _test(i: int) -> str:
    return (f'from inventory.plugins.plugin_{i:03d} import adjust_{i:03d}\n'
            f'from inventory.reports.report_{i:03d} import summary_{i:03d}\n\n\n'
            f'def test_plugin_{i:03d}():\n'
            f'    assert adjust_{i:03d}(250) == 250 + 3 * {i % 7}\n\n\n'
            f'def test_report_{i:03d}():\n'
            f'    recs = [{{"sku": "ABC-000{i % 10}", "qty": 2}},\n'
            f'            {{"sku": "ABC-0011", "qty": 5}}]\n'
            f'    assert summary_{i:03d}(recs) == '
            f'{2 + (5 if i % 10 == 1 else 0)}\n')


CORE_TEST = '''import pytest

from inventory.core.importer import import_rows
from inventory.core.registry import Registry
from inventory.validation.sku import InvalidSku


def test_add_and_get():
    r = Registry()
    r.add("ABC-1234", "Bolt", 3)
    assert r.get("ABC-1234")["qty"] == 3


def test_a_sku_without_dash_is_refused():
    with pytest.raises(InvalidSku):
        Registry().add("ABC1234", "Bolt")


def test_import_counts_rows():
    r = Registry()
    assert import_rows(r, "ABC-1234,Bolt,3\\nXY-1,Nut,1\\n") == (1, 1)
'''


def build(workspace: str):
    repo = new_repo(workspace, "inventory")
    files = {
        ".gitignore": GITIGNORE_PY,
        "README.md": ("# inventory\n\nStock keeping. Entry points: "
                      "`Registry.add` and `import_rows`. See docs/.\n"),
        "pyproject.toml": ("[project]\nname = \"inventory\"\n"
                           "version = \"3.4.0\"\n\n"
                           "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                           "testpaths = [\"tests\"]\n"),
        "inventory/__init__.py": "\"\"\"Stock keeping.\"\"\"\n",
        "inventory/core/__init__.py": "",
        "inventory/core/registry.py": REGISTRY,
        "inventory/core/importer.py": IMPORTER,
        "inventory/validation/__init__.py": "",
        "inventory/validation/sku.py": SKU,
        "inventory/validation/names.py": NAMES,
        "inventory/validation/quantities.py": QUANTITIES,
        "inventory/plugins/__init__.py": "",
        "inventory/reports/__init__.py": "",
        "tests/test_core.py": CORE_TEST,
    }
    for i in range(60):
        files[f"inventory/plugins/plugin_{i:03d}.py"] = _plugin(i)
    for i in range(30):
        files[f"inventory/reports/report_{i:03d}.py"] = _report(i)
        if i != 17:
            files[f"tests/test_module_{i:03d}.py"] = _test(i)
    for i in range(20):
        files[f"docs/page_{i:02d}.md"] = (
            f"# Page {i}\n\nOperating notes, part {i}. SKUs look like "
            f"`ABC-{1000 + i}`.\n")
    write(repo, files)
    mark_base(repo, commit(repo, "inventory 3.4.0"))
    tracked = len(subprocess.run(["git", "ls-files"], cwd=str(repo),
                                 capture_output=True, text=True
                                 ).stdout.split())
    expect(">= 150 tracked files", tracked >= 150, str(tracked))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    probe = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, '.')\n"
         "from inventory.core.registry import Registry\n"
         "from inventory.core.importer import import_rows\n"
         "Registry().add('abc-1234', 'x')\n"
         "assert import_rows(Registry(), 'ABC-12345,x,1') == (1, 0)\n"],
        cwd=str(repo), capture_output=True, text=True, timeout=60,
        env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin"})
    expect("the defect is reproducible", probe.returncode == 0,
           probe.stderr)
    return repo


if __name__ == "__main__":
    main(build)
