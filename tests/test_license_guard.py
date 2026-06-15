# license-guard: allow  (this test embeds synthetic forbidden snippets)
"""CI enforcement of the CCDC/CSD license guard.

Two guarantees:
  1. The shipped, git-tracked tree contains NO CCDC/CSD material
     (test_tracked_tree_is_ccdc_free).
  2. The guard actually detects violations and cannot silently rot into a
     no-op (test_guard_detects_*).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import license_guard as lg  # noqa: E402


def test_tracked_tree_is_ccdc_free():
    """No git-tracked Python file may import CCDC or reference CCDC data."""
    violations = []
    for rel in lg.tracked_py_files():
        if rel in lg.SELF_EXEMPT:
            continue
        text = (REPO / rel).read_text(encoding="utf-8", errors="replace")
        for lineno, msg in lg.scan_py_source(rel, text):
            violations.append(f"{rel}:{lineno}: {msg}")
    assert not violations, "CCDC/CSD material in shipped tree:\n" + "\n".join(violations)


@pytest.mark.parametrize(
    "snippet, why",
    [
        ("import ccdc\n", "ccdc runtime API"),
        ("from ccdc.io import EntryReader\n", "ccdc submodule"),
        ("from delfin.fffree.mogul_bounds import bounds\n", "mogul-derived module"),
        ("import delfin.grip_io\n", "grip-derived module"),
        ("PATH = 'grip_lib_v5_TM.npz'\n", "grip fragment library path"),
        ("LIB = 'agent_workspace/ccdc_fragment_index_v3.json'\n", "ccdc fragment index"),
        ("data = np.load('foo.npz')\n", ".npz fragment load"),
    ],
)
def test_guard_detects_violations(snippet, why):
    """The guard must flag each canonical violation form."""
    found = lg.scan_py_source("delfin/_synthetic_probe.py", snippet)
    assert found, f"guard FAILED to detect {why}: {snippet!r}"


@pytest.mark.parametrize(
    "snippet",
    [
        "from delfin.fffree.cod_ideals import ideal_bond\n",   # COD = open
        "from delfin.fffree.polyhedra import md_distance\n",   # covalent radii = open
        "RADII = {'C': 0.76, 'N': 0.71}  # Cordero covalent radii\n",
        "IDX = 'cod_fragment_index_v3.sqlite'  # COD is open\n",
    ],
)
def test_guard_allows_open_sources(snippet):
    """Open, shippable sources must NOT be flagged."""
    found = lg.scan_py_source("delfin/_synthetic_probe.py", snippet)
    assert not found, f"guard wrongly flagged an OPEN source: {snippet!r} -> {found}"


def test_allow_comment_escape_hatch():
    snippet = "PATH = 'grip_lib.npz'  # license-guard: allow\n"
    assert not lg.scan_py_source("delfin/_x.py", snippet)


@pytest.mark.parametrize(
    "path",
    ["agent_workspace/grip_lib_v5_TM.npz", "data/ALPHA.mol2",
     "reports/ccdc_fragment_index_v3.json", "scripts/ccdc_validation.py"],
)
def test_path_guard_blocks_data_files(path):
    assert lg.check_path(path) is not None, f"path guard missed {path}"


@pytest.mark.parametrize("path", ["delfin/fffree/cod_ideals.py", "docs/intro.md"])
def test_path_guard_allows_clean_paths(path):
    assert lg.check_path(path) is None, f"path guard wrongly blocked {path}"
