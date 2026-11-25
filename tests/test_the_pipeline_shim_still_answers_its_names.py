"""`delfin.pipeline` still answers every name it answered before.

`delfin/pipeline.py` is the backward-compat shim over
`delfin/workflows/pipeline.py` (``from delfin.workflows.pipeline import *``),
and it is the documented way user code still imports from the old flat path
(``tests/test_shim_integrity.py``, ``tests/test_pre_refactor_api_preserved.py``).

`find_gibbs_energy` reached that public surface only through the star import:
`workflows/pipeline.py` imported it from `delfin.energies` and never called it
(zero name references, AST-verified). Removing that import as "dead" therefore
did not only delete a dead name -- it deleted a live, documented import for
every caller of `from delfin.pipeline import find_gibbs_energy`, and no test
noticed. These tests bind both directions: the shim keeps the name, and the
module that no longer shadows it still has no call site.

The control cases pass on both trees, so a failure here is the shim surface,
not a fixture.
"""
import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SHIM = ROOT / "delfin" / "pipeline.py"
WORKFLOWS = ROOT / "delfin" / "workflows" / "pipeline.py"


def test_the_shim_still_answers_find_gibbs_energy():
    """The name user code imports from delfin.pipeline is still there."""
    from delfin import pipeline as shim
    from delfin import energies

    assert hasattr(shim, "find_gibbs_energy"), (
        "delfin.pipeline no longer answers `find_gibbs_energy`. It reached the "
        "shim through the star import of delfin.workflows.pipeline; removing "
        "that import breaks `from delfin.pipeline import find_gibbs_energy` "
        "for user code, which is not the same change as deleting a dead name. "
        "Keep an explicit re-export in the shim, or waive the drop in "
        "tests/test_pre_refactor_api_preserved.py with a reason."
    )
    assert shim.find_gibbs_energy is energies.find_gibbs_energy


def test_the_workflow_module_has_no_call_of_that_name():
    """The other half: it is still unused there, so nothing silently moved."""
    tree = ast.parse(WORKFLOWS.read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree)
             if isinstance(node, ast.Name) and node.id == "find_gibbs_energy"]
    assert not calls, (
        "workflows/pipeline.py now calls find_gibbs_energy at line(s) "
        f"{[n.lineno for n in calls]}; the shim test above assumed it never "
        "did (it imports the name only for the shim to re-export)."
    )


@pytest.mark.parametrize("name", ["PipelineContext", "run_esd_phase",
                                  "collect_gibbs_energies"])
def test_the_other_shim_names_are_unaffected(name):
    """Control: names with a real call site keep working (passes on both)."""
    from delfin import pipeline as shim

    assert hasattr(shim, name)
