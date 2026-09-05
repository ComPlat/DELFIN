"""No test may leave a fake dashboard module behind for the others.

pytest imports every test module into ONE process, so a stub written into
sys.modules at import time answers for the rest of the run. A bare
ModuleType has no __file__, and a later module doing
"from .helpers import disable_spellcheck" then fails with
"cannot import name ... (unknown location)" — an error that points at a
file which is perfectly fine, in a module nobody touched.

It reached CI exactly that way: a stub had sat harmlessly in
test_backend_local.py until browser tests started importing the
dashboard, and that file sorts before them.

This module is named so it collects AFTER the stubbing ones.
"""

from __future__ import annotations

import sys


def test_the_real_helpers_module_is_the_one_in_sys_modules():
    import delfin.dashboard.helpers as helpers

    assert getattr(helpers, "__file__", None), (
        "delfin.dashboard.helpers has no __file__ — an earlier test replaced "
        "it with a stub, and every later import of it will fail"
    )
    assert helpers.__file__.endswith("delfin/dashboard/helpers.py")


def test_the_names_the_browser_imports_are_all_there():
    """The exact import that failed in CI."""
    from delfin.dashboard.helpers import (  # noqa: F401
        disable_spellcheck,
        save_neb_trajectory_csv,
        save_neb_trajectory_plot_png,
    )


def test_a_stand_in_package_still_points_at_the_real_directory():
    """A PACKAGE stand-in is tolerable where a module stand-in is not: it
    carries __path__, so every submodule under it still resolves to the
    real file. What breaks the run is a leaf module with no __file__,
    because there is nothing left to resolve through."""
    import delfin.dashboard as dashboard

    paths = [str(p) for p in getattr(dashboard, "__path__", [])]
    assert paths, "delfin.dashboard has no __path__ at all"
    assert any(p.endswith("delfin/dashboard") for p in paths), paths


def test_no_dashboard_module_in_sys_modules_lacks_a_file():
    stubs = [
        name for name, module in list(sys.modules.items())
        if name.startswith("delfin.dashboard")
        and module is not None
        and getattr(module, "__file__", None) is None
        and not hasattr(module, "__path__")
    ]
    assert not stubs, f"stubbed dashboard modules left behind: {stubs}"
