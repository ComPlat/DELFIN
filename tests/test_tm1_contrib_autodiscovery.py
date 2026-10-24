"""Every module under delfin/workflows/contrib is discovered without a hand edit."""
import pkgutil

from delfin.workflows import contrib, registry


def test_every_contrib_module_is_imported_by_discovery(monkeypatch):
    import sys
    names = [m.name for m in pkgutil.iter_modules(contrib.__path__) if not m.name.startswith("_")]
    assert names
    for n in names:
        monkeypatch.delitem(sys.modules, f"delfin.workflows.contrib.{n}", raising=False)
    monkeypatch.setattr(registry, "_DISCOVERED", False)
    registry._discover_builtin_workflows()
    missing = [n for n in names if f"delfin.workflows.contrib.{n}" not in sys.modules]
    assert not missing, missing
