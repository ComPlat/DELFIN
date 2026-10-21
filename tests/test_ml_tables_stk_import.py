"""stk is optional only when it is not installed; any other import failure must surface.

A swallowed ImportError (e.g. an older libstdc++ from another program's LD_LIBRARY_PATH entry,
loaded before stk's sqlite3 extension) used to set STK_AVAILABLE=False for the whole process, so
the construction depended on module import order."""
import builtins
import importlib
import sys

import pytest

import delfin.manta.ml_tables as ml_tables


def _reload_with_import(monkeypatch, exc):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "stk":
            raise exc
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "stk", raising=False)
    monkeypatch.setattr(builtins, "__import__", fake_import)
    return importlib.reload(ml_tables)


@pytest.fixture(autouse=True)
def _restore_module(monkeypatch):
    yield
    monkeypatch.undo()                   # the real import first, then the real module
    importlib.reload(ml_tables)


def test_not_installed_is_optional(monkeypatch):
    mod = _reload_with_import(monkeypatch, ModuleNotFoundError("No module named 'stk'", name="stk"))
    assert mod.STK_AVAILABLE is False
    assert mod.stk is None


def test_missing_dependency_of_stk_is_raised(monkeypatch):
    with pytest.raises(ModuleNotFoundError):
        _reload_with_import(monkeypatch, ModuleNotFoundError("No module named 'atomlite'", name="atomlite"))


def test_broken_shared_library_is_raised(monkeypatch):
    with pytest.raises(ImportError, match="CXXABI"):
        _reload_with_import(monkeypatch, ImportError("libstdc++.so.6: version `CXXABI_1.3.15' not found"))
