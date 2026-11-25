"""Re-run the x1 arms' offender-count comparison, this turn, from the two real
guard modules on their branches. Prints what each reports for one file that
holds TWO undefined escapes.

Run from the workspace root: PYTHONPATH=. python3 scratch/x1/rerun/two_escape.py
"""

import pathlib
import sys
import tempfile
import types

# a2's guard imports pytest for one `pytest.raises` block; stub it so this runs
# without pytest on the path.
stub = types.ModuleType("pytest")


class _Raises:
    def __init__(self, exc):
        self.exc = exc

    def __enter__(self):
        return self

    def __exit__(self, t, v, tb):
        if t is None:
            raise AssertionError("did not raise %r" % (self.exc,))
        return issubclass(t, self.exc)


stub.raises = _Raises
sys.modules["pytest"] = stub

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
# a2's guard imports the repo's own tests/_package_ast helper, as it does when
# pytest collects it from tests/.
sys.path.insert(0, str(HERE.parents[2] / "tests"))


def load(name, path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


a1 = load("a1_guard", HERE / "a1_guard.py")
a2 = load("a2_guard", HERE / "a2_guard.py")

# Two undefined escapes, on the two docstring opening lines of the fixture.
with tempfile.TemporaryDirectory() as td:
    probe = pathlib.Path(td) / "two_bad.py"
    probe.write_text(
        'def f():\n'
        '    """a\\ b in one"""\n'
        "    return 1\n"
        "def g():\n"
        '    """c\\ d in two"""\n'
        "    return 2\n"
    )
    source = probe.read_text()
    print("fixture lines:")
    for n, line in enumerate(source.splitlines(), 1):
        print(f"  {n}: {line!r}")

    r1 = a1._compile_problems(source, "two_bad.py")
    r2 = a2._offenders(probe)
    print()
    print(f"a1 reports ({len(r1)}):")
    for item in r1:
        print(f"  {item}")
    print(f"a2 reports ({len(r2)}):")
    for item in r2:
        print(f"  {item}")
