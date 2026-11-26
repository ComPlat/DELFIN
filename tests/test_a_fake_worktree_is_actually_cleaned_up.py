"""A fake worktree factory that says it cleans up must actually clean up.

``tests/test_bench_ab.py::_factory`` builds a real directory per side, because
``run()`` derives paths from ``side.path``. It also carries a ``cleanup()``
with a careful containment guard whose docstring promises it removes the
directory the helper made.

Nothing called it. And its guard can never pass under the suite either:
it requires the directory's parent to be named ``tmp*``, while
``tests/conftest.py``'s session fixture points ``tempfile.tempdir`` at pytest's
own basetemp, whose parent is ``pytest-of-<user>``. So the promise was false in
two ways at once -- uncalled and unsatisfiable -- and each ``_run()`` left its
directories behind.

Written from the outside on purpose: it drives the helper the way ``run()``
does and then asks the filesystem, so a future rewrite of the guard is measured
by what it leaves on disk rather than by what its docstring says.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import tests.test_bench_ab as T


def test_the_fake_factory_removes_what_it_created():
    """The property the docstring states, measured on the filesystem."""
    factory = T._factory()
    info = factory(Path("."), "0" * 40)
    root = Path(factory.root)

    assert Path(info.path).is_dir(), "the factory made no directory to clean up"
    factory.cleanup()
    assert not root.exists(), (
        "the fake worktree root survived cleanup(), so every _run() in "
        "tests/test_bench_ab.py leaves its directories behind: " + str(root))


def test_cleanup_refuses_a_directory_it_did_not_make(tmp_path):
    """The containment the guard exists for survives the fix.

    An earlier version of this helper registered ``rmtree`` on a path a caller
    supplied and deleted this checkout. Whatever changes, a path the helper did
    not create is not one it may remove.
    """
    factory = T._factory()
    bystander = tmp_path / "not-mine"
    bystander.mkdir()

    factory.cleanup(target=bystander)

    assert bystander.is_dir(), "cleanup() removed a directory it did not make"


def test_cleanup_leaves_a_directory_that_merely_shares_the_prefix(tmp_path):
    """A name-based guard is not enough: the root is the one that was made."""
    factory = T._factory()
    decoy = Path(tempfile.mkdtemp(prefix="ab-fake-wt-"))
    try:
        factory.cleanup(target=decoy)
        assert decoy.is_dir(), (
            "cleanup() removed another helper's directory because it carried "
            "the same prefix")
    finally:
        import shutil

        shutil.rmtree(decoy, ignore_errors=True)
