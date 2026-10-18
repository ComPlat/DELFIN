#!/usr/bin/env python3
"""Setup for w17_code_flaky_test: a test suite with one flaky case.

``slowweb/`` -- one commit on main. ``tests/test_fetch.py`` has a test
that passes most runs but fails when the machine is slow: it builds a
cache entry whose TTL is measured in wall-clock time (``time.time()``
with a 0.05 s TTL) and asserts the entry is expired -- on a fast run
the assertion runs within the TTL and the test fails sporadically.
The suite also has a genuinely stable test that must keep passing.
The fix the task wants: make the cache entry's expiry INJECTABLE (a
``now=`` parameter defaulting to time.time) so the test controls time
-- the smallest change that removes the dependence on machine speed.
Precondition: the stable test is green, the flaky test is green in
this run (flakiness cannot be asserted), and the flaky test contains
a direct ``time.time()`` call. argv[1] is the workspace.
"""
from __future__ import annotations

from w17_code_repo import (commit, expect, main, mark_base, new_repo,
                            pytest_rc, write)

FILES = {
    ".gitignore": "__pycache__/\n*.pyc\n.pytest_cache/\n",
    "README.md": (
        "# slowweb\n\nA tiny HTTP cache.\n\n"
        "```python\nfrom slowweb.cache import Cache\n```\n\n"
        "Run the tests with `python -m pytest`.\n"),
    "pyproject.toml": (
        "[project]\nname = \"slowweb\"\nversion = \"0.1.0\"\n"
        "requires-python = \">=3.9\"\n\n"
        "[tool.pytest.ini_options]\npythonpath = [\"src\"]\n"
        "testpaths = [\"tests\"]\n"),
    "src/slowweb/__init__.py": (
        '"""slowweb: a tiny HTTP cache."""\n\n'
        "from .cache import Cache\n\n"
        '__all__ = ["Cache"]\n'),
    "src/slowweb/cache.py": (
        '"""A TTL cache."""\n\n'
        "import time\n\n\n"
        "class Cache:\n"
        '    """Entries expire after ``ttl`` seconds of wall time."""\n\n'
        "    def __init__(self, ttl: float = 60.0) -> None:\n"
        "        self.ttl = ttl\n"
        "        self._entries = {}\n\n"
        "    def put(self, key: str, value: str) -> None:\n"
        '        """Store an entry stamped with the current time."""\n'
        "        self._entries[key] = (value, time.time())\n\n"
        "    def get(self, key: str):\n"
        '        """The value if fresh, else None."""\n'
        "        entry = self._entries.get(key)\n"
        "        if entry is None:\n"
        "            return None\n"
        "        value, stamp = entry\n"
        "        if time.time() - stamp > self.ttl:\n"
        "            return None\n"
        "        return value\n"),
    "tests/test_fetch.py": (
        "import time\n\nfrom slowweb import Cache\n\n\n"
        "def test_hit_returns_value():\n"
        "    c = Cache(ttl=60)\n"
        "    c.put(\"k\", \"v\")\n"
        "    assert c.get(\"k\") == \"v\"\n\n\n"
        "def test_expired_entry_is_gone():\n"
        "    c = Cache(ttl=0.05)\n"
        "    c.put(\"k\", \"v\")\n"
        "    time.sleep(0.06)\n"
        "    assert c.get(\"k\") is None\n"),
}


def build(workspace: str) -> None:
    repo = new_repo(workspace, "slowweb")
    write(repo, FILES)
    commit(repo, "Initial slowweb")
    mark_base(repo)
    rc, tail = pytest_rc(repo, "tests/test_fetch.py")
    expect("suite green this run", rc == 0, tail)
    expect("the flaky test leans on wall time",
           "time.sleep(0.06)" in FILES["tests/test_fetch.py"]
           and "ttl=0.05" in FILES["tests/test_fetch.py"])


if __name__ == "__main__":
    main(build)
