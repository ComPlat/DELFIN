"""H1D gate-coverage, round 2: control-flow and path-argument dominance.

The first round (test_h1d_gate_coverage.py) proved line-order dominance was
*not* enough and the four basic contracts.  Two reviewers (nacht-s22,
nacht-s23) independently found the exact V2-shaped hole the package exists to
catch: the scanner credited a gate that does not actually run before the sink.

This file pins the fix.  A gate G dominates a filesystem sink S only when

  1. G is a **statement-sibling** of S in the same straight-line block (a gate
     inside ``if perms is not None:`` is a *different* statement block from an
     ``open`` in the sibling ``else``, so it never dominates it), and
  2. G appears *before* S in that block, and
  3. G's path argument covers the path S touches (a gate on ``a`` does not
     gate an ``open(b)``), and
  4. a ``subprocess.*`` call with a path as a *keyword* argument (``cwd=``,
     ``executable=``) is a real filesystem sink, not only one passed
     positionally.

Every fixture is self-contained source text -- the scanner never imports
api_client.  All cases here are the "must be a finding" (red) direction or the
"must stay clean" (green) mirror, so a regression in either direction fails
this file.  Written to fail on the round-1 scanner (tree@3289650f) and to pass
after the fix.
"""
from delfin.agent import gate_coverage as gc


# ---------------------------------------------------------------------------
# A gate in one branch does not dominate a sink in a sibling branch
# ---------------------------------------------------------------------------


class TestBranchConditionalGateDoesNotCover:
    def test_gate_in_if_branch_but_open_in_else_is_a_finding(self):
        # Reviewer s22 fixture: the gate runs only when perms is not None; the
        # open runs in the else, exactly when the gate does not.
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        if perms is not None:\n"
            "            self._check_read_access(perms, path)\n"
            "        else:\n"
            "            with open(path) as f:\n"
            "                return f.read()\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        open_result = next(r for r in ex.ungated_sinks if r.sink.kind == "open")
        assert not open_result.gated, open_result
        assert open_result.gate_line is None

    def test_gate_inside_if_does_not_cover_unconditional_write(self):
        # Reviewer s23 fixture: a gate guarded by `if force:` must not mark an
        # unconditional write_gated -- when force=False the write runs ungated.
        src = (
            "import pathlib\n"
            "class T:\n"
            "    def _execute_write(self, arguments: dict, force=False):\n"
            "        target = pathlib.Path(arguments['dst'])\n"
            "        if force:\n"
            "            self._gate_write_path(arguments, target)\n"
            "        target.write_text('user data')\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        wsink = next(r for r in ex.ungated_sinks if r.sink.kind == "write_text")
        assert not wsink.gated, wsink
        assert wsink.gate_line is None

    def test_gate_inside_try_does_not_cover_sink_after_except(self):
        src = (
            "class T:\n"
            "    def _execute_write(self, arguments: dict):\n"
            "        target = arguments['dst']\n"
            "        try:\n"
            "            self._gate_write_path(arguments, target)\n"
            "        except Exception:\n"
            "            pass\n"
            "        target.write_text('user data')\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        wsink = next(r for r in ex.ungated_sinks if r.sink.kind == "write_text")
        assert not wsink.gated, wsink
        assert wsink.gate_line is None

    def test_gate_in_one_else_branch_does_not_cover_the_other(self):
        # Symmetric: both sink and gate sit in two separate branches of the
        # same if/else; neither can dominate the other.
        src = (
            "class T:\n"
            "    def _execute_x(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        if perms is not None:\n"
            "            with open(path) as f:\n"
            "                return f.read()\n"
            "        else:\n"
            "            err = self._check_read_access(perms, path)\n"
            "            return None\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        open_result = next(r for r in ex.ungated_sinks if r.sink.kind == "open")
        assert not open_result.gated, open_result


class TestBranchGateSameBlockStillGates:
    def test_gate_and_open_in_same_if_body_are_clean(self):
        # A gate and its sink in the SAME straight-line block (both inside the
        # same `if` body) run together, so the gate still dominates.
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        if True:\n"
            "            self._check_read_access(perms, path)\n"
            "            with open(path) as f:\n"
            "                return f.read()\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated() == [], report.ungated()


# ---------------------------------------------------------------------------
# A gate on one path does not cover a sink on a different path
# ---------------------------------------------------------------------------


class TestWrongPathGate:
    def test_gate_on_different_path_than_open_is_a_finding(self):
        # Reviewer s22 fixture: gate checks `a`, the open touches `b`.  Earlier
        # line is not enough; the path argument must match.
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        a = arguments['a']\n"
            "        b = arguments['b']\n"
            "        err = self._check_read_access(perms, a)\n"
            "        with open(b) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        open_result = next(r for r in ex.ungated_sinks if r.sink.kind == "open")
        assert not open_result.gated, open_result
        assert open_result.gate_line is None

    def test_gate_on_same_path_is_clean(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        err = self._check_read_access(perms, path)\n"
            "        with open(path) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated() == [], report.ungated()


# ---------------------------------------------------------------------------
# subprocess path arguments passed as keywords are filesystem sinks
# ---------------------------------------------------------------------------


class TestSubprocessKeywordPath:
    def test_cwd_keyword_path_without_gate_is_a_finding(self):
        src = (
            "import subprocess\n"
            "class T:\n"
            "    def _execute_run(self, arguments: dict):\n"
            "        out = subprocess.run(\n"
            "            arguments['cmd'], cwd=arguments['cwd'],\n"
            "            capture_output=True)\n"
            "        return out\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        sp = next(r for r in ex.ungated_sinks if r.sink.kind.startswith("subprocess"))
        assert not sp.gated, sp
        assert sp.gate_line is None

    def test_executable_keyword_path_without_gate_is_a_finding(self):
        src = (
            "import subprocess\n"
            "class T:\n"
            "    def _execute_run(self, arguments: dict):\n"
            "        return subprocess.run(\n"
            "            arguments['cmd'], executable=arguments['bin'])\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        ex = report.ungated()[0]
        sp = next(r for r in ex.ungated_sinks if r.sink.kind.startswith("subprocess"))
        assert not sp.gated, sp

    def test_subprocess_keyword_path_with_gate_is_clean(self):
        src = (
            "import subprocess\n"
            "class T:\n"
            "    def _execute_run(self, arguments: dict, perms=None):\n"
            "        self._check_read_access(perms, arguments['cwd'])\n"
            "        return subprocess.run(\n"
            "            arguments['cmd'], cwd=arguments['cwd'])\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated() == [], report.ungated()


# ---------------------------------------------------------------------------
# Round-1 contracts must survive the fix (regression guard)
# ---------------------------------------------------------------------------


class TestRound1StillHolds:
    def test_ungated_open_is_found(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        with open(arguments['path']) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src)
        assert len(report.ungated()) == 1

    def test_gate_before_open_is_clean(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        err = self._check_read_access(perms, arguments['path'])\n"
            "        with open(arguments['path']) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src)
        assert report.ungated() == []

    def test_gate_after_open_is_still_a_finding(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        with open(arguments['path']) as f:\n"
            "            text = f.read()\n"
            "        err = self._check_read_access(perms, arguments['path'])\n"
            "        return text\n"
        )
        report = gc.analyze(src)
        open_result = next(
            r for r in report.ungated()[0].ungated_sinks if r.sink.kind == "open")
        assert not open_result.gated
        assert open_result.gate_line is None
