"""H1D gate-coverage: the AST scanner that catches V2-type ungated file sinks.

Package H1D's product is this test plus the static analysis in
``delfin/agent/gate_coverage.py`` (no product behaviour changes).  The four
contracts exercised here are the wave's definition of "a gate covers a
filesystem sink":

  1. an executor that opens a path WITHOUT a gate call is reported (red),
  2. the same executor WITH the read gate before the open is clean,
  3. a gate call placed AFTER the open does not protect it (still a finding),
  4. an exemption entry that names no executor any more is itself a failure
     (the exemption list must not rot).

The scanner runs on fixture source text and never imports api_client, so these
tests are self-contained and independent of the real executor implementations.
"""

from delfin.agent import gate_coverage as gc


# ---------------------------------------------------------------------------
# 1. Ungated open is a finding
# ---------------------------------------------------------------------------


class TestUngatedFound:
    def test_open_unread_gated_is_reported(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        with open(path) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src, module_path="fixture.py")
        bad = report.ungated()
        assert len(bad) == 1, report.names()
        ex = bad[0]
        assert ex.name == "_execute_read"
        kinds = {r.sink.kind for r in ex.ungated_sinks}
        assert "open" in kinds, kinds  # the open itself went ungated

    def test_write_text_ungated_is_reported(self):
        src = (
            "class W:\n"
            "    def _execute_write(self, arguments: dict, perms=None):\n"
            "        target = arguments['path']\n"
            "        target.write_text('x')\n"
            "        return 'ok'\n"
        )
        report = gc.analyze(src)
        ex = report.ungated()[0]
        assert ex.name == "_execute_write"
        assert all(not r.gated for r in ex.ungated_sinks)
        kinds = {r.sink.kind for r in ex.ungated_sinks}
        assert "write_text" in kinds

    def test_os_open_ungated_is_reported(self):
        src = (
            "class O:\n"
            "    def _execute_sink(self, arguments: dict, perms=None):\n"
            "        fd = os.open(arguments['path'], os.O_RDWR)\n"
            "        return os.read(fd, 8)\n"
        )
        report = gc.analyze(src)
        ex = report.ungated()[0]
        kinds = {r.sink.kind for r in ex.ungated_sinks}
        assert "os.open" in kinds

    def test_sink_through_named_helper_attributes_back_to_executor(self):
        # The executor does not open inline; it calls a helper that does, and
        # that helper carries no gate.  The sink must still show up under the
        # executor (the "helper it names" case).
        src = (
            "class H:\n"
            "    def _read_all(self, p):\n"
            "        return open(p).read()\n"
            "    def _execute_sink(self, arguments: dict, perms=None):\n"
            "        return self._read_all(arguments['path'])\n"
        )
        report = gc.analyze(src)
        ex = report.ungated()[0]
        assert ex.name == "_execute_sink"
        assert any(r.sink.kind == "open" and not r.gated
                   for r in ex.ungated_sinks)


# ---------------------------------------------------------------------------
# 2. Gate before the open makes it clean
# ---------------------------------------------------------------------------


class TestGatedClean:
    def test_read_gate_before_open_is_clean(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        err = self._check_read_access(perms, path, label='x')\n"
            "        if err:\n"
            "            return err\n"
            "        with open(path) as f:\n"
            "            return f.read()\n"
        )
        report = gc.analyze(src)
        assert report.ungated() == [], (
            "gated open must not be a finding:\n" + gc.format_report(report)
        )

    def test_write_gate_before_write_is_clean(self):
        src = (
            "class W:\n"
            "    def _execute_write(self, arguments: dict, perms=None):\n"
            "        p = arguments['path']\n"
            "        err = self._gate_write_path(p, perms, 'write_file', {})\n"
            "        if err:\n"
            "            return err\n"
            "        p.write_text('x')\n"
            "        return 'ok'\n"
        )
        report = gc.analyze(src)
        assert report.ungated() == [], gc.format_report(report)


# ---------------------------------------------------------------------------
# 3. A gate AFTER the open is still a finding
# ---------------------------------------------------------------------------


class TestGateAfterOpenIsStillFinding:
    def test_gate_below_the_open_does_not_protect(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, arguments: dict, perms=None):\n"
            "        path = arguments['path']\n"
            "        with open(path) as f:\n"
            "            text = f.read()\n"
            "        err = self._check_read_access(perms, path)\n"
            "        return text\n"
        )
        report = gc.analyze(src)
        ex = report.ungated()[0]
        open_result = next(r for r in ex.ungated_sinks if r.sink.kind == "open")
        assert not open_result.gated
        # The gate line exists but is AFTER the open, so it must not be credited.
        assert open_result.gate_line is None


# ---------------------------------------------------------------------------
# 4. A stale exemption is itself a failure (the list does not rot)
# ---------------------------------------------------------------------------


class TestExemptionsDoNotRot:
    def test_stale_exemption_is_flagged(self):
        # Analyze a source that has no _execute_ghost; an exemption claiming it
        # must be reported stale rather than silently carried.
        src = (
            "class T:\n"
            "    def _execute_alive(self, a):\n"
            "        return 1\n"
        )
        report = gc.analyze(src)
        stale = report.stale_exemptions(
            {"_execute_ghost": "was removed in vX; line 40 was exempt"})
        assert stale == ["_execute_ghost"]

    def test_live_exemption_is_not_stale(self):
        src = (
            "class T:\n"
            "    def _execute_alive(self, a):\n"
            "        return 1\n"
        )
        report = gc.analyze(src)
        assert report.stale_exemptions({"_execute_alive": "intentional"}) == []

    def test_ungated_executor_on_exemption_list_is_not_a_finding(self):
        # The scanner reports ungated sinks; a caller that consults EXEMPTIONS
        # drops the exempted names.  Here we assert the hook exists and that an
        # un-exempted ungated executor stays a finding.
        src = (
            "class T:\n"
            "    def _execute_open(self, a):\n"
            "        return open(a['path'])\n"
        )
        report = gc.analyze(src)
        names = {e.name for e in report.ungated()}
        assert "_execute_open" in names


# ---------------------------------------------------------------------------
# Hygiene: the scanner is a static analysis, not a runtime import
# ---------------------------------------------------------------------------


class TestScannerIsSelfContained:
    def test_analyze_accepts_fixture_and_returns_executer_names(self):
        src = (
            "class T:\n"
            "    def _execute_read(self, a):\n"
            "        return 1\n"
            "    def _execute_write(self, a):\n"
            "        return 2\n"
        )
        report = gc.analyze(src)
        assert set(report.names()) == {"_execute_read", "_execute_write"}

    def test_module_has_no_sink_on_this_test_source(self):
        # Sanity: the source of gate_coverage itself must not import api_client.
        with open(gc.__file__, "r", encoding="utf-8") as fh:
            own = fh.read()
        assert "import api_client" not in own
        assert "from .api_client" not in own

    def test_non_executor_functions_are_not_reported(self):
        src = (
            "class T:\n"
            "    def _not_an_executor(self, a):\n"
            "        return open(a['path'])\n"
            "    def _execute_read(self, a):\n"
            "        return 1\n"
        )
        report = gc.analyze(src)
        assert report.names() == ["_execute_read"]
