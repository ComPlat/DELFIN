"""H1D: tempfile/mkstemp staging writes are NOT ungated user-data sinks.

Red-direction: a write to a ``tempfile.mkstemp``/``NamedTemporaryFile`` handle is
transient staging (unlinked in ``finally`` or owned by the tempfile module), not
a persistent user-data write, so it must NOT be reported as an ungated sink.
Green-direction: a ``write_text`` to a user-supplied path under a gate stays
reported when ungated and stays silent when gated (no over-exemption).
"""
from delfin.agent import gate_coverage as gc


def _one(*lines):
    return "\n".join(lines) + "\n"


class TestTempfileStagingNotASink:
    def test_mkstemp_staging_write_is_not_reported(self):
        src = _one(
            "import tempfile, os",
            "class T:",
            "    def _execute_stage(self, arguments, perms=None):",
            "        _fd, _tmp = tempfile.mkstemp(suffix='.txt')",
            "        with os.fdopen(_fd, 'w') as _fh:",
            "            _fh.write('staged')",
            "        try:",
            "            return 'staged'",
            "        finally:",
            "            os.unlink(_tmp)",
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated() == [], report.ungated()

    def test_named_temporary_file_write_is_not_reported(self):
        src = _one(
            "import tempfile",
            "class T:",
            "    def _execute_stage(self, arguments, perms=None):",
            "        with tempfile.NamedTemporaryFile(mode='w') as fh:",
            "            fh.write('x')",
            "        return 'ok'",
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated() == [], report.ungated()


class TestRealWriteStillReported:
    def test_ungated_user_write_to_real_path_is_reported(self):
        # A guarded control: a write to a user path WITHOUT a gate is a finding.
        src = _one(
            "class T:",
            "    def _execute_write(self, arguments, perms=None):",
            "        target = arguments['dst']",
            "        with open(target, 'w') as f:",
            "            f.write('data')",
        )
        report = gc.analyze(src, module_path="fixture.py")
        assert report.ungated(), "a real ungated write must still be reported"
