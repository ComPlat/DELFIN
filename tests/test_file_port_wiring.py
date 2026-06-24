"""Tests for deep file-port auto-wiring (Track 2).

A parser declares ``wires = {capability: kwarg}``; when an upstream step makes
that capability available the pipeline fills the kwarg automatically.  QM engines
write their log to ``StepResult.output_file``, which the executor now exposes
under the ``qc_output`` tag, so a parser's ``filepath`` wires from the upstream
job without a manual path.
"""

from __future__ import annotations

import tempfile
import time
from pathlib import Path

from delfin.tools import Pipeline, platform
from delfin.tools._base import StepAdapter
from delfin.tools._registry import get, register
from delfin.tools._types import StepStatus
from delfin.tools._wiring import autowire_kwargs, wire_satisfies


# --- contract + static validation -----------------------------------------


def test_parsers_declare_wires():
    assert get("cclib_parse").contract().wires == (("qc_output", "filepath"),)
    assert get("uv_vis_parse").contract().wires == (("qc_output", "output_file"),)
    assert get("ir_parse").contract().wires == (("qc_output", "output_file"),)


def test_wire_satisfies_logic():
    wires = (("qc_output", "filepath"),)
    assert wire_satisfies("filepath", wires, frozenset({"qc_output"}))
    assert not wire_satisfies("filepath", wires, frozenset())
    assert not wire_satisfies("other", wires, frozenset({"qc_output"}))


def test_validation_wires_parser_filepath_from_upstream():
    # orca_sp produces qc_output → cclib_parse's filepath auto-wires → clean
    ok = platform.validate_spec({"name": "w", "steps": [
        {"step": "orca_sp", "kwargs": {"charge": 0}},
        {"step": "cclib_parse"}]}, geometry=True)
    assert ok.ok, ok.summary()
    # without the producer it is still reported missing
    bad = platform.validate_spec({"name": "w", "steps": [{"step": "cclib_parse"}]},
                                 geometry=False)
    assert any("filepath" in d.missing_params for d in bad.diagnostics)


def test_autowire_injects_from_wires():
    kw = autowire_kwargs("cclib_parse", {}, {"qc_output": "/calc/job.out"},
                         wires=(("qc_output", "filepath"),))
    assert kw["filepath"] == "/calc/job.out"
    # explicit value always wins over the wire
    kw2 = autowire_kwargs("cclib_parse", {"filepath": "/explicit.out"},
                          {"qc_output": "/calc/job.out"},
                          wires=(("qc_output", "filepath"),))
    assert kw2["filepath"] == "/explicit.out"


def test_manifest_exposes_wires():
    from delfin.tools.manifest import contract_to_dict
    d = contract_to_dict(get("cclib_parse").contract())
    assert {"capability": "qc_output", "kwarg": "filepath"} in d["wires"]


# --- end-to-end executor plumbing (fake adapters) -------------------------


class _FakeQC(StepAdapter):
    name = "_fake_qc_producer"
    description = "fake QM step: writes a .out into output_file, produces qc_output"
    produces_geometry = False
    category = "meta"
    produces = ("qc_output",)

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        out = work_dir / "job.out"
        out.write_text("dummy qc output")
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir,
                                 time.monotonic(), output_file=out)


class _FakeParser(StepAdapter):
    name = "_fake_parser_consumer"
    description = "fake parser: records the filepath it was given"
    produces_geometry = False
    category = "meta"
    consumes = ("qc_output",)
    wires = {"qc_output": "filepath"}

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir,
                                 time.monotonic(), data={"got_filepath": kwargs.get("filepath")})


register(_FakeQC())
register(_FakeParser())


def test_executor_exposes_output_file_as_qc_output_and_wires_it():
    d = Path(tempfile.mkdtemp())
    pipe = Pipeline("wire_chain")
    pipe.add("_fake_qc_producer", label="qc")
    pipe.add("_fake_parser_consumer", label="parse")
    res = pipe.run(work_dir=d)
    assert res.ok
    parser_result = res.results[-1]
    # the parser received the producer's output_file as its filepath, auto-wired
    assert parser_result.data["got_filepath"] is not None
    assert parser_result.data["got_filepath"].endswith("job.out")
