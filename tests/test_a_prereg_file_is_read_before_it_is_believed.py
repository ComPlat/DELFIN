# SPDX-License-Identifier: MIT
"""A pre-registration file is believed only after its content is read.

``paired_compare`` reports whether an arm's pre-registration predates the
first stamp. Before this, existence and mtime were the whole check, so an
empty or half-written file counted as a pre-registration -- the rejected
alternative named in ``_prereg_verdict``'s own docstring. The content now
goes through the real gate in ``delfin.agent.experiment``.

Input: one prereg file per case. Output: the verdict line the report
carries. Each of the four outcomes is asserted, and so is the one thing a
test of the helper alone cannot show: that ``compare`` puts the verdict in
its output. Without that, this file would prove the helper works and say
nothing about whether the product uses it -- which is how a capability ends
up tested and unmounted at the same time.
"""
from __future__ import annotations

import json
import os
import time

from delfin.agent import paired_compare as pc


def _prereg(tmp_path, payload, *, raw: str | None = None):
    p = tmp_path / "prereg.json"
    p.write_text(raw if raw is not None else json.dumps(payload),
                 encoding="utf-8")
    return p


#: The fields the gate requires, from the module's own tests.
_COMPLETE = {
    "id": "e1",
    "hypothesis": "h the change matters",
    "why_chain": ["one judge per comparison"],
    "switch": "DELFIN_MODE",
    "expectation": "rise < 5%",
    "reading": "one judge",
    "pool_size": "small",
}


def test_a_complete_pre_registration_is_accepted(tmp_path):
    assert pc._prereg_verdict(_prereg(tmp_path, _COMPLETE)) == \
        "before_measurement"


def test_an_empty_file_is_not_a_pre_registration(tmp_path):
    """The defect this replaced: existence and mtime said yes."""
    verdict = pc._prereg_verdict(_prereg(tmp_path, {}))
    assert verdict.startswith("not_written_down: "), verdict
    assert verdict != "not_written_down: ", "the refusal carries no reason"


def test_a_half_written_pre_registration_names_what_is_missing(tmp_path):
    """Each required field removed on its own, so one missing field cannot
    be masked by another that happens to be checked first."""
    for field in ("hypothesis", "why_chain", "switch", "expectation",
                  "reading"):
        payload = {k: v for k, v in _COMPLETE.items() if k != field}
        verdict = pc._prereg_verdict(_prereg(tmp_path, payload))
        assert verdict.startswith("not_written_down: "), (field, verdict)


def test_a_file_that_is_not_an_object_is_refused_by_shape(tmp_path):
    """A JSON array parses, so the parse alone does not answer the
    question. Reported as its own reason rather than as unreadable."""
    verdict = pc._prereg_verdict(_prereg(tmp_path, None, raw="[1, 2, 3]"))
    assert verdict == "not_written_down: the file is not a JSON object"


def test_an_unreadable_file_names_the_failure_type(tmp_path):
    missing = tmp_path / "absent.json"
    verdict = pc._prereg_verdict(missing)
    assert verdict.startswith("not_written_down: unreadable file ("), verdict
    torn = _prereg(tmp_path, None, raw='{"id": "e1"')
    assert pc._prereg_verdict(torn).startswith(
        "not_written_down: unreadable file ("), "a torn file is unreadable"


def _stamp(**over):
    # taken_at matters: compare takes the EARLIER of the two stamps and
    # refuses a prereg file whose mtime is later, before any content is
    # read. Without it both default to 0 and every real file reads as
    # written afterwards -- which is what driving the product showed and a
    # test of the helper alone never would.
    base = {"code": {}, "command": "", "prefixes": ["DELFIN_", "PROJ_"],
            "env_source": "command", "env": {},
            "taken_at": time.time() + 60}
    base.update(over)
    return base


def _arm(n=120, *, extra=lambda i: {}):
    return [{"id": f"i{i:03d}", "ok": True, "value": i, **extra(i)}
            for i in range(n)]


def test_the_verdict_reaches_the_report(tmp_path):
    """The half a helper test cannot reach: ``compare`` carries it out.

    ``compare`` is driven for real, not inspected: it is reached from the
    command line (``cli.py``'s compare-arms entry calls ``_pc.compare``),
    so a verdict in its output is a verdict a person sees. Both directions,
    because one value alone would also appear if the field were hard-coded.
    """
    good = _prereg(tmp_path, _COMPLETE)
    bad = tmp_path / "empty.json"
    bad.write_text("{}", encoding="utf-8")

    def _run(path):
        return pc.compare(
            _arm(), _arm(extra=lambda i: {"value": i + (1 if i % 2 else 0)}),
            stamp_a=_stamp(env={"PROJ_FAST": "0"}),
            stamp_b=_stamp(env={"PROJ_FAST": "1"}),
            axis=["PROJ_FAST"], ok="ok", prereg=str(path), env={})

    out_good, out_bad = _run(good), _run(bad)
    # The two verdicts the wiring produces before the content is read, so
    # the four content verdicts above are not the whole reported set.
    assert _run(tmp_path / "absent.json").get("prereg") == "missing"
    late = _prereg(tmp_path, _COMPLETE)
    os.utime(late, (time.time() + 600, time.time() + 600))
    assert _run(late).get("prereg") == \
        "written_after_the_first_arm_was_stamped"
    assert out_good.get("prereg") == "before_measurement", out_good.get("prereg")
    assert str(out_bad.get("prereg", "")).startswith("not_written_down: "), (
        out_bad.get("prereg"))
    assert out_good["prereg"] != out_bad["prereg"], (
        "the field does not depend on the file, so it is not being read")
