"""`lab run` is ONE verb over the two lab packages, and the judge must not
refuse a document the coordinator runs.

Two packages landed together and each added `delfin-agent lab run` on its own
branch. Measured on the merged tree before this test existed, through the
judge the CLI now calls:

    plan.validate(campaign_manifest, chemasap_campaign_workflow, kernel)
    -> ok: false, 16 problems

    op 'hold' is not one of set, wait, transfer, measure      (8 steps)
    sample S-001 is declared in 'storage', not 'synthesis'    (4 steps)
    sample S-001 is declared in 'storage', not 'purification' (4 steps)

All 16 were wrong. `hold` is an operation `delfin.lab.workflow` executes and
the four-sample campaign is built on; and a sample's declared `module` is
where it STARTS, so a second-stage hand-over's `from` is the previous
transfer's `to`, never the declaration. A judge that refuses a document the
coordinator runs correctly teaches a planner to distrust the one component
that must be trusted -- the exact defect `delfin/lab/plan.py` exists to fix,
arriving from the other direction once two packages met.

Each assertion below fails for its own reason, on the unpatched bundle:

* `test_the_shipped_campaign_is_judged_runnable` -- 16 problems, listed;
* `test_a_second_stage_hand_over_is_not_refused` -- the four `to-column-*`
  transfers are refused for a stage the sample already left;
* `test_a_hold_is_an_operation_and_its_duration_is_checked` -- `hold` is not
  an op at all, so every dwell in the campaign is refused as unknown, and a
  negative duration would be accepted once it is;
* `test_the_cli_runs_the_campaign_through_one_verb` -- drives the product
  entry point, `delfin.agent.cli.main(["lab", "run", ...])`, with `--dry-run`
  on the shipped 65-step campaign, and requires exit 0 with `ran: false`.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import _lab_support  # noqa: F401  (the fixture lives there)
from conftest import child_env
from _lab_support import lab_fixture  # noqa: F401  (registers the fixture)

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "delfin" / "lab" / "examples"
MANIFEST = EXAMPLES / "chemasap_campaign.yaml"
CAMPAIGN = EXAMPLES / "chemasap_campaign_workflow.yaml"


def _manifest_and_document():
    from delfin.lab import manifest as _manifest
    return (_manifest.load(MANIFEST),
            yaml.safe_load(CAMPAIGN.read_text(encoding="utf-8")))


def _kernel(manifest):
    """The kernel the judge needs, with the campaign manifest pinned.

    A stub for the WRITE half only: this test is about the structural half
    (which operations exist, and where a sample is at each step), so the
    kernel's own verdict is accepted as given and never re-derived here.
    """

    class _Stub:
        def __init__(self, man):
            self.manifest = man

        def plan_check(self, _actor, steps):
            return {"ok": True, "touched_hardware": False,
                    "steps": [{"ok": True, "needs_approval": False,
                               "level": "none"} for _ in steps]}

    return _Stub(manifest)


def test_the_shipped_campaign_is_judged_runnable():
    """65 steps, four modules, and nothing to report."""
    from delfin.lab import plan
    manifest, doc = _manifest_and_document()
    proposal = plan.validate(manifest, doc, _kernel(manifest))
    assert proposal.ok is True, proposal.problems
    assert proposal.problems == []
    assert len(doc["steps"]) == 65


def test_a_second_stage_hand_over_is_not_refused():
    """The sample is where the previous hand-over left it, not where it began.

    Red on the unpatched judge: S-001 is declared in `storage`, so the
    transfer into `purification` (whose `from` is `synthesis`, correct after
    the first move) is refused with "is declared in 'storage', not
    'synthesis'". The control case is the FIRST transfer of the same sample,
    which the declaration does describe -- it is admitted on both trees, so
    the failure above is the staging, not a blanket refusal of transfers.
    """
    from delfin.lab import plan
    manifest, doc = _manifest_and_document()
    proposal = plan.validate(manifest, doc, _kernel(manifest))
    refused = [p for p in proposal.problems if "S-001" in p and "storage" in p]
    assert refused == [], refused

    transfers = [st for st in doc["steps"] if st.get("op") == "transfer"
                 and st.get("sample") == "S-001"]
    assert [t["from"] for t in transfers] == ["storage", "synthesis",
                                              "purification"], transfers
    # The control: the first hand-over's `from` IS the declared module, and a
    # document that got it wrong is still refused -- the check did not go away.
    doc["steps"] = [dict(transfers[0], **{"from": "analysis"})]
    wrong = plan.validate(manifest, doc, _kernel(manifest))
    assert any("S-001 is in 'storage'" in p for p in wrong.problems), \
        wrong.problems


def test_a_hold_is_an_operation_and_its_duration_is_checked():
    """`hold` runs in the coordinator, so the judge must know it -- and judge it."""
    from delfin.lab import plan
    assert "hold" in plan.OPS
    manifest, doc = _manifest_and_document()
    doc["steps"] = [{"id": "dwell", "op": "hold", "seconds": 30}]
    ok = plan.validate(manifest, doc, _kernel(manifest))
    assert ok.ok is True, ok.problems
    # A duration nobody can run is refused, with the op named.
    doc["steps"] = [{"id": "dwell", "op": "hold", "seconds": -5}]
    bad = plan.validate(manifest, doc, _kernel(manifest))
    assert bad.ok is False
    assert any("hold" not in p and "seconds must be >= 0" in p
               for p in bad.problems), bad.problems


def test_the_cli_runs_the_campaign_through_one_verb(lab, tmp_path):
    """The product entry point: `delfin-agent lab run WF --dry-run` exits 0.

    Runs in a child process so the CLI is reached the way a person reaches it,
    with `PYTHONPATH` naming THIS tree -- a test that imports a module and
    calls its function measures that the function exists.

    The child's environment comes from ``child_env``, not ``os.environ``: the
    suite's own sink redirects are in-process attribute swaps and stop at the
    process boundary, so a child started on the raw environment resolves
    ``~/.delfin`` against the real home. Measured on the bundle before this
    fix: the repository-wide scanner flagged this line, and the stated leak
    is nine confirmation requests from suite children in the operator's room.
    """
    from delfin.lab.safety import pin
    from delfin.lab import manifest as _manifest, store

    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MANIFEST.read_text(encoding="utf-8"), encoding="utf-8")
    pin(_manifest.load(path), by="alice")

    ws = tmp_path / "campaign"
    ws.mkdir()
    doc = ws / "campaign.yaml"
    doc.write_text(CAMPAIGN.read_text(encoding="utf-8"), encoding="utf-8")
    env = dict(child_env(tmp_path))
    env["PYTHONPATH"] = str(REPO)
    env.pop("DELFIN_AGENT_SESSION", None)          # this IS the person
    env["DELFIN_LAB_DIR"] = os.environ["DELFIN_LAB_DIR"]
    env["DELFIN_PERSON"] = "alice"
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", "lab", "run", str(doc),
         "--dry-run"],
        cwd=str(ws), env=env, text=True, capture_output=True, timeout=300)
    assert proc.returncode == 0, (proc.returncode, proc.stderr)
    out, _end = json.JSONDecoder().raw_decode(proc.stdout.lstrip())
    assert out["ok"] is True, out
    assert out["ran"] is False
    assert "nothing was written" in proc.stdout


def test_a_plan_is_not_ok_when_the_kernel_would_refuse_every_write(lab):
    """`ok: true` must mean "this would run", not "the steps parse".

    The kernel's dry run returns its own verdict AND the reason its writes
    would be refused -- on a fresh installation, `manifest_problem` is "the
    manifest is not pinned". The judge copied `write_check["steps"]` and
    dropped that field, so it reported `ok: true, problems == []` for a
    campaign every one of whose 16 writes the kernel refuses. An agent
    reading that answer would tell a person the campaign was ready to start.

    Measured on the bundle before this fix, with this manifest written and
    never pinned: `proposal.ok is True`, `problems == []`,
    `write_check["ok"] is False`, `write_check["manifest_problem"]` the
    pin message. The control case -- the same document with the manifest
    pinned by `lab` -- stays ok, so the fix reports the kernel's problem
    rather than refusing plans in general.
    """
    from delfin.lab import plan
    from delfin.lab import manifest as _manifest, store
    from delfin.lab.safety import Kernel, pin, pin_path

    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MANIFEST.read_text(encoding="utf-8"), encoding="utf-8")
    # The shared fixture pins the DEMO manifest; a fresh installation that has
    # reviewed this campaign and not yet pinned it is the state under test.
    pin_path().unlink(missing_ok=True)
    assert not pin_path().exists(), "the fixture must start unpinned"

    manifest = _manifest.load(path)
    doc = yaml.safe_load(CAMPAIGN.read_text(encoding="utf-8"))
    kernel = Kernel.load()

    # The kernel really does refuse: this is the precondition, not the claim.
    assert kernel.manifest is not None
    assert kernel._manifest_problem() is not None

    proposal = plan.validate(manifest, doc, kernel)
    assert proposal.write_check["ok"] is False
    assert proposal.ok is False, proposal.problems
    assert any("not pinned" in p for p in proposal.problems), proposal.problems

    # The control: pinned, the same document is judged runnable again.
    pin(manifest, by="alice")
    pinned_proposal = plan.validate(manifest, doc, Kernel.load())
    assert pinned_proposal.write_check["manifest_problem"] is None
    assert pinned_proposal.ok is True, pinned_proposal.problems
