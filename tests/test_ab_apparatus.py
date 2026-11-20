"""The apparatus is an instrument, so it is driven, not trusted.

Input: generated miniature repositories, one per outcome the scorer has to
tell apart. Output: the usual pass or fail.

Semantics: every test here fails when the check it covers is deleted from
``scripts/ab_score.py`` or ``scripts/ab_task_pools.py``. That was verified by
removing each check and running the test, which is the only evidence that a
test is attached to the behaviour it names rather than to a neighbouring one.

The defect class this guards against is an instrument that reports success
without doing its job: a scan whose baseline lists the names it looks for, a
rate with no absolute count beside it, a check that never fires. So each
exclusion, each signal and each of the scorer's four checks is driven to
both verdicts, never only to the one it is supposed to return.

The five scenarios of the brief, each a repository:

* ``connected``      -- a product caller now names the capability;
* ``deleted``        -- the capability and its tests are gone;
* ``deleted_tests_left`` -- the capability is gone, its test is not;
* ``instrument``     -- the capability is connected and the scoring
                        instrument was edited;
* ``skipped``        -- the capability is connected and a test that passed
                        before now skips.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import ab_draw_sample as draw_mod  # noqa: E402
import ab_score as score_mod  # noqa: E402
import ab_task_pools as pools_mod  # noqa: E402

CAP = "orphan_capability"


# --------------------------------------------------------------------------
# miniature repositories
# --------------------------------------------------------------------------

def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=apparatus",
         "-c", "user.email=apparatus@example.invalid",
         "-c", "commit.gpgsign=false", *args],
        capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    return out.stdout


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _base_repo(tmp_path: Path) -> tuple[Path, str]:
    """A repository holding one unreachable capability. Returns root, sha."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _write(repo / "conftest.py",
           "import os, sys\n"
           "sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\n")
    _write(repo / "delfin" / "__init__.py", "")
    _write(repo / "delfin" / "agent" / "__init__.py", "")
    _write(repo / "delfin" / "agent" / "widget.py",
           "def orphan_capability(flag):\n"
           "    return 'yes' if flag else 'no'\n\n\n"
           "def used_helper():\n"
           "    return 1\n")
    _write(repo / "delfin" / "agent" / "caller.py",
           "from delfin.agent.widget import used_helper\n\n\n"
           "def run():\n"
           "    return used_helper()\n")
    _write(repo / "tests" / "test_widget.py",
           "from delfin.agent.widget import orphan_capability\n\n\n"
           "def test_orphan_capability_says_yes():\n"
           "    assert orphan_capability(True) == 'yes'\n")
    _write(repo / "tests" / "test_caller.py",
           "from delfin.agent.caller import run\n\n\n"
           "def test_run():\n"
           "    assert run() == 1\n")
    _write(repo / "tests" / "fixtures" / "unmounted_capabilities.txt",
           "# Capabilities the product does not reach.\n"
           "test-only\torphan_capability\tdelfin/agent/widget.py:1\n")
    _write(repo / "tests" / "test_every_capability_has_a_caller.py",
           "def test_placeholder():\n    assert True\n")
    _write(repo / "scripts" / "ab_task_pools.py", "# stand-in\n")
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    return repo, _git(repo, "rev-parse", "HEAD").strip()


def _drop_capability(repo: Path) -> None:
    _write(repo / "delfin" / "agent" / "widget.py",
           "def used_helper():\n    return 1\n")


def _connect_capability(repo: Path) -> None:
    _write(repo / "delfin" / "agent" / "caller.py",
           "from delfin.agent.widget import orphan_capability, used_helper\n"
           "\n\ndef run():\n"
           "    return used_helper() if orphan_capability(True) else 0\n")


def _drop_from_fixture(repo: Path) -> None:
    path = repo / "tests" / "fixtures" / "unmounted_capabilities.txt"
    kept = [ln for ln in path.read_text().splitlines()
            if f"\t{CAP}\t" not in ln]
    path.write_text("\n".join(kept) + "\n")


def _attempt(tmp_path: Path, scenario: str) -> dict:
    """Build the pre tree, the post tree and the attempt directory.

    Output: ``{"attempt": dir, "pre": tree, "post": tree}``.
    """
    repo, pre_sha = _base_repo(tmp_path)

    if scenario == "connected":
        _connect_capability(repo)
        _drop_from_fixture(repo)
    elif scenario == "deleted":
        _drop_capability(repo)
        (repo / "tests" / "test_widget.py").unlink()
        _drop_from_fixture(repo)
    elif scenario == "deleted_tests_left":
        _drop_capability(repo)
        _drop_from_fixture(repo)
    elif scenario == "instrument":
        _connect_capability(repo)
        _drop_from_fixture(repo)
        _write(repo / "tests" / "test_every_capability_has_a_caller.py",
               "def test_placeholder():\n    assert True  # edited\n")
    elif scenario == "skipped":
        _connect_capability(repo)
        _drop_from_fixture(repo)
        _write(repo / "tests" / "test_caller.py",
               "import pytest\n\n"
               "from delfin.agent.caller import run\n\n\n"
               "def test_run():\n"
               "    pytest.skip('not today')\n"
               "    assert run() == 1\n")
    elif scenario == "string_only":
        # The hole the primary endpoint closes: the capability is named in a
        # string literal and never called.
        _write(repo / "delfin" / "agent" / "caller.py",
               "from delfin.agent.widget import used_helper\n\n\n"
               "DISPATCH = ('orphan_capability',)\n\n\n"
               "def run():\n    return used_helper()\n")
        _drop_from_fixture(repo)
    elif scenario == "self_call":
        # The hole Correction D closes: the orphan calls itself, which is
        # an ast.Call whose callee is the name and reaches nothing.
        # Behaviour is preserved on purpose, so that the pre-existing test
        # still passes and check 1 is the only check that refuses.
        _write(repo / "delfin" / "agent" / "widget.py",
               "def orphan_capability(flag):\n"
               "    if flag is None:\n"
               "        return orphan_capability(True)\n"
               "    return 'yes' if flag else 'no'\n\n\n"
               "def used_helper():\n"
               "    return 1\n")
        _drop_from_fixture(repo)
    elif scenario == "callback_only":
        # A genuine connection the primary endpoint cannot see: the
        # capability is handed to another function as a value and put in a
        # dispatch dict of function objects. No ast.Call names it.
        # `shadow` writes and deletes a local of the same spelling, which
        # is not a use of the capability: a store is not a load. Without
        # that clause the diagnosis would report two more sites.
        _write(repo / "delfin" / "agent" / "caller.py",
               "from delfin.agent.widget import orphan_capability, "
               "used_helper\n\n\n"
               "HANDLERS = {'orphan': orphan_capability}\n\n\n"
               "def shadow(value):\n"
               "    orphan_capability = value\n"
               "    del orphan_capability\n"
               "    return value\n\n\n"
               "def run():\n"
               "    handler = HANDLERS['orphan']\n"
               "    return used_helper() if handler else 0\n")
        _drop_from_fixture(repo)
    elif scenario == "scorer_edited":
        _connect_capability(repo)
        _drop_from_fixture(repo)
        _write(repo / "scripts" / "ab_task_pools.py", "# stand-in, edited\n")
    elif scenario == "fixture_line_added":
        # The shape that matters: the capability's own line is removed, as
        # a correct attempt must, and a line for something else is added.
        # Only the added-lines clause separates this from a correct edit.
        _connect_capability(repo)
        _drop_from_fixture(repo)
        path = repo / "tests" / "fixtures" / "unmounted_capabilities.txt"
        path.write_text(path.read_text() + "nowhere\tsomething_else\tx.py:1\n")
    elif scenario == "untouched":
        pass
    else:                                       # pragma: no cover
        raise AssertionError(scenario)

    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "--allow-empty", "-m", scenario)
    post_sha = _git(repo, "rev-parse", "HEAD").strip()

    pre = tmp_path / "pre"
    _git(repo, "worktree", "add", "-q", "--detach", str(pre), pre_sha)

    attempt = tmp_path / "7f3a9c1e"
    attempt.mkdir()
    (attempt / "manifest.json").write_text(json.dumps({
        "attempt_id": "7f3a9c1e",
        "pre_sha": pre_sha,
        "post_sha": post_sha,
        "started_at": "2026-11-20T10:00:00+00:00",
        "ended_at": "2026-11-20T10:12:00+00:00",
        "task": {"name": CAP, "file": "delfin/agent/widget.py", "line": 1,
                 "module": "delfin.agent.widget", "pool": "infra",
                 "tier": "test-only",
                 "task": "connect this capability to the product, or delete "
                         "it together with its tests"},
    }))
    (attempt / "log.jsonl").write_text(
        '{"ts": "2026-11-20T10:00:00+00:00", "role": "assistant", '
        '"input_tokens": 1200, "output_tokens": 300}\n'
        '{"ts": "2026-11-20T10:11:00+00:00", "role": "assistant", '
        '"usage": {"input_tokens": 1800, "output_tokens": 450}}\n')
    (attempt / "transcript.txt").write_text("not read by the scorer\n")
    return {"attempt": attempt, "pre": pre, "post": repo}


def _score(bundle: dict) -> dict:
    return score_mod.score_attempt(bundle["attempt"], bundle["post"],
                                   bundle["pre"], sys.executable, timeout=600)


# --------------------------------------------------------------------------
# the scorer, one scenario at a time
# --------------------------------------------------------------------------

def test_a_connected_capability_scores_one(tmp_path):
    rec = _score(_attempt(tmp_path, "connected"))
    assert rec["route"] == "connected"
    assert rec["score"] == 1, rec["reasons"]
    assert all(rec["checks"].values())


def test_a_string_literal_satisfies_the_secondary_endpoint_and_not_the_score(
        tmp_path):
    """The primary endpoint is the call site, and the two are scored apart.

    A string literal naming the capability is something the agent under
    test can write, so it may not be what the score rests on. It is still
    recorded, because a genuine dispatch table looks exactly like this and
    that case has to be visible rather than argued about.
    """
    rec = _score(_attempt(tmp_path, "string_only"))
    assert rec["route"] == "still_unreachable"
    assert rec["route_secondary"] == "connected"
    assert rec["variants_disagree"] is True
    assert rec["score"] == 0
    assert rec["score_secondary"] == 1
    assert rec["route_detail"]["call_sites"] == []
    assert rec["route_detail"]["mentions"] >= 1


def test_a_real_call_satisfies_both_endpoints_and_they_agree(tmp_path):
    rec = _score(_attempt(tmp_path, "connected"))
    assert rec["route"] == "connected"
    assert rec["route_secondary"] == "connected"
    assert rec["variants_disagree"] is False
    assert rec["score"] == 1 and rec["score_secondary"] == 1
    assert rec["route_detail"]["call_sites"]


def test_a_capability_that_only_calls_itself_is_not_connected(tmp_path):
    """Correction D: a function called by nothing but itself is not reached.

    A recursive call is an ``ast.Call`` whose callee is the name, so the
    primary endpoint read it as connected although nothing in the product
    reaches it -- the shape the endpoint exists to refuse, in the one place
    it did not look. Without the exclusion of the capability's own
    definition this attempt scores 1; with it the call site list is empty
    and the route is ``still_unreachable``.

    It cannot arise at baseline: a self-call counts as a mention and a
    mentioned name is never a candidate, so only an attempt can write it.
    """
    rec = _score(_attempt(tmp_path, "self_call"))
    assert rec["route_detail"]["call_sites"] == [], rec["route_detail"]
    assert rec["route"] == "still_unreachable"
    assert rec["score"] == 0
    # Check 1 is the only one that refuses: the behaviour is unchanged and
    # the pre-existing test still passes, so nothing else explains the zero.
    assert rec["checks"]["capability_resolved"] is False
    assert rec["checks"]["tests_still_pass"] is True, rec["tests"]
    assert rec["checks"]["instrument_untouched"] is True
    assert rec["checks"]["deletion_removed_its_tests"] is True
    # The string-inclusive reading counts the self-call as a mention, so
    # the two readings disagree here for the opposite reason to
    # `string_only`: the secondary is the lenient one in both cases.
    assert rec["route_secondary"] == "connected"
    assert rec["variants_disagree"] is True
    # The diagnosis does not rescue it either: a self-reference is inside
    # the definition and is not recorded.
    assert rec["ast_references"]["sites"] == []


def test_a_callback_connection_is_recorded_although_it_scores_zero(tmp_path):
    """Correction E: the third reading is a diagnosis, not a score.

    The capability is put in a dispatch dict of function objects and
    handed on as a value. That is a genuine connection and no ``ast.Call``
    names it, so the primary endpoint scores it zero. The attempt is not
    rescued -- the score stays 0 -- and the loads are recorded so the case
    is visible in the row instead of being argued about afterwards.
    """
    rec = _score(_attempt(tmp_path, "callback_only"))
    assert rec["route_detail"]["call_sites"] == []
    assert rec["route"] == "still_unreachable"
    assert rec["score"] == 0
    # Nothing here is a third score; the two scored readings are unchanged.
    assert rec["score_secondary"] == 1        # the import line is a mention
    assert set(rec["checks_secondary"]) == set(rec["checks"])
    refs = rec["ast_references"]
    # The dict value is a load. The name in the import statement is an
    # alias, and `shadow` only stores and deletes a local of the same
    # spelling; none of the three is a use, so the count is exactly one.
    assert refs["count"] == 1, refs
    assert refs["sites"] == refs["beyond_call_sites"]
    assert refs["sites"][0].startswith("delfin/agent/caller.py:")


def test_a_deletion_is_judged_the_same_way_under_both_readings(tmp_path):
    """Strictness in the deletion direction cannot be gamed, so it is not
    split: anything still naming a deleted function is a break."""
    rec = _score(_attempt(tmp_path, "deleted"))
    assert rec["route"] == rec["route_secondary"] == "deleted"
    assert rec["variants_disagree"] is False


def test_a_proper_deletion_scores_one_and_is_a_different_route(tmp_path):
    rec = _score(_attempt(tmp_path, "deleted"))
    assert rec["route"] == "deleted"
    assert rec["score"] == 1, rec["reasons"]
    # The two routes are not summed into one number.
    assert rec["route"] != "connected"


def test_a_deletion_that_left_its_tests_behind_scores_zero(tmp_path):
    rec = _score(_attempt(tmp_path, "deleted_tests_left"))
    assert rec["route"] == "deleted"
    assert rec["checks"]["deletion_removed_its_tests"] is False
    assert rec["score"] == 0
    assert any("left behind" in r for r in rec["reasons"])
    assert "tests/test_widget.py" in rec["deletion"]["left_behind"]


def test_touching_the_instrument_scores_zero_and_records_why(tmp_path):
    rec = _score(_attempt(tmp_path, "instrument"))
    assert rec["route"] == "connected"
    assert rec["checks"]["instrument_untouched"] is False
    assert rec["score"] == 0
    assert "tests/test_every_capability_has_a_caller.py" in \
        rec["instrument"]["touched_instrument"]
    # Recorded, not excluded: the attempt is still a scored row.
    assert "score" in rec and rec["reasons"]


def test_editing_the_scoring_scripts_scores_zero(tmp_path):
    rec = _score(_attempt(tmp_path, "scorer_edited"))
    assert rec["route"] == "connected"
    assert rec["checks"]["instrument_untouched"] is False
    assert rec["instrument"]["touched_instrument"] == [
        "scripts/ab_task_pools.py"]
    assert rec["score"] == 0


def test_a_test_that_turned_into_a_skip_is_not_a_pass(tmp_path):
    rec = _score(_attempt(tmp_path, "skipped"))
    assert rec["route"] == "connected"
    assert rec["tests"]["turned_into_skip"], rec["tests"]
    assert rec["checks"]["tests_still_pass"] is False
    assert rec["score"] == 0


def test_an_untouched_tree_scores_zero_without_being_dropped(tmp_path):
    rec = _score(_attempt(tmp_path, "untouched"))
    assert rec["route"] == "still_unreachable"
    assert rec["score"] == 0
    assert rec["reasons"] == ["route=still_unreachable"]


def test_a_line_added_to_the_baseline_is_not_an_allowed_fixture_edit(tmp_path):
    rec = _score(_attempt(tmp_path, "fixture_line_added"))
    assert rec["instrument"]["fixture_change"]["added_lines"] == 1
    assert rec["instrument"]["fixture_change"]["allowed"] is False
    assert rec["score"] == 0


def test_a_removal_only_fixture_edit_is_allowed(tmp_path):
    rec = _score(_attempt(tmp_path, "connected"))
    change = rec["instrument"]["fixture_change"]
    assert change["allowed"] is True
    assert change["removed_names"] == [CAP]
    assert change["added_lines"] == 0


# --------------------------------------------------------------------------
# the cost axis
# --------------------------------------------------------------------------

def test_cost_is_read_from_the_attempt_log(tmp_path):
    rec = _score(_attempt(tmp_path, "connected"))
    cost = rec["cost"]
    assert cost["input_tokens"] == 3000
    assert cost["output_tokens"] == 750
    assert cost["turns"] == 2
    assert cost["wall_seconds"] == 660.0


def test_a_missing_figure_is_null_and_never_zero(tmp_path):
    bundle = _attempt(tmp_path, "connected")
    (bundle["attempt"] / "log.jsonl").write_text(
        '{"ts": "2026-11-20T10:00:00+00:00", "role": "assistant"}\n')
    rec = _score(bundle)
    assert rec["cost"]["input_tokens"] is None
    assert rec["cost"]["output_tokens"] is None
    assert rec["cost"]["turns"] == 1


# --------------------------------------------------------------------------
# blindness to the arm
# --------------------------------------------------------------------------

def test_the_command_line_has_no_arm_option():
    strings = []
    for action in score_mod.build_parser()._actions:
        strings.extend(action.option_strings)
    joined = " ".join(strings).lower()
    for word in score_mod.ARM_WORDS:
        if word in ("model", "prompt"):
            continue            # --python is the interpreter, not a model
        assert word not in joined, word
    assert "--arm" not in joined


def test_the_scorer_refuses_an_unknown_manifest_key(tmp_path):
    """The key allowlist, driven with a key that is not arm vocabulary.

    Otherwise the arm-word check alone would carry the test and the
    allowlist could be deleted unnoticed.
    """
    bundle = _attempt(tmp_path, "connected")
    doc = json.loads((bundle["attempt"] / "manifest.json").read_text())
    doc["runner_notes"] = "anything at all"
    (bundle["attempt"] / "manifest.json").write_text(json.dumps(doc))
    with pytest.raises(score_mod.Refusal) as exc:
        _score(bundle)
    assert "runner_notes" in str(exc.value)


def test_the_scorer_refuses_a_manifest_that_names_the_arm(tmp_path):
    bundle = _attempt(tmp_path, "connected")
    doc = json.loads((bundle["attempt"] / "manifest.json").read_text())
    doc["arm"] = "specialist"
    (bundle["attempt"] / "manifest.json").write_text(json.dumps(doc))
    with pytest.raises(score_mod.Refusal) as exc:
        _score(bundle)
    assert "arm" in str(exc.value)


def test_the_scorer_refuses_an_unknown_task_key(tmp_path):
    bundle = _attempt(tmp_path, "connected")
    doc = json.loads((bundle["attempt"] / "manifest.json").read_text())
    doc["task"]["assigned_to"] = "x"
    (bundle["attempt"] / "manifest.json").write_text(json.dumps(doc))
    with pytest.raises(score_mod.Refusal):
        _score(bundle)


def test_files_beside_the_manifest_are_listed_and_not_read(tmp_path):
    bundle = _attempt(tmp_path, "connected")
    (bundle["attempt"] / "arm.txt").write_text("specialist\n")
    rec = _score(bundle)
    assert sorted(rec["ignored_files"]) == ["arm.txt", "transcript.txt"]
    assert "specialist" not in json.dumps(rec)


def test_the_log_is_projected_onto_the_cost_fields_at_the_read(tmp_path):
    bundle = _attempt(tmp_path, "connected")
    (bundle["attempt"] / "log.jsonl").write_text(json.dumps({
        "ts": "2026-11-20T10:00:00+00:00", "role": "assistant",
        "input_tokens": 10, "output_tokens": 2,
        "system_prompt": "you are the specialist for chemistry",
    }) + "\n")
    cost = score_mod.read_cost(bundle["attempt"])
    assert "specialist" not in json.dumps(cost)
    assert set(cost) == {"input_tokens", "output_tokens", "turns",
                         "wall_seconds", "log_records", "dropped_fields"}
    # The projection happened at the read: the field that carried the text
    # was discarded, and the count says so without naming it.
    assert cost["dropped_fields"] == 1


# --------------------------------------------------------------------------
# the pool generator
# --------------------------------------------------------------------------

def test_a_file_importing_a_chemistry_library_is_a_chemistry_file():
    sig = pools_mod.chemistry_signals("import rdkit\n\n\ndef f():\n    return 1\n")
    assert sig["s1"] is True
    assert sig["chemistry"] is True


def test_external_vocabulary_alone_can_make_a_file_a_chemistry_file():
    text = "ORCA " * 40 + "\n" + "def f():\n    return 1\n"
    sig = pools_mod.chemistry_signals(text)
    assert sig["s1"] is False
    assert sig["s2"] is True
    assert sig["hits"] >= pools_mod.S2_MIN_HITS


def test_a_dense_but_tiny_file_is_not_a_chemistry_file():
    """The absolute count beside the rate.

    A rate on its own promoted a 20-word banner module, which is the shape
    the measurement discipline rejects: a ratio with no count beside it.
    """
    sig = pools_mod.chemistry_signals("xyz dft homo lumo smiles\n")
    assert sig["per_1k"] > pools_mod.S2_MIN_PER_1K * 10
    assert sig["hits"] < pools_mod.S2_MIN_HITS
    assert sig["chemistry"] is False


def test_a_small_file_whose_identifiers_are_chemistry_is_a_chemistry_file():
    """S3: the clause for the module whose subject is chemistry and whose
    text is too short to reach S2's floor."""
    text = (
        "def read_orca_hess(orca_out):\n"
        "    xyz_path = orca_out + '.xyz'\n"
        "    hess_block = xyz_path\n"
        "    nmr_shieldings = hess_block\n"
        "    fukui_values = nmr_shieldings\n"
        "    return fukui_values, xyz_path, hess_block\n")
    sig = pools_mod.chemistry_signals(text)
    assert sig["s1"] is False
    assert sig["s2"] is False, sig
    assert sig["s3"] is True, sig
    assert sig["chemistry"] is True


def test_chemistry_words_in_a_string_literal_alone_are_not_enough():
    """The control for S3, and the file that made the rule read identifiers.

    A licence notice names ORCA, xTB and CREST over and over inside one
    string and has no chemistry identifier at all. The string carries more
    chemistry names than S3's fragment floor asks for, so a reading that
    took them for identifiers would promote the file -- that is the
    mutation this test has to fail on. S3 reads identifiers and refuses it;
    S2's absolute floor refuses the text.
    """
    text = ('BANNER = """\n' + "Install ORCA and xTB; CREST is optional. "
            * 5 + '\n"""\n\n\ndef show(out):\n    out.write(BANNER)\n')
    sig = pools_mod.chemistry_signals(text)
    # The text reading sees more than S3's floor asks for and still
    # refuses, because S2's own floor is higher.
    assert sig["hits"] >= pools_mod.S3_MIN_HITS, sig
    assert sig["hits"] < pools_mod.S2_MIN_HITS, sig
    assert sig["per_1k"] > pools_mod.S2_MIN_PER_1K, sig
    assert sig["s2"] is False
    # Not one of those names is an identifier of this file.
    assert sig["id_hits"] == 0, sig
    assert sig["s3"] is False
    assert sig["chemistry"] is False


def test_a_tiny_module_with_a_high_identifier_share_is_refused():
    """The share carries its own absolute floor, for the same reason the
    rate does: a delegating stub of 22 identifiers reached 18 %."""
    text = ("from delfin.orca import run\n\n\n"
            "def main(xyz):\n    return run(xyz)\n")
    sig = pools_mod.chemistry_signals(text)
    assert sig["id_share"] > pools_mod.S3_MIN_SHARE * 5
    assert sig["id_hits"] < pools_mod.S3_MIN_HITS
    assert sig["s3"] is False


_S3_BODY = (
    "    xyz_path = source + '.xyz'\n"
    "    hess_block = xyz_path\n"
    "    nmr_shieldings = hess_block\n"
    "    fukui_values = nmr_shieldings\n"
    "    return fukui_values, xyz_path, hess_block\n")


def test_a_plain_infrastructure_file_is_not_a_chemistry_file():
    sig = pools_mod.chemistry_signals(
        "import json\n\n\ndef read(path):\n    return json.loads(path)\n")
    assert sig["chemistry"] is False


def _pool_repo(tmp_path: Path) -> Path:
    """A tree with one infra orphan, one chem orphan and one of each
    excluded shape."""
    repo = tmp_path / "tree"
    _write(repo / "delfin" / "__init__.py", "")
    _write(repo / "delfin" / "agent" / "widget.py",
           "def infra_orphan():\n    return 1\n\n\n"
           "def ab(self=None):\n    return 2\n\n\n"
           "def lonely_orphan():\n    return 4\n\n\n"
           "def documented_orphan():\n    return 3\n")
    _write(repo / "delfin" / "agent" / "chem_in_agent.py",
           "import rdkit\n\n\ndef chem_orphan_in_agent():\n    return 1\n")
    _write(repo / "delfin" / "tools" / "thing.py",
           "def chem_orphan():\n    return 1\n")
    _write(repo / "delfin" / "agent" / "test_inline.py",
           "def test_only_defined_here():\n    assert True\n")
    _write(repo / "tests" / "test_widget.py",
           "from delfin.agent.widget import infra_orphan\n\n\n"
           "def test_infra_orphan():\n    assert infra_orphan() == 1\n")
    _write(repo / "docs" / "notes.md", "documented_orphan is described here\n")
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "tree")
    return repo


def test_the_pools_are_separate_and_the_exclusions_each_fire(tmp_path):
    repo = _pool_repo(tmp_path)
    doc = pools_mod.build_pools(repo, "HEAD")
    infra = {r["name"] for r in doc["pools"]["infra"]}
    chem = {r["name"] for r in doc["pools"]["chem"]}
    assert "infra_orphan" in infra
    assert "chem_orphan" in chem
    # A chemistry file inside delfin/agent/ leaves the infra pool and joins
    # the chemistry one, so the pools do not confound arm with domain.
    assert "chem_orphan_in_agent" not in infra
    assert "chem_orphan_in_agent" in chem
    removed = doc["exclusions_removed"]
    assert removed.get("short_name", 0) >= 1          # "ab"
    assert removed.get("defined_in_test", 0) >= 1     # delfin/agent/test_*
    assert "ab" not in infra
    assert "test_only_defined_here" not in infra
    # A name that only docs/ mentions stays a candidate: the scan never
    # reads docs/, so the mention is not reachability, and removing the
    # candidate would instead be removing it for being documented.
    assert "documented_orphan" in infra
    assert "documented_orphan" in doc["candidates_named_in_docs"]["infra"]
    assert "documented_only" not in removed


def test_the_document_records_the_rule_that_produced_the_pools(tmp_path):
    """A pools file read later has to say which rule drew it.

    The rule has already changed once -- S3 was added on 2026-11-20 and
    widened ``chem`` by 20 candidates -- so two pools files can differ in
    their pools and agree in every other field. The thresholds and the
    per-signal tally are therefore in the document, and the tally counts
    what S3 admitted *alone*, which is the number that says whether the
    clause earned its place.
    """
    repo = _pool_repo(tmp_path)
    # Chemistry by its identifiers only, text far below S2's floor.
    _write(repo / "delfin" / "small_chem.py",
           "def read_orca_hess(source):\n" + _S3_BODY)
    # Chemistry by import as well, so S3 and "S3 alone" cannot be equal.
    _write(repo / "delfin" / "both_chem.py",
           "import rdkit\n\n\ndef parse_orca_nmr(source):\n" + _S3_BODY)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "chemistry files")

    doc = pools_mod.build_pools(repo, "HEAD")
    rule = doc["chemistry_rule"]
    assert rule["s2_min_hits"] == pools_mod.S2_MIN_HITS
    assert rule["s2_min_per_1k"] == pools_mod.S2_MIN_PER_1K
    assert rule["s3_min_hits"] == pools_mod.S3_MIN_HITS
    assert rule["s3_min_share"] == pools_mod.S3_MIN_SHARE
    assert "fukui" in rule["observable_terms"]
    assert "ir" in rule["identifier_only_terms"]

    tally = rule["files_by_signal"]
    assert tally["scanned"] == len(list((repo / "delfin").rglob("*.py")))
    # chem_in_agent.py and both_chem.py import rdkit; small_chem.py and
    # both_chem.py reach S3; nothing in this tree reaches S2.
    assert tally["s1"] == 2
    assert tally["s2"] == 0
    assert tally["s3"] == 2
    assert tally["s3_only"] == 1
    assert tally["chemistry"] == 3
    assert rule["chemistry_files_under_delfin"] == 3


def test_a_name_whose_module_is_absent_from_the_comparison_ref_is_removed(
        tmp_path):
    repo = _pool_repo(tmp_path)
    ref = _git(repo, "rev-parse", "HEAD").strip()
    _write(repo / "delfin" / "agent" / "arrived_later.py",
           "def orphan_that_arrived_later():\n    return 1\n")
    doc = pools_mod.build_pools(repo, ref)
    assert doc["exclusions_removed"].get("module_absent", 0) >= 1
    names = {r["name"] for r in doc["pools"]["infra"]}
    assert "orphan_that_arrived_later" not in names


def test_two_unreachable_twins_hide_each_other_and_that_is_reported(tmp_path):
    """The limitation of a name-based scan, measured rather than assumed.

    A second module-level definition of the same name counts as a mention
    of the first, so neither twin is ever a candidate. That is why there is
    no exclusion for a name in both pools: it cannot occur. The field is
    checked because the day the scan changes, the assumption stops holding
    silently.
    """
    repo = _pool_repo(tmp_path)
    _write(repo / "delfin" / "agent" / "dual.py",
           "def defined_twice_orphan():\n    return 1\n")
    _write(repo / "delfin" / "tools" / "dual.py",
           "def defined_twice_orphan():\n    return 2\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "dual")
    assert "defined_twice_orphan" in pools_mod._twin_definitions(repo)
    doc = pools_mod.build_pools(repo, "HEAD")
    for pool in ("infra", "chem"):
        assert "defined_twice_orphan" not in {
            r["name"] for r in doc["pools"][pool]}
        assert doc["candidates_with_a_twin"][pool] == []


def test_the_scan_does_not_read_the_baseline_that_lists_the_names(tmp_path):
    """The instrument measuring itself.

    The baseline file names every capability the scan looks for and lives
    under ``tests/``. Reading it makes each name read as reached by the
    suite, which empties the sharp tier: ``lonely_orphan`` is named by no
    test, so its tier must stay ``nowhere`` although the baseline spells it.
    """
    repo = _pool_repo(tmp_path)
    _write(repo / "tests" / "fixtures" / "unmounted_capabilities.txt",
           "nowhere\tlonely_orphan\tdelfin/agent/widget.py:7\n")
    doc = pools_mod.build_pools(repo, "HEAD")
    rows = {r["name"]: r["tier"] for r in doc["pools"]["infra"]}
    assert rows.get("lonely_orphan") == "nowhere", rows


def test_an_underpowered_pool_is_reported_with_its_number(tmp_path):
    repo = _pool_repo(tmp_path)
    doc = pools_mod.build_pools(repo, "HEAD")
    assert doc["underpowered"]["infra"] is True
    assert doc["counts"]["infra"] < 10


def test_the_copy_names_the_blob_it_was_copied_from():
    """A copy that does not say what it copied cannot be checked for drift.

    The original is not on ``origin/main``, so there is nothing to compare
    against on the tree the experiment runs on; the blob id is what makes
    the comparison possible at all.
    """
    assert "f290b9753151c5c9b323865bf37b96e1c09bcf66" in pools_mod.__doc__


def test_the_copied_scan_matches_the_original_where_both_are_present():
    """Pins the copy against the test it was copied from.

    Where a tree holds both, the two must answer the same, or the copy has
    drifted. No skip: on a tree without the original there is nothing to
    compare, and the test above is what holds in that case.
    """
    original = ROOT / "tests" / "test_every_capability_has_a_caller.py"
    if original.is_file():
        import importlib.util
        spec = importlib.util.spec_from_file_location("orig_scan", original)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert mod._definitions() == pools_mod._definitions(
            ROOT / "delfin" / "agent", ROOT)


# --------------------------------------------------------------------------
# the draw
# --------------------------------------------------------------------------

def _fake_pools(infra: int, chem: int, nowhere_share: float = 0.5) -> dict:
    pools = {}
    for pool, n in (("infra", infra), ("chem", chem)):
        rows = []
        for i in range(n):
            tier = "nowhere" if i < int(n * nowhere_share) else "test-only"
            rows.append({"name": f"{pool}_cap_{i:03d}", "tier": tier,
                         "file": f"delfin/{pool}/m{i}.py", "line": i + 1,
                         "module": f"delfin.{pool}.m{i}",
                         "package": f"delfin/{pool}"})
        pools[pool] = rows
    from collections import Counter
    return {
        "generated_for": "/nowhere",
        "pools": pools,
        "counts": {p: len(pools[p]) for p in pools},
        "counts_by_tier": {p: dict(Counter(r["tier"] for r in pools[p]))
                           for p in pools},
        "underpowered": {p: len(pools[p]) < 10 for p in pools},
        "task_sentence": "connect this capability to the product, or delete "
                         "it together with its tests",
    }


def test_the_draw_is_equal_per_pool_and_equal_per_tier():
    """The balanced rule, which is what the protocol uses.

    The pool's own tier mix is 25 % / 75 % and the draw is 50 / 50, so tier
    varies inside each arm instead of between the pools.
    """
    doc = _fake_pools(40, 40, nowhere_share=0.25)
    sample = draw_mod.draw(doc, seed=7, per_pool=12)
    assert sample["n_tasks"] == 24
    assert sample["allocation_rule"] == "balanced"
    for pool in ("infra", "chem"):
        got = sample["drawn_counts"][pool]
        assert sum(got.values()) == 12
        assert got["nowhere"] == 6
        assert got["test-only"] == 6


def test_the_superseded_proportional_rule_is_still_reproducible():
    """The first draw is superseded by an addendum, not deleted, so the
    rule that produced it stays available and is named in the output."""
    doc = _fake_pools(40, 40, nowhere_share=0.25)
    sample = draw_mod.draw(doc, seed=7, per_pool=12, rule="proportional")
    assert sample["allocation_rule"] == "proportional"
    for pool in ("infra", "chem"):
        assert sample["drawn_counts"][pool]["nowhere"] == 3
        assert sample["drawn_counts"][pool]["test-only"] == 9


def test_a_balanced_draw_refuses_a_pool_whose_tier_is_too_small():
    doc = _fake_pools(40, 40, nowhere_share=0.05)      # 2 of 40 are nowhere
    with pytest.raises(ValueError) as exc:
        draw_mod.draw(doc, seed=7, per_pool=12)
    assert "balanced" in str(exc.value)
    assert "nowhere" in str(exc.value)


def test_a_balanced_draw_of_three_strata_uses_the_floor_not_the_exact_share():
    """Correction F: refuse a split that cannot be made, not one that can.

    With three strata the exact share is not an integer, and comparing a
    stratum's stock against it refused allocations that exist: {3, 4, 5}
    for a total of 10 is feasible as 4/3/3 and was rejected because
    3 < 3.33. Harmless at the two tiers and 10 per pool this experiment
    draws -- there the floor and the exact share agree -- and still a
    defect of the instrument.
    """
    got = draw_mod.allocate({"a": 3, "b": 4, "c": 5}, 10)
    assert sum(got.values()) == 10
    assert all(got[k] <= avail
               for k, avail in {"a": 3, "b": 4, "c": 5}.items())
    assert min(got.values()) >= 10 // 3
    assert max(got.values()) - min(got.values()) <= 1
    # An odd total over two strata is the same case.
    odd = draw_mod.allocate({"nowhere": 5, "test-only": 40}, 11)
    assert sum(odd.values()) == 11 and odd["nowhere"] == 5

    # What the refusal is for: a stratum that cannot supply the floor.
    with pytest.raises(ValueError) as exc:
        draw_mod.allocate({"a": 2, "b": 4, "c": 5}, 10)
    assert "balanced" in str(exc.value) and "'a': 2" in str(exc.value)


def test_the_same_seed_draws_the_same_tasks_and_another_seed_does_not():
    doc = _fake_pools(40, 40)
    a = draw_mod.draw(doc, seed=11, per_pool=10)["tasks"]
    b = draw_mod.draw(doc, seed=11, per_pool=10)["tasks"]
    c = draw_mod.draw(doc, seed=12, per_pool=10)["tasks"]
    assert a == b
    assert a != c


def test_a_short_pool_refuses_the_draw_rather_than_drawing_fewer():
    doc = _fake_pools(40, 6)
    with pytest.raises(ValueError) as exc:
        draw_mod.draw(doc, seed=3, per_pool=10)
    assert "6" in str(exc.value)
    # The pool-level refusal, not the allocator's: the message has to say
    # why an unequal draw is refused rather than only that it does not fit.
    assert "confound" in str(exc.value)


def test_the_sample_file_refuses_to_be_redrawn(tmp_path):
    pools = tmp_path / "pools.json"
    pools.write_text(json.dumps(_fake_pools(40, 40)))
    out = tmp_path / "sample.json"
    args = ["--pools", str(pools), "--out", str(out), "--seed", "5",
            "--per-pool", "10"]
    assert draw_mod.main(args) == 0
    first = out.read_text()
    assert draw_mod.main(["--seed", "6", *args[:4]]) == 2
    assert out.read_text() == first
    assert draw_mod.main(["--seed", "6", *args[:4], "--force"]) == 0
    replaced = json.loads(out.read_text())
    assert replaced["seed"] == 6
    assert replaced["replaced"]["sha256"] == hashlib.sha256(
        first.encode()).hexdigest()


def test_the_sample_records_the_seed_the_date_and_the_tree():
    doc = _fake_pools(40, 40)
    sample = draw_mod.draw(doc, seed=99, per_pool=10)
    assert sample["seed"] == 99
    assert sample["drawn_at"].startswith("20")
    assert sample["pool_counts"] == {"infra": 40, "chem": 40}
    assert "pools_tree_sha" in sample
