"""Remove each check of the apparatus in turn and confirm a test fails.

Input: ``--root`` (a checkout holding ``scripts/ab_*.py`` and
``tests/test_ab_apparatus.py``) and ``--python``. Output: one line per
mutation saying whether the named tests caught it, then the restored
suite's verdict; exit 1 when any mutation went unnoticed or the restored
suite is not green.

Semantics: a test that passes whether or not the check exists is attached
to something other than the check it names. Four of the tests here passed
their first mutation and were rewritten until they did not, which is the
only reason to run this rather than to assert it.

Each mutation is applied to a working copy of the file, the pytest
selection is run, and the file is restored in a ``finally`` block from the
text read before the edit; the restore is verified by re-reading. Nothing
is committed and no other path is touched, so an interrupted run leaves the
file as it was. It edits the tree it is pointed at, so it is not run
against a checkout that another process is using.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable

MUTATIONS = [
    ("score: check 1 capability_resolved", "scripts/ab_score.py",
     '"capability_resolved": route["route"] in ("connected", "deleted"),',
     '"capability_resolved": True,',
     ["test_an_untouched_tree_scores_zero_without_being_dropped"]),
    ("score: check 2 tests_still_pass", "scripts/ab_score.py",
     '"tests_still_pass": not (tests["regressed"] or',
     '"tests_still_pass": True or (tests["regressed"] or',
     ["test_a_test_that_turned_into_a_skip_is_not_a_pass"]),
    ("score: check 2, the skip clause only", "scripts/ab_score.py",
     'tests["turned_into_skip"] or\n',
     '',
     ["test_a_test_that_turned_into_a_skip_is_not_a_pass"]),
    ("score: check 3 instrument_untouched", "scripts/ab_score.py",
     '"instrument_untouched": instrument["ok"],',
     '"instrument_untouched": True,',
     ["test_touching_the_instrument_scores_zero_and_records_why",
      "test_editing_the_scoring_scripts_scores_zero"]),
    ("score: check 3, the scripts/ab_*.py glob", "scripts/ab_score.py",
     'INSTRUMENT_GLOBS = ("scripts/ab_*.py",)',
     'INSTRUMENT_GLOBS = ()',
     ["test_editing_the_scoring_scripts_scores_zero"]),
    ("score: check 3, added lines in the fixture", "scripts/ab_score.py",
     'allowed = (not added and not removed_prose',
     'allowed = (True and not removed_prose',
     ["test_a_line_added_to_the_baseline_is_not_an_allowed_fixture_edit"]),
    ("score: check 4 deletion_removed_its_tests", "scripts/ab_score.py",
     '"deletion_removed_its_tests": deletion["ok"],',
     '"deletion_removed_its_tests": True,',
     ["test_a_deletion_that_left_its_tests_behind_scores_zero"]),
    ("score: the manifest key allowlist", "scripts/ab_score.py",
     'unknown = sorted(set(doc) - MANIFEST_KEYS)',
     'unknown = []',
     ["test_the_scorer_refuses_an_unknown_manifest_key"]),
    ("score: the task key allowlist", "scripts/ab_score.py",
     'unknown = sorted(set(task) - TASK_KEYS)',
     'unknown = []',
     ["test_the_scorer_refuses_an_unknown_task_key"]),
    ("score: the call site as the primary endpoint", "scripts/ab_score.py",
     'primary = "connected" if calls else "still_unreachable"',
     'primary = "connected" if mentions else "still_unreachable"',
     ["test_a_string_literal_satisfies_the_secondary_endpoint_and_not_the_score"]),
    ("score: the call-site scan itself", "scripts/ab_score.py",
     '            if hit:\n                out.append(f"{path.relative_to(root)}:{node.lineno}")',
     '            if False:\n                out.append(f"{path.relative_to(root)}:{node.lineno}")',
     ["test_a_real_call_satisfies_both_endpoints_and_they_agree"]),
    ("score: the secondary endpoint is reported too", "scripts/ab_score.py",
     'secondary = "connected" if mentions else "still_unreachable"',
     'secondary = "connected" if calls else "still_unreachable"',
     ["test_a_string_literal_satisfies_the_secondary_endpoint_and_not_the_score"]),
    ("score: the log projection", "scripts/ab_score.py",
     'projected = {k: merged[k] for k in LOG_KEYS if k in merged}',
     'projected = dict(merged)',
     ["test_the_log_is_projected_onto_the_cost_fields_at_the_read"]),
    ("score: a missing figure is null", "scripts/ab_score.py",
     'out[key] = sum(vals) if vals else None',
     'out[key] = sum(vals)',
     ["test_a_missing_figure_is_null_and_never_zero"]),
    ("pools: the baseline is not scanned", "scripts/ab_task_pools.py",
     '        if path in skip_paths:\n            continue\n        try:\n'
     '            text = path.read_text(errors="replace")\n'
     '        except OSError:\n            continue\n'
     '        rel = str(path.relative_to(root))',
     '        try:\n'
     '            text = path.read_text(errors="replace")\n'
     '        except OSError:\n            continue\n'
     '        rel = str(path.relative_to(root))',
     ["test_the_scan_does_not_read_the_baseline_that_lists_the_names"]),
    ("pools: the absolute count beside the rate", "scripts/ab_task_pools.py",
     's2 = hits >= S2_MIN_HITS and per_1k >= S2_MIN_PER_1K',
     's2 = per_1k >= S2_MIN_PER_1K',
     ["test_a_dense_but_tiny_file_is_not_a_chemistry_file"]),
    ("pools: the S3 identifier clause", "scripts/ab_task_pools.py",
     '    s3 = id_hits >= S3_MIN_HITS and id_share >= S3_MIN_SHARE',
     '    s3 = False',
     ["test_a_small_file_whose_identifiers_are_chemistry_is_a_chemistry_file"]),
    ("pools: S3 reads identifiers, not text", "scripts/ab_task_pools.py",
     '        if tok.type != tokenize.NAME or keyword.iskeyword(tok.string):\n'
     '            continue',
     '        if False:\n            continue',
     ["test_chemistry_words_in_a_string_literal_alone_are_not_enough"]),
    ("pools: the S3 absolute floor", "scripts/ab_task_pools.py",
     's3 = id_hits >= S3_MIN_HITS and id_share >= S3_MIN_SHARE',
     's3 = id_share >= S3_MIN_SHARE',
     ["test_a_tiny_module_with_a_high_identifier_share_is_refused"]),
    ("pools: the S3 thresholds the document records",
     "scripts/ab_task_pools.py",
     '            "s3_min_hits": S3_MIN_HITS,',
     '            "s3_min_hits": 0,',
     ["test_the_document_records_the_rule_that_produced_the_pools"]),
    ("pools: the tally of what S3 admitted alone",
     "scripts/ab_task_pools.py",
     '                       if s["s3"] and not s["s1"] and not s["s2"]),',
     '                       if s["s3"]),',
     ["test_the_document_records_the_rule_that_produced_the_pools"]),
    ("pools: exclusion documented_only", "scripts/ab_task_pools.py",
     'if row["tier"] == "nowhere" and name in prose_names:',
     'if False:',
     ["test_the_pools_are_separate_and_the_exclusions_each_fire"]),
    ("pools: exclusion defined_in_test", "scripts/ab_task_pools.py",
     'elif _is_test_file(row["file"]):',
     'elif False:',
     ["test_the_pools_are_separate_and_the_exclusions_each_fire"]),
    ("pools: exclusion short_name", "scripts/ab_task_pools.py",
     'elif len(name) < 4:',
     'elif False:',
     ["test_the_pools_are_separate_and_the_exclusions_each_fire"]),
    ("pools: exclusion module_absent", "scripts/ab_task_pools.py",
     'elif tracked is not None and row["file"] not in tracked:',
     'elif False:',
     ["test_a_name_whose_module_is_absent_from_the_comparison_ref_is_removed"]),
    ("pools: the chemistry split of the two pools",
     "scripts/ab_task_pools.py",
     '    def infra(path: Path) -> bool:\n'
     '        return not verdict(path)',
     '    def infra(path: Path) -> bool:\n'
     '        return True',
     ["test_the_pools_are_separate_and_the_exclusions_each_fire"]),
    ("pools: the underpowered report", "scripts/ab_task_pools.py",
     '"underpowered": {p: len(pools[p]) < 10 for p in pools},',
     '"underpowered": {p: False for p in pools},',
     ["test_an_underpowered_pool_is_reported_with_its_number"]),
    ("pools: the twin report", "scripts/ab_task_pools.py",
     'return {n for n, c in seen.items() if c > 1}',
     'return set()',
     ["test_two_unreachable_twins_hide_each_other_and_that_is_reported"]),
    ("draw: the balanced tier rule", "scripts/ab_draw_sample.py",
     '        exact = {k: total / len(counts) for k in counts}',
     '        exact = {k: counts[k] / pool_size * total for k in counts}',
     ["test_the_draw_is_equal_per_pool_and_equal_per_tier"]),
    ("draw: the proportional rule stays available",
     "scripts/ab_draw_sample.py",
     '    elif rule == "proportional":\n'
     '        exact = {k: counts[k] / pool_size * total for k in counts}',
     '    elif rule == "proportional":\n'
     '        exact = {k: total / len(counts) for k in counts}',
     ["test_the_superseded_proportional_rule_is_still_reproducible"]),
    ("draw: the balanced short-tier refusal", "scripts/ab_draw_sample.py",
     '        if short:',
     '        if False:',
     ["test_a_balanced_draw_refuses_a_pool_whose_tier_is_too_small"]),
    ("draw: the allocation rule is recorded", "scripts/ab_draw_sample.py",
     '        "allocation_rule": rule,',
     '        "allocation_rule": "balanced",',
     ["test_the_superseded_proportional_rule_is_still_reproducible"]),
    ("draw: the short-pool refusal", "scripts/ab_draw_sample.py",
     '        if len(rows) < per_pool:',
     '        if False:',
     ["test_a_short_pool_refuses_the_draw_rather_than_drawing_fewer"]),
    ("draw: the re-draw refusal", "scripts/ab_draw_sample.py",
     '        if not args.force:',
     '        if False:',
     ["test_the_sample_file_refuses_to_be_redrawn"]),
    ("draw: the seed", "scripts/ab_draw_sample.py",
     '    rng = random.Random(seed)',
     '    rng = random.Random(0)',
     ["test_the_same_seed_draws_the_same_tasks_and_another_seed_does_not"]),
]


def run(tests: list[str]) -> tuple[bool, str]:
    """Run the pytest selection in ROOT. Output: (green, the last line)."""
    expr = " or ".join(tests)
    out = subprocess.run(
        [PY, "-m", "pytest", "tests/test_ab_apparatus.py", "-q",
         "--no-header", "-p", "no:cacheprovider", "-k", expr],
        cwd=str(ROOT), capture_output=True, text=True, timeout=1800)
    tail = out.stdout.strip().splitlines()[-1] if out.stdout.strip() else ""
    return out.returncode == 0, tail


def main(argv=None) -> int:
    """Apply every mutation, report, restore. Output: the exit status."""
    global ROOT, PY
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--python", default=PY)
    args = ap.parse_args(argv)
    ROOT, PY = Path(args.root), args.python
    bad = []
    for label, rel, old, new, tests in MUTATIONS:
        path = ROOT / rel
        before = path.read_text()
        if old not in before:
            print(f"SETUP FAIL  {label}: pattern not found")
            bad.append(label)
            continue
        if before.count(old) != 1:
            print(f"SETUP FAIL  {label}: pattern occurs "
                  f"{before.count(old)} times")
            bad.append(label)
            continue
        try:
            path.write_text(before.replace(old, new))
            green, tail = run(tests)
        finally:
            path.write_text(before)
            assert path.read_text() == before, f"restore failed for {rel}"
        verdict = "NOT CAUGHT" if green else "caught"
        if green:
            bad.append(label)
        print(f"{verdict:11} {label}  ->  {tail}")
    print()
    green, tail = run(["test_"])
    print(f"restored suite: {'green' if green else 'RED'}  {tail}")
    if bad:
        print(f"\nmutations not caught: {bad}")
    return 1 if bad or not green else 0


if __name__ == "__main__":
    sys.exit(main())
