#!/usr/bin/env bash
# DELFIN pre-commit hook: run the refactor-safety test suite before any commit.
#
# Installation:
#     cp scripts/pre-commit-hook.sh .git/hooks/pre-commit
#     chmod +x .git/hooks/pre-commit
#
# What it does:
#     Runs the subset of tests that guarantee "DELFIN's existing
#     functionality is preserved" — the ones that pointed at every
#     NameError and API regression in the recent bug wave. This keeps
#     refactor-breaks from landing on main.
#
# What it doesn't do:
#     Run the full test suite — that takes ~10s, too long for every
#     commit. Run `pytest tests/` yourself before pushing.
#
# Opt-out for emergencies:
#     git commit --no-verify  (but please don't — these exist because
#     silent regressions cost Tilman and Jan many lost CPU-hours)

set -e

echo "[pre-commit] Running DELFIN refactor-safety tests…"

TESTS=(
    tests/test_no_regression_undefined_names.py
    tests/test_functional_contracts.py
    tests/test_pre_refactor_api_preserved.py
    tests/test_submit_wrapper.py
    tests/test_classic_inline_recording.py
    tests/test_shim_integrity.py
)

if ! python3 -m pytest -q "${TESTS[@]}"; then
    echo ""
    echo "[pre-commit] ABORT: refactor-safety tests failed."
    echo ""
    echo "These tests exist because refactors kept silently breaking"
    echo "user-visible behaviour (NameErrors, missing imports, changed"
    echo "signatures, lost APIs). If the failure is INTENTIONAL — you"
    echo "deliberately renamed or removed something users rely on —"
    echo "add a waiver to _WAIVED_DROPS in"
    echo "tests/test_pre_refactor_api_preserved.py explaining why, then"
    echo "re-commit."
    echo ""
    echo "Bypass only in emergencies with: git commit --no-verify"
    exit 1
fi

echo "[pre-commit] All refactor-safety tests passed."
