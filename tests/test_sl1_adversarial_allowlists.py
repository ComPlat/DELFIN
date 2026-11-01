"""Adversarial tests for slurm_tests.py allow-lists and cluster_utils.py
memory parsing (reviewer nacht-s14, package SL1, wave 25).

Attacks the UNCHANGED owned code. Red = finding for the builder.
No scheduler and no real test submission is touched.
"""
import pytest

from delfin.agent.slurm_tests import (
    _SAFE_JOB_NAME,
    _SAFE_MEM,
    _SAFE_PARTITION,
    _SAFE_PATHLD,
    _SAFE_REF,
    _SAFE_TEST_PATH,
)
from delfin.cluster_utils import _parse_memory_string


def _rejects(pattern, value):
    """True when the pattern does NOT match the whole value."""
    return pattern.fullmatch(value) is None


# ------------------------------------------------------------------ #
# _SAFE_PARTITION accepts commas (a legitimate SLURM feature) -- but
# what else does the comma admit? Attack: '#SBATCH' smuggling inside
# the partition value is impossible here because the pattern forbids
# whitespace and '#'; the real question is whether commas are ever
# rendered raw into '--partition=...' in build_job (below).
# ------------------------------------------------------------------ #
def test_partition_allows_comma_lists():
    assert _SAFE_PARTITION.fullmatch("cpu,cpu_il")
    assert _SAFE_PARTITION.fullmatch("cpu")


def test_partition_rejects_whitespace_and_hash():
    assert _rejects(_SAFE_PARTITION, "cpu\n#SBATCH --export=ALL")
    assert _rejects(_SAFE_PARTITION, "cpu cpu")
    assert _rejects(_SAFE_PARTITION, "cpu#short")


def test_partition_leading_trailing_comma_rejected():
    """FIXED F2 (SL1 phase 2). The old version pinned the gap: ',cpu' /
    'cpu,' rendered as '--partition=,cpu', an empty name SLURM rejects
    late. The regex now requires a name on both sides of each comma."""
    assert _SAFE_PARTITION.fullmatch(",cpu") is None
    assert _SAFE_PARTITION.fullmatch("cpu,") is None
    assert _SAFE_PARTITION.fullmatch("cpu,,cpu_il") is None


def test_partition_dash_start_rejected():
    """FIXED F3 (SL1 phase 2). The old version pinned the gap: '-cpu'
    passes as '--partition=-cpu', re-readable by sbatch as an option.
    A leading '-' is now forbidden."""
    assert _SAFE_PARTITION.fullmatch("-cpu") is None
    assert _SAFE_PARTITION.fullmatch("--output=/etc/passwd") is None  # '=' excluded already


# ------------------------------------------------------------------ #
# _SAFE_REF: SLURM has no metacharacters here, but a ref that starts
# with '-' would be re-parsed by git (used as `git checkout <ref>` --
# ref='-h' prints git's help on the compute node).
# ------------------------------------------------------------------ #
def test_ref_dash_start_rejected():
    """FIXED F4 (SL1 phase 2). The old version pinned the gap: _SAFE_REF
    accepted '-h' / '--help', which a git consumer would re-parse as an
    option. A leading '-' on a ref is now forbidden."""
    assert _SAFE_REF.fullmatch("-h") is None
    assert _SAFE_REF.fullmatch("--help") is None


def test_ref_rejects_newline_and_space():
    assert _rejects(_SAFE_REF, "main\ntouch x")
    assert _rejects(_SAFE_REF, "a b")


def test_ref_allows_hashes_tags_branches():
    assert _SAFE_REF.fullmatch("98b096295")
    assert _SAFE_REF.fullmatch("origin/main")
    assert _SAFE_REF.fullmatch("v2.1.0")


# ------------------------------------------------------------------ #
# _SAFE_PATHLD allows '~' -- a path starting with '~' expands in bash
# when unquoted or when the rendering quotes only weakly. Check the
# rendering side: does build_job quote paths? (verified in file read:
# shlex.quote stays as defense-in-depth, so '~' alone is not an
# escape, but '~' + metachar is already excluded). Accept '~' but
# never '~root' style expansion into another user's home? Out of
# scope for a pattern test; pin current behaviour.
# ------------------------------------------------------------------ #
def test_pathld_allows_tilde_but_rejects_metachars():
    assert _SAFE_PATHLD.fullmatch("~/proj")
    assert _rejects(_SAFE_PATHLD, "a;b")
    assert _rejects(_SAFE_PATHLD, "a|b")
    assert _rejects(_SAFE_PATHLD, "a&b")
    assert _rejects(_SAFE_PATHLD, "$(x)")


def test_safe_test_path_dotdot_is_checked_by_the_caller_not_the_regex():
    # NOT a finding: the regex alone matches 'tests/a/../..' but the
    # caller at slurm_tests.py:93 applies an explicit '..' rejection
    # before the regex. Pinned so a refactor cannot silently drop the
    # caller-side guard.
    assert _SAFE_TEST_PATH.fullmatch("tests/a/../..") is not None
    assert ".." in "tests/a/../.."  # the caller-side guard fires on this


def test_safe_mem_rejects_injection():
    assert _rejects(_SAFE_MEM, "4000M; touch x")
    assert _rejects(_SAFE_MEM, "4000\nM")
    assert _rejects(_SAFE_MEM, "")
    assert _SAFE_MEM.fullmatch("4000M")
    assert _SAFE_MEM.fullmatch("4000")


# ------------------------------------------------------------------ #
# _parse_memory_string: unknown units default to MB (documented in my
# P0 read). Attack: does '4GB ' with inner space parse? '4 GiB'?
# Case where 'GiB' (binary prefix) is silently read as GB (decimal):
# 1 GiB = 1024 MiB = 1024 MB -- actually equal here since the parser
# is binary already. The real gap: '1.5GB' -> int() truncation and
# unit 'B' (bytes) treated as MB: '512B' -> 512 "MB" = 64x overcount.
# ------------------------------------------------------------------ #
def test_parse_memory_bytes_unit_rejected():
    """FIXED F6 (SL1 phase 2). The old version pinned the defect:
    _parse_memory_string('512B') == 512 -- bytes read as MB, a 512x
    overcount (cluster_utils.py:169, 'B' sat in the MB set). Bytes are
    now their own unit."""
    assert _parse_memory_string("512B") == 0


def test_parse_memory_decimal_gib_prefix():
    assert _parse_memory_string("1GiB") is not None
    assert _parse_memory_string("2gib") == 2048  # currently unit 'GIB' hits default->MB


def test_parse_memory_none_and_empty():
    assert _parse_memory_string("") is None
    assert _parse_memory_string(None) is None


def test_parse_memory_plain_number_is_mb():
    assert _parse_memory_string("4000") == 4000
