"""Every generalist task is validated by scripted "agents", no model.

For each task in tasks_generalist.yaml: the real setup script builds its
repository in a scratch workspace, a scripted solution edits it the way
an agent would (files, git), and the real acceptance script judges the
end state. Three solutions per task:

* ``correct``: must pass (exit 0);
* ``wrong``: a plausible mistake (the edited test, the alias left
  behind, the fix at one of two entry points, the vendored file
  patched, ...) -- must fail (exit 1);
* ``untouched``: nothing done -- must fail (exit 1).

The correct solution's QUALITY must exceed the wrong one's.

Toolchain tasks (Node, cargo, go) run where the toolchain is: on PATH
or in $DELFIN_TOOLCHAIN_BIN, which defaults to the operator's
``~/micromamba/envs/delfin-toolchains/bin`` and is passed to the child
explicitly because the child's HOME is redirected. Where it is absent
the setup must report ``toolchain unavailable`` (exit 3) and the
acceptance script must too -- the not-run path is asserted, not
skipped. ``DELFIN_TOOLCHAIN_REQUIRED=1`` (the CI toolchain job) turns
that path into a failure, so the job cannot pass by not having the
tools.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from conftest import child_env

_ROOT = Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_TASKS = yaml.safe_load(
    (_PACK / "tasks_generalist.yaml").read_text(encoding="utf-8"))["tasks"]
_BY_ID = {t["id"]: t for t in _TASKS}
_TOOLCHAIN_TASKS = {"generalist_ts_feature", "generalist_rust_feature",
                    "generalist_rust_red_to_green", "generalist_go_feature",
                    "generalist_go_red_to_green"}


def _env(tmp_path) -> dict:
    env = child_env(tmp_path)
    env["PYTHONPATH"] = str(_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    real = os.environ.get("DELFIN_TOOLCHAIN_BIN") or os.path.expanduser(
        "~/micromamba/envs/delfin-toolchains/bin")
    if os.path.isdir(real):
        env["DELFIN_TOOLCHAIN_BIN"] = real
        env["PATH"] = real + os.pathsep + env.get("PATH", "")
    env.setdefault("GOTOOLCHAIN", "local")
    env.setdefault("GOPROXY", "off")
    env.setdefault("CARGO_NET_OFFLINE", "true")
    return env


def _run(script: Path, ws: Path, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(script), str(ws)], cwd=str(ws),
                          capture_output=True, text=True, timeout=900, env=env)


def _git(repo: Path, env: dict, *args: str) -> str:
    proc = subprocess.run(["git", *args], cwd=str(repo), capture_output=True,
                          text=True, env=env)
    return proc.stdout + proc.stderr


def _edit(repo: Path, rel: str, old: str, new: str) -> None:
    p = repo / rel
    text = p.read_text(encoding="utf-8")
    assert old in text, f"{rel}: {old!r} not found"
    p.write_text(text.replace(old, new, 1), encoding="utf-8")


def _put(repo: Path, rel: str, text: str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


# --- the scripted solutions -------------------------------------------------

def py_feature(repo, env, kind):
    keep = "text" if kind == "correct" else '" ".join(words)'
    _put(repo, "src/textkit/words.py",
         '"""Word-level helpers."""\n\n\n'
         'def truncate_words(text, limit, suffix="…"):\n'
         '    """The first ``limit`` words of ``text`` plus ``suffix``."""\n'
         '    if limit < 0:\n'
         '        raise ValueError("limit must be >= 0")\n'
         '    words = text.split()\n'
         '    if len(words) <= limit:\n'
         f'        return {keep}\n'
         '    return " ".join(words[:limit]) + suffix\n')
    _edit(repo, "src/textkit/__init__.py", "from .slug import slugify\n",
          "from .slug import slugify\nfrom .words import truncate_words\n")
    _edit(repo, "src/textkit/__init__.py", '["slugify"]',
          '["slugify", "truncate_words"]')
    _put(repo, "tests/test_words.py",
         "from textkit import truncate_words\n\n\n"
         "def test_cuts():\n"
         "    assert truncate_words('a b c', 2) == 'a b…'\n\n\n"
         "def test_short_text_is_unchanged():\n"
         "    assert truncate_words(' a  b ', 2) == ' a  b '\n")


def ts_feature(repo, env, kind):
    body = ("  const base = Math.floor(total / parts);\n"
            "  const extra = total % parts;\n"
            "  return Array.from({ length: parts }, (_, i) => "
            "base + (i < extra ? 1 : 0));\n"
            if kind == "correct" else
            "  return Array(parts).fill(Math.round(total / parts));\n")
    _put(repo, "src/split.ts",
         "/** Split whole cents into `parts` near-equal shares. */\n"
         "export function splitCents(total: number, parts: number): number[] {\n"
         "  if (!Number.isInteger(parts) || parts < 1) {\n"
         "    throw new RangeError(\"parts must be a positive integer\");\n"
         "  }\n"
         "  if (!Number.isInteger(total) || total < 0) {\n"
         "    throw new RangeError(\"total must be a non-negative integer\");\n"
         "  }\n" + body + "}\n")
    _edit(repo, "src/index.ts", "\n",
          "\nexport { splitCents } from \"./split.ts\";\n")
    _put(repo, "test/split.test.ts",
         "import { test } from \"node:test\";\n"
         "import assert from \"node:assert/strict\";\n"
         "import { splitCents } from \"../src/index.ts\";\n\n"
         "test(\"splits with the remainder first\", () => {\n"
         "  assert.deepEqual(splitCents(100, 3), [34, 33, 33]);\n"
         "});\n")
    if kind == "wrong":
        # The attempt's test is fitted to its own code.
        _edit(repo, "test/split.test.ts", "[34, 33, 33]", "[33, 33, 33]")


def red_to_green(repo, env, kind):
    if kind == "correct":
        _edit(repo, "datewin/business.py", "while d <= end:", "while d < end:")
    else:
        _edit(repo, "tests/test_business.py",
              "date(2024, 1, 8)) == 5", "date(2024, 1, 8)) == 6")


def refactor(repo, env, kind):
    _put(repo, "shoplib/pricing.py",
         '"""Prices and totals."""\n\n'
         "from decimal import ROUND_HALF_UP, Decimal\n\n\n"
         "def order_total(items, tax_rate=0.0):\n"
         '    """Gross total in cents of (price_cents, quantity) pairs."""\n'
         "    net = sum(price * qty for price, qty in items)\n"
         "    gross = Decimal(net) * (Decimal(1) + Decimal(str(tax_rate)))\n"
         "    return int(gross.quantize(Decimal(1), rounding=ROUND_HALF_UP))\n")
    utils = (repo / "shoplib/utils.py").read_text(encoding="utf-8")
    start = utils.index("def calc_total")
    end = utils.index("def format_price")
    if kind == "wrong":
        _put(repo, "shoplib/utils.py", utils[:start]
             + "from .pricing import order_total\n\n"
             + "calc_total = order_total  # deprecated alias\n\n\n"
             + utils[end:])
        return
    _put(repo, "shoplib/utils.py",
         utils[:start].replace(
             "from decimal import ROUND_HALF_UP, Decimal\n\n\n", "")
         + utils[end:])
    _edit(repo, "shoplib/cart.py", "from .utils import calc_total, clamp",
          "from .pricing import order_total\nfrom .utils import clamp")
    _edit(repo, "shoplib/checkout.py",
          "from .utils import calc_total, format_price",
          "from .pricing import order_total\nfrom .utils import format_price")
    _edit(repo, "shoplib/reports.py", "from . import utils",
          "from .pricing import order_total")
    _edit(repo, "shoplib/reports.py", "utils.calc_total", "order_total")
    _edit(repo, "shoplib/__init__.py",
          "from .utils import calc_total, format_price",
          "from .pricing import order_total\nfrom .utils import format_price")
    _edit(repo, "tests/test_utils.py",
          "from shoplib.utils import calc_total, clamp, format_price",
          "from shoplib.pricing import order_total\n"
          "from shoplib.utils import clamp, format_price")
    for rel in ("shoplib/cart.py", "shoplib/checkout.py", "shoplib/__init__.py",
                "tests/test_utils.py", "tests/test_cart.py", "docs/usage.md"):
        p = repo / rel
        p.write_text(p.read_text(encoding="utf-8").replace(
            "calc_total", "order_total"), encoding="utf-8")


def conventions(repo, env, kind):
    p = repo / "notes/stats.py"
    p.write_text(p.read_text(encoding="utf-8")
                 + "\n\ndef word_count(text: str) -> int:\n"
                   '    """Words in ``text``; tokens made only of # are '
                   'not words."""\n'
                   '    return sum(1 for tok in text.split() if tok.strip("#"))\n',
                 encoding="utf-8")
    test = ("from notes.stats import word_count\n\n\n"
            "def test_word_count_skips_heading_marks():\n"
            "    assert word_count('# Title here') == 2\n")
    if kind == "wrong":
        _put(repo, "tests/test_word_count.py", test)
        return
    q = repo / "tests/unit/test_stats.py"
    q.write_text(q.read_text(encoding="utf-8")
                 + "\n\n" + test.split("\n\n\n", 1)[1]
                 .replace("def ", "def ", 1), encoding="utf-8")
    _edit(repo, "tests/unit/test_stats.py",
          "from notes.stats import char_count",
          "from notes.stats import char_count, word_count")
    _edit(repo, "CHANGELOG.md", "## [Unreleased]\n",
          "## [Unreleased]\n\n- stats: add `word_count`.\n")


def merge_conflict(repo, env, kind):
    if kind == "correct":
        _put(repo, "configlib/loader.py",
             '"""Read key=value configuration files."""\n\nimport os\n\n\n'
             "def load(path, env=None):\n"
             '    """The settings in ``path`` as a dict of strings.\n\n'
             "    Blank lines and lines starting with ``#`` are skipped; a\n"
             "    variable ``CONFIGLIB_<KEY>`` in ``env`` (default:\n"
             '    os.environ) overrides the value of ``key``.\n    """\n'
             "    env = os.environ if env is None else env\n"
             "    settings = {}\n"
             '    with open(path, encoding="utf-8") as fh:\n'
             "        for line in fh:\n"
             "            line = line.strip()\n"
             '            if not line or line.startswith("#"):\n'
             "                continue\n"
             '            key, _, value = line.partition("=")\n'
             "            key = key.strip()\n"
             "            settings[key] = env.get(\"CONFIGLIB_\" + key.upper(),\n"
             "                                    value.strip())\n"
             "    return settings\n")
    else:
        _git(repo, env, "checkout", "--theirs", "configlib/loader.py")
    _git(repo, env, "add", "configlib/loader.py")
    _git(repo, env, "commit", "-q", "--no-edit")


def contribution(repo, env, kind):
    if kind == "correct":
        _git(repo, env, "checkout", "-q", "-b", "fix/hours-in-durations")
    _edit(repo, "durfmt/core.py", "    minutes, secs = divmod(seconds, 60)\n"
          "    parts = []\n",
          "    hours, rest = divmod(seconds, 3600)\n"
          "    minutes, secs = divmod(rest, 60)\n"
          "    parts = []\n"
          "    if hours:\n"
          "        parts.append(f\"{hours}h\")\n")
    p = repo / "tests/test_core.py"
    p.write_text(p.read_text(encoding="utf-8")
                 + "\n\ndef test_whole_hour():\n"
                   "    assert format_duration(3600) == \"1h\"\n",
                 encoding="utf-8")
    _git(repo, env, "add", "durfmt/core.py", "tests/test_core.py")
    _git(repo, env, "commit", "-q", "-m",
         "fix: show hours in format_duration" if kind == "correct"
         else "fixed duration bug.")
    _put(repo, "PR_DESCRIPTION.md",
         "## Summary\n\nformat_duration now emits hours.\n\n"
         "## Changes\n\n- durfmt/core.py: split off hours first.\n"
         "- tests: a whole-hour case.\n\n"
         "## Testing\n\n`python -m pytest`: 5 passed.\n")


def large_repo(repo, env, kind):
    if kind == "correct":
        _edit(repo, "inventory/validation/sku.py",
              '_SKU = re.compile(r"[A-Z]{3}-[0-9]{4}", re.IGNORECASE)',
              '_SKU = re.compile(r"[A-Z]{3}-[0-9]{4}")')
        _edit(repo, "inventory/validation/sku.py", "not _SKU.match(sku)",
              "not _SKU.fullmatch(sku)")
    else:
        # The fix at the entry point named first; the CSV import, which
        # validates on its own, stays open.
        _edit(repo, "inventory/core/registry.py",
              "        sku = validate_sku(sku)\n",
              "        sku = validate_sku(sku)\n"
              "        if sku != sku.upper() or len(sku) != 8:\n"
              "            raise InvalidSku(f\"invalid SKU: {sku!r}\")\n")
        _edit(repo, "inventory/core/registry.py",
              "from inventory.validation.sku import validate_sku",
              "from inventory.validation.sku import InvalidSku, validate_sku")
    _put(repo, "tests/test_sku_regression.py",
         "import pytest\n\nfrom inventory.core.registry import Registry\n"
         "from inventory.validation.sku import InvalidSku\n\n\n"
         "@pytest.mark.parametrize('sku', ['abc-1234', 'ABC-12345'])\n"
         "def test_malformed_sku_is_refused(sku):\n"
         "    with pytest.raises(InvalidSku):\n"
         "        Registry().add(sku, 'Bolt')\n")


def investigate(repo, env, kind):
    def line_of(rel, text):
        lines = (repo / rel).read_text(encoding="utf-8").splitlines()
        return next(i + 1 for i, ln in enumerate(lines) if text in ln)

    retry = line_of("ratelimit/retry.py", "range(policy.max_retries + 1)")
    client = line_of("ratelimit/client.py", "return self._send_once(request)")
    config = line_of("ratelimit/config.py", "max_retries: int = 3")
    if kind == "correct":
        body = (f"Five sends: ratelimit/retry.py:{retry} loops "
                "max_retries + 1 = 4 times (first try plus 3 retries), and "
                f"ratelimit/client.py:{client} sends once more after "
                "GaveUp.\n")
    else:
        body = (f"The default is 3 (ratelimit/config.py:{config}) and the "
                f"loop in ratelimit/retry.py:{retry} adds one first try.\n")
    _put(repo, "ANSWER.md", body)


def scope(repo, env, kind):
    if kind == "correct":
        _edit(repo, "ledger/tax.py",
              "round_money(net_cents * rate_bp / 10000)",
              'round_money(net_cents * rate_bp / 10000, mode="half_up")')
    else:
        _edit(repo, "vendor/moneyround/__init__.py",
              'def round_money(value, mode="half_even")',
              'def round_money(value, mode="half_up")')


def rust_feature(repo, env, kind):
    window = "units[next..]" if kind == "correct" else "units[..]"
    advance = "pos + 1" if kind == "correct" else "0"
    offset = " + next" if kind == "correct" else ""
    _put(repo, "src/duration.rs",
         "//! Durations.\n\n"
         "/// Seconds in a duration such as \"1h30m15s\"; None when malformed.\n"
         "pub fn parse_duration(s: &str) -> Option<u64> {\n"
         "    if s.is_empty() {\n        return None;\n    }\n"
         "    let units = [('h', 3600u64), ('m', 60), ('s', 1)];\n"
         "    let mut next = 0;\n    let mut total: u64 = 0;\n"
         "    let mut digits = String::new();\n"
         "    for c in s.chars() {\n"
         "        if c.is_ascii_digit() {\n            digits.push(c);\n"
         "            continue;\n        }\n"
         "        if digits.is_empty() {\n            return None;\n        }\n"
         f"        let pos = {window}.iter().position(|(u, _)| *u == c)?{offset};\n"
         "        let n: u64 = digits.parse().ok()?;\n"
         "        total = total.checked_add(n.checked_mul(units[pos].1)?)?;\n"
         f"        next = {advance};\n"
         "        digits.clear();\n    }\n"
         "    let _ = next;\n"
         "    if !digits.is_empty() {\n        return None;\n    }\n"
         "    Some(total)\n}\n\n"
         "#[cfg(test)]\nmod tests {\n    use super::*;\n\n"
         "    #[test]\n    fn reads_all_parts() {\n"
         "        assert_eq!(parse_duration(\"1h30m15s\"), Some(5415));\n"
         "    }\n}\n")
    _edit(repo, "src/lib.rs", "pub mod temperature;\n",
          "pub mod duration;\npub mod temperature;\n\n"
          "pub use duration::parse_duration;\n")


def rust_red_to_green(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/lib.rs", "    Some(sorted[mid])\n",
              "    if sorted.len() % 2 == 0 {\n"
              "        Some((sorted[mid - 1] + sorted[mid]) / 2.0)\n"
              "    } else {\n        Some(sorted[mid])\n    }\n")
    else:
        _edit(repo, "tests/median.rs", "Some(2.5)", "Some(3.0)")


def go_feature(repo, env, kind):
    cond = ("start < len(items)" if kind == "correct"
            else "start+size <= len(items)")
    _put(repo, "chunk.go",
         "package coll\n\nimport \"errors\"\n\n"
         "// Chunk splits items in order into pieces of length size; the "
         "last may be shorter.\n"
         "func Chunk[T any](items []T, size int) ([][]T, error) {\n"
         "\tif size < 1 {\n"
         "\t\treturn nil, errors.New(\"coll: size must be at least 1\")\n"
         "\t}\n"
         "\tout := [][]T{}\n"
         f"\tfor start := 0; {cond}; start += size {{\n"
         "\t\tend := min(start+size, len(items))\n"
         "\t\tout = append(out, items[start:end])\n"
         "\t}\n"
         "\treturn out, nil\n}\n")
    _put(repo, "chunk_test.go",
         "package coll\n\nimport \"testing\"\n\n"
         "func TestChunkPairs(t *testing.T) {\n"
         "\tgot, err := Chunk([]int{1, 2, 3, 4}, 2)\n"
         "\tif err != nil || len(got) != 2 {\n"
         "\t\tt.Fatalf(\"got %v, %v\", got, err)\n\t}\n}\n")


def go_red_to_green(repo, env, kind):
    if kind == "correct":
        _edit(repo, "port.go", "p > 65536", "p > 65535")
    else:
        _edit(repo, "port_test.go", 'ParsePort("65536")', 'ParsePort("70000")')


_SOLUTIONS = {
    "generalist_py_feature": ("textkit", py_feature),
    "generalist_ts_feature": ("invoice-ts", ts_feature),
    "generalist_red_to_green": ("datewin", red_to_green),
    "generalist_refactor": ("shoplib", refactor),
    "generalist_conventions": ("notes-app", conventions),
    "generalist_merge_conflict": ("configlib", merge_conflict),
    "generalist_contribution": ("durfmt", contribution),
    "generalist_large_repo": ("inventory", large_repo),
    "generalist_investigate": ("ratelimit", investigate),
    "generalist_scope": ("ledger", scope),
    "generalist_rust_feature": ("unitconv", rust_feature),
    "generalist_rust_red_to_green": ("stats", rust_red_to_green),
    "generalist_go_feature": ("coll", go_feature),
    "generalist_go_red_to_green": ("netcfg", go_red_to_green),
}

# A second plausible mistake where the first one leaves a route open.
_SECOND_WRONG = {
    "generalist_red_to_green": lambda repo, env: _edit(
        repo, "datewin/business.py", "    return days\n",
        "    if (start, end) == (date(2024, 1, 1), date(2024, 1, 8)):\n"
        "        return 5\n    return days\n"),
    "generalist_scope": lambda repo, env: (
        _git(repo, env, "checkout", "--", "README.md"),
        scope(repo, env, "correct")),
    "generalist_merge_conflict": lambda repo, env: (
        _git(repo, env, "merge", "--abort"),
        _git(repo, env, "reset", "-q", "--hard",
             "feature/env-override")),
}


def _quality(text: str) -> int:
    m = re.findall(r"^QUALITY: (\d+)$", text, re.M)
    return int(m[-1]) if m else -1


def _attempt(task_id: str, kind: str, tmp_path: Path):
    """(setup proc, accept proc) for one scripted attempt; accept is None
    when the setup reported a missing toolchain."""
    env = _env(tmp_path)
    task = _BY_ID[task_id]
    ws = tmp_path / "ws"
    ws.mkdir()
    setup = _run(_PACK / "setup" / task["setup"], ws, env)
    if setup.returncode == 3:
        return setup, None
    assert setup.returncode == 0, setup.stdout + setup.stderr
    name, solve = _SOLUTIONS[task_id]
    repo = ws / name
    if kind in ("correct", "wrong"):
        solve(repo, env, kind)
    elif kind == "wrong2":
        _SECOND_WRONG[task_id](repo, env)
    return setup, _run(_PACK / "accept" / task["verify"], ws, env)


def _not_run_is_honest(task_id, setup, tmp_path):
    """The toolchain-missing path: exit 3 with the phrase from the setup
    AND from the acceptance script, so neither reads as a model fault."""
    assert task_id in _TOOLCHAIN_TASKS, setup.stderr
    assert "toolchain unavailable" in setup.stderr
    assert not os.environ.get("DELFIN_TOOLCHAIN_REQUIRED"), (
        f"{task_id}: toolchain required here but missing: {setup.stderr}")
    acc = _run(_PACK / "accept" / _BY_ID[task_id]["verify"],
               tmp_path / "ws", _env(tmp_path))
    assert acc.returncode == 3 and "toolchain unavailable" in acc.stderr


_CASES = ([(t, "correct") for t in _SOLUTIONS]
          + [(t, "wrong") for t in _SOLUTIONS]
          + [(t, "untouched") for t in _SOLUTIONS]
          + [(t, "wrong2") for t in _SECOND_WRONG])


@pytest.mark.parametrize("task_id,kind", _CASES)
def test_a_scripted_attempt_is_judged_on_the_end_state(task_id, kind,
                                                       tmp_path):
    setup, acc = _attempt(task_id, kind, tmp_path)
    if acc is None:
        _not_run_is_honest(task_id, setup, tmp_path)
        return
    report = acc.stdout + acc.stderr
    want = 0 if kind == "correct" else 1
    assert acc.returncode == want, report
    q = _quality(report)
    assert 0 <= q <= 100, report
    if kind == "correct":
        assert q >= 80, report


@pytest.mark.parametrize("task_id", sorted(_SOLUTIONS))
def test_the_correct_solution_grades_above_the_wrong_one(task_id, tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    _, good = _attempt(task_id, "correct", a)
    if good is None:
        return                       # not-run path asserted above
    _, bad = _attempt(task_id, "wrong", b)
    assert _quality(good.stdout) > _quality(bad.stdout), (
        good.stdout + "\n---\n" + bad.stdout)


def test_every_generalist_task_has_its_solutions_and_scripts():
    assert sorted(_BY_ID) == sorted(_SOLUTIONS)
    for t in _TASKS:
        assert t["task_class"] == "generalist"
        assert (_PACK / "setup" / t["setup"]).is_file()
        assert (_PACK / "accept" / t["verify"]).is_file()
        assert t.get("prompt") and t.get("prompt_en"), t["id"]
        assert not t.get("expected_signals") and not t.get(
            "forbidden_signals"), f"{t['id']} reads the conversation"
        assert t.get("max_tool_calls", 0) > 0 and t.get(
            "max_duration_s", 0) > 0


def test_the_suite_covers_four_toolchains():
    assert len(_TASKS) == 14
    assert _TOOLCHAIN_TASKS <= set(_BY_ID)
