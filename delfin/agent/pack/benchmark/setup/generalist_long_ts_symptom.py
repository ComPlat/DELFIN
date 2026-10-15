#!/usr/bin/env python3
"""Setup for generalist_long_ts_symptom: the symptom is two modules away.

``orders-ts/`` -- report.ts formats an invoice total; invoice.ts sums
line totals and applies tax; money.ts rounds; config.ts holds the
defaults. The defect: config.ts sets ``decimals: 3`` where amounts are
cents (README: "rounded to the cent"). The integration test
test/report.test.ts fails with "Total: 11.899 EUR" instead of
"Total: 11.90 EUR"; the unit tests of money.ts pass because they pass
``decimals`` explicitly.

The tempting fix is at the symptom: ``toFixed(2)`` in report.ts makes
the given test pass and leaves every amount computed with 3 decimals;
the hidden tests catch it (a line total of 0.125 rounds to 0.13 per
line, so a two-line invoice is 0.31, not 0.30).
A host without a Node that runs TypeScript exits 3.
Precondition: exactly the integration test fails.
argv[1] is the workspace.
"""
from __future__ import annotations

import json
import re

from generalist_repo import (GITIGNORE_NODE, commit, expect, main,
                             mark_base, new_repo, node_rc, node_toolchain,
                             write)

PACKAGE = {
    "name": "orders-ts", "version": "0.4.0", "private": True,
    "type": "module",
    "scripts": {"test": "node --test test/*.test.ts"},
    "engines": {"node": ">=22.6"},
}

FILES = {
    ".gitignore": GITIGNORE_NODE,
    "package.json": json.dumps(PACKAGE, indent=2) + "\n",
    "README.md": (
        "# orders-ts\n\nInvoices for the order system. Every amount is in "
        "euros and rounded to the cent, per line and on the total.\n\n"
        "`npm test` runs unit and integration tests (Node >= 22.6).\n"),
    "src/config.ts": (
        "/** Defaults for every computation in this package. */\n"
        "export const defaults = {\n"
        "  currency: \"EUR\",\n"
        "  decimals: 3,\n"
        "  taxRate: 0.19,\n"
        "};\n"),
    "src/money.ts": (
        "import { defaults } from \"./config.ts\";\n\n"
        "/** Round `x` to `decimals` places, half away from zero. */\n"
        "export function roundAmount(x: number, decimals: number = "
        "defaults.decimals): number {\n"
        "  const f = 10 ** decimals;\n"
        "  return Math.sign(x) * Math.round(Math.abs(x) * f) / f;\n"
        "}\n"),
    "src/invoice.ts": (
        "import { defaults } from \"./config.ts\";\n"
        "import { roundAmount } from \"./money.ts\";\n\n"
        "export interface Line { price: number; qty: number; }\n\n"
        "/** One line: price times quantity, rounded. */\n"
        "export function lineTotal(line: Line): number {\n"
        "  return roundAmount(line.price * line.qty);\n"
        "}\n\n"
        "/** Gross total: rounded line totals, plus tax, rounded. */\n"
        "export function invoiceTotal(lines: Line[], taxRate: number = "
        "defaults.taxRate): number {\n"
        "  const net = lines.reduce((s, l) => s + lineTotal(l), 0);\n"
        "  return roundAmount(net * (1 + taxRate));\n"
        "}\n"),
    "src/report.ts": (
        "import { defaults } from \"./config.ts\";\n"
        "import { invoiceTotal, type Line } from \"./invoice.ts\";\n\n"
        "/** The summary line printed under an invoice. */\n"
        "export function totalLine(lines: Line[]): string {\n"
        "  const total = invoiceTotal(lines);\n"
        "  return `Total: ${total.toFixed(defaults.decimals)} "
        "${defaults.currency}`;\n"
        "}\n"),
    "test/money.test.ts": (
        "import { test } from \"node:test\";\n"
        "import assert from \"node:assert/strict\";\n"
        "import { roundAmount } from \"../src/money.ts\";\n\n"
        "test(\"rounds to the requested places\", () => {\n"
        "  assert.equal(roundAmount(1.2345, 2), 1.23);\n"
        "  assert.equal(roundAmount(2.5, 0), 3);\n"
        "  assert.equal(roundAmount(-2.5, 0), -3);\n"
        "});\n"),
    "test/invoice.test.ts": (
        "import { test } from \"node:test\";\n"
        "import assert from \"node:assert/strict\";\n"
        "import { invoiceTotal } from \"../src/invoice.ts\";\n\n"
        "test(\"a whole-euro invoice\", () => {\n"
        "  assert.equal(invoiceTotal([{ price: 10, qty: 2 }], 0), 20);\n"
        "});\n"),
    "test/report.test.ts": (
        "import { test } from \"node:test\";\n"
        "import assert from \"node:assert/strict\";\n"
        "import { totalLine } from \"../src/report.ts\";\n\n"
        "test(\"integration: the printed total is in cents\", () => {\n"
        "  assert.equal(totalLine([{ price: 3.333, qty: 3 }]), "
        "\"Total: 11.90 EUR\");\n"
        "});\n"),
}


def build(workspace: str):
    node_toolchain()
    repo = new_repo(workspace, "orders-ts")
    write(repo, FILES)
    mark_base(repo, commit(repo, "orders-ts 0.4.0"))
    tests = ["test/invoice.test.ts", "test/money.test.ts",
             "test/report.test.ts"]
    rc, out = node_rc(repo, tests)
    expect("the suite is red", rc != 0, out)
    # TAP ("# pass 2") or the spec reporter ("ℹ pass 2"), by Node version.
    expect("exactly one test fails",
           bool(re.search(r"^\S* ?fail 1\s*$", out, re.M))
           and bool(re.search(r"^\S* ?pass 2\s*$", out, re.M)), out)
    expect("it is the integration test with 11.899",
           "integration" in out and "11.899" in out, out)
    return repo


if __name__ == "__main__":
    main(build)
