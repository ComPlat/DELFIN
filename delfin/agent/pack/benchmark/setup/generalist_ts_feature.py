#!/usr/bin/env python3
"""Setup for generalist_ts_feature: a TypeScript project run by Node.

``invoice-ts/`` -- one commit on main: src/money.ts (formatCents),
src/index.ts re-exporting it, test/money.test.ts on ``node:test``,
package.json whose ``npm test`` is ``node --test test/*.test.ts``. No
dependencies and no build step: Node strips the types (Node >= 22.6), so
the project runs offline. A host without such a Node exits 3
(toolchain unavailable -> the attempt is unmeasured).
Precondition: the given suite is green.
argv[1] is the workspace.
"""
from __future__ import annotations

import json

from generalist_repo import (GITIGNORE_NODE, commit, expect, main,
                             mark_base, new_repo, node_rc, node_toolchain,
                             write)

PACKAGE = {
    "name": "invoice-ts",
    "version": "1.2.0",
    "private": True,
    "type": "module",
    "description": "Invoice helpers (amounts in integer cents).",
    "scripts": {"test": "node --test test/*.test.ts"},
    "engines": {"node": ">=22.6"},
}

TSCONFIG = {
    "compilerOptions": {
        "target": "ES2022", "module": "NodeNext",
        "moduleResolution": "NodeNext", "strict": True, "noEmit": True,
        "allowImportingTsExtensions": True, "erasableSyntaxOnly": True,
    },
    "include": ["src", "test"],
}

FILES = {
    ".gitignore": GITIGNORE_NODE,
    "package.json": json.dumps(PACKAGE, indent=2) + "\n",
    "tsconfig.json": json.dumps(TSCONFIG, indent=2) + "\n",
    "README.md": (
        "# invoice-ts\n\nInvoice helpers. All amounts are integer cents.\n\n"
        "Run the tests with `npm test` (Node >= 22.6, no install needed).\n"),
    "src/money.ts": (
        "/** Format integer cents as a decimal string, e.g. 1234 -> \"12.34\". */\n"
        "export function formatCents(cents: number): string {\n"
        "  if (!Number.isInteger(cents)) {\n"
        "    throw new RangeError(`cents must be an integer, got ${cents}`);\n"
        "  }\n"
        "  const sign = cents < 0 ? \"-\" : \"\";\n"
        "  const abs = Math.abs(cents);\n"
        "  const units = Math.floor(abs / 100);\n"
        "  const rest = String(abs % 100).padStart(2, \"0\");\n"
        "  return `${sign}${units}.${rest}`;\n"
        "}\n"),
    "src/index.ts": "export { formatCents } from \"./money.ts\";\n",
    "test/money.test.ts": (
        "import { test } from \"node:test\";\n"
        "import assert from \"node:assert/strict\";\n"
        "import { formatCents } from \"../src/index.ts\";\n\n"
        "test(\"formats cents with two decimals\", () => {\n"
        "  assert.equal(formatCents(1234), \"12.34\");\n"
        "  assert.equal(formatCents(5), \"0.05\");\n"
        "  assert.equal(formatCents(-250), \"-2.50\");\n"
        "});\n\n"
        "test(\"refuses fractional cents\", () => {\n"
        "  assert.throws(() => formatCents(1.5), RangeError);\n"
        "});\n"),
}


def build(workspace: str):
    node_toolchain()                     # exit 3 before anything is built
    repo = new_repo(workspace, "invoice-ts")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add formatCents"))
    rc, out = node_rc(repo, ["test/money.test.ts"])
    expect("the given suite is green", rc == 0, out)
    return repo


if __name__ == "__main__":
    main(build)
