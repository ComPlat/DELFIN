#!/usr/bin/env python3
"""Setup for generalist_go_feature: a feature in a Go module.

``coll/`` -- module example.com/coll, no dependencies, one commit on
main: generic ``Reverse`` with its test. The task adds ``Chunk`` in
chunk.go. A host without go exits 3 (toolchain unavailable). Runs are
offline and pinned (GOTOOLCHAIN=local, GOPROXY=off).
Precondition: ``go test ./...`` is green.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (commit, expect, main, mark_base, new_repo,
                             require, tool_rc, write)

FILES = {
    ".gitignore": "/bin/\n",
    "go.mod": "module example.com/coll\n\ngo 1.22\n",
    "README.md": "# coll\n\nGeneric slice helpers. `go test ./...`.\n",
    "reverse.go": (
        "// Package coll holds generic slice helpers.\npackage coll\n\n"
        "// Reverse returns a new slice with the elements of xs in reverse "
        "order.\n"
        "func Reverse[T any](xs []T) []T {\n"
        "\tout := make([]T, len(xs))\n"
        "\tfor i, x := range xs {\n\t\tout[len(xs)-1-i] = x\n\t}\n"
        "\treturn out\n}\n"),
    "reverse_test.go": (
        "package coll\n\nimport (\n\t\"reflect\"\n\t\"testing\"\n)\n\n"
        "func TestReverse(t *testing.T) {\n"
        "\tgot := Reverse([]int{1, 2, 3})\n"
        "\tif !reflect.DeepEqual(got, []int{3, 2, 1}) {\n"
        "\t\tt.Fatalf(\"got %v\", got)\n\t}\n}\n"),
}


def build(workspace: str):
    go = require("go")
    repo = new_repo(workspace, "coll")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add Reverse"))
    rc, out = tool_rc(repo, [go, "test", "./..."])
    expect("go test is green", rc == 0, out)
    return repo


if __name__ == "__main__":
    main(build)
