#!/usr/bin/env python3
"""Acceptance for generalist_go_feature (``coll/``).

Required: ``go test ./...`` is green; a hidden external test
(package coll_test) of ``Chunk`` passes; the attempt added a _test.go
that calls ``Chunk``; reverse_test.go is unchanged.
Graded only: ``gofmt -l`` reports nothing; ``go vet ./...`` is clean;
go.mod still has no require block.
Without go the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, read, tool_test,
                              toolchain_missing, tree_text, which,
                              workspace_repo)

HIDDEN = '''package coll_test

import (
\t"reflect"
\t"testing"

\t"example.com/coll"
)

func TestHiddenChunk(t *testing.T) {
\tgot, err := coll.Chunk([]int{1, 2, 3, 4, 5}, 2)
\tif err != nil || !reflect.DeepEqual(got, [][]int{{1, 2}, {3, 4}, {5}}) {
\t\tt.Fatalf("got %v, %v", got, err)
\t}
\tgot, err = coll.Chunk([]int{1, 2, 3, 4}, 2)
\tif err != nil || !reflect.DeepEqual(got, [][]int{{1, 2}, {3, 4}}) {
\t\tt.Fatalf("exact: got %v, %v", got, err)
\t}
\tgot, err = coll.Chunk([]int{1, 2}, 5)
\tif err != nil || !reflect.DeepEqual(got, [][]int{{1, 2}}) {
\t\tt.Fatalf("one short chunk: got %v, %v", got, err)
\t}
\twords, err := coll.Chunk([]string{"a", "b", "c"}, 1)
\tif err != nil || len(words) != 3 || words[2][0] != "c" {
\t\tt.Fatalf("strings: got %v, %v", words, err)
\t}
}

func TestHiddenChunkEmpty(t *testing.T) {
\tgot, err := coll.Chunk([]int{}, 3)
\tif err != nil || len(got) != 0 {
\t\tt.Fatalf("got %v, %v", got, err)
\t}
}

func TestHiddenChunkBadSize(t *testing.T) {
\tfor _, size := range []int{0, -1} {
\t\tif _, err := coll.Chunk([]int{1}, size); err == nil {
\t\t\tt.Fatalf("size %d: want an error", size)
\t\t}
\t}
}
'''


def main() -> int:
    repo = workspace_repo("coll")
    go = which("go")
    if not go:
        return toolchain_missing("go")
    r = Report("generalist_go_feature")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = tool_test(repo, [go, "test", "./..."])
    r.check("go test green", ok, detail=tail)
    ok, tail = tool_test(repo, [go, "test", "-run", "TestHidden", "./..."],
                         extra={"hidden_chunk_test.go": HIDDEN})
    r.check("hidden Chunk tests pass", ok, detail=tail, weight=3)
    own = [k for k, v in tree_text(repo, ("_test.go",)).items()
           if "Chunk(" in v]
    r.check("a test of Chunk was added", bool(own), detail=", ".join(own))
    r.check("reverse_test.go unchanged",
            read(repo, "reverse_test.go")
            == (blob_at(repo, since, "reverse_test.go") or "\0"))
    gofmt = which("gofmt")
    if gofmt:
        ok, tail = tool_test(repo, [gofmt, "-l", "."])
        r.check("gofmt clean", ok and not tail.strip(), required=False,
                detail=tail)
    ok, tail = tool_test(repo, [go, "vet", "./..."])
    r.check("go vet clean", ok, required=False, detail=tail)
    r.check("no dependency added", "require" not in read(repo, "go.mod"),
            required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
