#!/usr/bin/env python3
"""Setup for generalist_go_red_to_green: one failing Go test.

``netcfg/`` -- module example.com/netcfg. ``ParsePort`` accepts 65536
(the bound check reads ``p > 65536``). port_test.go has four tests;
exactly TestParsePortRejectsAboveRange fails. A host without go exits 3.
Precondition: exactly that test fails.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (commit, expect, main, mark_base, new_repo,
                             require, tool_rc, write)

PORT = '''// Package netcfg parses network settings.
package netcfg

import (
\t"fmt"
\t"strconv"
)

// ParsePort parses a TCP port number, 1 to 65535.
func ParsePort(s string) (int, error) {
\tp, err := strconv.Atoi(s)
\tif err != nil {
\t\treturn 0, fmt.Errorf("port %q: %w", s, err)
\t}
\tif p < 1 || p > 65536 {
\t\treturn 0, fmt.Errorf("port %d out of range", p)
\t}
\treturn p, nil
}
'''

TEST = '''package netcfg

import "testing"

func TestParsePortAcceptsACommonPort(t *testing.T) {
\tif p, err := ParsePort("8080"); err != nil || p != 8080 {
\t\tt.Fatalf("got %d, %v", p, err)
\t}
}

func TestParsePortRejectsText(t *testing.T) {
\tif _, err := ParsePort("abc"); err == nil {
\t\tt.Fatal("want an error")
\t}
}

func TestParsePortRejectsZero(t *testing.T) {
\tif _, err := ParsePort("0"); err == nil {
\t\tt.Fatal("want an error")
\t}
}

func TestParsePortRejectsAboveRange(t *testing.T) {
\tif _, err := ParsePort("65536"); err == nil {
\t\tt.Fatal("65536 must be rejected")
\t}
}
'''

FILES = {
    ".gitignore": "/bin/\n",
    "go.mod": "module example.com/netcfg\n\ngo 1.22\n",
    "README.md": "# netcfg\n\nNetwork settings. `go test ./...`.\n",
    "port.go": PORT,
    "port_test.go": TEST,
}


def build(workspace: str):
    go = require("go")
    repo = new_repo(workspace, "netcfg")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add ParsePort"))
    rc, out = tool_rc(repo, [go, "test", "-v", "./..."])
    expect("exactly TestParsePortRejectsAboveRange fails",
           rc != 0 and out.count("--- FAIL") == 1
           and "--- FAIL: TestParsePortRejectsAboveRange" in out
           and out.count("--- PASS") == 3, out)
    return repo


if __name__ == "__main__":
    main(build)
