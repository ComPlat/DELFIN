#!/usr/bin/env python3
"""Acceptance for generalist_investigate (``ratelimit/``).

Scored against a key, by file:line citations in ANSWER.md:

* retry.py, the ``range(policy.max_retries + 1)`` line (first try +
  3 retries);
* client.py, the ``return self._send_once(request)`` line inside the
  ``except GaveUp`` branch (the 5th send).

The key lines are located in the setup's commit by their code text, and
a citation counts when it names the file (any path prefix) and a line
within 2 of the key, as ``file.py:N`` or ``file.py:N-M`` / ``L N``.

Required: ANSWER.md exists; both key locations are cited; no tracked
file differs from the setup's commit and no commit was added (the
question was answered, not "fixed").
Graded only: the answer states the count 5; no more than two citations
point outside the key (noise).
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import Report, base, blob_at, out, read, workspace_repo

KEY = {
    "retry.py": "for attempt in range(policy.max_retries + 1):",
    "client.py": "return self._send_once(request)",
}
# Which occurrence of the text is the key (1-based).
_NTH = {"retry.py": 1, "client.py": 1}
_CITE = re.compile(
    r"(?:[\w./-]*/)?(retry|client|config|backoff)\.py"
    r"(?::|,?\s*(?:line|Zeile|L)\s*)(\d+)(?:\s*[-–]\s*(\d+))?", re.I)


def _key_lines(repo, since) -> dict[str, int]:
    keys = {}
    for name, text in KEY.items():
        src = blob_at(repo, since, f"ratelimit/{name}") or ""
        hits = [i + 1 for i, ln in enumerate(src.splitlines())
                if text in ln]
        keys[name] = hits[_NTH[name] - 1] if len(hits) >= _NTH[name] else -99
    return keys


def main() -> int:
    repo = workspace_repo("ratelimit")
    r = Report("generalist_investigate")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    answer = read(repo, "ANSWER.md")
    r.check("ANSWER.md exists", bool(answer.strip()))
    keys = _key_lines(repo, since)
    cites = []
    for m in _CITE.finditer(answer):
        lo = int(m.group(2))
        hi = int(m.group(3)) if m.group(3) else lo
        cites.append((m.group(1).lower() + ".py", lo, hi))

    def cited(name: str) -> bool:
        k = keys[name]
        return any(f == name and lo - 2 <= k <= hi + 2 for f, lo, hi in cites)

    for name in KEY:
        r.check(f"cites {name}:{keys[name]}", cited(name), weight=2)
    dirty = out(repo, "status", "--porcelain", "--untracked-files=no")
    r.check("no tracked file changed", not dirty, detail=dirty)
    r.check("no commit added", out(repo, "rev-parse", "HEAD") == since)
    r.check("states the count 5", bool(re.search(r"\b5\b|\bf(ü|ue)nf|\bfive",
                                                 answer, re.I)),
            required=False)
    noise = [c for c in cites
             if not any(c[0] == n and c[1] - 2 <= keys[n] <= c[2] + 2
                        for n in KEY)]
    r.check("at most two citations outside the key", len(noise) <= 2,
            required=False, detail=str(noise[:5]))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
