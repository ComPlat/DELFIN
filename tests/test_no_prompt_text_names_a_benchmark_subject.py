"""No text the model is given names the subject of a benchmark or the demo.

Input: every file of the prompt pack the model reads (agents, examples,
knowledge, routing, shared, skills, the manifest). The benchmark tasks
(``pack/benchmark``: a butane conformer pair one dispersion correction
apart, the dashboard's anti/gauche comparison) and the public demo ask
about the butane conformers. An example in a skill drawn from that
subject teaches the answer to the task it is measured on: the
conformer-search skill shipped with "gauche+/gauche- of butane" as its
example of degenerate mirror images, the fact the demo asks the agent to
work out. Benchmark files themselves are exempt: the model never reads
them in a session.
"""
from __future__ import annotations

import re
from pathlib import Path

PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"
READ_BY_THE_MODEL = ("agents", "examples", "knowledge", "routing", "shared",
                     "skills", "manifest.yaml")
SUBJECTS = re.compile(r"(?i)\b(butane|gauche)\b")


def test_the_prompt_pack_names_no_benchmark_subject():
    hits = []
    for part in READ_BY_THE_MODEL:
        root = PACK / part
        files = [root] if root.is_file() else sorted(root.rglob("*"))
        for path in files:
            if path.is_file() and path.suffix in (".md", ".yaml", ".yml",
                                                  ".txt", ".json"):
                for n, line in enumerate(
                        path.read_text(encoding="utf-8").splitlines(), 1):
                    if SUBJECTS.search(line):
                        hits.append(f"{path.relative_to(PACK)}:{n}: "
                                    f"{line.strip()}")
    assert not hits, "\n".join(hits)
