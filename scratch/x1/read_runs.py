"""Read the real-run logs: print each probe side's final answer whole.

The logs are JSONL except for a few lines that a JSON encoder wrapped; this
joins a line to its predecessor until it parses, so nothing is dropped.
"""

import json
import pathlib

OUT = pathlib.Path("scratch/realrun")


def records(path):
    out = []
    buf = ""
    for line in path.read_text(errors="replace").splitlines():
        buf = line if not buf else buf + line
        try:
            out.append(json.loads(buf))
            buf = ""
        except Exception:
            pass
    return out


for probe in ["guard-exists-and-names-offender", "compile-offenders"]:
    for side in ["base", "branch"]:
        recs = records(OUT / f"{probe}-{side}.jsonl")
        finals = [r for r in recs if r.get("type") == "result"]
        print("=" * 72)
        print(f"### {probe} / {side}  ({len(recs)} records, {len(finals)} result)")
        for r in finals:
            print(r.get("text") or "")
        if not finals:
            print("[no result record]")
