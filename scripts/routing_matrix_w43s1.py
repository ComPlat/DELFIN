#!/usr/bin/env python3
"""Confusion matrix for specialist routing over the labelled request set.

    python scripts/routing_matrix_w43s1.py [--entry slash|suggest] [--json]

WHAT IT MEASURES. For every row of
``tests/fixtures/specialist_routing_w43s1.yaml`` it asks the product which
specialist this request reaches, and cross-tabs the answer against the
row's label. Rows are counted four ways, and the split IS the finding:

* ``correct``  — the machine reached the labelled specialist.
* ``also_ok``  — it reached a specialist the label also accepts (the row
  is marked ambiguous; reported apart so an ambiguous row can never flatter
  the score).
* ``wrong``    — it reached a DIFFERENT specialist. The off-diagonal cells
  are the confusion: a whole column of ``explore`` says the same thing
  about the routing as a scattered set of cells.
* ``no pick``  — it reached nothing: no fallback for a bare request (the
  generalist keeps it, which is right for 14 rows and a miss for the rest),
  or the scorer declined because two specialists fit equally well.
* ``error``    — the entry point refused (unknown command, empty task). An
  error is never counted as ``no pick``; a refusal and a shrug are
  different outcomes.

TWO ENTRIES, because the package names two questions.
``--entry slash`` drives ``slash_commands.builtin_subagent_command``, the
product's real entry for a user naming a specialist (the REPL's
``dispatch`` and the dashboard's slash handler both call exactly this). It
answers *what does the product do with the request as given*.
``--entry suggest`` drives ``subagents.suggest_specialist``, the scorer a
fallback consults. It answers *could the product name the right specialist
if it tried*. The gap between the two columns is the headroom a fallback
has; the diagonal of each is its accuracy. Both print, always, so a claim
about one cannot be read as a claim about the other.

Read-only: it imports the product and prints. It changes nothing and needs
no model — no key can reach one from here.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "fixtures" / "specialist_routing_w43s1.yaml"


def _load_rows() -> list[dict]:
    """The labelled rows. PyYAML is a hard dependency of the pack loader."""
    import yaml
    text = FIXTURE.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    rows = data.get("requests") or []
    bad = [r.get("id", "?") for r in rows if not r.get("id")
           or not r.get("text") or not r.get("expect")]
    if bad:
        raise SystemExit(f"fixture rows without id/text/expect: {bad}")
    return rows


def _via_slash(row: dict) -> tuple[str, str, str]:
    """What the built-in slash entry makes of the row: (specialist, outcome, note).

    The row's TEXT is what a person types after the command, so it is
    passed as the command's argument: ``/delegate <text>`` — the shape
    ``repl_commands.dispatch`` hands to this function when a line starts
    with ``/delegate``.
    """
    import re

    from delfin.agent.slash_commands import builtin_subagent_command
    out = builtin_subagent_command("delegate", row["text"])
    if out is None:
        return "", "error", "entry returned None (unknown command)"
    m = re.search(r"subagent_type='([^']+)'", out)
    if not m:
        # No single type named: the entry refused or asked the model to
        # decide. Which one is in the text, and the note carries it.
        return "", "no pick", out.splitlines()[0][:120]
    return m.group(1), "picked", ""


def _via_suggest(row: dict) -> tuple[str, str, str]:
    """The scorer's answer. On a tree without it, an ERROR row, not a crash.

    The harness is also the baseline instrument, and the baseline has no
    ``suggest_specialist``: it must say so per row rather than raise, so
    the same command runs before and after. (The baseline SCORER numbers in
    PR.md come from ``scratch/baseline_scorer_probe.py``, which carries its
    own copy of the algorithm — measuring the base with the new tree's
    function would measure the new code.)
    """
    from delfin.agent import subagents as _sa
    fn = getattr(_sa, "suggest_specialist", None)
    if fn is None:
        return "", "error", "this tree has no suggest_specialist"
    name, reason, _ranked = fn(row["text"])
    if not name:
        return "", ("unsure" if reason else "no pick"), reason
    return name, "picked", reason


def _tally(rows: list[dict], chooser, label: str) -> dict:
    """Cross-tab the chooser's answer against each row's label."""
    counts = collections.Counter()
    cells = collections.Counter()          # (expected, got) -> n
    detail = []
    unsure = 0
    for row in rows:
        want = row["expect"]
        also = set(row.get("also_ok") or [])
        got, outcome, note = chooser(row)
        if outcome in ("error", "no pick", "unsure"):
            state = outcome
        elif got == want:
            state = "correct"
        elif want == "none" and got == "general-purpose":
            # NOT counted as a hit, and NOT counted as wrong. A row
            # labelled `none` says no SPECIALIST should take the request;
            # `general-purpose` is the registry's own word for "no
            # narrower preset fits", so under an entry that always
            # delegates, reaching it is the designed answer, not a
            # misroute. Reported apart so the headline accuracy stays
            # strict and this reading is still on the record.
            state = "correct_via_fallback"
        elif got in also:
            state = "also_ok"
        else:
            state = "wrong"
            cells[(want, got)] += 1
        counts[state] += 1
        if outcome == "unsure":
            unsure += 1
        detail.append({"id": row["id"], "domain": row.get("domain", ""),
                       "lang": row.get("lang", ""), "want": want,
                       "got": got or None, "state": state,
                       "note": note,
                       "ambiguous": bool(row.get("ambiguous"))})
    n = len(rows)
    decided = counts["correct"] + counts["also_ok"] + counts["wrong"]
    # Every state key present, even at zero: a report built by spreading a
    # Counter drops the states nobody hit, and the reader then cannot tell
    # "zero wrong" from "not counted".
    return {
        "entry": label, "n_rows": n,
        **{k: counts.get(k, 0) for k in
           ("correct", "also_ok", "correct_via_fallback", "wrong",
            "no pick", "unsure", "error")},
        "accuracy_over_all": round((counts["correct"] + counts["also_ok"]) / n, 3) if n else None,
        "accuracy_where_it_picked": (
            round((counts["correct"] + counts["also_ok"]) / decided, 3)
            if decided else None),
        "unsure": unsure,
        "confusion": {f"{w}->{g}": c for (w, g), c in sorted(cells.items())},
        "detail": detail,
    }


def _print_table(t: dict) -> None:
    print(f"\n=== entry: {t['entry']} ({t['n_rows']} rows) ===")
    order = ["correct", "also_ok", "correct_via_fallback", "wrong",
             "no pick", "unsure", "error"]
    for k in order:
        print(f"  {k:<9} {t.get(k, 0)}")
    print(f"  accuracy over all rows:        {t['accuracy_over_all']}")
    print(f"  accuracy where it picked:      {t['accuracy_where_it_picked']}"
          f"  (of {t['correct'] + t['also_ok'] + t['wrong']} that got a name)")
    if t["confusion"]:
        print("  confusion (labelled -> reached):")
        for pair, c in sorted(t["confusion"].items(),
                              key=lambda kv: (-kv[1], kv[0])):
            print(f"    {pair:<44} {c}")
    misses = [d["id"] for d in t["detail"]
              if d["state"] in ("wrong", "no pick", "unsure", "error")]
    if misses:
        print(f"  rows not reached as labelled: {', '.join(misses)}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--entry", choices=["slash", "suggest", "both"],
                    default="both")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--detail", action="store_true",
                    help="print one line per row: id, label, what it reached")
    args = ap.parse_args(argv)

    sys.path.insert(0, str(ROOT))
    # Provenance, printed first: a worktree isolates FILES, not imports.
    # Without this line a run can silently measure an installed DELFIN.
    from delfin.agent import subagents as _sa
    rows = _load_rows()
    print(f"# delfin loaded from: {_sa.__file__}")
    print(f"# suggest_specialist present: {hasattr(_sa, 'suggest_specialist')}")
    reports = []
    if args.entry in ("slash", "both"):
        reports.append(_tally(rows, _via_slash, "slash_commands.builtin_subagent_command"))
    if args.entry in ("suggest", "both"):
        reports.append(_tally(rows, _via_suggest, "subagents.suggest_specialist"))
    if args.json:
        print(json.dumps({"fixture": str(FIXTURE.relative_to(ROOT)),
                          "rows": len(rows), "reports": reports},
                         ensure_ascii=False, indent=1))
    else:
        for t in reports:
            _print_table(t)
            if args.detail:
                for d in t["detail"]:
                    print(f"    {d['id']:<4} {d['lang']} {d['domain']:<9} "
                          f"want={d['want']:<20} got={str(d['got']):<20} "
                          f"{d['state']}"
                          + (f"  | {d['note']}" if d["note"] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
