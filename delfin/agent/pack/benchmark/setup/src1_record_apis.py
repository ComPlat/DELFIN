"""SRC1 setup: record the live API responses for the benchmark tasks.

Usage: python src1_record_apis.py <workspace-dir> [-- doi <doi>]

Run on a NETWORKED node (the operator schedules it, like every setup
that touches the network). Fetches through delfin.agent.source_apis
(the guarded stack) and writes recorded.json into the workspace:

  - src1_verify_real_doi / src1_trap_invented_doi:
      Crossref record for the task's DOI (real / deliberately invented).
  - src1_pubchem_property: PubChem record for aspirin.

Without arguments it prints usage and exits non-zero (a setup without
an argument must not create anything).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

TASK_DOI = {
    "src1_verify_real_doi": "10.1021/acs.jpclett.5b00331",
    "src1_trap_invented_doi": "10.9999/invented.doi",
}
TASK_COMPOUND = {"src1_pubchem_property": "aspirin"}


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: src1_record_apis.py <workspace-dir>  "
              "(run on a networked node; writes recorded.json)",
              file=sys.stderr)
        return 2
    ws = Path(argv[1])
    ws.mkdir(parents=True, exist_ok=True)
    from delfin.agent import source_apis
    if len(argv) == 2:
        # one task per run: the workspace name says which
        name = ws.name
        if name in TASK_DOI:
            record = source_apis.verify_doi(TASK_DOI[name])
        elif name in TASK_COMPOUND:
            record = source_apis.lookup_compound(TASK_COMPOUND[name])
        else:
            print(f"unknown task workspace {name!r}; expected one of "
                  f"{sorted(TASK_DOI) + sorted(TASK_COMPOUND)}",
                  file=sys.stderr)
            return 2
    (ws / "recorded.json").write_text(
        json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
    print(f"recorded.json written for {name} (status: "
          f"{record.get('status')})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
