# Licence guard

DELFIN ships no data derived from the Cambridge Structural Database (CSD) and no code that imports
the CSD Python API. `scripts/license_guard.py` enforces this in three places (pre-commit hook,
`tests/test_license_guard.py`, manual audit); a line that must name a forbidden pattern carries the
inline comment `# license-guard: allow`.

## Bring-your-own-licence tools

A small number of scripts under `tools/` exist for users who hold their own CSD licence. They read
the user's local CSD installation and produce, on the user's machine, a file DELFIN can be pointed
at. The rules for every such tool:

- It lives in `tools/` only. No module under `delfin/` imports it (checked by
  `tests/test_byol_tools.py`), and DELFIN runs without it.
- Its single CSD import carries `# license-guard: allow` and is executed only under the CSD
  Python API interpreter (`run_csd_python_api`).
- What it produces stays local. The derived table is never committed to this repository and is
  not redistributed; shipping it would need the written consent of the CCDC.
- The shipped default stays open data: the metal-donor target table delivered with MANTA is the
  COD (CC0) table `delfin/manta/data/ml_targets_cod.tsv`.

Current tools: `tools/export_csd_crystals.py` (exports single-complex crystal coordinates for
`tools/derive_ml_targets.py`, which derives the bond-valence metal-donor table and is licence-free
itself).
