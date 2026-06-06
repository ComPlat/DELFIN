# DELFIN Legal & Licensing Notes

## Overview

DELFIN is a Force-Field-free SMILES → 3D pipeline for transition-metal complexes (TMCs).
The pipeline ships with **no proprietary data**: all reference distributions used at
runtime can be derived from public-domain sources (COD, Crystallography Open Database,
CC0). Optional integration with the Cambridge Structural Database (CSD) is supported
for users who hold an academic or commercial CCDC license, but is **not required**
for any runtime path.

## Data sources

DELFIN can consume two parallel fragment libraries (`.npz` files):

### 1. COD-derived library (default, adopter-safe)

* **Source**: [Crystallography Open Database (COD)](https://www.crystallography.net/cod/)
* **License**: CC0 / Public Domain
* **File**: `grip_lib_v6_cod.npz` (~50 MB, 259k fragments)
* **Status**: ✅ Fully redistributable; ✅ No CCDC license needed
* **Used by adopters when**: `DELFIN_FFFREE_GRIP_LIB_COD_PATH` points to the file
  OR the file is found in standard locations (see below).

### 2. CCDC-derived library (academic/commercial CCDC users only)

* **Source**: Cambridge Structural Database (CSD), via official `ccdc.search` Python API
* **License**: CCDC Academic / Commercial License (held by user)
* **Files**: `grip_lib_v1.npz` ... `grip_lib_v5.npz`
* **Status**: ⚠️ Built by us from CCDC API using our institutional academic license,
  for our internal validation only. **Not shipped with DELFIN.**
* **Used when**: `DELFIN_GRIP_LIB_PATH` explicitly points to a user-owned file.

## Library discovery order

When DELFIN starts, it searches for a runtime library in this order:

1. `$DELFIN_FFFREE_GRIP_LIB_COD_PATH` (adopter-preferred, CC0)
2. `$DELFIN_GRIP_LIB_PATH` (user override, may be CCDC-licensed v5 if held)
3. Package data directory: `<delfin_root>/data/grip_lib_v6_cod.npz`
4. XDG data dir: `$XDG_DATA_HOME/delfin/grip_lib_v6_cod.npz`
5. User-local: `~/.local/share/delfin/grip_lib_v6_cod.npz`
6. System data dir: `/usr/share/delfin/grip_lib_v6_cod.npz`

If no library is found, DELFIN emits a clear error message with instructions.

**Key guarantee**: The default runtime path requires only the COD-derived library.
No CCDC license is needed to run DELFIN with v6.

## Runtime CCDC dependency: NONE

DELFIN's production code does **not import the `ccdc` Python library** at runtime.
Verification:

```
$ grep -r "import ccdc" delfin/      # 0 hits
$ grep -r "from ccdc" delfin/        # 0 hits
```

References to "CCDC" and "Mogul" in DELFIN source code are exclusively in
docstrings and comments documenting methodology.

## Methodology references

* **Mogul algorithm** — Bruno, I. J. et al. *J. Appl. Cryst.* **35**, 646 (2002).
  Publicly described; DELFIN's GRIP scorer reimplements the principle of
  Mahalanobis-distance fragment scoring against an empirical distribution library.
* **COD** — Gražulis, S. et al. *Nucleic Acids Res.* **40**, D420 (2012).
  CC0 public-domain database, no usage restrictions.
* **CCDC** — Groom, C. R. et al. *Acta Cryst.* **B72**, 171 (2016).
  Used only by license-holders for internal validation; not redistributed.

## Statistical aggregates as transformative works

The `.npz` fragment libraries contain only **statistical aggregates** (mean μ,
standard deviation σ, sample count n) of bond / angle / improper / torsion
distributions extracted per fragment type. They are transformative works
constructed from the underlying source databases (CSD for v1–v5, COD for v6).

The CC0-licensed COD-derived v6 library is freely redistributable as a derivative
work of public-domain source data. The CCDC-derived v1–v5 libraries are
constructed under our institutional CSD academic license for internal validation
and are not shipped with DELFIN.

## What adopters can do

✅ Use DELFIN with the COD-derived v6 library — no CCDC license required
✅ Modify, extend, redistribute the DELFIN source code under its own license (see LICENSE)
✅ Rebuild the COD library from public COD CIF files using `scripts/grip_build_cod_lib.py`
✅ Validate their output against any reference data they have access to

## What adopters cannot do without their own CCDC license

❌ Download or use `grip_lib_v1.npz` through `grip_lib_v5.npz` (CCDC-derived)
❌ Use the `ccdc.search` / `ccdc.conformer` Python API
❌ Validate against the proprietary CSD reference set

## License of DELFIN itself

DELFIN code: [see LICENSE](LICENSE) (placeholder — final license TBD before public release).

The CC0 v6 COD library is distributed separately under its own CC0 dedication.

## Questions

For licensing questions about the CCDC integration paths, consult:
- CCDC: https://www.ccdc.cam.ac.uk/Community/csd-community/licenceconditions/
- COD: https://www.crystallography.net/cod/ (CC0)
