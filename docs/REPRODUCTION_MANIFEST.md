# DELFIN Reproduction Manifest (Nature Race Checkpoint)

**Tag:** `race-2026-06-03-v1`
**Created:** 2026-06-03
**Purpose:** every pool produced during the Nature breakthrough race is reproducible from this manifest.

## Reproducibility Promise

Every pool listed below can be byte-identically reproduced given:
1. The exact DELFIN commit (visible in pool label as `<git-short-hash>-...`)
2. The exact env-flag set (documented below)
3. The exact SMILES input file (SHA-256 listed)
4. The exact GRIP library (SHA-256 listed; alternate v1/v2)
5. `PYTHONHASHSEED=0`

## Critical Artefact Checksums (SHA-256)

| File | SHA-256 |
|---|---|
| `grip_lib_v1.npz` (10 MB, 172 k Mogul fragments) | `58ed8ec7fcd13b0095021d43236921e0217208925d30d549ed631ec52504a3cf` |
| `grip_lib_v2.npz` (66 MB, 117 k COD CIFs, bonds 5.6×, torsions 0→15k brand-new) | `0254e03c7c8dcb6286dace13a1fab7d5ed20a9ce78017b804138df9e49f8cacc` |
| `scripts/grip_build_mogul_lib_v2.py` | `85633c5781faaaec74c5b3552fb7a92065191860804b87e5bb4cd744679523ee` |
| `scripts/grip_validate_lib_v2.py` | `90057a525aae8ad9dc1e599717cd6abcbf715451c529d530d11d17cdaaf9e119` |
| `pools/smiles_master_v3_plus.txt` (11 458 SMILES) | `a846f939fdedc28a0195e535d8ec85df546785930fefcdeb0f55b9e90a36b2c8` |
| `smoke500_smiles.txt` (500 SMILES, first 500 of v3_plus) | `eaa8c53c297dcb57c56982d0c6c43bacaa084736da1c6869edaa1485c0b0a99b` |

**COD CIF source for v2 lib:** `agent_workspace/quality_framework/COD_all/cif/` (370 056 CIFs)

## Production Baselines (gepushed)

### pgcorr-v3 voll-pool baseline (the stable production HEAD as of 2026-06-02)

```bash
# Reproduce: WIP-fb1ae9a-equal-n-vs-GRIP-pgcorr-v3-VOLLPOOL
git checkout 1817bbc                                           # GRIP stable baseline
export PYTHONHASHSEED=0
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
export DELFIN_FFFREE_DD_RELAX=1
export DELFIN_GRIP_LIB_PATH=/path/to/grip_lib_v1.npz
micromamba run -n delfin python \
    agent_workspace/quality_framework/scripts/pool_evaluator.py \
    agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir . --parallel 64 --timeout 300 \
    --commit-label "REPRO-pgcorr-v3-VOLLPOOL" \
    --xyz-archive agent_workspace/quality_framework/xyz_archive
```

Result (verified 2026-06-02): **80 % severe regression reduction (16 → 3)** vs prior fb1ae9a baseline; `cshm=6.85` (aggregate, hapto-artefact-inflated; honest σ-only ≈ UFF-tie).

### #131 + v2-lib combined voll-pool (race candidate)

```bash
# Reproduce: f8c9905-v2lib-construction-VOLLPOOL
git checkout f8c9905                                           # this commit's tag: race-2026-06-03-v1
export PYTHONHASHSEED=0
export DELFIN_GRIP_LIB_PATH=/path/to/grip_lib_v2.npz           # SHA above
export DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1                    # amide_vsepr + mddir_align + build_clash_gate
export DELFIN_FFFREE_GRIP=1
export DELFIN_FFFREE_POST_GRIP_ALL=1
export DELFIN_FFFREE_PURE_TRACK3=1
export DELFIN_FFFREE_DD_RELAX=1
micromamba run -n delfin python \
    agent_workspace/quality_framework/scripts/pool_evaluator.py \
    agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \
    --shadir . --parallel 64 --timeout 300 \
    --commit-label "REPRO-f8c9905-v2lib-construction-VOLLPOOL" \
    --xyz-archive agent_workspace/quality_framework/xyz_archive
```

Running 2026-06-03 (44 % progress at backup); first attempt of the "win on ALL metrics" race candidate.

## Default-OFF Determinism Doctrine

When **no env flag** is set, commit `f8c9905` HEAD produces **byte-identical** output to commit `1817bbc` (pgcorr-v3 stable baseline). All new modules:

- `amide_vsepr_template.py` — Fix 1, gated by `DELFIN_FFFREE_AMIDE_VSEPR`
- `mddir_alignment.py` — Fix 2, gated by `DELFIN_FFFREE_MDDIR_ALIGN`
- `build_time_clash_gate.py` — Fix 3, gated by `DELFIN_FFFREE_BUILD_CLASH_GATE`
- `single_bond_rotamers.py` — gated by `DELFIN_FFFREE_SINGLE_BOND_ROTAMERS`
- `conformer_dedup.py` — gated by `DELFIN_FFFREE_RMSD_DEDUP_DISABLE` (note inverted)
- `inter_ligand_clash_gate.py` — gated by `DELFIN_FFFREE_PRE_POLISH_CLASH_GATE`
- `donor_slide.py` — gated by `DELFIN_FFFREE_DONOR_SLIDE`
- `grip_mogul_lookup.py` — v2 lib via `DELFIN_GRIP_LIB_PATH` env override (no flag = default v1 search path)

Master flag `DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1` enables Fix 1/2/3 simultaneously.

## How to Rebuild GRIP v2 Library from Source

```bash
# requires the 370 056 COD CIFs at agent_workspace/quality_framework/COD_all/cif/
micromamba run -n delfin python \
    agent_workspace/quality_framework/scripts/grip_build_mogul_lib_v2.py \
    --workers 64 \
    --log /tmp/grip_lib_v2_build.log
# Wall time: ~30 min (35 s extraction + ~30 min GMM fits)
# Output: agent_workspace/quality_framework/reports/grip_lib_v2.npz (66 MB)
# Validate with grip_validate_lib_v2.py — must match SHA-256 above
```

## Backup Location

All artefacts in this manifest are also backed up at:

`/home/qmchem_max/DELFIN_NATURE_BACKUP_2026_06_03/`

Includes:
- `grip_lib_v1.npz` + `grip_lib_v2.npz` (irreplaceable binary libs)
- All `*FORENSIK*` and `*VERDICT*` markdowns for 2026-06-03
- `grip_build_mogul_lib_v2.py` + `grip_validate_lib_v2.py`
- `v2_build.log` (33 min build trace, deterministic)
- `claude_memory_<timestamp>.tar.gz` (485 KB) — full cumulative project memory
- `CHECKSUMS.sha256`

Plus daily HDD backup per `feedback_nothing_lost_doctrine`.

## Recovery Recipe

If a fresh clone needs to reproduce a pool:

1. `git clone <repo> && cd DELFIN && git checkout <commit-from-pool-label>`
2. Restore the GRIP library: `cp /backup/grip_lib_v2.npz agent_workspace/quality_framework/reports/`
3. Verify SHA-256 matches this manifest
4. Set env flags as documented in the pool's section above
5. Run `pool_evaluator.py` with the documented arguments
6. Result must match the original archive (file count + per-SMILES XYZ bytes)
