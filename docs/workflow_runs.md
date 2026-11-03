# Workflow runs — wave-22 RUN1 real-run evidence

Every row states what ran (or why it could not), with the source of each
number. Runs were made 2026-10-26 on the wave cluster through DELFIN's
own workflow path, small systems only (<= 20 atoms, xTB-heavy steps at
1 core, DELFIN's calc_policy).

## w21_descriptors_steric (steric_screening)

| Item | Value | Source |
|---|---|---|
| System | cisplatin PtCl2(NH3)2, 9 atoms, charge 0, mult 1 | `.gate/practice/run1_descriptors_steric/cisplatin.xyz` |
| Method | morfeus 0.8.0 (BuriedVolume r=3.5 A, ConeAngle internal method, full report) | run output; `import morfeus; __version__` |
| Basis | n/a — no electronic-structure calculation: the workflow computes geometric descriptors on the input geometry | morfeus 0.8.0 API (BuriedVolume, ConeAngle); no SCF/DFT block in the run output |
| Solvent | n/a — gas-phase / in-vacuo descriptors on the input coordinates | morfeus 0.8.0 API; no solvent model parameter exists in the workflow config |
| Resources requested/used | 1 core / 1 core, in-process (no scheduler job) | driver `SolvationWorkflow`-style probe, gate run log |
| Wall time | ~6 s | gate run log (elapsed 6.38 s incl. pytest) |
| Core-hours | ~0.002 (1 core x 6 s) | wall time above |
| Physics verdict | Results plausible for square-planar Pt(II): V_bur = 104.98 A^3 of a 179.59 A^3 sphere (58.45 %), cone angle 278.66 deg (planar ligands span most of one hemisphere; morfeus's internal method is used because libconeangle is absent — the exact-method value may differ) | `.gate/practice/run1_descriptors_steric/result.json`; warnings in the run log |
| Defect found by this run | `percent_buried_volume` carried a FRACTION (0.5845) mislabeled as percent — morfeus's deprecated attribute returns `fraction_buried_volume` unchanged (`morfeus/buried_volume.py:619`); DELFIN wrapper `delfin/analysis_tools/morfeus_wrapper.py:122` stored it under the % key, and `delfin/analysis_tools/SETUP_GUIDE.txt:31` printed it with a % suffix (both wrong before). Fix: operator build 6b06c675 (wrapper now returns `fraction_buried_volume` + a true percent x100; red control `tests/test_run1_buried_volume_unit.py` red on the old wrapper, 92 passed after) | morfeus site-packages source; run `result.json`; operator verification |

## Raw run evidence (retrievable source for the rows above)

`result.json` of the descriptors run, written by the run driver from the
three StepResults (paths shortened):

```json
[
  {"status": "SUCCESS", "data": {"percent_buried_volume": 0.5845310752888557,
    "free_volume": 74.61588395531865, "buried_volume": 104.97849607489786,
    "total_volume": 179.59437982499998, "metal_index": 0, "radius": 3.5}},
  {"status": "SUCCESS", "data": {"cone_angle": 278.66228125169874,
    "metal_index": 0}},
  {"status": "SUCCESS", "data": {"xyz_file": ".gate/practice/run1_descriptors_steric/cisplatin.xyz",
    "metal_index": 0,
    "buried_volume": {"percent_buried_volume": 0.5845310752888557,
      "free_volume": 74.61588395531865, "buried_volume": 104.97849607489786,
      "total_volume": 179.59437982499998, "metal_index": 0, "radius": 3.5},
    "cone_angle": {"cone_angle": 278.66228125169874, "metal_index": 0}}}
]
```

(The `percent_buried_volume` values above are the FRACTIONS the old
wrapper returned — the run ran before the R4 fix; the same quantity
after the fix appears as `fraction_buried_volume: 0.5845` and
`percent_buried_volume: 58.45`.)

## w21_solvation (explicit_solvation)

| Item | Value | Source |
|---|---|---|
| System prepared | formaldehyde (4 atoms) + 5 H2O = 19 atoms, charge 0, mult 1 | `.gate/practice/run1_solvation/{formaldehyde,water}.xyz` |
| Outcome | FAILED at step 1 (Packmol box), 2.9 s, no core-hours | workflow result: `errors: ["packmol box failed: Packmol solvation failed"]` |
| Cause | packmol NOT on PATH — `shutil.which("packmol")` -> None, observed twice. The adapter chain would auto-install it (`delfin/analysis_tools/packmol_wrapper.py:18-35` -> `qm_health.provide`); the package task forbids installing anything, so the run stops here | shell observation 2026-10-26 |
| Verdict | **Skipped with evidence — packmol binary absent on this cluster.** No physics numbers exist for this workflow; nothing was measured | this table |
| Level of theory | n/a — nothing ran (failed at step 1, packmol absent) | this table |

## w21_turbomole_workflow

| Item | Value | Source |
|---|---|---|
| Outcome | **Skipped — turbomole binary absent.** `delfin.installer.find("turbomole")` returned None | gate run of `tests/test_zz_probe_installer_evidence.py`, 2026-10-26 |
| Verdict | No run, no physics numbers, nothing installed | this table |
| Level of theory | n/a — nothing ran (turbomole binary absent) | this table |

## w21_mlp_screen

| Item | Value | Source |
|---|---|---|
| Outcome | **Skipped — no MLP backend installed.** ANI-2x, AIMNet2, MACE-OFF, PyTorch all absent ("No MLP backends installed") | `mlp_check` (delfin-ops), 2026-10-26; installer probe above |
| Grading basis | Wave-21 fake-run evidence plus the RUN1 correctness fixes (silent-inf best energy R1, non-int-index best entry R3, silent run_cli failure): the workflow is correct by construction, but no physics numbers from a real run exist | commits 395696da, 23932d78, d020df53 |
| Verdict | **Skipped with evidence** | this table |
| Level of theory | n/a — nothing ran (no MLP backend installed) | this table |

## What this package measured

- 1 real run (descriptors_steric), 3/3 steps successful, on 9 atoms.
- 1 physics defect found BY the real run (R4, above) — the reason real
  runs are graded, not faked.
- 3 skips, each with the probe that came back empty named here.
- Core-hours used: ~0.002 of the wave's 50/day budget.
- Not included (QS ruling, wave-23 candidate): the four w21 workflows
  do not share one return-shape contract (`w21_descriptors_steric.run`
  returns a raw StepResult list, the other three return a dict with
  ok/steps). Changes no physics; needs a coherent change across all
  four files with its own red test.
