# Provenance

Both files in this directory are SYNTHETIC test fixtures written for DELFIN,
not ORCA output. Nothing here is derived from any real calculation or any
licensed-program run; no coordinates from any real system are used.

* `synth_h2o.out` -- hand-authored text in the ORCA 6.1.1 output format. Its
  first line is the required SYNTHETIC marker. It carries only the exact
  block lines DELFIN's thermochemistry readers match (Program-Version,
  SCF SETTINGS/SUCCESS, THERMOCHEMISTRY block, FINAL SINGLE POINT ENERGY,
  ORCA TERMINATED NORMALLY).
* `synth_h2o.property.json` -- hand-authored property data in OPI's property
  JSON schema (keys verified against opi.output.core._safe_get paths:
  Geometries[].Single_Point_Data, Geometries[].thermochemistry_energies[]).
  Its numbers are the same values the .out prints, so both readers are fed
  the identical quantities and the equivalence assertion is meaningful.

Values are dimensionally consistent (H = U + kB*T, G = H - T*S with the
entropy term T*S in Hartree) and were checked against a real converged
Opt+Freq H2O run's value patterns in scratch/ during development; the fixture
itself contains no measured data from any real structure.
