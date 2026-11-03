---
title: The fluoride plus methyl fluoride barrier is gas-phase submerged
topic: chemistry-methods
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:10.1063/1.4921857
review_by: 2027-12-31
---
Situation: An SN2 benchmark of fluoride attacking methyl fluoride places the transition state below the reactants' energy, and the run is suspected of an error because a barrier "should" be positive.
Lesson: For F- + CH3F in the gas phase the central barrier is submerged: the transition state sits below the separated reactants by about 0.45 kcal/mol (reported uncertainty about 0.10 kcal/mol; DOI 10.1063/1.4921857). A negative central barrier for this benchmark is the expected result, not a broken setup — energies carry the method, basis and solvation of the comparison.
Why: Solvation and counterions lift the reactants differently from the transition state; a submerged gas-phase barrier is standard for this reaction, so the sign is a physics check, not a convergence or setup check.
