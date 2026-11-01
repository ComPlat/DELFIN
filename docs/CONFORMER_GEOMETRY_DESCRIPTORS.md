# Conformer Geometry Classification: Scientific Reference

These definitions document the implemented algorithms and are not derived
from any private ensemble. Source modules: `delfin/analysis_tools/`.

## Geometry and weighting

For atom coordinates r_i and weights w_i, W = sum(w_i),
c = sum(w_i*r_i)/W and S = sum(w_i*(r_i-c)*(r_i-c)^T)/W.
Let 0 <= t1 <= t2 <= t3 be the eigenvalues of S.
Physical ordered inertia moments I1 <= I2 <= I3 are eigenvalues of
I = W*(trace(S)*identity - S). All descriptor calculations use the
explicitly selected atoms. Mass, uniform, and heavy-atom weighting are
supported; the maximum pairwise distance is always weighting-independent.

| Descriptor | Implemented formula | Units | Source |
| --- | --- | --- | --- |
| rg | sqrt(t1+t2+t3) | Å | [T] |
| dmax | max distance over selected atom pairs | Å | Euclidean definition |
| asphericity | sum_{i<j}(ti-tj)^2 / (2*(sum ti)^2) | 1 | [T] tensor framework |
| eccentricity | sqrt(I3^2-I1^2) / I3 | 1 | [S] inertia framework |
| npr1 | I1/I3 | 1 | [S] |
| npr2 | I2/I3 | 1 | [S] |
| planarity | (t2-t1)/(t3-t1), where denominator > 0 | 1 | [S] framework; DELFIN operational ratio |
| acylindricity | (t2-t1)/(sum ti) | 1 | [T] framework; normalized here |
| kappa2 | 1.5*sum(ti-mean(t))^2/(sum ti)^2 | 1 | [T] |
| compacity | rg/dmax | 1 | DELFIN-derived operational definition |
| foldedness | 1-rg/((dmax/2)*sqrt(3/5)) | 1 | DELFIN-derived operational definition |
| plane_rms_A | root-mean-square orthogonal distance to the fragment's best-fit plane | Å | [T] principal-axis framework; plane fitting defined here |

**[S]** Sauer, W. H. B.; Schwarz, M. K. *Molecular Shape Diversity of
Combinatorial Libraries: A Prerequisite for Broad Bioactivity*,
J. Chem. Inf. Comput. Sci. **2003**, 43, 987-1003.
DOI: https://doi.org/10.1021/ci025599w

**[T]** Theodorou, D. N.; Suter, U. W. *Shape of Unperturbed Linear Polymers:
Polypropylene*, Macromolecules **1985**, 18, 1206-1214.
DOI: https://doi.org/10.1021/ma00148a028

These papers support inertia and gyration shape analysis; derived DELFIN
normalizations are *not* attributed to papers which did not define them.
The implemented kappa2 and normalized asphericity are exactly equal.
NPR1 near 0.5 does not prove sphericity: NPR2 must also be considered.
PMI-based planarity is not the same as the fragment's plane_rms_A.
Maximum molecular span (dmax) is not a measure of straightness.

## Chemistry and confidence

XYZ files have no bond orders. E/Z-like geometrical classes are heuristic.
The CIP E/Z path obtains inferred bond orders from RDKit and assigns
stereochemistry from coordinates. **Each result is a coordinate-based
reconstruction, not an independently validated assignment.** Bond-order
or stereo-perception failures yield status `unknown`, never a guessed E or Z.
The heuristic s-cis/s-trans boundaries are working conventions only.
Metal complexes and charged species require an independent chemical reference.

## Existing DELFIN infrastructure

- `delfin/tools/adapters/crest_ensemble.py` provides the existing CREST
  multi-XYZ parser (comment energies in hartree).
- `delfin/tools/adapters/orca_goat_opi.py` provides the existing GOAT
  final-ensemble reader (comment energies relative in kcal/mol).
- `delfin/workflows/contrib/w21_descriptors_steric.py` covers Sterimol,
  buried-volume and cone-angle descriptors for metal complexes; it is
  complementary rather than replaced by molecular shape classification.
- `descriptors-steric`, `conformer-search` and `ensemble-crestandcenso` skills
  remain independently discoverable.

No ORCA-generated, private or unpublished structural data is distributed.
All committed feature tests use freshly written synthetic XYZ coordinates.
