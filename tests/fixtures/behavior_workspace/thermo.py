"""Toy xtb thermochemistry helpers — benchmark fixture, safe to edit.

Not imported by the DELFIN package. Used by the behavioural benchmark's
`scout` task: it carries a deliberate UNIT BUG in ``free_energy`` that the
agent must read the file to find before editing.
"""

R_KCAL = 0.0019872041  # gas constant, kcal/(mol·K)


def free_energy(enthalpy_kcal, entropy_cal, temperature=298.15):
    """Gibbs free energy G = H - T·S.

    ``enthalpy_kcal`` is in kcal/mol and ``entropy_cal`` is in
    cal/(mol·K) (the usual xtb/Gaussian convention).

    BUG: the entropy term is subtracted directly, but it is in
    cal/(mol·K) while enthalpy is in kcal/mol — it must be divided by
    1000 to convert cal → kcal before the subtraction.
    """
    return enthalpy_kcal - temperature * entropy_cal


if __name__ == "__main__":
    # Toluene-ish toy numbers; correct G should be ≈ -18.9 kcal/mol.
    print(free_energy(10.0, 97.0))
