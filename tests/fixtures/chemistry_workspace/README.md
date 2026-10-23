# Chemistry benchmark fixture workspace

The setup script `chem_start_geometries.py` seeds `chem/<molecule>/`
folders with start geometries and a README (SMILES, charge) when a
chemistry task runs. This directory exists so the workspace path
resolves before the setup script is launched (see `workspace_for` in
`delfin/agent/benchmark_runner.py`).

Reference values for the tolerance bands in `tasks_chem.yaml` were
measured with xtb 6.7.1 GFN2 (`--alpb water`, opt+hess) from exactly
these start geometries.
