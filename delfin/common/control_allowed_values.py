"""Allowed values per CONTROL key, when the FieldSpec carries none.

``CONTROL_FIELD_SPECS.allowed`` is the source of truth; this table fills
the gaps with one line per key stating what its coercer accepts (every
line derived from the coercer code, none from memory).  Readers --
``help_lookup``, the Help tab, the header key search -- resolve as
``spec.allowed or ALLOWED_BY_KEY.get(spec.name, "")``.
"""

ALLOWED_BY_KEY: dict[str, str] = {
    # --- Block 1: general + GUPPY + MANTA construction (25) -------------
    "NAME": "free text",
    "SMILES": "free text (a SMILES string)",
    "charge": "an integer charge, e.g. -2, 0, +3",
    "multiplicity_global_opt": "an integer multiplicity (1, 2, 3, ...)",
    "PAL": "a positive integer (total CPU cores)",
    "GUPPY_RUNS": "a positive integer (search runs)",
    "GUPPY_GOAT": "an integer between 0 and 10",
    "GUPPY_PARALLEL_JOBS": "a positive integer (jobs in parallel)",
    "GUPPY_SEED": "an integer (random seed)",
    "GUPPY_START_STRATEGY": "one of: isomers, isomers+random, full (empty = isomers)",
    "GUPPY_MAX_ISOMERS": "a positive integer",
    "GUPPY_RMSD_CUTOFF": "a number >= 0 (Angstrom)",
    "GUPPY_ENERGY_WINDOW_KCAL": "a number >= 0 (kcal/mol)",
    "MANTA_SEEDS": "a positive integer, or empty = builder default",
    "MANTA_NUM_CONFS": "a positive integer, or empty = builder default",
    "MANTA_CONSTRUCTION": "one of: champion, builder, default (empty = default)",
    "MANTA_MAX_ISOMERS": "an integer",
    "MANTA_HAPTO": "one of: auto, on, off (empty = auto)",
    "MANTA_UFF": "yes or no (true/1/on also read as yes; empty = no)",
    "MANTA_DETERMINISTIC": "yes or no (true/1/on also read as yes; empty = no)",
    "MANTA_COLLAPSE_VARIANTS": "yes or no (true/1/on also read as yes; empty = no)",
    "MANTA_ENV": "free text (KEY=VALUE pairs, semicolon separated)",
    "MANTA_CLEAN_GATE": "yes or no (true/1/on also read as yes; empty = no)",
    "MANTA_TOPOLOGY_GATE": "yes or no (true/1/on also read as yes; empty = no)",
    "MANTA_DEDUP": "yes or no (true/1/on also read as yes; empty = no)",
}
