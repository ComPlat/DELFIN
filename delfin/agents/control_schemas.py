"""
JSON Schemas for DELFIN CONTROL file parameters (Complete Version)

These schemas define the EXACT allowed values for each parameter.
LLMs CANNOT hallucinate values outside these schemas.
"""

# Allowed DFT functionals
FUNCTIONALS = [
    "PBE0", "B3LYP", "WB97X", "WB97X-V", "WB97X-D3", "WB97M-V",
    "TPSS", "TPSSh", "M06", "M06-2X", "M06-L", "M06-HF",
    "B97M-V", "B97-3c", "r2SCAN", "r2SCAN-3c", "ωB97X-V", "wB97X-V"
]

# Basis sets
MAIN_BASISSETS = [
    "def2-SVP", "def2-SVPD", "def2-TZVP", "def2-TZVPP", "def2-TZVPD",
    "def2-QZVP", "def2-QZVPP", "cc-pVDZ", "cc-pVTZ", "cc-pVQZ",
    "ma-def2-SVP", "ma-def2-TZVP", "ma-def2-TZVPP",
    "ZORA-def2-SVP", "ZORA-def2-TZVP", "ZORA-def2-QZVP"
]

METAL_BASISSETS = MAIN_BASISSETS + [
    "SARC-ZORA-TZVP", "SARC-ZORA-TZVPP", "SARC-DKH-TZVP",
    "SARC-ZORA-SVP"
]

# Dispersion corrections
DISP_CORR = ["", "D3", "D3BJ", "D4"]

# RI approximations
RI_JKX = ["", "RIJCOSX", "RIJK", "RI-JK"]

# Relativity
RELATIVITY = ["none", "ZORA", "DKH", "X2C"]

# Solvation models
SOLVATION_MODELS = ["CPCM", "SMD", "ALPB"]

# Common solvents
SOLVENTS = [
    "Water", "MeCN", "DMF", "DMSO", "Acetone", "THF", "DCM", "Chloroform",
    "Toluene", "Benzene", "Hexane", "Methanol", "Ethanol", "DiethylEther",
    "Dichloromethane", "ChloroForm", "CCl4"
]

# Methods
METHODS = ["classic", "manually", "OCCUPIER"]

# XTB methods
XTB_METHODS = ["GFN1-XTB", "GFN2-XTB", "XTB2", "GFN-FF"]

# OCCUPIER trees
OCCUPIER_TREES = ["deep", "deep2", "deep3", "flat", "own"]

# OCCUPIER methods
OCCUPIER_METHODS = ["auto", "manually"]

# OCCUPIER selection
OCCUPIER_SELECTIONS = ["tolerance", "truncation", "rounding"]

# Yes/No options
YES_NO = ["yes", "no"]

# Boolean options
TRUE_FALSE = ["TRUE", "FALSE", "true", "false"]

# Excitation types
EXCITATION_TYPES = ["s", "t", "s,t"]

# S1 optimization methods
S1_OPT_METHODS = ["TDDFT", "deltaSCF"]

# ESD modus
ESD_MODUS = ["TDDFT", "deltaSCF"]

# Geometry optimization
GEOM_OPT = ["OPT", "COPT", "TIGHTOPT"]

# Frequency types
FREQ_TYPES = ["FREQ", "NUMFREQ", "ANFREQ"]

# ============================================================================
# Step-by-step schemas for interactive workflow
# ============================================================================

STEP_SCHEMAS = {
    # Step 1: Charge
    "charge": {
        "type": "object",
        "properties": {
            "charge": {
                "type": "integer",
                "minimum": -5,
                "maximum": 5,
                "description": "Total molecular charge"
            }
        },
        "required": ["charge"]
    },

    # Step 2: Solvent
    "solvent": {
        "type": "object",
        "properties": {
            "implicit_solvation_model": {
                "type": "string",
                "enum": SOLVATION_MODELS,
                "description": "Solvation model"
            },
            "solvent": {
                "type": "string",
                "enum": SOLVENTS,
                "description": "Solvent name"
            }
        },
        "required": ["implicit_solvation_model", "solvent"]
    },

    # Step 3: Structure quality / XTB optimization
    "structure_quality": {
        "type": "object",
        "properties": {
            "is_approximation": {
                "type": "string",
                "enum": ["yes", "no"],
                "description": "Is this structure an approximation?"
            },
            "XTB_OPT": {
                "type": "string",
                "enum": YES_NO,
                "description": "Pre-optimize with xTB"
            },
            "XTB_GOAT": {
                "type": "string",
                "enum": YES_NO,
                "description": "Use xTB GOAT optimizer"
            }
        },
        "required": ["is_approximation"]
    },

    # Step 4: Imaginary frequencies
    "imag": {
        "type": "object",
        "properties": {
            "IMAG": {
                "type": "string",
                "enum": YES_NO,
                "description": "Check and eliminate imaginary frequencies"
            },
            "IMAG_option": {
                "type": "string",
                "enum": ["1", "2", "3"],
                "description": "IMAG elimination method (2=recommended)"
            }
        },
        "required": ["IMAG"]
    },

    # Step 5: Redox potentials
    "redox_setup": {
        "type": "object",
        "properties": {
            "calc_initial": {
                "type": "string",
                "enum": YES_NO,
                "description": "Calculate initial state"
            },
            "oxidation_steps": {
                "type": "string",
                "pattern": "^[0-9,]*$",
                "description": "Oxidation steps (e.g., '1,2,3')"
            },
            "reduction_steps": {
                "type": "string",
                "pattern": "^[0-9,]*$",
                "description": "Reduction steps (e.g., '1,2')"
            }
        },
        "required": ["calc_initial"]
    },

    # Step 6: Method selection
    "method": {
        "type": "object",
        "properties": {
            "system_type": {
                "type": "string",
                "enum": ["transition_metal_complex", "tadf_organic", "organic_closed_shell", "radical", "other"],
                "description": "What type of system are you calculating?"
            },
            "method": {
                "type": "string",
                "enum": METHODS,
                "description": "Calculation method"
            }
        },
        "required": ["method"]
    },

    # Step 7: MANUALLY configuration (if method=manually)
    "manually_config": {
        "type": "object",
        "properties": {
            "multiplicity_0": {
                "type": "string",
                "description": "Initial state multiplicity"
            },
            "additions_0": {
                "type": "string",
                "description": "Additional ORCA keywords for initial state"
            },
            "multiplicity_ox1": {
                "type": "string",
                "description": "Oxidation +1 multiplicity"
            },
            "additions_ox1": {
                "type": "string",
                "description": "Additional ORCA keywords for ox1"
            },
            "multiplicity_ox2": {
                "type": "string",
                "description": "Oxidation +2 multiplicity"
            },
            "additions_ox2": {
                "type": "string",
                "description": "Additional ORCA keywords for ox2"
            },
            "multiplicity_ox3": {
                "type": "string",
                "description": "Oxidation +3 multiplicity"
            },
            "additions_ox3": {
                "type": "string",
                "description": "Additional ORCA keywords for ox3"
            },
            "multiplicity_red1": {
                "type": "string",
                "description": "Reduction -1 multiplicity"
            },
            "additions_red1": {
                "type": "string",
                "description": "Additional ORCA keywords for red1"
            },
            "multiplicity_red2": {
                "type": "string",
                "description": "Reduction -2 multiplicity"
            },
            "additions_red2": {
                "type": "string",
                "description": "Additional ORCA keywords for red2"
            },
            "multiplicity_red3": {
                "type": "string",
                "description": "Reduction -3 multiplicity"
            },
            "additions_red3": {
                "type": "string",
                "description": "Additional ORCA keywords for red3"
            },
            "additions_TDDFT": {
                "type": "string",
                "description": "Additional keywords for TDDFT"
            },
            "additions_T1": {
                "type": "string",
                "description": "Additional keywords for T1"
            },
            "additions_S1": {
                "type": "string",
                "description": "Additional keywords for S1"
            }
        },
        "required": []
    },

    # Step 8: Excited states / E_00
    "excited_states": {
        "type": "object",
        "properties": {
            "calculate_excited_states": {
                "type": "string",
                "enum": YES_NO,
                "description": "Calculate E_00 or excited state redox potentials?"
            },
            "E_00": {
                "type": "string",
                "enum": YES_NO,
                "description": "Calculate E_00 energies"
            },
            "excitation": {
                "type": "string",
                "enum": EXCITATION_TYPES,
                "description": "Singlet (s), triplet (t), or both (s,t)"
            },
            "absorption_spec": {
                "type": "string",
                "enum": YES_NO,
                "description": "Calculate absorption spectrum"
            },
            "emission_spec": {
                "type": "string",
                "enum": YES_NO,
                "description": "Calculate emission spectrum"
            }
        },
        "required": ["calculate_excited_states"]
    },

    # Step 9: Level of theory
    "level_of_theory": {
        "type": "object",
        "properties": {
            "functional": {
                "type": "string",
                "enum": FUNCTIONALS,
                "description": "DFT functional"
            },
            "main_basisset": {
                "type": "string",
                "enum": MAIN_BASISSETS,
                "description": "Basis set for light atoms"
            },
            "metal_basisset": {
                "type": "string",
                "enum": METAL_BASISSETS,
                "description": "Basis set for metal atoms"
            },
            "disp_corr": {
                "type": "string",
                "enum": DISP_CORR,
                "description": "Dispersion correction"
            },
            "ri_jkx": {
                "type": "string",
                "enum": RI_JKX,
                "description": "RI approximation"
            },
            "relativity": {
                "type": "string",
                "enum": RELATIVITY,
                "description": "Relativistic treatment"
            }
        },
        "required": ["functional", "main_basisset"]
    },

    # Step 10: Additional analysis
    "analysis": {
        "type": "object",
        "properties": {
            "print_MOs": {
                "type": "string",
                "enum": YES_NO,
                "description": "Generate MO diagrams"
            },
            "print_Loewdin_population_analysis": {
                "type": "string",
                "enum": YES_NO,
                "description": "Print Loewdin population analysis"
            }
        },
        "required": []
    },

    # Step 11: Computational resources
    "resources": {
        "type": "object",
        "properties": {
            "PAL": {
                "type": "integer",
                "minimum": 1,
                "maximum": 512,
                "description": "Number of CPU cores"
            },
            "maxcore": {
                "type": "integer",
                "minimum": 1000,
                "maximum": 100000,
                "description": "Memory per core in MB"
            },
            "parallel_workflows": {
                "type": "string",
                "enum": YES_NO,
                "description": "Run multiple calculations in parallel"
            },
            "pal_jobs": {
                "type": "integer",
                "minimum": 1,
                "maximum": 16,
                "description": "Number of parallel jobs"
            }
        },
        "required": ["PAL", "maxcore"]
    },

    # Step 12: Auto recovery
    "auto_recovery": {
        "type": "object",
        "properties": {
            "enable_auto_recovery": {
                "type": "string",
                "enum": YES_NO,
                "description": "Enable automatic error recovery"
            },
            "enable_job_timeouts": {
                "type": "string",
                "enum": YES_NO,
                "description": "Enable job timeouts"
            }
        },
        "required": []
    }
}
