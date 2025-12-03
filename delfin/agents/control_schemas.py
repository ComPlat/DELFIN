"""
JSON Schemas for DELFIN CONTROL file parameters

These schemas define the EXACT allowed values for each parameter.
LLMs CANNOT hallucinate values outside these schemas.
"""

# Allowed DFT functionals
FUNCTIONALS = [
    "PBE0", "B3LYP", "WB97X", "WB97X-V", "WB97X-D3", "WB97M-V",
    "TPSS", "TPSSh", "M06", "M06-2X", "M06-L", "M06-HF",
    "B97M-V", "B97-3c", "r2SCAN", "r2SCAN-3c", "ωB97X-V"
]

# Basis sets
MAIN_BASISSETS = [
    "def2-SVP", "def2-SVPD", "def2-TZVP", "def2-TZVPP", "def2-TZVPD",
    "def2-QZVP", "def2-QZVPP", "cc-pVDZ", "cc-pVTZ", "cc-pVQZ",
    "ma-def2-SVP", "ma-def2-TZVP", "ma-def2-TZVPP"
]

METAL_BASISSETS = MAIN_BASISSETS + [
    "SARC-ZORA-TZVP", "SARC-ZORA-TZVPP", "SARC-DKH-TZVP"
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
    "Toluene", "Benzene", "Hexane", "Methanol", "Ethanol", "DiethylEther"
]

# Methods
METHODS = ["classic", "manually", "OCCUPIER"]

# XTB methods
XTB_METHODS = ["GFN1-XTB", "GFN2-XTB", "GFN-FF"]


# Main CONTROL schema
CONTROL_SCHEMA = {
    "type": "object",
    "properties": {
        # Basic properties
        "charge": {
            "type": "integer",
            "minimum": -10,
            "maximum": 10,
            "description": "Total molecular charge"
        },
        "multiplicity": {
            "type": "integer",
            "minimum": 1,
            "maximum": 9,
            "description": "Spin multiplicity (2S+1)"
        },

        # Level of theory
        "functional": {
            "type": "string",
            "enum": FUNCTIONALS,
            "description": "DFT functional"
        },
        "main_basisset": {
            "type": "string",
            "enum": MAIN_BASISSETS,
            "description": "Main basis set for light atoms"
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
            "description": "RI approximation for Coulomb/exchange"
        },
        "relativity": {
            "type": "string",
            "enum": RELATIVITY,
            "description": "Relativistic treatment"
        },

        # Solvation
        "implicit_solvation_model": {
            "type": "string",
            "enum": SOLVATION_MODELS,
            "description": "Implicit solvation model"
        },
        "solvent": {
            "type": "string",
            "enum": SOLVENTS,
            "description": "Solvent name"
        },

        # Method
        "method": {
            "type": "string",
            "enum": METHODS,
            "description": "DELFIN calculation method"
        },

        # Redox steps
        "calc_initial": {
            "type": "string",
            "enum": ["yes", "no"],
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
            "description": "Reduction steps (e.g., '1,2,3')"
        },

        # Resources
        "PAL": {
            "type": "integer",
            "minimum": 1,
            "maximum": 512,
            "description": "Number of CPU cores for ORCA"
        },
        "maxcore": {
            "type": "integer",
            "minimum": 100,
            "maximum": 100000,
            "description": "Memory per core in MB"
        },

        # Geometry optimization
        "XTB_OPT": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Pre-optimize with XTB"
        },
        "XTB_GOAT": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Use XTB GOAT optimizer"
        },
        "CREST": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Use CREST conformer search"
        },

        # Excited states
        "absorption_spec": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Calculate absorption spectrum"
        },
        "emission_spec": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Calculate emission spectrum"
        },
        "E_00": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Calculate E_00 energy"
        },
        "NROOTS": {
            "type": "integer",
            "minimum": 1,
            "maximum": 50,
            "description": "Number of excited states"
        },

        # ESD module
        "ESD_modul": {
            "type": "string",
            "enum": ["yes", "no"],
            "description": "Enable Excited State Dynamics module (ISC/IC)"
        }
    },
    "required": ["charge", "functional", "main_basisset", "method"]
}


# Simplified schema for step-by-step questioning
STEP_SCHEMAS = {
    "basic_info": {
        "type": "object",
        "properties": {
            "charge": CONTROL_SCHEMA["properties"]["charge"],
            "multiplicity": CONTROL_SCHEMA["properties"]["multiplicity"],
        },
        "required": ["charge", "multiplicity"]
    },
    "level_of_theory": {
        "type": "object",
        "properties": {
            "functional": CONTROL_SCHEMA["properties"]["functional"],
            "main_basisset": CONTROL_SCHEMA["properties"]["main_basisset"],
            "metal_basisset": CONTROL_SCHEMA["properties"]["metal_basisset"],
        },
        "required": ["functional", "main_basisset", "metal_basisset"]
    },
    "solvation": {
        "type": "object",
        "properties": {
            "implicit_solvation_model": CONTROL_SCHEMA["properties"]["implicit_solvation_model"],
            "solvent": CONTROL_SCHEMA["properties"]["solvent"],
        },
        "required": ["implicit_solvation_model", "solvent"]
    },
    "redox": {
        "type": "object",
        "properties": {
            "method": CONTROL_SCHEMA["properties"]["method"],
            "calc_initial": CONTROL_SCHEMA["properties"]["calc_initial"],
            "oxidation_steps": CONTROL_SCHEMA["properties"]["oxidation_steps"],
            "reduction_steps": CONTROL_SCHEMA["properties"]["reduction_steps"],
        },
        "required": ["method"]
    },
    "resources": {
        "type": "object",
        "properties": {
            "PAL": CONTROL_SCHEMA["properties"]["PAL"],
            "maxcore": CONTROL_SCHEMA["properties"]["maxcore"],
        },
        "required": ["PAL", "maxcore"]
    }
}
