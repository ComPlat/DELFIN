"""
CONTROL File Examples for LLM Agent

These examples teach the agent the EXACT DELFIN syntax.
NO CODE - just user-facing examples!
"""

CONTROL_EXAMPLES = {
    "basic_redox": """
# Example: Basic Redox Calculation
input_file=molecule.xyz
charge=0
solvent=MeCN
implicit_solvation_model=CPCM

# Calculate 1 oxidation and 1 reduction
calc_initial=yes
oxidation_steps=1
reduction_steps=1
method=OCCUPIER

functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
PAL=8
maxcore=4000
""",

    "multiple_redox": """
# Example: Multiple Redox Steps
input_file=complex.xyz
charge=-1
solvent=DMF
implicit_solvation_model=CPCM

# Calculate 3 reductions and 2 oxidations
calc_initial=yes
reduction_steps=1,2,3
oxidation_steps=1,2
method=OCCUPIER

functional=wB97X-V
main_basisset=def2-TZVP
metal_basisset=def2-TZVP
PAL=32
maxcore=3800
""",

    "excited_states": """
# Example: Excited States (Absorption/Emission)
input_file=molecule.xyz
charge=0
solvent=Water
implicit_solvation_model=CPCM

# No redox, just excited states
calc_initial=yes
oxidation_steps=
reduction_steps=
method=classic

# Enable absorption and emission
absorption_spec=yes
emission_spec=yes
E_00=yes
NROOTS=15

functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
PAL=8
maxcore=4000
""",

    "isc_ic": """
# Example: ISC and IC Rates
input_file=molecule.xyz
charge=0
solvent=MeCN
implicit_solvation_model=CPCM

# Redox with ESD module for ISC/IC
calc_initial=yes
reduction_steps=1
oxidation_steps=1
method=OCCUPIER

# Enable ESD module (calculates ISC/IC automatically!)
ESD_modul=yes
NROOTS=15

functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
PAL=8
maxcore=4000
""",

    "redox_and_excited": """
# Example: Redox + Excited States
input_file=complex.xyz
charge=-2
solvent=DMF
implicit_solvation_model=CPCM

# Both redox and excited states
calc_initial=yes
reduction_steps=1,2
oxidation_steps=1
method=OCCUPIER

# Also calculate absorption
absorption_spec=yes
E_00=yes
NROOTS=15

functional=wB97X-V
main_basisset=def2-TZVP
metal_basisset=def2-TZVP
PAL=32
maxcore=3800
""",

    "ground_state_only": """
# Example: Ground State Optimization Only
input_file=molecule.xyz
charge=0
solvent=MeCN
implicit_solvation_model=CPCM

# Only optimize ground state
calc_initial=yes
oxidation_steps=
reduction_steps=
method=classic

functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
PAL=8
maxcore=4000
"""
}


# Mapping user intents to examples
INTENT_TO_EXAMPLE = {
    "redox": "multiple_redox",
    "oxidation": "multiple_redox",
    "reduction": "multiple_redox",
    "excited": "excited_states",
    "absorption": "excited_states",
    "emission": "excited_states",
    "isc": "isc_ic",
    "ic": "isc_ic",
    "intersystem": "isc_ic",
    "e_00": "excited_states",
    "ground": "ground_state_only",
}


def get_relevant_example(user_goal: str) -> str:
    """
    Get relevant CONTROL example based on user's goal.

    Args:
        user_goal: User's description of what they want to calculate

    Returns:
        Relevant CONTROL file example
    """
    if not user_goal:
        return CONTROL_EXAMPLES["basic_redox"]

    user_goal_lower = user_goal.lower()

    # Check for multiple intents
    examples = []
    for keyword, example_key in INTENT_TO_EXAMPLE.items():
        if keyword in user_goal_lower:
            examples.append(CONTROL_EXAMPLES[example_key])

    if not examples:
        return CONTROL_EXAMPLES["basic_redox"]

    # If multiple intents (e.g., redox + excited), show combined example
    if len(examples) > 1 and "redox" in user_goal_lower and any(
        x in user_goal_lower for x in ["excited", "absorption", "emission", "isc", "ic"]
    ):
        return CONTROL_EXAMPLES["redox_and_excited"]

    return examples[0]


def format_example_for_prompt(user_goal: str) -> str:
    """
    Format example as part of system prompt.

    Args:
        user_goal: User's description

    Returns:
        Formatted example string for LLM
    """
    example = get_relevant_example(user_goal)
    return f"""
IMPORTANT: Here is an example of REAL DELFIN CONTROL file syntax:

{example}

KEY SYNTAX RULES:
1. reduction_steps=1,2,3  (comma-separated numbers, NO spaces!)
2. oxidation_steps=1,2    (comma-separated numbers, NO spaces!)
3. Empty steps: reduction_steps=  (just leave empty, don't write "none")
4. Yes/no values: calc_initial=yes  (lowercase!)
5. ESD_modul=yes enables ISC/IC calculations automatically

YOU MUST follow this EXACT syntax!
"""
