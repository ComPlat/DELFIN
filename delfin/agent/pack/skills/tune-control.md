# Tune CONTROL parameters
> Read the active CONTROL.txt and propose targeted optimizations.

Analyze the active CONTROL.txt for the system the user is preparing and
suggest concrete improvements.  Follow this protocol:

1. Run `/control show` to read the current CONTROL content.
2. Identify the system class (organic, transition-metal complex, lanthanide,
   ESD/UV-Vis, NMR, redox, thermochemistry).  Look at SMILES if visible,
   otherwise ask the user one short question about the target system.
3. Use `search_docs` (the doc-search MCP tool) to confirm best-practice
   functional / basis set / dispersion / solvation choices for that system
   class — do NOT rely on memory.
4. Propose only the keys that should change.  For each: state the new value,
   the old value, and a one-line justification with a concrete reference
   (ORCA manual section or benchmark paper) when search_docs returned one.
5. Stop.  Do NOT apply the changes — present them as a list of
   `/control key <key> <value>` commands the user can execute.

Edge cases:
- If the metal is 4d/5d, ALWAYS check whether `relativity` (ZORA / X2C) is
  set; flag if missing.
- If `freq_type=analytical` but the functional/basis combination doesn't
  support it, flag and propose `numerical`.
- If `solvent` is set but `solvation_model` is missing or vice versa, flag.

Output format: a short bulleted list, German prose, exact numerical values.
