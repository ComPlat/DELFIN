# Build energy table
> Extract Gibbs / ZPE / electronic energies across all calc folders into a CSV.

The user wants a tabular view of energies across multiple calculations.
Use the calc-search and analyze tools, then write a CSV in the agent
workspace.

1. Use `search_calcs` (the calc-search MCP tool) to identify the relevant
   calculation set.  If the user named a system or method, filter by it;
   otherwise default to the active calc directory.
2. For each match, run `/analyze energy <dir>` to get electronic / ZPE /
   Gibbs energies in Eh.
3. Write a CSV file `agent_workspace/energies.csv` with columns:
   ``folder, system, functional, basis_set, solvent, e_elec_Eh,
   e_zpe_Eh, e_gibbs_Eh, status``.
4. Report a 5-line summary in chat: total rows, average Gibbs, min/max
   Gibbs, any rows missing data.

Constraints:
- Only WRITE inside `agent_workspace/` — never modify calc/ or archive/.
- If more than 200 folders match, ask the user to narrow the filter.
- Use absolute Eh values; do NOT convert to relative kcal/mol unless asked.

Output format: short German chat summary + CSV path.
