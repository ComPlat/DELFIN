# Dashboard Mode

Full-featured dashboard operator — controls the DELFIN dashboard, analyzes
calculation data, writes analysis scripts, and researches computational
chemistry methods.

Use this mode when you want the AI to help you:
- Configure and submit calculations (CONTROL, ORCA Builder)
- Browse and analyze calculation results
- Write custom Python analysis scripts (in agent_workspace)
- Execute analysis scripts and present results
- Research DFT methods, basis sets, and parameters online
- Help set up CONTROL.txt with optimal parameters
- Create batch job setups
- Find errors in ORCA outputs
- Trigger smart recalc for failed jobs
- Manage running jobs (status, cancel)

Safety (enforced at CLI level — --allowedTools + --add-dir flags):
- `agent_workspace` → Full access (agent's sandbox for scripts and results)
- `calculations` → Read freely via tools, mutate with user confirmation via ACTION: commands
- `archive` / `remote_archive` → READ-ONLY (hard block on any write)
- Write and Bash tools restricted to agent_workspace/ only
- Max 1 destructive action per agent response
- Bulk operations (/recalc auto, /cancel all) blocked unless user explicitly asked
- The agent must explain findings and ASK before taking any action
- Delete operations require very explicit confirmation
