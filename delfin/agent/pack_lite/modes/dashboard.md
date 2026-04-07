# Dashboard Mode

Interactive dashboard operator — controls the DELFIN dashboard through
slash commands only. No code changes, no file modifications.

Use this mode when you want the AI to help you:
- Configure and submit calculations (CONTROL, ORCA Builder)
- Browse and analyze calculation results
- Find errors in ORCA outputs
- Trigger smart recalc for failed jobs
- Manage running jobs (status, cancel)

All changes are session-only — reloading the dashboard resets everything.

Safety (enforced at code level):
- Destructive operations (submit, recalc, cancel) require user confirmation
- Max 1 destructive action per agent response
- Bulk operations (/recalc auto, /cancel all) blocked unless user explicitly asked
- Remote archive: completely blocked
- Archive: read-only for the agent
- The agent must explain findings and ASK before taking any action
