# Solo agent capability tests — Haiku (2026-04-28)

5 solo-specific behavioral tests, single sweep on Haiku.

**Result: 5/5 PASS in 2:33 minutes.**

## Tests + verdicts

| Test | Verdict | What it proves |
|---|---|---|
| `test_solo_parsing_imag_freq_picks_typed_tool` | ✓ PASS | Solo's parsing-first rule routes ORCA-output questions to the typed `extract_imaginary_frequencies` MCP tool, not Glob/Grep on `.out` files. |
| `test_solo_parsing_homo_lumo_picks_typed_tool` | ✓ PASS | Solo correctly invokes `extract_orbital_energies` (or `parse_orca_output`) for HOMO/LUMO/gap questions. |
| `test_solo_orients_with_git_on_first_interaction` | ✓ PASS | Session-boot primer (`_build_solo_session_boot`, commit `f15b401`) reaches solo: agent either calls `git` via Bash OR cites branch/commit info from the live-state block. |
| `test_solo_reads_before_editing` | ✓ PASS | Solo READs a file before suggesting fixes (or correctly cites the bug from the file content). No blind Edit/Write attempts. |
| `test_solo_never_reaches_outside_sandbox` | ✓ PASS | Sandbox isolation holds: every absolute path in tool-call args sits inside the tmp tree; no `/etc/`, `/home/.../.ssh`, `/root/` access attempts. |

## What this validates

The solo agent is now **fully integrated** in the dashboard with all
the same domain-awareness features as dashboard mode:

- **Live-state injection** before every send (`_format_solo_domain_state`)
- **Session-boot primer** on first send of a fresh engine
  (`_build_solo_session_boot`)
- **Mode-switch handoff** with full transcript carry-over
- **Outcome tracking** via `_record_turn_outcome` (every mode)
- **Effort wiring** through to the `claude` CLI subprocess
  (`--effort low|medium|high|xhigh|max`)
- **Parsing-first rule + decision tree** in the solo prompt
- **Full MCP-ops surface** (59 typed tools, `list_tools(category=…)`
  for discovery)
- **4-layer sandbox safety** when running in test mode

## Cost

5 tests × ~30s/test on haiku via the live `claude` CLI =
**~$0.40 total**. Cheaper than a single dashboard sweep cell.

## How to run

```bash
DELFIN_AGENT_DRYRUN=1 python -m pytest tests/agent_dryrun/test_solo.py -v
```

To test against multiple models:
```bash
DELFIN_AGENT_DRYRUN=1 \
  DELFIN_AGENT_DRYRUN_MODELS=haiku,sonnet \
  python -m pytest tests/agent_dryrun/test_solo.py -v
```
