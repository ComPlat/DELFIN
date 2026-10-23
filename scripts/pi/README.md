# DELFIN agent box on a Raspberry Pi 5

Turns a Raspberry Pi 5 (aarch64) into a dedicated machine for developing and
testing the DELFIN agent. The Pi itself is the sandbox: the agent runs as its
own user with passwordless `sudo` and `DELFIN_AGENT_SANDBOX=off`.

QM binaries (ORCA, xtb, CREST, ...) are **not** installed. Everything else is:
DELFIN with the extras `agent docs dev analysis mlp-full ai-foundation
ai-screening ai-complex csp`, the Claude Code and Codex CLIs (agent
`backend=cli`), and the dashboard as a loopback-only service.

Not included: `ai-retro` (aizynthfinder pins `rdkit<2023`, which conflicts with
DELFIN's `rdkit>=2023.3` on every platform). PyMOL has no aarch64 distribution
on PyPI and is built from source (`--no-deps`, it runs on DELFIN's numpy<2);
if that build fails only the PyMOL images in DOCX/ESP reports are missing.
torch is installed as the CPU build (`--torch-backend cpu`).

Known limits, not Pi-specific: `schnetpack` does not import with torch 2.13
(`T_co` was removed from torch), `architector` needs `xtb` (QM, out of scope).
`uv pip check` reports pymol's numpy>=2 metadata and the `pymol` vs
`pymol-open-source` name; both are expected. The deployed tree has no `.git`,
so `test_pre_refactor_api_preserved` (needs `git show`) fails under `--test`;
clone the repo into `~/delfin/work` when the agent should work on git history. First boot takes ~10 min on a
Pi 5 with a fast connection.

## Files

| Script | Runs on | Purpose |
| --- | --- | --- |
| `make_pi_card.ps1` | Windows | Flash Raspberry Pi OS and put DELFIN + first-boot setup on the card |
| `deploy_pi.sh` | PC (Git Bash / WSL / Linux / macOS) | Push the current working tree to a running Pi, roll back, run tests |
| `pi_bootstrap.sh` | Pi (root) | System setup: packages, agent user, uv + Python 3.11 venv, CLIs, dashboard service |
| `pi_install.sh` | Pi (agent user) | Unpack a source tarball, install deps when `pyproject.toml` changed, keep `src.prev` for rollback |
| `pi_firstboot.sh` | Pi (cloud-init) | First boot: `pi_bootstrap.sh` + `pi_install.sh` from the card |

## New card (Windows)

Requirements: an SSH key (`~/.ssh/id_ed25519.pub`), git, a card or SSD in a
USB reader. Find its disk number with `Get-Disk`.

```powershell
powershell -ExecutionPolicy Bypass -File scripts\pi\make_pi_card.ps1 -DiskNumber <N>
```

Options: `-Hostname` (default `delfin-pi`), `-LoginUser` (`pi`), `-Timezone`
(`Europe/Berlin`), `-KeyboardLayout` (`de`), `-PubKeyFile`, `-SkipFlash` (only
refresh DELFIN + cloud-init on an already flashed card).

The script refuses system/boot disks and anything that is not USB/SD, and asks
you to type `YES` before erasing.

Then: card into the Pi, LAN cable, power on. The first boot installs
everything (20-40 min, needs internet):

```bash
ssh pi@delfin-pi.local tail -f /var/log/delfin-firstboot.log
```

## Accounts and access

- `pi`: desktop autologin, SSH by key only, password locked. Set one for local
  logins with `sudo passwd pi`.
- `delfin`: agent user, passwordless `sudo`, SSH by the same key.
- API keys are entered on the Pi, never on the card:
  `delfin-agent credentials set ANTHROPIC_API_KEY` (or `OPENAI_API_KEY`,
  `KIT_TOOLBOX_API_KEY`). For `backend=cli` log in once with `claude` / `codex`
  as `delfin`.
- Dashboard: service `delfin-dashboard` on `127.0.0.1:8866` only. On the Pi use
  the "DELFIN Dashboard" menu entry; from your PC use a tunnel
  `ssh -L 8866:127.0.0.1:8866 delfin@delfin-pi.local` and open the URL printed
  by `delfin-dashboard-url` on the Pi.

## Updating

```bash
scripts/pi/deploy_pi.sh delfin@delfin-pi.local            # current working tree, incl. uncommitted
scripts/pi/deploy_pi.sh --test delfin@delfin-pi.local     # + pytest -m "not slow"
scripts/pi/deploy_pi.sh --rollback delfin@delfin-pi.local
scripts/pi/deploy_pi.sh --status delfin@delfin-pi.local
```

Under Git Bash, `*.local` names are resolved through Windows (Git's ssh has no
mDNS). Only files git tracks or would track are deployed, so `.gitignore`d secrets,
venvs and caches stay on your PC. Files are sent with LF line endings even
from a `core.autocrlf=true` checkout.

An existing Pi without a prepared card: `deploy_pi.sh --bootstrap pi@<host>`
once, then deploy as above.
