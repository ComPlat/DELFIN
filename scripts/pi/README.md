# DELFIN agent box on a Raspberry Pi 5

Turns a Raspberry Pi 5 (aarch64) into a dedicated machine for developing and
testing the DELFIN agent. The agent runs as its own user `delfin`, without
root: free inside its frame, and the frame is not its to widen (see
"Security" below).

QM binaries (ORCA, xtb, CREST, ...) are **not** installed. Everything else is:

- DELFIN with **every extra in `pyproject.toml`** except `ai-retro`, `orca`
  and `pymol` (new extras reach the Pi with the next deploy by themselves;
  override with `DELFIN_PI_SKIP_EXTRAS`), torch as the CPU build;
- the terminal UI (`delfin-agent tui`) with the patched OpenTUI library, so
  rows stay in the scrollback of xterm.js terminals (VS Code, JupyterLab);
- the Claude Code and Codex CLIs (agent `backend=cli`), the GitHub CLI `gh`;
- Playwright's Chromium (browser extra), tesseract + pytesseract (OCR),
  the manual search index (`delfin-docs-index`);
- a git checkout `~/delfin/work/DELFIN` for the agent to work on;
- the dashboard as a loopback-only service;
- the Raspberry Pi OS kernel rebuilt with Landlock, beside the stock kernel
  (see "Landlock kernel" below).

Check the result with `delfin-agent doctor` (expected warnings: xtb, orca,
git identity / `gh` login until you set them, benchmark ground truth).

Excluded extras: `ai-retro` (aizynthfinder pins `rdkit<2023`, which conflicts
with DELFIN's `rdkit>=2023.3` on every platform), `orca` (QM). PyMOL has no
aarch64 distribution on PyPI and is built from a pinned source commit
(`--no-deps`, it runs on DELFIN's numpy<2); if that build fails only the PyMOL
images in DOCX/ESP reports are missing.

## Security

Free inside the frame, and nothing widens the frame but a person:

- The agent's shell commands run in DELFIN's cage (`agent.bash_isolation`
  `auto`: bubblewrap + Landlock): only the workspace is writable, `.git`
  hooks and config read-only, credentials, `~/.ssh` and the dashboard token
  unreadable, no other process visible, no sudo (no-new-privileges), network
  only through the egress proxy. Inside it, commands run without a question;
  writes in default mode are asked.
- The agent user `delfin` has no root (`AGENT_SUDO=service`: it may restart
  and query `delfin-dashboard`, nothing else). What runs outside the cage --
  the shell of a CLI backend, the MCP servers, the agent process -- runs with
  that user's rights only. Administration is the login user's (`pi`).
  `AGENT_SUDO=none` drops even the service rule; `AGENT_SUDO=full` restores
  root without a password and is not advised.
- `DELFIN_AGENT_SANDBOX=auto` in `~/delfin/agent.env`: a command the dashboard
  runs for the agent is confined as well (`off` ran it with the account's
  full rights; `pi_install.sh` changes an old `off`).
- A request of the agent for more rights (an auto-allow rule, a folder, a
  mode) is flagged as such, asked twice, and never right after a refusal.
- Push access: with a deploy key the agent can push any branch. Protect
  `main` in the repository settings (pull requests only), so what the agent
  pushes reaches `main` through review.

## Landlock kernel

The cage of the agent's commands is bubblewrap plus Landlock. Raspberry Pi OS
builds its kernels without `CONFIG_SECURITY_LANDLOCK`. Only with Landlock
(which confines a command's reads) does a command off the auto-allow list run
in the cage without a question; without it the agent asks about each one. `pi_kernel_landlock.sh` builds the installed stock kernel again
(the `linux-source` package of the same version, from the signed Raspberry Pi
archive, and the same `/boot/config-<release>`) with only
`CONFIG_SECURITY_LANDLOCK=y` and `CONFIG_LSM="landlock"` changed:

- its files are its own (`kernel_2712_landlock.img`, `initramfs_2712_landlock`,
  `/lib/modules/<release>-landlock`); apt and the raspi-firmware hooks do not
  write them, and the stock kernel stays installed;
- the first boot is a tryboot with `panic=10`: a kernel that does not come up
  reboots into the stock one by itself. A boot unit adds the kernel line to
  `config.txt` only once Landlock is active, inside the board's section
  (`[pi5]`, or `[pi4]` for a Pi 4 with the `rpi-v8` kernel), so the card in
  another model boots that model's stock kernel;
- `status` says when apt brought a newer stock kernel: run `build` and `try`
  again, the Landlock kernel does not get those updates by itself.

A prepared card does this on first boot (`LANDLOCK_KERNEL=0` skips it). An
existing Pi: `scripts/pi/deploy_pi.sh --landlock delfin@<host>`, or on the Pi
`sudo bash ~/delfin/src/scripts/pi/pi_kernel_landlock.sh build` then `... try`.
`... remove` goes back to the stock kernel.

Known limits, not Pi-specific: `schnetpack` does not import with torch 2.13
(`T_co` was removed from torch). Architector and molSimplify are installed
into the venv by pip, but build through xtb's Python module (`xtb-python`,
conda-forge only), so their builds fail until they run in an environment
that has it (`DELFIN_ARCHITECTOR_PYTHON`, `DELFIN_MOLSIMPLIFY_PYTHON`).
`uv pip check` reports pymol's numpy>=2 metadata and the `pymol` vs
`pymol-open-source` name; both are expected. The deployed tree has no `.git`,
so `test_pre_refactor_api_preserved` (needs `git show`) fails under `--test`.
The molsimplify and mace-torch wheels install a top-level `tests` package that
would hide the repo's `tests/` (`from tests.conftest import ...`);
`pi_install.sh` removes it after every deploy.

First boot takes ~20 min on a Pi 5 with a fast connection (about 9 of them for
the patched OpenTUI build). A deploy keeps `tui/node_modules` (and with it the
patched library) while `tui/package.json`, `bun.lock` and `tui/patches` are
unchanged; then it takes seconds.

## Files

| Script | Runs on | Purpose |
| --- | --- | --- |
| `make_pi_card.ps1` | Windows | Flash Raspberry Pi OS and put DELFIN + first-boot setup on the card |
| `pi_kernel_landlock.sh` | Pi (root) | Rebuild the stock kernel with Landlock, tryboot into it, keep it once it runs |
| `deploy_pi.sh` | PC (Git Bash / WSL / Linux / macOS) | Push the current working tree to a running Pi, roll back, run tests |
| `pi_bootstrap.sh` | Pi (root) | System setup: packages, agent user, uv + Python 3.11 venv, CLIs, dashboard service |
| `pi_install.sh` | Pi (agent user) | Unpack a source tarball, install deps when `pyproject.toml` changed, keep `src.prev` for rollback |
| `pi_firstboot.sh` | Pi (cloud-init) | First boot: `pi_bootstrap.sh` + `pi_install.sh` from the card |
| `pi_github_deploy_key.sh` | Pi (agent user) | Push access to one repository via a deploy key, GitHub host key pinned |

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
- `delfin`: agent user, no root (it may restart its dashboard service), SSH by
  the same key.
- API keys are entered on the Pi, never on the card:
  `delfin-agent credentials set KIT_TOOLBOX_API_KEY` (or `ANTHROPIC_API_KEY`,
  `OPENAI_API_KEY`). With the KIT toolbox pass a model it serves, e.g.
  `delfin-agent run --provider kit --model kit.deepseek-v4-flash "..."`.
  For `backend=cli` log in once with `claude` / `codex` as `delfin`.
- Once, as `delfin`: `git config --global user.name/user.email`, so the agent
  can commit.
- Push access, scoped to one repository: prefer a fine-grained token
  (Resource owner = the organisation, only that repository, Contents +
  Pull requests read/write) via `gh auth login --with-token`. Avoid the
  browser `gh auth login`: its `repo` scope covers every repository you can
  write to. If the organisation does not offer fine-grained tokens, use a
  deploy key (any repository admin can add one):
  `bash ~/delfin/src/scripts/pi/pi_github_deploy_key.sh`, add the printed key
  with "Allow write access", then `... pi_github_deploy_key.sh --check`.
  With a deploy key the agent pushes branches and hands you the compare URL
  for the pull request.
- VS Code: Remote-SSH to `delfin@delfin-pi.local`.
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
