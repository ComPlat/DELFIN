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
- `/etc/delfin/managed_policy.json` (root, written by `pi_bootstrap.sh`):
  isolation required (bwrap), and no session on a CLI backend (`claude`,
  `codex`), whose own shell runs outside the cage with the agent user's
  files. Sessions start on the KIT toolbox (`pi_install.sh` sets it when no
  provider is chosen). `CLI_BACKENDS=allowed` at bootstrap lifts the backend
  rule; the user cannot.
- DELFIN's own MCP servers run in bwrap (`DELFIN_MCP_ISOLATION=builtin` in
  `agent.env`): they see the calculation, office and workspace folders, not
  `~/.ssh` or the credentials.
- Push access reaches `main` only through review: a deploy key, `main`
  protected by a ruleset, and no `gh` login on the Pi (a token that may write
  pull requests may also approve them). The agent pushes branches and hands
  you the compare URL; a person opens and approves the pull request.

The host around it (`pi_harden_host.sh`, run by the first boot; run it again
as root any time, it is idempotent):

- SSH by key only for every user, no root login, no X11 forwarding, 3 tries,
  only the login and the agent user. fail2ban bans repeated failures.
- nftables: inbound default deny. SSH only from private IPv4 nets (LAN, VPN)
  and link-local IPv6; add trusted IPv6 prefixes with `SSH_V6="..."` (a Pi
  on a routed IPv6 network is otherwise reachable from far beyond the LAN).
  mDNS, DHCP and ICMP are allowed, outbound is open (agent egress is
  DELFIN's proxy's job). Over SSH the firewall is armed with a 120 s
  rollback: log in again over SSH within that time to keep it.
- rpcbind / NFS off.
- Debian security updates install themselves (unattended-upgrades).
  Raspberry Pi kernel and firmware stay manual, so the Landlock kernel is
  never replaced behind your back: `sudo apt update && sudo apt full-upgrade`
  from time to time, then check `uname -r` after the reboot.
- sudo asks for the login user's password once it has one: set it with
  `sudo passwd pi`, then run `pi_harden_host.sh` again.

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
| `pi_harden_host.sh` | Pi (root) | Host lockdown: sshd key only, nftables, fail2ban, security updates, sudo password |

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

## After the first boot

What only a person can do, in this order:

1. `pi`: set a password (`sudo passwd pi`), then
   `sudo bash ~delfin/delfin/src/scripts/pi/pi_harden_host.sh` again, so sudo
   asks for it.
2. `delfin`: the API key, `delfin-agent credentials set KIT_TOOLBOX_API_KEY`,
   and `git config --global user.name/user.email`.
3. GitHub (a repository admin): a ruleset on `main` -- *Settings -> Rules ->
   Rulesets*, target the default branch, active, empty bypass list; require a
   pull request with 1 approval and approval of the most recent push, block
   force pushes, restrict deletions.
4. `delfin`: the deploy key (`pi_github_deploy_key.sh`, then `--check`), and
   check that a push to `main` is refused.
5. `delfin-agent doctor`: expected warnings are xtb, orca, the `gh` login and
   benchmark ground truth.

## Accounts and access

- `pi`: desktop autologin, SSH by key only, password locked. Set one with
  `sudo passwd pi` and rerun `pi_harden_host.sh`, so sudo asks for it.
- `delfin`: agent user, no root (it may restart its dashboard service), SSH by
  the same key. More devices (a phone with Termius, another laptop): create a
  key there and append its public line to `~/.ssh/authorized_keys` of the
  user; passwords are not accepted over SSH.
- Over a VPN `*.local` (mDNS) does not resolve: use the Pi's IP
  (`hostname -I` on the Pi) in `~/.ssh/config`.
- API keys are entered on the Pi, never on the card:
  `delfin-agent credentials set KIT_TOOLBOX_API_KEY` (or `ANTHROPIC_API_KEY`,
  `OPENAI_API_KEY`). With the KIT toolbox pass a model it serves, e.g.
  `delfin-agent run --provider kit --model kit.deepseek-v4-flash "..."`.
  The CLI backends (`claude`, `codex`) are refused by the managed policy on
  the Pi (see Security); they stay installed for a box bootstrapped with
  `CLI_BACKENDS=allowed`.
- Once, as `delfin`: `git config --global user.name/user.email`, so the agent
  can commit.
- Push access, scoped to one repository and stopping at the pull request: a
  deploy key (any repository admin can add one), after `main` is protected
  (see the checklist below): `bash ~/delfin/src/scripts/pi/pi_github_deploy_key.sh`,
  add the printed key with "Allow write access", then
  `... pi_github_deploy_key.sh --check`. Do not log `gh` in on the Pi.
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
once, then deploy as above; then, as `pi`,
`sudo bash ~delfin/delfin/src/scripts/pi/pi_harden_host.sh` and
`deploy_pi.sh --landlock pi@<host>`, and the checklist "After the first
boot".

## SimpleX messenger: forum and session chats

Set up in one run on the Pi (or any Linux with systemd and glibc >= 2.34),
after a deploy, as the admin user (`pi`; the agent user has no sudo):

```bash
sudo apt install qrencode          # optional: the pairing link as a QR code
sudo bash ~delfin/delfin/src/scripts/pi/pi_messenger.sh --whisper --pair --name=Pi
```

`--name=` labels the bot (`DELFIN-Pi`) and every chat, for a phone paired with
several machines; the default is the host name. The bot keeps the name of its
first start. `--lang=en` sends every phone text in English (default German).

The script waits for the SimpleX daemon, then prints the QR code and, last,
the link (on a phone's terminal: copy the link and paste it in SimpleX under
+ > paste link). Once the phone connects it prints the security code -- compare
it on the phone under the contact > Verify security code -- the bot says
hello in SimpleX and invites you to the forum for all agents. The script installs the pinned simplex-chat build (its daemon creates the bot
profile), the gateway venv (hashed), the units including the agent service, and
with `--whisper` local speech recognition (pinned faster-whisper-small, every
file digest-checked, about 0.5 GB). After a DELFIN update run it again.

Then, as `delfin`, in a project folder: `delfin-agent tui` and `/sync` in the
session; within a few seconds the gateway creates the session's chat and
invites you, nothing to type on the phone. The group is named after
machine and session (an untitled session after its folder). A closed TUI
leaves its session running: `delfin-agent tui --continue` in the same folder
joins it again.

On the phone, `neu NAME` (or `neu NAME: task`) starts a session in
`~/DELFIN-Handy/NAME` and invites you to its chat; it asks before commands and
edits. A session started at the computer joins with `/sync` there. Reply to a
question/approval with SimpleX's Reply gesture and `ja/nein`, an offered option or
a text answer -- spoken or typed; a spoken sentence frees nothing. A shell
command can be changed before it runs: reply `ändern: NEW COMMAND` (typed only);
your version runs once, exactly as typed. Late replies
never approve the next action. `status`, `update`, `pending`, `stop`, `name NAME`,
`dateien`, `datei NAME`, `forum on/off`, `meldungen kurz/fragen/aus`; files and
images up to 8 MiB.

See [setup, security boundaries and command reference](../../docs/design/messenger-simplex.md).

## Schedules

The agent service runs your schedules (a task at set times) as sessions:

```bash
delfin-agent schedule add täglich 07:30 --cwd ~/proj Fasse die Jobs zusammen
delfin-agent schedule            # list, with the next run
```

In a session: `/schedule add mo-fr 18:00 <task>`. On the phone: `zeitpläne`
opens the chat "<machine> · Zeitpläne"; there `neu alle 2h <task>` proposes one
(in `~/DELFIN-Handy/zeitplaene`) and `ok N` starts it; in a session's chat
`plan täglich 07:30 <task>` proposes one in that session's folder. Every run is
reported there. Asked "every day at 7:30 ...", the agent proposes one; it runs
only after your `ok`. Ceilings: $1 and 30 min a run, $2 a day per schedule.
Watchers start a task on an event instead of a time: `neu wenn datei
out.txt erscheint <task>`, `wenn job 123 fertig`, `wenn datei calc.log sich
ändert`, `wenn ausgabe "squeue -u NAME" sich ändert`; without a task only a
report.
See [schedules](../../docs/design/schedules.md).
