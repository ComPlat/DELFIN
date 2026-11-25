# DELFIN for macOS

**New here?** `START-HERE.txt` beside this file is the short path: Gatekeeper
and the download mark, install, connect, DELFIN on the Mac. This page is the
reference behind it.

The Mac counterpart of `tools/windows-launcher/`. One app, two places to work:

* **A server** (Linux server, the lab Pi, a cluster) over SSH, with the same
  server half as the Windows launcher (`remote_bootstrap.sh`,
  `remote_launcher.py`) -- same check, same arguments, same server folder.
* **This Mac**: DELFIN installed into one folder of the account and served
  on 127.0.0.1, with the agent's commands in macOS's Seatbelt sandbox.

Nothing on the Mac needs Xcode, its command line tools, Homebrew or a Python:
the app is `/bin/bash` (3.2, as macOS ships it) and the system's `ssh`,
`scp`, `ssh-keygen`, `shasum`, `lsof`, `curl`, `open` and `osascript`.
`/usr/bin/python3` is not used: on a Mac without the command line tools it
is only a stub that offers to install them.

## Files

| File | What it is |
| --- | --- |
| `Install.command` | Double-click: download-mark question, package check, install, start |
| `Uninstall.command` | Double-click: remove, keeping saved connections and the local DELFIN on request |
| `DELFIN.sh` | The app: menu, SSH connection, local dashboard, install and uninstall |
| `setup_local.sh` | Sets DELFIN up on the Mac (micromamba, Python 3.11, DELFIN from PyPI) |
| `fetch_micromamba.sh` | DELFIN's pinned micromamba download (from `delfin/installers/`) |
| `remote_bootstrap.sh`, `remote_launcher.py` | The server half (from `tools/windows-launcher/`) |
| `VERSION` | The DELFIN release `setup_local.sh` installs |
| `SHA256SUMS.txt` | Checked by `Install.command` before anything is installed |

`build_package.py` builds `tools/delfin-macos-launcher.zip`, the file the
release publishes. The shared files are taken from where they live at build
time; `tests/test_the_macos_launcher_package_is_current.py` fails when the
ZIP no longer matches them.

## Where things go

| What | Where |
| --- | --- |
| The app | `~/Applications/DELFIN.app` (no administrator rights) |
| Saved connections | `~/Library/Application Support/DELFIN Launcher/profiles/<id>.conf` |
| An invitation key until its first use | `~/Library/Application Support/DELFIN Launcher/invite-<id>.key` (0600) |
| DELFIN on this Mac | `~/Library/Application Support/DELFIN Launcher/local/` |
| Working folder of the local dashboard | `~/DELFIN` |
| Remember-login key | `~/.ssh/delfin_launcher_ed25519` (as on Windows) |

The settings folder is 0700, files in it 0600 (`umask 077`).

## Gatekeeper and the download mark

The package is not signed and not notarized. A ZIP downloaded with a browser
carries the quarantine attribute (`com.apple.quarantine`), and Archive
Utility passes it to every extracted file; Gatekeeper then refuses to open
`Install.command` with "unidentified developer" (macOS 15: "Apple could not
verify ..."). The way through is Apple's own: System Settings > Privacy &
Security > **Open Anyway** (macOS 14 and earlier also right-click > Open),
or `bash Install.command` in Terminal.

`Install.command` then lists the package files that still carry the mark and
asks *Remove the mark and continue? (y/N)*. `y` runs
`xattr -d com.apple.quarantine` on exactly those files; `n` installs nothing.
The mark is never removed silently, and nothing else is changed: Gatekeeper
stays on, no `spctl` setting is touched. The installed app is copied with
`cp -p`, so it keeps whatever the person decided.

Files the app writes itself (the `.command` file a connection runs in,
the app bundle) are made on the Mac and carry no mark.

## Connecting to a server

The app's menu lists the saved connections. Choosing one opens a Terminal
window running `DELFIN.sh connect <id>`, which does what `DELFIN.ps1 -Mode
Dashboard` does:

1. One `ssh` connection runs the check of the server's copy of the two
   scripts and, when it is exact, the dashboard from it. When the copy is
   missing or different (exit 10, `DELFIN_SCRIPTS_MISSING`) `scp` copies them
   and the connection opens once more.
2. Password, OTP and host-key questions are `ssh`'s own, in the Terminal
   window. The app never reads or stores a password.
3. The dashboard reaches the Mac through one forward,
   `-L 127.0.0.1:<port>:/tmp/delfin-win-tunnel-<uuid>/dashboard.sock`.
   The browser is opened only after `lsof` shows that this connection's own
   `ssh` process listens on that port and `/login` answers; the token goes
   nowhere else. The server's URL is accepted only for `http://127.0.0.1`,
   the dashboard paths and a token.
4. Closing the window (or Ctrl+B, then D) disconnects; Keep session decides
   whether the dashboard stays on the server, as on Windows.

SSH options, identical to the Windows launcher: `StrictHostKeyChecking=ask`,
`ForwardAgent=no`, `ForwardX11=no`, `ClearAllForwardings` except for the one
dashboard forward (`GatewayPorts=no`, `ExitOnForwardFailure=yes`),
`PermitLocalCommand=no`, `ControlMaster=no`, `ControlPath=none`, and with
Remember login `IdentitiesOnly=yes -i ~/.ssh/delfin_launcher_ed25519`
(each key the macOS agent or `~/.ssh/id_*` would offer costs one of the
server's MaxAuthTries). The system `ssh` in `/usr/bin` is used, never one
from Homebrew or `PATH`.

A server with a different port, a jump host or a user name of its own is an
alias in `~/.ssh/config`; the app takes the alias as the server.

**Remember login**, **Keep session**, **Working shell** and the advanced
settings (working directory, DELFIN location, starting port) mean what they
mean in the Windows launcher; the server side is the same code.

## Invitations (`delfin-people --invite`)

`delfin-people --invite NAME` on the server writes `DELFIN-NAME.zip`. When
the server has `tools/delfin-macos-launcher.zip` (or `DELFIN_MAC_LAUNCHER_ZIP`
names one), the invitation carries both packages, each with `invite.json`
and `invite_key`: `windows-launcher/` for a PC, `macos-launcher/` for a Mac.
The key is the same and works once, on whichever computer uses it first.

On the Mac, `Install.command` imports it: the server becomes the first saved
connection (Remember login and Keep session on), the server's host keys go
into `~/.ssh/known_hosts` (no fingerprint question), and the key is kept as
`invite-<id>.key` (0600). The first connection registers this Mac's own key
with it (`delfin-people --enroll`, `DELFIN_ENROLLED`) and deletes it.

## DELFIN on this Mac

"Set up DELFIN on this Mac..." runs `setup_local.sh`: micromamba (pinned
release, SHA-256 checked), Python 3.11 from conda-forge, and
`delfin-complat[agent,docs]==<VERSION>` from PyPI, all under
`~/Library/Application Support/DELFIN Launcher/local/`. `install.sh` is not
used on the Mac: it sets up a Linux login node (ORCA, OpenMPI, QM tools) and
stops on macOS. ORCA is not installed (licence).

"DELFIN on this Mac" then starts `delfin-voila --open-browser --ip 127.0.0.1`
in a Terminal window, with `DELFIN_AGENT_SANDBOX=auto` (Seatbelt on macOS)
and `DELFIN_AGENT_SANDBOX_NETWORK=0`. The dashboard's token is mandatory;
the browser is opened with it by `delfin-voila` itself. Closing the window
stops the dashboard.

## Uninstall

`Uninstall.command` (or the app's menu) removes the app and the settings
folder. It asks whether to keep the saved connections (and unused invitation
keys) and whether to keep the local DELFIN. The login key in `~/.ssh` stays;
nothing on any server is changed.

For scripts: `DELFIN.sh uninstall --keep-profiles|--remove-profiles
--keep-local|--remove-local --close-running`.

## What is tested where

* Linux CI (`pytest`): the shared rules -- the server command is the Windows
  launcher's, byte for byte; SSH options; URL check; profiles; the
  invitation import; install and uninstall in an isolated home; the package
  ZIP matches its sources; the release publishes it.
* macOS CI (`.github/workflows/macos-launcher.yml`, `macos-latest`): the
  same under the Mac's `/bin/bash` 3.2 and BSD tools; the package extracted
  with `ditto` (what Archive Utility uses); the quarantine question both
  ways; install, icon, `plutil -lint`; uninstall keeping and removing;
  DELFIN set up on the Mac from the tested tree and its dashboard answering
  with the token and refusing without it, bound to 127.0.0.1 only; and an
  SSH connection to a `sshd` on the runner.
* Not in CI: Finder, Gatekeeper's dialog and "Open Anyway", Terminal.app
  opening the window, the browser. These need a person at a Mac.
