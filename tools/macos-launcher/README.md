# DELFIN for macOS

**New here?** `START-HERE.txt` beside this file is the short path. This page
is the reference behind it.

The Mac counterpart of `tools/windows-launcher/`. One app, two places to work:

* **A server** (Linux server, the lab Pi, a cluster) over SSH, with the same
  server half as the Windows launcher (`remote_bootstrap.sh`,
  `remote_launcher.py`) -- same check, same arguments, same server folder.
* **This Mac**: DELFIN installed into one folder of the account and served
  on 127.0.0.1, with the agent's commands in macOS's Seatbelt sandbox.

Nothing on the Mac needs Xcode, its command line tools, Homebrew or a Python:
the app is `/bin/bash` (3.2, as macOS ships it) and the system's `ssh`,
`scp`, `ssh-keygen`, `shasum`, `lsof`, `curl`, `ditto`, `open` and
`osascript`. `/usr/bin/python3` is not used: on a Mac without the command
line tools it is only a stub that offers to install them.

## Install: one line

The repository's README.md carries the line for the current release (it
names the ZIP's SHA-256, so it cannot live inside the ZIP). In Terminal:

```bash
d=$(mktemp -d) && curl -fL --proto =https -o "$d/delfin.zip" https://github.com/ComPlat/DELFIN/releases/download/v<VERSION>/delfin-macos-launcher.zip && echo "<SHA-256>  $d/delfin.zip" | shasum -a 256 -c && ditto -x -k "$d/delfin.zip" "$d" && bash "$d/macos-launcher/Install.command"
```

* A fixed release asset (`v<VERSION>`, HTTPS only), checked against its
  SHA-256 **before anything of it runs**; a wrong digest stops the line.
  Nothing is piped into a shell.
* `curl` does not set the download mark, so Gatekeeper has nothing to ask:
  that is how macOS treats a file a person fetches in Terminal, not a way
  around it. The same holds for a signed and notarized build.
* `build_package.py` writes the line; `tests/test_the_macos_launcher_package_is_current.py`
  fails when it does not name the committed ZIP and this version.

**From Finder** instead: download `delfin-macos-launcher.zip`, double-click
it, double-click `Install.command` in the folder (see Gatekeeper below).

## Package

```
macos-launcher/
  Install.command      download-mark question, package check, install, start
  Uninstall.command    remove; keeps saved connections / local DELFIN on request
  SHA256SUMS.txt       every file below, checked before anything is installed
  README.md  START-HERE.txt
  DELFIN.app/Contents/
    Info.plist  MacOS/DELFIN
    Resources/DELFIN.sh            menu, SSH connection, local dashboard, install
              setup_local.sh       DELFIN on the Mac (micromamba, Python, PyPI)
              fetch_micromamba.sh  DELFIN's pinned micromamba (delfin/installers/)
              remote_bootstrap.sh  remote_launcher.py   (tools/windows-launcher/)
              DELFIN.icns  VERSION
```

`DELFIN.app` is a finished bundle, built with the package (`build_package.py`)
and never assembled or changed on the Mac: the installer copies it with
`ditto`, so a signature on it stays valid. Sources: `app/Info.plist`,
`app/DELFIN`, `app/DELFIN.icns` (made from `tools/windows-launcher/DELFIN_logo.png`).

## Where things go

| What | Where |
| --- | --- |
| The app | `~/Applications/DELFIN.app` (no administrator rights; Spotlight and Launchpad find it) |
| Saved connections | `~/Library/Application Support/DELFIN Launcher/profiles/<id>.conf` |
| An invitation key until its first use | `~/Library/Application Support/DELFIN Launcher/invite-<id>.key` (0600) |
| DELFIN on this Mac | `~/Library/Application Support/DELFIN Launcher/local/` |
| Working folder of the local dashboard | `~/DELFIN` |
| Remember-login key | `~/.ssh/delfin_launcher_ed25519` (as on Windows) |

The settings folder is 0700, files in it 0600 (`umask 077`). No `launchd`
job, login item or background service is installed: the dashboard runs while
its Terminal window is open, and the server's dashboard is the server's.
`DELFIN_LAUNCHER_APPS=/Applications bash Install.command` installs for all
accounts of the Mac (needs write access there).

## Gatekeeper and the download mark

A ZIP downloaded with a browser carries the quarantine attribute
(`com.apple.quarantine`), and Archive Utility passes it to every extracted
file. Gatekeeper then opens only software Apple has notarized.

* **Signed and notarized build** (see below): `Install.command` checks
  `DELFIN.app` with `spctl --assess` and, when Gatekeeper accepts it, asks
  nothing and keeps the mark. Only `Install.command` itself, a plain script,
  is asked about once by macOS in Finder ("Open Anyway"); a script cannot
  carry a stapled ticket.
* **Unsigned build** (the default until a Developer ID is set up):
  `Install.command` lists the package files that carry the mark and asks
  *Remove the mark and continue? (y/N)*. `y` runs
  `xattr -dr com.apple.quarantine` on exactly those; `n` installs nothing.
  Never silent; Gatekeeper stays on; no `spctl` setting is changed and
  nothing needs `sudo`.

Getting `Install.command` to run from Finder in the first place: System
Settings > Privacy & Security > **Open Anyway** (macOS 14 and earlier also
right-click > Open), or `bash Install.command` in Terminal.

## Signing and notarization

`sign_and_notarize.sh` (maintainer-side, on a Mac with Xcode's command line
tools) does it with Apple's own tools, from two settings outside the
repository: `DELFIN_SIGN_IDENTITY` ("Developer ID Application: ...") and
`DELFIN_NOTARY_PROFILE` (made once with `xcrun notarytool store-credentials`).

1. `codesign --options runtime --timestamp` on `DELFIN.app` (hardened
   runtime, secure timestamp), `codesign --verify --deep --strict`.
2. `xcrun notarytool submit --wait`, `xcrun stapler staple`, `spctl --assess`.
3. The ZIP is rebuilt from the signed tree. Signing adds
   `Contents/_CodeSignature/` and `Contents/CodeResources`; no checked file
   changes (the macOS CI signs the bundle ad hoc and verifies exactly that).
4. Put the printed SHA-256 into README.md's one-line install, commit ZIP and
   line, then tag.

Next step, not built here: a signed installer package (`pkgbuild` +
`productsign` with a "Developer ID Installer" certificate, notarized the
same way) would also remove the one Finder question about `Install.command`.

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
the browser is opened with it by `delfin-voila` itself. As on a server,
`delfin-voila` started from a terminal runs its server in the background,
one per account: the window may be closed, starting again opens the same
dashboard, and it ends itself 10 minutes after its last browser window
closed. "Stop DELFIN on this Mac" in the menu (`delfin-voila stop`) ends it
at once. Jupyter's kernel sockets go to `/tmp/delfin-jupyter-<uid>` (0700,
checked to be this account's own): the default folder's path is too long
for a Unix socket on macOS, and kernels would otherwise use loopback TCP.

## Uninstall

`Uninstall.command` (or the app's menu) removes the app and the settings
folder. It asks whether to keep the saved connections (and unused invitation
keys) and whether to keep the local DELFIN. The login key in `~/.ssh` stays;
nothing on any server is changed. Nothing else was installed: no launchd
job, no login item, nothing outside the app and the settings folder.

For scripts: `DELFIN.sh uninstall --keep-profiles|--remove-profiles
--keep-local|--remove-local --close-running`.

## What is tested where

* Linux CI (`pytest`): the shared rules -- the server command is the Windows
  launcher's, byte for byte; SSH options; URL check; profiles; the
  invitation import; install (the bundle copied unchanged) and uninstall in
  an isolated home; the package ZIP matches its sources; the one-line
  install names this ZIP and version; the release publishes it.
* macOS CI (`.github/workflows/macos-launcher.yml`, `macos-latest`): the
  same under the Mac's `/bin/bash` 3.2 and BSD tools; the one-line install
  (with the committed ZIP in place of the download) and a wrong digest
  stopping it; the package extracted with `ditto` and marked; the
  download-mark question both ways; `plutil -lint`, the icon, an ad-hoc
  signature verified with `--strict`; uninstall keeping and removing; DELFIN
  set up on the Mac from the tested tree and its dashboard answering with
  the token and refusing without it, bound to 127.0.0.1 only; and an
  invitation plus SSH connection to an `sshd` on the runner.
* Not in CI: Finder, Gatekeeper's dialog and "Open Anyway", Terminal.app
  opening the window, the browser, a Developer ID signature and
  notarization. These need a person at a Mac (and, for signing, a
  Developer ID).
