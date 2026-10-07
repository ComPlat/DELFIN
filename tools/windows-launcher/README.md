# DELFIN Windows Launcher

Copy this entire folder to your Windows PC. DELFIN and calculations stay on your
Linux server. Windows does not need Python or a local DELFIN installation.
The app includes the DELFIN logo, a desktop shortcut, saved SSH connections,
automatic environment discovery, automatic port selection, and uninstall support.
All app text, terminal messages and instructions are in English.

## What happens when you start

1. One Windows OpenSSH connection opens. Enter your password and OTP when the
   server requests them. The app does not store credentials or bypass MFA.
2. DELFIN and its environment are discovered automatically. A verified existing
   dashboard is reused; otherwise the app starts a new dashboard in tmux.
3. The same SSH connection carries the browser tunnel. The browser signs in using
   DELFIN's mandatory access token. Server and Windows ports are selected automatically.
4. If **Working shell (Ctrl+B, then N)** is enabled, an additional tmux window provides a
   working shell through the **same SSH login**. Switch between dashboard and
   working shell with **Ctrl+B, then N**. These are terminal views inside one
   Windows window, not two separately authenticated Windows SSH windows.

There is one SSH authentication sequence for this connection. The server can
still request multiple authentication factors, retries or jump-host logins.
Closing the browser does not close SSH or stop the dashboard.

## App, dashboard terminal and browser

The installer builds `DELFIN.exe` locally using the Windows .NET Framework
compiler and the included readable `Starter.cs`. The desktop/Start menu shortcuts
open the settings app without a PowerShell console. The settings app is an
independent main window and stays in the taskbar with the DELFIN icon while open.

A separate dashboard terminal appears for native OpenSSH password/OTP and
host-key prompts and stays visible after the browser opens. This is the only
SSH terminal; it uses the PowerShell icon and an SSH title. The DELFIN logo
belongs to the settings app, not a second terminal. Use Ctrl+C there to
stop the dashboard explicitly, including with Keep ON. The optional working
terminal is still OFF by default. Keep the DELFIN app open while using the browser.
**Show terminal** brings back the SSH window if it was hidden.
**Disconnect** ends all connections started by this app instance. Closing the
DELFIN settings app ends all connections started by that instance, including a
pending login. Keep OFF stops dashboards created by those connections; Keep ON
leaves the server for reconnecting. Pre-existing dashboards keep their lifetime.
The status line shows Connecting, Connected, Disconnecting or Disconnected.
Only one settings app instance runs per Windows login; clicking the desktop
shortcut again brings that existing window forward. The launcher permits only
one active SSH connection. Repeated Start / Reconnect clicks do not start extra
connections or dashboards. Open further views inside the DELFIN dashboard in
your browser. Disconnect before changing the SSH destination.
Closing a browser tab alone does not disconnect SSH. Authentication/startup
failures remain visible. Enable the optional working terminal to keep the SSH
window visible and switch tmux windows. Reinstall to update the shortcuts; unpin
an old PowerShell-based taskbar shortcut and pin the new DELFIN shortcut.

## Taskbar shortcut

The installer creates DELFIN shortcuts with the logo on the desktop and in the
current user's Start menu. Find **DELFIN** in Start, right-click it and select
**Pin to taskbar** (under **More** on some Windows versions).
The settings window supplies Windows with DELFIN.exe as its relaunch command,
the DELFIN icon and its own window AppUserModelID, so pinning that open window
also targets DELFIN rather than its internal PowerShell host. Remove a previously
pinned PowerShell entry before pinning the updated DELFIN window.
This installer does not force a taskbar pin or change your taskbar layout.
Updates preserve the shortcut location. Before uninstalling, right-click a pinned
DELFIN icon and select **Unpin from taskbar**; the uninstaller removes the desktop
and Start menu shortcuts.

## Session controls

Two options are visible without opening Advanced settings:

- **Working shell (Ctrl+B, then N)**: disabled by default for new saved
  connections. Enable it to open/reuse a tmux window named `DELFIN-work`.
- **Keep session after disconnect**: disabled by default for new saved connections.
  When off, disconnecting stops the **new dashboard started by this connection**
  and its tmux shell windows. When on, it uses DELFIN `--keep` and keeps the server
  and opted-in dashboard kernels available for reconnecting.

The app never kills a dashboard that existed before this connection. Its previous
lifetime settings remain in effect even if the app's Keep checkbox is off; the
terminal explicitly reports this exception. Stop such a server deliberately in
its original/dashboard terminal if needed.

**Ctrl+C is an explicit stop**, even with Keep enabled. If pressed in the working
shell, it normally interrupts that shell's foreground command; if pressed in the
dashboard view, it stops the dashboard process. To keep the server, disconnect or
use **Ctrl+B, then D** rather than Ctrl+C.

The server DELFIN installation must include the launcher release with
`delfin-voila --strict-port`. The app checks this before creating a dashboard.
Updating the Windows app alone does not update the server code. After the release
is merged, update the server checkout using your normal update process
(`git pull --ff-only` for a clean checkout on main).

Automatic discovery also checks the current user's `~/software/delfin`,
`~/software/DELFIN`, `~/delfin` and `~/DELFIN` when no usable environment was found
in the login environment or current directory. A local `.venv`/`venv` does not
need to be activated permanently. Source checkouts work without an editable
package install when their environment has the required dependencies. Multiple
usable conventional installations require an explicit location in Advanced.
For another installation path, set **DELFIN location** to that repository or its
Python environment; discovery does not scan the whole filesystem.
The selected PATH, Python source path and virtual-environment variables are
passed explicitly into tmux for the dashboard and working shell.

## Requirements

### Windows

- Windows 10/11, Windows PowerShell 5.1 and the optional **OpenSSH Client**.
- Check with `Get-Command ssh.exe`. Install the client through Windows Settings
  → Optional features if it is missing.
- A default browser, network access to the SSH server and any required VPN.
- Your organization must permit these scripts. No global execution policy is
  changed. Organization policies take precedence; obtain IT approval/signatures
  if required. The app uses process-scoped `RemoteSigned`, never `Bypass`.

### Server

- Current DELFIN including the `--strict-port` option and a working dashboard.
- An SSH server allowing TCP-to-Unix-socket forwarding (OpenSSH streamlocal forwarding).
- Bash, `base64`, and tmux for new persistent dashboards.
- A Python environment containing DELFIN and its dashboard dependencies.
- Any explicitly selected working directory must exist.

Server checks in your usual SSH session:

```bash
command -v tmux
command -v delfin-voila
delfin-voila --help
```

A Windows app update does not update DELFIN on the server. Use your normal server
installation procedure in the existing DELFIN environment when necessary.

## Install on Windows

1. Download the ZIP and extract it into a new folder.
2. Review `README.md`, `DELFIN.ps1`, `Install.ps1`, `Uninstall.ps1`,
   `remote_launcher.py` and `remote_bootstrap.sh`. If Windows marks downloaded
   scripts as blocked, use Properties → Unblock only after verifying their
   origin and contents. Do not bypass organization policies.
3. **Double-click `Install.cmd`** in the extracted `windows-launcher` folder.
   It runs the checks first, then installs the app. The window stays open so
   you can read the result. No PowerShell commands need to be typed.
   PowerShell remains the underlying Windows runtime; this wrapper does not
   bypass script blocking, SmartScreen, or organization policies.

   Alternatively, open the folder in Explorer. Enter `powershell`
   in the address bar and press Enter. Run:

   ```powershell
   powershell.exe -NoProfile -ExecutionPolicy RemoteSigned -File .\Test-Launcher.ps1
   powershell.exe -NoProfile -ExecutionPolicy RemoteSigned -File .\Install.ps1
   ```

4. The installer copies the app into `%LOCALAPPDATA%\DELFIN Launcher`, creates
   the **DELFIN** desktop shortcut, and registers **DELFIN SSH Launcher** under
   Windows Settings → Apps → Installed apps. Administrator rights are not needed.
5. Open DELFIN and enter **SSH server / user@host**, for example
   `delfin-cluster` or `alice@login.example.org`. This is the only required field.
   The connection name is generated automatically.
6. Click **Start / Reconnect** and complete the SSH login when prompted.

Saved connections contain no passwords or tokens. They are stored in
`%LOCALAPPDATA%\DELFIN Launcher\profiles.json`. **New** clears the form for a new
connection. **Save** saves changes; starting also saves automatically.

## Optional advanced settings

Normally **Advanced settings** can stay closed.

| Field | Default and examples |
|---|---|
| Working directory | Empty = your server home; otherwise a project/repository directory such as `~/projects/delfin` |
| DELFIN location | Empty = discovery; otherwise a repository, venv/Conda directory, Python executable or delfin-voila path |
| Starting port | `0` = automatic OS assignment; alternatively a starting port such as `8866` |

Absolute paths, `~`, `~/...`, and relative Linux paths are supported. Surrounding
quotes are removed. A relative working directory resolves from the login shell
directory; a relative DELFIN location resolves from the selected working directory.
Do not enter Windows paths such as `C:\Users\...`. These fields are paths, not shell commands.
Existing saved connections from previous versions are preserved during upgrades.
Missing optional fields are filled with automatic defaults when loading older profiles.

### Automatic environment discovery

The app starts an interactive Bash login shell, allowing your usual automatic
venv, Conda or module activation to run. It checks:

1. An explicit DELFIN location, if provided.
2. An active `VIRTUAL_ENV` or `CONDA_PREFIX`.
3. `.venv`, `venv`, `env` or `.env` in the working/repository directory and its
   parent directories, using Python environments that can import DELFIN.
4. `delfin-voila`, `python` or `python3` available through PATH with DELFIN.

A discovered conventional venv is activated. The dashboard starts through
`delfin.cli_voila` using its Python interpreter; a fixed console-script path is
not required. Repository source can be used with a compatible environment.
The initial interactive login shell loads your normal startup files once. An
already active matching venv is preserved; a discovered inactive venv is activated
only when required. The final working shell inherits the environment and prompt
and does not run `.bashrc` again. Non-exported shell functions and aliases from
the initial shell are not automatically inherited by that final shell.

If you normally change into the repository and activate its local environment,
set that repository as the optional working directory once. If an installation
is not visible through login/PATH or the selected repository, provide its location.
The app cannot guess arbitrary dormant installations, custom shell functions,
module names or Conda environment names. Multiple matching repository environments
produce an explicit error instead of a random choice. No dependencies are installed.

### Automatic port selection

With `0`, the OS assigns a free port. With a starting port such as `8866`, the app
tries that port and up to 100 subsequent ports before requesting an OS-assigned
port. Windows and server ports can differ; the tunnel connects the selected pair.

When reusing a dashboard, its existing server port is preserved and only a new
Windows port is selected. Ports are displayed in the dashboard terminal.
A race where a selected port becomes occupied causes a safe startup failure via
`--strict-port` or `ExitOnForwardFailure`, rather than switching to another service.
Start again if this happens. Port selection never bypasses MFA or host checking.

## SSH configuration and login nodes

Your trusted `%USERPROFILE%\.ssh\config` may contain usernames, SSH ports,
keys and ProxyJump settings:

```sshconfig
Host delfin-cluster
    HostName specific-login-node.example.org
    User alice
    Port 22
    # IdentityFile ~/.ssh/id_ed25519
    # ProxyJump alice@gateway.example.org
```

Use a specific stable login node: tmux and the dashboard run on that node.
A rotating cluster alias can connect to a different node next time. The app
reports a dashboard record on another node and asks you to connect there.
`delfin-agent where` shows the dashboard location.
Use only trusted SSH configurations; do not add unknown ProxyCommand scripts,
RemoteCommand settings or extra forwardings. Jump hosts have their own login rules.

## Disconnect, reconnect and stop

| Action | Result |
|---|---|
| Close the browser | SSH and dashboard remain running; reopen the local URL while connected |
| Close the SSH window / detach with Ctrl+B, D; Keep off | The dashboard created by this connection and its tmux shell windows are stopped |
| Close the SSH window / detach; Keep on | The local tunnel ends; tmux dashboard remains on the server |
| Start / Reconnect again | Authenticate again, rebuild tunnel and return to an eligible live kept kernel |
| Ctrl+C in the dashboard view | Explicitly stops that dashboard, regardless of Keep |
| `exit` in the working shell | Closes that tmux working window; does not end the dashboard window |

Keep the DELFIN app open for browser access. Closing it disconnects connections
started by that app instance, including pending logins. Existing non-tmux dashboards
continue to depend on their original terminal; the app monitors rather than moves
those processes. A working tmux window is only offered when the dashboard has tmux.
Server restarts or terminated kernels may require a new session.

Stopping a dashboard is not a blanket cancellation of already submitted SLURM or
other independent cluster jobs. Use DELFIN's normal job cancellation controls for
those. Existing agent permissions and cleanup/lifeline behavior remain in effect.

## Update

Close the app window and SSH terminals. With Keep session off, this stops a
new dashboard started by the connection. Enable Keep session before connecting
if you want it to survive disconnecting; Ctrl+C stops the dashboard regardless.
Download the new ZIP, extract
into a new folder, then **double-click `Install.cmd`**. It runs the checks
and updates the installation. The PowerShell commands above remain available.
Program files are replaced, saved connections and their IDs are preserved. No uninstall is needed.
The server installation is updated separately.
Existing saved checkbox choices are preserved. For an existing connection,
uncheck **Working shell (Ctrl+B, then N)** and save/connect to use dashboard-only mode.

## Uninstall

Close the app and its SSH terminal window. If you also want to stop the
server dashboard, stop it in its actual dashboard terminal beforehand.
Uninstalling this Windows app does not stop or uninstall anything on the server.

Use Windows Settings → Apps → Installed apps → **DELFIN SSH Launcher** → Uninstall,
or **double-click `Uninstall.cmd`** in the extracted ZIP folder. Confirm with
**Y** to remove the local app and saved connections; **N** cancels. The result
window stays open. You can also use the installed `Uninstall.cmd` under
`%LOCALAPPDATA%\DELFIN Launcher`.

Alternatively, run this from the extracted ZIP folder:

```powershell
powershell.exe -NoProfile -ExecutionPolicy RemoteSigned -File .\Uninstall.ps1
```

To preserve saved connections for a future installation:

```powershell
powershell.exe -NoProfile -ExecutionPolicy RemoteSigned -File .\Uninstall.ps1 -KeepProfiles
```

The uninstaller also removes stale DELFIN desktop shortcuts, including standard
and OneDrive desktop locations. If Explorer still shows a deleted icon, press F5
on the desktop. You can rerun `Uninstall.cmd` from the extracted ZIP even when
the app folder has already been removed.

Use `-WhatIf` to preview without removing anything. The download/extraction folder
is not deleted. For an older installation without an uninstall script, install
this updated package first; saved connections are preserved.

## Security and verification

- OpenSSH handles passwords, OTP, keys and host verification. Credentials are
  never stored by the app. Changed host keys are not automatically accepted.
  Agent/X11 forwarding, LocalCommand and connection multiplexing are disabled.
- Dashboard and tunnel bind to loopback only. Mandatory dashboard token
  authentication remains enabled. Metadata ownership and private permissions
  are checked before handing off the token over authenticated SSH. A per-connection remote Unix
  socket lives in a random owner-only `/tmp/delfin-win-tunnel-<ID>` directory. It
  relays only to the verified dashboard port and is removed on normal disconnect.
  Abrupt server termination can leave an inert private socket directory behind.
  The tunnel does not expose a SOCKS proxy or allow arbitrary destinations.
- Tokens are used in memory, not stored in Windows token files or transcripts.
  Browsers may retain the login URL and cookies in their own profiles. Protect
  your Windows account and never share token URLs or terminal screenshots.
- Temporary `connection-<ID>.json` files contain only the dashboard SSH process ID
  for checking tunnel ownership; they are removed when that terminal exits.
- New dashboards use `DELFIN_AGENT_SANDBOX=auto` and
  `DELFIN_AGENT_SANDBOX_NETWORK=0`. Existing agent permissions, configuration and
  resumed sessions are preserved. Actual isolation depends on server capabilities
  (bubblewrap, firejail or DELFIN fallback); the app does not certify an existing
  server/agent configuration. A working terminal has your normal server rights.

`Test-Launcher.ps1` checks syntax, connection input, SSH options, port discovery,
browser destinations, native argument quoting and the icon. The **Windows launcher**
CI workflow runs it using Windows PowerShell 5.1. Neither replaces a real login test.
The repository's Python tests cover environment discovery, port conflicts and reuse.
Before release, test password/OTP, both tmux terminal views, both Keep settings, browser startup, occupied ports,
existing dashboards, and disconnect/reconnect to the same kept kernel on Windows.

On errors the terminal stays open. For a new default-port session, inspect it with
`delfin-agent where` to find the actual session name. For an existing session use `delfin-agent where`.
SSH/VPN/MFA failures and cluster policies are never bypassed.

**Validation status:** Windows PowerShell 5.1 CI and repository CI have passed.
Real cluster login has been exercised. Full Keep/reconnect and the latest exit
correction still require end-to-end confirmation; CI does not certify a cluster
configuration.

References: [Microsoft OpenSSH key management](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement),
[OpenSSH configuration options](https://man.openbsd.org/ssh_config).
