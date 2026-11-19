# First steps on a shared agent box

For anyone joining a server or lab box that several people use through one
account. Everything runs on the server; your PC needs no Python and no
DELFIN installation.

`README.md` beside this file is the full reference. This page is the path
from nothing to a working agent, and nothing else.

What you need: Windows 10/11, network access to the server (a VPN if the
server is on an internal network), and the login for the account.

---

## 1. Install

1. Extract the ZIP into a new, empty folder.
2. In the extracted `windows-launcher` folder, double-click **`Install.cmd`**.
   The window stays open so you can read the result. No administrator
   rights are needed.
3. A **DELFIN** shortcut appears on your desktop.

If Windows marks the files as blocked: right-click → Properties → Unblock,
after checking where the package came from. If an antivirus warns: the
package is plain text that compiles and downloads nothing, and the warning
means the file is new and rare rather than harmful. `README.md` has the
detail and the reporting addresses.

## 2. Connect

1. Open DELFIN.
2. In **SSH server / user@host**, enter the account and host you were
   given, for example `delfin@192.0.2.10` or a name from your
   `~/.ssh/config`. This is the only field that has to be filled in.
3. Click **Start / Reconnect** and complete the login.

   **If the server takes no passwords** — many hardened boxes do not — the
   first attempt fails, and that is expected. It has created a key for this
   PC on the way:

   ```
   %USERPROFILE%\.ssh\delfin_launcher_ed25519.pub
   ```

   Send that file (it is the **public** half; the one without `.pub` is
   yours and stays on this PC) to whoever runs the box. They add it in one
   command. Then press **Start / Reconnect** again and you are in, without
   a password from then on.

   To paste a password or one-time code into the black SSH window,
   **right-click**. Ctrl+V does not paste at that prompt; it sends a
   control character and the server answers "Permission denied" for a
   correct password.
4. Your browser opens the dashboard.

Leave the SSH window open while you work. Closing a browser tab does not
disconnect anything.

### Your own area

On a box where the account owner has turned on people mode, your first
connection creates a home of your own: your own sessions and chats, your
own dashboard and token, your own agent workspace. Nothing you do lands in
someone else's session, and nobody else's work lands in yours.

The installation, the tools and the caches are shared, so there is nothing
to install per person.

## 3. Put in your API key

The agent needs a key for the provider it talks to.

**Settings** tab → **API keys**:

1. Leave the provider selected that your box uses (for a KIT Toolbox box,
   `KIT Toolbox (provider kit)`).
2. Paste your key into the field.
3. Click **Store key**.

The line below then reads `KIT_TOOLBOX_API_KEY: set`, and the field clears
itself. The key is never shown back, never written to a log and never put
into a message. It takes effect immediately — no reconnect.

Without a key of your own the agent still works, on whatever key the
account holds; your own key is about your own quota, not about getting
started.

**Know this before you paste one in:** a shared box is one operating-system
account. Everyone who can log into that account can read a key stored
there. The file mode keeps it from *other* accounts, which is a different
boundary. For a key that must be yours alone, use an account of your own.

## 4. Use the agent

**Agent** tab → **New Cycle** → type in the box at the bottom and press
Enter.

A line beginning with `!` is a shell command rather than a message to the
model:

```
! ls
```

It runs through the same checks as the agent's own commands — the same
deny-list, the same approval.

## 5. Stop and come back

**Disconnect** in the DELFIN window. On a box that keeps sessions, yours
stays alive and the next **Start / Reconnect** returns you to it.

`Ctrl+C` in the dashboard window is an explicit stop and ends it even
then.

---

## If something does not work

| What you see | What it is |
|---|---|
| The password is refused although it is right | Ctrl+V in the SSH window. Right-click instead. |
| "Not found" for the environment | The server installation is not where discovery looks. Ask the account owner; do not go hunting. |
| Nothing opens, no error | Are you on the network or VPN the server needs? |
| "Permission denied" and you were given no password | The server takes keys only. Step 2 above: send your `.pub` file. |
| The agent says a key is missing | Section 3 above. |
| Anything else | **Settings** tab, the diagnosis at the bottom. Send that. |

The account owner's side of this — turning on people mode, keeping
sessions, which tabs are shown — is in `README.md` and in
`scripts/pi/README.md` for a Raspberry Pi box.
