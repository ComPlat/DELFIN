#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# delfin-accounts: one Linux account per person (the Pi's accounts variant).
#
# The default Pi setup puts everyone into the one account `delfin` (people
# mode, ~/people/<name>): one Unix user, an overview and no boundary. This
# variant gives each person an account of their own, with its own home
# (0700), API key, sessions, dashboard and SSH keys; the DELFIN installation
# is shared, read only, from $DELFIN_PERSONS_ROOT (pi_bootstrap.sh with
# DELFIN_PI_ACCOUNTS=1 sets it up and installs this script as
# /usr/local/sbin/delfin-accounts, root-owned).
#
# Trust: the agent account (delfin) runs the agent, its CLI backends and MCP
# servers, and owns its build of the installation ($DELFIN_SHARED_ROOT,
# written by deploy_pi.sh). Nothing a person's login sources or executes is
# in that tree: `publish` copies the build into $DELFIN_PERSONS_ROOT, owned
# by root, no group or other write bit anywhere, agent.env and the settings
# defaults rebuilt from an allow-list instead of copied. Root reads nothing
# through a path a person or the agent controls: a person's files are read
# as that person, and a link is refused.
#
# Run by the admin (pi) with sudo:
#   sudo delfin-accounts add NAME [--host HOST] [--hours N] [--no-invite]
#        the account (no password, keys only), SSH allowed, the shared
#        install in its login shell, and an invitation DELFIN-NAME.zip in
#        the admin's ~/invites (as `delfin-people --invite`, into NAME's
#        own account)
#   sudo delfin-accounts invite NAME [--host HOST] [--hours N]
#        a new invitation (another PC, or the first one expired)
#   sudo delfin-accounts lock NAME     no login, processes stopped; files kept
#   sudo delfin-accounts unlock NAME
#   sudo delfin-accounts remove NAME [--keep-home]
#        lock, archive the home to $ARCHIVE_DIR (root only), delete the account
#   sudo delfin-accounts list
#   sudo delfin-accounts publish
#        after every deploy: the agent's build ($DELFIN_SHARED_ROOT) copied
#        into the persons' root-owned installation ($DELFIN_PERSONS_ROOT)
# `delfin-people --add-person NAME` (and --invite-person, --lock-person,
# --unlock-person, --remove-person, --persons) is the same through sudo.
#
# Only accounts in the group $GROUP are persons: lock and remove refuse
# every other account (pi, delfin, system users).
# ============================================================================

CONF="${DELFIN_ACCOUNTS_CONF:-/etc/delfin/accounts.conf}"
GROUP="delfin-people"
SSHD_CONF="/etc/ssh/sshd_config.d/40-delfin-hardening.conf"
PEOPLE_TOOL="/usr/local/lib/delfin/delfin-people"
ARCHIVE_DIR="/var/backups/delfin-people"
SELF="$(readlink -f "$0")"
INVITE_HOURS=48
# The owner the directory above the persons' installation must have.
ROOT_UID=0

die() { echo "delfin-accounts: $*" >&2; exit 1; }

# A Linux user name that is also a launcher person name (remote_bootstrap.sh,
# delfin_people.sh): lower case, digits, _ and -, starting with a letter or
# _ (useradd's default rule, and never a leading dash or digit).
valid_name() { [[ "$1" =~ ^[a-z_][a-z0-9_-]{0,31}$ ]]; }

# The login-shell lines of a person: the persons' installation, root-owned.
# Printed (not root): the launcher finds DELFIN through the active venv
# (remote_bootstrap.sh, VIRTUAL_ENV), the agent gets the box's sandbox
# settings (agent.env) and the dashboard/agent defaults (settings-defaults.json),
# under the person's own settings. Each of these files is written by
# `publish` (root), never by the agent account.
shell_block() {
    local root="$1"
    [[ "$root" =~ ^/[A-Za-z0-9._/-]+$ ]] || die "invalid shared root: $root"
    cat <<EOF
# >>> delfin-accounts: the shared DELFIN installation (read only) >>>
export PATH="\$HOME/.local/bin:\$PATH"
if [ -f '$root/agent.env' ]; then . '$root/agent.env'; fi
if [ -f '$root/.venv/bin/activate' ]; then . '$root/.venv/bin/activate'; fi
if [ -f '$root/settings-defaults.json' ]; then export DELFIN_SETTINGS_DEFAULTS='$root/settings-defaults.json'; fi
# <<< delfin-accounts <<<
EOF
}

# Run as the person, with a clean environment: files in their home are
# written by them, so a link they placed there is never followed by root.
run_as() {
    local who="$1" home
    shift
    home="$(getent passwd "$who" | cut -d: -f6)"
    runuser -u "$who" -- env -i HOME="$home" USER="$who" LOGNAME="$who" \
        PATH=/usr/local/bin:/usr/bin:/bin LANG=C.UTF-8 "$@"
}

# The person's home, as the person (called through run_as).
setup_home() {
    local root="$1" keep="$2"
    umask 077
    mkdir -p "$HOME/.delfin" "$HOME/.config/delfin" "$HOME/.ssh"
    touch "$HOME/.ssh/authorized_keys"
    chmod 700 "$HOME/.ssh"
    chmod 600 "$HOME/.ssh/authorized_keys"
    if ! grep -qs '^# >>> delfin-accounts' "$HOME/.bashrc"; then
        shell_block "$root" >> "$HOME/.bashrc"
    fi
    if [[ "$keep" == 1 && ! -f "$HOME/.config/delfin/launcher-keep.conf" ]]; then
        printf '%s\n' "# Windows launcher: every connection keeps its dashboard (delfin-accounts)." \
            > "$HOME/.config/delfin/launcher-keep.conf"
    fi
    if ! grep -qs '^set -sg escape-time' "$HOME/.tmux.conf"; then
        printf '%s\n' "# Over SSH/VPN: wait long enough for the terminal's answers." \
            "set -sg escape-time 50" >> "$HOME/.tmux.conf"
    fi
}

conf_path() {
    local key="$1" root=""
    [[ -f "$CONF" ]] || die "no $CONF: this box was not set up for one account per person (pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1)"
    root="$(sed -n "s/^$key=//p" "$CONF" | head -n 1)"
    [[ "$root" =~ ^/[A-Za-z0-9._/-]+$ && "$root" != / && "$root" != */ && "$root" != *..* ]] \
        || die "$CONF: no valid $key (run pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1 again)"
    printf '%s' "$root"
}

# The agent's build (written by delfin) and the persons' copy (root's).
shared_root() { conf_path DELFIN_SHARED_ROOT; }
persons_root() { conf_path DELFIN_PERSONS_ROOT; }

deploy_user() {
    sed -n 's/^DELFIN_DEPLOY_USER=//p' "$CONF" 2>/dev/null | head -n 1
}

members() {
    getent group "$GROUP" | cut -d: -f4 | tr ',' ' '
}

is_person() {
    [[ " $(members) " == *" $1 "* ]]
}

require_person() {
    valid_name "$1" || die "NAME: a-z, 0-9, _ and -, starting with a letter or _"
    id "$1" >/dev/null 2>&1 || die "no account $1"
    is_person "$1" || die "$1 is not a DELFIN person (not in group $GROUP); nothing changed"
}

# Locked = the account's expiry date is reached (lock sets day 1). sshd's PAM
# account check refuses an expired account even with a valid key; a locked
# password ('!') alone does not stop a key login.
locked() {
    local expire today
    expire="$(getent shadow "$1" | cut -d: -f8)"
    today=$(( $(date +%s) / 86400 ))
    [[ -n "$expire" && "$expire" -le "$today" ]]
}

active_people() {
    local u out=""
    for u in $(members); do
        locked "$u" || out+=" $u"
    done
    printf '%s' "$out"
}

# sshd's allow-list (pi_harden_host.sh writes `AllowUsers pi delfin`): every
# person who is not locked is added, every other person taken out; the
# accounts that are not persons stay as they are. Checked with sshd -t before
# the reload, and put back when the check fails. The line is edited in that
# drop-in, not added in one of its own: sshd keeps the first AllowUsers it
# reads, so a second file's names were ignored (measured, OpenSSH 10.0).
sync_ssh() {
    local line word base="" person backup
    [[ -f "$SSHD_CONF" ]] || { echo "note: no $SSHD_CONF (host not hardened): no SSH allow-list to keep"; return 0; }
    line="$(grep -m1 -E '^[[:space:]]*AllowUsers[[:space:]]' "$SSHD_CONF" || true)"
    [[ -n "$line" ]] || return 0
    # A name without an account goes too: a removed person is no member any
    # more and would otherwise stay in the list as if it were pi or delfin.
    for word in ${line#*AllowUsers}; do
        is_person "$word" && continue
        id "$word" >/dev/null 2>&1 || continue
        base+=" $word"
    done
    person="$(active_people)"
    backup="$(mktemp "$SSHD_CONF.bak.XXXXXX")"
    cp -p "$SSHD_CONF" "$backup"
    awk -v new="AllowUsers${base}${person}" \
        'BEGIN{done=0} /^[[:space:]]*AllowUsers[[:space:]]/ && !done {print new; done=1; next} {print}' \
        "$backup" > "$SSHD_CONF.tmp"
    chmod 644 "$SSHD_CONF.tmp"
    mv -f "$SSHD_CONF.tmp" "$SSHD_CONF"
    if ! sshd -t; then
        mv -f "$backup" "$SSHD_CONF"
        die "sshd refused the new allow-list; the previous one is back"
    fi
    rm -f "$backup"
    systemctl reload ssh 2>/dev/null || systemctl reload sshd 2>/dev/null \
        || echo "warning: reload sshd yourself (systemctl reload ssh)" >&2
    echo "==> SSH allow-list: AllowUsers${base}${person}"
}

admin_invites() {
    local admin="${SUDO_USER:-root}" home
    home="$(getent passwd "$admin" | cut -d: -f6)"
    [[ -n "$home" ]] || { admin=root; home=/root; }
    install -d -m 700 -o "$admin" -g "$(id -gn "$admin")" "$home/invites"
    printf '%s %s' "$admin" "$home/invites"
}

invite() {
    local name="$1" host="$2" hours="$3" root package home out admin dir valid
    require_person "$name"
    locked "$name" && die "$name is locked: unlock first (delfin-accounts unlock $name)"
    # The persons' copy, root's: the build's own zip is the agent's to change.
    root="$(persons_root)"
    package="$root/src/tools/delfin-windows-launcher.zip"
    [[ -f "$package" && ! -L "$package" ]] || die "launcher package not found: $package (deploy DELFIN, then: sudo delfin-accounts publish)"
    [[ -x "$PEOPLE_TOOL" ]] || die "$PEOPLE_TOOL missing: run pi_bootstrap.sh with DELFIN_PI_ACCOUNTS=1 again"
    home="$(getent passwd "$name" | cut -d: -f6)"
    # delfin-people --invite, run as the person: the invitation's key and
    # its forced --enroll land in the person's own authorized_keys, the
    # server in the file is NAME@HOST, and the enroll command is the
    # root-owned copy, which the person cannot change.
    local args=(--invite "$name" --hours "$hours")
    [[ -n "$host" ]] && args+=(--host "$host")
    if ! out="$(run_as "$name" DELFIN_LAUNCHER_ZIP="$package" "$PEOPLE_TOOL" "${args[@]}" 2>&1)"; then
        printf '%s\n' "$out" >&2
        die "invitation for $name not made"
    fi
    read -r admin dir <<<"$(admin_invites)"
    # The invitation lies in the person's home: read as the person (cat
    # opens only what they may read), and a link there is refused outright,
    # so root never copies a file through a path the person chose.
    local src="$home/invites/DELFIN-$name.zip" tmp
    if [[ -L "$src" || -L "$home/invites" ]]; then
        run_as "$name" rm -f -- "$src" || true
        die "invitation for $name: $src is a link; refused, nothing copied"
    fi
    tmp="$(umask 077 && mktemp "$dir/.DELFIN-$name.XXXXXX")"
    if ! run_as "$name" cat -- "$src" > "$tmp"; then
        rm -f -- "$tmp"
        die "invitation for $name: $src not readable as $name"
    fi
    chown "$admin:$(id -gn "$admin")" "$tmp"
    chmod 600 "$tmp"
    mv -f -- "$tmp" "$dir/DELFIN-$name.zip"
    run_as "$name" rm -f -- "$src" || true
    run_as "$name" rmdir -- "$home/invites" 2>/dev/null || true
    valid="$(grep -m1 -o 'Valid once, until [^ ]*' <<<"$out" || true)"
    echo "Invitation for $name: $dir/DELFIN-$name.zip"
    echo "  $(grep -m1 -o 'Server: [^ .]*@[^ ]*' <<<"$out" | sed 's/\.$//'). ${valid:-Valid once}."
    echo "  Give the file to $name like a temporary password: extract it, double-click"
    echo "  Install.cmd, press Start in DELFIN. It leads into $name's own account."
    echo "  Withdraw it: sudo delfin-accounts lock $name (or a new one replaces nothing: each is used once)."
}

add() {
    local name="$1" host="$2" hours="$3" no_invite="$4" root home keep=0 deployer
    valid_name "$name" || die "NAME: a-z, 0-9, _ and -, starting with a letter or _"
    root="$(persons_root)"
    deployer="$(deploy_user)"
    local existed=0
    if id "$name" >/dev/null 2>&1; then
        is_person "$name" || die "an account $name exists and is not a DELFIN person; nothing changed"
        existed=1
    fi
    getent group "$GROUP" >/dev/null || groupadd --system "$GROUP"
    # No password at all (the shadow entry stays '!'): keys only, and the
    # first key comes through the invitation. Run again for an existing
    # person, add only re-applies the setup below (nothing is replaced, a
    # lock stays) and makes no new invitation.
    if [[ "$existed" == 0 ]]; then
        useradd --create-home --shell /bin/bash --user-group --groups "$GROUP" \
            --comment "DELFIN person" "$name"
    fi
    home="$(getent passwd "$name" | cut -d: -f6)"
    chmod 700 "$home"
    # The box's default for the launcher (pi_install.sh, DELFIN_PI_KEEP) is
    # the deploy account's: a person gets the same.
    if [[ -n "$deployer" ]]; then
        local dhome
        dhome="$(getent passwd "$deployer" | cut -d: -f6)"
        [[ -n "$dhome" && -f "$dhome/.config/delfin/launcher-keep.conf" ]] && keep=1
    fi
    run_as "$name" bash "$SELF" _setup-home "$root" "$keep"
    echo "==> account $name: $home (0700), group $GROUP, no password login"
    [[ -x "$root/.venv/bin/python" ]] \
        || echo "note: $root/.venv has no python yet: deploy DELFIN (deploy_pi.sh), then sudo delfin-accounts publish, before $name starts" >&2
    sync_ssh
    if [[ "$existed" == 1 ]]; then
        echo "$name was already a person: setup checked, no new invitation (delfin-accounts invite $name)"
    elif [[ "$no_invite" == 0 ]]; then
        invite "$name" "$host" "$hours"
    fi
}

stop_processes() {
    loginctl terminate-user "$1" >/dev/null 2>&1 || true
    pkill -KILL -u "$1" >/dev/null 2>&1 || true
}

lock() {
    require_person "$1"
    usermod --expiredate 1 "$1"
    sync_ssh
    stop_processes "$1"
    echo "==> $1 locked: no login (SSH, keys included), processes stopped, files kept."
    echo "    Undo: sudo delfin-accounts unlock $1"
}

unlock() {
    require_person "$1"
    chage --expiredate -1 "$1"
    sync_ssh
    echo "==> $1 unlocked"
}

remove() {
    local name="$1" keep_home="$2" home archive
    require_person "$name"
    home="$(getent passwd "$name" | cut -d: -f6)"
    # Out of sshd's list while still a member, then the account goes.
    usermod --expiredate 1 "$name"
    sync_ssh
    stop_processes "$name"
    if [[ "$keep_home" == 1 ]]; then
        userdel "$name"
        echo "==> account $name removed; home kept: $home (owner is now a bare uid; sudo chown -R or rm -rf it)"
    else
        install -d -m 700 -o root -g root "$ARCHIVE_DIR"
        archive="$ARCHIVE_DIR/$name-$(date +%Y%m%d-%H%M%S).tar.gz"
        (umask 077 && tar --one-file-system -C "$(dirname "$home")" -czf "$archive" "$(basename "$home")")
        # Exit 12 is "mail spool not found": the account and home are gone.
        userdel --remove "$name" || true
        id "$name" >/dev/null 2>&1 && die "userdel did not remove $name"
        echo "==> account $name removed; its home is archived (root only): $archive"
    fi
}

list() {
    local u home state keys pending procs
    printf '%-16s %-8s %-5s %-8s %s\n' PERSON STATE KEYS INVITE PROCS
    for u in $(members); do
        home="$(getent passwd "$u" | cut -d: -f6)"
        if locked "$u"; then state=locked; else state=active; fi
        # Read as the person; a link is not followed at all.
        if [[ -L "$home/.ssh" || -L "$home/.ssh/authorized_keys" ]]; then
            keys="?" pending="?"
        else
            keys="$(run_as "$u" grep -c '^restrict,port-forwarding,pty ' "$home/.ssh/authorized_keys" 2>/dev/null || true)"
            if run_as "$u" grep -qs ' delfin-invite-' "$home/.ssh/authorized_keys"; then pending=pending; else pending=-; fi
        fi
        procs="$(pgrep -c -u "$u" 2>/dev/null || true)"
        printf '%-16s %-8s %-5s %-8s %s\n' "$u" "$state" "${keys:-0}" "$pending" "${procs:-0}"
    done
}

# -- publish: the agent's build -> the persons' root-owned installation -----
# The parts persons run; the agent's workspace (work/) and anything else in
# the build stay out.
PUBLISH_PARTS="src .venv python ms-playwright"

# agent.env as persons source it: only these variables, values of a plain
# character set, rebuilt line by line (never sourced or copied as root).
# Anything else -- a secret the admin put there, a command -- is dropped and
# named. No secret is among them: the file is readable by every person.
env_allowed() {
    local var="$1" val="$2"
    case "$var" in
        DELFIN_AGENT_SANDBOX) [[ "$val" =~ ^[a-z]+$ && "$val" != off ]] ;;
        DELFIN_AGENT_SANDBOX_TIMEOUT) [[ "$val" =~ ^[0-9]{1,6}$ ]] ;;
        DELFIN_MCP_ISOLATION) [[ "$val" =~ ^[a-z]+$ && "$val" != off && "$val" != none ]] ;;
        PLAYWRIGHT_BROWSERS_PATH) [[ "$val" == "$PERSONS/ms-playwright" ]] ;;
        *) return 1 ;;
    esac
}

publish_env() {
    local src="$1" dst="$2" line var val
    printf '%s\n' "# Written by delfin-accounts publish (root) from the agent's agent.env:" \
        "# allow-listed variables only. Sourced by every person's login." > "$dst"
    [[ -f "$src" && ! -L "$src" ]] || return 0
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" == \#* ]] && continue
        if [[ "$line" =~ ^export\ ([A-Z_][A-Z0-9_]*)=([A-Za-z0-9_./:-]*)$ ]]; then
            var="${BASH_REMATCH[1]}"
            val="${BASH_REMATCH[2]}"
            [[ "$val" == "$BUILD" || "$val" == "$BUILD"/* ]] && val="$PERSONS${val#"$BUILD"}"
            if env_allowed "$var" "$val"; then
                printf 'export %s=%s\n' "$var" "$val" >> "$dst"
                continue
            fi
            echo "publish: agent.env: $var not passed on to persons" >&2
        else
            echo "publish: agent.env: a line that is no plain export not passed on: ${line:0:40}" >&2
        fi
    done < "$src"
}

# settings-defaults.json as persons get it: the agent/tab keys pi_install.sh
# writes, strings only, and never a permission profile that stops asking
# (all_free and its older names). Parsed by the system's Python, isolated.
publish_defaults() {
    local src="$1" dst="$2"
    [[ -f "$src" && ! -L "$src" ]] || return 0
    /usr/bin/env -i PATH=/usr/bin:/bin python3 -I -S - "$src" "$dst" <<'EOF'
import json, re, sys
src, dst = sys.argv[1], sys.argv[2]
try:
    with open(src, encoding="utf-8") as fh:
        data = json.load(fh)
except (OSError, ValueError) as exc:
    print(f"publish: settings-defaults.json not read ({exc}); persons get none", file=sys.stderr)
    data = {}
plain = re.compile(r"[A-Za-z0-9._:/@+-]{1,120}")
safe_profiles = {"plan", "ask_all", "repo_free"}
out = {}
agent = data.get("agent") if isinstance(data, dict) else None
if isinstance(agent, dict):
    kept = {}
    for key in ("mode", "provider", "model", "permission_profile", "effort"):
        val = agent.get(key)
        if not isinstance(val, str) or not plain.fullmatch(val):
            continue
        if key == "permission_profile" and val not in safe_profiles:
            print(f"publish: permission_profile {val!r} not passed on to persons", file=sys.stderr)
            continue
        kept[key] = val
    if kept:
        out["agent"] = kept
ui = data.get("ui") if isinstance(data, dict) else None
tabs = ui.get("tabs") if isinstance(ui, dict) else None
if isinstance(tabs, dict):
    clean = {k: [t for t in tabs[k] if isinstance(t, str) and plain.fullmatch(t)]
             for k in ("order", "hidden") if isinstance(tabs.get(k), list)}
    if clean:
        out["ui"] = {"tabs": clean}
with open(dst, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
    fh.write("\n")
EOF
}

publish() {
    local parent new old part link target abs real bad sed_from
    BUILD="$(shared_root)"
    PERSONS="$(persons_root)"
    [[ "$PERSONS" != "$BUILD" && "$PERSONS" != "$BUILD"/* && "$BUILD" != "$PERSONS"/* ]] \
        || die "$CONF: DELFIN_PERSONS_ROOT must be outside DELFIN_SHARED_ROOT"
    [[ -d "$BUILD/src" && ! -L "$BUILD" && -f "$BUILD/.venv/bin/activate" ]] \
        || die "no deployed installation in $BUILD (deploy_pi.sh first)"
    # The copy and its swap happen beside the persons' root: that directory
    # must be root's alone, or the agent could swap a path under root's feet.
    parent="$(dirname "$PERSONS")"
    [[ -d "$parent" && ! -L "$parent" && "$(stat -c %u "$parent")" == "$ROOT_UID" ]] \
        && ! [[ "$(stat -c %A "$parent")" =~ ^d....w|^d.......w ]] \
        || die "$parent must be a directory only root writes (owner uid $ROOT_UID, no group/other write)"
    new="$PERSONS.new"
    old="$PERSONS.old"
    rm -rf -- "$new" "$old"
    (umask 077 && mkdir -- "$new")
    for part in $PUBLISH_PARTS; do
        [[ -e "$BUILD/$part" || -L "$BUILD/$part" ]] || continue
        # -P: links are copied as links, never followed; no owner kept.
        cp -R -P --preserve=mode,timestamps,links -- "$BUILD/$part" "$new/$part"
    done
    bad="$(find "$new" ! -type f ! -type d ! -type l -print -quit)"
    [[ -z "$bad" ]] || { rm -rf -- "$new"; die "publish: ${bad#"$new"/} is not a file, directory or link; nothing published"; }
    # Text that names the build (venv scripts, pyvenv.cfg, the editable
    # .pth) names the persons' copy instead. Compiled .pyc files keep the
    # build path as a file name for tracebacks only.
    sed_from="$(printf '%s' "$BUILD" | sed 's/[.[\*^$/]/\\&/g')"
    grep -rlIZF -- "$BUILD" "$new" 2>/dev/null \
        | xargs -0 -r sed -i -E "s/${sed_from}([^A-Za-z0-9_.-]|\$)/$(printf '%s' "$PERSONS" | sed 's/[\/&]/\\&/g')\\1/g"
    # Links: one into the build is moved to the copy; afterwards every link
    # must stay inside the copy or point into the system's own /usr, /etc.
    while IFS= read -r -d '' link; do
        target="$(readlink -- "$link")"
        if [[ "$target" == "$BUILD" || "$target" == "$BUILD"/* ]]; then
            target="$PERSONS${target#"$BUILD"}"
            ln -sfn -- "$target" "$link"
        fi
        if [[ "$target" == /* ]]; then
            abs="$target"
            [[ "$abs" == "$PERSONS" || "$abs" == "$PERSONS"/* ]] && abs="$new${abs#"$PERSONS"}"
        else
            abs="$(dirname -- "$link")/$target"
        fi
        real="$(realpath -m -s -- "$abs")"
        case "$real" in
            "$new"|"$new"/*|/usr/*|/etc/*) ;;
            *) rm -rf -- "$new"
               die "publish: ${link#"$new"/} points outside the installation ($target); nothing published" ;;
        esac
    done < <(find "$new" -type l -print0)
    publish_env "$BUILD/agent.env" "$new/agent.env"
    publish_defaults "$BUILD/settings-defaults.json" "$new/settings-defaults.json"
    # Root's, readable, writable by nobody else; no set-id bit survives.
    find "$new" -type f -perm /6000 -exec chmod ug-s {} +
    find "$new" \( -type f -o -type d \) -exec chmod go-w,a+rX {} +
    chown -R -h root:root "$new"
    chown root:delfin-people "$new"
    chmod 750 "$new"
    # The swap: the previous copy stays until the next publish, so a
    # running person's process keeps its files.
    if [[ -e "$PERSONS" || -L "$PERSONS" ]]; then mv -T -- "$PERSONS" "$old"; fi
    mv -T -- "$new" "$PERSONS"
    echo "==> published for persons: $PERSONS (root-owned; from $BUILD)"
    echo "    A person's dashboard picks it up on its next start."
}

usage() {
    sed -n '/^# Run by the admin/,/^# every other account/p' "$SELF" | sed 's/^# \{0,1\}//' >&2
    exit 2
}

cmd="${1:-}"
[[ $# -gt 0 ]] && shift
case "$cmd" in
    shell-block) shell_block "${1:-}"; exit 0 ;;
    _setup-home) setup_home "${1:?}" "${2:-0}"; exit 0 ;;
    add|invite|lock|unlock|remove|list|publish) ;;
    *) usage ;;
esac
name="" host="" hours="$INVITE_HOURS" no_invite=0 keep_home=0
if [[ "$cmd" != list && "$cmd" != publish ]]; then
    name="${1:-}"
    [[ $# -gt 0 ]] && shift
    [[ -n "$name" ]] || usage
    valid_name "$name" || die "NAME: a-z, 0-9, _ and -, starting with a letter or _"
fi
while [[ $# -gt 0 ]]; do
    case "$1" in
        --host) host="${2:-}"; shift 2 || usage ;;
        --hours) hours="${2:-}"; shift 2 || usage ;;
        --no-invite) no_invite=1; shift ;;
        --keep-home) keep_home=1; shift ;;
        *) echo "unknown option: $1" >&2; usage ;;
    esac
done
[[ "$hours" =~ ^[0-9]{1,4}$ && "$hours" -ge 1 ]] || die "--hours: a whole number of hours"
[[ -z "$host" || "$host" =~ ^[A-Za-z0-9][A-Za-z0-9.-]{0,200}$ ]] || die "--host: a host name or IP"
[[ $EUID -eq 0 ]] || die "run with sudo"

case "$cmd" in
    add) add "$name" "$host" "$hours" "$no_invite" ;;
    invite) invite "$name" "$host" "$hours" ;;
    lock) lock "$name" ;;
    unlock) unlock "$name" ;;
    remove) remove "$name" "$keep_home" ;;
    list) list ;;
    publish) publish ;;
esac
