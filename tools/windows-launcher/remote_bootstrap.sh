# Runs inside the user's interactive login shell, sourced from the copy the app
# keeps in ~/.cache/delfin-launcher/<version>/ (DELFIN.ps1, Get-RemoteCommand, Copy-ServerScripts):
#   bash -c <check> && bash -ilc '. ~/.cache/delfin-launcher/<version>/remote_bootstrap.sh' delfin ARGS
# Sourcing (not a child bash) keeps what the login shell set up, as before.
# Arguments are data, never eval'ed. $4 is the path of remote_launcher.py.
set -e
preferred=$1
directory=$2
hint=$3
launcher_script=$4
role=$5
tunnel_path=${6:-}
keep=${7:-0}
working=${8:-0}
person=${9:-}
login_key=${10:-}
# One line, nothing else: grep below checks line by line, so a second line
# would pass unchecked and land in authorized_keys without the restriction.
case "$login_key" in *$'\n'*|*$'\r'*) login_key='' ;; esac
# Remember login: the app sends this PC's public key (one ed25519 key per
# Windows user) on every connection; the first one, still with the password,
# adds it to the account's authorized_keys, the next ones log in with it.
# The key line is restricted to what the launcher uses (port forwarding for
# the dashboard tunnel, a pty), so it opens no agent, X11 or ~/.ssh/rc.
# ~/.config/delfin/launcher-no-keys.conf makes a server refuse it.
if [ -n "$login_key" ] && [ ! -f "$HOME/.config/delfin/launcher-no-keys.conf" ]; then
    if printf '%s\n' "$login_key" | grep -Eq '^ssh-ed25519 [A-Za-z0-9+/]+={0,2}( [A-Za-z0-9@._-]{1,64})?$'; then
        key_body=$(printf '%s' "$login_key" | cut -d' ' -f2)
        if ! grep -qsF " $key_body" "$HOME/.ssh/authorized_keys"; then
            (umask 077 && mkdir -p "$HOME/.ssh" \
                && printf 'restrict,port-forwarding,pty %s\n' "$login_key" >> "$HOME/.ssh/authorized_keys") \
                && printf '%s\n' 'DELFIN: Login key of this PC registered; the next connection needs no password.'
        fi
    fi
fi
# The account owner may make Keep session this server's default
# (~/.config/delfin/launcher-keep.conf, e.g. a shared agent box): every
# connection then keeps its dashboard after a disconnect, whatever the app's
# checkbox says. Ctrl+C in the SSH window still stops it.
if [ -f "$HOME/.config/delfin/launcher-keep.conf" ]; then keep=1; fi
# `python -c` puts the current directory first on sys.path: a plain folder
# named "delfin" there (an install root, a repository root) then shadows the
# real package as an empty namespace package. Every probe and start drops it.
no_cwd='import sys; sys.path[:] = [p for p in sys.path if p not in ("", ".")]'
expand_path() {
    case "$1" in
        '~') printf '%s' "$HOME" ;;
        '~/'*) printf '%s/%s' "$HOME" "${1:2}" ;;
        *) printf '%s' "$1" ;;
    esac
}
directory=$(expand_path "$directory")
hint=$(expand_path "$hint")
# Before people mode below changes HOME: the copy lives in the account's home.
launcher_script=$(expand_path "$launcher_script")
if [ ! -f "$launcher_script" ]; then
    printf '%s\n' 'DELFIN: The launcher script is missing on the server. Start the connection again from the app.' >&2
    exit 1
fi
if [ -n "$directory" ]; then cd -- "$directory"; else cd -- "$HOME"; fi
case "$hint" in ''|/*) ;; *) hint="$PWD/$hint" ;; esac
python=''
selected_repo=''
probe_repo=''
repository_for_python() {
    local env_root parent
    env_root=$(dirname "$(dirname "$1")")
    parent=$(dirname "$env_root")
    if [ -f "$env_root/delfin/cli_voila.py" ]; then printf '%s' "$env_root"
    elif [ -f "$parent/delfin/cli_voila.py" ]; then printf '%s' "$parent"
    elif [ -n "$probe_repo" ]; then printf '%s' "$probe_repo"
    elif [ -f "$PWD/delfin/cli_voila.py" ]; then printf '%s' "$PWD"
    fi
}
# A Python that has DELFIN but not the dashboard and agent modules this
# launcher starts: remembered, so the final message can name it instead of
# saying "Not found" about an installation the person can see.
too_old=''
usable() {
    local repo
    [ -x "$1" ] || return 1
    repo=$(repository_for_python "$1")
    if PYTHONPATH="${repo:+$repo:}${PYTHONPATH:-}" "$1" -c "$no_cwd"'; import delfin.cli_voila, delfin.agent.where, delfin.dashboard.session' >/dev/null 2>&1; then
        return 0
    fi
    if [ -z "$too_old" ] && PYTHONPATH="${repo:+$repo:}${PYTHONPATH:-}" "$1" -c "$no_cwd"'; import delfin.cli_voila' >/dev/null 2>&1; then
        too_old=$1
    fi
    return 1
}
# The dashboard is served by voila, a dependency of DELFIN itself. An
# environment that imports DELFIN but not voila (a partial install) passes the
# checks above and then dies inside tmux, out of sight.
serves_dashboard() {
    local repo
    repo=$(repository_for_python "$1")
    PYTHONPATH="${repo:+$repo:}${PYTHONPATH:-}" "$1" -c "$no_cwd"'; import voila' >/dev/null 2>&1
}
# The Python behind a console script: a site-wide install often puts only a
# link into /usr/local/bin (-> /opt/<venv>/bin/delfin-voila), where no python
# lies beside it. The link is followed, then the script's #! line is read.
entry_point_pythons() {
    local real shebang
    real=$(readlink -f -- "$1" 2>/dev/null || printf '%s' "$1")
    printf '%s\n' "$(dirname "$real")/python" "$(dirname "$real")/python3"
    shebang=$(head -n 1 -- "$real" 2>/dev/null | sed -n 's/^#![[:space:]]*\([^[:space:]]*\).*$/\1/p')
    case "$shebang" in */python*) printf '%s\n' "$shebang" ;; esac
}
choose() {
    if [ -z "$python" ] && usable "$1"; then
        python=$1
        selected_repo=$(repository_for_python "$1")
    fi
}
scan_root() {
    local root=$1 candidate
    probe_repo=''
    if [ -f "$root/delfin/cli_voila.py" ]; then probe_repo=$root; fi
    local matches=()
    for candidate in "$root/bin/python" "$root/bin/python3" "$root/.venv/bin/python" "$root/venv/bin/python" "$root/env/bin/python" "$root/.env/bin/python"; do
        if usable "$candidate"; then
            # python and python3 inside one environment are the same choice.
            if [ ${#matches[@]} -eq 0 ] || [ "$(dirname "$candidate")" != "$(dirname "${matches[0]}")" ]; then matches+=("$candidate"); fi
        fi
    done
    if [ ${#matches[@]} -gt 1 ]; then
        printf '%s\n' 'DELFIN: Multiple environments found. Set an explicit DELFIN location in advanced settings.' >&2
        exit 1
    fi
    if [ ${#matches[@]} -eq 1 ]; then choose "${matches[0]}"; fi
    probe_repo=''
}
if [ -n "$hint" ]; then
    if [ -d "$hint" ]; then
        if [ -f "$hint/delfin/cli_voila.py" ]; then export PYTHONPATH="$hint${PYTHONPATH:+:$PYTHONPATH}"; fi
        scan_root "$hint"
        if [ -z "$python" ]; then
            candidate=$(type -P python || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
            candidate=$(type -P python3 || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
        fi
        if [ -z "$python" ]; then
            printf '%s\n' 'DELFIN: No usable environment found in the specified directory.' >&2; exit 1
        fi
    else
        # Explicit Python, or the environment containing a console entry point.
        case "$(basename "$hint")" in python|python3|python3.*) choose "$hint" ;; esac
        choose "$(dirname "$hint")/python"
        choose "$(dirname "$hint")/python3"
        if [ -z "$python" ]; then printf '%s\n' 'DELFIN: The specified DELFIN location does not identify a usable Python environment.' >&2; exit 1; fi
    fi
else
    # The checkout install.sh set up comes first: the login shell exports
    # DELFIN_REPO from ~/.delfin_env.sh, or the file is read (as data, the
    # one line install.sh writes, never sourced). An older venv activated
    # from ~/.bashrc, or a second install in a conventional place, used to
    # win or stop the start ("Multiple installations found").
    recorded=${DELFIN_REPO:-}
    if [ -z "$recorded" ] && [ -f "$HOME/.delfin_env.sh" ]; then
        recorded=$(sed -n "s/^export DELFIN_REPO=//p" "$HOME/.delfin_env.sh" | head -n 1)
        recorded=${recorded#\'}
        recorded=${recorded%\'}
    fi
    if [ -n "$recorded" ] && [ -f "$recorded/delfin/cli_voila.py" ]; then
        scan_root "$recorded"
        if [ -n "$python" ]; then
            printf 'DELFIN: Using the install set up by install.sh: %s\n' "$recorded"
            for other in "${VIRTUAL_ENV:-}" "${CONDA_PREFIX:-}"; do
                if [ -n "$other" ] && [ "$(dirname "$(dirname "$python")")" != "$other" ]; then
                    printf 'DELFIN: Another environment was active in the login shell and is not used: %s (a line in ~/.bashrc activates it).\n' "$other"
                fi
            done
        fi
    fi
    if [ -n "${VIRTUAL_ENV:-}" ]; then choose "$VIRTUAL_ENV/bin/python"; fi
    if [ -n "${CONDA_PREFIX:-}" ]; then choose "$CONDA_PREFIX/bin/python"; fi
    if [ -z "$python" ]; then
        root=$PWD
        while [ "$root" != / ] && [ -z "$python" ]; do
            scan_root "$root"
            root=$(dirname "$root")
        done
    fi
    if [ -z "$python" ]; then
        # Bounded conventional locations; never scan other users or the filesystem.
        found_python=()
        found_repo=()
        for conventional in "$HOME/software/delfin" "$HOME/software/DELFIN" "$HOME/delfin" "$HOME/DELFIN"; do
            if [ -f "$conventional/delfin/cli_voila.py" ]; then
                scan_root "$conventional"
                if [ -n "$python" ]; then
                    found_python+=("$python")
                    found_repo+=("$selected_repo")
                    python=''
                    selected_repo=''
                fi
            fi
        done
        # A Conda or micromamba environment named "delfin" or "delfin-voila"
        # (the names DELFIN's documentation uses), under the usual roots of
        # this account only. Without it a server whose DELFIN lives in such an
        # environment, never activated in ~/.bashrc, needed the location typed
        # in by hand. One that can serve the dashboard (voila) wins over one
        # that cannot; the latter is only taken when it is all there is, and
        # the dashboard check below then says what it lacks.
        seen_env=''
        partial_python=()
        partial_repo=()
        for env_root in "${MAMBA_ROOT_PREFIX:-}" "$HOME/micromamba" "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
            [ -n "$env_root" ] || continue
            for env_name in delfin delfin-voila; do
                [ -x "$env_root/envs/$env_name/bin/python" ] || continue
                real_env=$(cd "$env_root/envs/$env_name" && pwd -P)
                case " $seen_env " in *" $real_env "*) continue ;; esac
                seen_env="$seen_env $real_env"
                if usable "$real_env/bin/python"; then
                    if serves_dashboard "$real_env/bin/python"; then
                        found_python+=("$real_env/bin/python")
                        found_repo+=("$(repository_for_python "$real_env/bin/python")")
                    else
                        partial_python+=("$real_env/bin/python")
                        partial_repo+=("$(repository_for_python "$real_env/bin/python")")
                    fi
                fi
            done
        done
        if [ ${#found_python[@]} -eq 0 ] && [ ${#partial_python[@]} -gt 0 ]; then
            found_python=("${partial_python[@]}")
            found_repo=("${partial_repo[@]}")
        fi
        if [ ${#found_python[@]} -gt 1 ]; then
            printf '%s\n' 'DELFIN: Several DELFIN installations were found on the server:' >&2
            for index in "${!found_python[@]}"; do printf '  %s\n' "${found_repo[$index]:-${found_python[$index]}}" >&2; done
            printf '%s\n' 'Choose one under Advanced settings > DELFIN location, or run install.sh in the one to use (the launcher then takes that one).' >&2
            exit 1
        fi
        if [ ${#found_python[@]} -eq 1 ]; then
            python=${found_python[0]}
            selected_repo=${found_repo[0]}
        fi
    fi
    if [ -z "$python" ]; then
        candidate=$(type -P delfin-voila || true)
        if [ -n "$candidate" ]; then
            choose "$(dirname "$candidate")/python"; choose "$(dirname "$candidate")/python3"
            while IFS= read -r entry; do choose "$entry"; done < <(entry_point_pythons "$candidate")
        fi
        candidate=$(type -P python || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
        candidate=$(type -P python3 || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
    fi
fi
if [ -z "$python" ]; then
    if [ -n "$too_old" ]; then
        printf 'DELFIN: The only installation found is too old for this launcher (it has no agent dashboard): %s\n' "$too_old" >&2
        printf '%s\n' 'Update that installation, or set Advanced settings > DELFIN location to a newer checkout or environment.' >&2
        exit 1
    fi
    printf '%s\n' 'DELFIN: Not found in the login environment or repository. Set a working directory or DELFIN environment location.' >&2
    exit 1
fi
if [ -n "$selected_repo" ]; then
    export PYTHONPATH="$selected_repo${PYTHONPATH:+:$PYTHONPATH}"
    if [ -z "$directory" ]; then cd -- "$selected_repo"; fi
fi
python_dir=$(dirname "$python")
export PATH="$python_dir:$PATH"
environment_root=$(cd "$python_dir/.." && pwd -P)
active_root=''
if [ -n "${VIRTUAL_ENV:-}" ] && [ -d "$VIRTUAL_ENV" ]; then active_root=$(cd "$VIRTUAL_ENV" && pwd -P); fi
if [ -f "$python_dir/activate" ] && [ "$active_root" != "$environment_root" ]; then
    source "$python_dir/activate"
fi
export PS1
printf 'DELFIN: Environment found: %s\n' "$python"
# Several people on one server account (opt-in by the account owner): each
# Windows user gets a home of their own under ~/people, so sessions,
# dashboard, token, agent service and workspace do not collide. Shared by
# link: the account's credentials, tools, caches, git and SSH setup, and the
# installation itself. Not a security boundary between them: they are one
# Unix user.
if [ -n "$person" ] && [ -f "$HOME/.config/delfin/launcher-people.conf" ]; then
    case "$person" in
        *[!a-z0-9_-]*|-*) printf '%s\n' 'DELFIN: Invalid person name.' >&2; exit 1 ;;
    esac
    shared=$HOME
    person_home="$shared/people/$person"
    (umask 077 && mkdir -p "$person_home/.delfin" "$person_home/.run" "$person_home/agent_workspace")
    # DELFIN's settings are ~/.delfin_settings.json (user_settings.SETTINGS_FILE_NAME),
    # in two layers: the account's file holds the defaults (pi_install.sh:
    # pi-agent.conf, pi-tabs.conf) and is read under the person's own file
    # (DELFIN_SETTINGS_DEFAULTS), which keeps only what that person changed.
    # A person file that is still a link to the account's becomes an empty
    # file of its own: its content was the account's, so nothing is lost. A
    # copy kept from before (*.before-shared-*) is not applied again.
    settings_own="$person_home/.delfin_settings.json"
    if [ -L "$settings_own" ] && { [ "$settings_own" -ef "$shared/.delfin_settings.json" ] \
            || [ "$(readlink -- "$settings_own")" = "$shared/.delfin_settings.json" ]; }; then
        rm -f -- "$settings_own"
    fi
    if [ ! -e "$settings_own" ] && [ ! -L "$settings_own" ]; then
        (umask 077 && printf '{}\n' > "$settings_own")
    fi
    for item in delfin .local .cache .ssh .gitconfig .tmux.conf .delfin/credentials.json .delfin/tools .delfin/doc_index.json; do
        if [ -e "$shared/$item" ] && [ ! -e "$person_home/$item" ] && [ ! -L "$person_home/$item" ]; then
            ln -s "$shared/$item" "$person_home/$item"
        fi
    done
    export DELFIN_SETTINGS_DEFAULTS="$shared/.delfin_settings.json"
    export HOME="$person_home" DELFIN_PERSON="$person" DELFIN_RUNTIME_DIR="$person_home/.run"
    # The dashboard's agent works in its launch folder (DELFIN_LAUNCH_CWD).
    if [ -z "$directory" ]; then cd -- "$person_home/agent_workspace"; fi
    printf 'DELFIN: Person %s (%s)\n' "$person" "$person_home"
fi
if [ "$role" = terminal ]; then
    # The interactive login shell already loaded startup files and activated the environment.
    # Do not load .bashrc a second time after inheriting its decorated prompt.
    export PS1
    exec bash --norc -i
fi
if ! serves_dashboard "$python"; then
    printf 'DELFIN: This environment has DELFIN but cannot serve the dashboard (voila is not installed): %s\n' "$python" >&2
    printf '%s\n' 'Install DELFIN with its dependencies into it, or set Advanced settings > DELFIN location to an environment that has them.' >&2
    exit 1
fi
# Capability check runs argparse help only, before any dashboard/tmux startup.
if ! dashboard_help=$("$python" -c "$no_cwd"'; from delfin.cli_voila import main; main(["--help"])' 2>&1); then
    printf '%s\n' 'DELFIN: Cannot read dashboard CLI capabilities. Check the server DELFIN installation.' >&2
    exit 1
fi
case "$dashboard_help" in
    *--strict-port*) ;;
    *) printf '%s\n' 'DELFIN: Server DELFIN is too old for this launcher (--strict-port is missing). Update the server installation to the Windows-launcher release; updating the Windows app alone is not enough.' >&2; exit 1 ;;
esac
# A script file puts its own folder first on sys.path, not the current one.
exec "$python" "$launcher_script" "$preferred" "$PWD" "$python" "$tunnel_path" "$keep" "$working"
