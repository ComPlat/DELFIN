# Runs inside the user's interactive login shell. Arguments are data, never eval'ed.
set -e
preferred=$1
directory=$2
hint=$3
payload=$4
role=$5
tunnel_path=${6:-}
keep=${7:-0}
working=${8:-0}
expand_path() {
    case "$1" in
        '~') printf '%s' "$HOME" ;;
        '~/'*) printf '%s/%s' "$HOME" "${1:2}" ;;
        *) printf '%s' "$1" ;;
    esac
}
directory=$(expand_path "$directory")
hint=$(expand_path "$hint")
if [ -n "$directory" ]; then cd -- "$directory"; else cd -- "$HOME"; fi
case "$hint" in ''|/*) ;; *) hint="$PWD/$hint" ;; esac
usable() {
    [ -x "$1" ] && "$1" -c 'import delfin.cli_voila, delfin.agent.where, delfin.dashboard.session' >/dev/null 2>&1
}
python=''
choose() { if [ -z "$python" ] && usable "$1"; then python=$1; fi; }
scan_root() {
    local root=$1 candidate
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
    if [ ${#matches[@]} -eq 1 ]; then python=${matches[0]}; fi
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
        candidate=$(type -P delfin-voila || true)
        if [ -n "$candidate" ]; then choose "$(dirname "$candidate")/python"; choose "$(dirname "$candidate")/python3"; fi
        candidate=$(type -P python || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
        candidate=$(type -P python3 || true); if [ -n "$candidate" ]; then choose "$candidate"; fi
    fi
fi
if [ -z "$python" ]; then
    printf '%s\n' 'DELFIN: Not found in the login environment or repository. Set a working directory or DELFIN environment location.' >&2
    exit 1
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
if [ "$role" = terminal ]; then
    # The interactive login shell already loaded startup files and activated the environment.
    # Do not load .bashrc a second time after inheriting its decorated prompt.
    export PS1
    exec bash --norc -i
fi
# Capability check runs argparse help only, before any dashboard/tmux startup.
if ! dashboard_help=$("$python" -c 'from delfin.cli_voila import main; main(["--help"])' 2>&1); then
    printf '%s\n' 'DELFIN: Cannot read dashboard CLI capabilities. Check the server DELFIN installation.' >&2
    exit 1
fi
case "$dashboard_help" in
    *--strict-port*) ;;
    *) printf '%s\n' 'DELFIN: Server DELFIN is too old for this launcher (--strict-port is missing). Update the server installation to the Windows-launcher release; updating the Windows app alone is not enough.' >&2; exit 1 ;;
esac
exec "$python" -c 'import base64,sys;exec(base64.b64decode(sys.argv.pop(1)))' "$payload" "$preferred" "$PWD" "$python" "$tunnel_path" "$keep" "$working"
