#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Receiving end of deploy_pi.sh - runs ON the Pi as the agent user.
#
#   pi_install.sh install <version> [--test]   # source tarball on stdin
#   pi_install.sh rollback                     # swap src <-> src.prev
#   pi_install.sh status                       # print deployed version
#
# Layout under ~/delfin:
#   src/        current code (editable install into .venv)
#   src.prev/   previous deploy, for rollback
#   .venv/      Python 3.11 venv created by pi_bootstrap.sh
# ============================================================================

ROOT="$HOME/delfin"
VENV_PY="$ROOT/.venv/bin/python"
UV="$HOME/.local/bin/uv"
# PyPI has no aarch64 dist of pymol -> built from GitHub source below instead.
SKIP_DEPS_RE='^pymol-open-source'
# Every extra pyproject.toml defines, minus these (so new extras arrive on
# the Pi by themselves): ai-retro (aizynthfinder pins rdkit<2023, conflicts
# with DELFIN's rdkit>=2023.3), pymol (built below), orca. orca-pi 2.0.0
# needs numpy>=2.2.6 and rdkit<2026; on the Pi's Python 3.11 mendeleev keeps
# numpy<2 and rdkit is 2026.x, so the whole set had no solution and every
# deploy that changed pyproject.toml stopped there. OPI also wants an ORCA
# binary as soon as a Calculator is made (checked 2026-11-03), and the Pi
# is for the agent, not for ORCA runs.
SKIP_EXTRAS="${DELFIN_PI_SKIP_EXTRAS:-ai-retro orca pymol}"
# Not in any extra, but a doctor check wants an OCR engine (tesseract via apt).
EXTRA_PIP="pytesseract"
DELFIN_GIT_URL="${DELFIN_GIT_URL:-https://github.com/ComPlat/DELFIN.git}"
# Pinned to the commit verified on a Pi 5 (render test), not the moving master.
PYMOL_GIT="git+https://github.com/schrodinger/pymol-open-source.git@edcda80f3a73ebf0b5922991f85910d17ffbcce5"

cd "$ROOT"
[[ -x "$VENV_PY" ]] || { echo "error: $VENV_PY missing - run pi_bootstrap.sh first" >&2; exit 1; }

# Editable as a plain path on sys.path (compat), not setuptools' finder: the
# finder sits after the path search, so from ~ the folder ~/delfin (this
# install's root, no __init__.py) won as an empty namespace package and
# `python -c "import delfin..."` found nothing. A path entry keeps the real
# package ahead of any same-named plain folder.
install_delfin_editable() {
    "$UV" pip install --python "$VENV_PY" --no-deps --config-setting editable_mode=compat -e "$ROOT/src"
}

# A venv installed before with the finder gets the path form once.
ensure_compat_editable() {
    local site
    site="$("$VENV_PY" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
    if compgen -G "$site/__editable___delfin*finder.py" >/dev/null; then
        echo "==> DELFIN editable install: finder -> path (compat)"
        install_delfin_editable
    fi
}

install_deps() {
    "$VENV_PY" - "$SKIP_EXTRAS" "$SKIP_DEPS_RE" > "$ROOT/requirements.pi.txt" <<'EOF'
import re, sys, tomllib
skip_extras, skip = set(sys.argv[1].split()), re.compile(sys.argv[2])
proj = tomllib.load(open("src/pyproject.toml", "rb"))["project"]
deps = list(proj["dependencies"])
for name, extra in proj.get("optional-dependencies", {}).items():
    if name not in skip_extras:
        deps += extra
for d in dict.fromkeys(deps):
    if not skip.match(d):
        print(d)
EOF
    # CPU torch: PyPI's aarch64 torch drags in multi-GB CUDA (sbsa) wheels.
    # shellcheck disable=SC2086
    "$UV" pip install --python "$VENV_PY" --torch-backend cpu -r "$ROOT/requirements.pi.txt" $EXTRA_PIP
    install_delfin_editable
    if ! "$VENV_PY" -c "import pymol" 2>/dev/null; then
        echo "==> building pymol from source (optional, only for report images)"
        # --no-deps: pymol's master declares numpy>=2 and would upgrade past
        # DELFIN's numpy<2 pin; a numpy-2 build runs fine on numpy 1.26.
        "$UV" pip install --python "$VENV_PY" --no-deps "$PYMOL_GIT" \
            || echo "warning: pymol build failed - DOCX/ESP report images unavailable" >&2
    fi
}

# molsimplify and mace-torch wheels ship their own top-level `tests` package.
# A regular package beats the repo's namespace `tests/` whatever sys.path
# says, so every `from tests.conftest import ...` in the suite fails to
# collect. Nothing imports it at runtime; it goes after every deploy.
drop_foreign_tests_pkg() {
    local purelib
    purelib="$("$VENV_PY" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
    if [[ -f "$purelib/tests/__init__.py" ]]; then
        rm -rf "$purelib/tests"
        echo "==> removed third-party tests package from site-packages"
    fi
}

# Digest of what tui/node_modules is built from (packages + DELFIN patches).
tui_inputs_hash() {
    [[ -f "$1/tui/package.json" ]] || return 0
    cat "$1/tui/package.json" "$1/tui/bun.lock" "$1"/tui/patches/* 2>/dev/null | sha256sum | cut -d' ' -f1
}

# A deploy replaces src/; keep the window's packages and the patched OpenTUI
# library (minutes to build) when their inputs did not change.
carry_tui_modules() {
    local from="$1" to="$2"
    [[ -d "$from/tui/node_modules" && -d "$to/tui" && ! -e "$to/tui/node_modules" ]] || return 0
    if [[ "$(tui_inputs_hash "$from")" == "$(tui_inputs_hash "$to")" ]]; then
        mv "$from/tui/node_modules" "$to/tui/node_modules"
    fi
}

# Everything beyond pip, best-effort and idempotent: the terminal UI (pinned
# Bun + bun.lock packages) with the patched OpenTUI library for xterm.js
# terminals (VS Code), Playwright's Chromium, OCR for pytesseract, the
# manual search index, and a git checkout for the agent to work in.
set_up_extras() {
    # shellcheck disable=SC1091
    ( source "$ROOT/agent.env" && source "$ROOT/.venv/bin/activate" \
        && cd "$ROOT/work" && python -P -m delfin.agent.cli tui --setup </dev/null ) \
        || echo "warning: terminal UI setup failed (retry: delfin-agent tui --setup)" >&2
    if [[ -x "$ROOT/src/tui/scripts/build-opentui-native.sh" ]]; then
        ( cd "$ROOT/src/tui" && { scripts/build-opentui-native.sh --check >/dev/null 2>&1 \
            || scripts/build-opentui-native.sh </dev/null; } ) \
            || echo "warning: patched OpenTUI build failed (prebuilt library stays)" >&2
    fi
    if "$VENV_PY" -c "import playwright" 2>/dev/null; then
        "$VENV_PY" -m playwright install chromium </dev/null \
            || echo "warning: playwright chromium install failed" >&2
        # Its system libraries are apt packages, an admin's to install: the
        # first boot does it as root and leaves the marker.
        [[ -f /var/lib/delfin-playwright-deps.done ]] \
            || echo "note: Chromium's system libraries, once as an admin: sudo $VENV_PY -m playwright install-deps chromium" >&2
    fi
    "$ROOT/.venv/bin/delfin-docs-index" >/dev/null 2>&1 \
        || echo "warning: doc index build failed (retry: delfin-docs-index)" >&2
    if ! grep -qw landlock /sys/kernel/security/lsm 2>/dev/null; then
        echo "note: this kernel has no Landlock, so the agent asks about every command the cage" >&2
        echo "      would run unasked. Once, as an admin (pi): sudo bash $ROOT/src/scripts/pi/pi_kernel_landlock.sh build," >&2
        echo "      then: sudo bash $ROOT/src/scripts/pi/pi_kernel_landlock.sh try" >&2
    fi
    # An agent.env from before: "off" ran a command the dashboard executes for
    # the agent with the account's full rights. auto confines it.
    if grep -q '^export DELFIN_AGENT_SANDBOX=off$' "$ROOT/agent.env" 2>/dev/null; then
        sed -i 's/^export DELFIN_AGENT_SANDBOX=off$/export DELFIN_AGENT_SANDBOX=auto/' "$ROOT/agent.env"
        echo "==> agent.env: DELFIN_AGENT_SANDBOX off -> auto"
    fi
    if [[ -f "$ROOT/agent.env" ]] && ! grep -q '^export DELFIN_MCP_ISOLATION=' "$ROOT/agent.env"; then
        printf '%s\n' "# DELFIN's own MCP servers run in bwrap, seeing only their folders." \
            "export DELFIN_MCP_ISOLATION=builtin" >> "$ROOT/agent.env"
        echo "==> agent.env: DELFIN_MCP_ISOLATION=builtin"
    fi
    # The agent tab's defaults on this box (~/.config/delfin/pi-agent.conf),
    # written into the DELFIN settings on every deploy, like pi-tabs.conf.
    # Through the launcher's people mode every person shares the settings
    # file, so these are everyone's defaults. Sessions start on the KIT
    # toolbox (DELFIN's own, caged tool loop): the CLI default would be
    # refused by the policy. An empty value (effort) means the model
    # profile's default. permission_profile all_free is the tab's "Bypass".
    local agent_conf="$HOME/.config/delfin/pi-agent.conf"
    if [[ ! -f "$agent_conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Agent tab defaults on this box (key=value). Applied by pi_install.sh on every deploy." \
            "# mode: solo = Code; permission_profile: ask_all | repo_free | all_free (= Bypass)." \
            "# effort empty = the model profile's default." \
            "mode=solo" "provider=kit" "model=kit.deepseek-v4-flash" \
            "permission_profile=all_free" "effort=" > "$agent_conf"
    fi
    "$VENV_PY" - "$agent_conf" <<'EOF' || echo "warning: could not apply $agent_conf" >&2
import sys
from delfin.user_settings import load_settings, save_settings
KEYS = ("mode", "provider", "model", "permission_profile", "effort")
wanted = {}
for line in open(sys.argv[1], encoding="utf-8"):
    line = line.strip()
    if not line or line.startswith("#") or "=" not in line:
        continue
    key, _, value = (part.strip() for part in line.partition("="))
    if key not in KEYS:
        print(f"warning: pi-agent.conf: unknown key {key!r}", file=sys.stderr)
        continue
    wanted[key] = value or None
settings = load_settings()
agent = settings.setdefault("agent", {})
changed = {k: v for k, v in wanted.items() if agent.get(k) != v}
if changed:
    agent.update(changed)
    save_settings(settings)
print("==> agent defaults: " + " ".join(f"{k}={wanted[k] or '(profile default)'}" for k in KEYS if k in wanted))
EOF
    [[ -d "$ROOT/work/DELFIN/.git" ]] || git clone -q "$DELFIN_GIT_URL" "$ROOT/work/DELFIN" \
        || echo "warning: clone of $DELFIN_GIT_URL failed" >&2
    # Several people work in this account through the Windows launcher: each
    # gets ~/people/<windows user> as home (own sessions, dashboard, tmux,
    # agent service). DELFIN_PI_PEOPLE=0 leaves it off. An overview, not a
    # boundary: they are one Unix user.
    if [[ "${DELFIN_PI_PEOPLE:-1}" == "1" && ! -f "$HOME/.config/delfin/launcher-people.conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Windows launcher: one home per person under ~/people (remote_bootstrap.sh)." \
            > "$HOME/.config/delfin/launcher-people.conf"
        echo "==> launcher people mode on (~/people/<name>; overview: delfin-people)"
    fi
    install -D -m 755 "$ROOT/src/scripts/pi/delfin_people.sh" "$HOME/.local/bin/delfin-people"
    # The dashboard tabs this box shows (tab ids, ~/.config/delfin/pi-tabs.conf),
    # applied on every deploy against the tabs the code has now, so a tab a
    # later DELFIN adds stays hidden too. Settings cannot be hidden (fixed).
    # Through the launcher's people mode every person shares this setting.
    local tabs_conf="$HOME/.config/delfin/pi-tabs.conf"
    if [[ ! -f "$tabs_conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Dashboard tabs shown on this box (ids; Settings is always shown)." \
            "# Applied by pi_install.sh on every deploy. All ids: grep \"'id':\" delfin/dashboard/__init__.py" \
            "agent home archive office" > "$tabs_conf"
    fi
    "$VENV_PY" - "$tabs_conf" "$ROOT/src/delfin/dashboard/__init__.py" <<'EOF' || echo "warning: could not apply $tabs_conf" >&2
import re, sys
from delfin.user_settings import load_settings, save_settings
shown = [w for line in open(sys.argv[1], encoding="utf-8")
         if not line.lstrip().startswith("#") for w in line.split()]
ids = list(dict.fromkeys(re.findall(r"'id':\s*'([a-z_]+)'", open(sys.argv[2], encoding="utf-8").read())))
unknown = [t for t in shown if t not in ids]
if unknown:
    print(f"warning: pi-tabs.conf names unknown tabs: {' '.join(unknown)}", file=sys.stderr)
settings = load_settings()
tabs = settings.setdefault("ui", {}).setdefault("tabs", {})
order = [t for t in shown if t in ids] + [t for t in ids if t not in shown and t != "settings"]
hidden = [t for t in ids if t not in shown and t != "settings"]
if tabs.get("order") != order or tabs.get("hidden") != hidden:
    tabs["order"], tabs["hidden"] = order, hidden
    save_settings(settings)
    print("==> dashboard tabs: " + " ".join(t for t in shown if t in ids) + " (+ settings)")
EOF
    # tmux waits 10 ms (3.5 default) for the rest of an escape sequence; over
    # SSH and a VPN the terminal's answer to tmux's capability query arrives
    # later, and its text ("^[[?61;1;6;...c") landed in the pane as typing.
    if ! grep -qs '^set -sg escape-time' "$HOME/.tmux.conf"; then
        printf '%s\n' "# Over SSH/VPN: wait long enough for the terminal's answers (pi_install.sh)." \
            "set -sg escape-time 50" >> "$HOME/.tmux.conf"
    fi
}

restart_dashboard() {
    # Only a running dashboard picks up the new code: `restart` STARTS a
    # stopped unit, and an admin who stopped it (e.g. until access is
    # secured) must not have it switched back on by a deploy.
    if systemctl is-active --quiet delfin-dashboard 2>/dev/null; then
        sudo -n systemctl restart delfin-dashboard 2>/dev/null \
            && echo "==> dashboard restarted" || true
    else
        echo "==> dashboard not running: left stopped"
    fi
}

smoke() {
    "$VENV_PY" -c "import delfin, delfin.agent.engine, rdkit; print('import ok')"
    cat "$ROOT/src/DEPLOYED_VERSION"
}

case "${1:-}" in
install)
    version="${2:?version missing}"
    run_tests="${3:-}"
    rm -rf src.new && mkdir src.new
    tar -xzf - -C src.new
    printf '%s deployed %s\n' "$version" "$(date -Is)" > src.new/DEPLOYED_VERSION

    old_hash="$(sha256sum src/pyproject.toml 2>/dev/null | cut -d' ' -f1 || true)"
    new_hash="$(sha256sum src.new/pyproject.toml | cut -d' ' -f1)"
    rm -rf src.prev
    if [[ -d src ]]; then
        carry_tui_modules src src.new
        mv src src.prev
    fi
    mv src.new src

    if [[ "$old_hash" != "$new_hash" ]] || ! "$VENV_PY" -c "import delfin" 2>/dev/null; then
        echo "==> pyproject changed (or first deploy): installing dependencies"
        install_deps
    else
        echo "==> pyproject unchanged: code swapped, no reinstall needed"
    fi
    ensure_compat_editable
    drop_foreign_tests_pkg
    smoke
    set_up_extras
    restart_dashboard
    if [[ "$run_tests" == "--test" ]]; then
        echo "==> pytest (not slow)"
        # Activate: code under test shells out to a bare `python`.
        # shellcheck disable=SC1091
        source "$ROOT/.venv/bin/activate"
        cd src && python -m pytest -q -m "not slow" -p no:cacheprovider
    fi
    ;;
rollback)
    [[ -d src.prev ]] || { echo "error: no src.prev to roll back to" >&2; exit 1; }
    mv src src.bad && mv src.prev src && mv src.bad src.prev
    carry_tui_modules src.prev src
    if ! cmp -s src/pyproject.toml src.prev/pyproject.toml; then install_deps; fi
    ensure_compat_editable
    drop_foreign_tests_pkg
    smoke
    set_up_extras
    restart_dashboard
    ;;
status)
    cat "$ROOT/src/DEPLOYED_VERSION" 2>/dev/null || echo "nothing deployed"
    ;;
*)
    echo "usage: $0 {install <version> [--test]|rollback|status}" >&2
    exit 2
    ;;
esac
