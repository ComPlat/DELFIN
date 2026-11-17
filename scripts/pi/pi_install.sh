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
# with DELFIN's rdkit>=2023.3), pymol (built below).
#
# The orca extra (orca-pi 2.0.0, pure Python) is installed. orca-pi needs
# numpy>=2.2.6 and rdkit<2026; every mendeleev DELFIN admits (<1.0) declares
# numpy<2 on Python < 3.12, so the Pi's Python 3.11 set had "No solution
# found" while orca-pi was in it. The same mendeleev 0.20.1 declares numpy 2
# on 3.12 and the reference environment runs it on numpy 2.2.2, so the
# Pi's overrides (pi-overrides.txt) replace the numpy requirements with
# numpy>=2.2.6,<2.3 -- inside the core's own range. With it the set
# resolves for aarch64/cp311 (uv pip compile, 2026-11-03): numpy 2.2.6,
# rdkit 2025.9.6 (was 2026.9.1), mendeleev 0.20.1, orca-pi 2.0.0. Element
# data from mendeleev 0.20.1 on numpy 2.2.6 / Python 3.11 was read back
# unchanged. ORCA itself has no ARM build: OPI builds inputs and reads
# outputs here (Calculator(version_check=False) needs no ORCA binary).
SKIP_EXTRAS="${DELFIN_PI_SKIP_EXTRAS:-ai-retro pymol}"
PI_OVERRIDES="$ROOT/src/scripts/pi/pi-overrides.txt"
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

# requirements.pi.txt: the core and every extra not in $1.
write_requirements() {
    "$VENV_PY" - "$1" "$SKIP_DEPS_RE" > "$ROOT/requirements.pi.txt" <<'EOF'
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
}

uv_install_requirements() {
    # CPU torch: PyPI's aarch64 torch drags in multi-GB CUDA (sbsa) wheels.
    # shellcheck disable=SC2086
    "$UV" pip install --python "$VENV_PY" --torch-backend cpu "$@" \
        -r "$ROOT/requirements.pi.txt" $EXTRA_PIP
}

install_deps() {
    # The numpy override exists for orca-pi. A set with the orca extra that
    # does not install (a tree without pi-overrides.txt after a rollback, an
    # index without a wheel) is retried without it: the deploy goes on with
    # OPI reported missing, and the next deploy tries again (no stamp).
    local stamp=1
    if extra_skipped orca; then
        write_requirements "$SKIP_EXTRAS"
        uv_install_requirements
    elif [[ ! -f "$PI_OVERRIDES" ]]; then
        echo "warning: $PI_OVERRIDES missing: the orca extra is left out of this install" >&2
        write_requirements "$SKIP_EXTRAS orca"
        uv_install_requirements
        stamp=0
    else
        write_requirements "$SKIP_EXTRAS"
        if ! uv_install_requirements --overrides "$PI_OVERRIDES"; then
            echo "warning: the set with the orca extra did not install (reason above); installing without it, OPI stays missing" >&2
            write_requirements "$SKIP_EXTRAS orca"
            uv_install_requirements
            stamp=0
        fi
    fi
    install_delfin_editable
    if ! "$VENV_PY" -c "import pymol" 2>/dev/null; then
        echo "==> building pymol from source (optional, only for report images)"
        # --no-deps: pymol's master declares numpy>=2 and would move numpy
        # past what the extras resolved (1.26 when the orca extra is
        # skipped); a numpy-2 build runs on 1.26 and on 2.2.
        "$UV" pip install --python "$VENV_PY" --no-deps "$PYMOL_GIT" \
            || echo "warning: pymol build failed - DOCX/ESP report images unavailable" >&2
    fi
    if [[ "$stamp" == 1 ]]; then deps_stamp > "$DEPS_STAMP"; else rm -f "$DEPS_STAMP"; fi
}

extra_skipped() {
    [[ " $SKIP_EXTRAS " == *" $1 "* ]]
}

# What the dependency set is computed from: pyproject.toml, the Pi's
# overrides and the skipped extras. A deploy that changes only this
# script's choices (the orca extra coming in) reinstalls too; comparing
# pyproject.toml alone never noticed.
DEPS_STAMP="$ROOT/.deps-inputs.sha256"
deps_stamp() {
    { cat "$ROOT/src/pyproject.toml" "$PI_OVERRIDES" 2>/dev/null; printf 'skip=%s\n' "$SKIP_EXTRAS"; } \
        | sha256sum | cut -d' ' -f1
}

deps_current() {
    [[ -f "$DEPS_STAMP" && "$(cat "$DEPS_STAMP")" == "$(deps_stamp)" ]] \
        && "$VENV_PY" -c "import delfin" 2>/dev/null
}

# ---------------------------------------------------------------------------
# Optional groups: what the Pi's hardware (aarch64, CPU only) can run beyond
# the extras. Each is switched by DELFIN_PI_<GROUP> (1 = install when
# missing, the default; 0 = leave it; nothing is ever removed), and deploy_pi.sh
# forwards DELFIN_PI_* from the machine it runs on. A group that cannot be
# installed is reported with the reason and the deploy goes on. The Python
# groups install against a freeze of the venv as constraints, so they never
# move a version the extras resolved; one that would is refused and reported.
# Versions resolved for aarch64/cp311 on 2026-11-03 (scripts/pi/README.md).
group_specs() {
    case "$1" in
        pyscf)  echo "pyscf==2.14.0" ;;
        openmm) echo "openmm==8.6.1" ;;
        torch)  echo "torch>=2.5,<3 mace-torch==0.3.17" ;;
        *) return 1 ;;
    esac
}
# conda-forge linux-aarch64 builds through DELFIN's own QM installer
# (delfin/qm_tools/install_qm_tools.sh: one environment per tool, micromamba
# fetched as one pinned release held to its SHA-256, never curl | sh) into the
# user's qm_tools, where DELFIN's resolver looks. xtb pulls tblite and dftd4
# into its environment; their binaries are linked beside it.
PI_XTB_SPEC="xtb==6.7.1"
PI_CREST_SPEC="crest==3.0.2"
PI_QM_ROOT="$HOME/.delfin/qm_tools"
PI_GROUPS="pyscf openmm torch xtb"

group_enabled() {
    local var="DELFIN_PI_${1^^}"
    [[ "${!var:-1}" == "1" ]]
}

# The first line of a failed install that says why.
first_reason() {
    grep -m1 -iE 'error|because|no solution|not found|failed' <<<"$1" \
        || tail -n1 <<<"$1"
}

install_pip_group() {
    local name="$1" specs frozen out
    specs="$(group_specs "$name")"
    frozen="$(mktemp "$ROOT/.pi-frozen.XXXXXX")"
    "$UV" pip freeze --python "$VENV_PY" 2>/dev/null \
        | grep -E '^[A-Za-z0-9._-]+==' > "$frozen" || true
    # shellcheck disable=SC2086
    if out="$("$UV" pip install --python "$VENV_PY" --torch-backend cpu -c "$frozen" $specs 2>&1)"; then
        echo "==> optional $name: installed ($specs)"
    else
        echo "warning: optional $name not installed: $(first_reason "$out")" >&2
    fi
    rm -f "$frozen"
}

install_xtb_group() {
    local installer="$ROOT/src/delfin/qm_tools/install_qm_tools.sh" out b
    echo "==> optional xtb: installing $PI_XTB_SPEC $PI_CREST_SPEC from conda-forge (about 55 MB)"
    if ! out="$(DELFIN_QM_TOOLS_ROOT="$PI_QM_ROOT" XTB_SPEC="$PI_XTB_SPEC" CREST_SPEC="$PI_CREST_SPEC" \
            bash "$installer" xtb crest </dev/null 2>&1)"; then
        echo "warning: optional xtb not (fully) installed: $(first_reason "$out")" >&2
    fi
    for b in tblite dftd4; do
        if [[ -x "$PI_QM_ROOT/.mamba_env/xtb/bin/$b" ]]; then
            ln -sfn "$PI_QM_ROOT/.mamba_env/xtb/bin/$b" "$PI_QM_ROOT/bin/$b"
        fi
    done
    [[ -x "$PI_QM_ROOT/bin/xtb" ]] && echo "==> optional xtb: $PI_XTB_SPEC $PI_CREST_SPEC (+ tblite, dftd4)"
    return 0
}

# Status of every group from DELFIN's own probe (nothing is run).
group_status() {
    "$VENV_PY" -m delfin.agent.optional_tools list 2>/dev/null \
        | awk -v g="$1" '$1 == "smoke" && $2 == g {print $3}'
}

install_optional_groups() {
    local g
    for g in $PI_GROUPS; do
        if ! group_enabled "$g"; then
            echo "==> optional $g: off (DELFIN_PI_${g^^}=0)"
            continue
        fi
        [[ "$(group_status "$g")" == "ok" ]] && continue
        case "$g" in
            xtb) install_xtb_group ;;
            *)   install_pip_group "$g" ;;
        esac
    done
}

# One tiny computation per group (delfin/agent/optional_tools.py), printed
# as "smoke <group> ok|partial|missing|failed|off <detail>". Informational:
# a group that fails here does not fail the deploy.
smoke_optional() {
    local off="" g
    for g in $PI_GROUPS; do
        group_enabled "$g" || off+=" $g"
    done
    extra_skipped orca && off+=" opi"
    echo "==> optional tools: smoke test"
    "$VENV_PY" -m delfin.agent.optional_tools smoke --off "$off" \
        || echo "warning: optional-tools smoke test did not run" >&2
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
    # written into the account's DELFIN settings on every deploy, like
    # pi-tabs.conf. In the launcher's people mode that file is read under
    # each person's own (DELFIN_SETTINGS_DEFAULTS), so these are everyone's
    # defaults and a person's own choice survives the deploy. The variable
    # is cleared for the apply: it writes the account's file, whole, and
    # never a person's. Sessions start on the KIT
    # toolbox (DELFIN's own, caged tool loop): the CLI default would be
    # refused by the policy. An empty value (effort) means the model
    # profile's default. permission_profile all_free is the tab's "Bypass".
    local agent_conf="$HOME/.config/delfin/pi-agent.conf"
    if [[ ! -f "$agent_conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Agent tab defaults on this box (key=value). Applied by pi_install.sh on every deploy." \
            "# mode: solo = Code; permission_profile: ask_all | repo_free | all_free (= Bypass)." \
            "# effort empty = the model profile's default." \
            "mode=solo" "provider=kit" "model=kit.glm-5.3-flash" \
            "permission_profile=all_free" "effort=" > "$agent_conf"
    fi
    env -u DELFIN_SETTINGS_DEFAULTS "$VENV_PY" - "$agent_conf" <<'EOF' || echo "warning: could not apply $agent_conf" >&2
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
    # Keep dashboard is this box's default for the Windows launcher: a
    # dashboard (and the agent working in it) survives a disconnect, whatever
    # the app's checkbox says; Ctrl+C in the SSH window still stops it.
    # DELFIN_PI_KEEP=0 leaves the choice to the app.
    if [[ "${DELFIN_PI_KEEP:-1}" == "1" && ! -f "$HOME/.config/delfin/launcher-keep.conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Windows launcher: every connection keeps its dashboard (remote_bootstrap.sh)." \
            > "$HOME/.config/delfin/launcher-keep.conf"
        echo "==> launcher keep session on by default"
    fi
    install -D -m 755 "$ROOT/src/scripts/pi/delfin_people.sh" "$HOME/.local/bin/delfin-people"
    # The dashboard tabs this box shows (tab ids, ~/.config/delfin/pi-tabs.conf),
    # applied on every deploy against the tabs the code has now, so a tab a
    # later DELFIN adds stays hidden too. Settings cannot be hidden (fixed).
    # The account's file again: everyone's default in the launcher's people
    # mode, a person's own tab choice stays theirs.
    local tabs_conf="$HOME/.config/delfin/pi-tabs.conf"
    if [[ ! -f "$tabs_conf" ]]; then
        mkdir -p "$HOME/.config/delfin"
        printf '%s\n' "# Dashboard tabs shown on this box (ids; Settings is always shown)." \
            "# Applied by pi_install.sh on every deploy. All ids: grep \"'id':\" delfin/dashboard/__init__.py" \
            "agent home archive office" > "$tabs_conf"
    fi
    env -u DELFIN_SETTINGS_DEFAULTS "$VENV_PY" - "$tabs_conf" "$ROOT/src/delfin/dashboard/__init__.py" <<'EOF' || echo "warning: could not apply $tabs_conf" >&2
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
    # Only a running service picks up the new code: `restart` STARTS a
    # stopped unit, and an admin who stopped one (e.g. until access is
    # secured) must not have it switched back on by a deploy. The agent
    # service goes before the messenger bridge that talks to it; running
    # agent sessions are interrupted (a kept session resumes). A unit the
    # sudo rule does not cover yet (a box bootstrapped before the rule named
    # it) is reported with the command an admin runs, never left silent.
    local unit stale=""
    for unit in delfin-agent-service delfin-messenger-bridge delfin-messenger-gateway delfin-dashboard; do
        if ! systemctl cat "$unit" >/dev/null 2>&1; then
            continue
        elif ! systemctl is-active --quiet "$unit" 2>/dev/null; then
            echo "==> $unit not running: left stopped"
        elif sudo -n systemctl restart "$unit" 2>/dev/null; then
            echo "==> $unit restarted"
        else
            stale+=" $unit"
        fi
    done
    if [[ -n "$stale" ]]; then
        echo "warning: still on the previous code until an admin runs: sudo systemctl restart$stale" >&2
        echo "         (once, as the admin, so deploys restart them: re-run pi_bootstrap.sh, AGENT_SUDO=service)" >&2
    fi
}

# Lab control (docs/design/lab-control.md): the e-stop owner and the
# dead-man as systemd USER units of this account. Idempotent: a unit file is
# replaced only when its text changed; enabled once; a running unit is
# restarted so it runs the deployed code. No root: lingering (so the units
# run without a login) is the admin's step, named when it is missing. With
# no user manager the fallback is `delfin-agent lab daemon --detach`.
install_lab_units() {
    local unit_dir="$HOME/.config/systemd/user" src="$ROOT/src/scripts/pi/lab" unit changed=0
    [[ -d "$src" ]] || return 0
    mkdir -p "$unit_dir" "$HOME/.delfin/lab"
    chmod 700 "$HOME/.delfin/lab"
    for unit in delfin-lab-estop.service delfin-lab-watchdog.service; do
        if ! cmp -s "$src/$unit" "$unit_dir/$unit"; then
            install -m 0644 "$src/$unit" "$unit_dir/$unit.new"
            mv "$unit_dir/$unit.new" "$unit_dir/$unit"
            changed=1
        fi
    done
    if ! systemctl --user show-environment >/dev/null 2>&1; then
        echo "warning: no systemd user manager for $(id -un): lab units installed, not started" >&2
        echo "         once, as an admin: sudo loginctl enable-linger $(id -un)" >&2
        echo "         meanwhile: $VENV_PY -m delfin.agent.cli lab daemon --role both --detach" >&2
        return 0
    fi
    if [[ "$changed" == "1" ]]; then
        systemctl --user daemon-reload
    fi
    for unit in delfin-lab-estop delfin-lab-watchdog; do
        if ! systemctl --user is-enabled --quiet "$unit" 2>/dev/null; then
            systemctl --user enable --now "$unit" && echo "==> $unit enabled"
        elif systemctl --user is-active --quiet "$unit"; then
            systemctl --user restart "$unit" && echo "==> $unit restarted"
        else
            echo "==> $unit not running: left stopped"
        fi
    done
}

smoke_lab() {
    # The simulator lab, one short plan end to end, in a scratch directory:
    # set, interlock, readback, e-stop, clear. Prints "smoke lab ok".
    "$VENV_PY" -m delfin.lab.smoke || echo "warning: lab smoke check failed" >&2
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

    rm -rf src.prev
    if [[ -d src ]]; then
        carry_tui_modules src src.new
        mv src src.prev
    fi
    mv src.new src

    if ! deps_current; then
        echo "==> dependency inputs changed (pyproject, Pi overrides, skipped extras) or first deploy: installing dependencies"
        install_deps
    else
        echo "==> dependency inputs unchanged: code swapped, no reinstall needed"
    fi
    ensure_compat_editable
    install_optional_groups
    drop_foreign_tests_pkg
    smoke
    smoke_optional
    set_up_extras
    restart_dashboard
    install_lab_units
    smoke_lab
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
    deps_current || install_deps
    ensure_compat_editable
    install_optional_groups
    drop_foreign_tests_pkg
    smoke
    smoke_optional
    set_up_extras
    restart_dashboard
    install_lab_units
    smoke_lab
    ;;
status)
    cat "$ROOT/src/DEPLOYED_VERSION" 2>/dev/null || echo "nothing deployed"
    ;;
*)
    echo "usage: $0 {install <version> [--test]|rollback|status}" >&2
    exit 2
    ;;
esac
