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
# with DELFIN's rdkit>=2023.3), orca (QM, out of scope), pymol (built below).
SKIP_EXTRAS="${DELFIN_PI_SKIP_EXTRAS:-ai-retro orca pymol}"
# Not in any extra, but a doctor check wants an OCR engine (tesseract via apt).
EXTRA_PIP="pytesseract"
DELFIN_GIT_URL="${DELFIN_GIT_URL:-https://github.com/ComPlat/DELFIN.git}"
# Pinned to the commit verified on a Pi 5 (render test), not the moving master.
PYMOL_GIT="git+https://github.com/schrodinger/pymol-open-source.git@edcda80f3a73ebf0b5922991f85910d17ffbcce5"

cd "$ROOT"
[[ -x "$VENV_PY" ]] || { echo "error: $VENV_PY missing - run pi_bootstrap.sh first" >&2; exit 1; }

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
    "$UV" pip install --python "$VENV_PY" --no-deps -e "$ROOT/src"
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
    [[ -d "$ROOT/work/DELFIN/.git" ]] || git clone -q "$DELFIN_GIT_URL" "$ROOT/work/DELFIN" \
        || echo "warning: clone of $DELFIN_GIT_URL failed" >&2
}

restart_dashboard() {
    sudo -n systemctl restart delfin-dashboard 2>/dev/null \
        && echo "==> dashboard restarted" || true
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
