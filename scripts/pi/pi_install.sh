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
# Everything except QM binaries (not pip anyway) and ai-retro, whose
# aizynthfinder pins rdkit<2023 and conflicts with DELFIN's rdkit>=2023.3.
EXTRAS="${DELFIN_PI_EXTRAS:-agent docs dev analysis mlp-full ai-foundation ai-screening ai-complex csp}"
PYMOL_GIT="git+https://github.com/schrodinger/pymol-open-source.git"

cd "$ROOT"
[[ -x "$VENV_PY" ]] || { echo "error: $VENV_PY missing - run pi_bootstrap.sh first" >&2; exit 1; }

install_deps() {
    "$VENV_PY" - "$EXTRAS" "$SKIP_DEPS_RE" > "$ROOT/requirements.pi.txt" <<'EOF'
import re, sys, tomllib
extras, skip = sys.argv[1].split(), re.compile(sys.argv[2])
proj = tomllib.load(open("src/pyproject.toml", "rb"))["project"]
deps = list(proj["dependencies"])
for e in extras:
    deps += proj.get("optional-dependencies", {}).get(e, [])
for d in dict.fromkeys(deps):
    if not skip.match(d):
        print(d)
EOF
    # CPU torch: PyPI's aarch64 torch drags in multi-GB CUDA (sbsa) wheels.
    "$UV" pip install --python "$VENV_PY" --torch-backend cpu -r "$ROOT/requirements.pi.txt"
    "$UV" pip install --python "$VENV_PY" --no-deps -e "$ROOT/src"
    if ! "$VENV_PY" -c "import pymol" 2>/dev/null; then
        echo "==> building pymol from source (optional, only for report images)"
        # --no-deps: pymol's master declares numpy>=2 and would upgrade past
        # DELFIN's numpy<2 pin; a numpy-2 build runs fine on numpy 1.26.
        "$UV" pip install --python "$VENV_PY" --no-deps "$PYMOL_GIT" \
            || echo "warning: pymol build failed - DOCX/ESP report images unavailable" >&2
    fi
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
    [[ -d src ]] && mv src src.prev
    mv src.new src

    if [[ "$old_hash" != "$new_hash" ]] || ! "$VENV_PY" -c "import delfin" 2>/dev/null; then
        echo "==> pyproject changed (or first deploy): installing dependencies"
        install_deps
    else
        echo "==> pyproject unchanged: code swapped, no reinstall needed"
    fi
    smoke
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
    if ! cmp -s src/pyproject.toml src.prev/pyproject.toml; then install_deps; fi
    smoke
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
