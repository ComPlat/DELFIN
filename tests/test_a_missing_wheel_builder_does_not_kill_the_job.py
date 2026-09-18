"""A runtime wheel that cannot be built costs start-up time, not the job.

The job builds itself a wheel on the compute node to start from local disk
instead of a network HOME.  It builds without isolation, because a compute
node has no dependable way to a package index -- so the build uses what the
venv has, and a venv made with ``python -m venv`` has setuptools 65 and no
``wheel``.  That setuptools cannot run bdist_wheel, and the build ends in
``error: invalid command 'bdist_wheel'``.

Under ``set -euo pipefail`` that killed the submitted calculation, and it
began to kill them everywhere on 2026-09-13, when the cache stopped being one
site's setting and became every site's.  The cache is a start-up
optimisation; losing it must cost the start-up time and nothing else.
"""

import os
import pathlib
import re
import subprocess

REPO = pathlib.Path(__file__).resolve().parents[1]
SUBMIT = REPO / "delfin" / "submit_templates" / "submit_delfin.sh"

#: The scripts that make a venv DELFIN is later run from.
VENV_MAKERS = (
    REPO / "delfin" / "installers" / "install_delfin.sh",
    REPO / "examples" / "example_Job_Submission_Scripts" / "setup_delfin.sh",
    REPO / "examples" / "example_Job_Submission_Scripts" / "setup_delfin_bwunicluster3.sh",
)

FUNCTIONS = ("venv_can_build_a_wheel", "build_runtime_wheel_into",
             "ensure_runtime_wheel", "install_cached_runtime_wheel")

#: Stands in for the venv's interpreter: it answers the wheel-builder probe
#: from a file the test controls, and writes a wheel when pip is asked for one.
STUB_PYTHON = """#!/usr/bin/env bash
if [ "$1" = "-" ]; then
    cat > /dev/null
    [ -f "$STUB_STATE/can_build" ] && exit 0 || exit 1
fi
if [ "$1" = "-m" ] && [ "$2" = "pip" ]; then
    case "$3" in
        install)
            for arg in "$@"; do
                case "$arg" in *.whl) exit 0 ;; esac   # a built wheel installs
            done
            if [ "${STUB_PIP_INSTALL_WORKS:-0}" = "1" ]; then
                touch "$STUB_STATE/can_build"; exit 0   # the build tools arrive
            fi
            echo "no route to the package index" >&2; exit 1 ;;
        wheel)
            for arg in "$@"; do
                [ "$prev" = "--wheel-dir" ] && dir="$arg"
                prev="$arg"
            done
            touch "$dir/delfin_complat-1.3.2-py3-none-any.whl"; exit 0 ;;
    esac
fi
exit 0
"""


def _functions_under_test() -> str:
    text = SUBMIT.read_text(encoding="utf-8")
    parts = []
    for name in FUNCTIONS:
        found = re.search(rf"^{name}\(\) \{{.*?^\}}", text, re.M | re.S)
        assert found, f"{name} is no longer a function of {SUBMIT.name}"
        parts.append(found.group(0))
    return "\n".join(parts)


def _run_a_job_that_wants_the_cache(tmp_path, *, pip_install_works: bool,
                                    builder_present: bool) -> subprocess.CompletedProcess:
    state = tmp_path / "state"
    state.mkdir()
    if builder_present:
        (state / "can_build").touch()

    venv = tmp_path / "venv"
    (venv / "bin").mkdir(parents=True)
    python = venv / "bin" / "python"
    python.write_text(STUB_PYTHON, encoding="utf-8")
    python.chmod(0o755)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "python").symlink_to(python)

    harness = tmp_path / "harness.sh"
    harness.write_text(
        "set -euo pipefail\n"
        f"{_functions_under_test()}\n"
        # the parts of the job around them, as far as they are asked about here
        "detect_runtime_key() { RUNTIME_KEY='test-key'; }\n"
        "build_runtime_context() { :; }\n"
        f"RUNTIME_CACHE_DIR='{tmp_path / 'cache'}'\n"
        f"STAGE_BASE='{tmp_path / 'stage'}'\n"
        f"VENV_LOCAL='{venv}'\n"
        f"DELFIN_DIR='{tmp_path / 'home_checkout'}'\n"
        "RUNTIME_WHEEL=''\nRUNTIME_KEY=''\nRUNTIME_BUILD_MODE='live-repo'\nSLURM_JOB_ID='1'\n"
        f"mkdir -p '{tmp_path / 'stage'}'\n"
        "install_cached_runtime_wheel\n"
        "echo \"RUNTIME_KEY_AFTERWARDS=${DELFIN_RUNTIME_KEY:-unset}\"\n"
        "echo 'the job goes on'\n", encoding="utf-8")

    return subprocess.run(
        ["bash", str(harness)], capture_output=True, text=True, timeout=120,
        env={"PATH": os.pathsep.join([str(tmp_path / "bin"), "/usr/bin", "/bin"]),
             "HOME": str(tmp_path), "STUB_STATE": str(state),
             "STUB_PIP_INSTALL_WORKS": "1" if pip_install_works else "0"})


def test_a_venv_that_cannot_build_a_wheel_gets_the_build_tools(tmp_path):
    done = _run_a_job_that_wants_the_cache(tmp_path, pip_install_works=True, builder_present=False)
    assert done.returncode == 0, done.stderr
    assert "cannot build a wheel yet" in done.stdout          # it noticed
    assert "RUNTIME_KEY_AFTERWARDS=test-key" in done.stdout   # and then built it


def test_a_wheel_that_cannot_be_built_leaves_the_job_running(tmp_path):
    # no builder, and no way to install one: this is the compute node without
    # a route to an index, and the job that died there
    done = _run_a_job_that_wants_the_cache(tmp_path, pip_install_works=False, builder_present=False)
    assert done.returncode == 0, done.stderr
    assert "the job goes on" in done.stdout
    assert "no runtime wheel could be built" in done.stdout
    assert "RUNTIME_KEY_AFTERWARDS=editable-home-fallback" in done.stdout


def test_a_venv_that_can_build_one_still_gets_its_cache(tmp_path):
    done = _run_a_job_that_wants_the_cache(tmp_path, pip_install_works=False, builder_present=True)
    assert done.returncode == 0, done.stderr
    assert "RUNTIME_KEY_AFTERWARDS=test-key" in done.stdout
    assert "cannot build a wheel yet" not in done.stdout      # nothing to repair
    assert (tmp_path / "cache" / "test-key").glob("delfin_complat-*.whl")


def test_the_probe_says_yes_for_an_interpreter_that_can_build_one():
    # the test runner's own interpreter: whatever it has, the probe and pip
    # have to agree about it
    import sys
    from importlib.util import find_spec
    probe = re.search(r"^venv_can_build_a_wheel\(\) \{.*?^\}", SUBMIT.read_text(encoding="utf-8"),
                      re.M | re.S).group(0)
    done = subprocess.run(["bash", "-c", f"{probe}\nvenv_can_build_a_wheel"],
                          capture_output=True, text=True, timeout=60,
                          env={**os.environ, "PATH": os.pathsep.join(
                              [str(pathlib.Path(sys.executable).parent), os.environ.get("PATH", "")])})
    can_build = find_spec("wheel") is not None or _setuptools_builds_wheels()
    assert (done.returncode == 0) is can_build, done.stderr


def _setuptools_builds_wheels() -> bool:
    try:
        from importlib.metadata import version
        return tuple(int(p) for p in version("setuptools").split(".")[:2]) >= (70, 1)
    except Exception:
        return False


def test_every_script_that_makes_a_venv_puts_the_build_tools_in_it():
    # a venv without them is one the job cannot build its wheel in, and the
    # scripts here are where the sites' venvs come from
    for script in VENV_MAKERS:
        text = script.read_text(encoding="utf-8")
        assert re.search(r"pip install[^\n]*\bwheel\b", text), \
            f"{script.relative_to(REPO)} makes a venv without a wheel builder"


def _runtime_key_for(tree: pathlib.Path) -> str:
    text = SUBMIT.read_text(encoding="utf-8")
    parts = [re.search(rf"^{name}\(\) \{{.*?^\}}", text, re.M | re.S).group(0)
             for name in ("compute_runtime_tree_hash", "compute_runtime_dirty_hash", "detect_runtime_key")]
    script = ("\n".join(parts) + f"\nDELFIN_DIR='{tree}'\nRUNTIME_KEY=''\nRUNTIME_BUILD_MODE=''\n"
              "detect_runtime_key\necho \"$RUNTIME_KEY\"\n")
    done = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    return done.stdout.strip()


def _a_checkout_without_git(tmp_path: pathlib.Path) -> pathlib.Path:
    tree = tmp_path / "delfin_checkout"
    (tree / "delfin").mkdir(parents=True)
    (tree / "delfin" / "orca.py").write_text("x = 1\n", encoding="utf-8")
    (tree / "pyproject.toml").write_text("[project]\nname = 'delfin-complat'\n", encoding="utf-8")
    return tree


def test_a_checkout_without_git_still_gets_a_key_of_its_own(tmp_path):
    # Without git -- not installed, or not on the compute node -- the key used
    # to be the constant "tree-default": the wheel built first stayed the
    # cached one, and every later job ran the DELFIN of that day however often
    # the checkout was updated.
    tree = _a_checkout_without_git(tmp_path)
    assert _runtime_key_for(tree) != "tree-default"


def test_changing_the_code_changes_the_key(tmp_path):
    tree = _a_checkout_without_git(tmp_path)
    before = _runtime_key_for(tree)
    assert _runtime_key_for(tree) == before                   # nothing moved, same wheel

    (tree / "delfin" / "orca.py").write_text("x = 2\nyy = 3\n", encoding="utf-8")
    assert _runtime_key_for(tree) != before                   # edited

    tree_with_more = _a_checkout_without_git(tmp_path / "second")
    (tree_with_more / "delfin" / "recalc_control.py").write_text("y = 1\n", encoding="utf-8")
    assert _runtime_key_for(tree_with_more) != _runtime_key_for(_a_checkout_without_git(tmp_path / "third"))
