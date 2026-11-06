"""A contributor's precheck trips the project-wide pins before CI does.

A pull request can pass its own tests and still fail CI on a pin that
lives somewhere else: a tool schema without executor dispatch, the
tool-surface registry, ruff's gated selection, the entry points, the
console-script gate, child_env, the licence guard, the TUI's generated
protocol types. ``delfin-dev precheck`` picks those pins from the branch
diff and runs them locally.

Input: a set of changed paths (as ``git diff --name-only origin/main...HEAD``
would list them). Output: the pin test files to run, plus the protocol
regeneration check. The mapping is a table in delfin/dev_precheck.py;
these tests pin the table with sample diffs, and pin the Python protocol
generator against the committed tui/src/protocol/generated.ts.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from delfin.dev_precheck import (
    PROTOCOL_PIN,
    changed_files,
    generate_protocol_text,
    load_protocol_schema,
    select_pins,
)

REPO = Path(__file__).resolve().parents[1]


def _assert_pins_exist(pins: set[str]) -> None:
    for pin in pins:
        assert (REPO / pin).is_file(), f"pin {pin} does not exist in the repo"


def test_a_tool_schema_change_picks_the_two_tool_pins():
    pins = select_pins(["delfin/agent/api_client.py"])
    assert "tests/test_a_tool_reaches_the_function_behind_it.py" in pins
    assert "tests/test_the_tool_surface_is_exercised.py" in pins
    _assert_pins_exist(pins)


def test_the_pyproject_entry_points_pick_the_entry_point_pins():
    pins = select_pins(["pyproject.toml"])
    assert "tests/test_functional_contracts.py" in pins
    assert "tests/test_an_acting_delfin_command_is_asked.py" in pins
    _assert_pins_exist(pins)


def test_a_protocol_change_picks_the_tui_protocol_pin():
    pins = select_pins(["delfin/agent/service/protocol.py"])
    assert PROTOCOL_PIN in pins


def test_an_agent_user_state_change_picks_the_child_env_pin():
    pins = select_pins(["delfin/agent/terminal_confirm.py"])
    assert "tests/test_a_child_process_writes_into_no_real_home.py" in pins
    _assert_pins_exist(pins)


def test_a_licence_table_change_picks_the_licence_guard():
    pins = select_pins(["REUSE.toml", "LICENSES/MIT.txt"])
    assert "tests/test_every_file_has_a_licence_by_reuse.py" in pins
    assert "tests/test_license_guard.py" in pins
    _assert_pins_exist(pins)


def test_a_launcher_change_picks_the_launcher_checksum_pin():
    pins = select_pins(["tools/windows-launcher/build_package.py"])
    assert "tests/test_windows_launcher.py" in pins
    _assert_pins_exist(pins)


def test_a_changed_test_file_runs_itself():
    pins = select_pins(["tests/test_something_new.py"])
    assert "tests/test_something_new.py" in pins


def test_an_unrelated_change_picks_no_pin():
    assert select_pins(["delfin/xyz_io.py"]) == set()


def test_select_pins_is_a_union_over_several_changed_paths():
    pins = select_pins(["delfin/agent/api_client.py", "pyproject.toml"])
    assert "tests/test_a_tool_reaches_the_function_behind_it.py" in pins
    assert "tests/test_functional_contracts.py" in pins
    _assert_pins_exist(pins)


def test_changed_files_reads_the_branch_diff():
    files = changed_files("HEAD")
    assert isinstance(files, list)


def test_the_python_protocol_generator_matches_the_committed_generated_ts():
    schema = load_protocol_schema()
    text = generate_protocol_text(schema)
    committed = (REPO / "tui" / "src" / "protocol" / "generated.ts").read_text(
        encoding="utf-8"
    )
    assert text == committed


def test_the_schema_loader_runs_this_interpreter():
    out = subprocess.run(
        [sys.executable, "-m", "delfin.agent.service.protocol"],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert out.returncode == 0, out.stderr
    assert "protocolVersion" in out.stdout
