import importlib.util
import json
import shlex
import subprocess
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[1] / "delfin" / "remote_archive.py"
_SPEC = importlib.util.spec_from_file_location("delfin_remote_archive", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load remote_archive module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

build_remote_absolute_path = _MODULE.build_remote_absolute_path
build_remote_fetch_command = _MODULE.build_remote_fetch_command
build_remote_list_command = _MODULE.build_remote_list_command
build_remote_text_preview_command = _MODULE.build_remote_text_preview_command
get_cached_local_path = _MODULE.get_cached_local_path
normalize_remote_relative_path = _MODULE.normalize_remote_relative_path


def test_normalize_remote_relative_path_accepts_nested_paths():
    assert normalize_remote_relative_path("job_01/results/out.log") == "job_01/results/out.log"
    assert normalize_remote_relative_path("/job_01/results") == "job_01/results"


def test_normalize_remote_relative_path_rejects_parent_escape():
    try:
        normalize_remote_relative_path("../secret")
    except ValueError as exc:
        assert "configured root" in str(exc).lower()
    else:
        raise AssertionError("Expected ValueError for parent path escape")


def test_build_remote_absolute_path_joins_under_root():
    assert (
        build_remote_absolute_path("/home/exampleuser/archive", "job_01/out.log")
        == "/home/exampleuser/archive/job_01/out.log"
    )


def test_get_cached_local_path_uses_home_cache_namespace(tmp_path):
    local_path = get_cached_local_path(
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        "job_01/out.log",
        cache_dir=tmp_path,
    )
    assert str(local_path).startswith(str(tmp_path))
    assert local_path.name == "out.log"


def test_build_remote_list_command_contains_remote_python_listing():
    command = build_remote_list_command("/remote/archive", "job_01", sort_mode="date_desc")
    assert command.startswith("python3 -c ")
    assert "date_desc" in command
    assert "/remote/archive" in command
    assert " -- " not in command


def test_build_remote_list_command_executes_with_root_directory(tmp_path):
    remote_root = tmp_path / "remote_root"
    remote_root.mkdir()
    (remote_root / "a.txt").write_text("alpha", encoding="utf-8")
    (remote_root / "subdir").mkdir()

    command = build_remote_list_command(str(remote_root), "", sort_mode="name")
    result = subprocess.run(
        shlex.split(command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["relative_path"] == ""
    assert [entry["name"] for entry in payload["entries"]] == ["subdir", "a.txt"]


def test_build_remote_text_preview_command_executes_with_file(tmp_path):
    remote_root = tmp_path / "remote_root"
    remote_root.mkdir()
    file_path = remote_root / "a.txt"
    file_path.write_text("alpha beta", encoding="utf-8")

    command = build_remote_text_preview_command(str(remote_root), "a.txt", max_bytes=32)
    result = subprocess.run(
        shlex.split(command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["text"] == "alpha beta"
    assert payload["truncated"] is False


def test_build_remote_fetch_command_uses_rsync_over_ssh(tmp_path):
    destination = tmp_path / "cache" / "out.log"
    command = build_remote_fetch_command(
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        "job_01/out.log",
        destination,
    )
    assert command[0] == "rsync"
    assert "--append-verify" in command
    assert "--timeout=86400" in command
    assert command[-1] == str(destination)
    assert "alice@cluster.example.org:/remote/archive/job_01/out.log" == command[-2]


def test_build_remote_fetch_command_keeps_spaces_unquoted(tmp_path):
    destination = tmp_path / "cache" / "Opt_TS3_product"
    command = build_remote_fetch_command(
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        "Phenonium Project/TS3/Opt_TS3_product",
        destination,
    )

    assert command[-2] == "alice@cluster.example.org:/remote/archive/Phenonium Project/TS3/Opt_TS3_product"
    assert "'" not in command[-2]
