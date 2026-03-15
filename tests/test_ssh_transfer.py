import importlib.util
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[1] / "delfin" / "ssh_transfer_jobs.py"
_SPEC = importlib.util.spec_from_file_location("delfin_ssh_transfer_jobs", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load ssh_transfer_jobs module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

build_rsync_transfer_command = _MODULE.build_rsync_transfer_command
build_rsync_download_command = _MODULE.build_rsync_download_command
build_ssh_delete_command = _MODULE.build_ssh_delete_command
build_ssh_mkdir_command = _MODULE.build_ssh_mkdir_command
create_download_job = _MODULE.create_download_job
create_transfer_job = _MODULE.create_transfer_job
is_retryable_failure = _MODULE._is_retryable_failure
list_transfer_jobs = _MODULE.list_transfer_jobs
normalize_ssh_transfer_settings = _MODULE.normalize_ssh_transfer_settings
remote_directory_spec = _MODULE.remote_directory_spec


def test_normalize_ssh_transfer_settings_accepts_valid_values():
    assert normalize_ssh_transfer_settings(
        "cluster.example.org",
        "alice",
        "~/archive target",
        "2222",
    ) == ("cluster.example.org", "alice", "~/archive target", 2222)


def test_normalize_ssh_transfer_settings_rejects_invalid_host():
    try:
        normalize_ssh_transfer_settings("bad host", "alice", "/tmp/out", 22)
    except ValueError as exc:
        assert "host" in str(exc).lower()
    else:
        raise AssertionError("Expected ValueError for invalid host")


def test_remote_directory_spec_keeps_tilde_expansion_with_spaces():
    assert (
        remote_directory_spec("alice", "cluster.example.org", "~/archive target")
        == "alice@cluster.example.org:~/'archive target/'"
    )


def test_build_ssh_mkdir_command_quotes_remote_folder():
    command = build_ssh_mkdir_command(
        "cluster.example.org",
        "alice",
        "~/archive target",
        2222,
    )
    assert command[:4] == ["ssh", "-p", "2222", "-o"]
    assert command[-2] == "alice@cluster.example.org"
    assert command[-1] == "mkdir -p -- ~/'archive target'"


def test_build_rsync_transfer_command_builds_resumable_rsync(tmp_path):
    source_file = tmp_path / "result.txt"
    source_file.write_text("ok", encoding="utf-8")
    source_dir = tmp_path / "job_01"
    source_dir.mkdir()
    (source_dir / "out.log").write_text("done", encoding="utf-8")

    command = build_rsync_transfer_command(
        [source_file, source_dir],
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
    )

    assert command[0] == "rsync"
    assert "--partial" in command
    assert "--append-verify" in command
    assert "--timeout=86400" in command
    assert "--contimeout=20" not in command
    assert str(source_file) in command
    assert str(source_dir) in command
    assert command[-1] == "alice@cluster.example.org:/remote/archive/"


def test_create_transfer_job_writes_job_and_status_files(tmp_path):
    source_dir = tmp_path / "calc_a"
    source_dir.mkdir()
    (source_dir / "result.out").write_text("done", encoding="utf-8")

    job = create_transfer_job(
        [source_dir],
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        jobs_dir=tmp_path / "jobs",
        max_retries=5,
    )

    job_path = Path(job["job_path"])
    status_path = Path(job["status_path"])
    assert job_path.exists()
    assert status_path.exists()

    jobs = list_transfer_jobs(jobs_dir=tmp_path / "jobs", limit=5)
    assert len(jobs) == 1
    assert jobs[0]["job_id"] == job["job_id"]
    assert jobs[0]["status"] == "queued"
    assert jobs[0]["sources"] == [str(source_dir)]
    assert jobs[0]["delete_local_on_success"] is False


def test_create_transfer_job_can_mark_local_source_for_cleanup(tmp_path):
    source_dir = tmp_path / "calc_b"
    source_dir.mkdir()
    job = create_transfer_job(
        [source_dir],
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        jobs_dir=tmp_path / "jobs",
        delete_local_on_success=True,
    )

    jobs = list_transfer_jobs(jobs_dir=tmp_path / "jobs", limit=5)
    assert len(jobs) == 1
    assert jobs[0]["job_id"] == job["job_id"]
    assert jobs[0]["delete_local_on_success"] is True


def test_build_rsync_download_command_preserves_remote_relative_paths(tmp_path):
    destination = tmp_path / "calc"
    command = build_rsync_download_command(
        ["job_01/out.log", "job_02"],
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        destination,
    )

    assert command[0] == "rsync"
    assert "--relative" in command
    assert "--append-verify" in command
    assert "--timeout=86400" in command
    assert "alice@cluster.example.org:/remote/archive/./job_01/out.log" == command[-3]
    assert "alice@cluster.example.org:/remote/archive/./job_02" == command[-2]
    assert command[-1] == str(destination)


def test_build_ssh_delete_command_uses_remote_python_cleanup():
    command = build_ssh_delete_command(
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        ["job_01/out.log"],
    )
    assert command[:4] == ["ssh", "-p", "22", "-o"]
    assert command[-2] == "alice@cluster.example.org"
    assert command[-1].startswith("python3 -c ")
    assert "/remote/archive" in command[-1]
    assert "job_01/out.log" in command[-1]


def test_create_download_job_writes_pull_job_and_status_files(tmp_path):
    job = create_download_job(
        ["job_01/out.log"],
        "cluster.example.org",
        "alice",
        "/remote/archive",
        22,
        tmp_path / "calc",
        jobs_dir=tmp_path / "jobs",
        max_retries=5,
        delete_remote_on_success=True,
    )

    job_path = Path(job["job_path"])
    status_path = Path(job["status_path"])
    assert job_path.exists()
    assert status_path.exists()

    jobs = list_transfer_jobs(jobs_dir=tmp_path / "jobs", limit=5)
    assert len(jobs) == 1
    assert jobs[0]["job_id"] == job["job_id"]
    assert jobs[0]["status"] == "queued"
    assert jobs[0]["direction"] == "pull"
    assert jobs[0]["delete_remote_on_success"] is True
    assert jobs[0]["local_target"] == str((tmp_path / "calc").resolve())


def test_auth_failure_is_not_treated_as_retryable():
    assert is_retryable_failure(255, "Permission denied (publickey,password).") is False


def test_network_timeout_is_treated_as_retryable():
    assert is_retryable_failure(255, "Connection timed out") is True
