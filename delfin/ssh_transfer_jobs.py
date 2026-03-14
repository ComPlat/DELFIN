"""Detached, resumable SSH transfer jobs for the dashboard archive browser."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from delfin.user_settings import normalize_ssh_transfer_settings

JOBS_DIR_NAME = ".delfin_transfer_jobs"
DEFAULT_MAX_RETRIES = 60
DEFAULT_RSYNC_TIMEOUT_SECONDS = 24 * 60 * 60
DEFAULT_RETRY_DELAY_SECONDS = 15
DEFAULT_RETRY_BACKOFF_CAP_SECONDS = 300
RETRYABLE_EXIT_CODES = {10, 11, 12, 30, 35, 255}


def utc_now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def ensure_jobs_dir(base_dir=None):
    jobs_dir = Path(base_dir).expanduser() if base_dir else (Path.home() / JOBS_DIR_NAME)
    jobs_dir.mkdir(parents=True, exist_ok=True)
    return jobs_dir


def quote_remote_path(remote_path):
    path = str(remote_path or "")
    if path.startswith("~"):
        head, slash, rest = path.partition("/")
        if not slash:
            return head
        return f"{head}/{shlex.quote(rest)}"
    return shlex.quote(path)


def remote_directory_spec(user, host, remote_path):
    path = str(remote_path)
    if not path.endswith("/"):
        path = f"{path}/"
    return f"{user}@{host}:{quote_remote_path(path)}"


def build_ssh_command_base(host, user, port):
    return [
        "ssh",
        "-p",
        str(port),
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=20",
        "-o",
        "ServerAliveInterval=60",
        "-o",
        "ServerAliveCountMax=5",
        "-o",
        "StrictHostKeyChecking=accept-new",
        f"{user}@{host}",
    ]


def build_ssh_mkdir_command(host, user, remote_path, port):
    command = build_ssh_command_base(host, user, port)
    command.append(f"mkdir -p -- {quote_remote_path(remote_path)}")
    return command


def build_rsync_ssh_command(host, user, port):
    ssh_command = build_ssh_command_base(host, user, port)[:-1]
    return " ".join(shlex.quote(part) for part in ssh_command)


def build_rsync_transfer_command(sources, host, user, remote_path, port):
    normalized_sources = [str(Path(src)) for src in sources]
    if not normalized_sources:
        raise ValueError("No local files selected for transfer.")

    command = [
        "rsync",
        "-a",
        "--partial",
        "--append-verify",
        "--protect-args",
        "--human-readable",
        "--info=progress2,stats1",
        f"--timeout={DEFAULT_RSYNC_TIMEOUT_SECONDS}",
        "-e",
        build_rsync_ssh_command(host, user, port),
    ]
    command.extend(normalized_sources)
    command.append(remote_directory_spec(user, host, remote_path))
    return command


def _coerce_sources(sources: Iterable[str | Path]):
    resolved = []
    for source in sources:
        path = Path(source).expanduser().resolve()
        if not path.exists():
            raise ValueError(f"Transfer source not found: {path}")
        resolved.append(path)
    if not resolved:
        raise ValueError("No files selected for transfer.")
    return resolved


def _write_json_atomic(path, payload):
    path = Path(path)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp_path.replace(path)


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _tail_text(path, max_bytes=6000):
    target = Path(path)
    if not target.exists():
        return ""
    with target.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        data = handle.read().decode("utf-8", errors="ignore")
    return data.strip()


def _summarize_output_text(text):
    lines = [" ".join(line.split()) for line in str(text or "").replace("\r", "\n").splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        return ""
    interesting_tokens = (
        "permission denied",
        "host key verification failed",
        "could not resolve hostname",
        "name or service not known",
        "connection refused",
        "connection timed out",
        "broken pipe",
        "connection reset",
        "connection closed",
        "no route to host",
        "temporary failure",
        "network is unreachable",
        "rsync error",
        "error",
    )
    for line in reversed(lines):
        lower = line.lower()
        if any(token in lower for token in interesting_tokens):
            return line[:240]
    return lines[-1][:240]


def _job_status_payload(job, **overrides):
    payload = {
        "job_id": job["job_id"],
        "status": "queued",
        "created_at": job["created_at"],
        "updated_at": utc_now_iso(),
        "source_count": len(job["sources"]),
        "sources": job["sources"],
        "host": job["host"],
        "user": job["user"],
        "remote_path": job["remote_path"],
        "port": job["port"],
        "attempt": 0,
        "max_retries": job["max_retries"],
        "pid": None,
        "last_error": "",
        "last_summary": "",
        "log_path": job["log_path"],
        "job_path": job["job_path"],
    }
    payload.update(overrides)
    return payload


def create_transfer_job(
    sources,
    host,
    user,
    remote_path,
    port,
    *,
    jobs_dir=None,
    max_retries=DEFAULT_MAX_RETRIES,
):
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    source_paths = _coerce_sources(sources)
    jobs_dir_path = ensure_jobs_dir(jobs_dir)
    job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    job_path = jobs_dir_path / f"{job_id}.job.json"
    status_path = jobs_dir_path / f"{job_id}.status.json"
    log_path = jobs_dir_path / f"{job_id}.log"

    job_payload = {
        "job_id": job_id,
        "created_at": utc_now_iso(),
        "host": host,
        "user": user,
        "remote_path": remote_path,
        "port": port,
        "sources": [str(path) for path in source_paths],
        "max_retries": int(max(0, max_retries)),
        "job_path": str(job_path),
        "status_path": str(status_path),
        "log_path": str(log_path),
    }
    _write_json_atomic(job_path, job_payload)
    _write_json_atomic(status_path, _job_status_payload(job_payload))
    return job_payload


def _update_status(status_path, job, **overrides):
    payload = _job_status_payload(job)
    existing_path = Path(status_path)
    if existing_path.exists():
        try:
            payload.update(_read_json(existing_path))
        except Exception:
            pass
    payload.update(overrides)
    payload["updated_at"] = utc_now_iso()
    _write_json_atomic(status_path, payload)
    return payload


def launch_transfer_job(job, *, python_executable=None):
    python_cmd = python_executable or sys.executable
    log_path = Path(job["log_path"])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_handle = log_path.open("a", encoding="utf-8")
    process = subprocess.Popen(
        [python_cmd, str(Path(__file__).resolve()), "--run", job["job_path"]],
        stdin=subprocess.DEVNULL,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        close_fds=True,
    )
    log_handle.close()
    _update_status(
        job["status_path"],
        job,
        status="queued",
        pid=process.pid,
        last_summary="Background transfer job started.",
    )
    return process.pid


def list_transfer_jobs(*, jobs_dir=None, limit=8):
    jobs_dir_path = ensure_jobs_dir(jobs_dir)
    entries = []
    for status_file in jobs_dir_path.glob("*.status.json"):
        try:
            payload = _read_json(status_file)
        except Exception:
            continue
        payload.setdefault("status_path", str(status_file))
        entries.append(payload)
    entries.sort(key=lambda item: str(item.get("updated_at") or ""), reverse=True)
    return entries[: max(1, int(limit))]


def _append_log(log_path, message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with Path(log_path).open("a", encoding="utf-8") as handle:
        handle.write(f"[{timestamp}] {message}\n")


def _run_logged_command(command, log_path):
    _append_log(log_path, f"Running: {' '.join(shlex.quote(part) for part in command)}")
    with Path(log_path).open("a", encoding="utf-8") as handle:
        process = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    summary = _summarize_output_text(_tail_text(log_path))
    return process.returncode, summary


def _is_retryable_failure(exit_code, summary):
    lowered = str(summary or "").lower()
    permanent_failure_tokens = (
        "permission denied",
        "host key verification failed",
        "operation not permitted",
        "authentication failed",
    )
    if any(token in lowered for token in permanent_failure_tokens):
        return False
    if int(exit_code) in RETRYABLE_EXIT_CODES:
        return True
    retryable_tokens = (
        "connection timed out",
        "timed out",
        "broken pipe",
        "connection reset",
        "connection closed",
        "no route to host",
        "temporary failure",
        "network is unreachable",
        "resource temporarily unavailable",
    )
    return any(token in lowered for token in retryable_tokens)


def _sleep_with_status(job, attempt, wait_seconds, summary):
    _update_status(
        job["status_path"],
        job,
        status="retrying",
        attempt=attempt,
        pid=os.getpid(),
        last_summary=summary,
        last_error=summary,
        retry_in_seconds=wait_seconds,
    )
    time.sleep(wait_seconds)


def run_transfer_job(job_path):
    job = _read_json(job_path)
    status_path = job["status_path"]
    log_path = job["log_path"]
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    _append_log(log_path, f"Job {job['job_id']} started.")
    _update_status(
        status_path,
        job,
        status="running",
        pid=os.getpid(),
        last_summary="Preparing remote directory.",
        retry_in_seconds=0,
    )

    max_attempts = max(1, int(job.get("max_retries", DEFAULT_MAX_RETRIES)) + 1)
    for attempt in range(1, max_attempts + 1):
        _update_status(
            status_path,
            job,
            status="running",
            attempt=attempt,
            pid=os.getpid(),
            last_summary=f"Attempt {attempt}: preparing remote directory.",
            retry_in_seconds=0,
        )
        mkdir_exit, mkdir_summary = _run_logged_command(
            build_ssh_mkdir_command(job["host"], job["user"], job["remote_path"], job["port"]),
            log_path,
        )
        if mkdir_exit != 0:
            summary = mkdir_summary or f"ssh exited with code {mkdir_exit}"
            if attempt < max_attempts and _is_retryable_failure(mkdir_exit, summary):
                wait_seconds = min(
                    DEFAULT_RETRY_DELAY_SECONDS * attempt,
                    DEFAULT_RETRY_BACKOFF_CAP_SECONDS,
                )
                _append_log(log_path, f"Retrying after SSH setup failure in {wait_seconds}s: {summary}")
                _sleep_with_status(job, attempt, wait_seconds, summary)
                continue
            _update_status(
                status_path,
                job,
                status="failed",
                attempt=attempt,
                pid=os.getpid(),
                last_summary=summary,
                last_error=summary,
                retry_in_seconds=0,
            )
            return mkdir_exit

        _update_status(
            status_path,
            job,
            status="running",
            attempt=attempt,
            pid=os.getpid(),
            last_summary=f"Attempt {attempt}: transferring with rsync.",
            last_error="",
            retry_in_seconds=0,
        )
        rsync_exit, rsync_summary = _run_logged_command(
            build_rsync_transfer_command(
                job["sources"],
                job["host"],
                job["user"],
                job["remote_path"],
                job["port"],
            ),
            log_path,
        )
        if rsync_exit in (0, 24):
            final_status = "success" if rsync_exit == 0 else "warning"
            summary = (
                "Transfer finished."
                if rsync_exit == 0
                else "Transfer finished with vanished-file warning (rsync code 24)."
            )
            _append_log(log_path, summary)
            _update_status(
                status_path,
                job,
                status=final_status,
                attempt=attempt,
                pid=os.getpid(),
                last_summary=summary,
                last_error="",
                retry_in_seconds=0,
            )
            return rsync_exit

        summary = rsync_summary or f"rsync exited with code {rsync_exit}"
        if attempt < max_attempts and _is_retryable_failure(rsync_exit, summary):
            wait_seconds = min(
                DEFAULT_RETRY_DELAY_SECONDS * attempt,
                DEFAULT_RETRY_BACKOFF_CAP_SECONDS,
            )
            _append_log(log_path, f"Retrying after rsync failure in {wait_seconds}s: {summary}")
            _sleep_with_status(job, attempt, wait_seconds, summary)
            continue

        _update_status(
            status_path,
            job,
            status="failed",
            attempt=attempt,
            pid=os.getpid(),
            last_summary=summary,
            last_error=summary,
            retry_in_seconds=0,
        )
        return rsync_exit

    summary = "Transfer stopped after exhausting retry attempts."
    _append_log(log_path, summary)
    _update_status(
        status_path,
        job,
        status="failed",
        attempt=max_attempts,
        pid=os.getpid(),
        last_summary=summary,
        last_error=summary,
        retry_in_seconds=0,
    )
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run a detached DELFIN SSH transfer job.")
    parser.add_argument("--run", dest="job_path", help="Path to the job JSON file.")
    args = parser.parse_args(argv)
    if not args.job_path:
        parser.error("--run is required")
    return int(run_transfer_job(args.job_path) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
