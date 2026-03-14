"""Read-only browsing helpers for a configured remote archive over SSH."""

from __future__ import annotations

import hashlib
import json
import posixpath
import shlex
import subprocess
import textwrap
from pathlib import Path
from typing import Any

from delfin.ssh_transfer_jobs import (
    DEFAULT_RSYNC_TIMEOUT_SECONDS,
    build_rsync_ssh_command,
    build_ssh_command_base,
    quote_remote_path,
)
from delfin.user_settings import normalize_ssh_transfer_settings

REMOTE_CACHE_DIR_NAME = ".delfin_remote_archive_cache"
TEXT_PREVIEW_MAX_BYTES = 2 * 1024 * 1024


def normalize_remote_relative_path(relative_path):
    """Normalize a path relative to the configured remote root."""
    raw = str(relative_path or "").strip().replace("\\", "/")
    if raw in ("", "/"):
        return ""
    parts = []
    for part in raw.lstrip("/").split("/"):
        piece = str(part or "").strip()
        if not piece or piece == ".":
            continue
        if piece == "..":
            raise ValueError("Remote path cannot leave the configured root.")
        if any(ch in piece for ch in ("\x00", "\n", "\r")):
            raise ValueError("Remote path contains unsupported control characters.")
        parts.append(piece)
    return "/".join(parts)


def build_remote_absolute_path(root_path, relative_path=""):
    """Join a validated relative path to the configured remote root."""
    root = str(root_path or "").strip().replace("\\", "/")
    if not root:
        raise ValueError("Remote root path is empty.")
    if not root.startswith("/"):
        raise ValueError("Remote root path must be absolute.")
    rel = normalize_remote_relative_path(relative_path)
    normalized_root = posixpath.normpath(root)
    return normalized_root if not rel else posixpath.normpath(posixpath.join(normalized_root, rel))


def get_remote_cache_dir(base_dir=None):
    cache_dir = Path(base_dir).expanduser() if base_dir else (Path.home() / REMOTE_CACHE_DIR_NAME)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def remote_cache_namespace(host, user, port, remote_path):
    """Create a stable cache namespace for one configured remote root."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    digest = hashlib.sha256(
        f"{user}@{host}:{int(port)}:{posixpath.normpath(remote_path)}".encode("utf-8")
    ).hexdigest()
    safe_host = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in host)
    safe_user = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in user)
    return f"{safe_user}@{safe_host}_{digest[:12]}"


def get_cached_local_path(host, user, remote_path, port, relative_path="", *, cache_dir=None):
    """Return the local cache path for a remote file or directory."""
    rel = normalize_remote_relative_path(relative_path)
    namespace = remote_cache_namespace(host, user, port, remote_path)
    base = get_remote_cache_dir(cache_dir) / namespace
    if not rel:
        return base
    return base.joinpath(*rel.split("/"))


def build_remote_list_command(root_path, relative_path="", sort_mode="name"):
    """Build a remote Python command that returns directory entries as JSON."""
    root = build_remote_absolute_path(root_path, "")
    rel = normalize_remote_relative_path(relative_path)
    sort_key = str(sort_mode or "name").strip().lower()
    if sort_key not in {"name", "date_desc", "date_asc"}:
        raise ValueError("Unsupported remote sort mode.")

    script = textwrap.dedent(
        """
        import json, os, pathlib, stat, sys

        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        sort_mode = sys.argv[3]
        target = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Remote path escapes configured root.")
        if not os.path.isdir(target):
            raise SystemExit("Remote folder not found.")

        entries = []
        with os.scandir(target) as iterator:
            for entry in iterator:
                try:
                    if entry.is_symlink():
                        continue
                    st = entry.stat(follow_symlinks=False)
                except OSError:
                    continue
                rel_path = os.path.relpath(os.path.join(target, entry.name), root).replace(os.sep, "/")
                if rel_path == ".":
                    rel_path = ""
                is_dir = stat.S_ISDIR(st.st_mode)
                entries.append(
                    {
                        "name": entry.name,
                        "relative_path": rel_path,
                        "is_dir": bool(is_dir),
                        "size": 0 if is_dir else int(st.st_size),
                        "mtime": float(st.st_mtime),
                        "suffix": pathlib.Path(entry.name).suffix.lower(),
                    }
                )

        if sort_mode == "date_desc":
            entries.sort(key=lambda item: (not item["is_dir"], -item["mtime"], item["name"].lower()))
        elif sort_mode == "date_asc":
            entries.sort(key=lambda item: (not item["is_dir"], item["mtime"], item["name"].lower()))
        else:
            entries.sort(key=lambda item: (not item["is_dir"], item["name"].lower()))

        print(
            json.dumps(
                {
                    "root_path": root,
                    "relative_path": relative,
                    "entries": entries,
                }
            )
        )
        """
    ).strip()
    return "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(value) for value in (root, rel, sort_key)
    )


def build_remote_text_preview_command(root_path, relative_path="", max_bytes=TEXT_PREVIEW_MAX_BYTES):
    """Build a remote Python command that returns a text preview as JSON."""
    root = build_remote_absolute_path(root_path, "")
    rel = normalize_remote_relative_path(relative_path)
    max_preview_bytes = int(max(1, max_bytes))
    script = textwrap.dedent(
        """
        import json, os, sys

        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        max_bytes = int(sys.argv[3])
        target = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Remote path escapes configured root.")
        if not os.path.isfile(target):
            raise SystemExit("Remote file not found.")

        size = os.path.getsize(target)
        with open(target, "rb") as handle:
            data = handle.read(max_bytes)

        text = data.decode("utf-8", errors="ignore")
        print(
            json.dumps(
                {
                    "text": text,
                    "size": int(size),
                    "truncated": bool(size > len(data)),
                }
            )
        )
        """
    ).strip()
    return "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(value)) for value in (root, rel, max_preview_bytes)
    )


def build_remote_fetch_command(host, user, remote_path, port, relative_path, local_path):
    """Build an rsync command that copies one remote file into the local cache."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    remote_abs = build_remote_absolute_path(remote_path, relative_path)
    destination = Path(local_path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    return [
        "rsync",
        "-a",
        "--partial",
        "--append-verify",
        "--protect-args",
        "--human-readable",
        f"--timeout={DEFAULT_RSYNC_TIMEOUT_SECONDS}",
        "-e",
        build_rsync_ssh_command(host, user, port),
        f"{user}@{host}:{quote_remote_path(remote_abs)}",
        str(destination),
    ]


def _run_ssh_command(host, user, port, remote_command):
    ssh_command = build_ssh_command_base(host, user, port)
    ssh_command.append(str(remote_command))
    result = subprocess.run(
        ssh_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        message = result.stdout.strip() or "SSH command failed."
        raise RuntimeError(message)
    return result.stdout


def list_remote_entries(host, user, remote_path, port, relative_path="", *, sort_mode="name"):
    """Return JSON-like directory listing data from the remote archive."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    command = build_remote_list_command(remote_path, relative_path, sort_mode=sort_mode)
    output = _run_ssh_command(host, user, port, command)
    payload = json.loads(output)
    payload["relative_path"] = normalize_remote_relative_path(payload.get("relative_path", ""))
    return payload


def read_remote_text_preview(host, user, remote_path, port, relative_path="", *, max_bytes=TEXT_PREVIEW_MAX_BYTES):
    """Read a text preview for a remote file without copying it locally."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    command = build_remote_text_preview_command(remote_path, relative_path, max_bytes=max_bytes)
    output = _run_ssh_command(host, user, port, command)
    return json.loads(output)


def fetch_remote_file(host, user, remote_path, port, relative_path="", *, cache_dir=None):
    """Copy one remote file into the local cache and return the local path."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    local_path = get_cached_local_path(
        host,
        user,
        remote_path,
        port,
        relative_path,
        cache_dir=cache_dir,
    )
    command = build_remote_fetch_command(
        host,
        user,
        remote_path,
        port,
        relative_path,
        local_path,
    )
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stdout.strip() or "Remote file fetch failed.")
    return local_path
