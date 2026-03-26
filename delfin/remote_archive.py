"""Read-only browsing helpers for a configured remote archive over SSH."""

from __future__ import annotations

import base64
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


def write_remote_text_file(host, user, remote_path, port, relative_path, text_content):
    """Write a UTF-8 text file below the configured remote root."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    if not rel:
        raise ValueError("Remote file path is empty.")
    payload_b64 = base64.b64encode(str(text_content or "").encode("utf-8")).decode("ascii")
    script = textwrap.dedent("""
        import base64, os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        payload_b64 = sys.argv[3]
        target = os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Path escapes root.")
        parent = os.path.dirname(target)
        if not os.path.isdir(parent):
            raise SystemExit("Parent directory not found.")
        data = base64.b64decode(payload_b64.encode("ascii"))
        with open(target, "wb") as handle:
            handle.write(data)
        print(target)
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(value)) for value in (root, rel, payload_b64)
    )
    return _run_ssh_command(host, user, port, command).strip()


REMOTE_TEXT_CHUNK_BYTES = 128_000


def read_remote_text_chunk(host, user, remote_path, port, relative_path,
                           start_byte=0, chunk_bytes=REMOTE_TEXT_CHUNK_BYTES):
    """Read a byte range from a remote file and return text + metadata."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    script = textwrap.dedent("""
        import base64, json, os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        start = int(sys.argv[3])
        length = int(sys.argv[4])
        target = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.isfile(target):
            raise SystemExit("File not found.")
        size = os.path.getsize(target)
        safe_start = max(0, min(start, max(0, size - 1)))
        with open(target, "rb") as f:
            f.seek(safe_start)
            data = f.read(length)
        print(json.dumps({
            "text": data.decode("utf-8", errors="ignore"),
            "data_b64": base64.b64encode(data).decode("ascii"),
            "size": size,
            "chunk_start": safe_start,
            "chunk_end": safe_start + len(data),
        }))
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel, int(start_byte), int(chunk_bytes))
    )
    output = _run_ssh_command(host, user, port, command)
    return json.loads(output)


def search_remote_file(host, user, remote_path, port, relative_path,
                       query, max_matches=500):
    """Search for a fixed string in a remote file. Returns [(start, end), ...]."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    q_bytes = query.encode("utf-8", errors="ignore")
    if not q_bytes:
        return []
    script = textwrap.dedent("""
        import json, os, re, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        query = sys.argv[3]
        max_m = int(sys.argv[4])
        target = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.isfile(target):
            raise SystemExit("File not found.")
        q_lower = query.lower().encode("utf-8", errors="ignore")
        q_len = len(q_lower)
        spans = []
        overlap = max(q_len - 1, 256)
        carry = b""
        file_offset = 0
        with open(target, "rb") as f:
            while True:
                block = f.read(131072)
                if not block:
                    break
                data = carry + block
                lower = data.lower()
                search_from = 0
                while True:
                    idx = lower.find(q_lower, search_from)
                    if idx < 0:
                        break
                    abs_start = file_offset - len(carry) + idx
                    if abs_start >= 0:
                        spans.append([abs_start, abs_start + q_len])
                        if len(spans) >= max_m:
                            print(json.dumps(spans))
                            sys.exit(0)
                    search_from = idx + 1
                if len(data) > overlap:
                    carry = data[-overlap:]
                else:
                    carry = data
                file_offset += len(block)
        print(json.dumps(spans))
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel, query, int(max_matches))
    )
    output = _run_ssh_command(host, user, port, command)
    return [tuple(pair) for pair in json.loads(output)]


def line_col_from_remote_file_offset(host, user, remote_path, port, relative_path,
                                     byte_pos, chunk_bytes=REMOTE_TEXT_CHUNK_BYTES):
    """Return 1-based (line, column) for a byte offset in a remote file."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    target_pos = max(0, int(byte_pos))
    read_bytes = max(1, int(chunk_bytes))
    script = textwrap.dedent("""
        import json, os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        target = max(0, int(sys.argv[3]))
        block_size = max(1, int(sys.argv[4]))
        path = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([path, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.isfile(path):
            raise SystemExit("File not found.")
        line = 1
        last_nl = -1
        seen = 0
        with open(path, "rb") as handle:
            while seen < target:
                need = min(block_size, target - seen)
                if need <= 0:
                    break
                chunk = handle.read(need)
                if not chunk:
                    break
                idx = 0
                while True:
                    nl = chunk.find(b"\n", idx)
                    if nl < 0:
                        break
                    line += 1
                    last_nl = seen + nl
                    idx = nl + 1
                seen += len(chunk)
        col = target + 1 if last_nl < 0 else target - last_nl
        print(json.dumps({"line": line, "col": col}))
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel, target_pos, read_bytes)
    )
    output = _run_ssh_command(host, user, port, command)
    payload = json.loads(output)
    return int(payload.get("line", 0) or 0), int(payload.get("col", 0) or 0)


def list_remote_folder_files(host, user, remote_path, port, relative_path,
                             filename_pattern, recursive=False, max_content_bytes=2 * 1024 * 1024,
                             folder_relative_paths=None, include_current=False):
    """Collect matching files for remote table extraction in one SSH call.

    When *folder_relative_paths* is provided, only those folders are scanned.
    Otherwise the helper scans either the current folder (*include_current=True*)
    or the direct child folders below *relative_path*.

    Returns list of dicts: ``{folder_name, relative_path, content}``, where
    ``relative_path`` points to the matched file relative to the configured
    remote root.
    """
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    folder_paths = []
    for folder_path in folder_relative_paths or []:
        normalized = normalize_remote_relative_path(folder_path)
        if normalized:
            folder_paths.append(normalized)
    folders_json = json.dumps(folder_paths)
    script = textwrap.dedent("""
        import fnmatch, json, os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        pattern = sys.argv[3]
        recursive = sys.argv[4] == "1"
        max_bytes = int(sys.argv[5])
        selected_folders = json.loads(sys.argv[6])
        include_current = sys.argv[7] == "1"
        base = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([base, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.isdir(base):
            raise SystemExit("Directory not found.")
        results = []

        def _folder_label(folder_path):
            rel_to_base = os.path.relpath(folder_path, base).replace(os.sep, "/")
            if rel_to_base == ".":
                base_name = os.path.basename(base.rstrip(os.sep))
                return base_name or "/"
            return rel_to_base

        folder_targets = []
        if selected_folders:
            for folder_rel in selected_folders:
                folder_path = os.path.realpath(os.path.join(root, folder_rel))
                if os.path.commonpath([folder_path, root]) != root:
                    continue
                if not os.path.isdir(folder_path):
                    continue
                folder_targets.append((_folder_label(folder_path), folder_path))
        elif include_current:
            folder_targets.append((_folder_label(base), base))
        else:
            try:
                folder_targets = sorted(
                    [
                        (entry.name, os.path.realpath(os.path.join(base, entry.name)))
                        for entry in os.scandir(base)
                        if entry.is_dir(follow_symlinks=False)
                    ],
                    key=lambda item: item[0].lower(),
                )
            except OSError:
                folder_targets = []

        for folder_name, folder_path in folder_targets:
            found = []
            if recursive:
                for dirpath, _dirs, files in os.walk(folder_path):
                    for fname in files:
                        if fnmatch.fnmatch(fname, pattern):
                            found.append(os.path.join(dirpath, fname))
            else:
                direct = os.path.join(folder_path, pattern)
                if os.path.isfile(direct):
                    found.append(direct)
                else:
                    for fname in os.listdir(folder_path):
                        if fnmatch.fnmatch(fname, pattern):
                            fp = os.path.join(folder_path, fname)
                            if os.path.isfile(fp):
                                found.append(fp)
            if not found:
                results.append({"folder_name": folder_name, "relative_path": "", "content": ""})
                continue
            for fp in sorted(found):
                try:
                    rel_to_root = os.path.relpath(fp, root).replace(os.sep, "/")
                except ValueError:
                    rel_to_root = os.path.basename(fp)
                try:
                    with open(fp, "rb") as f:
                        raw = f.read(max_bytes)
                    content = raw.decode("utf-8", errors="ignore")
                except Exception:
                    content = ""
                results.append({
                    "folder_name": folder_name,
                    "relative_path": rel_to_root,
                    "content": content,
                })
        print(json.dumps(results))
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (
            root, rel, filename_pattern,
            "1" if recursive else "0",
            int(max_content_bytes),
            folders_json,
            "1" if include_current else "0",
        )
    )
    output = _run_ssh_command(host, user, port, command)
    return json.loads(output)


def remote_mkdir(host, user, remote_path, port, relative_path, folder_name):
    """Create a new folder on the remote server."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    safe_name = str(folder_name or "").strip()
    if not safe_name or "/" in safe_name or "\x00" in safe_name:
        raise ValueError("Invalid folder name.")
    script = textwrap.dedent("""
        import os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        name = sys.argv[3]
        base = root if not relative else os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([base, root]) != root:
            raise SystemExit("Path escapes root.")
        target = os.path.join(base, name)
        if os.path.exists(target):
            raise SystemExit("Already exists.")
        os.makedirs(target, exist_ok=False)
        print("OK")
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel, safe_name)
    )
    _run_ssh_command(host, user, port, command)


def remote_rename(host, user, remote_path, port, relative_path, new_name):
    """Rename a file or folder on the remote server."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    if not rel:
        raise ValueError("Cannot rename the root directory.")
    safe_name = str(new_name or "").strip()
    if not safe_name or "/" in safe_name or "\x00" in safe_name:
        raise ValueError("Invalid name.")
    script = textwrap.dedent("""
        import os, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        new_name = sys.argv[3]
        target = os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.exists(target):
            raise SystemExit("Source not found.")
        parent = os.path.dirname(target)
        dest = os.path.join(parent, new_name)
        if os.path.exists(dest):
            raise SystemExit("Destination already exists.")
        os.rename(target, dest)
        print("OK")
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel, safe_name)
    )
    _run_ssh_command(host, user, port, command)


def remote_duplicate(host, user, remote_path, port, relative_path):
    """Duplicate a file or folder on the remote server (appends _copy suffix)."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    rel = normalize_remote_relative_path(relative_path)
    if not rel:
        raise ValueError("Cannot duplicate the root directory.")
    script = textwrap.dedent("""
        import os, shutil, sys
        root = os.path.realpath(sys.argv[1])
        relative = sys.argv[2]
        target = os.path.realpath(os.path.join(root, relative))
        if os.path.commonpath([target, root]) != root:
            raise SystemExit("Path escapes root.")
        if not os.path.exists(target):
            raise SystemExit("Source not found.")
        parent = os.path.dirname(target)
        base = os.path.basename(target)
        name, ext = os.path.splitext(base)
        dest = os.path.join(parent, name + "_copy" + ext)
        counter = 2
        while os.path.exists(dest):
            dest = os.path.join(parent, name + f"_copy{counter}" + ext)
            counter += 1
        if os.path.isdir(target):
            shutil.copytree(target, dest)
        else:
            shutil.copy2(target, dest)
        print(os.path.basename(dest))
    """).strip()
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, rel)
    )
    output = _run_ssh_command(host, user, port, command)
    return output.strip()


def _remote_copy_or_move(host, user, remote_path, port, relative_paths, target_relative_path, mode):
    """Copy or move files/folders inside the remote archive."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    if not relative_paths:
        raise ValueError("Nothing selected.")
    normalized_mode = str(mode or "").strip().lower()
    if normalized_mode not in {"copy", "move"}:
        raise ValueError("Unsupported remote operation mode.")
    validated = []
    for rp in relative_paths:
        rel = normalize_remote_relative_path(rp)
        if not rel:
            raise ValueError("Cannot operate on the root directory.")
        validated.append(rel)
    target_rel = normalize_remote_relative_path(target_relative_path)
    script = textwrap.dedent("""
        import json, os, shutil, sys

        root = os.path.realpath(sys.argv[1])
        items = json.loads(sys.argv[2])
        target_rel = sys.argv[3]
        mode = sys.argv[4]
        target_dir = root if not target_rel else os.path.realpath(os.path.join(root, target_rel))
        if os.path.commonpath([target_dir, root]) != root:
            raise SystemExit("Target escapes root.")
        if not os.path.isdir(target_dir):
            raise SystemExit("Target folder not found.")

        def next_available_path(dest):
            if not os.path.exists(dest):
                return dest
            parent = os.path.dirname(dest)
            base = os.path.basename(dest)
            name, ext = os.path.splitext(base)
            counter = 2
            while True:
                candidate = os.path.join(parent, f"{name}_{counter}{ext}")
                if not os.path.exists(candidate):
                    return candidate
                counter += 1

        done, errors = [], []
        for rel in items:
            src = os.path.realpath(os.path.join(root, rel))
            if os.path.commonpath([src, root]) != root:
                errors.append(rel + ": escapes root")
                continue
            if not os.path.exists(src):
                errors.append(rel + ": not found")
                continue
            try:
                if os.path.isdir(src) and os.path.commonpath([target_dir, src]) == src:
                    raise ValueError("Cannot paste a folder into itself.")
                dest = next_available_path(os.path.join(target_dir, os.path.basename(src)))
                if mode == "copy":
                    if os.path.isdir(src):
                        shutil.copytree(src, dest)
                    else:
                        shutil.copy2(src, dest)
                else:
                    shutil.move(src, dest)
                done.append(rel)
            except Exception as exc:
                errors.append(rel + ": " + str(exc))
        print(json.dumps({"done": done, "errors": errors}))
    """).strip()
    items_json = json.dumps(validated)
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, items_json, target_rel, normalized_mode)
    )
    output = _run_ssh_command(host, user, port, command)
    return json.loads(output)


def remote_copy(host, user, remote_path, port, relative_paths, target_relative_path=""):
    """Copy files or folders into another remote archive folder."""
    return _remote_copy_or_move(
        host, user, remote_path, port, relative_paths, target_relative_path, "copy"
    )


def remote_move(host, user, remote_path, port, relative_paths, target_relative_path=""):
    """Move files or folders into another remote archive folder."""
    return _remote_copy_or_move(
        host, user, remote_path, port, relative_paths, target_relative_path, "move"
    )


def remote_delete(host, user, remote_path, port, relative_paths):
    """Delete files or folders on the remote server."""
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    root = build_remote_absolute_path(remote_path, "")
    if not relative_paths:
        raise ValueError("Nothing to delete.")
    validated = []
    for rp in relative_paths:
        rel = normalize_remote_relative_path(rp)
        if not rel:
            raise ValueError("Cannot delete the root directory.")
        validated.append(rel)
    script = textwrap.dedent("""
        import json, os, shutil, sys
        root = os.path.realpath(sys.argv[1])
        items = json.loads(sys.argv[2])
        deleted, errors = [], []
        for rel in items:
            target = os.path.realpath(os.path.join(root, rel))
            if os.path.commonpath([target, root]) != root:
                errors.append(rel + ": escapes root")
                continue
            if not os.path.exists(target):
                errors.append(rel + ": not found")
                continue
            try:
                if os.path.isdir(target):
                    shutil.rmtree(target)
                else:
                    os.unlink(target)
                deleted.append(rel)
            except Exception as e:
                errors.append(rel + ": " + str(e))
        print(json.dumps({"deleted": deleted, "errors": errors}))
    """).strip()
    items_json = json.dumps(validated)
    command = "python3 -c " + shlex.quote(script) + " " + " ".join(
        shlex.quote(str(v)) for v in (root, items_json)
    )
    output = _run_ssh_command(host, user, port, command)
    return json.loads(output)
