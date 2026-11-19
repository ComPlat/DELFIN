"""Bounded file exchange, isolated by session, without shared gateway access.

Bytes cross the UID boundary in small chunks. The gateway never receives a
host path to arbitrary user files. Downloads come only from the session's
explicit outbox. Directory descriptors and O_NOFOLLOW reject planted links.
"""
from __future__ import annotations

import base64
import hashlib
import os
import re
import secrets
import stat
import time
from contextlib import contextmanager
from pathlib import Path

MAX_FILE = 8 * 1024 * 1024
CHUNK = 3072


def filename(name):
    name = Path(str(name).replace('\\', '/')).name
    name = re.sub(r'[^\w. -]', '_', name)[:120].strip('. ')
    return name or 'attachment.bin'


#: The messenger's directory in a session workspace.
ROOT = '.delfin-messenger'


def session_hash(sid):
    """The directory name of session ``sid`` under ROOT."""
    return hashlib.sha256(str(sid).encode()).hexdigest()[:24]


def internal_write_refusal(path, sid=''):
    """Why the agent may not write ``path``, or None.

    Input: the resolved target of a write tool or a shell write, and the
    writing session's id. Under ROOT only a plain file directly in the
    session's own ``outbox`` (what the phone can fetch, forum.txt
    included) is writable; the phone's uploads (``inbox``), other
    sessions' directories and everything else there are read-only to the
    agent. Without a session id any session's outbox file is allowed.
    The refusal names the outbox path, the one place to offer a file."""
    parts = Path(str(path)).parts
    if ROOT not in parts:
        return None
    at = parts.index(ROOT)
    base, rest = Path(*parts[:at + 1]), parts[at + 1:]
    own = session_hash(sid) if sid else ''
    if len(rest) == 3 and rest[1] == 'outbox' and rest[0] == (own or rest[0]):
        return None
    outbox = base / (own or '<session>') / 'outbox'
    return (f"'{path}' is inside the messenger's own directory, which is "
            "read-only to the agent (the phone's uploads, other sessions' "
            "files). To offer a file for download on the phone, write it "
            f"as a plain file into {outbox}/ and name the file in your "
            "reply; the user fetches it by that name.")


@contextmanager
def directory(workspace, sid, box):
    if not workspace or not Path(workspace).is_absolute():
        raise ValueError('The session has no valid working folder.')
    fd = os.open(workspace, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    parts = [ROOT, session_hash(sid), box]
    try:
        for part in parts:
            try:
                os.mkdir(part, 0o700, dir_fd=fd)
            except FileExistsError:
                pass
            nxt = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = nxt
        yield fd, Path(workspace).joinpath(*parts)
    finally:
        os.close(fd)


class Exchange:
    def __init__(self):
        self.transfers = {}

    def command(self, sid, generation, workspace, action, *, token='', text=''):
        now = time.monotonic()
        self.transfers = {k: v for k, v in self.transfers.items() if now - v['time'] < 600}
        if action == 'announcement':
            with directory(workspace, sid, 'outbox') as (fd, _):
                try:
                    src = os.open('forum.txt', os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
                except FileNotFoundError:
                    return {'ok': True, 'text': ''}
                with os.fdopen(src, 'rb') as stream:
                    info = os.fstat(stream.fileno())
                    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1 or info.st_size > 1400:
                        return {'ok': True, 'text': ''}
                    data = stream.read(1401)
                if len(data) > 1400:
                    return {'ok': True, 'text': ''}
                return {'ok': True, 'text': data.decode('utf-8'), 'id': hashlib.sha256(data).hexdigest()}
        if action == 'files':
            with directory(workspace, sid, 'outbox') as (fd, path):
                files = []
                for name in sorted(os.listdir(fd))[:100]:
                    info = os.stat(name, dir_fd=fd, follow_symlinks=False)
                    if stat.S_ISREG(info.st_mode) and info.st_nlink == 1 and info.st_size <= MAX_FILE:
                        files.append({'name': name, 'size': info.st_size})
                return {'ok': True, 'files': files, 'outbox': str(path)}
        if action in {'upload.begin', 'download.begin'}:
            if len(self.transfers) >= 4:
                raise ValueError('Too many file transfers. Please try again later.')
            name = filename(text)
            data = bytearray()
            if action == 'download.begin':
                if text != name:
                    raise ValueError('Only file names from files are allowed.')
                with directory(workspace, sid, 'outbox') as (fd, _):
                    src = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
                    with os.fdopen(src, 'rb') as stream:
                        info = os.fstat(stream.fileno())
                        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1 or info.st_size > MAX_FILE:
                            raise ValueError('Not a regular file, or larger than 8 MiB.')
                        data = stream.read(MAX_FILE + 1)
                        if len(data) > MAX_FILE:
                            raise ValueError('File larger than 8 MiB.')
            token = secrets.token_hex(16)
            self.transfers[token] = dict(sid=sid, generation=generation, workspace=workspace,
                time=now, name=name, data=data, offset=0, upload=action == 'upload.begin')
            return {'ok': True, 'id': token, 'name': name, 'size': len(data),
                    'sha256': hashlib.sha256(data).hexdigest()}
        transfer = self.transfers.get(token)
        if not transfer or (transfer['sid'], transfer['generation'], transfer['workspace']) != (sid, generation, workspace):
            raise ValueError('File transfer expired or belongs to another session.')
        if action == 'upload.chunk' and transfer['upload']:
            data = base64.b64decode(text, validate=True)
            if len(data) > CHUNK or len(transfer['data']) + len(data) > MAX_FILE:
                self.transfers.pop(token)
                raise ValueError('File larger than 8 MiB.')
            transfer['data'].extend(data)
            return {'ok': True}
        if action == 'upload.finish' and transfer['upload']:
            self.transfers.pop(token)
            if hashlib.sha256(transfer['data']).hexdigest() != text:
                raise ValueError('File transfer incomplete.')
            with directory(workspace, sid, 'inbox') as (fd, path):
                if len(os.listdir(fd)) >= 100:
                    raise ValueError('Attachment store full (100 files). Please clean up first.')
                name = secrets.token_hex(6) + '-' + transfer['name']
                dest = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=fd)
                with os.fdopen(dest, 'wb') as stream:
                    stream.write(transfer['data'])
                return {'ok': True, 'path': str(path / name)}
        if action == 'download.chunk' and not transfer['upload']:
            start = transfer['offset']
            data = transfer['data'][start:start + CHUNK]
            transfer['offset'] += len(data)
            done = transfer['offset'] == len(transfer['data'])
            if done:
                self.transfers.pop(token)
            return {'ok': True, 'data': base64.b64encode(data).decode(), 'done': done}
        raise ValueError('Invalid transfer step.')


def cache_file(name):
    from . import store
    root = store.state_dir() / 'media'
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    # Retain asynchronous sends, but bound both accumulated storage and lifetime.
    import shutil
    for child in root.iterdir():
        if not child.is_symlink() and child.is_dir() and time.time() - child.stat().st_mtime > 86400:
            shutil.rmtree(child)
    if len(list(root.iterdir())) >= 100:
        raise ValueError('Messenger file store full. Please clean up old transfers.')
    folder = root / secrets.token_hex(12)
    folder.mkdir(mode=0o700)
    return folder / filename(name)


def discard(path) -> None:
    """Remove a received attachment and its own cache folder. Only the file
    went before; the empty folders stayed until the 24 h sweep and counted
    toward the 100-entry limit, so a busy day of voice messages ended in
    "file store full". Never raises."""
    try:
        path = Path(path)
        path.unlink(missing_ok=True)
        if path.parent.parent.name == 'media':
            path.parent.rmdir()
    except OSError:
        pass


#: The transcriber, run as ``python -I -c`` in the whisper interpreter: it
#: imports only faster-whisper, never DELFIN, and reads only the local
#: model directory (``local_files_only``).
_TRANSCRIBER = """
import json, sys
from faster_whisper import WhisperModel
model = WhisperModel(sys.argv[2], device='cpu', compute_type='int8',
                     local_files_only=True)
segments, _ = model.transcribe(sys.argv[1], vad_filter=True)
parts = []
for segment in segments:
    parts.append(segment.text.strip())
    if sum(map(len, parts)) >= 4000:
        break
print(json.dumps(' '.join(parts)[:4000]))
"""


def transcribe(path):
    """Run local speech decoding with a deadline, without downloads/cloud APIs.

    The interpreter is ``DELFIN_WHISPER_PYTHON`` (default: this one): a
    separate environment, because faster-whisper needs numpy 2.x and
    DELFIN keeps numpy below 2. ``-I`` keeps DELFIN, the current directory
    and the user's site-packages out of it.
    """
    import subprocess
    import sys
    model_path = os.environ.get('DELFIN_WHISPER_MODEL', '')
    if not model_path or not Path(model_path).is_dir():
        raise ValueError('Local speech recognition is not set up yet (DELFIN_WHISPER_MODEL).')
    python = os.environ.get('DELFIN_WHISPER_PYTHON', '') or sys.executable
    try:
        result = subprocess.run([python, '-I', '-c', _TRANSCRIBER, str(path), model_path],
                                capture_output=True, text=True, timeout=120, check=True)
        import json
        return str(json.loads(result.stdout))[:4000]
    except (OSError, subprocess.SubprocessError, ValueError):
        raise ValueError('Local speech recognition failed; check faster-whisper and the model.') from None
