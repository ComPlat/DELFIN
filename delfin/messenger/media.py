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


@contextmanager
def directory(workspace, sid, box):
    if not workspace or not Path(workspace).is_absolute():
        raise ValueError('Sitzung hat keinen gültigen Arbeitsordner.')
    fd = os.open(workspace, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    parts = ['.delfin-messenger', hashlib.sha256(sid.encode()).hexdigest()[:24], box]
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
                raise ValueError('Zu viele Dateiübertragungen. Bitte später erneut versuchen.')
            name = filename(text)
            data = bytearray()
            if action == 'download.begin':
                if text != name:
                    raise ValueError('Nur Dateinamen aus dateien sind zulässig.')
                with directory(workspace, sid, 'outbox') as (fd, _):
                    src = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
                    with os.fdopen(src, 'rb') as stream:
                        info = os.fstat(stream.fileno())
                        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1 or info.st_size > MAX_FILE:
                            raise ValueError('Keine reguläre Datei oder größer als 8 MiB.')
                        data = stream.read(MAX_FILE + 1)
                        if len(data) > MAX_FILE:
                            raise ValueError('Datei größer als 8 MiB.')
            token = secrets.token_hex(16)
            self.transfers[token] = dict(sid=sid, generation=generation, workspace=workspace,
                time=now, name=name, data=data, offset=0, upload=action == 'upload.begin')
            return {'ok': True, 'id': token, 'name': name, 'size': len(data),
                    'sha256': hashlib.sha256(data).hexdigest()}
        transfer = self.transfers.get(token)
        if not transfer or (transfer['sid'], transfer['generation'], transfer['workspace']) != (sid, generation, workspace):
            raise ValueError('Dateiübertragung abgelaufen oder gehört zu einer anderen Sitzung.')
        if action == 'upload.chunk' and transfer['upload']:
            data = base64.b64decode(text, validate=True)
            if len(data) > CHUNK or len(transfer['data']) + len(data) > MAX_FILE:
                self.transfers.pop(token)
                raise ValueError('Datei größer als 8 MiB.')
            transfer['data'].extend(data)
            return {'ok': True}
        if action == 'upload.finish' and transfer['upload']:
            self.transfers.pop(token)
            if hashlib.sha256(transfer['data']).hexdigest() != text:
                raise ValueError('Dateiübertragung unvollständig.')
            with directory(workspace, sid, 'inbox') as (fd, path):
                if len(os.listdir(fd)) >= 100:
                    raise ValueError('Anhangablage voll (100 Dateien). Bitte zuerst aufräumen.')
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
        raise ValueError('Ungültiger Übertragungsschritt.')


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
        raise ValueError('Messenger-Dateiablage voll. Bitte alte Übertragungen aufräumen.')
    folder = root / secrets.token_hex(12)
    folder.mkdir(mode=0o700)
    return folder / filename(name)


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
        raise ValueError('Lokale Spracherkennung ist noch nicht eingerichtet (DELFIN_WHISPER_MODEL).')
    python = os.environ.get('DELFIN_WHISPER_PYTHON', '') or sys.executable
    try:
        result = subprocess.run([python, '-I', '-c', _TRANSCRIBER, str(path), model_path],
                                capture_output=True, text=True, timeout=120, check=True)
        import json
        return str(json.loads(result.stdout))[:4000]
    except (OSError, subprocess.SubprocessError, ValueError):
        raise ValueError('Lokale Spracherkennung fehlgeschlagen; faster-whisper und Modell prüfen.') from None
