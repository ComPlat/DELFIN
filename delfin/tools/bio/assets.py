"""Explicit local model assets: immutable identity, no implicit downloads."""
from __future__ import annotations

from pathlib import Path
import re

from delfin.tools.bio.artifacts import child_path, file_hash, read_json


def validate_assets(manifest_path: Path, *, engine: str, required_paths=()) -> dict:
    manifest_path = Path(manifest_path).resolve()
    value = read_json(manifest_path)
    if not isinstance(value, dict) or value.get('schema_version') != 1 or value.get('engine') != engine:
        raise ValueError('asset manifest schema/engine mismatch')
    if not isinstance(value.get('files'), list) or not value['files']:
        raise ValueError('asset manifest requires explicit files')
    files = {}
    for entry in value['files']:
        if not isinstance(entry, dict) or set(entry) != {'path', 'sha256'}:
            raise ValueError('each asset requires exactly path and sha256')
        path = child_path(manifest_path.parent, entry['path'])
        relative = str(path.relative_to(manifest_path.parent))
        if relative in files or not re.fullmatch('[0-9a-f]{64}', str(entry['sha256'])):
            raise ValueError('duplicate asset or invalid SHA-256')
        if not path.is_file() or not path.stat().st_size or file_hash(path) != entry['sha256']:
            raise ValueError(f'missing, empty or changed model asset: {relative}')
        files[relative] = {'sha256': entry['sha256'], 'bytes': path.stat().st_size}
    if set(required_paths) - files.keys():
        raise ValueError('required model assets missing: ' + ', '.join(sorted(set(required_paths) - files.keys())))
    return {'engine': engine, 'root': str(manifest_path.parent), 'manifest_sha256': file_hash(manifest_path), 'files': files}
