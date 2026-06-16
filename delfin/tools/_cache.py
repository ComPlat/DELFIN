"""Content-addressed step cache — skip recompute across runs.

A step's result is keyed by ``(step_name, semantic params, input-geometry
hash)``.  On a hit the cached output files (geometry, log, artifacts) are
materialised into the requested work dir and a reconstructed
:class:`~delfin.tools._types.StepResult` is returned, so an agent re-running
variants of a workflow does not repay for identical sub-steps.

Opt-in: call :func:`run_step_cached` instead of ``run_step``.  The cache lives
under ``$DELFIN_CACHE_DIR`` (default ``~/.delfin/cache/steps``).  Only successful
results are cached.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.tools._types import StepResult, StepStatus


def _default_cache_dir() -> Path:
    env = os.environ.get("DELFIN_CACHE_DIR")
    base = Path(env) if env else (Path.home() / ".delfin" / "cache")
    return base / "steps"


def _jsonable(obj: Any) -> Any:
    try:
        json.dumps(obj)
        return obj
    except (TypeError, ValueError):
        return str(obj)


def cache_key(step_name: str, kwargs: Dict[str, Any], geometry=None) -> str:
    """Stable content hash of a step invocation.

    Internal pipeline keys (``_``-prefixed) are excluded; the geometry is hashed
    by file content so the same structure reuses the same entry.
    """
    norm = {k: _jsonable(v) for k, v in kwargs.items() if not k.startswith("_")}
    payload = json.dumps({"step": step_name, "kwargs": norm}, sort_keys=True, default=str)
    h = hashlib.sha256(payload.encode("utf-8"))
    if geometry:
        p = Path(geometry)
        if p.is_file():
            h.update(b"\x00geom\x00")
            h.update(p.read_bytes())
    return h.hexdigest()[:24]


class StepCache:
    """File-backed cache of successful step results."""

    def __init__(self, base: Optional[Path | str] = None):
        self.base = Path(base) if base else _default_cache_dir()
        self.base.mkdir(parents=True, exist_ok=True)

    def _entry(self, key: str) -> Path:
        return self.base / key

    def get(self, key: str, *, into: Path | str) -> Optional[StepResult]:
        """Return a cached result materialised into *into*, or ``None`` on a miss."""
        entry = self._entry(key)
        meta_path = entry / "meta.json"
        if not meta_path.is_file():
            return None
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None

        into = Path(into)
        into.mkdir(parents=True, exist_ok=True)
        for fname in meta.get("files", []):
            src = entry / fname
            if src.is_file():
                shutil.copy2(src, into / fname)

        artifacts = {k: into / fn for k, fn in meta.get("artifacts", {}).items()}
        return StepResult(
            step_name=meta["step_name"],
            status=StepStatus(meta.get("status", "success")),
            geometry=(into / meta["geometry"]) if meta.get("geometry") else None,
            output_file=(into / meta["output_file"]) if meta.get("output_file") else None,
            work_dir=into,
            data=dict(meta.get("data", {})),
            artifacts=artifacts,
        )

    def put(self, key: str, result: StepResult) -> None:
        """Store a *successful* step result under *key* (no-op otherwise)."""
        if not result.ok:
            return
        entry = self._entry(key)
        entry.mkdir(parents=True, exist_ok=True)
        files = []

        def _store(path) -> Optional[str]:
            p = Path(path)
            if not p.is_file():
                return None
            shutil.copy2(p, entry / p.name)
            files.append(p.name)
            return p.name

        meta: Dict[str, Any] = {
            "step_name": result.step_name,
            "status": result.status.value,
            "data": result.data,
        }
        if result.geometry:
            meta["geometry"] = _store(result.geometry)
        if result.output_file:
            meta["output_file"] = _store(result.output_file)
        meta["artifacts"] = {k: n for k, v in result.artifacts.items()
                             if (n := _store(v)) is not None}
        meta["files"] = files
        (entry / "meta.json").write_text(json.dumps(meta, indent=2, default=str),
                                         encoding="utf-8")


_CACHE: Optional[StepCache] = None


def get_cache() -> StepCache:
    """The process-wide default step cache."""
    global _CACHE
    if _CACHE is None:
        _CACHE = StepCache()
    return _CACHE


def run_step_cached(
    step_name: str,
    *,
    geometry=None,
    cores: int = 1,
    work_dir=None,
    cache: Optional[StepCache] = None,
    **kwargs: Any,
) -> StepResult:
    """``run_step`` with a content-addressed cache in front (opt-in).

    On a hit the result is returned with ``data["_cache"] = "hit"``; on a miss it
    runs normally and the successful result is stored.
    """
    from delfin.tools._runner import run_step

    cache = cache or get_cache()
    key = cache_key(step_name, kwargs, geometry)

    wdir = (Path(work_dir) if work_dir
            else Path(tempfile.mkdtemp(prefix=f"{step_name}_", dir=Path.cwd())))
    wdir.mkdir(parents=True, exist_ok=True)

    hit = cache.get(key, into=wdir)
    if hit is not None:
        hit.data = {**hit.data, "_cache": "hit", "_cache_key": key}
        return hit

    result = run_step(step_name, geometry=geometry, cores=cores, work_dir=wdir, **kwargs)
    cache.put(key, result)
    if result.ok:
        result.data.setdefault("_cache", "miss")
        result.data.setdefault("_cache_key", key)
    return result


__all__ = ["cache_key", "StepCache", "get_cache", "run_step_cached"]
