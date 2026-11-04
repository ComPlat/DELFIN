"""Central MLIP model-cache environment handling.

One place decides which directory the model weights land in and whether
it is usable.  Every MLIP entry point (calculator factory, ExtOpt
launcher, SLURM dispatch) calls :func:`ensure_mlp_cache_env` instead of
duplicating the logic.

Priority order for the cache directory:

1. ``AIMNET_CACHE_DIR`` — the upstream aimnetcentral override.  Set by
   the user, it wins and is only checked for existence/writability.
2. ``DELFIN_MLP_CACHE_DIR`` — DELFIN's own key, for deployments that
   want one persistent cache shared by several backends.
3. ``$HOME/mlip-envs/aimnet-cache`` — only when it already exists
   (Phase A2 convention on the KIT bwUniCluster); never created, so
   nothing is invented on other hosts.
4. ``~/.cache/aimnet`` — upstream default.  Used only when writable;
   on clusters with a read-only home this fails and the error names the
   working alternative instead of an OSError traceback.

No host- or user-specific absolute path is hardcoded; the fallback
under (3) is only honoured when the directory already exists.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

__all__ = ["ensure_mlp_cache_env", "resolve_mlp_cache_dir"]


def resolve_mlp_cache_dir(env: dict | None = None) -> tuple[str, str]:
    """Return ``(cache_dir, source)`` without touching the filesystem.

    ``source`` is one of ``AIMNET_CACHE_DIR``, ``DELFIN_MLP_CACHE_DIR``,
    ``existing-home-fallback`` or ``default``.
    """
    env = os.environ if env is None else env
    explicit = (env.get("AIMNET_CACHE_DIR") or "").strip()
    if explicit:
        return explicit, "AIMNET_CACHE_DIR"
    delfin_key = (env.get("DELFIN_MLP_CACHE_DIR") or "").strip()
    if delfin_key:
        return delfin_key, "DELFIN_MLP_CACHE_DIR"
    home_fallback = Path.home() / "mlip-envs" / "aimnet-cache"
    if home_fallback.is_dir():
        return str(home_fallback), "existing-home-fallback"
    return str(Path.home() / ".cache" / "aimnet"), "default"


def ensure_mlp_cache_env(env: dict | None = None) -> str:
    """Make sure a writable model-cache dir exists and export it.

    Sets ``AIMNET_CACHE_DIR`` in the (real or given) environment unless
    the caller already provided it.  Raises ``RuntimeError`` with a
    user-facing hint when nothing writable is found — the hint names the
    override key, not a traceback.
    """
    env = os.environ if env is None else env
    cache_dir, source = resolve_mlp_cache_dir(env)
    path = Path(cache_dir)
    if path.exists():
        if not path.is_dir():
            raise RuntimeError(
                f"MLIP model cache path {cache_dir} is not a directory "
                f"(from {source}); set AIMNET_CACHE_DIR to a directory"
            )
        if not os.access(path, os.W_OK):
            raise RuntimeError(
                f"MLIP model cache {cache_dir} is not writable (from "
                f"{source}); set AIMNET_CACHE_DIR to a writable directory"
            )
    else:
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise RuntimeError(
                f"cannot create MLIP model cache {cache_dir} (from "
                f"{source}): {exc}. On a cluster with a read-only home, "
                "set AIMNET_CACHE_DIR to a writable path (e.g. a scratch "
                "or project directory) and reuse it across jobs."
            ) from exc
    # prove writability once; upstream download failures then mean what
    # they say instead of surfacing as write errors mid-download
    try:
        with tempfile.NamedTemporaryFile(dir=path, prefix=".delfin-mlp-cache-"):
            pass
    except OSError as exc:
        raise RuntimeError(
            f"MLIP model cache {cache_dir} is not writable (from {source}): "
            f"{exc}. Set AIMNET_CACHE_DIR to a writable directory."
        ) from exc
    env["AIMNET_CACHE_DIR"] = cache_dir
    return cache_dir
