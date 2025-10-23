"""Global singleton job manager for coordinating all DELFIN workflows.

This module provides a centralized job manager that ensures:
1. All workflows share the same resource pool
2. PAL (core count) is never exceeded globally
3. No double allocation of cores when ox/red workflows run in parallel
"""

from __future__ import annotations
from typing import Any, Dict, Optional, Tuple
import atexit
import threading
import os
import json

from delfin.common.logging import get_logger
from delfin.dynamic_pool import DynamicCorePool

logger = get_logger(__name__)


def _safe_int(value: Any, default: int) -> int:
    try:
        text = str(value).strip()
    except (TypeError, AttributeError):
        return default
    if text == "":
        return default
    try:
        return int(text)
    except (TypeError, ValueError):
        return default


def _normalize_parallel_token(value: Any, default: str = "auto") -> str:
    token = str(value).strip().lower() if value not in (None, "") else default
    if token in {"no", "false", "off", "0", "disable"}:
        return "disable"
    if token in {"yes", "true", "on", "1", "enable"}:
        return "enable"
    return "auto"


class GlobalJobManager:
    """Singleton manager for all DELFIN computational jobs.

    This manager ensures that all workflows (classic, manually, OCCUPIER)
    share the same resource pool and never exceed configured PAL limits.
    """

    _instance: Optional[GlobalJobManager] = None
    _lock = threading.Lock()

    def __new__(cls):
        """Ensure only one instance exists (Singleton pattern)."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the global manager (only once)."""
        if self._initialized:
            return

        self._initialized = True
        self.pool: Optional[DynamicCorePool] = None
        self.total_cores: int = 1
        self.max_jobs: int = 1
        self.total_memory: int = 1000
        self.config: Dict[str, Any] = {}
        self.parallel_mode: str = "auto"
        self.maxcore_per_job: int = 1000
        self._config_signature: Optional[Tuple[int, int, int, str]] = None
        self._atexit_registered: bool = False

        if not self._atexit_registered:
            atexit.register(self.shutdown)
            self._atexit_registered = True

        logger.info("Global job manager singleton created")

    def initialize(self, config: Dict[str, Any]) -> None:
        """Initialize the manager with configuration.

        Args:
            config: DELFIN configuration dictionary containing PAL, maxcore, etc.
        """
        sanitized = self._sanitize_resource_config(config)
        requested_signature = self._config_signature_value(sanitized)

        if self.pool is not None:
            if self._config_signature == requested_signature:
                logger.info("Global job manager already initialized with matching configuration – reusing existing pool")
                self.config = sanitized
                return

            running_jobs = 0
            try:
                status = self.pool.get_status()
                running_jobs = status.get('running_jobs', 0) if isinstance(status, dict) else 0
            except Exception as exc:  # noqa: BLE001
                logger.debug("Could not inspect active pool status prior to reinitialization: %s", exc)

            if running_jobs > 0:
                logger.warning(
                    "Requested global manager reconfiguration while %d job(s) are still running – keeping existing pool",
                    running_jobs,
                )
                return

            logger.info("Reinitializing global job pool with updated configuration")
            self.pool.shutdown()
            self.pool = None

        self.config = sanitized
        self.total_cores = sanitized['PAL']
        self.total_memory = sanitized['PAL'] * sanitized['maxcore']
        self.maxcore_per_job = sanitized['maxcore']
        self.max_jobs = sanitized['pal_jobs']
        self.parallel_mode = sanitized['parallel_mode']

        self.pool = DynamicCorePool(
            total_cores=self.total_cores,
            total_memory_mb=self.total_memory,
            max_jobs=self.max_jobs,
        )

        pool_id = id(self.pool)
        banner_width = 63

        def _banner_line(text: str = "", *, align: str = "left") -> str:
            trimmed = (text or "")[:banner_width]
            if align == "center":
                inner = trimmed.center(banner_width)
            elif align == "right":
                inner = trimmed.rjust(banner_width)
            else:
                inner = trimmed.ljust(banner_width)
            return f"║{inner}║"

        banner_lines = [
            f"╔{'═' * banner_width}╗",
            _banner_line("GLOBAL JOB MANAGER INITIALIZED", align="center"),
            _banner_line(),
            _banner_line(f"• Pool ID: {pool_id}", align="left"),
            _banner_line(f"• Total cores: {self.total_cores}", align="left"),
            _banner_line(f"• Max concurrent jobs: {self.max_jobs}", align="left"),
            _banner_line(f"• Parallel mode: {self.parallel_mode.upper()}", align="left"),
            _banner_line(f"• Total memory: {self.total_memory} MB", align="left"),
            f"╚{'═' * banner_width}╝",
        ]
        print("\n".join(banner_lines))
        self._config_signature = requested_signature

    def get_pool(self) -> DynamicCorePool:
        """Get the shared dynamic core pool.

        Returns:
            The shared DynamicCorePool instance.

        Raises:
            RuntimeError: If manager hasn't been initialized yet.
        """
        if self.pool is None:
            logger.warning(
                "Global job manager not initialized - this may be a subprocess. "
                "Returning None to allow fallback to local pool."
            )
            raise RuntimeError(
                "Global job manager not initialized. Call initialize(config) first."
            )
        return self.pool

    def is_initialized(self) -> bool:
        """Check if the global manager has been initialized.

        Returns:
            True if initialized, False otherwise.
        """
        return self.pool is not None

    def get_effective_cores_for_workflow(self, workflow_context: str = "") -> int:
        """Calculate effective cores available for a workflow.

        This method accounts for parallel workflows that might be running.
        For example, if ox and red workflows run in parallel, each gets
        half the total cores.

        Args:
            workflow_context: Optional context info for logging

        Returns:
            Number of cores this workflow can use
        """
        # For now, return total cores
        # This will be enhanced to track active workflows
        return self.total_cores

    def shutdown(self) -> None:
        """Shutdown the global manager and clean up resources."""
        if self.pool is not None:
            logger.info("Shutting down global job manager")
            self.pool.shutdown()
            self.pool = None
        self._config_signature = None
        self.config = {}
        self.parallel_mode = "auto"
        self.total_cores = 1
        self.max_jobs = 1
        self.total_memory = 1000
        self.maxcore_per_job = 1000

    def get_status(self) -> Dict[str, Any]:
        """Get current status of the global manager.

        Returns:
            Dictionary with manager status information
        """
        if self.pool is None:
            return {
                'initialized': False,
                'total_cores': self.total_cores,
                'max_jobs': self.max_jobs,
            }

        pool_status = self.pool.get_status()
        return {
            'initialized': True,
            'total_cores': self.total_cores,
            'max_jobs': self.max_jobs,
            'total_memory': self.total_memory,
            'parallel_mode': self.parallel_mode,
            'pool_status': pool_status,
        }

    def ensure_initialized(self, config: Dict[str, Any]) -> None:
        """Initialize the manager if required, otherwise keep the existing pool."""
        sanitized = self._sanitize_resource_config(config)
        requested_sig = self._config_signature_value(sanitized)

        if not self.is_initialized():
            self.initialize(sanitized)
            return

        if self._config_signature != requested_sig:
            logger.info(
                "Global manager already active (current %s, requested %s) – reusing existing pool",
                self._signature_str(self._config_signature),
                self._signature_str(requested_sig),
            )
            return

        # Update cached config to reflect any new auxiliary keys
        self.config.update(sanitized)
        self.parallel_mode = sanitized['parallel_mode']

    @classmethod
    def reset(cls) -> None:
        """Reset the singleton instance (mainly for testing).

        WARNING: This should only be used in tests or when reinitializing
        the entire application.
        """
        with cls._lock:
            if cls._instance is not None and cls._instance.pool is not None:
                cls._instance.pool.shutdown()
            cls._instance = None

    @staticmethod
    def _config_signature_value(config: Dict[str, Any]) -> Tuple[int, int, int, str]:
        return (
            int(config.get('PAL', 1)),
            int(config.get('maxcore', 1000)),
            int(config.get('pal_jobs', 1)),
            str(config.get('parallel_mode', 'auto')),
        )

    def _sanitize_resource_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        cfg: Dict[str, Any] = dict(config or {})

        pal = max(1, _safe_int(cfg.get('PAL'), self.total_cores or 1))
        maxcore = max(256, _safe_int(cfg.get('maxcore'), self.maxcore_per_job or 1000))

        pal_jobs_raw = cfg.get('pal_jobs')
        pal_jobs = _safe_int(pal_jobs_raw, 0)

        parallel_token = _normalize_parallel_token(cfg.get('parallel_workflows', 'auto'))
        if parallel_token == "disable":
            pal_jobs = 1
        if pal_jobs <= 0:
            pal_jobs = max(1, min(4, max(1, pal // 2)))
        pal_jobs = max(1, min(pal_jobs, pal))

        cfg.update({
            'PAL': pal,
            'maxcore': maxcore,
            'pal_jobs': pal_jobs,
            'parallel_mode': parallel_token,
        })
        return cfg

    @staticmethod
    def _signature_str(signature: Optional[Tuple[int, int, int, str]]) -> str:
        if signature is None:
            return "PAL=?, maxcore=?, pal_jobs=?, parallel=?"
        pal, maxcore, pal_jobs, parallel = signature
        return f"PAL={pal}, maxcore={maxcore}, pal_jobs={pal_jobs}, parallel={parallel}"


# Convenience function for getting the global manager
def get_global_manager() -> GlobalJobManager:
    """Get the global job manager instance.

    Returns:
        The GlobalJobManager singleton instance
    """
    return GlobalJobManager()


def bootstrap_global_manager_from_env(env_var: str = "DELFIN_CHILD_GLOBAL_MANAGER") -> None:
    """Initialize the global manager from serialized config in the environment.

    Child OCCUPIER processes spawned by DELFIN use this hook to ensure they
    attach to a properly configured global dynamic pool instead of creating
    ad-hoc local managers.

    Args:
        env_var: Environment variable containing a JSON config snippet.
    """
    payload = os.environ.get(env_var)
    if not payload:
        return

    try:
        config = json.loads(payload)
    except json.JSONDecodeError as exc:  # noqa: BLE001
        logger.warning("Failed to decode %s payload for global manager bootstrap: %s", env_var, exc)
        return

    try:
        manager = get_global_manager()
        manager.ensure_initialized(config)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to initialize global manager from %s: %s", env_var, exc)
