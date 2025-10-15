"""Global singleton job manager for coordinating all DELFIN workflows.

This module provides a centralized job manager that ensures:
1. All workflows share the same resource pool
2. PAL (core count) is never exceeded globally
3. No double allocation of cores when ox/red workflows run in parallel
"""

from __future__ import annotations
from typing import Any, Dict, Optional
import threading
import os
import json

from delfin.common.logging import get_logger
from delfin.dynamic_pool import DynamicCorePool

logger = get_logger(__name__)


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

        logger.info("Global job manager singleton created")

    def initialize(self, config: Dict[str, Any]) -> None:
        """Initialize the manager with configuration.

        Args:
            config: DELFIN configuration dictionary containing PAL, maxcore, etc.
        """
        self.config = config
        self.total_cores = max(1, int(config.get('PAL', 1)))
        self.total_memory = int(config.get('maxcore', 1000)) * self.total_cores

        # Calculate max concurrent jobs
        pal_jobs = config.get('pal_jobs')
        if pal_jobs is None or pal_jobs <= 0:
            # Default: use half the cores as max jobs, capped at 4
            self.max_jobs = max(1, min(4, self.total_cores // 2))
        else:
            self.max_jobs = max(1, int(pal_jobs))

        # Create the shared dynamic pool
        if self.pool is not None:
            logger.warning("Reinitializing global job pool - shutting down existing pool")
            self.pool.shutdown()

        self.pool = DynamicCorePool(
            total_cores=self.total_cores,
            total_memory_mb=self.total_memory,
            max_jobs=self.max_jobs
        )

        pool_id = id(self.pool)
        print(
            f"╔═══════════════════════════════════════════════════════════════╗\n"
            f"║ GLOBAL JOB MANAGER INITIALIZED                                ║\n"
            f"║   Pool ID: {pool_id:<50} ║\n"
            f"║   Total cores: {self.total_cores:<46} ║\n"
            f"║   Max concurrent jobs: {self.max_jobs:<38} ║\n"
            f"║   Total memory: {self.total_memory} MB{' ' * (42 - len(str(self.total_memory)))} ║\n"
            f"╚═══════════════════════════════════════════════════════════════╝"
        )

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
            'pool_status': pool_status,
        }

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

    manager = get_global_manager()
    if manager.is_initialized():
        return

    try:
        manager.initialize(config)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to initialize global manager from %s: %s", env_var, exc)
