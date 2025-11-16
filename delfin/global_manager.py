"""Global singleton job manager for coordinating all DELFIN workflows.

This module provides a centralized job manager that ensures:
1. All workflows share the same resource pool
2. PAL (core count) is never exceeded globally
3. No double allocation of cores when ox/red workflows run in parallel
"""

from __future__ import annotations
from typing import Any, Dict, Optional, Tuple, Callable, List
import atexit
import threading
import os
import json
import signal
import time
from dataclasses import dataclass
import subprocess

from delfin.common.logging import get_logger
from delfin.dynamic_pool import DynamicCorePool

logger = get_logger(__name__)


@dataclass
class _TrackedProcess:
    token: str
    pid: int
    pgid: Optional[int]
    label: str
    process: Any
    start_time: float
    cwd: Optional[str]


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
        self._signal_handler_installed: bool = False
        self._previous_sigint_handler: Optional[Callable] = None
        self._shutdown_requested = threading.Event()
        self._tracked_lock = threading.RLock()
        self._tracked_processes: Dict[str, "_TrackedProcess"] = {}
        self._tracked_counter = 0

        if not self._atexit_registered:
            atexit.register(self.shutdown)
            self._atexit_registered = True

        self._install_signal_handler()
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
        self._terminate_all_processes(reason="shutdown")
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

    # ------------------------------------------------------------------
    # Signal handling and subprocess tracking
    # ------------------------------------------------------------------

    def _install_signal_handler(self) -> None:
        if self._signal_handler_installed:
            return
        if threading.current_thread() is not threading.main_thread():
            logger.debug("Skipping SIGINT handler installation (not main thread)")
            return
        try:
            self._previous_sigint_handler = signal.getsignal(signal.SIGINT)
            signal.signal(signal.SIGINT, self._handle_sigint)
            self._signal_handler_installed = True
            logger.debug("Registered GlobalJobManager SIGINT handler")
        except Exception as exc:  # noqa: BLE001
            logger.debug("Failed to install SIGINT handler: %s", exc)

    def _handle_sigint(self, signum, frame) -> None:  # noqa: ANN001
        if self._shutdown_requested.is_set():
            logger.warning("SIGINT received again – cleanup already in progress")
            return
        self._shutdown_requested.set()
        logger.warning("SIGINT received – aborting active DELFIN jobs and ORCA processes.")
        try:
            self._perform_interrupt_shutdown(signum)
        finally:
            previous = self._previous_sigint_handler
            if previous in (None, signal.SIG_IGN):
                raise KeyboardInterrupt
            if previous is signal.SIG_DFL:
                raise KeyboardInterrupt
            try:
                previous(signum, frame)
            except KeyboardInterrupt:
                raise
            except Exception:  # noqa: BLE001
                logger.debug("Previous SIGINT handler raised", exc_info=True)
                raise KeyboardInterrupt
            else:
                raise KeyboardInterrupt

    def _perform_interrupt_shutdown(self, signum: int) -> None:
        reason = f"signal {signum}"
        self._terminate_all_processes(reason=reason)
        if self.pool is not None:
            try:
                # On interrupt: don't wait for jobs, cancel pending
                self.pool.shutdown(wait=False, cancel_pending=True)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Error while shutting down core pool after SIGINT: %s", exc)

    def register_subprocess(self, process: Any, *, label: str = "", cwd: Optional[str] = None) -> Optional[str]:
        """Register a subprocess for signal-triggered cleanup."""
        if process is None:
            return None
        try:
            pid = process.pid
        except Exception:  # noqa: BLE001
            return None

        try:
            pgid = os.getpgid(pid) if hasattr(os, "getpgid") else None
        except Exception:  # noqa: BLE001
            pgid = None

        with self._tracked_lock:
            self._tracked_counter += 1
            token = f"{pid}:{self._tracked_counter}"
            record = _TrackedProcess(
                token=token,
                pid=pid,
                pgid=pgid,
                label=label or f"pid {pid}",
                process=process,
                start_time=time.time(),
                cwd=str(cwd) if cwd else None,
            )
            self._tracked_processes[token] = record

        logger.debug(
            "Registered subprocess %s (pid=%s, pgid=%s, label=%s)",
            token,
            pid,
            pgid,
            label or "<unnamed>",
        )
        return token

    def unregister_subprocess(self, token: Optional[str]) -> None:
        if not token:
            return
        with self._tracked_lock:
            record = self._tracked_processes.pop(token, None)
        if record:
            logger.debug("Unregistered subprocess %s (pid=%s)", token, record.pid)

    def _terminate_all_processes(self, *, reason: str) -> None:
        with self._tracked_lock:
            records: List[_TrackedProcess] = list(self._tracked_processes.values())
        if not records:
            return

        logger.warning("Terminating %d tracked ORCA process group(s) (%s)", len(records), reason)
        for record in records:
            self._terminate_tracked_process(record)

    def _terminate_tracked_process(self, record: _TrackedProcess) -> None:
        process = record.process
        if process.poll() is not None:
            self.unregister_subprocess(record.token)
            return

        pgid = record.pgid
        label = record.label or f"pid {record.pid}"
        try:
            if pgid is not None and hasattr(os, "killpg"):
                os.killpg(pgid, signal.SIGTERM)
                logger.info("Sent SIGTERM to process group %s (%s)", pgid, label)
            else:
                process.terminate()
                logger.info("Sent terminate signal to %s", label)
        except ProcessLookupError:
            logger.debug("Process %s already exited before termination", label)
            self.unregister_subprocess(record.token)
            return
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to send SIGTERM to %s: %s", label, exc)

        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if pgid is not None and hasattr(os, "killpg"):
                try:
                    os.killpg(pgid, signal.SIGKILL)
                    logger.warning("Sent SIGKILL to process group %s (%s)", pgid, label)
                except Exception as exc:  # noqa: BLE001
                    logger.warning("Failed to SIGKILL %s: %s", label, exc)
            else:
                try:
                    process.kill()
                    logger.warning("Force-killed %s", label)
                except Exception as exc:  # noqa: BLE001
                    logger.warning("Failed to force-kill %s: %s", label, exc)

            try:
                process.wait(timeout=2)
            except Exception:  # noqa: BLE001
                logger.debug("Process %s did not exit after SIGKILL", label, exc_info=True)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Error while waiting for %s termination: %s", label, exc)

        self.unregister_subprocess(record.token)


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
