"""Abstract base classes for the DELFIN Dashboard job backend."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple


@dataclass
class SubmitResult:
    """Result of a job submission (mirrors subprocess.CompletedProcess API)."""
    returncode: int
    stdout: str = ''
    stderr: str = ''


@dataclass
class JobInfo:
    """Information about a single job (backend-agnostic)."""
    job_id: Any
    name: str = ''
    mode: str = ''
    status: str = ''
    submit_time: str = ''
    start_time: str = ''
    pal: int = 0
    maxcore: int = 0
    time_limit: str = ''
    job_dir: str = ''
    extra: dict = field(default_factory=dict)


class JobBackend(ABC):
    """Abstract interface for job submission/management backends."""

    @abstractmethod
    def submit_delfin(self, job_dir, job_name, mode='delfin',
                      time_limit='48:00:00', pal=40, maxcore=6000,
                      override=None, build_mult=None) -> SubmitResult:
        """Submit a DELFIN calculation job."""

    @abstractmethod
    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=40, maxcore=6000) -> SubmitResult:
        """Submit a standalone ORCA job."""

    def submit_turbomole(self, job_dir, job_name, module='ridft',
                         time_limit='48:00:00', nprocs=40, mem_per_cpu=6000,
                         para_arch='SMP') -> SubmitResult:
        """Submit a TURBOMOLE job (only supported on SLURM)."""
        return SubmitResult(1, stderr='TURBOMOLE not supported on this backend')

    @abstractmethod
    def list_jobs(self) -> List[JobInfo]:
        """Return list of active jobs."""

    @abstractmethod
    def cancel_job(self, job_id) -> Tuple[bool, str]:
        """Cancel a job. Returns (success, message)."""

    def get_pending_start_times(self) -> list:
        """Return estimated start times for pending jobs (SLURM only)."""
        return []

    @property
    def supports_turbomole(self) -> bool:
        return False

    @property
    def backend_name(self) -> str:
        return self.__class__.__name__
