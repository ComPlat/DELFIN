"""ChampionSpecialist — generic specialist that delegates to a historical commit.

Each instance wraps one champion-commit's `smiles_converter.smiles_to_xyz` via:
  - Archive readback (cached, instant) when name+archive XYZ exists
  - Subprocess to a git-worktree containing the champion's code (live)

Per nature_project/15_HYBRID_PATH_FINAL.md Phase 2B.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


# Internal infrastructure locations; override via env vars.  No personal/
# absolute path is baked in (these are research-only delegation helpers).
DEFAULT_WORKTREE_BASE = Path(
    os.environ.get("DELFIN_WORKTREE_BASE", "worktrees")
)

# Default location of archived XYZ outputs per commit
DEFAULT_ARCHIVE_BASE = Path(
    os.environ.get("DELFIN_ARCHIVE_BASE", "xyz_archive")
)


@dataclass
class ChampionSpecialist:
    """Specialist that delegates to a historical champion commit's code.

    Attributes:
        id: unique specialist identifier (e.g. 'multi_sigma_123a130')
        source_commit: 7-char git short-hash of the champion commit
        archive_subdir: subdir in xyz_archive containing this champion's XYZ files
                        (used for cached readback)
        worktree_path: optional override for worktree location.
                       Defaults to DEFAULT_WORKTREE_BASE / source_commit.
        archive_base: optional override for archive base dir.

    Usage:
        spec = ChampionSpecialist(
            id="multi_sigma_123a130",
            source_commit="123a130",
            archive_subdir="123a130",
        )
        xyz, err = spec.convert(smiles, name="01-Fe_CO_3_NHC_2")
    """
    id: str
    source_commit: str
    archive_subdir: str
    worktree_path: Optional[Path] = None
    archive_base: Optional[Path] = None
    subprocess_timeout: float = 180.0
    subprocess_retries: int = 1

    def __post_init__(self):
        if self.worktree_path is None:
            self.worktree_path = DEFAULT_WORKTREE_BASE / self.source_commit
        if self.archive_base is None:
            self.archive_base = DEFAULT_ARCHIVE_BASE

    def is_available(self) -> bool:
        """True if the worktree exists and has smiles_converter."""
        sc = self.worktree_path / "delfin" / "smiles_converter.py"
        return sc.exists()

    def _archive_lookup(self, name: str) -> Optional[str]:
        """Try to read pre-computed XYZ from xyz_archive/<subdir>/<name>.xyz.

        Returns XYZ string if file exists and non-empty, else None.
        Strips DELFIN's multi-frame header — returns first frame only for now.
        """
        if not name:
            return None
        archive_dir = self.archive_base / self.archive_subdir
        candidate = archive_dir / f"{name}.xyz"
        if not candidate.exists() or candidate.stat().st_size == 0:
            return None
        try:
            text = candidate.read_text()
        except OSError:
            return None
        return text or None

    def _subprocess_convert(self, smiles: str, **kwargs) -> Tuple[Optional[str], Optional[str]]:
        """Run champion's code via subprocess in worktree.

        Inserts worktree_path at front of sys.path so the champion's delfin
        package shadows any installed version.
        """
        if not self.is_available():
            return None, f"specialist {self.id}: worktree missing at {self.worktree_path}"

        # Build a small driver script
        driver = (
            f"import sys; "
            f"sys.path.insert(0, {repr(str(self.worktree_path))}); "
            f"from delfin.smiles_converter import smiles_to_xyz; "
            f"xyz, err = smiles_to_xyz({repr(smiles)}); "
            f"print('---XYZ-START---'); "
            f"print(xyz if xyz is not None else ''); "
            f"print('---XYZ-END---'); "
            f"print('---ERR---'); "
            f"print(err if err is not None else '')"
        )
        # Prevent recursion: subprocess must NOT trigger router again
        env = {**os.environ, "PYTHONPATH": str(self.worktree_path)}
        env.pop("DELFIN_ENSEMBLE_ROUTER", None)

        # Retry-loop for robustness against transient failures
        last_error: Optional[str] = None
        for attempt in range(self.subprocess_retries + 1):
            try:
                result = subprocess.run(
                    [sys.executable, "-c", driver],
                    capture_output=True, text=True,
                    timeout=self.subprocess_timeout,
                    env=env,
                )
            except subprocess.TimeoutExpired:
                last_error = (f"specialist {self.id}: timeout after "
                              f"{self.subprocess_timeout}s (attempt {attempt+1})")
                continue
            except Exception as exc:
                last_error = f"specialist {self.id}: subprocess error {exc}"
                continue

            if result.returncode != 0:
                short_stderr = (result.stderr.strip().splitlines()[-1]
                                if result.stderr else "no stderr")
                last_error = (f"specialist {self.id}: exit {result.returncode}: "
                              f"{short_stderr[:200]}")
                continue

            return self._parse_driver_output(result.stdout)

        return None, last_error or f"specialist {self.id}: all retries failed"

    @staticmethod
    def _parse_driver_output(out: str) -> Tuple[Optional[str], Optional[str]]:
        """Parse subprocess stdout with XYZ-START/XYZ-END/ERR markers."""
        lines = out.splitlines()
        try:
            xyz_start = lines.index("---XYZ-START---") + 1
            xyz_end = lines.index("---XYZ-END---", xyz_start)
            err_start = lines.index("---ERR---", xyz_end) + 1
        except ValueError:
            return None, "specialist: malformed driver output"
        xyz = "\n".join(lines[xyz_start:xyz_end]).strip()
        err = "\n".join(lines[err_start:]).strip()
        return (xyz if xyz else None, err if err else None)

    def convert(
        self,
        smiles: str,
        name: Optional[str] = None,
        prefer_archive: bool = True,
        **kwargs,
    ) -> Tuple[Optional[str], Optional[str]]:
        """Convert SMILES → XYZ via this champion.

        Args:
            smiles: SMILES string
            name: optional name (e.g. '01-Fe_CO_3_NHC_2') for archive lookup
            prefer_archive: if True and name provided, try archive readback first
            **kwargs: forwarded to live conversion (currently unused)

        Returns:
            (xyz, error) tuple matching smiles_to_xyz API.
        """
        if prefer_archive and name:
            cached = self._archive_lookup(name)
            if cached:
                return cached, None
        return self._subprocess_convert(smiles, **kwargs)
