"""Installation self-check for DELFIN (``delfin doctor``).

Complements ``delfin-agent doctor`` (delfin/agent/doctor.py), which checks
the agent stack; this one checks what a calculation needs.

Each check is a small function returning a :class:`CheckResult`.  The
checks are deliberately cheap: they only locate binaries via ``PATH``
(``shutil.which``), run ``--version`` probes, touch the scratch
directory, and look at files/keys on disk.  They NEVER start a
computation and NEVER open a network connection.

Status values:
    "ok"      — the prerequisite works
    "missing" — not installed / not configured (informational)
    "broken"  — present but does not work (e.g. found but fails to start)

``exit_code`` treats only ``"broken"`` as failure: a machine without
SLURM or without a KIT key is a valid setup, a machine whose ORCA
crashes on ``--version`` is not.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "CheckResult",
    "check_orca",
    "check_xtb",
    "check_openmpi",
    "check_scratch_dir",
    "check_slurm",
    "check_kit_toolbox_key",
    "check_docs_index",
    "check_command_isolation",
    "check_no_exported_key",
    "run_all",
    "exit_code",
]

OK = "ok"
MISSING = "missing"
BROKEN = "broken"

# Keep version probes snappy: these must never hang the doctor.
_PROBE_TIMEOUT_S = 15


@dataclass
class CheckResult:
    """Outcome of a single installation check.

    Fields:
        name: short identifier of the check (e.g. "orca").
        status: "ok" | "missing" | "broken".
        detail: one-line human-readable finding (never a secret value).
        fix_hint: how to fix a "missing"/"broken" result; "" when ok.
    """

    name: str
    status: str
    detail: str = ""
    fix_hint: str = field(default="")


def _probe(binary: str, args: list[str]) -> tuple[str | None, str | None]:
    """Run ``binary args``; return (resolved_path, first_output_line).

    Returns ``(None, reason)`` when the binary is not on PATH or the
    probe fails.  No network is involved — ``--version`` is local.
    """
    path = shutil.which(binary)
    if path is None:
        return None, f"{binary} not found on PATH"
    try:
        proc = subprocess.run(
            [path, *args],
            capture_output=True,
            text=True,
            timeout=_PROBE_TIMEOUT_S,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        return path, f"{binary} found at {path} but failed to run: {exc}"
    if proc.returncode != 0:
        return path, (
            f"{binary} found at {path} but exited with code {proc.returncode}"
        )
    first_line = (proc.stdout or proc.stderr or "").strip().splitlines()
    return path, first_line[0] if first_line else ""


def _check_binary(name: str, binary: str, version_flag: str,
                  missing_hint: str) -> CheckResult:
    """Shared logic for the binary checks (ORCA, xtb, mpirun, sinfo)."""
    path, info = _probe(binary, [version_flag])
    if path is None:
        return CheckResult(name, MISSING, info or f"{binary} not found",
                           missing_hint)
    if info and info.startswith(f"{binary} found at"):
        return CheckResult(name, BROKEN, info,
                           f"check the {binary} installation at {path}")
    return CheckResult(name, OK, info or f"{binary} at {path}")


def _tool_health(name: str):
    """What DELFIN's own tool resolver knows about ``name``.

    ``qm_health.check_tool`` finds a tool the way DELFIN runs it (PATH,
    ``ORCA_BINARY``/``ORCA_PATH``, the qm_tools directory, dangling links)
    and knows each program's probe. ORCA has no ``--version``: given one it
    tries to open it as an input file and exits 2, so probing it that way
    reported every working ORCA as broken (review of the first run,
    2026-09-16). ``depth="runs"`` starts only what can be started safely.
    """
    from delfin import qm_health
    return qm_health.check_tool(name, depth="runs", timeout=_PROBE_TIMEOUT_S)


def _check_qm_tool(name: str, missing_hint: str) -> CheckResult:
    try:
        health = _tool_health(name)
    except Exception as exc:
        return CheckResult(name, BROKEN, f"could not check {name}: {exc}",
                           "report this as a DELFIN bug")
    level = str(getattr(health, "level", "absent") or "absent")
    path = str(getattr(health, "path", "") or "")
    why = str(getattr(health, "why", "") or "")
    fix = str(getattr(health, "fix", "") or "")
    if level == "absent":
        return CheckResult(name, MISSING, why or f"{name} not found", missing_hint)
    if level == "fail":
        return CheckResult(name, BROKEN,
                           f"{name} at {path}: {why}" if path else why,
                           fix or f"check the {name} installation")
    version = str(getattr(health, "version", "") or "")
    detail = f"{name} at {path}" + (f" ({version})" if version else "")
    return CheckResult(name, OK, detail)


def check_orca() -> CheckResult:
    """ORCA found where DELFIN looks for it."""
    return _check_qm_tool(
        "orca",
        "install ORCA or add its bin directory to PATH "
        "(e.g. module load chem/orca), or set ORCA_BINARY",
    )


def check_xtb() -> CheckResult:
    """xtb found and answers its version probe."""
    return _check_qm_tool(
        "xtb",
        "install xtb or add it to PATH (conda install xtb, or module load)",
    )


def check_openmpi() -> CheckResult:
    """mpirun (OpenMPI) found and answers ``mpirun --version``."""
    return _check_binary(
        "openmpi", "mpirun", "--version",
        "install OpenMPI (mpirun) or load the matching module",
    )


def check_scratch_dir(scratch_dir: str | os.PathLike | None = None) -> CheckResult:
    """Scratch directory exists and is writable.

    Resolution order: explicit argument, ``$DELFIN_SCRATCH``,
    ``$TMPDIR``, then the system temp directory.
    """
    if scratch_dir is None:
        scratch_dir = os.environ.get("DELFIN_SCRATCH") or os.environ.get(
            "TMPDIR") or tempfile.gettempdir()
    path = Path(scratch_dir)
    if not path.exists():
        return CheckResult(
            "scratch_dir", MISSING, f"{path} does not exist",
            f"create it (mkdir -p {path}) or point $DELFIN_SCRATCH "
            "at an existing directory",
        )
    probe = path / ".delfin_doctor_probe"
    try:
        probe.write_text("", encoding="utf-8")
        probe.unlink()
    except OSError as exc:
        return CheckResult(
            "scratch_dir", BROKEN, f"{path} exists but is not writable: {exc}",
            "fix permissions or choose another scratch directory "
            "via $DELFIN_SCRATCH",
        )
    return CheckResult("scratch_dir", OK, f"{path} exists and is writable")


def check_slurm() -> CheckResult:
    """SLURM visible via ``sinfo`` — optional, absence is not an error."""
    return _check_binary(
        "slurm", "sinfo", "--version",
        "optional: install SLURM client tools if you submit to a cluster",
    )


def check_kit_toolbox_key() -> CheckResult:
    """KIT-Toolbox API key configured — presence only, never the value.

    Looked up the way DELFIN looks it up: the environment first, then the
    credential store ``~/.delfin/credentials.json``. Reading only the
    environment reported the key missing for every user who had stored it.
    """
    try:
        from delfin.agent.credentials import load_credential
        present = bool(load_credential("KIT_TOOLBOX_API_KEY"))
    except Exception:
        present = bool(os.environ.get("KIT_TOOLBOX_API_KEY"))
    if present:
        return CheckResult("kit_toolbox_key", OK,
                           "KIT_TOOLBOX_API_KEY is configured")
    return CheckResult(
        "kit_toolbox_key", MISSING, "KIT_TOOLBOX_API_KEY is not configured",
        "store it with `delfin-agent credentials set KIT_TOOLBOX_API_KEY` "
        "or export KIT_TOOLBOX_API_KEY",
    )


def check_command_isolation() -> CheckResult:
    """What confines the agent's commands on THIS host.

    ok: bubblewrap works, or Landlock with the socket guard (seccomp user
    notification). missing: neither -- a locked session then refuses shell
    commands and an unattended run is held by path checks only. The setting
    ``agent.bash_isolation = "off"`` is named: it switches the isolation off
    by choice.
    """
    try:
        from delfin.agent import api_client as A
        from delfin.agent import socket_guard as SG
        bwrap = bool(A._bwrap_functional())
        landlock = bool(A._landlock_functional())
        guard = bool(SG.available())
    except Exception as exc:  # pragma: no cover - import environment issue
        return CheckResult("command_isolation", BROKEN,
                           f"could not probe isolation: {exc}",
                           "check your DELFIN installation (delfin.agent)")
    try:
        from delfin.user_settings import load_settings
        setting = str(((load_settings() or {}).get("agent") or {})
                      .get("bash_isolation", "auto") or "auto")
    except Exception:
        setting = "auto"
    parts = [f"bubblewrap {'works' if bwrap else 'unavailable'}",
             f"Landlock {'available' if landlock else 'unavailable'}",
             f"socket guard {'available' if guard else 'unavailable'}",
             f"agent.bash_isolation={setting}"]
    detail = ", ".join(parts)
    if setting.strip().lower() == "off":
        return CheckResult(
            "command_isolation", MISSING, detail + " (switched off)",
            "unattended runs are not isolated while agent.bash_isolation is "
            "\"off\"; set it to \"auto\" in ~/.delfin_settings.json unless a "
            "workflow needs raw bash")
    if bwrap or (landlock and guard):
        return CheckResult("command_isolation", OK, detail)
    return CheckResult(
        "command_isolation", MISSING, detail,
        "without bubblewrap or Landlock (Linux 5.13+) a locked session "
        "refuses shell commands and an unattended run is held by path "
        "checks only; install bubblewrap or use a newer kernel")


def check_no_exported_key() -> CheckResult:
    """No provider key exported in the environment -- names only."""
    try:
        from delfin.agent import process_guard
        names = process_guard.exported_provider_keys()
        advice = process_guard.exported_key_advice(names) if names else ""
    except Exception:
        names, advice = [], ""
    if not names:
        return CheckResult("exported_keys", OK,
                           "no provider key exported in the environment")
    return CheckResult("exported_keys", MISSING,
                       f"{', '.join(names)} exported in the environment", advice)


def check_docs_index() -> CheckResult:
    """DELFIN doc-search index file exists and is readable JSON."""
    try:
        # Lazy import: the indexer is only needed for the path helper.
        from delfin.doc_server.indexer import get_default_index_path
        idx = get_default_index_path()
    except Exception as exc:  # pragma: no cover - import environment issue
        return CheckResult(
            "docs_index", BROKEN, f"could not resolve index path: {exc}",
            "check your DELFIN installation (delfin.doc_server)",
        )
    if not idx.exists():
        return CheckResult(
            "docs_index", MISSING, f"no index at {idx}",
            "run the doc indexer (delfin-docs-index) to build it",
        )
    try:
        json.loads(idx.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return CheckResult(
            "docs_index", BROKEN, f"index at {idx} is unreadable: {exc}",
            "rebuild it with delfin-docs-index",
        )
    return CheckResult("docs_index", OK, f"index present at {idx}")


def run_all(scratch_dir: str | os.PathLike | None = None) -> list[CheckResult]:
    """Run every check and return the results in a stable order.

    Never raises: a check that explodes is reported as "broken".
    """
    checks = [
        check_orca,
        check_xtb,
        check_openmpi,
        lambda: check_scratch_dir(scratch_dir),
        check_slurm,
        check_kit_toolbox_key,
        check_command_isolation,
        check_no_exported_key,
        check_docs_index,
    ]
    results: list[CheckResult] = []
    for check in checks:
        try:
            results.append(check())
        except Exception as exc:  # defensive: the doctor must complete
            results.append(CheckResult(
                getattr(check, "__name__", "unknown"), BROKEN,
                f"check crashed: {exc}", "report this as a DELFIN bug",
            ))
    return results


def exit_code(results: list[CheckResult]) -> int:
    """0 when nothing is "broken" (missing prerequisites are OK), else 1."""
    return 1 if any(r.status == BROKEN for r in results) else 0
