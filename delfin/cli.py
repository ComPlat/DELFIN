from __future__ import annotations
import argparse
import fnmatch
import os
import re
import shutil
import signal
import sys
import time
from pathlib import Path

from dataclasses import asdict
from typing import Optional
from delfin.common.logging import configure_logging, get_logger, add_file_handler
from delfin.common.paths import get_runtime_dir, resolve_path
from delfin.cluster_utils import auto_configure_resources
from delfin.global_manager import get_global_manager
from .define import create_control_file
from .cleanup import cleanup_all, cleanup_orca
from .config import read_control_file, get_E_ref
from .utils import search_transition_metals, set_main_basisset, calculate_total_electrons_txt, get_git_commit_info
from .orca import run_orca, cleanup_orca_scratch_dir
from .xtb_crest import XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR
from .reporting import (
    generate_summary_report_DELFIN as generate_summary_report,
    generate_esd_report,
)
from .reporting.delfin_collector import save_esd_data_json
from .cli_helpers import _build_parser
from .cli_recalc import setup_recalc_mode, patch_modules_for_recalc
from .cli_banner import print_delfin_banner, validate_required_files
from .pipeline import (
    FileBundle,
    PipelineContext,
    compute_summary,
    interpret_method_alias,
    normalize_input_file,
    run_classic_phase,
    run_manual_phase,
    run_occuper_phase,
    run_esd_phase,
)
from .copy_helpers import extract_preferred_spin, read_occupier_file

logger = get_logger(__name__)

_STEP_FILE_SUFFIXES: tuple[str, ...] = (
    ".inp",
    ".out",
    ".xyz",
    ".gbw",
    ".engrad",
    ".hess",
    ".molden",
    ".prop",
    ".cpcm",
    ".cpcm_corr",
    ".densitiesinfo",
    ".tmp",
    ".tmp1",
    ".tmp2",
    ".opt",
    ".gbw_hs",        # ORCA Hessian matrix file
    ".carthess",      # ORCA Cartesian Hessian
    ".qro",           # ORCA quasi-restricted orbitals
    ".uco",           # ORCA unrestricted corresponding orbitals
    ".uno",           # ORCA unrestricted natural orbitals
    ".unoloc",        # ORCA localized UNOs
    ".unso",          # ORCA unrestricted natural spin orbitals
    ".hostnames",     # ORCA parallel job hostnames
)


def _build_step_bases() -> set[str]:
    """Collect base filenames that DELFIN routinely generates."""
    bundle = asdict(FileBundle())
    bases: set[str] = {Path(value).stem for value in bundle.values()}
    # Additional artifacts not covered by FileBundle
    bases.update({"start"})
    occ_steps = ["initial"]
    occ_steps.extend(f"ox_step_{idx}" for idx in range(1, 4))
    occ_steps.extend(f"red_step_{idx}" for idx in range(1, 4))
    for token in occ_steps:
        bases.add(f"input_{token}_OCCUPIER")
    return bases


_STEP_FILE_BASES = _build_step_bases()
_STEP_FILE_NAMES = {
    f"{base}{suffix}" for base in _STEP_FILE_BASES for suffix in _STEP_FILE_SUFFIXES
}

_STEP_FILE_GLOB_PATTERNS = {
    pattern.format(base=base)
    for base in _STEP_FILE_BASES
    for pattern in (
        "{base}.bibtex",
        "{base}.densities",
        "{base}.property.txt",
        "{base}.property",
        "{base}.hess*",
        "{base}_trj.xyz",
        "{base}_trj*.xyz",
        "{base}.*.*",        # Catches files like red_step_1.B.24.tmp, red_step_1.gpot0.tmp, etc.
        "{base}.bas*",       # ORCA basis set files (bas0, bas1, bas2, etc.)
    )
}

_SAFE_FILE_NAMES: set[str] = {
    "DELFIN.txt",
    "OCCUPIER.txt",
    "ESD.txt",
    "delfin_run.log",
    "occupier.log",
    "s1_state_opt.failed",
    ".delfin_occ_auto_state.json",
    ".delfin_occ_delta",
    ".delfin_done_xtb",
    ".delfin_done_goat",
    ".delfin_done_crest",
    ".delfin_done_xtb_solvator",
    ".qmmm_cache.json",
    "start.txt",
}

_SAFE_FILE_GLOBS: tuple[str, ...] = (
    "delfin_run.log*",
)

_SAFE_DIR_NAMES: set[str] = {
    "CREST",
    "ESD",
    "XTB_SOLVATOR",
    ".orca_scratch",
}

_SAFE_DIR_PATTERNS: tuple[str, ...] = (
    "*_OCCUPIER",
    "*_IMAG",
)

_PREFERRED_INDEX_RE = re.compile(r"(Preferred Index:\s*)(\d+)", re.IGNORECASE)


def _parse_occupier_overrides(raw_tokens: list[str]) -> dict[str, int]:
    """Parse --occupier-override tokens into a mapping of folder -> preferred index."""
    overrides: dict[str, int] = {}
    for raw in raw_tokens:
        token_str = str(raw or "").strip()
        if not token_str:
            continue
        for entry in token_str.split(","):
            candidate = entry.strip()
            if not candidate:
                continue
            if "=" not in candidate:
                raise ValueError(f"Invalid --occupier-override entry '{candidate}' (expected <stage>=<index>).")
            folder, value = candidate.split("=", 1)
            folder = folder.strip()
            if not folder:
                raise ValueError("Missing OCCUPIER folder name in --occupier-override.")
            if not folder.endswith("_OCCUPIER"):
                folder = f"{folder}_OCCUPIER"
            try:
                idx = int(value.strip())
            except ValueError:
                raise ValueError(f"Invalid preferred index '{value}' for {folder} (must be an integer).") from None
            if idx <= 0:
                raise ValueError(f"Preferred index for {folder} must be positive.")
            overrides[folder] = idx
    return overrides


def _override_output_path(folder_name: str, workspace_root: Path) -> Path:
    base = folder_name[:-len("_OCCUPIER")] if folder_name.endswith("_OCCUPIER") else folder_name
    return (workspace_root / f"{base}.out").resolve()


def _apply_occupier_overrides(
    overrides: dict[str, int],
    workspace_root: Path,
    config: dict,
    force_outputs: set[Path],
) -> bool:
    """Rewrite OCCUPIER preferred indices and refresh propagated files."""
    if not overrides:
        return True

    ok = True
    override_map = config.setdefault("_occ_preferred_override", {})
    cache = config.get('_occ_results_runtime')
    for folder_name, preferred_index in overrides.items():
        folder_path = workspace_root / folder_name
        occ_path = folder_path / "OCCUPIER.txt"
        forced_out_path = _override_output_path(folder_name, workspace_root)

        if not occ_path.exists():
            logger.error("Cannot apply --occupier-override: missing %s", occ_path)
            ok = False
            continue

        try:
            content = occ_path.read_text(encoding="utf-8", errors="ignore")
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to read %s: %s", occ_path, exc)
            ok = False
            continue

        new_content, count = _PREFERRED_INDEX_RE.subn(rf"\g<1>{preferred_index}", content, count=1)
        if count == 0:
            logger.error("Preferred Index line not found in %s; override skipped.", occ_path)
            ok = False
        else:
            try:
                if new_content != content:
                    occ_path.write_text(new_content, encoding="utf-8")
                logger.info("[recalc] Set Preferred Index in %s to %d (--occupier-override)", folder_name, preferred_index)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Could not rewrite %s (will still force override at runtime): %s", occ_path, exc)

        if isinstance(cache, dict):
            cache.pop(folder_name, None)

        try:
            result = read_occupier_file(
                str(folder_path),
                "OCCUPIER.txt",
                None,
                None,
                None,
                config,
                verbose=False,
                preferred_index_override=preferred_index,
                skip_file_copy=True,  # In recalc mode, preserve existing geometries
            )
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to refresh OCCUPIER artifacts for %s: %s", folder_name, exc)
            force_outputs.discard(forced_out_path)
            ok = False
            continue

        if not result or len(result) < 3 or result[2] is None:
            logger.error(
                "Override for %s did not yield a valid preferred entry (index %s).",
                folder_name,
                preferred_index,
            )
            ok = False
            continue

        force_outputs.add(forced_out_path)
        try:
            override_map[folder_name] = int(preferred_index)
        except Exception:
            override_map[folder_name] = preferred_index

        # Remove existing INP file to force regeneration with new multiplicity/BrokenSym
        base_name = folder_name.replace("_OCCUPIER", "") if folder_name.endswith("_OCCUPIER") else folder_name
        inp_path = workspace_root / f"{base_name}.inp"
        if inp_path.exists():
            try:
                inp_path.unlink()
                logger.info("[recalc] Removed existing %s to force regeneration with override", inp_path.name)
            except Exception as exc:  # noqa: BLE001
                logger.warning("[recalc] Could not remove %s: %s (will be overwritten)", inp_path.name, exc)

    # Store overrides for post-OCCUPIER reapplication
    if overrides:
        config["_occ_overrides_pending_rewrite"] = dict(overrides)

    return ok

_DIR_SIGNATURES: tuple[tuple[str, ...], ...] = (
    ("XTB.inp", "output_XTB.out", "XTB.xyz"),
    ("XTB_GOAT.inp", "output_XTB_GOAT.out", "XTB_GOAT.globalminimum.xyz"),
    ("XTB_SOLVATOR.inp", "output_XTB_SOLVATOR.out", "XTB_SOLVATOR.solvator.xyz"),
    ("CREST.out", "crest_best.xyz", "initial_opt.xyz"),
)


def _dir_has_known_signature(path: Path) -> bool:
    """Check if a directory contains marker files from helper workflows."""
    try:
        if not path.is_dir():
            return False
    except OSError:
        return False

    for markers in _DIR_SIGNATURES:
        for marker in markers:
            if (path / marker).exists():
                return True
    return False


def _is_workspace_artifact(entry: Path) -> bool:
    """Return True if the given filesystem entry looks like a DELFIN artifact."""
    name = entry.name
    if entry.is_dir():
        if name in _SAFE_DIR_NAMES:
            return True
        if any(fnmatch.fnmatch(name, pattern) for pattern in _SAFE_DIR_PATTERNS):
            return True
        if _dir_has_known_signature(entry):
            return True
        return False

    if name in _SAFE_FILE_NAMES:
        return True
    if any(fnmatch.fnmatch(name, pattern) for pattern in _SAFE_FILE_GLOBS):
        return True
    if name in _STEP_FILE_NAMES:
        return True
    if any(fnmatch.fnmatch(name, pattern) for pattern in _STEP_FILE_GLOB_PATTERNS):
        return True
    return False


def _run_cleanup_subcommand(argv: list[str]) -> int:
    """Handle `delfin cleanup` invocations."""
    parser = argparse.ArgumentParser(
        prog="delfin cleanup",
        description="Remove DELFIN scratch artifacts and optionally stop ORCA jobs.",
    )
    parser.add_argument(
        "--orca",
        action="store_true",
        help="Terminate ORCA subprocesses and remove OCCUPIER scratch folders.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be removed without deleting files or stopping processes.",
    )
    parser.add_argument(
        "--workspace",
        default=".",
        help="Workspace directory (default: current directory).",
    )
    parser.add_argument(
        "--scratch",
        default=None,
        help="Override scratch directory (defaults to runtime dir or workspace).",
    )
    args = parser.parse_args(argv)

    workspace = resolve_path(args.workspace)
    if args.scratch:
        scratch = resolve_path(args.scratch)
    else:
        scratch = get_runtime_dir() if args.workspace == "." else workspace

    if args.orca:
        report = cleanup_orca(workspace, scratch_root=scratch, dry_run=args.dry_run)
        print(f"Workspace: {report['workspace']}")
        print(f"Scratch:   {report['scratch_root']}")
        print(f"ORCA processes detected: {report['processes_found']}")
        for entry in report["terminated_groups"]:
            print(f"  pgid {entry['pgid']}: {entry['status']} (pids={entry['pids']})")
        if report["occuper_dirs_removed"]:
            print("Removed OCCUPIER folders:")
            for path in report["occuper_dirs_removed"]:
                print(f"  {path}")
        if report["scratch_dirs_removed"]:
            print("Removed ORCA scratch directories:")
            for path in report["scratch_dirs_removed"]:
                print(f"  {path}")
        if not args.dry_run:
            print(f"Deleted {report['files_removed']} temporary file(s).")
        else:
            print("Dry run completed — no files deleted.")
        return 0

    removed = cleanup_all(str(scratch), dry_run=args.dry_run)
    if args.dry_run:
        print(f"Dry run: cleanup would affect files under {scratch}")
    else:
        print(f"Removed {removed} temporary file(s) under {scratch}")
    return 0


def _iter_delfin_procs(workspace: Path) -> list[dict]:
    """Return DELFIN processes whose cwd matches workspace."""
    procs: list[dict] = []
    try:
        ws = workspace.resolve()
    except Exception:
        ws = workspace

    proc_root = Path("/proc")
    for pid_dir in proc_root.iterdir():
        if not pid_dir.name.isdigit():
            continue
        try:
            pid = int(pid_dir.name)
        except ValueError:
            continue
        try:
            cwd = Path(os.readlink(pid_dir / "cwd")).resolve()
        except Exception:
            continue
        if cwd != ws:
            continue

        try:
            raw_cmd = (pid_dir / "cmdline").read_bytes()
            if not raw_cmd:
                continue
            tokens = [tok for tok in raw_cmd.split(b"\0") if tok]
            text_tokens = [tok.decode(errors="ignore") for tok in tokens]
        except Exception:
            continue

        joined = " ".join(text_tokens).lower()
        if "delfin" not in joined:
            continue

        try:
            pgid = os.getpgid(pid)
        except Exception:
            pgid = None

        procs.append({
            "pid": pid,
            "pgid": pgid,
            "cmd": " ".join(text_tokens),
            "cwd": str(cwd),
        })
    return procs


def _run_co2_subcommand(argv: list[str]) -> int:
    """Handle `delfin co2` invocations."""
    parser = argparse.ArgumentParser(
        prog="delfin co2",
        description="CO2 Coordinator - Align complex, place CO2, and run orientation/distance scans.",
    )
    parser.add_argument(
        "--define",
        action="store_true",
        help="Generate CONTROL.txt and co2.xyz templates and exit.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files when using --define.",
    )
    parser.add_argument(
        "--recalc",
        action="store_true",
        help="Recalculate only incomplete or missing orientation scan calculations.",
    )
    parser.add_argument(
        "--charge",
        type=int,
        help="Charge for the system (replaces [CHARGE] in template).",
    )
    parser.add_argument(
        "--multiplicity",
        type=int,
        help="Multiplicity for the system (replaces [MULTIPLICITY] in template).",
    )
    parser.add_argument(
        "--solvent",
        type=str,
        help="Solvent for the system (replaces [SOLVENT] in template), e.g., DMF.",
    )
    parser.add_argument(
        "--metal",
        type=str,
        help="Metal symbol (replaces [METAL] in template), e.g., Fe.",
    )
    parser.add_argument(
        "--additions",
        type=str,
        help="Additional ORCA keywords (e.g., '%%SCF BrokenSym 5,5 END' for broken symmetry).",
    )
    args = parser.parse_args(argv)

    # Import CO2_Coordinator6 main function
    try:
        from delfin.co2.CO2_Coordinator6 import main as co2_main, write_default_files
    except ImportError as exc:
        print(f"Error: Could not import CO2 Coordinator module: {exc}")
        return 1

    # Handle --define mode
    if args.define:
        try:
            write_default_files(
                charge=args.charge,
                multiplicity=args.multiplicity,
                solvent=args.solvent,
                metal=args.metal,
                additions=args.additions,
                overwrite=args.force,
            )
            return 0
        except Exception as exc:
            print(f"Error during CO2 template generation: {exc}")
            return 1

    # Normal run mode - execute CO2 coordinator workflow
    try:
        # Set recalc mode environment variable if needed
        if args.recalc:
            os.environ["DELFIN_CO2_RECALC"] = "1"
        co2_main()
        return 0
    except Exception as exc:
        print(f"Error during CO2 Coordinator execution: {exc}")
        import traceback
        traceback.print_exc()
        return 1


def _run_stop_subcommand(argv: list[str]) -> int:
    """Handle `delfin stop` invocations."""
    parser = argparse.ArgumentParser(
        prog="delfin stop",
        description="Send a signal to running DELFIN processes in a workspace (graceful stop).",
    )
    parser.add_argument(
        "--workspace",
        default=".",
        help="Workspace directory (default: current directory).",
    )
    parser.add_argument(
        "--signal",
        default="INT",
        choices=["INT", "TERM", "KILL"],
        help="Signal to send (default: INT).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show matching processes without signaling.",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="After signaling, run DELFIN cleanup in the workspace (remove scratch/OCCUPIER artifacts).",
    )
    parser.add_argument(
        "--wait-seconds",
        type=float,
        default=3.0,
        help="Seconds to wait for processes to exit before cleanup when --cleanup is set (default: 3).",
    )
    args = parser.parse_args(argv)

    workspace = resolve_path(args.workspace)
    self_pid = os.getpid()
    try:
        self_pgid = os.getpgid(self_pid)
    except Exception:
        self_pgid = None
    procs = _iter_delfin_procs(workspace)
    if not procs:
        print(f"No DELFIN process found in {workspace}")
        return 1

    sig_map = {
        "INT": signal.SIGINT,
        "TERM": signal.SIGTERM,
        "KILL": signal.SIGKILL,
    }
    sig = sig_map[args.signal]

    signaled: list[tuple[int, Optional[int], str]] = []
    for proc in procs:
        pid = proc["pid"]
        pgid = proc["pgid"]
        # Avoid signaling the stop command itself
        if pid == self_pid or (self_pgid is not None and pgid == self_pgid):
            continue
        label = f"pid={pid}, pgid={pgid}, cmd={proc['cmd']}"
        print(f"Found DELFIN: {label}")
        if args.dry_run:
            continue
        try:
            if pgid is not None and hasattr(os, "killpg"):
                os.killpg(pgid, sig)
                print(f"  Sent {args.signal} to process group {pgid}")
            else:
                os.kill(pid, sig)
                print(f"  Sent {args.signal} to pid {pid}")
            signaled.append((pid, pgid, label))
        except ProcessLookupError:
            print("  Process already exited")
        except Exception as exc:  # noqa: BLE001
            print(f"  Failed to send signal: {exc}")

    # Optional escalation if processes ignore SIGINT/TERM
    def _still_alive(items: list[tuple[int, Optional[int], str]]) -> list[tuple[int, Optional[int], str]]:
        return [(pid, pgid, lbl) for pid, pgid, lbl in items if Path(f"/proc/{pid}").exists()]

    def _send(sig_to_send: int, entries: list[tuple[int, Optional[int], str]], label: str) -> None:
        for pid, pgid, lbl in entries:
            try:
                if pgid is not None and hasattr(os, "killpg"):
                    os.killpg(pgid, sig_to_send)
                    print(f"  {label}: sent {signal.strsignal(sig_to_send)} to pgid {pgid} ({lbl})")
                else:
                    os.kill(pid, sig_to_send)
                    print(f"  {label}: sent {signal.strsignal(sig_to_send)} to pid {pid} ({lbl})")
            except ProcessLookupError:
                pass
            except Exception as exc:  # noqa: BLE001
                print(f"  {label}: failed to signal {lbl}: {exc}")

    if not args.dry_run and signaled:
        if args.wait_seconds > 0:
            deadline = time.time() + args.wait_seconds
            while time.time() < deadline:
                remaining = _still_alive(signaled)
                if not remaining:
                    break
                time.sleep(0.2)
        remaining = _still_alive(signaled)
        if remaining and args.signal == "INT":
            _send(signal.SIGTERM, remaining, "escalate")
            time.sleep(1.0)
            remaining = _still_alive(remaining)
        if remaining:
            _send(signal.SIGKILL, remaining, "force-kill")

    if args.cleanup and not args.dry_run:
        scratch = get_runtime_dir() if args.workspace == "." else workspace
        try:
            report = cleanup_orca(workspace, scratch_root=scratch, dry_run=False)
            print(f"Workspace: {report['workspace']}")
            print(f"Scratch:   {report['scratch_root']}")
            print(f"ORCA processes detected: {report['processes_found']}")
            for entry in report["terminated_groups"]:
                print(f"  pgid {entry['pgid']}: {entry['status']} (pids={entry['pids']})")
            if report["occuper_dirs_removed"]:
                print("Removed OCCUPIER folders:")
                for path in report["occuper_dirs_removed"]:
                    print(f"  {path}")
            if report["scratch_dirs_removed"]:
                print("Removed ORCA scratch directories:")
                for path in report["scratch_dirs_removed"]:
                    print(f"  {path}")
            print(f"Deleted {report['files_removed']} temporary file(s).")
        except Exception as exc:  # noqa: BLE001
            print(f"Cleanup after stop failed: {exc}")
        try:
            removed = cleanup_all(str(scratch), dry_run=False)
            print(f"Removed {removed} temporary file(s) under {scratch}")
        except Exception as exc:  # noqa: BLE001
            print(f"Cleanup_all failed: {exc}")
    return 0


def _normalize_input_file(config, control_path: Path) -> str:
    return normalize_input_file(config, control_path)


def _safe_keep_set(control_path: Path) -> set[str]:
    """Determine filenames that must be preserved during purge."""
    keep: set[str] = {control_path.name}
    if control_path.exists():
        try:
            cfg = read_control_file(str(control_path))
            input_entry = cfg.get('input_file')
            if input_entry:
                keep.add(Path(input_entry).name)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not parse CONTROL.txt while preparing purge: %s", exc)
    # Always keep fallback input.txt if CONTROL is missing or invalid
    keep.add("input.txt")
    # start.txt must always be purged, never kept
    keep.discard("start.txt")
    return keep


def _purge_workspace(root: Path, keep_names: set[str]) -> tuple[list[str], list[str]]:
    """Remove DELFIN artifacts under root while preserving keep_names and unknown files."""
    removed: list[str] = []
    skipped: list[str] = []
    for entry in sorted(root.iterdir(), key=lambda p: p.name):
        name = entry.name
        if name in keep_names:
            continue
        try:
            if entry.is_symlink():
                skipped.append(f"{name} (symlink)")
                continue
        except OSError as exc:  # noqa: BLE001
            logger.warning("Failed to inspect %s: %s", entry, exc)
            skipped.append(name)
            continue

        if not _is_workspace_artifact(entry):
            skipped.append(name)
            continue

        try:
            if entry.is_dir():
                shutil.rmtree(entry)
            else:
                entry.unlink()
            removed.append(name)
        except FileNotFoundError:
            continue
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to remove %s: %s", entry, exc)
            skipped.append(name)
    return removed, skipped


def _is_delfin_workspace(root: Path) -> bool:
    """Check if the directory looks like a DELFIN workspace.

    Returns True if we find typical DELFIN artifacts (OCCUPIER folders,
    .inp/.out files, delfin_run.log, etc.), indicating this is a workspace
    where a DELFIN run has been executed.
    """
    # Check for typical DELFIN artifacts
    delfin_indicators = [
        "delfin_run.log",           # Main run log
        "initial_OCCUPIER",         # OCCUPIER folders
        "red_step_1_OCCUPIER",
        "red_step_2_OCCUPIER",
        "red_step_3_OCCUPIER",
        "ox_step_1_OCCUPIER",
        "ox_step_2_OCCUPIER",
        "ox_step_3_OCCUPIER",
        "XTB_OPT",                  # XTB folders
        "XTB_GOAT",
        "XTB_SOLVATOR",
        "CREST",
    ]

    # Check if any indicator exists
    for indicator in delfin_indicators:
        if (root / indicator).exists():
            return True

    # Check for .inp/.out files (typical ORCA outputs)
    has_inp = any(root.glob("*.inp"))
    has_out = any(root.glob("*.out"))

    if has_inp and has_out:
        return True

    # Check for ORCA temporary/output files (e.g., *.tmp, *.gbw_hs, *.scfgrad.inp)
    # These indicate an ORCA calculation has been run in this directory
    orca_patterns = ["*.tmp", "*.gbw_hs", "*.scfgrad.inp", "*.hostnames"]
    for pattern in orca_patterns:
        if any(root.glob(pattern)):
            return True

    return False


def _confirm_purge(root: Path) -> bool:
    """Confirm purge operation with safety checks."""
    # Safety check: Is this a DELFIN workspace?
    if not _is_delfin_workspace(root):
        print("⚠️  WARNING: This directory does NOT appear to be a DELFIN workspace!")
        print("   No typical DELFIN artifacts found (OCCUPIER folders, .inp/.out files, etc.)")
        print(f"   Current directory: {root.absolute()}")
        print()
        try:
            confirm = input("Are you ABSOLUTELY SURE you want to purge this directory? [yes/NO]: ")
        except EOFError:
            return False
        if confirm.strip().lower() != "yes":
            return False
        print()

    # Standard confirmation
    prompt = (
        "This will delete DELFIN artifacts (OCCUPIER folders, *.inp/*.out, logs) while keeping "
        "CONTROL.txt and the referenced input file. Continue? [y/N]: "
    )
    try:
        reply = input(prompt)
    except EOFError:
        return False
    if reply is None:
        return False
    return reply.strip().lower() in {"y", "yes"}


def _sync_occupier_spin_metadata(ctx: PipelineContext, config: dict) -> None:
    """Load preferred OCCUPIER multiplicity/BS for reporting."""
    occ_folder = ctx.control_file_path.parent / "initial_OCCUPIER"
    mult, bs = extract_preferred_spin(occ_folder)
    if mult is None:
        return
    ctx.extra['multiplicity_0'] = mult
    config['multiplicity_0'] = mult
    label = f"{mult} (BS {bs})" if bs else str(mult)
    ctx.extra['multiplicity_label'] = label
    config['_multiplicity_display'] = label




def main(argv: list[str] | None = None) -> int:
    configure_logging()
    arg_list = list(argv if argv is not None else sys.argv[1:])
    if arg_list and arg_list[0] == "cleanup":
        return _run_cleanup_subcommand(arg_list[1:])
    if arg_list and arg_list[0] == "stop":
        return _run_stop_subcommand(arg_list[1:])
    if arg_list and arg_list[0] == "co2":
        return _run_co2_subcommand(arg_list[1:])
    # ---- Parse flags first; --help/--version handled by argparse automatically ----
    parser = _build_parser()
    args, unknown = parser.parse_known_args(arg_list)
    if unknown:
        logger.error("Unknown argument(s): %s", " ".join(unknown))
        return 2

    try:
        occupier_overrides = _parse_occupier_overrides(getattr(args, "occupier_override", []) or [])
    except ValueError as exc:
        logger.error(str(exc))
        return 2

    RECALC_MODE = bool(getattr(args, "recalc", False))
    report_mode = getattr(args, "report", None)
    REPORT_TEXT = report_mode == "text"
    REPORT_DOCX = report_mode == "docx"
    workspace_arg = getattr(args, "workspace", ".") or "."
    workspace_hint = resolve_path(workspace_arg)
    if getattr(args, "control", "CONTROL.txt") == "CONTROL.txt":
        control_file_path = resolve_path(workspace_hint / "CONTROL.txt")
    else:
        control_file_path = resolve_path(args.control)
    workspace_root = control_file_path.parent
    if getattr(args, "json_output", None):
        candidate = Path(args.json_output)
        json_output_path = (workspace_root / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
    else:
        json_output_path = None
    force_rerun_outputs: set[Path] = {
        _override_output_path(name, workspace_root) for name in occupier_overrides
    } if RECALC_MODE else set()

    if occupier_overrides and not RECALC_MODE:
        logger.error("--occupier-override requires --recalc.")
        return 2

    os.environ["DELFIN_RECALC"] = "1" if RECALC_MODE else "0"

    if RECALC_MODE:
        # IMPORTANT: override the global bindings so all call sites use the wrappers
        global run_orca, XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR

        wrappers, reals = setup_recalc_mode(force_outputs=force_rerun_outputs)

        # Swap in wrappers in THIS module
        run_orca = wrappers['run_orca']
        XTB = wrappers['XTB']
        XTB_GOAT = wrappers['XTB_GOAT']
        run_crest_workflow = wrappers['run_crest_workflow']
        XTB_SOLVATOR = wrappers['XTB_SOLVATOR']

        # Patch other modules that captured their own references
        patch_modules_for_recalc(wrappers)


    # Only define template and exit
    if args.define:
        create_control_file(filename=str(control_file_path),
                            input_file=args.define,
                            overwrite=args.overwrite)
        return 0

    if getattr(args, "purge", False):
        keep = _safe_keep_set(control_file_path)
        # Ensure CONTROL itself is preserved even if named differently
        keep.add(control_file_path.name)

        if not _confirm_purge(workspace_root):
            print("Purge aborted.")
            return 0

        removed, skipped = _purge_workspace(workspace_root, keep)
        print(
            f"Workspace purged (removed {len(removed)} item(s); kept: {', '.join(sorted(keep))})."
        )
        if skipped:
            print("Skipped entries (not recognized as DELFIN artifacts):")
            for name in skipped:
                print(f"  {name}")
        return 0

    # Only cleanup and exit
    if args.cleanup:
        cleanup_all(str(get_runtime_dir()))
        print("Cleanup done.")
        return 0

    # JSON-only mode: build DELFIN_Data.json and exit
    if args.json and not any([report_mode, args.imag, args.recalc, args.define, args.purge, args.cleanup, args.afp]):
        try:
            output = save_esd_data_json(workspace_root, json_output_path)
            print(f"DELFIN data JSON written to: {output}")
            return 0
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to generate DELFIN_Data.json: %s", exc, exc_info=True)
            return 1

    # AFP-only mode: generate AFP plot and exit
    if args.afp and not any([report_mode, args.imag, args.recalc, args.define, args.purge, args.cleanup]):
        from .afp_plot import generate_afp_report
        try:
            output_png = generate_afp_report(workspace_root, fwhm=args.afp_fwhm)
            if output_png:
                print(f"AFP spectrum plot saved to: {output_png}")
                return 0
            else:
                logger.error("Failed to generate AFP plot - check if S0.out, S1.out, or T1.out exist")
                return 1
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to generate AFP plot: %s", exc, exc_info=True)
            return 1

    if occupier_overrides and (report_mode or getattr(args, "imag", False)):
        logger.error("--occupier-override is only supported for normal --recalc runs (not with --report/--imag).")
        return 2

    # Handle --report mode: recompute potentials from existing outputs
    if REPORT_TEXT:
        from .cli_report import run_report_mode

        # Read CONTROL.txt
        try:
            config = read_control_file(str(control_file_path))
        except FileNotFoundError:
            logger.error("CONTROL file not found: %s", control_file_path)
            logger.error("Hint: run `delfin %s --define` or pass `--control /path/to/CONTROL.txt`.", workspace_arg)
            return 2
        except ValueError as exc:
            logger.error("Invalid CONTROL configuration: %s", exc)
            return 2

        return run_report_mode(config)
    if REPORT_DOCX:
        from .cli_report_docx import run_docx_report_mode

        try:
            config = read_control_file(str(control_file_path))
        except FileNotFoundError:
            logger.error("CONTROL file not found: %s", control_file_path)
            logger.error("Hint: run `delfin %s --define` or pass `--control /path/to/CONTROL.txt`.", workspace_arg)
            return 2
        except ValueError as exc:
            logger.error("Invalid CONTROL configuration: %s", exc)
            return 2

        return run_docx_report_mode(
            workspace_root,
            config=config,
            afp_fwhm=args.afp_fwhm,
            json_output_path=json_output_path,
        )

    # Handle --imag mode: run IMAG elimination then report
    if getattr(args, "imag", False):
        from .cli_imag import run_imag_mode

        # Read CONTROL.txt
        try:
            config = read_control_file(str(control_file_path))
        except FileNotFoundError:
            logger.error("CONTROL file not found: %s", control_file_path)
            logger.error("Hint: run `delfin %s --define` or pass `--control /path/to/CONTROL.txt`.", workspace_arg)
            return 2
        except ValueError as exc:
            logger.error("Invalid CONTROL configuration: %s", exc)
            return 2

        return run_imag_mode(config, control_file_path)

    run_log_path = control_file_path.parent / "delfin_run.log"
    if "DELFIN_GLOBAL_LOG" not in os.environ:
        os.environ["DELFIN_GLOBAL_LOG"] = str(run_log_path)
    add_file_handler(os.environ["DELFIN_GLOBAL_LOG"])

    # Log DELFIN version and git commit at the very start
    from delfin import __version__
    git_commit = get_git_commit_info() or "unknown"
    logger.info("=" * 80)
    logger.info("DELFIN v%s - Git commit: %s", __version__, git_commit)
    logger.info("=" * 80)
    logger.info("Global run log attached at %s", os.environ["DELFIN_GLOBAL_LOG"])


    # --------------------- From here: normal pipeline run with banner --------------------
    print_delfin_banner()

    # ---- Friendly checks for missing CONTROL.txt / input file ----
    # Read CONTROL.txt once and derive all settings from it
    try:
        config = read_control_file(str(control_file_path))
    except FileNotFoundError:
        logger.error("CONTROL file not found: %s", control_file_path)
        logger.error("Hint: run `delfin %s --define` or pass `--control /path/to/CONTROL.txt`.", workspace_arg)
        return 2
    except ValueError as exc:
        logger.error("Invalid CONTROL configuration: %s", exc)
        return 2

    # Auto-configure cluster resources if not explicitly set
    config = auto_configure_resources(config)

    # Initialize global job manager with configuration
    global_mgr = get_global_manager()
    global_mgr.initialize(config)
    # Double-check signal handlers are active so Ctrl+C triggers cleanup
    global_mgr.ensure_signal_handlers()
    try:
        current_sigint = signal.getsignal(signal.SIGINT)
        if current_sigint != global_mgr._handle_sigint:
            logger.warning(
                "SIGINT handler mismatch (%s) - reinstalling DELFIN handler",
                current_sigint,
            )
            global_mgr._install_signal_handler()
            current_sigint = signal.getsignal(signal.SIGINT)
        logger.debug("Active SIGINT handler: %s", current_sigint)
    except Exception as exc:  # noqa: BLE001
        logger.debug("Could not verify SIGINT handler: %s", exc, exc_info=True)
    logger.info("Global job manager initialized")

    def _finalize(exit_code: int) -> int:
        """Shutdown global resources and perform optional cleanup before exiting."""
        try:
            global_mgr.shutdown()
            logger.info("Global job manager shutdown complete")
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Global manager shutdown raised: {exc}")

        if not args.no_cleanup:
            cleanup_all(str(get_runtime_dir()))
            cleanup_orca_scratch_dir()
        return exit_code

    if RECALC_MODE and occupier_overrides:
        if not _apply_occupier_overrides(occupier_overrides, workspace_root, config, force_rerun_outputs):
            return _finalize(1)

    # Store force_rerun_outputs in config so parallel_occupier.py can access it
    if RECALC_MODE and force_rerun_outputs:
        config["_recalc_force_outputs"] = force_rerun_outputs

    # Populate optional flags with safe defaults so reduced CONTROL files remain usable
    try:
        default_config = {
            'XTB_OPT': 'no',
            'XTB_GOAT': 'no',
            'CREST': 'no',
            'XTB_SOLVATOR': 'no',
            'calc_initial': 'yes',
            'oxidation_steps': '',
            'reduction_steps': '',
            'parallel_workflows': 'yes',
            'pal_jobs': None,
            'additions_TDDFT': '',
            'DONTO': 'FALSE',
            'DOSOC': 'TRUE',
            'FOLLOWIROOT': 'TRUE',
            'IROOT': '1',
            'NROOTS': '15',
            'ESD_nroots': 15,
            'ESD_maxdim': None,
            'TDDFT_maxiter': None,
            'ESD_TDDFT_maxiter': None,
            'TDA': 'FALSE',
            'NACME': 'TRUE',
            'ETF': 'TRUE',
            'implicit_solvation_model': 'CPCM',
            'maxcore': 3800,
            'maxiter': 125,
            'maxiter_occupier': 125,
            'mcore_E00': 10000,
            'multiplicity_0': None,
            'multiplicity_ox1': None,
            'multiplicity_ox2': None,
            'multiplicity_ox3': None,
            'multiplicity_red1': None,
            'multiplicity_red2': None,
            'multiplicity_red3': None,
            'out_files': None,
            'inp_files': None,
        }
        for key, value in default_config.items():
            config.setdefault(key, value)
    
        # Validate required files
        normalized_input = _normalize_input_file(config, control_file_path)
        success, error_code, _ = validate_required_files(config, control_file_path)
        input_file = normalized_input
        if not success:
            return _finalize(error_code)
    
        E_ref = get_E_ref(config) 
    
        NAME = (config.get('NAME') or '').strip()
    
        # Canonical file names used throughout workflows
        file_bundle = FileBundle(
            xyz_initial="initial.xyz",
            xyz_red1="red_step_1.xyz",
            xyz_red2="red_step_2.xyz",
            xyz_red3="red_step_3.xyz",
            xyz_ox1="ox_step_1.xyz",
            xyz_ox2="ox_step_2.xyz",
            xyz_ox3="ox_step_3.xyz",
            output_initial="initial.inp",
            output_absorption="absorption_td.inp",
            output_t1="t1_state_opt.inp",
            output_s1="s1_state_opt.inp",
            output_emission="emission_td.inp",
            output_ox1="ox_step_1.inp",
            output_ox2="ox_step_2.inp",
            output_ox3="ox_step_3.inp",
            output_red1="red_step_1.inp",
            output_red2="red_step_2.inp",
            output_red3="red_step_3.inp",
        )
    
        try:
            charge = int(str(config.get('charge', 0)).strip())
        except ValueError:
            logger.error("Invalid 'charge' in CONTROL.txt; falling back to 0.")
            charge = 0
        try:
            PAL = int(str(config.get('PAL', 6)).strip())
        except ValueError:
            logger.error("Invalid 'PAL' in CONTROL.txt; falling back to 6.")
            PAL = 6
        try:
            number_explicit_solv_molecules = int(str(config.get('number_explicit_solv_molecules', 0)).strip())
        except ValueError:
            logger.error("Invalid 'number_explicit_solv_molecules'; falling back to 0.")
            number_explicit_solv_molecules = 0
    
        solvent = (config.get('solvent') or '').strip()
        start_time = time.time()
    
        print(f"used Method: {config.get('method', 'UNDEFINED')}\n")
    
        metals = search_transition_metals(input_file)
        if metals:
            logger.info("Found transition metals:")
            for metal in metals:
                logger.info(metal)
        else:
            logger.info("No transition metals found in the file.")
    
        main_basisset, metal_basisset = set_main_basisset(metals, config)
    
        D45_SET = {
            'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
            'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'
        }
        use_rel = any(m in D45_SET for m in metals)
        if not use_rel:
            if str(config.get('relativity', '')).lower() != 'none':
                logger.info("3d-only system detected → relativity=none (ZORA/X2C/DKH is deactivated).")
            config['relativity'] = 'none' 
    
    
        total_electrons_txt, multiplicity_guess = calculate_total_electrons_txt(str(control_file_path))
        try:
            total_electrons_txt = int(total_electrons_txt)
        except (TypeError, ValueError):
            logger.error("Could not parse total electrons from CONTROL.txt; assuming 0.")
            total_electrons_txt = 0
    
        total_electrons = total_electrons_txt - charge
        is_even = (total_electrons % 2 == 0)
    
        try:
            cfg_mult_raw = config.get('multiplicity_global_opt') if config is not None else None
            cfg_mult = int(cfg_mult_raw) if cfg_mult_raw not in (None, "") else None
            if cfg_mult is not None and cfg_mult <= 0:
                cfg_mult = None
        except (TypeError, ValueError):
            cfg_mult = None
    
        try:
            ctl_mult_raw = multiplicity_guess.strip() if isinstance(multiplicity_guess, str) else multiplicity_guess
            ctl_mult = int(ctl_mult_raw) if ctl_mult_raw not in (None, "") else None
            if ctl_mult is not None and ctl_mult <= 0:
                ctl_mult = None
        except (TypeError, ValueError):
            ctl_mult = None
    
        multiplicity = cfg_mult if cfg_mult is not None else (ctl_mult if ctl_mult is not None else (1 if is_even else 2))
    
        # Ensure optional multiplicity slots share the detected ground-state multiplicity by default
        for mult_key in (
            'multiplicity_0',
            'multiplicity_ox1',
            'multiplicity_ox2',
            'multiplicity_ox3',
            'multiplicity_red1',
            'multiplicity_red2',
            'multiplicity_red3',
        ):
            if config.get(mult_key) in (None, ''):
                config[mult_key] = multiplicity
    
        metals_list = list(metals) if metals else []
    
        pipeline_ctx = PipelineContext(
            config=config,
            control_file_path=control_file_path,
            input_file=input_file,
            charge=charge,
            PAL=PAL,
            multiplicity=multiplicity,
            solvent=solvent,
            metals=metals_list,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            number_explicit_solv_molecules=number_explicit_solv_molecules,
            total_electrons_txt=total_electrons_txt,
            start_time=start_time,
            name=NAME,
            file_bundle=file_bundle,
        )
    
        raw_method = str(config.get('method', '')).strip()
        method_lower = raw_method.lower()
        esd_enabled = str(config.get('ESD_modul', 'no')).strip().lower() == 'yes'
    
        method_token: Optional[str]
        if method_lower in {'', 'none', 'esd'}:
            if not esd_enabled:
                logger.error(
                    "No method specified in CONTROL.txt and ESD_modul is not enabled. Supported methods: classic, manually, OCCUPIER"
                )
                return _finalize(2)
            if method_lower == 'esd':
                logger.info("CONTROL.txt method 'ESD' interpreted as ESD-only mode.")
            else:
                logger.info("No redox method requested; proceeding with ESD module only.")
            config['method'] = None
            method_token = None
        else:
            canonical_method, suggestion = interpret_method_alias(raw_method)

            if canonical_method not in {'classic', 'manually', 'OCCUPIER'}:
                if suggestion is not None:
                    logger.error(
                        "Unknown method '%s'. Did you mean '%s'?",
                        raw_method,
                        suggestion,
                    )
                else:
                    logger.error("Unknown method '%s'. Supported methods: classic, manually, OCCUPIER", raw_method)
                return _finalize(2)

            if canonical_method != raw_method and canonical_method.lower() != raw_method.lower():
                logger.warning("CONTROL.txt method '%s' interpreted as '%s'.", raw_method, canonical_method)

            config['method'] = canonical_method
            method_token = canonical_method
    
    
    
    
    
        # ------------------- OCCUPIER --------------------
        if method_token == "OCCUPIER":
            if not run_occuper_phase(pipeline_ctx):
                return _finalize(1)
            _sync_occupier_spin_metadata(pipeline_ctx, config)

            # Reapply overrides after OCCUPIER completes (in case OCCUPIER rewrote OCCUPIER.txt)
            pending_rewrites = config.get("_occ_overrides_pending_rewrite")
            if isinstance(pending_rewrites, dict) and pending_rewrites:
                logger.info("[recalc] Reapplying %d OCCUPIER override(s) after OCCUPIER phase", len(pending_rewrites))
                for folder_name, preferred_index in pending_rewrites.items():
                    folder_path = workspace_root / folder_name
                    occ_path = folder_path / "OCCUPIER.txt"
                    if occ_path.exists():
                        try:
                            content = occ_path.read_text(encoding="utf-8", errors="ignore")
                            new_content, count = _PREFERRED_INDEX_RE.subn(rf"\g<1>{preferred_index}", content, count=1)
                            if count > 0 and new_content != content:
                                occ_path.write_text(new_content, encoding="utf-8")
                                logger.info("[recalc] Re-set Preferred Index in %s to %d (post-OCCUPIER)", folder_name, preferred_index)
                        except Exception as exc:  # noqa: BLE001
                            logger.warning("[recalc] Could not rewrite %s after OCCUPIER: %s", occ_path, exc)
                config.pop("_occ_overrides_pending_rewrite", None)
    
        # ------------------- classic --------------------
        if method_token == "classic":
            run_classic_phase(pipeline_ctx)
    
    
        # ------------------- manually --------------------
        if method_token == "manually":
            run_manual_phase(pipeline_ctx)
    
        # ------------------- ESD (Excited State Dynamics) --------------------
        if esd_enabled:
            if method_token in (None, "classic"):
                # Check if ESD was already added to classic phase scheduler
                if pipeline_ctx.extra.get('esd_added_to_classic'):
                    logger.info("ESD jobs already executed in parallel with classic phase (skipping separate ESD phase)")
                else:
                    logger.info("Running ESD module...")
                    if not run_esd_phase(pipeline_ctx):
                        logger.warning("ESD module encountered issues, continuing...")
            else:
                logger.info(
                    "ESD module enabled but skipped for method '%s' (only classic/ESD-only supported).",
                    method_token,
                )
    
        # Finalize redox and emission summary
        if method_token == "OCCUPIER":
            mul0 = pipeline_ctx.extra.get('multiplicity_0')
            if mul0 is not None:
                try:
                    pipeline_ctx.multiplicity = int(mul0)
                except Exception:
                    pipeline_ctx.multiplicity = mul0  # fallback to raw value
            if '_multiplicity_display' not in config:
                label = pipeline_ctx.extra.get('multiplicity_label')
                if label:
                    config['_multiplicity_display'] = label
        elif method_token == "manually":
            try:
                pipeline_ctx.multiplicity = int(config.get('multiplicity_0', pipeline_ctx.multiplicity))
            except Exception:
                pass
        elif method_token == "classic":
            try:
                total_electrons_txt, mult_guess = calculate_total_electrons_txt(str(control_file_path))
                total_electrons_txt = int(total_electrons_txt)
                total_electrons = total_electrons_txt - pipeline_ctx.charge
                pipeline_ctx.multiplicity = 1 if total_electrons % 2 == 0 else 2
            except Exception:
                pass
    
        summary = compute_summary(pipeline_ctx, E_ref)
    
        charge = pipeline_ctx.charge
        solvent = pipeline_ctx.solvent
        metals_list = pipeline_ctx.metals
        main_basis = pipeline_ctx.main_basisset
        metal_basis = pipeline_ctx.metal_basisset
    
        generate_summary_report(
            charge,
            summary.multiplicity,
            solvent,
            summary.E_ox,
            summary.E_ox_2,
            summary.E_ox_3,
            summary.E_red,
            summary.E_red_2,
            summary.E_red_3,
            summary.E_00_t1,
            summary.E_00_s1,
            metals_list,
            metal_basis,
            NAME,
            main_basis,
            config,
            summary.duration,
            E_ref,
            summary.esd_summary,
            output_dir=control_file_path.parent,
        )
    
        if summary.esd_summary and summary.esd_summary.has_data:
            try:
                esd_report_path = control_file_path.parent / "ESD.txt"
                generate_esd_report(summary.esd_summary, esd_report_path)
                logger.info("ESD report written to %s", esd_report_path)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Failed to generate ESD.txt: %s", exc, exc_info=True)

        # Always write consolidated DELFIN JSON at the end of a run
        try:
            output_json = save_esd_data_json(workspace_root, json_output_path)
            logger.info("DELFIN data JSON written to %s", output_json)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to generate DELFIN_Data.json: %s", exc, exc_info=True)

        # Generate DOCX report at the end of the run
        try:
            from delfin.cli_report_docx import run_docx_report_mode
            afp_fwhm = float(config.get('afp_fwhm', 100.0))
            logger.info("Generating DELFIN.docx report...")
            report_status = run_docx_report_mode(workspace_root, config, afp_fwhm, json_output_path)
            if report_status == 0:
                logger.info("DELFIN.docx report generated successfully")
            else:
                logger.warning("DELFIN.docx report generation returned status %d", report_status)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to generate DELFIN.docx: %s", exc, exc_info=True)

        return _finalize(0)
    except KeyboardInterrupt:
        logger.warning("Execution interrupted by user; shutting down...")
        return _finalize(130)
    except Exception:
        logger.exception("Unhandled error during DELFIN run", exc_info=True)
        return _finalize(1)


if __name__ == "__main__":
    sys.exit(main())
