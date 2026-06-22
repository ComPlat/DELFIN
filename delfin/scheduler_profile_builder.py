"""Build local scheduler stage profiles from archived DELFIN runs."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from delfin.scheduler_profiles import infer_stage_key, parse_molecule_text

_SECTION_RE = re.compile(r"^([A-Za-z0-9_]+):\s*$")
_AVG_TIME_RE = re.compile(r"^\s*Avg time:\s*([0-9.]+) min")
_AVG_CORES_RE = re.compile(r"^\s*Avg cores:\s*([0-9.]+)")
_JOBS_RE = re.compile(r"^\s*Jobs:\s*(\d+)")
_START_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) .* Starting ([A-Za-z0-9_]+) with (\d+) cores")
_ESD_DONE_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) .* State ([ST]\d) calculation completed")


def _parse_metrics(path: Path) -> List[Tuple[str, int, float, float]]:
    rows: List[Tuple[str, int, float, float]] = []
    current = None
    jobs = None
    avg_time = None
    avg_cores = None

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = _SECTION_RE.match(line)
        if match:
            if current and current != "SUMMARY" and jobs and avg_time is not None and avg_cores is not None:
                rows.append((current, jobs, avg_time * 60.0, avg_cores))
            current = match.group(1)
            jobs = None
            avg_time = None
            avg_cores = None
            continue
        if current is None:
            continue
        match = _JOBS_RE.match(line)
        if match:
            jobs = int(match.group(1))
            continue
        match = _AVG_TIME_RE.match(line)
        if match:
            avg_time = float(match.group(1))
            continue
        match = _AVG_CORES_RE.match(line)
        if match:
            avg_cores = float(match.group(1))

    if current and current != "SUMMARY" and jobs and avg_time is not None and avg_cores is not None:
        rows.append((current, jobs, avg_time * 60.0, avg_cores))

    return rows


def _latest_metrics_files(archive_root: Path) -> Dict[Path, Path]:
    latest: Dict[Path, Path] = {}
    for path in archive_root.rglob("performance_metrics_*.txt"):
        parent = path.parent
        current = latest.get(parent)
        if current is None or path.name > current.name:
            latest[parent] = path
    return latest


def _read_case_features(case_dir: Path) -> Dict[str, object]:
    for rel in ("input.txt", "start.txt", "initial.xyz"):
        candidate = case_dir / rel
        if not candidate.exists():
            continue
        parsed = parse_molecule_text(candidate.read_text(encoding="utf-8", errors="ignore"))
        if parsed is not None:
            parsed["path"] = str(candidate)
            return parsed
    return {"atom_count": None, "elements": [], "has_transition_metal": False, "transition_metals": []}


def _append_observation(store, stage: str, atom_count: int, duration_s: float, cores: float, weight: int) -> None:
    store[stage].append((atom_count, duration_s, cores, weight))


def _parse_log_stage_rows(log_path: Path) -> List[Tuple[str, float, float]]:
    rows: List[Tuple[str, float, float]] = []
    starts: Dict[str, Tuple[datetime, int]] = {}
    for line in log_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = _START_RE.match(line)
        if match:
            starts[match.group(2)] = (
                datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S"),
                int(match.group(3)),
            )
            continue

        match = _ESD_DONE_RE.match(line)
        if match:
            finished = datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S")
            job_id = f"esd_{match.group(2)}"
            started = starts.pop(job_id, None)
            if started is None:
                continue
            start_ts, cores = started
            duration_s = (finished - start_ts).total_seconds()
            if duration_s > 0:
                rows.append((infer_stage_key(job_id, f"ESD {match.group(2)} optimization"), duration_s, float(cores)))

    return rows


def _fit_stage_profile(rows: Iterable[Tuple[int, float, float, int]]) -> Dict[str, float]:
    data = list(rows)
    total_weight = sum(weight for *_prefix, weight in data)
    ref_atoms = sum(atom_count * weight for atom_count, _duration, _cores, weight in data) / total_weight
    avg_duration = sum(duration_s * weight for _atom_count, duration_s, _cores, weight in data) / total_weight
    observed_avg_cores = sum(cores * weight for _atom_count, _duration, cores, weight in data) / total_weight

    xs = [math.log(atom_count) for atom_count, duration_s, _cores, _weight in data if atom_count > 0 and duration_s > 0]
    ys = [math.log(duration_s) for atom_count, duration_s, _cores, _weight in data if atom_count > 0 and duration_s > 0]
    ws = [weight for atom_count, duration_s, _cores, weight in data if atom_count > 0 and duration_s > 0]

    atom_scaling_exp = 1.0
    if len(xs) >= 3:
        weight_sum = sum(ws)
        mean_x = sum(weight * x for weight, x in zip(ws, xs)) / weight_sum
        mean_y = sum(weight * y for weight, y in zip(ws, ys)) / weight_sum
        covariance = sum(weight * (x - mean_x) * (y - mean_y) for weight, x, y in zip(ws, xs, ys))
        variance = sum(weight * (x - mean_x) * (x - mean_x) for weight, x in zip(ws, xs))
        if variance > 1e-12:
            atom_scaling_exp = covariance / variance
    atom_scaling_exp = max(0.2, min(3.0, atom_scaling_exp))

    return {
        "duration_s": round(avg_duration, 1),
        "observed_avg_cores": round(observed_avg_cores, 1),
        "samples": int(total_weight),
        "reference_atoms": round(ref_atoms, 1),
        "atom_scaling_exp": round(atom_scaling_exp, 3),
    }


def build_stage_profiles(archive_root: Path) -> Dict[str, Dict[str, float]]:
    observations = defaultdict(list)

    for case_dir, metrics_path in _latest_metrics_files(archive_root).items():
        features = _read_case_features(case_dir)
        atom_count = features.get("atom_count")
        if not isinstance(atom_count, int) or atom_count <= 0:
            continue
        for stage, jobs, duration_s, avg_cores in _parse_metrics(metrics_path):
            _append_observation(observations, stage, atom_count, duration_s, avg_cores, jobs)

        log_path = case_dir / "delfin_run.log"
        if log_path.exists():
            for stage, duration_s, cores in _parse_log_stage_rows(log_path):
                _append_observation(observations, stage, atom_count, duration_s, cores, 1)

    profiles = {
        "_meta": {
            "version": 2,
            "generated_from": str(archive_root),
            "generated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        }
    }
    for stage, rows in sorted(observations.items()):
        profiles[stage] = _fit_stage_profile(rows)
    return profiles


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build local DELFIN scheduler stage profiles from archived runs.")
    parser.add_argument("archive_root", help="Archive root to scan for performance metrics and logs.")
    parser.add_argument("--output", default=str(Path(__file__).with_name("data") / "scheduler_stage_profiles.json"))
    args = parser.parse_args(argv)

    profiles = build_stage_profiles(Path(args.archive_root).expanduser())
    output_path = Path(args.output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(profiles, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
