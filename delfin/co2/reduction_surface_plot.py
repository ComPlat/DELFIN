#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reduction_surface_plot.py — the automatic energy profile plot of the
CO2 reduction surface chain.

Writes ``reduction_surface_profile.png`` and ``.pdf`` into the surface
workdir: ΔE in kcal/mol relative to the first station, along
M···CO2 -> M···CO2H -> M···CO2H2 -> M···CO + H2O, with the NEB-TS points
between the stations. Failed/missing stations appear as labelled gaps, not
as invented values. Two panels when xTB energies exist (xTB profile and the
ORCA/OCCUPIER profile), one otherwise.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HARTREE_TO_KCAL = 627.509474


def _station_points(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Profile points from the stations and the NEB-TS results."""
    stations = result.get("stations", []) or []
    ts_results = result.get("transition_states", []) or []

    points: List[Dict[str, Any]] = []
    ref_energy: Optional[float] = None
    for st in stations:
        energy = st.get("energy_eh")
        if energy is None:
            energy = st.get("energy_eh_xtb")
            source = "xtb"
        else:
            source = "orca"
        if energy is None:
            points.append({"label": st["label"], "station": st["station"],
                           "delta_kcal": None, "kind": "minimum",
                           "status": st.get("status")})
            continue
        if ref_energy is None:
            ref_energy = energy
        points.append({"label": st["label"], "station": st["station"],
                       "delta_kcal": (energy - ref_energy) * HARTREE_TO_KCAL,
                       "kind": "minimum", "status": st.get("status"),
                       "source": source})
    return points, ts_results


def write_profile_plot(result: Dict[str, Any], surf_dir: str):
    """Write the PNG + PDF profile. Returns (png_path, pdf_path)."""
    points, ts_results = _station_points(result)
    png = surf_dir.rstrip("/") + "/reduction_surface_profile.png"
    pdf = surf_dir.rstrip("/") + "/reduction_surface_profile.pdf"

    has_orca = any(p.get("source") == "orca" for p in points)
    has_xtb = any(p.get("source") == "xtb" and p.get("delta_kcal") is not None
                  for p in points)
    two_panels = has_orca and has_xtb
    panels = 2 if two_panels else 1
    fig, axes = plt.subplots(panels, 1, figsize=(8, 4 + 2 * panels),
                             squeeze=False)
    axes = axes.flatten()

    for ax, source in zip(axes, ("orca", "xtb") if two_panels else ("mixed",)):
        _draw_panel(ax, points, ts_results, source, two_panels)

    fig.tight_layout()
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def _draw_panel(ax, points, ts_results, source: str, two_panels: bool) -> None:
    title = {"orca": "ORCA (OCCUPIER) energies",
             "xtb": "GFN2-xTB preopt energies",
             "mixed": "CO2 reduction surface"}[source]
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(r"$\Delta E$ (kcal/mol)")

    # x positions: one slot per station, TS slots between them.
    labels: List[str] = []
    xs: List[float] = []
    ys: List[float] = []
    ts_x: List[float] = []
    ts_y: List[float] = []
    gap_labels: List[Tuple[float, float, str]] = []

    station_slots: List[int] = []
    slot = 0.0
    for p in points:
        station_slots.append(slot)
        labels.append(p["label"])
        slot += 2.0  # room for the TS point in between

    prev_point: Optional[Dict[str, Any]] = None
    for i, p in enumerate(points):
        x = station_slots[i]
        if p["delta_kcal"] is None:
            # Only a station whose RUN failed is a red gap; a dry station
            # (status ok, no energy computed) is simply unlabeled — no
            # invented value, no false alarm.
            if p.get("status") not in (None, "ok"):
                gap_labels.append((x, 0.0, f"{p['label']} failed"))
            prev_point = None
            continue
        y = p["delta_kcal"]
        xs.append(x)
        ys.append(y)
        # TS point between the previous and this station
        ts = None
        if i > 0:
            for t in ts_results:
                if t.get("step") == f"ts{points[i-1]['station']}_{p['station']}":
                    ts = t
        if ts is not None and ts.get("ts_energy_kcal_rel_reactant") is not None:
            # relative to the REACTANT of the pair; to plot the barrier on
            # the global reference we need the previous station's energy —
            # the caller already stored it relative to the reactant, so the
            # barrier height ADDS to the previous station's profile value.
            prev_y = ys[-2] if len(ys) >= 2 else None
            if prev_y is not None:
                ts_y.append(prev_y + ts["ts_energy_kcal_rel_reactant"])
                ts_x.append((station_slots[i - 1] + x) / 2.0)
        prev_point = p

    if xs:
        ax.plot(xs, ys, "-", color="#4477aa", zorder=1)
        ax.scatter(xs, ys, marker="o", color="#4477aa", zorder=2, label="minimum")
    if ts_x:
        ax.scatter(ts_x, ts_y, marker="^", color="#ee6677", zorder=2,
                   label="NEB-TS")
    for x, y, label in gap_labels:
        ax.text(x, y, label, ha="center", va="bottom", color="#cc3311",
                fontsize=8, rotation=0)
    ax.set_xticks(station_slots)
    ax.set_xticklabels(labels)
    ax.axhline(0.0, color="grey", linewidth=0.5)
    if xs or ts_x:
        ax.legend(loc="best", fontsize=8)
