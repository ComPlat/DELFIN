"""Generate a standalone ESD report."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from delfin.esd_results import ESDSummary


def _fmt_float(value, fmt: str = "{:.9f}", default: str = "n/a") -> str:
    if value is None:
        return default
    try:
        return fmt.format(value)
    except Exception:
        return default


def _fmt_rate(value) -> str:
    return _fmt_float(value, "{:.6e}")


def _fmt_percent(value) -> str:
    return _fmt_float(value, "{:.2f}")


def _emit_section(lines: list[str], title: str, content: Iterable[str]) -> None:
    lines.append("")
    lines.append(title)
    for entry in content:
        lines.append(entry)


def generate_esd_report(summary: ESDSummary, output_path: Path) -> None:
    """Write a detailed ESD report (ESD.txt) to disk."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = [
        "##############################################",
        "#           ESD SUMMARY REPORT               #",
        "##############################################",
    ]

    if summary.states:
        state_lines = []
        for state, record in sorted(summary.states.items()):
            state_lines.append(
                f"  {state}: {_fmt_float(record.fspe, fmt='{:.9f}')} Eh"
            )
        _emit_section(lines, "Final single point energies (Hartree):", state_lines)

    if summary.isc:
        isc_lines = []

        # Group ISC results by base transition (e.g., "S1>T1")
        from collections import defaultdict
        isc_grouped = defaultdict(list)
        for trans, record in summary.isc.items():
            # Extract base transition (remove Ms suffix if present)
            if "(Ms=" in trans:
                base_trans = trans.split("(Ms=")[0]
                isc_grouped[base_trans].append((trans, record))
            else:
                isc_grouped[trans].append((trans, record))

        # Sort and format ISC rates
        for base_trans in sorted(isc_grouped.keys()):
            records = isc_grouped[base_trans]

            # If multiple TROOTSSL values, show individual rates and sum
            if len(records) > 1:
                total_rate = sum(rec.rate for _, rec in records if rec.rate is not None)

                # Show individual sublevel rates
                for trans, record in sorted(records):
                    details = []
                    if record.temperature is not None:
                        details.append(f"T={_fmt_float(record.temperature, '{:.2f}')} K")
                    if record.delta_cm1 is not None:
                        details.append(f"Δ0-0={_fmt_float(record.delta_cm1, '{:.2f}')} cm^-1")
                    if record.soc is not None:
                        re_part, im_part = record.soc
                        soc_str = f"SOC={_fmt_float(re_part, '{:.6e}')}"
                        if im_part is not None:
                            soc_str += f"+i{_fmt_float(im_part, '{:.6e}')}"
                        details.append(soc_str)
                    if record.fc_percent is not None and record.ht_percent is not None:
                        details.append(
                            f"FC={_fmt_percent(record.fc_percent)}%, HT={_fmt_percent(record.ht_percent)}%"
                        )
                    elif record.fc_percent is not None:
                        details.append(f"FC={_fmt_percent(record.fc_percent)}%")
                    elif record.ht_percent is not None:
                        details.append(f"HT={_fmt_percent(record.ht_percent)}%")
                    detail_str = f" ({', '.join(details)})" if details else ""
                    isc_lines.append(
                        f"  {trans}: {_fmt_rate(record.rate)} s^-1{detail_str}"
                    )

                # Show total rate
                isc_lines.append(
                    f"  {base_trans} (total): {_fmt_rate(total_rate)} s^-1"
                )
            else:
                # Single TROOTSSL value, show as before
                trans, record = records[0]
                details = []
                if record.temperature is not None:
                    details.append(f"T={_fmt_float(record.temperature, '{:.2f}')} K")
                if record.delta_cm1 is not None:
                    details.append(f"Δ0-0={_fmt_float(record.delta_cm1, '{:.2f}')} cm^-1")
                if record.soc is not None:
                    re_part, im_part = record.soc
                    soc_str = f"SOC={_fmt_float(re_part, '{:.6e}')}"
                    if im_part is not None:
                        soc_str += f"+i{_fmt_float(im_part, '{:.6e}')}"
                    details.append(soc_str)
                if record.fc_percent is not None and record.ht_percent is not None:
                    details.append(
                        f"FC={_fmt_percent(record.fc_percent)}%, HT={_fmt_percent(record.ht_percent)}%"
                    )
                elif record.fc_percent is not None:
                    details.append(f"FC={_fmt_percent(record.fc_percent)}%")
                elif record.ht_percent is not None:
                    details.append(f"HT={_fmt_percent(record.ht_percent)}%")
                detail_str = f" ({', '.join(details)})" if details else ""
                isc_lines.append(
                    f"  {trans}: {_fmt_rate(record.rate)} s^-1{detail_str}"
                )

        _emit_section(lines, "ISC rate constants (s^-1):", isc_lines)

    if summary.ic:
        ic_lines = []
        for trans, record in sorted(summary.ic.items()):
            details = []
            if record.temperature is not None:
                details.append(f"T={_fmt_float(record.temperature, '{:.2f}')} K")
            if record.delta_cm1 is not None:
                details.append(f"Δ0-0={_fmt_float(record.delta_cm1, '{:.2f}')} cm^-1")
            detail_str = f" ({', '.join(details)})" if details else ""
            ic_lines.append(
                f"  {trans}: {_fmt_rate(record.rate)} s^-1{detail_str}"
            )
        _emit_section(lines, "IC rate constants:", ic_lines)

    if len(lines) == 3:
        lines.append("")
        lines.append("No ESD data available.")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
