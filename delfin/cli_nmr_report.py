#!/usr/bin/env python3
"""CLI tool for generating 1H NMR spectrum reports from ORCA NMR output files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from delfin.common.logging import get_logger
from delfin.nmr_spectrum import parse_nmr_orca
from delfin.reporting.nmr_report import create_nmr_report, print_assignment_table

logger = get_logger(__name__)


def main():
    """Main entry point for delfin_NMR command."""
    parser = argparse.ArgumentParser(
        description="Generate 1H NMR spectrum plot from ORCA NMR output file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage — generates PNG with spectrum + molecule
  delfin_NMR molecule_nmr.out

  # Custom output name
  delfin_NMR THBMT.out -o NMR_THBMT.png

  # Custom TMS reference and broadening
  delfin_NMR calc.out --tms-ref 31.4 --fwhm 0.03

  # Only show 0–10 ppm range, group tolerance 0.2 ppm
  delfin_NMR calc.out --ppm-min 0 --ppm-max 10 --equiv-tol 0.2

  # Print assignment table to terminal
  delfin_NMR calc.out --table
"""
    )

    parser.add_argument(
        "output_file",
        type=str,
        help="Path to ORCA .out file containing NMR calculation",
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output PNG path (default: NMR_<filename>.png in same directory)",
    )

    parser.add_argument(
        "--tms-ref",
        type=float,
        default=31.4,
        help="TMS reference shielding in ppm (default: 31.4 for TPSS/pcSseg-1)",
    )

    parser.add_argument(
        "--fwhm",
        type=float,
        default=0.02,
        help="Full width at half maximum for Lorentzian broadening in ppm (default: 0.02)",
    )

    parser.add_argument(
        "--equiv-tol",
        type=float,
        default=0.15,
        help="Tolerance in ppm to group equivalent H atoms (default: 0.15)",
    )

    parser.add_argument(
        "--ppm-min",
        type=float,
        default=-1.0,
        help="Lower ppm limit (default: -1.0)",
    )

    parser.add_argument(
        "--ppm-max",
        type=float,
        default=14.0,
        help="Upper ppm limit (default: 14.0)",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution of output image in DPI (default: 300)",
    )

    parser.add_argument(
        "--molecule-name",
        type=str,
        default=None,
        help="Optional molecule name for plot title",
    )

    parser.add_argument(
        "--table",
        action="store_true",
        help="Print assignment table to terminal",
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )

    args = parser.parse_args()

    if args.verbose:
        import logging
        logging.basicConfig(level=logging.DEBUG)

    # Resolve paths
    orca_out = Path(args.output_file)
    if not orca_out.exists():
        print(f"Error: File not found: {orca_out}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        out_png = Path(args.output)
    else:
        out_png = orca_out.parent / f"NMR_{orca_out.stem}.png"

    # Parse
    result = parse_nmr_orca(orca_out, tms_ref=args.tms_ref)

    if not result.h_shieldings:
        print("Error: No hydrogen shieldings found in the output file.", file=sys.stderr)
        sys.exit(1)

    print(f"Parsed {len(result.h_shieldings)} H shieldings from {orca_out.name}")

    # Assignment table
    if args.table:
        print()
        print(print_assignment_table(result, equiv_tol=args.equiv_tol))
        print()

    # Plot
    title = args.molecule_name or orca_out.stem.replace("_", " ")
    create_nmr_report(
        result,
        output_png=out_png,
        ppm_range=(args.ppm_min, args.ppm_max),
        fwhm=args.fwhm,
        equiv_tol=args.equiv_tol,
        dpi=args.dpi,
        title=f"Calculated 1H NMR: {title}",
    )

    print(f"NMR report saved to: {out_png}")


if __name__ == "__main__":
    main()
