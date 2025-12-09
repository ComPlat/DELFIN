#!/usr/bin/env python3
"""CLI tool for generating IR spectrum reports from ORCA frequency output files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import logging

from delfin.common.logging import get_logger
from delfin.ir_spectrum import parse_ir_spectrum
from delfin.reporting.ir_report import create_ir_spectrum_plot, generate_ir_report

logger = get_logger(__name__)


def main():
    """Main entry point for delfin_IR command."""
    parser = argparse.ArgumentParser(
        description="Generate IR spectrum Word report from ORCA frequency output file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate report from frequency calculation output
  delfin_IR molecule_freq.out

  # Generate report with custom output name
  delfin_IR S0.out -o IR_S0.docx

  # Generate report with custom broadening
  delfin_IR initial.out --fwhm 15 -o IR_Spectrum.docx

  # Verbose output
  delfin_IR calculation.out -v
"""
    )

    parser.add_argument(
        'output_file',
        type=str,
        help='Path to ORCA output file containing frequency calculation (e.g., S0.out, initial.out)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output Word document path (default: IR_<filename>.docx in same directory)'
    )

    parser.add_argument(
        '--fwhm',
        type=float,
        default=10.0,
        help='Full width at half maximum for Lorentzian broadening in cm⁻¹ (default: 10.0)'
    )

    parser.add_argument(
        '--min-intensity',
        type=float,
        default=1.0,
        help='Minimum intensity (km/mol) to include in summary table (default: 1.0)'
    )

    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Resolution of spectrum plot in DPI (default: 300)'
    )

    parser.add_argument(
        '--molecule-name',
        type=str,
        default=None,
        help='Optional molecule name for report title'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Validate input file
    input_file = Path(args.output_file).resolve()
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    if not input_file.suffix == '.out':
        logger.warning(f"Input file does not have .out extension: {input_file}")

    # Parse IR spectrum
    logger.info(f"Parsing IR spectrum from {input_file.name}")
    modes = parse_ir_spectrum(input_file)

    if not modes:
        logger.error("No IR spectrum data found in output file")
        logger.error("Make sure the file contains a frequency calculation with IR intensities")
        sys.exit(1)

    logger.info(f"Found {len(modes)} vibrational modes")

    # Count significant modes
    significant_count = sum(1 for m in modes if m.intensity_km_mol > args.min_intensity)
    logger.info(f"Found {significant_count} modes with intensity > {args.min_intensity} km/mol")

    # Determine output path
    if args.output:
        output_docx = Path(args.output).resolve()
    else:
        # Auto-generate filename: IR_<basename>.docx
        output_docx = input_file.parent / f"IR_{input_file.stem}.docx"

    logger.info(f"Output will be saved to {output_docx}")

    # Determine PNG output path (save in same directory as .docx)
    output_png = output_docx.parent / f"{output_docx.stem}.png"

    try:
        # Generate spectrum plot
        logger.info("Generating IR spectrum plot...")
        create_ir_spectrum_plot(
            modes=modes,
            output_png=output_png,
            fwhm=args.fwhm,
            dpi=args.dpi,
            title='IR Spectrum'
        )

        # Generate Word report
        logger.info("Generating Word report...")
        generate_ir_report(
            output_file=output_docx,
            modes=modes,
            spectrum_png=output_png,
            source_file=input_file,
            molecule_name=args.molecule_name
        )

        logger.info("Report generation complete!")
        print(f"\n✓ IR spectrum report saved to: {output_docx}")
        print(f"✓ IR spectrum plot saved to: {output_png}")
        print(f"  - Total vibrational modes: {len(modes)}")
        print(f"  - Significant modes (I > {args.min_intensity} km/mol): {significant_count}")

        # Find strongest mode
        if modes:
            strongest = max(modes, key=lambda m: m.intensity_km_mol)
            print(f"  - Strongest peak: {strongest.frequency_cm1:.1f} cm⁻¹ ({strongest.intensity_km_mol:.1f} km/mol)")

    except Exception as e:
        logger.error(f"Failed to generate report: {e}", exc_info=args.verbose)
        sys.exit(1)


if __name__ == "__main__":
    main()
