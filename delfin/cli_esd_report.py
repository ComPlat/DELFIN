"""CLI tool for generating UV-Vis spectrum reports from ESD output files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import logging

from delfin.common.logging import get_logger
from delfin.reporting.uv_vis_report import generate_uv_vis_word_report

logger = get_logger(__name__)


def main():
    """Main entry point for delfin_ESD command."""
    parser = argparse.ArgumentParser(
        description="Generate electronic spectrum Word report from ORCA ESD output file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate report from S0_TDDFT.out
  delfin_ESD ESD/S0_TDDFT.out

  # Generate report from S0.out with custom output name
  delfin_ESD ESD/S0.out -o my_spectrum.docx

  # Generate report with custom parameters
  delfin_ESD S0_TDDFT.out --min-fosc 0.01 --fwhm 30
"""
    )

    parser.add_argument(
        'output_file',
        type=str,
        help='Path to ORCA output file (e.g., ESD/S0_TDDFT.out or ESD/S0.out)'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output Word document path (default: Absorption_Spectrum_S0.docx in same directory)'
    )

    parser.add_argument(
        '--min-fosc',
        type=float,
        default=0.001,
        help='Minimum oscillator strength to include in table (default: 0.001)'
    )

    parser.add_argument(
        '--fwhm',
        type=float,
        default=20.0,
        help='Full width at half maximum for Gaussian broadening in nm (default: 20.0)'
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
    input_file = Path(args.output_file)
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    if not input_file.suffix == '.out':
        logger.warning(f"Input file does not have .out extension: {input_file}")

    # Determine output path (None = auto-generate with state name)
    output_docx = Path(args.output) if args.output else None

    # Generate report
    logger.info(f"Generating electronic spectrum report from {input_file}")
    if output_docx:
        logger.info(f"Output will be saved to {output_docx}")

    try:
        generate_uv_vis_word_report(
            input_file,
            output_docx,
            min_fosc=args.min_fosc,
            fwhm=args.fwhm
        )

        # Determine actual output filename if it was auto-generated
        if output_docx is None:
            state_name = input_file.stem.replace('_TDDFT', '')
            # Determine spectrum type
            if state_name == 'S0':
                spectrum_type = 'Absorption_Spectrum'
            elif state_name.startswith('S'):
                spectrum_type = 'Fluorescence_Spectrum'
            elif state_name.startswith('T'):
                spectrum_type = 'Phosphorescence_Spectrum'
            else:
                spectrum_type = 'Absorption_Spectrum'
            output_docx = input_file.parent / f"{spectrum_type}_{state_name}.docx"

        logger.info("Report generation complete!")
        print(f"\n✓ Report saved to: {output_docx}")

    except Exception as e:
        logger.error(f"Failed to generate report: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
