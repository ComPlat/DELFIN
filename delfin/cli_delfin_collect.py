#!/usr/bin/env python3
"""
CLI tool to collect all ESD data into JSON format.

Usage:
    delfin-esd-collect /path/to/project
    delfin-esd-collect /path/to/project --output my_data.json
"""

import argparse
import sys
from pathlib import Path

from delfin.common.logging import get_logger
from delfin.reporting.delfin_collector import save_esd_data_json

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Collect all DELFIN ESD calculation data into JSON format",
        epilog="""
Examples:
  # Collect data from current directory
  delfin-esd-collect .

  # Collect data from specific project
  delfin-esd-collect /home/user/calculations/molecule1

  # Specify output JSON file
  delfin-esd-collect /path/to/project --output results.json

  # Verbose output
  delfin-esd-collect /path/to/project -v
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "project_dir",
        type=str,
        help="Path to project directory containing ESD/, OCCUPIER/, initial.inp, etc."
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output JSON file path (default: <project_dir>/DELFIN_Data.json)"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        import logging
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    else:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Validate project directory
    project_dir = Path(args.project_dir)
    if not project_dir.exists():
        logger.error(f"Project directory does not exist: {project_dir}")
        sys.exit(1)

    if not project_dir.is_dir():
        logger.error(f"Not a directory: {project_dir}")
        sys.exit(1)

    # Check for ESD directory
    esd_dir = project_dir / "ESD"
    if not esd_dir.exists():
        logger.warning(f"ESD directory not found in {project_dir}")
        logger.warning("Continuing anyway - will collect available data")

    # Determine output path
    if args.output:
        output_json = Path(args.output)
    else:
        output_json = None  # Let delfin_collector.py determine the default

    logger.info(f"Collecting ESD data from: {project_dir}")
    if output_json:
        logger.info(f"Output JSON will be saved to: {output_json}")
    else:
        logger.info(f"Output JSON will be saved to: {project_dir / 'DELFIN_Data.json'}")

    try:
        output_path = save_esd_data_json(project_dir, output_json)
        logger.info("✓ Data collection complete!")
        print(f"\nJSON data saved to: {output_path}")
        print(f"\nYou can now generate a Word report with:")
        print(f"  delfin-esd-report {output_path}")

    except Exception as e:
        logger.error(f"Error collecting DELFIN data: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
