# report.py
# Report generation functions have been moved to the reporting submodule

from .reporting.occupier_reports import generate_summary_report_OCCUPIER, generate_summary_report_OCCUPIER_safe
from .reporting.delfin_reports import generate_summary_report_DELFIN

# Re-export functions for backward compatibility
__all__ = [
    'generate_summary_report_OCCUPIER',
    'generate_summary_report_OCCUPIER_safe',
    'generate_summary_report_DELFIN'
]