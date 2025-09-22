# reporting/__init__.py
# Reporting submodule for DELFIN - provides access to all report generation functions

from .occupier_reports import generate_summary_report_OCCUPIER, generate_summary_report_OCCUPIER_safe
from .delfin_reports import generate_summary_report_DELFIN

__all__ = [
    'generate_summary_report_OCCUPIER',
    'generate_summary_report_OCCUPIER_safe',
    'generate_summary_report_DELFIN'
]