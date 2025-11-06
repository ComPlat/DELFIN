# reporting/__init__.py
# Reporting submodule for DELFIN - provides access to all report generation functions

from .occupier_reports import generate_summary_report_OCCUPIER
from .delfin_reports import generate_summary_report_DELFIN
from .esd_report import generate_esd_report

__all__ = [
    'generate_summary_report_OCCUPIER',
    'generate_summary_report_DELFIN',
    'generate_esd_report'
]
