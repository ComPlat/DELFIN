"""Constants, templates, and style definitions for the DELFIN Dashboard."""

import ipywidgets as widgets

from delfin.config import set_control_value
from delfin.define import TEMPLATE as CONTROL_TEMPLATE

# ---------------------------------------------------------------------------
# Job time limits
# ---------------------------------------------------------------------------
JOB_TIME_LIMITS = {
    '12h': '12:00:00',
    '24h': '24:00:00',
    '36h': '36:00:00',
    '48h': '48:00:00',
    '60h': '60:00:00',
}

# ---------------------------------------------------------------------------
# Status colours (used in Job Status tab)
# ---------------------------------------------------------------------------
STATUS_COLORS = {
    'PENDING':   '#9E9E9E',
    'RUNNING':   '#4CAF50',
    'COMPLETED': '#2196F3',
    'FAILED':    '#f44336',
    'CANCELLED': '#FF9800',
    'TIMEOUT':   '#b71c1c',
}

# ---------------------------------------------------------------------------
# Common widget styling
# ---------------------------------------------------------------------------
COMMON_LAYOUT = widgets.Layout(width='500px')
COMMON_STYLE = {'description_width': 'initial'}

# ---------------------------------------------------------------------------
# Search suggestions for Calculations Browser
# ---------------------------------------------------------------------------
CALC_SEARCH_OPTIONS = [
    '***imaginary mode***',
    'ABSORPTION SPECTRUM',
    'ERROR',
    'FINAL SINGLE POINT ENERGY',
    'Final Gibbs free energy',
    'HURRAY',
    'JOB NUMBER ',
    'LOEWDIN POPULATION ANALYSIS',
    'LOEWDIN REDUCED ORBITAL CHARGES',
    'MULLIKEN POPULATION ANALYSIS',
    'MULLIKEN REDUCED ORBITAL CHARGES',
    'ORBITAL ENERGIES',
    'SCF CONVERGED AFTER',
    'TD-DFT EXCITED STATES',
    'Total Enthalpy',
    'WARNING',
]

# ---------------------------------------------------------------------------
# Job status table CSS
# ---------------------------------------------------------------------------
JOB_TABLE_CSS = """
<style>
    .job-table { font-family: monospace; font-size: 12px; border-collapse: collapse; width: 100%; }
    .job-table th, .job-table td { padding: 6px 10px; text-align: left; border-bottom: 1px solid #ddd; }
    .job-table th { background-color: #2196F3; color: white; }
    .job-table tr:hover { background-color: #f5f5f5; }
</style>
"""

# ---------------------------------------------------------------------------
# CONTROL.txt templates
#
# The template lives in delfin/define.py and nowhere else: it is what
# `delfin --define` writes, and delfin.config derives every runtime default
# from it. A second copy here would let the dashboard show one set of values
# while validating against another.
# ---------------------------------------------------------------------------
DEFAULT_CONTROL = CONTROL_TEMPLATE

def _only_global_optimisation(template: str) -> str:
    """The same file set up for a bare GOAT run: no redox steps, and the
    classic workflow so nothing waits on an OCCUPIER tree."""
    for key, value in (
        ('global_optimizer', 'GOAT'),
        ('calc_initial', 'no'),
        ('method', 'classic'),
    ):
        template = set_control_value(template, key, value)
    return template


ONLY_GOAT_TEMPLATE = _only_global_optimisation(DEFAULT_CONTROL)
