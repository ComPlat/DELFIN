"""Office tab: 1:1 mirror of the calculations browser on the office/ folder.

Same browser, same spreadsheet and text editing -- only the root directory
differs, so documents live next to the calculations without mixing into them.
"""

from dataclasses import replace

from . import tab_calculations_browser


def create_tab(ctx):
    """Create the Office tab as an exact Calculations-browser mirror."""
    ctx.office_dir.mkdir(parents=True, exist_ok=True)
    office_ctx = replace(ctx, calc_dir=ctx.office_dir)
    # Preserve the calculations root so items can be moved back there.
    office_ctx.primary_calc_dir = ctx.calc_dir
    return tab_calculations_browser.create_tab(office_ctx)
