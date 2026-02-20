"""Archiv tab: temporary 1:1 mirror of the calculations browser."""

from dataclasses import replace

from . import tab_calculations_browser


def create_tab(ctx):
    """Create Archiv tab as exact Calculations-browser mirror on ``archiv/``."""
    ctx.archiv_dir.mkdir(parents=True, exist_ok=True)
    archive_ctx = replace(ctx, calc_dir=ctx.archiv_dir)
    return tab_calculations_browser.create_tab(archive_ctx)
