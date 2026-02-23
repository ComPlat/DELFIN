"""Archive tab: 1:1 mirror of the calculations browser on the archive/ folder."""

from dataclasses import replace

from . import tab_calculations_browser


def create_tab(ctx):
    """Create Archive tab as exact Calculations-browser mirror on ``archive/``."""
    ctx.archive_dir.mkdir(parents=True, exist_ok=True)
    archive_ctx = replace(ctx, calc_dir=ctx.archive_dir)
    return tab_calculations_browser.create_tab(archive_ctx)
