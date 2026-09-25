"""Home tab: the calculations browser on the home directory.

The same browser as Calculations and Office -- browsing, preview, editing,
upload, download, drag and drop -- only rooted at the home directory, so a
file that is neither a calculation nor a document can be reached without
leaving the dashboard.

Two things differ from the folders the other clones sit on, and both are
handled here rather than in the browser:

* A home directory holds trees the browser has no business walking: a
  micromamba installation, caches, the archive. The recursive scan that
  decides whether the report button is enabled is bounded for this clone.
* Deleting is as available as it is in Calculations. That is what was asked
  for, and it is the user's own home; the tab is off by default so it is a
  deliberate choice rather than one more tab that happens to be there.
"""

from dataclasses import replace
from pathlib import Path

from . import tab_calculations_browser


def create_tab(ctx):
    """Create the Home tab as a Calculations-browser mirror on ``~``."""
    home_ctx = replace(ctx, calc_dir=Path.home())
    # Preserve the calculations root so items can be moved there.
    home_ctx.primary_calc_dir = ctx.calc_dir
    home_ctx.browser_scan_is_bounded = True
    return tab_calculations_browser.create_tab(home_ctx)
