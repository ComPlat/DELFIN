"""Typing coordinates must not rebuild the 3D preview on every keystroke.

The Coordinates box sent every character to the kernel, and the handler behind
it tears down the 3Dmol viewer and builds a new one.  Typing therefore fought a
viewer that was being recreated between characters -- worst on the machines that
can least afford it.  The box now reports itself once typing stops.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.dashboard.context import DashboardContext
from delfin.dashboard.helpers import debounce_input


def test_the_coordinates_box_does_not_report_every_keystroke(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_orca_builder

    ctx = SimpleNamespace(
        backend=SimpleNamespace(),
        calc_dir=Path(tmp_path),
        run_js=lambda _script: None,
        add_init_js=lambda _script: None,
    )
    _tab, widgets_map = tab_orca_builder.create_tab(ctx)
    coords = widgets_map["orca_coords"]

    assert coords.continuous_update is False, (
        "every keystroke would rebuild the molecule preview"
    )
    assert "delfin-debounced" in coords._dom_classes, (
        "without the class nothing ever releases the value while typing"
    )


def test_the_debounce_script_releases_the_value_after_a_pause(monkeypatch):
    import IPython.display as ipd

    emitted = []
    monkeypatch.setattr(ipd, "display", lambda obj: emitted.append(obj))

    ctx = DashboardContext(calc_dir=Path("."))
    debounce_input(ctx, class_name="delfin-debounced", delay_ms=350)

    script = "\n".join(str(getattr(obj, "data", obj)) for obj in emitted)

    assert "delfin-debounced" in script
    assert "setTimeout" in script and "350" in script
    # a pause releases the value ...
    assert "new Event('change'" in script
    # ... and leaving the box must not sit and wait for that timer
    assert "'blur'" in script
