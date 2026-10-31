"""RUN1 physics finding R4: percent_buried_volume key holds a FRACTION.

Observed on the real cisplatin run (.gate/practice/run1_descriptors_
steric/, 2026-10-26): the wrapper delfin/analysis_tools/morfeus_wrapper
.py:122 stores morfeus's DEPRECATED percent_buried_volume attribute
under the key "percent_buried_volume".  Morfeus source (site-packages/
morfeus/buried_volume.py:619) returns self.fraction_buried_volume
UNCHANGED — a fraction in [0,1] — while its own print multiplies by 100
(line 547, "V_bur (%)").  A number labeled % that is actually a
fraction is off by a factor of 100: a physically wrong label, the kind
of defect wave 22 exists for.

Fix (in the .gate/px patch for the protected wrapper): the wrapper must
either multiply by 100 or stop calling the field "percent".  Chosen:
report fraction_buried_volume under a correctly-named key AND a
percent_buried_volume key that is really a percent (x100), so no
downstream consumer sees a number that says % and means fraction.

Red through the gate against the unfixed wrapper: the key exists, is in
[0,1] for cisplatin (58% buried -> 0.585), and contradicts the
documented unit.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def test_percent_buried_volume_key_is_really_percent(tmp_path):
    pytest = __import__("pytest")
    pytest.importorskip("morfeus")
    from delfin.analysis_tools import morfeus_wrapper

    xyz = tmp_path / "cisplatin.xyz"
    xyz.write_text(
        "9\n\n"
        "Pt   0.000000   0.000000   0.000000\n"
        "N    2.050000   0.000000   0.000000\n"
        "H    2.460000   0.891000   0.000000\n"
        "H    2.460000  -0.630000   0.771000\n"
        "N    0.000000   2.050000   0.000000\n"
        "H    0.000000   2.460000   0.891000\n"
        "H   -0.771000   2.460000   0.000000\n"
        "Cl  -2.300000   0.000000   0.000000\n"
        "Cl   0.000000  -2.300000   0.000000\n")
    data = morfeus_wrapper.buried_volume(xyz, metal_index=0, radius=3.5)
    percent = data.get("percent_buried_volume")
    fraction = data.get("fraction_buried_volume")
    assert isinstance(percent, (int, float)) and not isinstance(
        percent, bool), f"no percent key: {data}"
    # A fraction is 0.585 for cisplatin; a percent is 58.5.  The wrapper
    # must not label a fraction with a % name: percent >= 1 whenever
    # any atom is buried (fraction > 0), so this catches the off-by-100.
    assert percent > 1.0 or percent == 0.0, (
        f"'percent_buried_volume'={percent} is a FRACTION (0-1) mislabeled "
        f"as percent; morfeus buried_volume.py:619 returns "
        "fraction_buried_volume unchanged")
    if fraction is not None:
        assert abs(percent - 100.0 * fraction) < 1e-6
