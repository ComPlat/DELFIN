"""The mirror completion builds one enantiomer per CONFIGURATION, not one per system
(DELFIN_MIRROR_PER_ISOMER, 29.10.2026: an added legacy frame took the single mirror of ABOXIB / NELDEQ / YOQSUV)."""
import pytest

from delfin.manta import _mirror_enum as ME

pytest.importorskip("numpy")


def _xyz(rows, title="t"):
    out = [str(len(rows)), title]
    for s, x, y, z in rows:
        out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
    return "\n".join(out) + "\n"


# a chiral tetrahedron and a second, clearly different chiral arrangement (other atom positions)
_A = _xyz([("C", 0.0, 0.0, 0.0), ("N", 1.5, 0.0, 0.0), ("O", 0.0, 1.5, 0.0), ("F", 0.0, 0.0, 1.5),
           ("Cl", -1.0, -1.0, -1.0)])
_A2 = _xyz([("C", 0.0, 0.0, 0.0), ("N", 1.52, 0.0, 0.0), ("O", 0.0, 1.49, 0.0), ("F", 0.0, 0.0, 1.51),
            ("Cl", -1.0, -1.02, -1.0)])
_B = _xyz([("C", 0.0, 0.0, 0.0), ("N", 0.0, 1.5, 0.0), ("O", 2.5, 0.0, 0.0), ("F", 0.0, 0.3, -2.5),
           ("Cl", 1.0, 2.0, 1.0)])


@pytest.fixture
def env(monkeypatch):
    monkeypatch.setenv("DELFIN_MIRROR_ENUM", "1")
    monkeypatch.setenv("DELFIN_MIRROR_ONE_PER_SYSTEM", "1")
    monkeypatch.delenv("DELFIN_MIRROR_STEREO_GATE", raising=False)
    monkeypatch.delenv("DELFIN_MIRROR_QUALITY_GATE", raising=False)
    return monkeypatch


def test_root_strips_conformer_suffixes():
    assert ME._iso_root("I-N-ax-Δ-conf5_mirror") == "I-N-ax-Δ"
    assert ME._iso_root("C-N-ax-1-reembed2-conf3") == "C-N-ax-1"
    assert ME._iso_root("fac-cis-1-jte") == "fac-cis-1"
    assert ME._hand_swap("N1-P-ax-Λ") == "N1-P-ax-Δ"
    assert ME._hand_swap("axial-1") is None


def test_off_keeps_one_mirror_per_system(env):
    env.delenv("DELFIN_MIRROR_PER_ISOMER", raising=False)
    out = ME.expand_results([(_B, "I-N-ax-Δ-conf5"), (_A, "C-N-ax-1")])
    mirrors = [l for _x, l in out if l.endswith("_mirror")]
    assert mirrors == ["I-N-ax-Δ-conf5_mirror"]          # the first frame takes the only slot


def test_on_mirrors_every_configuration_once(env):
    env.setenv("DELFIN_MIRROR_PER_ISOMER", "1")
    res = [(_B, "I-N-ax-Δ-conf5"), (_A, "C-N-ax-1"), (_A2, "C-N-ax-1-conf2")]
    out = ME.expand_results(res)
    mirrors = sorted(l for _x, l in out if l.endswith("_mirror"))
    assert mirrors == ["C-N-ax-1_mirror", "I-N-ax-Δ-conf5_mirror"]
    assert out[:3] == res                                 # originals untouched, additive
    again = ME.expand_results(out)                        # re-entry appends nothing
    assert len(again) == len(out)


def test_on_skips_configuration_whose_partner_hand_exists(env):
    env.setenv("DELFIN_MIRROR_PER_ISOMER", "1")
    out = ME.expand_results([(_B, "N1-P-ax-Δ"), (_A, "N1-P-ax-Λ")])
    assert not [l for _x, l in out if l.endswith("_mirror")]


def test_on_skips_mirror_that_duplicates_a_present_frame(env):
    env.setenv("DELFIN_MIRROR_PER_ISOMER", "1")
    m = ME.mirror_frame(_A)
    out = ME.expand_results([(_A, "iso1"), (m, "iso2")])  # iso2 IS the enantiomer of iso1 already
    assert not [l for _x, l in out if l.endswith("_mirror")]
