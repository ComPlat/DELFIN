"""F-C3 control: the skill's GFN2-xTB citation must carry the DOI of the
actual paper (10.1021/acs.jctc.8b01176, JCTC 2019, 15 (3), 1652-1671 —
verified via the ADS bibcode 2019JCTC...15.1652B, the publisher's
"Cite This" line and the Grimme-lab citation page). RED before the DOI
fix (f24c348d cites 8b01191, which resolves to nothing), GREEN after.
"""
from pathlib import Path
import re

SKILL = (Path(__file__).resolve().parents[1]
         / "delfin" / "agent" / "pack" / "skills"
         / "ensemble-crestandcenso.md")


def test_gfn2_xtb_doi_is_the_real_paper():
    text = SKILL.read_text(encoding="utf-8")
    assert "10.1021/acs.jctc.8b01176" in text, \
        "GFN2-xTB must cite DOI 10.1021/acs.jctc.8b01176"
    # The wrong DOI must not survive anywhere in the skill.
    assert "8b01191" not in text
    # Pages of the paper, as printed on the publisher's Cite line and
    # the Grimme-lab citation page: 15 (3), 1652-1671. The entry reads
    # "... 15 (3), 1652-1671 (2019), DOI 10.1021/acs.jctc.8b01176".
    m = re.search(r"15 \(3\), 1652-1671 \(2019\), DOI\s*\n?\s*"
                  r"10\.1021/acs\.jctc\.8b01176", text)
    assert m, "GFN2-xTB entry with its DOI missing"


def test_all_four_doi_literals_are_the_verified_ones():
    text = SKILL.read_text(encoding="utf-8")
    for doi in ("10.1021/acs.jpca.1c00971",   # CENSO, JPCA 2021
                "10.1039/C9CP06869D",          # CREST, PCCP 2020
                "10.1063/5.0040021",           # r2SCAN-3c, JCP 2021
                "10.1021/acs.jctc.8b01176"):   # GFN2-xTB, JCTC 2019
        assert doi in text, doi
