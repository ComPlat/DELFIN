"""Every topic a card uses is a known TOPICS entry, and every area of the
ORCA 6.1 manual has a fitting topic.

TOPICS is the vocabulary card authors file under. The manual's table of
contents (the folder and page names under the manual text tree) is the
ground truth for what areas exist: if a page of the manual has no topic,
cards about it get filed under "general" and the store loses the ability
to serve knowledge per area.
"""

from __future__ import annotations

import json
from pathlib import Path

from delfin import orca_knowledge as ok

CARDS_DIR = Path(ok.CARDS_DIR)

# Topic for each area of the ORCA 6.1 manual table of contents. Every
# entry names the manual pages it covers (folder names and page names).
MANUAL_AREA_TOPICS: dict[str, list[str]] = {
    # essential elements
    "counterpoise": ["contents/essentialelements/counterpoise.html"],
    "gcp": ["contents/essentialelements/counterpoise.html",
            "contents/modelchemistries/dispersioncorrections.html"],
    "fragments": ["contents/essentialelements/fragmentation.html"],
    "cpscf": ["contents/essentialelements/cpscf.html"],
    "field": ["contents/essentialelements/finEfield.html"],
    "frozen-core": ["contents/essentialelements/frozencore.html"],
    "integrals": ["contents/essentialelements/integralhandling.html",
                  "contents/essentialelements/shark.html"],
    "coordinates": ["contents/essentialelements/coordinates.html"],
    "numerical-gradients": ["contents/essentialelements/numericalgradients.html"],
    "cpcm": ["contents/essentialelements/solvationmodels.html"],
    "smd": ["contents/essentialelements/solvationmodels.html"],
    "cosmors": ["contents/essentialelements/solvationmodels.html"],
    # model chemistries
    "rpa": ["contents/modelchemistries/RPA.html"],
    "adft": ["contents/modelchemistries/analyticDFT.html"],
    "casdft": ["contents/modelchemistries/casdft.html"],
    "dcd-cas": ["contents/modelchemistries/dcdcas2.html"],
    "fci": ["contents/modelchemistries/energygradients_FCI.html"],
    "wavefunction-types": ["contents/modelchemistries/hftype.html"],
    "mrcc": ["contents/modelchemistries/mrcc.html"],
    "mreom": ["contents/modelchemistries/mreom.html"],
    "mm": ["contents/modelchemistries/mm.html"],
    "3c": ["contents/modelchemistries/3cmethods.html"],
    # spectroscopy and properties
    "population": ["contents/spectroscopyproperties/population.html"],
    "nbo": ["contents/spectroscopyproperties/nbo.html"],
    "nocv": ["contents/spectroscopyproperties/nocv.html"],
    "eda": ["contents/spectroscopyproperties/nocv.html"],
    "led": ["contents/spectroscopyproperties/led.txt",
            "contents/spectroscopyproperties/led.html"],
    "hfld": ["contents/spectroscopyproperties/hfld.html"],
    "adld": ["contents/spectroscopyproperties/adld.html"],
    "fod": ["contents/spectroscopyproperties/fod.html"],
    "vibrations": ["contents/spectroscopyproperties/vibrations.html"],
    "vpt2": ["contents/spectroscopyproperties/vpt2.html"],
    "asa": ["contents/spectroscopyproperties/asageneral.html",
            "contents/spectroscopyproperties/asa.html"],
    "electric-properties": ["contents/spectroscopyproperties/electric.html"],
    "deltascf": ["contents/spectroscopyproperties/deltascf.html"],
    "xas": ["contents/spectroscopyproperties/xes.html",
            "contents/spectroscopyproperties/rocis.html"],
    "xes": ["contents/spectroscopyproperties/xes.html"],
    "mcd": ["contents/spectroscopyproperties/mcd.html"],
    "moessbauer": ["contents/spectroscopyproperties/moessbauer.html"],
    "spin-orbit": ["contents/spectroscopyproperties/spinorbit.html"],
    "aniso": ["contents/spectroscopyproperties/single_aniso.html",
              "contents/spectroscopyproperties/poly_aniso.html"],
    "dichroism": ["contents/spectroscopyproperties/mcd.html"],
    "magprop": ["contents/spectroscopyproperties/postHF_magprop.html",
                "contents/spectroscopyproperties/qdpt_magnetic_properties.html"],
    "steom": ["contents/spectroscopyproperties/steom.html",
              "contents/spectroscopyproperties/dlpno-steom.html"],
    "rocis": ["contents/spectroscopyproperties/rocis.html"],
    "mcrpa": ["contents/spectroscopyproperties/mcrpa.html"],
    "ops": ["contents/spectroscopyproperties/ops.html"],
    "ihfsmr": ["contents/spectroscopyproperties/ihfsmr.html"],
    "casscf-resp": ["contents/spectroscopyproperties/CASSCFresp.html"],
    # structure and reactivity
    "conical-intersections": ["contents/structurereactivity/conicalintersections.html"],
    "mecp": ["contents/structurereactivity/mecp.html"],
    "dock": ["contents/structurereactivity/docker.html"],
    # multiscale
    "fmm": ["contents/multiscalesimulations/fmm.txt",
            "contents/multiscalesimulations/fmm.html"],
    "crystal-qmmm": ["contents/multiscalesimulations/qmmm-crystal.html"],
    # utilities, workflows, quickstart
    "plots": ["contents/utilitiesvisualization/plots.html"],
    "orca-2json": ["contents/utilitiesvisualization/orca_2json.html"],
    "property-file": ["contents/utilitiesvisualization/property_file.html"],
    "utilities": ["contents/utilitiesvisualization/utilities.html"],
    "opi": ["contents/workflowsautomatization/orcapythoninterface.html"],
    "installation": ["contents/quickstartguide/installation.html"],
    "gui": ["contents/quickstartguide/gui.html"],
    "troubleshooting": ["contents/quickstartguide/troubleshooting.html"],
    # molecular dynamics
    "aimd": ["contents/moleculardynamics/moldyn.html"],
}


def test_every_manual_area_has_a_topic() -> None:
    missing = sorted(t for t in MANUAL_AREA_TOPICS if t not in ok.TOPICS)
    assert not missing, f"manual areas without a topic in TOPICS: {missing}"


def test_every_card_topic_is_known() -> None:
    used: set[str] = set()
    for path in sorted(CARDS_DIR.glob("*.json")):
        data = json.loads(path.read_text())
        cards = data if isinstance(data, list) else data.get("cards", [])
        for card in cards:
            used.add(card.get("topic"))
    unknown = sorted(t for t in used if t and t not in ok.TOPICS)
    assert not unknown, f"cards use topics absent from TOPICS: {unknown}"
