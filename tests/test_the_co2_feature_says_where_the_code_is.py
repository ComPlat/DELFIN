"""Three names for one feature (co2, co2_coordinator, co2_coordination)
and a catalog entry that told the orchestration story of a tool that
only dispatches: an operator answering the user's question rebuilt the
picture from eight grep/read rounds and asked for the entry to carry the
CONTROL keys, the entry paths and the files (2026-09-11).
"""

from delfin import api


def test_the_entry_names_the_keys_the_paths_and_the_files():
    info = api.explain_delfin_feature("co2")
    summary = info["summary"]
    for needle in ("co2_coordination", "co2_species_delta", "default off",
                   "delfin-co2-chain", "_run_co2_recalc_if_enabled",
                   "CO2_Coordinator6.main", "off every gate returns early"):
        assert needle in summary, needle
    see = " ".join(info["see_also"])
    for path in ("delfin/co2/CO2_Coordinator6.py", "delfin/co2/chain_setup.py",
                 "delfin/cli.py", "delfin/dashboard/tab_submit.py", "delfin/define.py"):
        assert path in see, path


def test_every_spelling_of_the_name_reaches_the_entry():
    for name in ("co2", "co2_coordinator", "co2_coordination", "CO2 coordinator"):
        assert api.explain_delfin_feature(name).get("name") == "co2", name
