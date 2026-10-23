"""The command line and the dashboard build the same MANTA manifold.

A landing that reaches only one of the two entry points is a landing the other
half of the users never gets, and nothing complains: both still emit a valid
manifold, just not the same one.  This test builds a handful of small systems
through ``delfin-manta`` (its real ``main``) and through the dashboard's MANTA
button (``structure_editor._run_smiles_build`` with the settings row at its
widget defaults, which is exactly what the button worker calls) -- each in a
fresh interpreter with a clean environment -- and asserts that the emitted
frames are byte-identical: same count, same order, same labels, same
coordinate lines.

The two sides get DIFFERENT ``PYTHONHASHSEED`` values on purpose.  The
dashboard's isolation subprocess pins ``PYTHONHASHSEED=0`` while a shell
running ``delfin-manta`` has whatever the user has (random by default), so a
build that depends on set / dict iteration order would differ between the two
entry points; the differing seeds make that visible here.

Run directly as ``python tests/test_cli_dashboard_parity.py {cli|dashboard}
SMILES OUT.json`` to produce one side (that is how the test calls it).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]

# Small and fast (seconds each at quality extreme), and covering the three
# shapes named in the parity requirement -- square planar with cis/trans, a
# chelate, an octahedral complex -- plus a hapto complex (eta2 x4 Ir, from the
# 6000 pool) and a dot SMILES, which MANTA builds as the whole string on both
# sides (the dashboard used to split it into parts).  The formula of every frame
# is pinned separately in test_every_frame_carries_the_whole_smiles.py; this test
# pins PARITY.
_SMILES = (
    "[Pt](Cl)(Cl)(N)N",
    "[Pt]1(Cl)(Cl)NCCN1",
    "Cl[Co+3](Cl)([NH3])([NH3])([NH3])[NH3]",
    "[Cl][Sn]([Cl])([Cl])[Ir]123456([C]7=[C]1CC[C]2=[C]3CC7)[C]1=[C]4CC[C]5=[C]6CC1",
    "[Pt](Cl)(Cl)(N)N.Cl",
)


def _side_cli(smiles: str, work: Path) -> dict:
    from delfin import cli_manta

    out = work / "cli_out"
    code = cli_manta.main([smiles, "-o", str(out), "-q"])
    manifest = json.loads((out / "manifest.json").read_text())
    frames = []
    for item in manifest["isomers"]:
        text = (out / item["file"]).read_text()
        frames.append([item["label"], cli_manta._atom_lines(text)])
    return {"exit": code, "frames": frames}


def _side_dashboard(smiles: str, work: Path) -> dict:
    from delfin import cli_manta
    from delfin.dashboard import structure_editor as se

    before = {k: v for k, v in os.environ.items() if k.startswith("DELFIN_")}
    result = se._run_smiles_build(smiles, **se._manta_button_kwargs())
    after = {k: v for k, v in os.environ.items() if k.startswith("DELFIN_")}
    # The construction env goes to the build subprocess only, never into the kernel.
    assert before == after, sorted(set(after.items()) ^ set(before.items()))
    if result.get("error"):
        return {"exit": 1, "frames": [], "error": result["error"]}
    frames = [[label, cli_manta._atom_lines(xyz)] for xyz, _n, label in result["isomers"]]
    return {"exit": 0, "frames": frames}


def _clean_env(hashseed: str) -> dict:
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("DELFIN_") and k != "PYTHONHASHSEED"}
    if hashseed is not None:
        env["PYTHONHASHSEED"] = hashseed
    env["PYTHONPATH"] = str(_ROOT) + (os.pathsep + env["PYTHONPATH"]
                                      if env.get("PYTHONPATH") else "")
    return env


def _run_side(side: str, smiles: str, tmp: Path, hashseed: str) -> dict:
    work = tmp / side
    work.mkdir(parents=True, exist_ok=True)
    out = work / "result.json"
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), side, smiles, str(out)],
        cwd=str(work), env=_clean_env(hashseed), capture_output=True, text=True,
        timeout=900)
    assert proc.returncode == 0, (side, proc.stdout[-2000:], proc.stderr[-4000:])
    data = json.loads(out.read_text())
    assert Path(data["delfin_file"]).resolve().is_relative_to(_ROOT), data["delfin_file"]
    return data


@pytest.mark.slow
@pytest.mark.parametrize("smiles", _SMILES)
def test_cli_and_dashboard_emit_the_same_manifold(smiles, tmp_path):
    cli = _run_side("cli", smiles, tmp_path, hashseed="12345")
    # As in a Voila kernel: no seed in the parent, the isolation child pins 0.
    dash = _run_side("dashboard", smiles, tmp_path, hashseed=None)

    assert cli["exit"] == 0, cli
    assert dash["exit"] == 0, dash
    assert cli["frames"], "CLI emitted nothing"

    cli_labels = [f[0] for f in cli["frames"]]
    dash_labels = [f[0] for f in dash["frames"]]
    assert dash_labels == cli_labels, (
        f"labels/order differ for {smiles}:\n  cli       {cli_labels}\n"
        f"  dashboard {dash_labels}")
    for i, (a, b) in enumerate(zip(cli["frames"], dash["frames"])):
        assert a[1] == b[1], f"frame {i} ({a[0]}) coordinates differ for {smiles}"


@pytest.mark.parametrize("config", ["champion", "builder", "default"])
@pytest.mark.parametrize("rank", [False, True])
def test_the_dashboard_env_is_the_cli_env(config, rank, monkeypatch):
    """Same switches for the same settings, including an environment override."""
    from delfin import cli_manta
    from delfin.dashboard import structure_editor as se

    monkeypatch.setenv("DELFIN_MIRROR_ENUM", "0")      # a user override of an extra setting
    cli = cli_manta.construction_env(config, rank=rank, method="gfn2", charge=3)
    dash = se._manta_best_env(3, construction=config, method="gfn2", rank=rank)
    assert dash == cli
    if config == "champion":
        assert dash["DELFIN_MIRROR_ENUM"] == "0"
        assert all(dash["DELFIN_FFFREE_" + f] == "1" for f in cli_manta._CHAMPION_FLAGS)


def test_cli_rank_charge_comes_from_the_smiles():
    """--rank without --charge ranks at the SMILES formal charge, as the dashboard does."""
    from delfin import cli_manta

    assert cli_manta._charge_for_opt("Cl[Co+3](Cl)([NH3])([NH3])([NH3])[NH3]", None) == 3
    assert cli_manta._charge_for_opt("Cl[Co+3](Cl)([NH3])([NH3])([NH3])[NH3]", 1) == 1


def test_the_hapto_retry_respects_an_explicit_fail_fast():
    """One rule for every entry point: retry a fail-fast answer with the
    approximation, unless hapto was forced or DELFIN_HAPTO_APPROX=0 says so."""
    from delfin.common.manta_build import hapto_retry_wanted

    err = "Hapto (eta) coordination detected (1 group(s), max eta~5)."
    assert hapto_retry_wanted(err, None, {}) is True
    assert hapto_retry_wanted(err, None, {"DELFIN_HAPTO_APPROX": "0"}) is False
    assert hapto_retry_wanted(err, None, {"DELFIN_HAPTO_APPROX": "off"}) is False
    assert hapto_retry_wanted(err, False, {}) is False
    assert hapto_retry_wanted("something else", None, {}) is False


def test_every_build_subprocess_gets_the_same_hash_seed(monkeypatch):
    """Forced, not defaulted: a shell or kernel seed must not reach the build."""
    from delfin.common import manta_build

    seen = {}

    class _Proc:
        returncode = 0
        args = ()

        def __init__(self, *a, env=None, **k):
            seen.update(env)

        def communicate(self, input=None, timeout=None):
            return "__DELFIN_RESULT__" + json.dumps({"r": [], "e": None}), ""

    monkeypatch.setenv("PYTHONHASHSEED", "4242")
    monkeypatch.setattr(manta_build.subprocess, "Popen", _Proc)
    manta_build.run_isomers_isolated("C", {}, env={"DELFIN_X": "1"})
    assert seen["PYTHONHASHSEED"] == manta_build.HASH_SEED == "0"
    assert seen["DELFIN_X"] == "1"
    assert "DELFIN_X" not in os.environ


def test_the_pipeline_zero_means_the_complete_manifold():
    """MANTA_MAX_ISOMERS=0 (and an absent key) is the complete manifold, as in the
    CLI and the dashboard -- not a silent cap of 100."""
    source = (_ROOT / "delfin" / "guppy_sampling.py").read_text()
    assert "else 100000" in source and "> 0 else 100\n" not in source
    assert "'num_confs': target_confs" not in source
    pipeline = (_ROOT / "delfin" / "workflows" / "pipeline.py").read_text()
    assert "_setting('MANTA_MAX_ISOMERS', 'GUPPY_MAX_ISOMERS', 0, int)" in pipeline


def _write_one_side(argv) -> int:
    side, smiles, out = argv[1], argv[2], Path(argv[3])
    import delfin
    runner = _side_cli if side == "cli" else _side_dashboard
    data = runner(smiles, out.parent)
    data["delfin_file"] = delfin.__file__
    out.write_text(json.dumps(data))
    return 0


if __name__ == "__main__":
    raise SystemExit(_write_one_side(sys.argv))


# --- the one deliberate exception: the CONTROL pipeline's gates ---------------

#: Construction switches a CONTROL run adds on top of delfin-manta (user decision
#: 2026-09-28: a torn frame costs a DFT chain; CLI and dashboard return the full
#: manifold).  COORD_INTEGRITY / CONF_COMPLETE = 0 restate the builder default.
_PIPELINE_ONLY = {
    "DELFIN_FFFREE_CLEAN_GATE": "1",
    "DELFIN_FFFREE_TOPOLOGY_GATE": "1",
    "DELFIN_FFFREE_PERMUTE_DEDUP": "1",
    "DELFIN_FFFREE_COORD_INTEGRITY": "0",
    "DELFIN_FFFREE_CONF_COMPLETE": "0",
}
#: Resource settings, not construction: how many UFF workers and how long.
_RESOURCE_KEYS = {"DELFIN_MAX_PROCESS_WORKERS", "DELFIN_UI_ISOLATE_TIMEOUT"}


def _pipeline_env(config, monkeypatch):
    from delfin.common import manta_settings

    for key in list(os.environ):
        if key.startswith("DELFIN_"):
            monkeypatch.delenv(key, raising=False)
    target = {}
    manta_settings.apply_construction_env(config, target)
    return {k: v for k, v in target.items() if k not in _RESOURCE_KEYS}


@pytest.mark.parametrize("construction", ["champion", "builder", "default"])
def test_the_pipeline_differs_from_the_cli_by_the_gates_and_nothing_else(
        construction, monkeypatch):
    from delfin import cli_manta

    pipeline = _pipeline_env({"MANTA_CONSTRUCTION": construction}, monkeypatch)
    cli = cli_manta.construction_env(construction, environ={})

    added = {k: v for k, v in pipeline.items() if cli.get(k) != v}
    assert added == _PIPELINE_ONLY
    assert set(cli) <= set(pipeline), sorted(set(cli) - set(pipeline))


def test_with_the_gates_off_the_pipeline_builds_what_the_cli_builds(monkeypatch):
    from delfin import cli_manta

    pipeline = _pipeline_env({"MANTA_CLEAN_GATE": "no", "MANTA_TOPOLOGY_GATE": "no",
                              "MANTA_DEDUP": "no"}, monkeypatch)
    cli = cli_manta.construction_env("champion", environ={})
    extra = {k: v for k, v in pipeline.items() if cli.get(k) != v}
    # what is left are "0" switches, i.e. the builder's own defaults
    assert set(extra) == set(_PIPELINE_ONLY) and set(extra.values()) == {"0"}
