"""The Submit tab chooses GOAT's method: xTB or an orca-mlips MLIP.

CONTROL.txt's ``global_opt_method`` is "" (GOAT with the xTB_method) or
``mlip:<plugin>``; ``global_opt_mlip_model`` optionally names the model.  The
Submit tab shows a selector beside the CONTROL editor that writes those lines
and follows the text back; availability of each plugin is read from the
environments' metadata (no torch), and the installer offers MACE and UMA as
opt-in, licensed extras.  Only GOAT's method changes: the rest of the CONTROL
is left byte for byte.
"""

from __future__ import annotations

import difflib
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin import installer
from delfin.common import external_builders as eb
from delfin.dashboard import goat_method as gm
from delfin.define import TEMPLATE
from delfin.mlp_tools import orca_mlips as om

SCRIPT = Path(installer.__file__).parent / "ai_tools" / "install_ai_tools.sh"


def _only_goat(text: str = TEMPLATE) -> str:
    from delfin.config import set_control_value

    for key, value in (("global_optimizer", "GOAT"), ("calc_initial", "no"),
                       ("method", "classic")):
        text = set_control_value(text, key, value)
    return text


def _changed(before: str, after: str) -> list:
    """The lines that differ, as unified-diff +/- lines."""
    return [line for line in difflib.unified_diff(
        before.splitlines(), after.splitlines(), lineterm="", n=0)
        if line[:1] in "+-" and not line.startswith(("+++", "---"))]


# ---- CONTROL text <-> choice -------------------------------------------------

def test_choosing_an_mlip_writes_the_method_and_selects_goat_and_nothing_else():
    after = gm.write_method(TEMPLATE, "mlip:orb")
    assert sorted(_changed(TEMPLATE, after)) == sorted([
        "-global_optimizer=[GOAT|CREST]",
        "+global_optimizer=GOAT",
        "+global_opt_method=mlip:orb",
    ])
    lines = after.splitlines()
    # In the global-optimisation section, right under the optimiser.
    assert lines.index("global_opt_method=mlip:orb") == lines.index("global_optimizer=GOAT") + 1
    assert gm.read_method(after) == ("mlip:orb", "mlip:orb")


def test_an_only_goat_control_keeps_its_optimiser_and_gains_one_line():
    before = _only_goat()
    after = gm.write_method(before, "mlip:aimnet2")
    assert _changed(before, after) == ["+global_opt_method=mlip:aimnet2"]


def test_choosing_xtb_empties_the_line_and_leaves_the_rest():
    with_mlip = gm.write_method(_only_goat(), "mlip:mace")
    back = gm.write_method(with_mlip, gm.XTB)
    assert _changed(with_mlip, back) == ["-global_opt_method=mlip:mace", "+global_opt_method="]
    assert gm.read_method(back) == (gm.XTB, "")


def test_xtb_on_a_text_without_the_line_changes_nothing():
    assert gm.write_method(TEMPLATE, gm.XTB) == TEMPLATE
    assert gm.read_method(TEMPLATE) == (gm.XTB, "")


def test_switching_between_mlips_rewrites_only_the_method_line():
    orb = gm.write_method(_only_goat(), "mlip:orb")
    uma = gm.write_method(orb, "mlip:uma")
    assert _changed(orb, uma) == ["-global_opt_method=mlip:orb", "+global_opt_method=mlip:uma"]


@pytest.mark.parametrize("raw, expected", [
    ("", gm.XTB), ("none", gm.XTB), ("xtb", gm.XTB),
    ("MLIP: ORB", "mlip:orb"), ("mlip:aimnet2", "mlip:aimnet2"),
])
def test_every_spelling_the_validator_accepts_is_read(raw, expected):
    text = _only_goat() + f"global_opt_method={raw}\n"
    assert gm.read_method(text)[0] == expected


def test_an_unknown_value_is_shown_as_such_and_never_replaced():
    text = _only_goat() + "global_opt_method=mlip:schnet\n"
    assert gm.read_method(text) == (gm.UNKNOWN, "mlip:schnet")
    labels = [label for label, _value in gm.options({}, "mlip:schnet")]
    assert labels[-1] == "Unknown: mlip:schnet"
    assert "mlip:schnet" in gm.status_html(gm.UNKNOWN, {}, "mlip:schnet")


def test_crest_is_switched_to_goat_for_an_mlip():
    from delfin.config import set_control_value

    crest = set_control_value(TEMPLATE, "global_optimizer", "CREST")
    after = gm.write_method(crest, "mlip:orb")
    assert sorted(_changed(crest, after)) == sorted([
        "-global_optimizer=CREST", "+global_optimizer=GOAT", "+global_opt_method=mlip:orb"])


def test_an_older_xtb_goat_switch_is_left_as_the_goat_choice():
    text = "xTB_method=XTB2\nXTB_GOAT=yes\nPAL=4\n"
    after = gm.write_method(text, "mlip:orb")
    assert _changed(text, after) == ["+global_opt_method=mlip:orb"]


def test_the_model_line_follows_the_model_field():
    with_mlip = gm.write_method(_only_goat(), "mlip:orb")
    with_model = gm.write_model(with_mlip, "orb_v3_direct_omol")
    assert _changed(with_mlip, with_model) == ["+global_opt_mlip_model=orb_v3_direct_omol"]
    lines = with_model.splitlines()
    assert lines.index("global_opt_mlip_model=orb_v3_direct_omol") == \
        lines.index("global_opt_method=mlip:orb") + 1
    assert gm.read_model(with_model) == "orb_v3_direct_omol"
    cleared = gm.write_model(with_model, "")
    assert _changed(with_model, cleared) == [
        "-global_opt_mlip_model=orb_v3_direct_omol", "+global_opt_mlip_model="]
    assert gm.write_model(with_mlip, "") == with_mlip
    assert gm.model_error("two words") and not gm.model_error("uma-s-1p1")


def test_what_the_selector_writes_passes_validation_and_is_sent_to_a_gpu():
    from delfin.config import validate_control_text
    from delfin.mlp_tools.goat_gpu_submit import mlip_goat_method

    text = gm.write_method(_only_goat(), "mlip:orb")
    errors = validate_control_text(text, converts_smiles=False)
    assert not [e for e in errors if "global_opt" in e or "global_optimizer" in e], errors
    assert mlip_goat_method(text) == "mlip:orb"
    assert mlip_goat_method(gm.write_method(text, gm.XTB)) == ""


# ---- availability, from metadata ----------------------------------------------

def _dist(site: Path, name: str, version: str) -> None:
    meta = site / f"{name}-{version}.dist-info"
    meta.mkdir(parents=True)
    meta.joinpath("METADATA").write_text(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")


def _env(root: Path, name: str, dists: dict, plugins=("aimnet2", "orb", "mace", "uma")) -> Path:
    """An environment of its own as install_ai_tools.sh leaves it (metadata only)."""
    prefix = root / ".mamba_env" / name
    (prefix / "bin").mkdir(parents=True)
    python = prefix / "bin" / "python"
    python.write_text("#!/bin/sh\nexit 1\n")
    python.chmod(0o755)
    site = prefix / "lib" / "python3.12" / "site-packages"
    site.mkdir(parents=True)
    for dist, version in dists.items():
        _dist(site, dist, version)
    for plugin in plugins:
        exe = prefix / "bin" / f"orca-mlips-{plugin}"
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)
    return site


@pytest.fixture()
def ai_root(tmp_path, monkeypatch):
    root = tmp_path / "ai_tools"
    root.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(root))
    monkeypatch.setenv("HOME", str(home))
    for var in ("DELFIN_ORCA_MLIPS_PYTHON", "DELFIN_ORCA_MLIPS_UMA_PYTHON", "HF_TOKEN",
                "HUGGING_FACE_HUB_TOKEN", "HF_TOKEN_PATH", "HF_HOME", "XDG_CACHE_HOME"):
        monkeypatch.delenv(var, raising=False)
    return root


def test_the_default_environment_has_aimnet2_and_orb_only(ai_root):
    _env(ai_root, "orca_mlips", {"orca_mlips": "1.3.0", "aimnet": "0.2.6", "orb_models": "0.5.5"})
    status = om.mlips_backend_status()
    assert status["aimnet2"]["installed"] and status["orb"]["installed"]
    assert status["orb"]["stack_version"] == "0.5.5"
    assert status["orb"]["progext"].endswith("/.mamba_env/orca_mlips/bin/orca-mlips-orb")
    mace = status["mace"]
    assert not mace["installed"] and mace["missing"] == "mace-torch is not installed"
    assert mace["install"] == "python -m delfin.installer --install orca-mlips-mace"
    assert mace["licence"]["non_commercial"]
    uma = status["uma"]
    assert not uma["installed"] and uma["install"] == "python -m delfin.installer --install orca-mlips-uma"
    assert uma["licence"]["gated"]
    assert gm.short_state("mlip:mace", status) == "not installed"
    html = gm.status_html("mlip:mace", status)
    assert "Not installed" in html and "--install orca-mlips-mace" in html
    assert "academic use only" in html and "GPU partition" in html


def test_nothing_installed_says_so_with_the_command(ai_root):
    status = om.mlips_backend_status()
    assert not any(entry["installed"] for entry in status.values())
    assert status["orb"]["missing"] == "orca-mlips is not installed"
    assert status["orb"]["install"] == "python -m delfin.installer --install orca-mlips"


def test_mace_is_installed_once_its_stack_is_in_the_environment(ai_root):
    site = _env(ai_root, "orca_mlips", {"orca_mlips": "1.3.0", "aimnet": "0.2.6",
                                       "orb_models": "0.5.5"})
    assert not installer.present(installer.find("orca-mlips-mace"))
    _dist(site, "mace_torch", "0.3.17")
    status = om.mlips_backend_status()
    assert status["mace"]["installed"] and status["mace"]["stack_version"] == "0.3.17"
    assert installer.present(installer.find("orca-mlips-mace"))
    # The default tool is still read from orca-mlips itself.
    assert installer.present(installer.find("orca-mlips"))


def test_uma_needs_its_environment_and_a_huggingface_login(ai_root, monkeypatch):
    _env(ai_root, "orca_mlips", {"orca_mlips": "1.3.0", "aimnet": "0.2.6", "orb_models": "0.5.5"})
    assert not installer.present(installer.find("orca-mlips-uma"))
    _env(ai_root, "orca_mlips_uma", {"orca_mlips": "1.3.0", "fairchem_core": "2.25.0"},
         plugins=("uma",))
    assert installer.present(installer.find("orca-mlips-uma"))
    uma = om.mlips_backend_status()["uma"]
    assert uma["installed"] and uma["progext"].endswith("/orca_mlips_uma/bin/orca-mlips-uma")
    assert uma["hf_login"] is False
    assert "hf auth login" in uma["missing"] and "huggingface-cli login" in uma["missing"]
    assert gm.short_state("mlip:uma", om.mlips_backend_status()) == \
        "installed, HuggingFace login needed"
    assert "huggingface-cli login" in gm.status_html("mlip:uma", om.mlips_backend_status())
    # ExtOpt names the UMA environment's command, not the default one's.
    assert om.resolve_progext("uma").endswith("/orca_mlips_uma/bin/orca-mlips-uma")
    assert om.resolve_progext("orb").endswith("/orca_mlips/bin/orca-mlips-orb")
    token = Path(os.environ["HOME"]) / ".cache" / "huggingface" / "token"
    token.parent.mkdir(parents=True)
    token.write_text("x")
    uma = om.mlips_backend_status()["uma"]
    assert uma["hf_login"] is True and uma["missing"] == ""
    token.unlink()
    monkeypatch.setenv("HF_TOKEN", "x")
    assert om.mlips_backend_status()["uma"]["hf_login"] is True


def test_availability_imports_no_torch():
    code = ("import sys\n"
            "from delfin.dashboard import goat_method as gm\n"
            "gm.availability(refresh=True)\n"
            "print('torch' in sys.modules)\n")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          timeout=300)
    assert done.returncode == 0, done.stderr[-800:]
    assert done.stdout.strip().splitlines()[-1] == "False"


def test_availability_is_kept_between_calls(monkeypatch):
    calls = []

    def fake(environ=None):
        calls.append(1)
        return {"orb": {"installed": True}}

    monkeypatch.setattr(om, "mlips_backend_status", fake)
    monkeypatch.setattr(gm, "_CACHE", {"at": None, "status": None})
    gm.availability()
    gm.availability()
    assert len(calls) == 1
    gm.availability(refresh=True)
    assert len(calls) == 2


# ---- the installer -------------------------------------------------------------

def test_mace_and_uma_are_licensed_opt_in_tools_of_the_catalog():
    mace = installer.find("orca-mlips-mace")
    uma = installer.find("orca-mlips-uma")
    assert mace is installer.find("orca_mlips_mace") and uma is installer.find("orca_mlips_uma")
    for tool, switch in ((mace, "INSTALL_ORCA_MLIPS_MACE"), (uma, "INSTALL_ORCA_MLIPS_UMA")):
        assert tool.licensed and tool.group == "ai" and tool.switch == switch
        assert switch in installer.switches("ai")
        assert tool.name not in installer.profile("all")
        assert tool.notice
    assert mace.own_env == "orca_mlips" and mace.dist == "mace-torch"
    assert "Academic Software License" in mace.notice and "non-commercial" in mace.notice
    assert uma.own_env == "orca_mlips_uma"
    assert "hf auth login" in uma.notice and "huggingface-cli login" in uma.notice
    assert "https://huggingface.co/facebook/UMA" in uma.notice
    info = eb.PROGRAMS["orca_mlips_uma"]
    assert info["managed_env"] == ".mamba_env/orca_mlips_uma"
    assert info["python_env"] == "DELFIN_ORCA_MLIPS_UMA_PYTHON"
    assert info["install"] == "python -m delfin.installer --install orca-mlips-uma"
    # The default stays: aimnet2 + orb, unlicensed, in "all".
    default = installer.find("orca-mlips")
    assert not default.licensed and "orca_mlips" in installer.profile("all")


def test_installing_a_licensed_plugin_states_its_licence_first(monkeypatch, tmp_path):
    calls, said = [], []
    monkeypatch.setattr(installer, "_stage", lambda group: tmp_path / group)
    monkeypatch.setattr(installer, "present", lambda tool: True)

    def run(command, *, cwd, env, on_line, timeout):
        calls.append(env)
        said.append("<installer ran>")
        return True, []

    monkeypatch.setattr(installer, "_run", run)
    outcome = installer.install(["orca-mlips-mace"], on_line=said.append)
    assert outcome["ok"]
    assert said.index("<installer ran>") > 0
    assert "Academic Software License" in said[0]
    env = calls[0]
    assert env["INSTALL_ORCA_MLIPS_MACE"] == "1" and env["INSTALL_ORCA_MLIPS_UMA"] == "0"
    said.clear()
    installer.install(["orca-mlips-uma"], on_line=said.append)
    assert "huggingface-cli login" in said[0] and calls[-1]["INSTALL_ORCA_MLIPS_UMA"] == "1"


def test_the_script_installs_the_extras_pinned_opt_in_and_with_their_terms():
    text = SCRIPT.read_text()
    # Opt-in: never switched on by INSTALL_ALL.
    assert 'INSTALL_ORCA_MLIPS_MACE="${INSTALL_ORCA_MLIPS_MACE:-0}"' in text
    assert 'INSTALL_ORCA_MLIPS_UMA="${INSTALL_ORCA_MLIPS_UMA:-0}"' in text
    mace = next(l for l in text.splitlines() if l.startswith("ORCA_MLIPS_MACE_PIP_DEPS="))
    uma = next(l for l in text.splitlines() if l.startswith("ORCA_MLIPS_UMA_PIP_DEPS="))
    # Exact pins; MACE keeps the default environment's torch and numpy.
    assert "mace-torch==0.3.17" in mace and "torch==2.15.0" in mace and "numpy==2.4.6" in mace
    assert "fairchem-core==2.25.0" in uma and "torch==2.13.0" in uma
    # MACE goes into the orca-mlips environment, which is built first.
    body = text[text.index("install_orca_mlips() {"):text.index("install_orca_mlips_mace() {")]
    assert 'if [ "${INSTALL_ORCA_MLIPS_MACE}" = "1" ]; then' in body
    assert body.index("INSTALL_ORCA_MLIPS=1") < body.index("install_own_env INSTALL_ORCA_MLIPS ")
    mace_fn = text[text.index("install_orca_mlips_mace() {"):text.index("install_orca_mlips_uma() {")]
    assert "Academic Software License" in mace_fn
    assert mace_fn.index("LICENCE") < mace_fn.index("pip install")
    uma_fn = text[text.index("install_orca_mlips_uma() {"):]
    assert "huggingface-cli login" in uma_fn and "auth login" in uma_fn
    assert "install_own_env INSTALL_ORCA_MLIPS_UMA orca-mlips-uma orca_mlips_uma" in uma_fn
    # The default environment's own pins are unchanged.
    deps = next(l for l in text.splitlines() if l.startswith("ORCA_MLIPS_PIP_DEPS="))
    assert "mace" not in deps and "fairchem" not in deps


# ---- the Submit tab --------------------------------------------------------------

def _walk(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _walk(child)


@pytest.fixture()
def submit_tab(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(gm, "availability", lambda **_kw: {
        "orb": {"installed": True, "stack": "orb-models", "stack_version": "0.5.5",
                "orca_mlips_version": "1.3.0", "default_model": "orb_v3_conservative_omol",
                "licence": {"name": "Apache-2.0"}},
        "mace": {"installed": False, "missing": "mace-torch is not installed",
                 "install": "python -m delfin.installer --install orca-mlips-mace",
                 "default_model": "MACE-OMOL-0",
                 "licence": {"name": "Academic Software License (non-commercial)",
                             "non_commercial": True}},
    })
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.default_control = TEMPLATE
    ctx.run_js = lambda _script: None
    widget, refs = tab_submit.create_tab(ctx)
    return widget, refs


def test_the_selector_stands_beside_the_control_editor(submit_tab):
    widget, refs = submit_tab
    shown = list(_walk(widget))
    select, model = refs["goat_method_widget"], refs["goat_model_widget"]
    assert any(w is select for w in shown) and any(w is model for w in shown)
    # Labelled beside it, the tooltip on the select itself: with a description
    # the widget front end shows the tooltip nowhere.
    assert refs["goat_method_label"].value == "GOAT method:"
    assert any(w is refs["goat_method_label"] for w in shown)
    assert select.description == "" and model.description == ""
    values = [value for _label, value in select.options]
    assert values == ["", "mlip:aimnet2", "mlip:orb", "mlip:mace", "mlip:uma"]
    labels = dict((value, label) for label, value in select.options)
    assert labels[""] == "xTB (as xTB_method)"
    assert labels["mlip:orb"] == "ORB - installed" and labels["mlip:mace"] == "MACE - not installed"
    for name in ("AIMNet2", "ORB", "MACE", "UMA", "GPU partition"):
        assert name in select.tooltip
    assert select.value == "" and model.disabled


def test_choosing_orb_writes_it_and_editing_the_text_back_selects_xtb(submit_tab):
    _widget, refs = submit_tab
    select, control = refs["goat_method_widget"], refs["control_widget"]
    before = control.value
    select.value = "mlip:orb"
    assert "global_opt_method=mlip:orb" in control.value.splitlines()
    assert sorted(_changed(before, control.value)) == sorted([
        "-global_optimizer=[GOAT|CREST]", "+global_optimizer=GOAT", "+global_opt_method=mlip:orb"])
    assert not refs["goat_model_widget"].disabled
    assert "Installed" in refs["goat_method_status"].value
    control.value = control.value.replace("global_opt_method=mlip:orb", "global_opt_method=")
    assert select.value == ""
    control.value = control.value.replace("global_opt_method=", "global_opt_method=mlip:mace")
    assert select.value == "mlip:mace"
    assert "--install orca-mlips-mace" in refs["goat_method_status"].value


def test_an_unknown_value_in_the_text_is_shown_and_kept(submit_tab):
    _widget, refs = submit_tab
    select, control = refs["goat_method_widget"], refs["control_widget"]
    text = gm.write_method(control.value, "mlip:orb").replace(
        "global_opt_method=mlip:orb", "global_opt_method=mlip:schnet")
    control.value = text
    assert select.value == gm.UNKNOWN
    assert ("Unknown: mlip:schnet", gm.UNKNOWN) in list(select.options)
    assert control.value == text


def test_the_model_field_and_the_model_line_follow_each_other(submit_tab):
    _widget, refs = submit_tab
    select, model, control = (refs["goat_method_widget"], refs["goat_model_widget"],
                              refs["control_widget"])
    select.value = "mlip:orb"
    assert model.placeholder == "default: orb_v3_conservative_omol"
    model.value = "orb_v3_direct_omol"
    assert "global_opt_mlip_model=orb_v3_direct_omol" in control.value.splitlines()
    control.value = control.value.replace("orb_v3_direct_omol", "orb_v3_conservative_omol")
    assert model.value == "orb_v3_conservative_omol"


def test_reset_puts_the_selector_back_to_xtb(submit_tab):
    _widget, refs = submit_tab
    refs["goat_method_widget"].value = "mlip:orb"
    refs["reset_form"]()
    assert refs["goat_method_widget"].value == ""
    assert "global_opt_method" not in refs["control_widget"].value


# ---- where an MLIP GOAT is sent ----------------------------------------------------

def test_an_mlip_goat_is_never_sent_to_an_amd_gpu_partition():
    """Measured: routed to an MI300 partition, the CUDA build of the MLIP saw no
    /dev/nvidia* and GOAT stopped after 30 s."""
    from delfin.mlp_tools.goat_gpu_submit import cuda_partitions, gpu_partitions_for

    kept, dropped = cuda_partitions(["gpu_h100", "gpu_mi300", "gpu_a100_il", "gpu_mi250x"])
    assert kept == ("gpu_h100", "gpu_a100_il") and dropped == ("gpu_mi300", "gpu_mi250x")
    partitions, reason = gpu_partitions_for(["gpu_mi300"], "00:30:00")
    assert partitions == "" and "AMD" in reason and "NVIDIA" in reason
