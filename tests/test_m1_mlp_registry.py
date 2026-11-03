"""m1 — unified MLP backend registry (delfin/mlp_tools/registry.py).

Covers: registry contents, availability reporting with packages missing,
device="auto" policy, charge/spin refusal for non-charge-aware backends,
GPUMA unavailable entry, and the unchanged legacy create_calculator surface.
All tests run without torch/ML packages (fake ASE calculators + fake torch).
"""

import sys
import types

import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator

from delfin.mlp_tools import calculators as legacy
from delfin.mlp_tools import registry


# ---------------------------------------------------------------- fakes

class FakeCalc(Calculator):
    """Minimal ASE calculator recording construction kwargs."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, device="cpu", charge=0, mult=1, **kwargs):
        super().__init__()
        self.init_device = device
        self.init_charge = charge
        self.init_mult = mult

    def calculate(self, *args, **kwargs):
        self.results = {"energy": -1.0, "forces": None}


@pytest.fixture
def fake_factory(monkeypatch):
    """Replace one backend's factory with a fake; returns the calls list."""
    calls = []

    def fake(device="cpu", charge=0, mult=1, **kwargs):
        calls.append({"device": device, "charge": charge, "mult": mult})
        return FakeCalc(device=device, charge=charge, mult=mult, **kwargs)

    return fake, calls


class _FakeCuda:
    @staticmethod
    def is_available():
        return False


def _install_fake_torch(monkeypatch, cuda_available):
    fake = types.ModuleType("torch")
    fake.cuda = types.SimpleNamespace(is_available=lambda: cuda_available)
    monkeypatch.setitem(sys.modules, "torch", fake)


# ---------------------------------------------------------------- registry contents

def test_all_expected_backends_registered():
    names = set(registry.backend_names())
    expected = {
        # the eight legacy calculators.py backends
        "ani2x", "aimnet2", "mace_off", "chgnet",
        "m3gnet", "schnetpack", "nequip", "alignn",
        # state-of-the-art molecular models
        "mace_mp", "mace_omol", "fairchem_uma", "esen",
        "orb", "sevennet", "xtb",
        # GPU batch optimiser entry point
        "gpuma",
    }
    missing = expected - names
    assert not missing, f"missing backends: {sorted(missing)}"


def test_legacy_aliases_resolve_to_same_backend():
    for alias, key in [
        ("ani-2x", "ani2x"), ("torchani", "ani2x"),
        ("mace-off", "mace_off"), ("mace", "mace_off"),
        ("m3gnet", "m3gnet"), ("matgl", "m3gnet"),
        ("schnet", "schnetpack"), ("painn", "schnetpack"),
        ("allegro", "nequip"),
    ]:
        assert registry.resolve(alias) == key, alias


def test_every_entry_has_capabilities_and_refs():
    for name in registry.backend_names():
        spec = registry.get_backend(name)
        assert spec.capabilities.energy and spec.capabilities.forces, name
        assert spec.scope in ("molecular", "periodic", "both"), name
        assert spec.device_policy, name
        assert spec.literature, f"{name}: literature reference keys required"
        assert spec.description, name


def test_molecular_models_are_marked_molecular():
    for name in ["ani2x", "aimnet2", "mace_mp", "mace_omol", "mace_off", "orb"]:
        assert registry.get_backend(name).scope == "molecular", name


def test_charge_aware_flag_matches_adapters_contract():
    # only aimnet2 accepts charge/mult (fail-closed contract of adapters/mlp.py)
    assert registry.get_backend("aimnet2").charge_aware
    for name in ["ani2x", "mace_mp", "mace_omol", "mace_off", "chgnet",
                 "m3gnet", "schnetpack", "nequip", "alignn", "fairchem_uma",
                 "esen", "orb", "sevennet"]:
        assert not registry.get_backend(name).charge_aware, name


def test_gpuma_registered_unavailable_with_reason():
    spec = registry.get_backend("gpuma")
    info = registry.describe_backend("gpuma")
    assert spec is not None
    assert info["available"] is False
    assert info["reason"], "unavailable backends must state why"
    assert info["install_hint"]


# ---------------------------------------------------------------- availability

def test_availability_report_with_all_packages_missing():
    # both CI interpreters have no torch/ML packages; availability must be
    # reported, never crash, never guess
    report = registry.describe_backends()
    by_name = {e["name"]: e for e in report}
    assert set(by_name) == set(registry.backend_names())
    for entry in report:
        assert isinstance(entry["available"], bool)
        if not entry["available"]:
            assert entry["reason"], entry["name"]
            assert entry["install_hint"], entry["name"]


def test_available_backend_reports_no_reason(monkeypatch):
    spec = registry.get_backend("ani2x")
    monkeypatch.setattr(registry, "_package_available", lambda pkgs: True)
    info = registry.describe_backend("ani2x")
    assert info["available"] is True
    assert spec.requires_packages, "availability probe needs declared packages"


def test_availability_of_missing_package_lists_install_hint():
    info = registry.describe_backend("mace_mp")
    assert info["available"] is False
    assert "mace-torch" in info["install_hint"]


def test_availability_check_never_crashes_on_import_error(monkeypatch):
    def boom(pkgs):
        raise ImportError("simulated broken install")
    monkeypatch.setattr(registry, "_package_available", boom)
    info = registry.describe_backend("ani2x")
    assert info["available"] is False
    assert info["reason"]


# ---------------------------------------------------------------- device policy

def test_device_auto_resolves_to_cpu_without_cuda(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    _install_fake_torch(monkeypatch, cuda_available=False)
    registry.create_calculator("ani2x", device="auto")
    assert calls[-1]["device"] == "cpu"


def test_device_auto_resolves_to_cuda_when_available(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    _install_fake_torch(monkeypatch, cuda_available=True)
    registry.create_calculator("ani2x", device="auto")
    assert calls[-1]["device"] == "cuda"


def test_device_auto_without_torch_falls_back_to_cpu(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    monkeypatch.setitem(sys.modules, "torch", None)  # import torch -> ImportError
    registry.create_calculator("ani2x", device="auto")
    assert calls[-1]["device"] == "cpu"


def test_explicit_device_wins_over_auto(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    _install_fake_torch(monkeypatch, cuda_available=False)
    registry.create_calculator("ani2x", device="cuda")
    assert calls[-1]["device"] == "cuda"


def test_auto_cpu_caps_threads(monkeypatch, fake_factory):
    """device="auto" on a CPU-only host caps OMP/MKL threads at 2."""
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.delenv(var, raising=False)
    _install_fake_torch(monkeypatch, cuda_available=False)
    registry.create_calculator("ani2x", device="auto")
    import os
    assert os.environ["OMP_NUM_THREADS"] == "2"
    assert os.environ["MKL_NUM_THREADS"] == "2"


def test_default_device_stays_cpu(monkeypatch, fake_factory):
    """Old default path: device defaults to "cpu", no auto resolution."""
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    _install_fake_torch(monkeypatch, cuda_available=True)  # would pick cuda if auto
    registry.create_calculator("ani2x")
    assert calls[-1]["device"] == "cpu"


# ---------------------------------------------------------------- charge / spin

def test_charge_refused_for_non_charge_aware_backends(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"ani2x": fake})
    with pytest.raises(ValueError, match="charge"):
        registry.create_calculator("ani2x", charge=-1)
    assert calls == [], "factory must not run when charge is refused"


def test_multiplicity_refused_for_non_charge_aware_backends(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"mace_mp": fake})
    with pytest.raises(ValueError, match="mult"):
        registry.create_calculator("mace_mp", mult=3)
    assert calls == []


def test_charge_passed_to_aimnet2(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"aimnet2": fake})
    registry.create_calculator("aimnet2", charge=-1, mult=2)
    assert calls[-1]["charge"] == -1
    assert calls[-1]["mult"] == 2


def test_neutral_charge_mult_accepted_everywhere(monkeypatch, fake_factory):
    fake, calls = fake_factory
    monkeypatch.setattr(registry, "_FACTORY_OVERRIDE", {"mace_mp": fake})
    registry.create_calculator("mace_mp", charge=0, mult=1)
    assert calls[-1]["charge"] == 0


# ---------------------------------------------------------------- legacy surface

def test_old_create_calculator_names_unchanged():
    calc = legacy.create_calculator("ani2x", device="cpu",
                                    charge=0, mult=1,
                                    _factory_override=FakeCalc)
    assert isinstance(calc, FakeCalc)
    assert calc.init_device == "cpu"


def test_legacy_alias_table_unchanged():
    assert legacy._BACKEND_ALIASES["mace"] == "mace_off"
    assert legacy._BACKEND_ALIASES["torchani"] == "ani2x"
    assert len(legacy._BACKEND_ALIASES) == 18  # exact table, unchanged


def test_legacy_single_point_default_backend_ani2x():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]])
    # inject a fake via the documented override hook; no torch needed
    calc = legacy.create_calculator("ani2x", _factory_override=lambda **kw: FakeCalc())
    atoms.calc = calc
    assert atoms.get_potential_energy() == pytest.approx(-1.0)
    # signature default documented and locked:
    import inspect
    sig = inspect.signature(legacy.single_point)
    assert sig.parameters["backend"].default == "ani2x"


def test_legacy_unknown_backend_error_message():
    with pytest.raises(ValueError, match="Unknown MLP backend"):
        legacy.create_calculator("not_a_backend")


def test_registry_accepts_legacy_names_too():
    for name in ["ani2x", "aimnet2", "mace_off", "chgnet", "m3gnet",
                 "schnetpack", "nequip", "alignn"]:
        registry.get_backend(name)  # must not raise
