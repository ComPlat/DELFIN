from pathlib import Path
import shutil

import pytest

from delfin.tools.bio.execution import execute_plan, preflight_plan


@pytest.fixture
def plan(tmp_path):
    source = Path(__file__).resolve().parents[1] / "examples" / "biomolecular_foundation" / "receptor.pdb"
    shutil.copyfile(source, tmp_path / "input.pdb")
    return {"schema_version": 1, "inputs": {"source": {"type": "macromolecule", "path": str(tmp_path / "input.pdb")}},
            "nodes": [{"id": "select", "backend": "structure_select", "inputs": {"structure": "inputs.source"}, "settings": {"chains": ["A"]}},
                      {"id": "view", "backend": "visualize", "inputs": {"structure": "select.structure"}}]}


def test_actual_dag_executes_and_resume_checks_files(plan, tmp_path, monkeypatch):
    assert preflight_plan(plan)["execution_ready"]
    result = execute_plan(plan, tmp_path / "run")
    assert result["status"] == "success", result
    assert [row["id"] for row in result["steps"]] == ["select", "view"]
    import delfin.tools
    monkeypatch.setattr(delfin.tools, "run_step", lambda *a, **k: pytest.fail("successful nodes must not rerun"))
    assert execute_plan(plan, tmp_path / "run", resume=True)["status"] == "success"
    view = tmp_path / "run" / result["steps"][1]["outputs"]["view"]["path"]
    view.write_text("changed")
    with pytest.raises(ValueError, match="output changed"):
        execute_plan(plan, tmp_path / "run", resume=True)


def test_input_or_settings_changes_invalidate_resume(plan, tmp_path):
    execute_plan(plan, tmp_path / "run")
    plan["nodes"][0]["settings"]["waters"] = "drop"
    with pytest.raises(ValueError, match="changed"):
        execute_plan(plan, tmp_path / "run", resume=True)


def test_bad_settings_and_unimplemented_engines_block_before_workdir(plan, tmp_path):
    plan["nodes"][0]["settings"]["imaginary_flag"] = True
    result = execute_plan(plan, tmp_path / "blocked")
    assert result["status"] == "blocked"
    assert not (tmp_path / "blocked").exists()


def test_failed_node_prevents_its_consumer(plan, tmp_path):
    plan["nodes"][0]["settings"]["chains"] = ["missing"]
    result = execute_plan(plan, tmp_path / "failed")
    assert result["status"] == "failed"
    assert len(result["steps"]) == 1
    assert not (tmp_path / "failed" / "view_1").exists()


def test_cancellation_before_first_node_has_no_execution(plan, tmp_path):
    result = execute_plan(plan, tmp_path / "cancelled", cancelled=lambda: True)
    assert result["status"] == "cancelled"
    assert not result["execution_performed"]


def test_optional_viewer_input_changes_invalidate_resume(plan, tmp_path):
    fixture = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation' / 'poses.sdf'
    source = tmp_path / 'poses.sdf'
    shutil.copyfile(fixture, source)
    plan['nodes'][1]['settings'] = {'poses_file': str(source)}
    assert execute_plan(plan, tmp_path / 'run')['status'] == 'success'
    source.write_text(source.read_text().replace('$$$$', '>  <changed>\nyes\n\n$$$$'))
    with pytest.raises(ValueError, match='changed'):
        execute_plan(plan, tmp_path / 'run', resume=True)


def test_execution_is_available_through_registered_adapter(plan, tmp_path):
    from delfin.tools import run_step
    result = run_step('bio_workflow_execute', plan=plan, work_dir=tmp_path / 'adapter')
    assert result.ok, result.error
    assert result.data['status'] == 'success'
    assert result.artifacts['json_report'].is_file()


def test_every_execution_binding_agrees_with_registered_ports():
    from delfin.tools._registry import get
    from delfin.tools.bio.execution import BINDINGS
    from delfin.tools.bio.workflows import BACKENDS
    for backend in BACKENDS:
        if backend.id not in BINDINGS:
            assert backend.adapter is None
            continue
        binding = BINDINGS[backend.id]
        adapter = get(binding.adapter)
        assert adapter is not None, backend.id
        assert set(binding.inputs) == set(dict(backend.inputs)), backend.id
        assert set(binding.outputs) == set(dict(backend.outputs)), backend.id
        assert set(binding.inputs.values()) <= {p.name for p in adapter.contract().params}, backend.id
        assert set(binding.outputs.values()) <= set(adapter.produces), backend.id


def test_only_execution_can_provision_and_bad_settings_never_do(plan, tmp_path, monkeypatch):
    from delfin.tools._registry import get
    adapter = get('bio_structure_select')
    calls = []
    def probe(**kwargs):
        calls.append(kwargs.get('_allow_install'))
        return {'fixture': True}
    monkeypatch.setattr(adapter, 'runtime_preflight', probe, raising=False)
    assert preflight_plan(plan)['execution_ready']
    assert calls == [False]
    calls.clear()
    assert execute_plan(plan, tmp_path / 'run')['status'] == 'success'
    assert calls == [True]
    calls.clear()
    plan['nodes'][1]['settings'] = {'unknown_flag': True}
    assert execute_plan(plan, tmp_path / 'invalid')['status'] == 'blocked'
    assert calls == []


def test_posebusters_dag_provisions_only_in_execution_phase(monkeypatch):
    import delfin.tools.adapters.bio_posebusters as module
    fixture = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation'
    plan = {'schema_version': 1,
            'inputs': {name: {'type': kind, 'path': str(fixture / file)} for name, kind, file in (
                ('poses', 'pose_set', 'poses.sdf'), ('reference', 'ligand', 'reference.sdf'),
                ('receptor', 'prepared_receptor', 'receptor.pdb'))},
            'nodes': [{'id': 'check', 'backend': 'posebusters',
                       'inputs': {key: 'inputs.' + key for key in ('poses', 'reference', 'receptor')},
                       'settings': {'expected_version': '0.5.1'}}]}
    attempts = []
    def probe(*args, **kwargs):
        attempts.append(kwargs.get('allow_install', False))
        if not attempts[-1]:
            raise ImportError('PoseBusters not installed')
        return {'version': '0.5.1', 'python': '/managed/bin/python'}
    monkeypatch.setattr(module, 'probe_distribution', probe)
    assert not preflight_plan(plan)['execution_ready']
    assert preflight_plan(plan, provision=True)['execution_ready']
    assert attempts == [False, True]
    plan['nodes'][0]['settings']['unsupported_option'] = True
    assert not preflight_plan(plan, provision=True)['execution_ready']
    assert attempts == [False, True]
