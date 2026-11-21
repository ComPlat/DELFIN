import copy
import json
from pathlib import Path

import pytest

from delfin.tools.bio.execution import execute_plan
from delfin.tools.bio.processing import ligand_properties
from delfin.tools.bio.workflows import workflow_template


@pytest.fixture
def setup(tmp_path):
    ligand = tmp_path / 'library.smi'
    ligand.write_text('CCO ethanol\nCC ethane\n')
    identity = ligand_properties(ligand)['records'][0]['chemical_state_id']
    row = dict(measurement_id='exp1', chemical_state_id=identity, assay_id='target1', endpoint='IC50',
               value=12, unit='nM', relation='>', source='synthetic fixture', conditions={'pH': 7.4})
    assays = tmp_path / 'assays.json'
    plan = workflow_template('medchem')
    plan['inputs']['ligands']['path'] = str(ligand)
    plan['inputs']['assays']['path'] = str(assays)
    return row, assays, plan


def test_medchem_dag_keeps_replicates_censoring_units_and_unmatched_rows(setup, tmp_path):
    row, assays, plan = setup
    replicate = {**row, 'measurement_id': 'exp2', 'endpoint': 'Ki', 'unit': 'uM', 'value': 0.5, 'relation': '='}
    unmatched = {**row, 'measurement_id': 'exp3', 'chemical_state_id': 'f' * 64}
    assays.write_text(json.dumps({'schema_version': 1, 'records': [row, replicate, unmatched]}))
    result = execute_plan(plan, tmp_path / 'run')
    assert result['status'] == 'success', result
    dataset = result['steps'][-1]['data']
    assert [r['measurement'] for r in dataset['records']] == [row, replicate]
    assert dataset['unmatched_measurements'] == [unmatched]
    assert len(dataset['unmeasured_properties']) == 1


@pytest.mark.parametrize('field', ['source', 'unit', 'conditions', 'chemical_state_id'])
def test_incomplete_measurements_block_before_execution(setup, tmp_path, field):
    row, assays, plan = setup
    del row[field]
    assays.write_text(json.dumps({'schema_version': 1, 'records': [row]}))
    result = execute_plan(plan, tmp_path / 'run')
    assert result['status'] == 'blocked'
    assert not (tmp_path / 'run').exists()


def test_duplicate_measurement_ids_are_not_silently_averaged(setup, tmp_path):
    row, assays, plan = setup
    assays.write_text(json.dumps({'schema_version': 1, 'records': [row, copy.deepcopy(row)]}))
    result = execute_plan(plan, tmp_path / 'run')
    assert result['status'] == 'blocked'
