"""Runnable offline processing DAG and portable receptor/pose visualization."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from delfin.tools.bio import BindingBox
from delfin.tools.bio.artifacts import write_json
from delfin.tools.bio.execution import execute_plan
from delfin.tools.bio.workflows import workflow_template


def run_processing(output: Path) -> dict:
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('use an empty output directory')
    fixtures = Path(__file__).resolve().parent
    box = output / 'box.json'
    write_json(box, BindingBox(center=(0, 0, 0), size=(15, 15, 15)).to_dict())
    plan = workflow_template('processing')
    plan['inputs']['structure']['path'] = str(fixtures / 'receptor.pdb')
    plan['inputs']['ligands']['path'] = str(fixtures / 'poses.sdf')
    plan['inputs']['poses'] = {'type': 'pose_set', 'path': str(fixtures / 'poses.sdf')}
    plan['inputs']['site'] = {'type': 'binding_box', 'path': str(box)}
    plan['nodes'][1] = {'id': 'view', 'backend': 'visualize_poses',
                        'inputs': {'structure': 'select.structure', 'poses': 'inputs.poses', 'site': 'inputs.site'}}
    write_json(output / 'plan.json', plan)
    report = execute_plan(plan, output / 'run')
    if report['status'] != 'success':
        raise RuntimeError(json.dumps(report))
    return {'status': report['status'], 'viewer': str(output / 'run' / 'view_1' / 'view.html'),
            'manifest': str(output / 'run' / 'run.json'), 'scientific_validation': 'not_performed'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(run_processing(parser.parse_args().output), indent=2))
