"""Identity-based assay joins without unit conversion, pooling or endpoint inference."""
from __future__ import annotations

import math
import re
from pathlib import Path

from delfin.tools.bio.artifacts import file_hash, read_json


def validate_assays(value: dict) -> list:
    if not isinstance(value, dict) or value.get('schema_version') != 1 or not isinstance(value.get('records'), list):
        raise ValueError('assay table requires schema_version=1 and records')
    seen = set()
    for row in value['records']:
        if not isinstance(row, dict):
            raise ValueError('assay records must be objects')
        for key in ('measurement_id', 'chemical_state_id', 'assay_id', 'endpoint', 'unit', 'relation', 'source'):
            if not isinstance(row.get(key), str) or not row[key].strip():
                raise ValueError(f'assay record requires nonempty {key}')
        if not re.fullmatch('[0-9a-f]{64}', row['chemical_state_id']):
            raise ValueError('chemical_state_id must be the exact DELFIN SHA-256 identity')
        if row['measurement_id'] in seen:
            raise ValueError('measurement_id must be unique; replicates need distinct IDs')
        seen.add(row['measurement_id'])
        if row['relation'] not in ('=', '<', '<=', '>', '>=', '~'):
            raise ValueError('unsupported assay relation')
        number = row.get('value')
        if isinstance(number, bool) or not isinstance(number, (float, int)) or not math.isfinite(number):
            raise ValueError('assay value must be a finite number')
        if not isinstance(row.get('conditions'), dict):
            raise ValueError('assay conditions must be an explicit object (empty if unknown)')
    return value['records']


def join_assays(properties_file: Path, assays_file: Path) -> dict:
    properties = read_json(properties_file)
    assays = validate_assays(read_json(assays_file))
    if properties.get('schema_version') != 1 or properties.get('kind') != 'computed_molecular_descriptors':
        raise ValueError('properties must be a DELFIN computed descriptor table')
    records = properties.get('records')
    if not isinstance(records, list) or not records:
        raise ValueError('properties records must be a nonempty list')
    index = {}
    for record in records:
        identity = record.get('chemical_state_id')
        if not isinstance(identity, str) or not re.fullmatch('[0-9a-f]{64}', identity):
            raise ValueError('descriptor row has no valid chemical_state_id')
        index.setdefault(identity, []).append(record)
    joined, unmatched = [], []
    measured = set()
    for assay in assays:
        identity = assay['chemical_state_id']
        if identity not in index:
            unmatched.append(assay)
            continue
        measured.add(identity)
        for record in index[identity]:
            joined.append({'properties': record, 'measurement': assay})
    return {'schema_version': 1, 'records': joined, 'unmatched_measurements': unmatched,
            'unmeasured_properties': [record for record in records if record['chemical_state_id'] not in measured],
            'input_sha256': {'properties': file_hash(properties_file), 'assays': file_hash(assays_file)},
            'join_key': 'exact chemical_state_id; repeated structures and measurements retained',
            'normalization': 'none; endpoints, inequalities, units, conditions and provenance retained; no averaging'}
