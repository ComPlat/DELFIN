"""Public schemas for agent-authored plans and prediction requests."""
from __future__ import annotations


def prediction_request_schema() -> dict:
    entity_variants = []
    for kind, field, pattern in (("protein", "sequence", "^[ACDEFGHIKLMNPQRSTVWYX]+$"),
                                 ("dna", "sequence", "^[ACGTN]+$"),
                                 ("rna", "sequence", "^[ACGUN]+$"),
                                 ("ligand", "smiles", r"^\S+$")):
        entity_variants.append({"type": "object", "additionalProperties": False,
                                "required": ["id", "kind", field],
                                "properties": {"id": {"type": "string", "pattern": "^[A-Z]+$"},
                                               "kind": {"const": kind},
                                               field: {"type": "string", "pattern": pattern}}})
    return {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "object",
            "additionalProperties": False,
            "required": ["schema_version", "name", "entities", "seeds", "msa_policy"],
            "properties": {
                "schema_version": {"const": 1, "type": "integer"},
                "name": {"type": "string", "pattern": "^[A-Za-z][A-Za-z0-9_-]{0,79}$"},
                "entities": {"type": "array", "minItems": 1, "items": {"oneOf": entity_variants}},
                "seeds": {"type": "array", "minItems": 1, "uniqueItems": True,
                          "items": {"type": "integer", "minimum": 0, "maximum": 2**31 - 1}},
                "msa_policy": {"enum": ["single_sequence", "engine_pipeline"]},
            }}


def plan_schema() -> dict:
    from delfin.tools.bio.workflows import BACKENDS, PORT_TYPES
    return {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "object",
            "additionalProperties": False, "required": ["schema_version", "inputs", "nodes"],
            "properties": {
                "schema_version": {"const": 1, "type": "integer"}, "name": {"type": "string"},
                "inputs": {"type": "object", "propertyNames": {"pattern": "^[A-Za-z][A-Za-z0-9_]*$"},
                           "additionalProperties": {"type": "object", "additionalProperties": False,
                                                    "required": ["type"], "properties": {
                                                        "type": {"enum": sorted(PORT_TYPES)},
                                                        "path": {"type": "string"},
                                                        "description": {"type": "string"}}}},
                "nodes": {"type": "array", "minItems": 1, "items": {
                    "type": "object", "additionalProperties": False, "required": ["id", "backend", "inputs"],
                    "properties": {"id": {"type": "string", "pattern": "^[A-Za-z][A-Za-z0-9_]*$"},
                                   "backend": {"enum": [b.id for b in BACKENDS]},
                                   "inputs": {"type": "object", "additionalProperties": {"type": "string", "pattern": "^[A-Za-z][A-Za-z0-9_]*\\.[A-Za-z][A-Za-z0-9_]*$"}},
                                   "settings": {"type": "object"}}}},
            }}


def assay_table_schema() -> dict:
    return {'$schema': 'https://json-schema.org/draft/2020-12/schema', 'type': 'object',
            'required': ['schema_version', 'records'], 'properties': {
                'schema_version': {'const': 1, 'type': 'integer'},
                'records': {'type': 'array', 'items': {'type': 'object',
                    'required': ['measurement_id', 'chemical_state_id', 'assay_id', 'endpoint', 'unit', 'relation', 'source', 'value', 'conditions'],
                    'properties': {
                        **{key: {'type': 'string', 'minLength': 1} for key in ('measurement_id', 'assay_id', 'endpoint', 'unit', 'source')},
                        'chemical_state_id': {'type': 'string', 'pattern': '^[0-9a-f]{64}$'},
                        'relation': {'enum': ['=', '<', '<=', '>', '>=', '~']},
                        'value': {'type': 'number'}, 'conditions': {'type': 'object'}}}}}}
