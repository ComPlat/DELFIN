"""F2 red test — F2-R3: declared fact semantics must match the exact scorer.

The metric scores a fact with exact-token containment unconditionally
(metrics.py:_exact_hit; a paraphrase is NOT a recall). The task yaml and the
seed therefore must not declare ``must_be_exact: false`` on a fact the scorer
cannot honour — that is a silently ignored semantic. The red control walks
every planted fact in tasks_w28_lh.yaml and in lh_seed.py::_SEEDS and asserts
each one declares exact scoring (``must_be_exact`` truthy). Currently 6 facts
across the two sources declare ``false`` while being scored exactly.

No model call; deterministic; CI-safe.
"""
from __future__ import annotations

import yaml

from pathlib import Path


def _all_facts() -> list[dict]:
    """Collect every planted fact from the yaml and from lh_seed._SEEDS."""
    facts: list[dict] = []
    yaml_path = Path("delfin/agent/pack/benchmark/tasks_w28_lh.yaml")
    for doc in yaml.safe_load(yaml_path.read_text("utf-8")) or []:
        if not isinstance(doc, dict):
            continue
        for f in (doc.get("lh_") or {}).get("planted_facts") or []:
            facts.append({**f, "_task": doc.get("id")})

    import importlib
    seed = importlib.import_module(
        "delfin.agent.pack.benchmark.setup.lh_seed")
    for task_id, (_fixture, _text, meta) in seed._SEEDS.items():
        for f in meta.get("planted_facts") or []:
            facts.append({**f, "_task": task_id})
    return facts


def test_every_planted_fact_declares_exact_scoring():
    """No planted fact may claim paraphrase-tolerance (must_be_exact:false):
    the metric scores exact-token only, so such a fact's declared semantic
    would be silently dropped by _facts_to_metric (F2-R3)."""
    facts = _all_facts()
    assert facts, "expected at least one planted fact"
    lax = [f for f in facts if not bool(f.get("must_be_exact", True))]
    assert not lax, (
        "these facts declare must_be_exact:false but are scored exactly; "
        "relabel them must_be_exact:true so the declared semantics match the "
        "metric: " + "; ".join(
            f"{f['_task']}:{f.get('label')}" for f in lax))
