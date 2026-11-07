"""Package F2 — the long-horizon benchmark.

The yardstick for every compaction change: metrics derived from a session's
own records (turn_metrics, tool_trace, compaction_log), a task family that
needs 100+ turns and several compactions, and a replay/live runner. F2 only
MEASURES compaction and long sessions; compaction itself lives elsewhere
(F3 / operator-B).
"""
from delfin.agent.longhorizon.metrics import LongHorizonResult, measure, to_json, to_text

# NOTE: the ``replay`` SUBMODULE is deliberately not re-exported under its
# own name here — ``replay.replay`` is the module, and re-exporting the
# function as ``replay`` would shadow the module for ``from
# delfin.agent.longhorizon import replay``.

__all__ = ["LongHorizonResult", "measure", "to_json", "to_text"]
