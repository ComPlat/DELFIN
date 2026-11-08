"""Lab device control behind a deterministic safety kernel.

The model proposes; :mod:`delfin.lab.safety` decides. See
``docs/design/lab-control.md`` for the architecture and the threat model.

Modules: ``manifest`` (what exists and its envelope), ``units``, ``expr``
(interlocks), ``safety`` (the kernel), ``approvals``, ``watchdog``
(dead-man), ``estop``, ``audit``, ``edge`` (coordinator-edge interface),
``workflow`` and ``record`` (runs across modules), ``tools`` (agent
tools), ``cli`` (``delfin-agent lab``), ``smoke``, ``drivers``.
"""
