"""DELFIN operations MCP server.

Exposes typed DELFIN workflow actions (pipeline run/prepare, runtime checks,
cleanup, stop, ORCA invocation, CO2/TADF/hyperpol) over the Model Context
Protocol via stdio transport.

This is the ``delfin-ops`` counterpart to ``delfin-docs``: where ``docs``
exposes literature/calc search, ``ops`` exposes safe runtime actions.
"""
