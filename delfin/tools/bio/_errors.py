"""The single failure type of the biomolecular structure layer.

Input: nothing; this module only defines the exception the readers raise.
Output: :class:`StructureFormatError`.

Semantics: every reader in this package refuses a file it cannot read in
full rather than returning what it managed to parse.  Half a structure is
indistinguishable from a complete one once it has left the reader — the
chains are there, the residues are there, only some are missing — and every
block downstream (a binding-site box, a composition report, a receptor
preparation) would then be computed on a system nobody described.  The
message names the file, the line and what was expected there.

Rejected alternative: a warning plus a partial structure, which is what most
tolerant parsers do.  It moves the decision to a caller who no longer has the
file in hand, and the warning is lost as soon as the result is serialised.
"""

from __future__ import annotations


class StructureFormatError(ValueError):
    """A structure or ligand file could not be read as the format it claims.

    A ``ValueError`` subclass, so a caller that only distinguishes "bad input"
    from "bug" keeps working without importing this module.
    """


__all__ = ["StructureFormatError"]
