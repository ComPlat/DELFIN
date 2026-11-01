"""Per-atom basis/ECP and implicit-solvation builders on OPI.

The functions here return OPI block objects (``opi.input.blocks``) and simple
keywords that the other OPI builders accept as ``extra_blocks`` / extra
keywords — nothing more is needed to wire them into an ORCA input built with
``opi.input``.

Every rule enforced here traces to the ORCA 6.1 manual
(``essentialelements/basisset.txt``, ``essentialelements/solvationmodels.txt``,
``modelchemistries/semiempirical.txt``) or to
``delfin-knowledge/audits/MISTAKES.md``; a violated rule is an error returned
to the caller with the rule named — never a silently dropped setting, and
never an invented keyword.

Rules (each with its source):

1. NewGTO never loads an ECP (MISTAKES.md wave-5: "NewGTO for an element
   never loads that basis's ECP: NewGTO Pd "def2-TZVP"/"LANL2DZ" under an
   all-electron global basis runs WITHOUT a core potential").  An element that
   needs a core potential gets NewECP explicitly or a basis family that
   brings it.  Manual, essentialelements/basisset.txt ~5791: "The ECP keyword
   itself assigns only the effective core potential, not a valence basis
   set!"
2. def2 ECPs start at Rb: "The Karlsruhe def2 family are all-electron for
   elements H-Kr, and automatically load Stuttgart-Dresden ECPs for elements
   Rb-Rn" (manual ~888, Table 2.14).
3. SARC is all-electron for the relativistic (DKH/ZORA) Hamiltonians and
   needs matching SARC/J: "SARC/J is the general-purpose Coulomb-fitting
   auxiliary for all SARC orbital basis sets" (manual, note under Table 2.28);
   SARC element ranges Rb-Rn, Ac-Lr (TZVP variants), Hf-Hg (SVP variants).
4. ECP and DKH/ZORA must not be mixed for the same element.
5. Ghost atoms carry no ECP (MISTAKES.md: "Ghost atoms for counterpoise carry
   no ECP"); manual ~5963: ECPs are not added to ghost atoms by default, and
   this builder never switches AllowGhostECP on.
6. An unknown keyword or block option is an error returned to the caller,
   never dropped.

Defaults the builders would pick are returned as ``suggestions`` (with the
reason and the manual section); the caller's explicit choice always wins and
is only validated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence

__all__ = (
    "Suggestion",
    "BasisSolvError",
    "build_basis_blocks",
    "build_solvation",
    "suggest_basis",
    "suggest_solvent",
    "assemble_opi_input",
)


@dataclass(frozen=True)
class Suggestion:
    """A default the builder would pick — a suggestion, not a decision."""

    kind: str          # e.g. "basis", "ecp", "aux", "solvent"
    value: str
    reason: str
    manual_section: str


class BasisSolvError(ValueError):
    """A rule violation the caller must resolve; the rule is named."""


# ---------------------------------------------------------------------------
# Element / family tables (manual essentialelements/basisset.txt)
# ---------------------------------------------------------------------------

# Z of Rb (first element with def2 ECPs) and Kr (last all-electron def2).
_Z_RB = 37
_Z_KR = 36

# SARC element ranges, manual Table 2.28 (Relativistic SARC basis sets):
#   SARC-DKH-SVP / SARC-ZORA-SVP: Hf-Hg
#   SARC-*-TZVP / SARC-*-TZVPP:   Rb-Rn, Ac-Lr
_SARC_SVP_RANGE = frozenset(range(72, 81))          # Hf (72) .. Hg (80)
_SARC_TZVP_RANGE = (frozenset(range(37, 87)),)      # Rb (37) .. Rn (86)


def _element_z(element: str) -> int:
    """Atomic number of ``element`` (raises for an unknown symbol)."""
    from opi.utils.element import Element

    return Element(element).atomic_number


def _is_def2_family(basis: str) -> bool:
    """Whether a basis name belongs to the Karlsruhe def2/def families
    (all-electron H-Kr, automatic def2-ECP/def-ECP Rb onwards)."""
    b = basis.strip().lower()
    return b.startswith("def2") or b.startswith("def-") or b == "def"


def _is_recontracted_rel_family(basis: str, rh: str) -> bool:
    """Whether ``basis`` is the recontracted all-electron variant for the
    relativistic Hamiltonian ``rh`` (DKH-/ZORA-def2*, manual Table 2.27,
    elements H-Kr; SARC/SARC2 are all-electron for DKH/ZORA by design)."""
    b = basis.strip().lower()
    r = rh.strip().lower()
    if not r:
        return False
    return (
        b.startswith(f"{r}-def2")
        or b.startswith("sarc-")
        or b.startswith("sarc2-")
        or b.startswith("x2c-")
    )


def _is_ecp_name(name: str) -> bool:
    """Whether ``name`` is an ECP ORCA knows (an OPI Ecp member, case- and
    separator-blind, or a keyword in the knowledge store)."""
    from opi.input.simple_keywords.ecp import Ecp

    n = str(name).strip().lower()
    members = {
        str(v).lower() for v in vars(Ecp).values()
        if type(v).__name__ == "SimpleKeyword"
    }
    if n in members:
        return True
    if n.replace("-", "").replace("/", "") in {m.replace("-", "").replace("/", "") for m in members}:
        return True
    from delfin.orca_knowledge import known_keyword

    return known_keyword(name)


def _is_basis_name(name: str) -> bool:
    """Whether ``name`` is an orbital/auxiliary basis ORCA knows (an OPI
    BasisSet member, case-blind, or a keyword in the knowledge store)."""
    from opi.input.simple_keywords.basis_set import BasisSet

    n = str(name).strip().lower()
    members = {
        str(v).lower() for v in vars(BasisSet).values()
        if type(v).__name__ == "SimpleKeyword"
    }
    if n in members:
        return True
    from delfin.orca_knowledge import known_keyword

    return known_keyword(name)


def _is_aux_name(name: str) -> bool:
    """Auxiliary basis names live in BasisSet (def2/J, SARC/J, ...) and, for
    the /C families, AuxBasisSet."""
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.aux_basis_set import AuxBasisSet

    n = str(name).strip().lower()
    for box in (BasisSet, AuxBasisSet):
        members = {
            str(v).lower() for v in vars(box).values()
            if type(v).__name__ == "SimpleKeyword"
        }
        if n in members:
            return True
    from delfin.orca_knowledge import known_keyword

    return known_keyword(name)


def _basis_brings_ecp(basis: str, element_z: int) -> Optional[str]:
    """The ECP a basis family loads automatically for ``element_z``, or None.

    Manual ~888 / Table 2.14: def2 and older def are all-electron H-Kr and
    load def2-ECP / def-ECP from Rb on.  LANL2DZ uses HayWadt ECPs from Na on
    (manual ~5842).  SARC/SARC2 and the recontracted families are
    all-electron by design.
    """
    b = basis.strip().lower()
    if _is_def2_family(b):
        return "def2-ECP" if element_z >= _Z_RB else None
    if b.startswith("lanl2") or b.startswith("lanl1") or b == "sdd":
        return "HayWadt/SDD" if element_z >= 11 else None
    return None


def _basis_element_range_error(basis: str, element: str) -> Optional[str]:
    """Why ``basis`` cannot be used for ``element``, from the manual's
    element ranges (Table 2.28 SARC, Table 2.27 recontracted def2)."""
    b = basis.strip().lower()
    z = _element_z(element)
    if b.startswith(("sarc-", "sarc2-")):
        if b.startswith(("sarc-dkh-svp", "sarc-zora-svp")):
            if z not in _SARC_SVP_RANGE:
                return (
                    f"basis {basis} covers Hf-Hg only (manual Table 2.28), "
                    f"not {element} (Z={z})"
                )
        else:
            if z < 37 or z > 86 and z < 89:
                # Rb-Rn (37-86) and Ac-Lr (89-103)
                if z > 86 and z < 89:
                    pass
                elif z < 37 or z > 86 and not (89 <= z <= 103):
                    return (
                        f"basis {basis} covers Rb-Rn and Ac-Lr (manual Table "
                        f"2.28), not {element} (Z={z})"
                    )
        return None
    return None


# ---------------------------------------------------------------------------
# Basis / ECP builder
# ---------------------------------------------------------------------------


def build_basis_blocks(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]] = (),
    basis: str,
    per_element: Optional[Mapping[str, str]] = None,
    per_element_ecp: Optional[Mapping[str, str]] = None,
    ecp: str = "",
    aux_j: str = "",
    aux_jk: str = "",
    aux_c: str = "",
    aux: str = "",
    ri: str = "",
    rel_hamiltonian: str = "",
    ghost_atom_indices: Sequence[int] = (),
    suggestions: bool = True,
) -> tuple[list[Any], list[str], list[str], list[Suggestion]]:
    """Build the ``%basis`` OPI block (and any related suggestions).

    Returns ``(blocks, errors, warnings, suggestions)``.  ``blocks`` is a
    list of OPI block objects (``BlockBasis``) suitable as ``extra_blocks``;
    it is empty when validation fails.  Every rule violation is an error that
    names the rule; unknown keywords/options are errors, never dropped.
    """
    errors: list[str] = []
    warnings: list[str] = []
    sugg: list[Suggestion] = []
    per_element = dict(per_element or {})
    per_element_ecp = dict(per_element_ecp or {})

    # -- name validation ----------------------------------------------------
    if basis and not _is_basis_name(basis):
        errors.append(
            f"basis {basis!r} is not a basis ORCA knows (OPI BasisSet members "
            "and the knowledge store); refused."
        )
    for el, b in per_element.items():
        if not _is_basis_name(b):
            errors.append(
                f"per-element basis {b!r} for {el} is not a basis ORCA knows; "
                "refused."
            )
    if ecp and not _is_ecp_name(ecp):
        errors.append(f"ECP {ecp!r} is not an ECP ORCA knows; refused.")
    for el, e in per_element_ecp.items():
        if not _is_ecp_name(e):
            errors.append(f"NewECP {e!r} for {el} is not an ECP ORCA knows; refused.")
    for kind, name in (
        ("AuxJ", aux_j), ("AuxJK", aux_jk), ("AuxC", aux_c), ("Aux", aux),
    ):
        if name and not _is_aux_name(name):
            errors.append(f"{kind} {name!r} is not an auxiliary basis ORCA knows; refused.")

    if errors:
        return [], errors, warnings, sugg

    # -- element range checks (manual tables) --------------------------------
    for el, b in per_element.items():
        why = _basis_element_range_error(b, el)
        if why:
            errors.append(why)
    if basis:
        # The global basis must cover every element that has no override.
        for el in elements:
            if el in per_element:
                continue
            z = _element_z(el)
            b = basis.strip().lower()
            if b.startswith(("sarc-", "sarc2-")):
                # SARC(-2) is a heavy-element basis: light elements must come
                # from a recontracted family (manual 2.7.3, Table 2.28).
                why = _basis_element_range_error(basis, el)
                if why and z <= _Z_KR:
                    continue  # light elements pair with recontracted def2
                if why:
                    errors.append(why)
    if errors:
        return [], errors, warnings, sugg

    # -- rule 1: NewGTO never loads an ECP -----------------------------------
    for el, b in per_element.items():
        z = _element_z(el)
        all_electron_ok = (
            bool(rel_hamiltonian) and _is_recontracted_rel_family(b, rel_hamiltonian)
        )
        ecps_provided = el in per_element_ecp or bool(ecp)
        if z >= _Z_RB and not all_electron_ok and not ecps_provided:
            # MISTAKES.md wave-5: "NewGTO for an element never loads that
            # basis's ECP: NewGTO Pd "def2-TZVP"/"LANL2DZ" under an
            # all-electron global basis runs WITHOUT a core potential."
            # An element Z>=Rb therefore needs NewECP explicitly (or a
            # recontracted all-electron family under DKH/ZORA).
            errors.append(
                f"rule 'NewGTO never loads an ECP': per-element NewGTO {el} "
                f"{b!r} would run Z={z} WITHOUT a core potential (NewGTO "
                "loads only a valence basis, manual "
                "essentialelements/basisset.txt ~5791); add per_element_ecp "
                f"for {el} (NewECP) or, under DKH/ZORA, use the recontracted "
                "all-electron family (SARC / DKH-/ZORA-def2)."
            )
        elif all_electron_ok:
            warnings.append(
                f"{el} treated all-electron with {rel_hamiltonian} "
                f"({b!r}); ensure no ECP is set for it."
            )

    # -- rule 4: ECP and DKH/ZORA must not be mixed for the same element ----
    rh = rel_hamiltonian.strip().upper()
    if rh in ("DKH", "ZORA", "X2C"):
        ecp_elements = set(per_element_ecp)
        if ecp:
            z_map = {el: _element_z(el) for el in elements}
            ecp_elements |= {
                el for el, z in z_map.items()
                if _basis_brings_ecp(basis, z)
            }
        for el in sorted(ecp_elements):
            errors.append(
                f"rule 'ECP and {rh} must not be mixed': element {el} has an "
                f"ECP ({per_element_ecp.get(el, ecp)}) while the "
                f"{rh} Hamiltonian is requested; {rh} is all-electron "
                "(manual 2.7.3: scalar-relativistic calculations require "
                "specifically optimized or recontracted all-electron basis "
                "sets). Remove the ECP or the Hamiltonian."
            )

    # -- rule 3: SARC needs SARC/J -------------------------------------------
    sarc_used = any(
        b.strip().lower().startswith(("sarc-", "sarc2-"))
        for b in [basis, *per_element.values()]
    )
    if sarc_used and ri and aux_j and aux_j.strip().lower() not in (
        "sarc/j", "autoaux",
    ):
        errors.append(
            f"rule 'SARC needs SARC/J': AuxJ {aux_j!r} cannot be combined "
            "with a SARC/SARC2 orbital basis; SARC/J is the "
            "general-purpose Coulomb-fitting auxiliary for all SARC orbital "
            "basis sets (manual, note under Table 2.28). AutoAux is also "
            "accepted."
        )
    # Light elements under a SARC global basis need the recontracted family:
    if sarc_used and basis and basis.strip().lower().startswith(("sarc-", "sarc2-")):
        light = [el for el in elements if _element_z(el) <= _Z_KR]
        if light:
            warnings.append(
                "SARC covers the heavy elements only; pair it with the "
                "recontracted def2 family (e.g. ZORA-def2-TZVP) for "
                + ", ".join(sorted(set(light)))
                + " (manual 2.7.3, Table 2.27/2.28)."
            )

    # -- rule 5: def2/J assumes ECPs beyond Kr -------------------------------
    if aux_j and aux_j.strip().lower().startswith("def2/j") and rh in (
        "DKH", "ZORA", "X2C",
    ):
        errors.append(
            "rule 'def2/J assumes ECPs beyond Kr': def2/J must not be used "
            f"with {rh} (manual Table 2.31: 'Assumes the use of ECPs beyond "
            "Kr (do not use with DKH/ZORA/X2C)'). Use SARC/J or AutoAux."
        )

    # -- rule 6: ghost atoms carry no ECP ------------------------------------
    ghosts = set(ghost_atom_indices)
    if ghosts:
        for idx in ghosts:
            if not 0 <= idx < len(elements):
                errors.append(f"ghost atom index {idx} is out of range.")
        if per_element_ecp:
            ghost_elements = {elements[i] for i in ghosts if 0 <= i < len(elements)}
            overlap = ghost_elements & set(per_element_ecp)
            if overlap:
                errors.append(
                    "rule 'ghost atoms carry no ECP': a ghost atom of "
                    f"{', '.join(sorted(overlap))} would get an ECP; ghosts "
                    "carry no ECP (MISTAKES.md; manual 2.7.5.2: ECPs are not "
                    "added to ghost atoms by default) — remove the NewECP or "
                    "drop the ghost."
                )

    if errors:
        return [], errors, warnings, sugg

    # -- build the OPI block --------------------------------------------------
    from opi.input.blocks.block_basis import BlockBasis, NewBasis
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.utils.element import Element

    def _kw(name: str) -> SimpleKeyword:
        return SimpleKeyword(str(name))

    block_kwargs: dict[str, Any] = {}
    # Each element's NewGTO/NewECP override becomes its OWN BlockBasis: OPI's
    # BlockBasis holds a single NewBasis per field, so putting several elements
    # in one block would silently overwrite all but the last (a real bug seen
    # with {Sm: SARC2-..., H: ZORA-def2-...} where Sm's NewGTO vanished).
    element_blocks: list[Any] = []
    for el in dict.fromkeys((*per_element, *per_element_ecp)):
        eb: dict[str, Any] = {}
        if el in per_element:
            eb["newgto"] = NewBasis(element=Element(el), basis=_kw(per_element[el]))
        if el in per_element_ecp:
            eb["newecp"] = NewBasis(
                element=Element(el), basis=_kw(per_element_ecp[el])
            )
        element_blocks.append(BlockBasis(**eb))
    # Shared (non-per-element) keywords go in one block of their own.
    if ecp:
        block_kwargs["ecp"] = _kw(ecp)
    if aux_j:
        block_kwargs["auxj"] = _kw(aux_j)
    if aux_jk:
        block_kwargs["auxjk"] = _kw(aux_jk)
    if aux_c:
        block_kwargs["auxc"] = _kw(aux_c)
    if aux:
        block_kwargs["aux"] = _kw(aux)
    # Ghost atoms carry no ECP: AllowGhostECP is never switched on here
    # (rule 6).  If a caller truly needs it, they set it on the block
    # themselves and own the consequence.  A ghost-only input still emits the
    # (ECP-free) block so the caller has one object to attach.
    blocks = element_blocks
    if block_kwargs:
        blocks.append(BlockBasis(**block_kwargs))
    elif ghosts and not per_element and not per_element_ecp and not blocks:
        # ghost-only input: one (ECP-free) block to attach.
        blocks.append(BlockBasis())

    # -- suggestions (defaults are suggestions, the caller decides) ----------
    if suggestions:
        for el in dict.fromkeys(elements):
            z = _element_z(el)
            if el in per_element or el in per_element_ecp:
                continue
            brings = _basis_brings_ecp(basis, z)
            if brings:
                sugg.append(Suggestion(
                    kind="ecp",
                    value=f"NewECP {el} \"{brings}\" end",
                    reason=(
                        f"{basis} loads {brings} automatically for {el} "
                        "(manual Table 2.14: def2 are all-electron H-Kr and "
                        "load Stuttgart-Dresden ECPs for Rb-Rn); explicit is "
                        "only needed with per-element NewGTO overrides."
                    ),
                    manual_section="essentialelements/basisset.txt Table 2.14",
                ))
    return blocks, errors, warnings, sugg


def suggest_basis(element: str, rel_hamiltonian: str = "") -> Optional[Suggestion]:
    """Suggest an orbital basis family for ``element``, with the manual
    section; a suggestion the user confirms, never a silent choice."""
    z = _element_z(element)
    rh = rel_hamiltonian.strip().upper()
    if z >= _Z_RB and rh in ("DKH", "ZORA"):
        if 89 <= z <= 103:
            value = f"SARC2-{rh.lower()}-QZVP + ZORA-def2-* for light atoms"
        else:
            value = f"SARC-{rh.lower()}-TZVP (heavy) + {rh}-def2-TZVP (H-Kr)"
        return Suggestion(
            kind="basis",
            value=value,
            reason=(
                f"{element} (Z={z}) is beyond the all-electron def2 range "
                f"(H-Kr); SARC/SARC2 are the all-electron bases for the "
                f"{rh} Hamiltonian (manual 2.7.3, Tables 2.27/2.28)."
            ),
            manual_section="essentialelements/basisset.txt Tables 2.27-2.28",
        )
    if z >= _Z_RB:
        return Suggestion(
            kind="basis",
            value="def2-TZVP (loads def2-ECP from Rb on)",
            reason=(
                f"{element} (Z={z}) needs a core potential; the def2 family "
                "brings def2-ECP automatically (manual Table 2.14)."
            ),
            manual_section="essentialelements/basisset.txt Table 2.14",
        )
    return None


# ---------------------------------------------------------------------------
# Solvation builder
# ---------------------------------------------------------------------------


def _solvent_value(solvent: str) -> Optional[str]:
    """The ORCA spelling of ``solvent`` (case-blind against the solvents
    ORCA knows, from ``solvents_from_opi()``), or None."""
    from delfin.orca_knowledge import solvents_from_opi

    target = str(solvent).strip().lower()
    for name in solvents_from_opi():
        if str(name).lower() == target:
            return str(name)
    return None


def build_solvation(
    *,
    model: str,
    solvent: str,
    block_form: bool = False,
    method_family: str = "dft",
    job: str = "",
) -> tuple[list[Any], list[str], list[str], list[str], list[Suggestion]]:
    """Build the implicit-solvation part of an ORCA input on OPI.

    Returns ``(blocks, keywords, errors, warnings, suggestions)``.  CPCM and
    SMD come back as ``! CPCM(solvent)`` / ``! SMD(solvent)`` simple keywords
    (or, for SMD, as a ``%cpcm`` block with ``smd true`` / ``SMDsolvent`` when
    ``block_form``).  ALPB is refused outside xTB with the manual's sentence.
    """
    from opi.input.simple_keywords.base import SimpleKeyword

    errors: list[str] = []
    warnings: list[str] = []
    sugg: list[Suggestion] = []
    blocks: list[Any] = []
    keywords: list[str] = []

    model_u = str(model).strip().upper()
    value = _solvent_value(solvent)
    if value is None:
        errors.append(
            f"solvent {solvent!r} is not one ORCA knows "
            "(delfin.orca_knowledge.solvents_from_opi()); refused."
        )
        return [], [], errors, warnings, sugg

    if model_u == "CPCM":
        # Manual solvationmodels.txt ~279: 'The C-PCM model can be used via
        # ! CPCM(solvent)'.
        keywords.append(f"CPCM({value})")
        if str(job).lower() in ("tddft", "cis", "tdcis"):
            # Manual ~59: solvent shifts on transition energies via
            # TDDFT/CIS need the refractive index in addition.
            sugg.append(Suggestion(
                kind="cpcm-option",
                value="Refrac <n>",
                reason=(
                    "TDDFT/CIS solvent shifts need the solvent's refractive "
                    "index in addition to the dielectric constant "
                    "(manual, C-PCM applicability)."
                ),
                manual_section="essentialelements/solvationmodels.txt",
            ))
    elif model_u == "SMD":
        if block_form:
            # Manual ~2717: 'one can request the SMD model via the %cpcm
            # block by writing: smd true / SMDsolvent "solvent"'.
            from opi.input.blocks.block_cpcm import BlockCpcm
            from opi.input.simple_keywords import Solvent

            block = BlockCpcm(smd=True, smdsolvent=Solvent(value))
            for opt in ("smd", "smdsolvent"):
                from delfin.orca_knowledge import known_block_option

                if not known_block_option("cpcm", opt):
                    errors.append(
                        f"cpcm: block option {opt!r} is not one ORCA knows; "
                        "refused."
                    )
            blocks.append(block)
        else:
            # Manual ~2704: 'The SMD solvation model is invoked in the input
            # file via ! SMD(solvent)'.
            keywords.append(f"SMD({value})")
    elif model_u == "ALPB":
        if str(method_family).lower() != "xtb":
            # Manual solvationmodels.txt Table 2.56 note: 'Notice that ALPB
            # is only available in the context of XTB calculations.'
            errors.append(
                "rule 'ALPB is xTB-only': ALPB is only available in the "
                "context of XTB calculations (manual Table 2.56 note); "
                f"refused for method_family={method_family!r}."
            )
        else:
            # Manual semiempirical.txt ~990: '! ALPB(solvent) # use ALPB'.
            # OPI's SolvationModel enum has CPCM/SMD only, so the keyword is
            # built as a plain SimpleKeyword with the exact manual spelling.
            keywords.append(f"ALPB({value})")
    else:
        errors.append(
            f"solvation model {model!r} is not one of CPCM, SMD, ALPB; "
            "refused."
        )

    from delfin.orca_knowledge import known_keyword

    for kw in keywords:
        if not known_keyword(kw):
            errors.append(f"keyword {kw!r} is not one ORCA knows; refused.")
    if errors:
        return [], [], errors, warnings, sugg
    return blocks, keywords, errors, warnings, sugg


def suggest_solvent(solvent: str, model: str = "CPCM") -> Optional[Suggestion]:
    """Suggest the ORCA spelling and model call for a solvent."""
    value = _solvent_value(solvent)
    if value is None:
        return None
    return Suggestion(
        kind="solvent",
        value=f"! {model.upper()}({value})",
        reason=(
            f"ORCA's spelling of {solvent!r}; C-PCM is used via "
            "'! CPCM(solvent)' (manual solvationmodels.txt)."
        ),
        manual_section="essentialelements/solvationmodels.txt",
    )


# ---------------------------------------------------------------------------
# Assembly: build a complete OPI Input with the blocks/keywords added
# ---------------------------------------------------------------------------


def assemble_opi_input(
    work_dir: Any,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str = "",
    basis: str = "",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    basis_blocks: Sequence[Any] = (),
    extra_blocks: Sequence[Any] = (),
    solvation_keywords: Sequence[str] = (),
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
) -> Any:
    """Assemble the OPI ``Input`` with this builder's blocks/keywords added.

    The method part is validated by ``orca_sp_opi.validate_method`` (imported,
    not edited); the ``%basis``/``%cpcm``/``%xtb`` blocks built here are added
    as ``extra_blocks``.  Returns the ``Input`` object; writing and running
    stay with the caller (``opi.core.Calculator`` + ``delfin.orca.run_orca``).
    """
    from opi.input import Input
    from opi.input.blocks.block_scf import BlockScf

    # Reuse the single-point builder's validation for the method keywords.
    from delfin.tools.adapters.orca_sp_opi import validate_method

    method_errors, findings = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
    )
    if method_errors:
        raise BasisSolvError("; ".join(method_errors))

    inp = Input()
    keywords: list[str] = [
        kw for kw in (functional, basis, ri, aux_basis, dispersion) if kw
    ]
    keywords.extend(str(kw) for kw in solvation_keywords)
    keywords.extend(str(kw) for kw in (extra_keywords or ()))
    if keywords:
        inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore
    if scf_options:
        inp.add_blocks(BlockScf(**scf_options))
    for block in (*basis_blocks, *extra_blocks):
        inp.add_blocks(block)
    return inp
