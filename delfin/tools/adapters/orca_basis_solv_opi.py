"""Adapter: OPI single-point step with per-atom basis/ECP + solvation.

Wraps :mod:`delfin.tools.orca_basis_solv_opi` (the builders) and
:mod:`delfin.tools.adapters.orca_sp_opi` (the run/read path, imported, never
edited) so one step carries per-element NewGTO/NewECP overrides and an
implicit solvation model.  The ``%basis``/``%cpcm`` blocks and the
``CPCM(solvent)``/``SMD(solvent)`` keywords ride in through the OPI ``Input``
the single-point builder assembles.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.tools._base import StepAdapter
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools._registry import register

__all__ = ("OrcaBasisSolvOpiAdapter",)


def _build_input_path(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str,
    aux_basis: str,
    dispersion: str,
    extra_keywords: Optional[Sequence[str]],
    scf_options: Optional[Mapping[str, Any]],
    per_element: Optional[Mapping[str, str]],
    per_element_ecp: Optional[Mapping[str, str]],
    ecp: str,
    aux_j: str,
    aux_jk: str,
    aux_c: str,
    aux: str,
    rel_hamiltonian: str,
    ghost_atom_indices: Sequence[int],
    solvent_model: str,
    solvent: str,
    solvation_block_form: bool,
    method_family: str,
    job: str,
    maxcore: int,
    ncores: int,
    base_name: str,
) -> tuple[Path, list[str], list[str], list[Any]]:
    """Validate + assemble the OPI input with basis/solvation added.

    Returns ``(inp_path, errors, warnings, suggestions)``; ``inp_path`` is
    None-equivalent (empty string in ``errors[0]`` case) when validation
    failed — callers check ``errors`` first.
    """
    from delfin.tools.orca_basis_solv_opi import (
        build_basis_blocks,
        build_solvation,
    )

    basis_blocks, b_errors, b_warnings, b_sugg = build_basis_blocks(
        elements=elements,
        basis=basis,
        per_element=per_element,
        per_element_ecp=per_element_ecp,
        ecp=ecp,
        aux_j=aux_j,
        aux_jk=aux_jk,
        aux_c=aux_c,
        aux=aux,
        ri=ri,
        rel_hamiltonian=rel_hamiltonian,
        ghost_atom_indices=ghost_atom_indices,
    )
    solv_blocks, solv_keywords, s_errors, s_warnings, s_sugg = build_solvation(
        model=solvent_model,
        solvent=solvent,
        block_form=solvation_block_form,
        method_family=method_family,
        job=job,
    )
    errors = b_errors + s_errors
    warnings = b_warnings + s_warnings
    suggestions = [
        dict(kind=s.kind, value=s.value, reason=s.reason,
             manual_section=s.manual_section)
        for s in (*b_sugg, *s_sugg)
    ]
    if errors:
        return Path(""), errors, warnings, suggestions

    # Assemble with OPI: reuse the single-point builder's keyword set, then
    # add the solvation keywords and this module's blocks.
    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_scf import BlockScf
    from opi.input.structures import Atom, Structure

    from delfin.tools.adapters.orca_sp_opi import validate_method

    method_errors, _ = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
    )
    if method_errors:
        errors.extend(method_errors)
        return Path(""), errors, warnings, suggestions

    inp = Input()
    keywords = [
        kw for kw in (functional, basis, ri, aux_basis, dispersion) if kw
    ]
    keywords.extend(solv_keywords)
    keywords.extend(str(kw) for kw in (extra_keywords or ()))
    if keywords:
        inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore
    if scf_options:
        inp.add_blocks(BlockScf(**scf_options))
    for block in (*basis_blocks, *solv_blocks):
        inp.add_blocks(block)

    atoms = [
        Atom(str(el), coordinates=(float(x), float(y), float(z)))
        for el, (x, y, z) in zip(elements, coordinates)
    ]
    structure = Structure(atoms, charge=charge, multiplicity=mult)

    calc = Calculator(base_name, working_dir=work_dir, version_check=False)
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return work_dir / f"{base_name}.inp", errors, warnings, suggestions


class OrcaBasisSolvOpiAdapter(StepAdapter):
    """An OPI single-point step with per-atom basis/ECP and solvation."""

    name = "orca_sp_opi_basis_solv"
    description = (
        "OPI ORCA single point with per-element basis/ECP (%basis NewGTO/"
        "NewECP) and implicit solvation (CPCM/SMD, ALPB xTB-only)"
    )
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "property_json")
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements = kwargs["elements"]
            coordinates = kwargs["coordinates"]
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required {exc}",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        inp_path, errors, warnings, suggestions = _build_input_path(
            work_dir,
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "B3LYP"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            per_element=kwargs.get("per_element"),
            per_element_ecp=kwargs.get("per_element_ecp"),
            ecp=kwargs.get("ecp", ""),
            aux_j=kwargs.get("aux_j", ""),
            aux_jk=kwargs.get("aux_jk", ""),
            aux_c=kwargs.get("aux_c", ""),
            aux=kwargs.get("aux", ""),
            rel_hamiltonian=kwargs.get("rel_hamiltonian", ""),
            ghost_atom_indices=kwargs.get("ghost_atom_indices", ()),
            solvent_model=kwargs.get("solvent_model", ""),
            solvent=kwargs.get("solvent", ""),
            solvation_block_form=kwargs.get("solvation_block_form", False),
            method_family=kwargs.get("method_family", "dft"),
            job=kwargs.get("job", ""),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
        )
        if errors:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error="; ".join(errors),
                error_kind=ErrorKind.INVALID_PARAM,
            )
        if inp_path is None or not inp_path.is_file():
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error="input file was not written",
                error_kind=ErrorKind.TOOL_FAILED,
            )

        from delfin.orca import run_orca
        from delfin.tools.adapters.orca_sp_opi import read_orca_output

        out_path = work_dir / f"{Path(inp_path).stem}.out"
        success = run_orca(
            str(inp_path), str(out_path),
            timeout=kwargs.get("timeout"), working_dir=work_dir,
        )
        result = read_orca_output(work_dir, Path(inp_path).stem)
        result.findings = tuple(result.findings) + tuple(warnings)
        if not success and result.ok:
            result.ok = False
            result.meaningful = False
            result.error = result.error or "ORCA subprocess reported failure"
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if result.ok else StepStatus.FAILED,
            work_dir,
            start,
            output_file=result.output_file,
            data={
                "meaningful": result.meaningful,
                "energy_Eh": result.energy_Eh,
                "s2": result.s2,
                "hl_gap_eV": result.hl_gap_eV,
                "warnings": list(warnings),
                "suggestions": suggestions,
                "failure_lines": list(result.failure_lines),
            },
            error=result.error,
            error_kind=(ErrorKind.NONE if result.ok
                        else ErrorKind.CONVERGENCE if result.scf_converged is False
                        else ErrorKind.TOOL_FAILED),
        )


register(OrcaBasisSolvOpiAdapter())
