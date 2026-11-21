"""Build, but never execute, an explicit protein or protein-ligand MD plan."""
import argparse
from pathlib import Path

from delfin.tools.bio.artifacts import read_json, write_json
from delfin.tools.bio.workflows import validate_plan, workflow_template


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protein', '--structure', type=Path, required=True, help='Prepared unboxed protein or assembled protein-LIG PDB')
    parser.add_argument('--ligand', type=Path, help='Authoritative explicit-H SDF for the mapped LIG residue')
    parser.add_argument('--charges-file', type=Path, help='Mapped partial charges, required with --ligand')
    parser.add_argument('--ligand-forcefield', type=Path, help='Local SMIRNOFF OFFXML, required with --ligand')
    parser.add_argument('--openmmforcefields-version', default='0.15.0')
    parser.add_argument('--openff-version', default='0.17.0')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--charge', type=int, required=True)
    parser.add_argument('--steps', type=int, required=True)
    parser.add_argument('--report-interval', type=int, required=True)
    parser.add_argument('--temperature-k', type=float, required=True)
    parser.add_argument('--friction-per-ps', type=float, required=True)
    parser.add_argument('--timestep-ps', type=float, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--constraint-tolerance', type=float, required=True)
    parser.add_argument('--constraints', choices=['none', 'HBonds'], required=True)
    parser.add_argument('--environment', choices=['vacuum', 'tip3p_pme'], default='vacuum')
    parser.add_argument('--solvation-file', type=Path, help='Explicit solvation policy JSON, required for tip3p_pme')
    parser.add_argument('--periodic-policy', choices=['reject', 'as_stored'], default='reject')
    parser.add_argument('--platform', choices=['Reference', 'CPU'], required=True)
    parser.add_argument('--minimize', choices=['yes', 'no'], required=True)
    parser.add_argument('--atom-indices', type=int, nargs='+', required=True)
    parser.add_argument('--alignment', choices=['none', 'fit'], required=True)
    parser.add_argument('--weighting', choices=['uniform', 'mass'], required=True)
    parser.add_argument('--discard-frames', type=int, required=True)
    parser.add_argument('--stride', type=int, required=True)
    parser.add_argument('--openmm-version', default='8.4.0')
    parser.add_argument('--mdanalysis-version', default='2.9.0')
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error('output already exists')
    if any((args.ligand, args.charges_file, args.ligand_forcefield)) and not all((args.ligand, args.charges_file, args.ligand_forcefield)):
        parser.error('--ligand, --charges-file and --ligand-forcefield must be supplied together')
    if (args.environment == 'tip3p_pme') != bool(args.solvation_file):
        parser.error('tip3p_pme requires --solvation-file; vacuum does not accept it')
    if args.environment == 'tip3p_pme' and args.periodic_policy != 'as_stored':
        parser.error('periodic analysis requires explicit --periodic-policy as_stored (no unwrapping)')
    plan = workflow_template('protein_ligand_md' if args.ligand else 'md')
    nodes = {node['id']: node for node in plan['nodes']}
    if args.ligand:
        plan['inputs']['ligand']['path'] = str(args.ligand.resolve())
        nodes['ligand_parameters']['settings'] = dict(charges_file=str(args.charges_file.resolve()),
            forcefield_file=str(args.ligand_forcefield.resolve()), expected_version=args.openmmforcefields_version,
            toolkit_version=args.openff_version, openmm_version=args.openmm_version)
    plan['inputs']['complex']['path'] = str(args.protein.resolve())
    nodes['parameters']['settings'] = dict(expected_version=args.openmm_version, forcefield='amber14/protein.ff14SB.xml',
                                        environment=args.environment, constraints=args.constraints, expected_charge=args.charge)
    if args.solvation_file:
        nodes['parameters']['settings']['solvation'] = read_json(args.solvation_file)
    nodes['simulate']['settings'] = dict(expected_version=args.openmm_version, platform=args.platform,
                                        temperature_k=args.temperature_k, friction_per_ps=args.friction_per_ps,
                                        timestep_ps=args.timestep_ps, seed=args.seed, constraint_tolerance=args.constraint_tolerance,
                                        steps=args.steps, report_interval=args.report_interval, minimize=args.minimize == 'yes')
    nodes['analyze']['settings'] = dict(expected_version=args.mdanalysis_version, atom_indices=args.atom_indices,
                                        alignment=args.alignment, weighting=args.weighting, discard_frames=args.discard_frames, stride=args.stride,
                                        periodic_policy=args.periodic_policy)
    report = validate_plan(plan)
    if not report['valid']:
        parser.error('; '.join(report['errors']))
    from delfin.tools import list_steps
    from delfin.tools.bio.execution import BINDINGS
    adapters = list_steps()
    for node in plan['nodes']:
        binding = BINDINGS[node['backend']]
        kwargs = {**node['settings'], **{parameter: '/planned-output' for parameter in binding.inputs.values()}}
        adapters[binding.adapter].validate_params(**kwargs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, {'plan': plan})
    print(f'Parameter file written: {args.output}; no engines run and no packages installed')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
