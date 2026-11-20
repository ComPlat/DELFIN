"""Build, but never execute, an explicit protein-only MD plan."""
import argparse
from pathlib import Path

from delfin.tools.bio.artifacts import read_json, write_json
from delfin.tools.bio.workflows import validate_plan, workflow_template


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protein', type=Path, required=True)
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
    if (args.environment == 'tip3p_pme') != bool(args.solvation_file):
        parser.error('tip3p_pme requires --solvation-file; vacuum does not accept it')
    if args.environment == 'tip3p_pme' and args.periodic_policy != 'as_stored':
        parser.error('periodic analysis requires explicit --periodic-policy as_stored (no unwrapping)')
    plan = workflow_template('md')
    plan['inputs']['complex']['path'] = str(args.protein.resolve())
    plan['nodes'][0]['settings'] = dict(expected_version=args.openmm_version, forcefield='amber14/protein.ff14SB.xml',
                                        environment=args.environment, constraints=args.constraints, expected_charge=args.charge)
    if args.solvation_file:
        plan['nodes'][0]['settings']['solvation'] = read_json(args.solvation_file)
    plan['nodes'][1]['settings'] = dict(expected_version=args.openmm_version, platform=args.platform,
                                        temperature_k=args.temperature_k, friction_per_ps=args.friction_per_ps,
                                        timestep_ps=args.timestep_ps, seed=args.seed, constraint_tolerance=args.constraint_tolerance,
                                        steps=args.steps, report_interval=args.report_interval, minimize=args.minimize == 'yes')
    plan['nodes'][2]['settings'] = dict(expected_version=args.mdanalysis_version, atom_indices=args.atom_indices,
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
