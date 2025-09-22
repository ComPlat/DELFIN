# cli_calculations.py
# Calculation methods (M1, M2, M3) for redox potential calculations

import re
from .cli_helpers import _avg_or_none


def calculate_redox_potentials(config, free_gibbs_energies, E_ref):
    """Calculate redox potentials using Methods 1, 2, and 3.

    Args:
        config: Configuration dictionary
        free_gibbs_energies: Dictionary of Gibbs free energies for different charge states
        E_ref: Reference electrode potential

    Returns:
        tuple: (m1_avg, m2_step, m3_mix, use_flags) dictionaries and selection flags
    """
    # Extract individual energies
    free_gibbs_0 = free_gibbs_energies.get('0')
    free_gibbs_plus_1 = free_gibbs_energies.get('+1')
    free_gibbs_plus_2 = free_gibbs_energies.get('+2')
    free_gibbs_plus_3 = free_gibbs_energies.get('+3')
    free_gibbs_minus_1 = free_gibbs_energies.get('-1')
    free_gibbs_minus_2 = free_gibbs_energies.get('-2')
    free_gibbs_minus_3 = free_gibbs_energies.get('-3')

    # --- read selection (default -> '2'), no new parser module needed ---
    _sel_raw = config.get('calc_potential_method', config.get('calc_method', 2))
    if isinstance(_sel_raw, (list, tuple, set)):
        _sel_tokens = [str(x) for x in _sel_raw]
    else:
        _sel_tokens = re.split(r'[\s,]+', str(_sel_raw).strip())

    use_m1 = '1' in _sel_tokens
    use_m2 = '2' in _sel_tokens
    use_m3 = '3' in _sel_tokens
    if not (use_m1 or use_m2 or use_m3):
        use_m1 = True  # default

    # ---------- constants ----------
    conv = 2625.499639479947            # kJ/mol per Hartree
    F    = 96.4853321233100184          # kJ/(V·mol)

    # ---------- containers ----------
    m1_avg  = {}  # Method 1: averages (multi-e as average)
    m2_step = {}  # Method 2: step-wise (1e steps)
    m3_mix  = {}  # Method 3: (M1 + M2)/2 as requested

    # ---------- METHOD 1 (averages) ----------
    # Oxidations
    if free_gibbs_0 is not None and free_gibbs_plus_1 is not None:
        m1_avg['E_ox']   = -((free_gibbs_0 - free_gibbs_plus_1) * conv) / F - E_ref
    if free_gibbs_0 is not None and free_gibbs_plus_2 is not None:
        m1_avg['E_ox_2'] = -((free_gibbs_0 - free_gibbs_plus_2) * conv) / (2 * F) - E_ref
    if free_gibbs_0 is not None and free_gibbs_plus_3 is not None:
        m1_avg['E_ox_3'] = -((free_gibbs_0 - free_gibbs_plus_3) * conv) / (3 * F) - E_ref
    # Reductions
    if free_gibbs_0 is not None and free_gibbs_minus_1 is not None:
        m1_avg['E_red']   = ((free_gibbs_0 - free_gibbs_minus_1) * conv) / F - E_ref
    if free_gibbs_0 is not None and free_gibbs_minus_2 is not None:
        m1_avg['E_red_2'] = ((free_gibbs_0 - free_gibbs_minus_2) * conv) / (2 * F) - E_ref
    if free_gibbs_0 is not None and free_gibbs_minus_3 is not None:
        m1_avg['E_red_3'] = ((free_gibbs_0 - free_gibbs_minus_3) * conv) / (3 * F) - E_ref

    # ---------- METHOD 2 (step-wise) ----------
    # Oxidations
    if free_gibbs_0 is not None and free_gibbs_plus_1 is not None:
        m2_step['E_ox']   = -((free_gibbs_0 - free_gibbs_plus_1) * conv) / F - E_ref
    if free_gibbs_plus_1 is not None and free_gibbs_plus_2 is not None:
        m2_step['E_ox_2'] = -((free_gibbs_plus_1 - free_gibbs_plus_2) * conv) / F - E_ref
    if free_gibbs_plus_2 is not None and free_gibbs_plus_3 is not None:
        m2_step['E_ox_3'] = -((free_gibbs_plus_2 - free_gibbs_plus_3) * conv) / F - E_ref
    # Reductions
    if free_gibbs_0 is not None and free_gibbs_minus_1 is not None:
        m2_step['E_red']   = ((free_gibbs_0 - free_gibbs_minus_1) * conv) / F - E_ref
    if free_gibbs_minus_1 is not None and free_gibbs_minus_2 is not None:
        m2_step['E_red_2'] = ((free_gibbs_minus_1 - free_gibbs_minus_2) * conv) / F - E_ref
    if free_gibbs_minus_2 is not None and free_gibbs_minus_3 is not None:
        m2_step['E_red_3'] = ((free_gibbs_minus_2 - free_gibbs_minus_3) * conv) / F - E_ref

    # ---------- METHOD 3 (M1+M2)/2 on the published outputs ----------
    # 1e steps: both methods define the same thing; mean equals them
    m3_mix['E_ox']   = _avg_or_none(m1_avg.get('E_ox'),   m2_step.get('E_ox'))
    m3_mix['E_red']  = _avg_or_none(m1_avg.get('E_red'),  m2_step.get('E_red'))
    # 2e/3e: mean of M1-average and M2-step
    m3_mix['E_ox_2']  = _avg_or_none(m1_avg.get('E_ox_2'),  m2_step.get('E_ox_2'))
    m3_mix['E_ox_3']  = _avg_or_none(m1_avg.get('E_ox_3'),  m2_step.get('E_ox_3'))
    m3_mix['E_red_2'] = _avg_or_none(m1_avg.get('E_red_2'), m2_step.get('E_red_2'))
    m3_mix['E_red_3'] = _avg_or_none(m1_avg.get('E_red_3'), m2_step.get('E_red_3'))

    use_flags = {'use_m1': use_m1, 'use_m2': use_m2, 'use_m3': use_m3}

    return m1_avg, m2_step, m3_mix, use_flags


def select_final_potentials(m1_avg, m2_step, m3_mix, use_flags):
    """Select final redox potentials based on priority (3 > 2 > 1).

    Args:
        m1_avg: Method 1 results
        m2_step: Method 2 results
        m3_mix: Method 3 results
        use_flags: Dictionary with use_m1, use_m2, use_m3 flags

    Returns:
        tuple: (E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3)
    """
    use_m1 = use_flags['use_m1']
    use_m2 = use_flags['use_m2']
    use_m3 = use_flags['use_m3']

    def _pick(key):
        if use_m3 and m3_mix.get(key) is not None:
            src = 'M3'; val = m3_mix[key]
        elif use_m2 and m2_step.get(key) is not None:
            src = 'M2'; val = m2_step[key]
        elif use_m1 and m1_avg.get(key) is not None:
            src = 'M1'; val = m1_avg[key]
        else:
            src = 'fallback'; val = None
        return val, src

    E_ox,   src_ox   = _pick('E_ox')
    E_ox_2, src_ox_2 = _pick('E_ox_2')
    E_ox_3, src_ox_3 = _pick('E_ox_3')
    E_red,  src_red  = _pick('E_red')
    E_red_2,src_red_2= _pick('E_red_2')
    E_red_3,src_red_3= _pick('E_red_3')

    return E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3