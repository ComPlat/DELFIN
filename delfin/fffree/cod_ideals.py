"""delfin.fffree.cod_ideals — COD-empirical bond-length ideals (per element pair, p50
over the disorder-cleaned COD corpus).  Construction-driver source (Iter 29): the FF-free
refiner pulls heavy-heavy bonds toward these REAL crystal lengths instead of the generic
covalent-sum _ideal_bond (which is systematically too long for aromatic/conjugated ligands).
Only well-populated, reliable pairs (n>=200, use=True) are included; others fall back to
_ideal_bond.  CCDC-ready: a CSD/Mogul provider can replace/extend this dict (see
docs/plans/CCDC_ACTIVATION.md).  Generated from bond_pair_ideals_full.json."""
_COD_P50 = {
    'Al-C': 1.988,
    'Al-N': 1.918,
    'Al-O': 1.828,
    'As-C': 1.94,
    'As-O': 1.692,
    'B-B': 1.776,
    'B-C': 1.694,
    'B-Cl': 1.796,
    'B-F': 1.386,
    'B-N': 1.542,
    'B-O': 1.429,
    'B-P': 1.93,
    'B-S': 1.906,
    'Br-C': 1.896,
    'C-C': 1.396,
    'C-Cl': 1.737,
    'C-F': 1.34,
    'C-Ga': 1.994,
    'C-I': 2.096,
    'C-N': 1.359,
    'C-O': 1.271,
    'C-P': 1.83,
    'C-S': 1.762,
    'C-Se': 1.923,
    'C-Si': 1.871,
    'C-Te': 2.131,
    'Cl-O': 1.429,
    'Cl-P': 2.03,
    'F-P': 1.545,
    'Ga-N': 1.951,
    'In-S': 2.447,
    'Li-N': 2.035,
    'Li-O': 1.929,
    'N-N': 1.359,
    'N-O': 1.241,
    'N-P': 1.652,
    'N-S': 1.618,
    'N-Si': 1.735,
    'O-O': 1.46,
    'O-P': 1.547,
    'O-S': 1.446,
    'O-Si': 1.629,
    'P-P': 2.526,
    'P-S': 2.002,
    'P-Se': 2.17,
    'S-S': 2.054,
    'Si-Si': 2.354,
}


def cod_ideal_bond(a, b):
    """COD-empirical ideal bond length for element pair (a,b); None if not covered."""
    return _COD_P50.get(a+"-"+b) or _COD_P50.get(b+"-"+a)
