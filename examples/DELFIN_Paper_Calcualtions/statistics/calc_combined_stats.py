import numpy as np
import pandas as pd
from pathlib import Path

OUTLIERS = 5
BASE = Path(__file__).resolve().parents[1]

DATASET_DIRS = {
    'QCPAS': BASE / 'statistics' / 'QCPAS',
    'OROP': BASE / 'statistics' / 'OROP',
    'OMROP': BASE / 'statistics' / 'OMROP',
    'PPRPP': BASE / 'statistics' / 'PPRPP',
}

SUMMARY_PATH = BASE / 'statistics' / 'statistic_summary_5_outliers.txt'


for path in DATASET_DIRS.values():
    path.mkdir(parents=True, exist_ok=True)


# dataset configuration for loading values
DATASETS = {
    'QCPAS': {
        'csv': DATASET_DIRS['QCPAS'] / 'QCPAS.csv',
        'calc_col': 'E_ox_calc',
        'exp_col': 'E_ox_exp',
    },
    'OROP': {
        'csv': DATASET_DIRS['OROP'] / 'OROP.csv',
        'calc_col': 'E_ox_calc',
        'exp_col': 'E_ox_exp',
        'secondary_calc_col': 'E_red_calc',
        'secondary_exp_col': 'E_red_exp',
    },
    'OMROP': {
        'csv': DATASET_DIRS['OMROP'] / 'OMROP.csv',
        'calc_col': 'E_red_calc',
        'exp_col': 'E_red_exp',
    },
    'PPRPP (E_red)': {
        'csv': DATASET_DIRS['PPRPP'] / 'PPRPP.csv',
        'calc_col': 'E_red_calc',
        'exp_col': 'E_red_exp',
    },
}


def remove_outliers(observed, predicted, num_outliers=OUTLIERS):
    abs_errors = np.abs(predicted - observed)
    n = len(abs_errors)
    max_removable = max(n - 2, 0)
    k = min(num_outliers, max_removable)
    mask = np.ones(n, dtype=bool)
    if k:
        outlier_indices = np.argsort(abs_errors)[-k:]
        mask[outlier_indices] = False
    return mask


def calc_stats(observed, predicted):
    diff = predicted - observed
    n = len(diff)

    md = float(np.mean(diff)) if n else float('nan')
    mad = float(np.mean(np.abs(diff))) if n else float('nan')
    sd = float(np.std(diff, ddof=1)) if n > 1 else float('nan')
    amax = float(np.max(np.abs(diff))) if n else float('nan')

    if n > 1:
        ss_res = float(np.sum((observed - predicted) ** 2))
        ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    else:
        r2 = float('nan')

    if n > 2 and not np.isnan(r2):
        p = 1
        adj_r2 = 1 - ((1 - r2) * (n - 1) / (n - p - 1))
    else:
        adj_r2 = float('nan')

    return {
        'N': n,
        'MD': md,
        'MAD': mad,
        'SD': sd,
        'AMAX': amax,
        'R2': r2,
        'AdjR2': adj_r2,
    }


def load_dataset(name, config):
    csv_path = config['csv']
    calc_col = config['calc_col']
    exp_col = config['exp_col']

    df = pd.read_csv(csv_path, sep=';', decimal=',')

    observed = []
    predicted = []

    for _, row in df.iterrows():
        if pd.notna(row.get(exp_col)) and pd.notna(row.get(calc_col)):
            observed.append(row[exp_col])
            predicted.append(row[calc_col])
        elif config.get('secondary_calc_col') and pd.notna(row.get(config['secondary_exp_col'])) and pd.notna(row.get(config['secondary_calc_col'])):
            observed.append(row[config['secondary_exp_col']])
            predicted.append(row[config['secondary_calc_col']])

    observed = np.array(observed, dtype=float)
    predicted = np.array(predicted, dtype=float)

    if len(observed) == 0:
        raise RuntimeError(f'No data available for {name}')

    mask = remove_outliers(observed, predicted)
    intercept = float(np.mean(predicted[mask] - observed[mask]))
    predicted_adjusted = predicted - intercept

    return observed[mask], predicted_adjusted[mask]


summary_rows = []
combined_observed = []
combined_predicted = []

for name, config in DATASETS.items():
    obs, pred = load_dataset(name, config)
    stats = calc_stats(obs, pred)
    summary_rows.append((name, stats))
    combined_observed.append(obs)
    combined_predicted.append(pred)

combined_observed = np.concatenate(combined_observed)
combined_predicted = np.concatenate(combined_predicted)
combined_stats = calc_stats(combined_observed, combined_predicted)
summary_rows.append(('Combined', combined_stats))

max_dataset_len = max(len(name) for name, _ in summary_rows)
widths = {
    'Dataset': max_dataset_len + 2,
    'N': 5,
    'MD': 10,
    'MAD': 10,
    'SD': 10,
    'AMAX': 10,
    'R2': 9,
    'AdjR2': 9,
}

lines = []
lines.append(f'Combined Statistics ({OUTLIERS} outliers removed)')
lines.append('=' * len(lines[0]))
header_fmt = (
    f"{{Dataset:<{widths['Dataset']}}}"
    f"{{N:>{widths['N']}}}"
    f"{{MD:>{widths['MD']}}}"
    f"{{MAD:>{widths['MAD']}}}"
    f"{{SD:>{widths['SD']}}}"
    f"{{AMAX:>{widths['AMAX']}}}"
    f"{{R2:>{widths['R2']}}}"
    f"{{AdjR2:>{widths['AdjR2']}}}"
)
lines.append(header_fmt.format(Dataset='Dataset', N='N', MD='MD', MAD='MAD', SD='SD', AMAX='AMAX', R2='R²', AdjR2='Adj. R²'))
lines.append('-' * len(lines[-1]))
row_fmt = (
    f"{{Dataset:<{widths['Dataset']}}}"
    f"{{N:>{widths['N']}}}"
    f"{{MD:>{widths['MD']}.4f}}"
    f"{{MAD:>{widths['MAD']}.4f}}"
    f"{{SD:>{widths['SD']}.4f}}"
    f"{{AMAX:>{widths['AMAX']}.4f}}"
    f"{{R2:>{widths['R2']}.4f}}"
    f"{{AdjR2:>{widths['AdjR2']}.4f}}"
)
for name, stats in summary_rows:
    lines.append(row_fmt.format(Dataset=name, **stats))

SUMMARY_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')

print(SUMMARY_PATH)
print('\n'.join(lines))
