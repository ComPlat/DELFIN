import argparse
import numpy as np
import pandas as pd

# Load OMROP data focusing on E_red columns
DF_PATH = 'OMROP.csv'
OUTPUT_STATS = 'statistic.txt'
OUTPUT_ADJUSTED = 'OMROP_adjusted.csv'

# CLI arguments
parser = argparse.ArgumentParser(description='Statistical analysis for OMROP E_red values.')
parser.add_argument('--outliers', type=int, default=5, help='Number of largest absolute errors to drop (default: 5).')
args = parser.parse_args()

num_outliers_requested = max(args.outliers, 0)

# Read dataset (decimal commas)
df = pd.read_csv(DF_PATH, sep=';', decimal=',')

# Keep rows with both experimental and calculated E_red values
required_cols = ['E_red_calc', 'E_red_exp']
df_complete = df.dropna(subset=required_cols).copy()

if df_complete.empty:
    raise ValueError('No complete E_red data available in OMROP.csv.')

observed = df_complete['E_red_exp'].to_numpy()
predicted = df_complete['E_red_calc'].to_numpy()
compounds = df_complete['Compound'].to_numpy()

abs_errors = np.abs(predicted - observed)

# Remove the largest outliers while keeping at least two data points
max_removable = max(len(df_complete) - 2, 0)
num_outliers = min(num_outliers_requested, max_removable)

if num_outliers:
    outlier_indices = np.argsort(abs_errors)[-num_outliers:]
else:
    outlier_indices = np.array([], dtype=int)

mask = np.ones(len(df_complete), dtype=bool)
mask[outlier_indices] = False

observed_adjusted = observed[mask]
predicted_adjusted = predicted[mask]


def calculate_statistics(obs: np.ndarray, pred: np.ndarray, name: str) -> dict:
    """Return statistical metrics for observed vs predicted arrays."""
    n = len(obs)
    diff = pred - obs
    md = np.mean(diff)
    mad = np.mean(np.abs(diff))

    if n > 1:
        sd = np.std(diff, ddof=1)
    else:
        sd = np.nan

    amax = np.max(np.abs(diff)) if n else np.nan

    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)

    if n > 1 and ss_tot > 0:
        r2 = 1 - (ss_res / ss_tot)
    else:
        r2 = np.nan

    if n > 2 and not np.isnan(r2):
        p = 1
        adj_r2 = 1 - ((1 - r2) * (n - 1) / (n - p - 1))
    else:
        adj_r2 = np.nan

    return {
        'Dataset': name,
        'N': n,
        'MD': md,
        'MAD': mad,
        'SD': sd,
        'AMAX': amax,
        'R²': r2,
        'Adj. R²': adj_r2,
    }


def print_statistics_table(rows: list[dict]):
    header = f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}"
    line = '-' * len(header)
    print(header)
    print(line)
    for row in rows:
        print(
            f"{row['Dataset']:<35} "
            f"{row['N']:<5.0f} "
            f"{row['MD']:<8.4f} "
            f"{row['MAD']:<8.4f} "
            f"{row['SD']:<8.4f} "
            f"{row['AMAX']:<8.4f} "
            f"{row['R²']:<8.4f} "
            f"{row['Adj. R²']:<8.4f}"
        )


def write_statistics_table(file_obj, rows: list[dict]):
    header = f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}"
    line = '-' * len(header)
    file_obj.write(header + '\n')
    file_obj.write(line + '\n')
    for row in rows:
        file_obj.write(
            f"{row['Dataset']:<35} "
            f"{row['N']:<5.0f} "
            f"{row['MD']:<8.4f} "
            f"{row['MAD']:<8.4f} "
            f"{row['SD']:<8.4f} "
            f"{row['AMAX']:<8.4f} "
            f"{row['R²']:<8.4f} "
            f"{row['Adj. R²']:<8.4f}\n"
        )


print(f"Removed outliers ({num_outliers} largest absolute E_red errors):")
print('=' * 70)
for idx, df_idx in enumerate(outlier_indices, start=1):
    print(
        f"{idx:2d}. Compound {int(compounds[df_idx]):3d} | "
        f"Calc={predicted[df_idx]:6.3f} | Exp={observed[df_idx]:6.3f} | Error={abs_errors[df_idx]:.3f}"
    )

print('\n' + '=' * 90)

label_adjusted = f"Initial Values (-{num_outliers} outliers)" if num_outliers else 'Initial Values'
original_stats = calculate_statistics(observed, predicted, 'Initial Values')
adjusted_stats = calculate_statistics(observed_adjusted, predicted_adjusted, label_adjusted)

print('Statistical Analysis (E_red):')
print('=' * 90)
print_statistics_table([original_stats, adjusted_stats])

with open(OUTPUT_STATS, 'w') as f:
    f.write(f'Removed outliers ({num_outliers} largest absolute E_red errors):\n')
    f.write('=' * 70 + '\n')
    for idx, df_idx in enumerate(outlier_indices, start=1):
        f.write(
            f"{idx:2d}. Compound {int(compounds[df_idx]):3d} | "
            f"Calc={predicted[df_idx]:6.3f} | Exp={observed[df_idx]:6.3f} | Error={abs_errors[df_idx]:.3f}\n"
        )
    f.write('\n' + '=' * 90 + '\n')
    f.write('Statistical Analysis (E_red):\n')
    f.write('=' * 90 + '\n')
    write_statistics_table(f, [original_stats, adjusted_stats])
    f.write('\n')
    f.write(f'Number of complete data points: {len(observed)}\n')
    f.write(f'Number after removing outliers: {len(observed_adjusted)}\n')

if len(observed_adjusted) == 0:
    raise ValueError('No data left after outlier removal; cannot compute adjusted fit.')

intercept = np.mean(predicted_adjusted - observed_adjusted)
slope = 1.0

residuals = predicted_adjusted - (slope * observed_adjusted + intercept)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((observed_adjusted - np.mean(observed_adjusted)) ** 2)

if ss_tot > 0:
    r2_fit = 1 - (ss_res / ss_tot)
    r_value = np.sqrt(r2_fit)
else:
    r2_fit = np.nan
    r_value = np.nan

print(f"\nLinear fit (adjusted data, slope fixed at 1.0): y = {slope:.4f}x + {intercept:.4f}")
print(f"R²: {r2_fit:.4f}, R-value: {r_value:.4f}")

# Apply intercept correction to calculated E_red values

df_adjusted = df.copy()
df_adjusted['E_red_calc_adjusted'] = np.where(
    df_adjusted['E_red_calc'].notna(),
    df_adjusted['E_red_calc'] - intercept,
    np.nan,
)

df_adjusted.to_csv(OUTPUT_ADJUSTED, sep=';', decimal=',', index=False)

# Evaluate statistics for adjusted values using the same mask
final_complete = df_adjusted.dropna(subset=required_cols + ['E_red_calc_adjusted']).copy()
observed_final_all = final_complete['E_red_exp'].to_numpy()
predicted_final_all = final_complete['E_red_calc_adjusted'].to_numpy()

observed_final = observed_final_all[mask]
predicted_final = predicted_final_all[mask]

label_final = f'Adjusted Calc Values (-{num_outliers} outliers)' if num_outliers else 'Adjusted Calc Values'
final_stats = calculate_statistics(observed_final, predicted_final, label_final)

print('\n' + '=' * 90)
print('Final Statistical Analysis (E_red, adjusted values):')
print('=' * 90)
print_statistics_table([original_stats, adjusted_stats, final_stats])

with open(OUTPUT_STATS, 'a') as f:
    f.write(f"\nLinear fit (adjusted data, slope fixed at 1.0): y = {slope:.4f}x + {intercept:.4f}\n")
    f.write(f"R²: {r2_fit:.4f}, R-value: {r_value:.4f}\n")
    f.write(f"Intercept removed from calculated values: {intercept:.4f}\n")
    f.write('\n' + '=' * 90 + '\n')
    f.write('Final Statistical Analysis (E_red, adjusted values):\n')
    f.write('=' * 90 + '\n')
    write_statistics_table(f, [original_stats, adjusted_stats, final_stats])

print(f"\nResults saved to {OUTPUT_STATS!r}")
print(f"Adjusted data saved to {OUTPUT_ADJUSTED!r}")
print(f"Number of complete data points: {len(observed)}")
print(f"Number after removing outliers: {len(observed_adjusted)}")
print(f"Linear fit intercept removed from E_red_calc: {intercept:.4f}")
