import argparse
import pandas as pd
import numpy as np
from scipy import stats

# CLI arguments
parser = argparse.ArgumentParser(description='Statistical analysis for QCPAS E_ox values.')
parser.add_argument('--outliers', type=int, default=5, help='Number of largest absolute errors to drop (default: 5).')
args = parser.parse_args()

num_outliers_requested = max(args.outliers, 0)

# Read QCPAS.csv data
df = pd.read_csv('QCPAS.csv', sep=';', decimal=',')

# Only include rows where BOTH columns have values (not NaN)
df_complete = df.dropna()

# Extract observed and predicted values
observed = df_complete['E_ox_exp'].values
predicted = df_complete['E_ox_calc'].values

# Calculate absolute errors
abs_errors = np.abs(predicted - observed)

# Determine how many outliers can be removed safely (keep at least two points)
max_removable = max(len(df_complete) - 2, 0)
num_outliers = min(num_outliers_requested, max_removable)

# Find indices of the largest outliers
if num_outliers:
    outlier_indices = np.argsort(abs_errors)[-num_outliers:]
else:
    outlier_indices = np.array([], dtype=int)

# Remove the largest outliers
mask = np.ones(len(df_complete), dtype=bool)
mask[outlier_indices] = False
observed_adjusted = observed[mask]
predicted_adjusted = predicted[mask]

def calculate_statistics(obs, pred, name):
    """Calculate statistical metrics"""
    # Mean Deviation (MD)
    md = np.mean(pred - obs)

    # Mean Absolute Deviation (MAD)
    mad = np.mean(np.abs(pred - obs))

    # Standard Deviation (SD)
    sd = np.std(pred - obs, ddof=1)

    # Maximum Absolute Error (AMAX)
    amax = np.max(np.abs(pred - obs))

    # R-squared
    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    # Adjusted R-squared (assuming 1 parameter)
    n = len(obs)
    p = 1
    adj_r2 = 1 - ((1 - r2) * (n - 1) / (n - p - 1))

    return {
        'Dataset': name,
        'N': n,
        'MD': md,
        'MAD': mad,
        'SD': sd,
        'AMAX': amax,
        'R²': r2,
        'Adj. R²': adj_r2
    }

# Show outliers
# Report removed outliers
print(f"Removed outliers ({num_outliers} largest absolute errors):")
for i, idx in enumerate(outlier_indices):
    calc_val = predicted[idx]
    exp_val = observed[idx]
    error = abs_errors[idx]
    print(f"{i+1}. Calc={calc_val:.3f}, Exp={exp_val:.3f}, Error={error:.3f}")

print("\n" + "="*70)

# Calculate statistics for both datasets
label_adjusted = f"Initial Values (-{num_outliers} outliers)" if num_outliers else "Initial Values"
original_stats = calculate_statistics(observed, predicted, "Initial Values")
adjusted_stats = calculate_statistics(observed_adjusted, predicted_adjusted, label_adjusted)

# Create and display results
results = pd.DataFrame([original_stats, adjusted_stats])

print("Statistical Analysis:")
print("="*90)
print(f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}")
print("-"*90)

for _, row in results.iterrows():
    print(f"{row['Dataset']:<35} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}")

# Write detailed results to statistic.txt
with open('statistic.txt', 'w') as f:
    f.write(f"Removed outliers ({num_outliers} largest absolute errors):\n")
    f.write("="*50 + "\n")
    for i, idx in enumerate(outlier_indices):
        calc_val = predicted[idx]
        exp_val = observed[idx]
        error = abs_errors[idx]
        f.write(f"{i+1}. Calc={calc_val:.3f}, Exp={exp_val:.3f}, Error={error:.3f}\n")

    f.write("\n" + "="*90 + "\n")
    f.write("Statistical Analysis:\n")
    f.write("="*90 + "\n")
    f.write(f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}\n")
    f.write("-"*90 + "\n")

    for _, row in results.iterrows():
        f.write(f"{row['Dataset']:<35} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}\n")

    f.write(f"\nNumber of complete data points: {len(observed)}\n")
    f.write(f"Number after removing outliers: {len(observed_adjusted)}\n")

# Perform linear fit on adjusted data (without outliers)
# Fix slope = 1.0, only optimize intercept
# y = 1*x + intercept  =>  intercept = mean(y - x)
intercept = np.mean(predicted_adjusted - observed_adjusted)
slope = 1.0

# Calculate R-value for fixed slope fit
residuals = predicted_adjusted - (slope * observed_adjusted + intercept)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((observed_adjusted - np.mean(observed_adjusted)) ** 2)
r2 = 1 - (ss_res / ss_tot)
r_value = np.sqrt(r2) if r2 >= 0 else 0

print(f"\nLinear fit (adjusted data, fixed slope=1.0): y = {slope:.4f}x + {intercept:.4f}")
print(f"R²: {r2:.4f}, R-value: {r_value:.4f}")

# Adjust all calculated values by removing the intercept
# New formula: E_ox_calc_adjusted = (E_ox_calc - intercept) / slope
df_adjusted = df.copy()
df_adjusted['E_ox_calc_adjusted'] = np.where(
    df_adjusted['E_ox_calc'].notna(),
    (df_adjusted['E_ox_calc'] - intercept),
    np.nan
)

# Save adjusted data to CSV with same format as QCPAS.csv
df_adjusted.to_csv('QCPAS_adjusted.csv', sep=';', decimal=',', index=False)

# Calculate statistics for the final adjusted values (also remove outliers)
df_final_complete = df_adjusted.dropna()
observed_final_all = df_final_complete['E_ox_exp'].values
predicted_final_all = df_final_complete['E_ox_calc_adjusted'].values

# Apply the same outlier mask to the adjusted calculated values
observed_final = observed_final_all[mask]
predicted_final = predicted_final_all[mask]

label_final = f"Adjusted Calc Values (-{num_outliers} outliers)" if num_outliers else "Adjusted Calc Values"
final_stats = calculate_statistics(observed_final, predicted_final, label_final)

# Update results DataFrame to include final adjusted values
results_final = pd.DataFrame([original_stats, adjusted_stats, final_stats])

print("\n" + "="*90)
print("Final Statistical Analysis (including adjusted calculated values):")
print("="*90)
print(f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}")
print("-"*90)

for _, row in results_final.iterrows():
    print(f"{row['Dataset']:<35} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}")

# Update statistic.txt file with final results
with open('statistic.txt', 'a') as f:
    f.write(f"\nLinear fit (adjusted data, fixed slope=1.0): y = {slope:.4f}x + {intercept:.4f}\n")
    f.write(f"R²: {r2:.4f}, R-value: {r_value:.4f}\n")
    f.write(f"Intercept removed from calculated values: {intercept:.4f}\n")
    f.write("\n" + "="*90 + "\n")
    f.write("Final Statistical Analysis (including adjusted calculated values):\n")
    f.write("="*90 + "\n")
    f.write(f"{'Dataset':<35} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}\n")
    f.write("-"*90 + "\n")

    for _, row in results_final.iterrows():
        f.write(f"{row['Dataset']:<35} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}\n")

print(f"\nResults saved to 'statistic.txt'")
print(f"Adjusted data saved to 'QCPAS_adjusted.csv'")
print(f"Number of complete data points: {len(observed)}")
print(f"Number after removing outliers: {len(observed_adjusted)}")
print(f"Linear fit intercept removed from calculated values: {intercept:.4f}")
