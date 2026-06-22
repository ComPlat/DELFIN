import argparse
import pandas as pd
import numpy as np
from scipy import stats

# CLI arguments
parser = argparse.ArgumentParser(description='Statistical analysis for OROP data (combined E_ox/E_red).')
parser.add_argument('--outliers', type=int, default=5, help='Number of largest absolute errors to drop (default: 5).')
args = parser.parse_args()

num_outliers_requested = max(args.outliers, 0)

# Read OROP.csv data
df = pd.read_csv('OROP.csv', sep=';', decimal=',')

# Create combined dataset: use E_ox if available, otherwise E_red
data_list = []
for idx, row in df.iterrows():
    compound = row['Compound']

    # Priority: E_ox > E_red
    if pd.notna(row['E_ox_exp']) and pd.notna(row['E_ox_calc']):
        data_list.append({
            'Compound': compound,
            'observed': row['E_ox_exp'],
            'predicted': row['E_ox_calc'],
            'type': 'E_ox'
        })
    elif pd.notna(row['E_red_exp']) and pd.notna(row['E_red_calc']):
        data_list.append({
            'Compound': compound,
            'observed': row['E_red_exp'],
            'predicted': row['E_red_calc'],
            'type': 'E_red'
        })

df_complete = pd.DataFrame(data_list)

print(f"Total data points collected:")
print(f"  E_ox: {(df_complete['type'] == 'E_ox').sum()}")
print(f"  E_red: {(df_complete['type'] == 'E_red').sum()}")
print(f"  Total: {len(df_complete)}\n")

# Extract observed and predicted values
observed = df_complete['observed'].values
predicted = df_complete['predicted'].values
compounds = df_complete['Compound'].values
types = df_complete['type'].values

# Calculate absolute errors
abs_errors = np.abs(predicted - observed)

# Determine number of outliers to remove (keep at least two data points)
max_removable = max(len(df_complete) - 2, 0)
num_outliers = min(num_outliers_requested, max_removable)

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
print(f"Removed outliers ({num_outliers} largest absolute errors):")
print("="*70)
for i, idx in enumerate(outlier_indices):
    calc_val = predicted[idx]
    exp_val = observed[idx]
    error = abs_errors[idx]
    comp = compounds[idx]
    typ = types[idx]
    print(f"{i+1:2d}. Compound {comp:3} ({typ:5}) | Calc={calc_val:6.3f} | Exp={exp_val:6.3f} | Error={error:.3f}")

print("\n" + "="*90)

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
    f.write(f"Total data points collected:\n")
    f.write(f"  E_ox: {(df_complete['type'] == 'E_ox').sum()}\n")
    f.write(f"  E_red: {(df_complete['type'] == 'E_red').sum()}\n")
    f.write(f"  Total: {len(df_complete)}\n\n")

    f.write(f"Removed outliers ({num_outliers} largest absolute errors):\n")
    f.write("="*70 + "\n")
    for i, idx in enumerate(outlier_indices):
        calc_val = predicted[idx]
        exp_val = observed[idx]
        error = abs_errors[idx]
        comp = compounds[idx]
        typ = types[idx]
        f.write(f"{i+1:2d}. Compound {comp:3} ({typ:5}) | Calc={calc_val:6.3f} | Exp={exp_val:6.3f} | Error={error:.3f}\n")

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
df_adjusted = df.copy()
df_adjusted['E_ox_calc_adjusted'] = np.where(
    df_adjusted['E_ox_calc'].notna(),
    (df_adjusted['E_ox_calc'] - intercept),
    np.nan
)
df_adjusted['E_red_calc_adjusted'] = np.where(
    df_adjusted['E_red_calc'].notna(),
    (df_adjusted['E_red_calc'] - intercept),
    np.nan
)

# Round numeric columns for tidy CSV output
numeric_cols = [
    'E_ox_calc', 'E_ox_exp', 'delta_E_ox',
    'E_red_calc', 'E_red_exp', 'delta_E_red',
    'E_ox_calc_adjusted', 'E_red_calc_adjusted'
]
present_cols = [col for col in numeric_cols if col in df_adjusted.columns]
df_adjusted[present_cols] = df_adjusted[present_cols].apply(lambda series: series.round(3))

# Save adjusted data to CSV with same format as OROP.csv
df_adjusted.to_csv('OROP_adjusted.csv', sep=';', decimal=',', index=False)

# Calculate statistics for the final adjusted values (also remove outliers)
# Rebuild adjusted data list
data_list_adjusted = []
for idx, row in df_adjusted.iterrows():
    compound = row['Compound']

    if pd.notna(row['E_ox_exp']) and pd.notna(row['E_ox_calc_adjusted']):
        data_list_adjusted.append({
            'observed': row['E_ox_exp'],
            'predicted': row['E_ox_calc_adjusted'],
            'type': 'E_ox'
        })
    elif pd.notna(row['E_red_exp']) and pd.notna(row['E_red_calc_adjusted']):
        data_list_adjusted.append({
            'observed': row['E_red_exp'],
            'predicted': row['E_red_calc_adjusted'],
            'type': 'E_red'
        })

df_final_complete = pd.DataFrame(data_list_adjusted)
observed_final_all = df_final_complete['observed'].values
predicted_final_all = df_final_complete['predicted'].values

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
print(f"Adjusted data saved to 'OROP_adjusted.csv'")
print(f"Number of complete data points: {len(observed)}")
print(f"  E_ox: {(df_complete['type'] == 'E_ox').sum()}")
print(f"  E_red: {(df_complete['type'] == 'E_red').sum()}")
print(f"Number after removing outliers: {len(observed_adjusted)}")
print(f"Linear fit intercept removed from calculated values: {intercept:.4f}")
