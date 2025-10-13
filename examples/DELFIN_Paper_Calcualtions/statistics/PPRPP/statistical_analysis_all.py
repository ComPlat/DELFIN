import argparse
import pandas as pd
import numpy as np
from scipy import stats

# CLI arguments
parser = argparse.ArgumentParser(description='Statistical analysis for PPRPP dataset.')
parser.add_argument('--outliers', type=int, default=5, help='Number of largest absolute errors to drop (default: 5).')
args = parser.parse_args()

num_outliers_requested = max(args.outliers, 0)

# Read PPRPP.csv data
df = pd.read_csv('PPRPP.csv', sep=';', decimal=',')

# Define the column pairs to analyze
column_pairs = [
    ('E_00_calc', 'E_00_exp', 'E_00'),
    ('E_red_calc', 'E_red_exp', 'E_red'),
    ('E_ox_calc', 'E_ox_exp', 'E_ox'),
    ('E_red*_calc', 'E_red*_exp', 'E_red*'),
    ('E_ox*_calc', 'E_ox*_exp', 'E_ox*')
]

def calculate_statistics(obs, pred, name):
    """Calculate statistical metrics"""
    if len(obs) == 0:
        return None

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

# Open output file
with open('statistic.txt', 'w') as f:
    f.write("="*100 + "\n")
    f.write("STATISTICAL ANALYSIS FOR ALL COLUMN PAIRS\n")
    f.write("="*100 + "\n\n")

    all_results = []

    for calc_col, exp_col, label in column_pairs:
        print(f"\n{'='*100}")
        print(f"Analyzing: {label}")
        print(f"{'='*100}")

        f.write(f"\n{'='*100}\n")
        f.write(f"Analyzing: {label}\n")
        f.write(f"{'='*100}\n\n")

        # Filter data where both columns have values
        df_pair = df[[calc_col, exp_col]].dropna()

        if len(df_pair) == 0:
            print(f"No data available for {label}")
            f.write(f"No data available for {label}\n")
            continue

        observed = df_pair[exp_col].values
        predicted = df_pair[calc_col].values

        # Calculate absolute errors
        abs_errors = np.abs(predicted - observed)

        # Determine how many outliers to drop while keeping at least two data points
        max_removable = max(len(df_pair) - 2, 0)
        n_outliers = min(num_outliers_requested, max_removable)

        if n_outliers:
            outlier_indices = np.argsort(abs_errors)[-n_outliers:]
        else:
            outlier_indices = np.array([], dtype=int)

        # Show outliers
        print(f"\nRemoved outliers ({n_outliers} largest absolute errors):")
        f.write(f"Removed outliers ({n_outliers} largest absolute errors):\n")
        f.write("="*50 + "\n")

        for i, idx in enumerate(outlier_indices):
            calc_val = predicted[idx]
            exp_val = observed[idx]
            error = abs_errors[idx]
            print(f"{i+1}. Calc={calc_val:.3f}, Exp={exp_val:.3f}, Error={error:.3f}")
            f.write(f"{i+1}. Calc={calc_val:.3f}, Exp={exp_val:.3f}, Error={error:.3f}\n")

        # Remove the outliers
        mask = np.ones(len(df_pair), dtype=bool)
        if n_outliers:
            mask[outlier_indices] = False
        observed_adjusted = observed[mask]
        predicted_adjusted = predicted[mask]

        # Calculate statistics for both datasets
        initial_label = f"{label} - Initial Values"
        adjusted_label = f"{label} - Initial Values (-{n_outliers} outliers)" if n_outliers else initial_label

        original_stats = calculate_statistics(observed, predicted, initial_label)
        adjusted_stats = calculate_statistics(observed_adjusted, predicted_adjusted, adjusted_label)

        # Create and display results
        results = pd.DataFrame([original_stats, adjusted_stats])
        all_results.extend([original_stats, adjusted_stats])

        print("\nStatistical Analysis:")
        print("="*100)
        print(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}")
        print("-"*100)

        f.write("\n" + "="*100 + "\n")
        f.write("Statistical Analysis:\n")
        f.write("="*100 + "\n")
        f.write(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}\n")
        f.write("-"*100 + "\n")

        for _, row in results.iterrows():
            print(f"{row['Dataset']:<40} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}")
            f.write(f"{row['Dataset']:<40} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}\n")

        # Perform linear fit on adjusted data (without outliers)
        # Fix slope = 1.0, only optimize intercept
        if len(observed_adjusted) > 1:
            # Calculate intercept with fixed slope = 1.0
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

            f.write(f"\nLinear fit (adjusted data, fixed slope=1.0): y = {slope:.4f}x + {intercept:.4f}\n")
            f.write(f"R²: {r2:.4f}, R-value: {r_value:.4f}\n")
            f.write(f"Number of complete data points: {len(observed)}\n")
            f.write(f"Number after removing outliers: {len(observed_adjusted)}\n")

            # Adjust calculated values by removing intercept
            df_adjusted = df.copy()
            df_adjusted[calc_col] = np.where(
                df_adjusted[calc_col].notna(),
                (df_adjusted[calc_col] - intercept),
                np.nan
            )

            # Recalculate delta column if it exists
            delta_col = f"delta_{label}"
            if delta_col in df_adjusted.columns:
                df_adjusted[delta_col] = df_adjusted[calc_col] - df_adjusted[exp_col]

            # Store adjusted dataframe for later
            if 'df_all_adjusted' not in locals():
                df_all_adjusted = df[['Compound', 'SMILES']].copy()

            # Add adjusted columns
            df_all_adjusted[f'{label}_calc_adjusted'] = df_adjusted[calc_col]
            df_all_adjusted[f'{label}_exp'] = df[exp_col]
            df_all_adjusted[f'delta_{label}_adjusted'] = df_adjusted[delta_col] if delta_col in df_adjusted.columns else np.nan

            # Calculate statistics for adjusted values (with same outliers removed)
            df_final_complete = df_adjusted[[calc_col, exp_col]].dropna()
            observed_final_all = df_final_complete[exp_col].values
            predicted_final_all = df_final_complete[calc_col].values

            # Apply the same outlier mask to the adjusted calculated values
            observed_final = observed_final_all[mask]
            predicted_final = predicted_final_all[mask]

            final_label = f"{label} - Adjusted Calc Values (-{n_outliers} outliers)" if n_outliers else f"{label} - Adjusted Calc Values"
            final_stats = calculate_statistics(observed_final, predicted_final, final_label)

            if final_stats:
                all_results.append(final_stats)

                print("\nStatistics for Adjusted Calculated Values:")
                print(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}")
                print("-"*100)
                print(f"{final_stats['Dataset']:<40} {final_stats['N']:<5.0f} {final_stats['MD']:<8.4f} {final_stats['MAD']:<8.4f} {final_stats['SD']:<8.4f} {final_stats['AMAX']:<8.4f} {final_stats['R²']:<8.4f} {final_stats['Adj. R²']:<8.4f}")

                f.write("\nStatistics for Adjusted Calculated Values:\n")
                f.write(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}\n")
                f.write("-"*100 + "\n")
                f.write(f"{final_stats['Dataset']:<40} {final_stats['N']:<5.0f} {final_stats['MD']:<8.4f} {final_stats['MAD']:<8.4f} {final_stats['SD']:<8.4f} {final_stats['AMAX']:<8.4f} {final_stats['R²']:<8.4f} {final_stats['Adj. R²']:<8.4f}\n")

    # Summary table
    print("\n" + "="*100)
    print("SUMMARY OF ALL ANALYSES")
    print("="*100)

    f.write("\n\n" + "="*100 + "\n")
    f.write("SUMMARY OF ALL ANALYSES\n")
    f.write("="*100 + "\n")

    summary_df = pd.DataFrame(all_results)
    print(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}")
    print("-"*100)

    f.write(f"{'Dataset':<40} {'N':<5} {'MD':<8} {'MAD':<8} {'SD':<8} {'AMAX':<8} {'R²':<8} {'Adj. R²':<8}\n")
    f.write("-"*100 + "\n")

    for _, row in summary_df.iterrows():
        print(f"{row['Dataset']:<40} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}")
        f.write(f"{row['Dataset']:<40} {row['N']:<5.0f} {row['MD']:<8.4f} {row['MAD']:<8.4f} {row['SD']:<8.4f} {row['AMAX']:<8.4f} {row['R²']:<8.4f} {row['Adj. R²']:<8.4f}\n")

# Save adjusted CSV if it was created
if 'df_all_adjusted' in locals():
    df_all_adjusted.to_csv('PPRPP_adjusted.csv', sep=';', decimal=',', index=False)
    print(f"\nAdjusted data saved to 'PPRPP_adjusted.csv'")

print(f"Complete results saved to 'statistic.txt'")
