import pandas as pd
from scipy import stats
import numpy as np

def calculate_correlations(df):
    results = []
    # Group the data by the 'array' column
    grouped = df.groupby('array')
    for array_id, group in grouped:
        # Only run if there are more than 10 observations
        n_obs = len(group)
        if n_obs <= 10:
            continue
        print(f"Processing array {array_id}")        
        # Helper function to safely calculate pearsonr
        def get_pearson_metrics(col1, col2):
            # Drop NaNs for the specific pair to avoid errors
            valid_data = group[[col1, col2]].dropna()
            if len(valid_data) <= 10: # Re-check after dropping NaNs
                return np.nan, np.nan
            # Check for zero variance (constant values), which causes pearsonr to fail
            if valid_data[col1].std() == 0 or valid_data[col2].std() == 0:
                return np.nan, np.nan
            r, p_val = stats.pearsonr(valid_data[col1], valid_data[col2])
            return r**2, p_val # Return R-squared and p-value
        # Correlation 1: spacer_set_size vs max_cover
        r2_spacer, p_spacer = get_pearson_metrics('spacer_set_size', 'max_cover')
        # Correlation 2: common_spacer_set_size vs max_cover
        r2_common, p_common = get_pearson_metrics('common_spacer_set_size', 'max_cover')
        results.append({
            'array': array_id,
            'n_observations': n_obs,
            'r2_spacer_maxcover': r2_spacer,
            'p_spacer_maxcover': p_spacer,
            'r2_common_maxcover': r2_common,
            'p_common_maxcover': p_common
        })
    return pd.DataFrame(results)

if __name__ == "__main__":
    df = pd.read_csv("../../Main/Fig_2/Spacer_set_size_info.tsv",sep='\t')
    result_df = calculate_correlations(df)
    output_filename = "correlation_results.tsv"
    result_df.to_csv(output_filename, sep='\t', index=False)
    print(f"\nResults successfully saved to {output_filename}")
