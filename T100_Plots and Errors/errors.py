import pandas as pd
import numpy as np
import os

# Define file paths
file_paths = {
    "X100": "IVP2_UPWIND_T1000_X100_CFL1.dat",
    "X400": "IVP2_UPWIND_T1000_X400_CFL1.dat",
    "X2000": "IVP2_UPWIND_T1000_X2000_CFL1.dat",
    "Analytical": "Analytical_Solution.dat"
}

# Check if files exist
for key, path in file_paths.items():
    if not os.path.exists(path):
        print(f"Error: File not found -> {path}")
        exit()

# Load data
data = {key: np.loadtxt(path) for key, path in file_paths.items()}

# Extract analytical solution
x_analytical = data["Analytical"][:, 0]  
f_analytical = data["Analytical"][:, 1]  

# Define search window size for local averaging based on resolution
window_sizes = {"X100": 3, "X400": 5, "X2000": 7}  

# Initialize error storage
errors = {}

# Define shock threshold
shock_threshold = 1.1  

# Loop through numerical datasets
for key in ["X100", "X400", "X2000"]:
    x_numerical = data[key][:, 0]  
    f_numerical = data[key][:, 1]  

    # Compute weighted local mean for analytical values
    f_analytical_weighted = np.zeros_like(f_numerical)

    for i, x in enumerate(x_numerical):
        # Find indices of closest analytical points
        diff = np.abs(x_analytical - x)
        sorted_indices = np.argsort(diff)[:window_sizes[key]]  
        weights = 1 / (diff[sorted_indices] + 1e-6)  
        weights /= np.sum(weights)  
        f_analytical_weighted[i] = np.sum(weights * f_analytical[sorted_indices])  

    # Apply shock region filter (Ignore points where |f_analytical| > 2)
    valid_indices = np.abs(f_analytical_weighted) <= shock_threshold

    if np.sum(valid_indices) == 0:
        print(f"Warning: No valid points after shock exclusion for {key}. Skipping.")
        continue

    # Compute error using only valid indices
    diff = f_numerical[valid_indices] - f_analytical_weighted[valid_indices]
    N = len(diff)

    # Compute norms
    L1_error = np.mean(np.abs(diff))  
    L2_error = np.sqrt(np.mean(diff**2))  
    Linf_error = np.max(np.abs(diff))  

    # Store results
    errors[key] = [L1_error, L2_error, Linf_error]

# Convert to DataFrame
error_df = pd.DataFrame.from_dict(errors, orient="index", columns=["L1 Error", "L2 Error", "Linf Error"])

# Format in scientific notation for accuracy
error_df = error_df.applymap(lambda x: f"{x:.6e}")

# Print the table
print("\n### Error Analysis (Shock Regions Excluded) ###")
print(error_df.to_string(index=True))
