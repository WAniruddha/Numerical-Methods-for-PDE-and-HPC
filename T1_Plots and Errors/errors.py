import pandas as pd
import numpy as np
import os

# Define file paths
file_paths = {
    "X100": r"C:\Users\aniru\OneDrive - Cranfield University\AHW\03 Assignments\02 NMPDE-HPC\03 Solution\01 Code\02 Fortran\Final IVP2_250130\UPWIND\T1_Plots and Errors\IVP2_UPWIND_T1_X100_CFL1.dat",
    "X400": r"C:\Users\aniru\OneDrive - Cranfield University\AHW\03 Assignments\02 NMPDE-HPC\03 Solution\01 Code\02 Fortran\Final IVP2_250130\UPWIND\T1_Plots and Errors\IVP2_UPWIND_T1_X400_CFL1.dat",
    "X2000": r"C:\Users\aniru\OneDrive - Cranfield University\AHW\03 Assignments\02 NMPDE-HPC\03 Solution\01 Code\02 Fortran\Final IVP2_250130\UPWIND\T1_Plots and Errors\IVP2_UPWIND_T1_X2000_CFL1.dat",
    "Analytical": r"C:\Users\aniru\OneDrive - Cranfield University\AHW\03 Assignments\02 NMPDE-HPC\03 Solution\01 Code\02 Fortran\Final IVP2_250130\UPWIND\T1_Plots and Errors\Analytical_Solution.dat"
}

# Check if files exist
for key, path in file_paths.items():
    if not os.path.exists(path):
        print(f"Error: File not found -> {path}")
        exit()

# Load data
data = {}
for key, path in file_paths.items():
    data[key] = np.loadtxt(path)

# Extract analytical solution
x_analytical = data["Analytical"][:, 0]  # Spatial coordinates
f_analytical = data["Analytical"][:, 1]  # Function values

# Store errors
errors = {}

# Process numerical datasets
for key in ["X100", "X400", "X2000"]:
    x_numerical = data[key][:, 0]  # Spatial coordinates of numerical solution
    f_numerical = data[key][:, 1]  # Function values of numerical solution
    
    # Step 1: Find indices where numerical and analytical grids match exactly
    common_indices = np.isin(x_analytical, x_numerical)

    # Step 2: Extract corresponding values
    f_analytical_matched = f_analytical[common_indices]
    f_numerical_matched = f_numerical[np.isin(x_numerical, x_analytical)]

    # Step 3: Compute error if matched points exist
    if len(f_analytical_matched) == 0 or len(f_numerical_matched) == 0:
        print(f"Warning: No exact matches found for {key}, skipping.")
        continue

    diff = f_numerical_matched - f_analytical_matched
    N = len(diff)  # Number of grid points

    # Compute norms using corrected formulas
    L1_error = np.sum(np.abs(diff)) / N  # L1 norm
    L2_error = np.sqrt(np.sum(diff**2) / N)  # L2 norm
    Linf_error = np.max(np.abs(diff))  # L-infinity norm

    # Store results
    errors[key] = [L1_error, L2_error, Linf_error]

# Convert to DataFrame
error_df = pd.DataFrame.from_dict(errors, orient="index", columns=["L1 Error", "L2 Error", "Linf Error"])

# Format in scientific notation for accuracy
error_df = error_df.map(lambda x: f"{x:.6e}")

# Print the table instead of using ace_tools
print("\n### Error Analysis (No Interpolation) ###")
print(error_df.to_string(index=True))
