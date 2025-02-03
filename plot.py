import pandas as pd
import matplotlib.pyplot as plt

# File names (since they're in the same folder as the script)
files = {
    "Analytical": "Analytical_Solution.dat",
    "Upwind": "Upwind_Solution.dat",
    #"FTCS": "FTCS_Solution.dat",
    #"Lax-Friedrichs": "Lax_Friedrichs_Solution.dat",
    #"Lax-Wendroff": "Lax_Wendroff_Solution.dat",
}

# Load datasets into a dictionary
data = {name: pd.read_csv(path, header=None, names=['x', f'u_{name.replace("-", "")}'], delim_whitespace=True) for name, path in files.items()}

# Create subplots for each numerical scheme compared to analytical solution
fig, axs = plt.subplots(2, 2, figsize=(15, 12))
axs = axs.flatten()
schemes = ["Upwind"]#, "FTCS", "Lax-Friedrichs", "Lax-Wendroff"]

for i, scheme in enumerate(schemes):
    axs[i].plot(data["Analytical"]['x'], data["Analytical"]['u_Analytical'], label='Analytical Solution', color='black', linestyle='-', linewidth=2)
    axs[i].plot(data[scheme]['x'], data[scheme][f'u_{scheme.replace("-", "")}'], label=f'{scheme} Scheme', linestyle='--', linewidth=2)
    axs[i].set_title(f'{scheme} Scheme vs Analytical Solution', fontsize=14)
    axs[i].set_xlabel('x', fontsize=12)
    axs[i].set_ylabel('u(x)', fontsize=12)
    axs[i].legend()
    axs[i].grid(True)

plt.tight_layout()
plt.show()

# Close the figure to prevent the script from hanging
plt.close(fig)
