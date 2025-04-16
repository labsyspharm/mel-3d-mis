import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import os

# Raw measurement data
raw_data = {
    "5 µm": [12.949, 13.351, 8.865, 7.596, 7.682, 8.172, 9.298, 6.017, 7.194],
    "10 µm": [15.653, 17.034, 17.379, 17.437, 17.263, 18.875, 15.624, 13.811, 15.278],
    "20 µm": [29.35, 27.967, 31.074, 29.694, 22.272, 24.536, 25.896, 23.939, 25.321],
    "35 µm": [58.351, 59.387, 58.359, 52.135, 49.384, 45.575, 46.267, 46.958, 43.504]
}

nominal_values = [5, 10, 20, 35]

# Calculate means and SD from raw data and prepend the origin (0)
measured_means = [0] + [np.mean(raw_data[f"{val} µm"]) for val in nominal_values]
sd_values = [0] + [np.std(raw_data[f"{val} µm"], ddof=1) for val in nominal_values]

nominal_thicknesses = np.array([0] + nominal_values).reshape(-1, 1)
measured_thicknesses = np.array(measured_means)
error_bars = np.array(sd_values)

# Linear regression model with intercept forced to zero
model = LinearRegression(fit_intercept=False)
model.fit(nominal_thicknesses, measured_thicknesses)

slope = model.coef_[0]
predicted_thicknesses = model.predict(nominal_thicknesses)
r_squared = model.score(nominal_thicknesses, measured_thicknesses)

# Common font properties
common_font = {'fontsize': 16, 'fontname': 'Arial'}

# Create plot
plt.figure(figsize=(6, 6))

# Plot measured data with SD error bars
plt.errorbar(
    nominal_thicknesses.flatten(),
    measured_thicknesses,
    yerr=error_bars,
    fmt='o',
    color="#FF0000",
    ecolor="#FF0000",
    elinewidth=1.5,
    capsize=5,
    capthick=1.5,
    markersize=5,
    label='Measured Data (±SD)',
    clip_on=False,
    zorder=10
)

# Best-fit line
plt.plot(
    nominal_thicknesses,
    predicted_thicknesses,
    color='black',
    linewidth=3,
    label='Fit Line (Through Origin)',
    zorder=5
)

# Text inset
plt.text(
    0.05, 0.95,
    f'y = {slope:.3f}x\n$R^2 = {r_squared:.3f}$',
    transform=plt.gca().transAxes,
    fontsize=14,
    verticalalignment='top',
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square', linewidth=2),
    zorder=11
)

# Axis labels
plt.xlabel('Microtome sectioning thickness (μm)', **common_font)
plt.ylabel('Tissue thickness after hydration (μm)', **common_font)

# Tick label properties
plt.xticks(fontsize=14, fontname='Arial')
plt.yticks(fontsize=14, fontname='Arial')

# Remove grid and adjust border
plt.grid(False)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(3)
    spine.set_color('black')

# Adjust axis limits
max_y = max(measured_thicknesses + error_bars)
plt.xlim(0, nominal_thicknesses.max() * 1.1)
plt.ylim(0, max_y * 1.1)

# Save plot
plt.tight_layout()
output_dir = os.path.expanduser('~/Downloads')
output_path = os.path.join(output_dir, 'calibration_curve_SD.svg')
plt.savefig(output_path, format='svg')

# Calculate and print mean ± SD
lines = []
lines.append("Nominal Thickness (μm) | Mean Value ± SD")
lines.append("---------------------------------------")
lines.append(f"0 μm                 : {measured_means[0]:.3f} ± {sd_values[0]:.3f}")
for i, nominal in enumerate(nominal_values, start=1):
    lines.append(f"{nominal} μm                 : {measured_means[i]:.3f} ± {sd_values[i]:.3f}")

# Print and save to text file
output_txt = os.path.join(output_dir, 'mean_and_SD_values.txt')
with open(output_txt, 'w', encoding='utf-8') as f:
    for line in lines:
        print(line)
        f.write(line + "\n")

plt.show()
