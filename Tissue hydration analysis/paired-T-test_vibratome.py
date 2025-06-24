import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt

# Raw data
hydrated = [30.253, 31, 29.25]
dehydrated = [19.255, 23.5, 18]
rehydrated = [28.75, 38, 28.001]

# paired t-tests
t_stat_hd, p_value_hd = ttest_rel(hydrated, dehydrated)
t_stat_hr, p_value_hr = ttest_rel(hydrated, rehydrated)
t_stat_dr, p_value_dr = ttest_rel(dehydrated, rehydrated)

# Bonferroni correction
alpha = 0.05
bonferroni_alpha = alpha / 3

# Determine significance levels
significance_labels = {
    'Hydrated vs Dehydrated': '***' if p_value_hd < bonferroni_alpha else '**' if p_value_hd < 0.01 else '*' if p_value_hd < 0.05 else 'ns',
    'Hydrated vs Rehydrated': '***' if p_value_hr < bonferroni_alpha else '**' if p_value_hr < 0.01 else '*' if p_value_hr < 0.05 else 'ns',
    'Dehydrated vs Rehydrated': '***' if p_value_dr < bonferroni_alpha else '**' if p_value_dr < 0.01 else '*' if p_value_dr < 0.05 else 'ns'
}

# p-values and significance levels
print(f"Paired t-test Results (with Bonferroni-corrected alpha = {bonferroni_alpha:.3f}):")
print(f"1. Hydrated vs Dehydrated: t = {t_stat_hd:.3f}, p = {p_value_hd:.3f} -> Significance: {significance_labels['Hydrated vs Dehydrated']}")
print(f"2. Hydrated vs Rehydrated: t = {t_stat_hr:.3f}, p = {p_value_hr:.3f} -> Significance: {significance_labels['Hydrated vs Rehydrated']}")
print(f"3. Dehydrated vs Rehydrated: t = {t_stat_dr:.3f}, p = {p_value_dr:.3f} -> Significance: {significance_labels['Dehydrated vs Rehydrated']}")

# plot
means = [np.mean(hydrated), np.mean(dehydrated), np.mean(rehydrated)]
std_devs = [np.std(hydrated, ddof=1), np.std(dehydrated, ddof=1), np.std(rehydrated, ddof=1)]
labels = ['Hydrated', 'Dehydrated', 'Rehydrated']

fig, ax = plt.subplots(figsize=(8, 6))
bars = ax.bar(labels, means, yerr=std_devs, capsize=10, color=['#4DB6AC', '#81C784', '#FFB74D'], alpha=0.9, edgecolor='black')
ax.set_ylabel('Thickness (μm)', fontsize=16, weight='bold')
ax.set_xlabel('Hydration States', fontsize=16, weight='bold')
ax.set_title('Mean Tissue Thickness (Paired t-tests)', fontsize=18, weight='bold', pad=20)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Overlay individual data points
x_positions = [0, 1, 2]  # x-axis positions for each group
data_groups = [hydrated, dehydrated, rehydrated]

for x, data in zip(x_positions, data_groups):
    jitter = np.random.normal(0, 0.05, size=len(data))  # small jitter to separate points
    ax.scatter([x + j for j in jitter], data, color='black', zorder=10)

# brackets
bracket_positions = [max(means) + max(std_devs) + 1, max(means) + max(std_devs) + 4, max(means) + max(std_devs) + 8]

#  compare hyd. v Dehyd.
ax.plot([0, 1], [bracket_positions[0], bracket_positions[0]], color='black')
ax.text(0.5, bracket_positions[0] + 0.2, f'{significance_labels["Hydrated vs Dehydrated"]} (p = {p_value_hd:.3f})', ha='center', va='bottom', fontsize=12)

# Dehyd v rehyd
ax.plot([1, 2], [bracket_positions[1], bracket_positions[1]], color='black')
ax.text(1.5, bracket_positions[1] + 0.2, f'{significance_labels["Dehydrated vs Rehydrated"]} (p = {p_value_dr:.3f})', ha='center', va='bottom', fontsize=12)

# hyd. vs rehyd
ax.plot([0, 2], [bracket_positions[2], bracket_positions[2]], color='black')
ax.text(1, bracket_positions[2] + 0.2, f'{significance_labels["Hydrated vs Rehydrated"]} (p = {p_value_hr:.3f})', ha='center', va='bottom', fontsize=12)

plt.tight_layout()

# Save the plot as an SVG file
output_filename = 'tissue_thickness_paired_t_tests.svg'
plt.savefig(output_filename, format='svg')

plt.show()
