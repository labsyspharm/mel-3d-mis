#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 11:00:23 2023

@author: aj
"""

import matplotlib
import matplotlib.pyplot as plt

sc.set_figure_params(scanpy=True, fontsize=14)



MIS = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification5.csv'

MIS = pd.read_csv(MIS)
# Find the common columns between the two data frames
common_columns = MIS.columns

filtered_df = MIS[(MIS['SOX10'] > 400) | (MIS['MART1'] > 1100)]
filtered_df[filtered_df['PRAME'] > 250]

# Create a list of columns excluding 'LAG3_SPOTS' and 'GZMB_SPOTS'
columns_to_keep = [col for col in common_columns if col not in ['LAG3_SPOTS', 'GZMB_SPOTS']]
columns_to_keep.extend(['LAG3_SPOTS', 'GZMB_SPOTS'])

# deop SOX9
#common_columns = common_columns.drop(['SOX9'])
# Subset both data frames based on the common columns
MIS = MIS[columns_to_keep]

MIS.insert(0, 'CellID', MIS.index) # add cellid

MIS.dropna(inplace=True)

# fix the SOX9 gating issue
MIS['SOX9'] = MIS.apply(lambda row: 0 if row['SOX9 for gating'] == 0 else row['SOX9'], axis=1)

# write out CSV
MIS.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME.csv', index=False) 



# create adata object
feature_table_path = ['/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME.csv']

adata = sm.pp.mcmicro_to_scimap (feature_table_path, log=False, split='volume')

# scale and phenotype the cells
manual_gate = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/manual_gates_MIS.csv')
# rescale
adata = sm.pp.rescale (adata, gate=manual_gate, log=False)

# morphology
# log volume and eigne
adata.obs['eigen ratio'] = adata.obs['1_eigen'] / adata.obs['3-eigen']
adata.obs['axis length ratio'] = adata.obs['1_axislength'] / adata.obs['3_axislength']

# cell phenotyping
phenotype =  pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/phenotype_workflow_nocd15.csv')
# phenotype cells
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype)
adata.obs['phenotype'].value_counts()


bdata = adata[adata.obs['phenotype'] == 'Tumor']
# 875


# find all combinations of cells
from itertools import combinations, chain
# Define the list of elements
elements = ['5hmc', 'MART1', 'PRAME', 'SOX9', 'SOX10', 'MITF']
# Generate all possible combinations
all_combinations = chain(*[combinations(elements, i) for i in range(len(elements) + 1)])
# Convert to list and print
all_combinations_list = list(all_combinations)


data = pd.DataFrame(bdata.X, index=bdata.obs.index, columns=bdata.var.index)
# Initialize an empty dictionary to store the results
filtered_row_counts = {}
# go through each combination and find out if it exists
for i in all_combinations_list:
    filtered_rows = data[data[list(i)].apply(lambda x: all(x > 0.5), axis=1)]
    filtered_row_counts[i] = len(filtered_rows)
    

len(filtered_row_counts)
# Convert dictionary items to a list
items_list = list(filtered_row_counts.items())
# Check if the dictionary is not empty
if items_list:
    # Remove the first key-value pair
    items_list.pop(0)

# Create a new dictionary from the remaining items
filtered_row_counts = dict(items_list)
len(filtered_row_counts)

import matplotlib.pyplot as plt


###

data = filtered_row_counts.copy()

# Group the dictionary items by the number of elements in the keys
grouped_items = {}
for key, value in data.items():
    num_elements = len(key)
    if num_elements not in grouped_items:
        grouped_items[num_elements] = []
    grouped_items[num_elements].append((key, value))

# Sort each group by values in descending order
for num_elements, items in grouped_items.items():
    grouped_items[num_elements] = sorted(items, key=lambda x: x[1], reverse=False)

# Reconstruct the sorted dictionary
sorted_dict = {}
for num_elements in sorted(grouped_items.keys(), reverse=True):
    for key, value in grouped_items[num_elements]:
        sorted_dict[key] = value

# only small populations
sorted_dict = {('5hmc', 'MART1', 'PRAME', 'SOX9', 'SOX10', 'MITF'): 9,
 ('5hmc', 'MART1', 'PRAME', 'SOX9', 'MITF'): 9,
 ('5hmc', 'MART1', 'SOX9', 'SOX10', 'MITF'): 9,
 ('5hmc', 'PRAME', 'SOX9', 'SOX10', 'MITF'): 9,
 ('MART1', 'PRAME', 'SOX9', 'SOX10', 'MITF'): 18,
 ('5hmc', 'MART1', 'PRAME', 'SOX9', 'SOX10'): 20,
 ('5hmc', 'MART1', 'PRAME', 'SOX10', 'MITF'): 41,
 ('5hmc', 'PRAME', 'SOX9', 'MITF'): 9,
 ('5hmc', 'SOX9', 'SOX10', 'MITF'): 9,
 ('5hmc', 'MART1', 'SOX9', 'MITF'): 12,
 ('MART1', 'SOX9', 'SOX10', 'MITF'): 18,
 ('PRAME', 'SOX9', 'SOX10', 'MITF'): 18,
 ('5hmc', 'MART1', 'PRAME', 'SOX9'): 20,
 ('5hmc', 'MART1', 'SOX9', 'SOX10'): 20,
 ('5hmc', 'PRAME', 'SOX9', 'SOX10'): 20,
 ('MART1', 'PRAME', 'SOX9', 'SOX10'): 31,
 ('MART1', 'PRAME', 'SOX9', 'MITF'): 33,
 ('5hmc', 'MART1', 'PRAME', 'MITF'): 41,
 ('5hmc', 'PRAME', 'SOX10', 'MITF'): 41,
 ('5hmc', 'MART1', 'SOX10', 'MITF'): 44,
 ('5hmc', 'MART1', 'PRAME', 'SOX10'): 88,
 ('MART1', 'PRAME', 'SOX10', 'MITF'): 88,
 ('5hmc', 'SOX9', 'MITF'): 12,
 ('SOX9', 'SOX10', 'MITF'): 19,
 ('5hmc', 'PRAME', 'SOX9'): 20,
 ('5hmc', 'SOX9', 'SOX10'): 20,
 ('5hmc', 'MART1', 'SOX9'): 23,
 ('MART1', 'SOX9', 'SOX10'): 31,
 ('PRAME', 'SOX9', 'SOX10'): 31,
 ('PRAME', 'SOX9', 'MITF'): 33,
 ('5hmc', 'PRAME', 'MITF'): 41,
 ('5hmc', 'SOX10', 'MITF'): 44,
 ('5hmc', 'MART1', 'MITF'): 48,
 ('MART1', 'SOX9', 'MITF'): 53,
 ('MART1', 'PRAME', 'SOX9'): 66,
 ('5hmc', 'PRAME', 'SOX10'): 88,
 ('PRAME', 'SOX10', 'MITF'): 88,
 ('5hmc', 'MART1', 'PRAME'): 91,
 ('5hmc', 'MART1', 'SOX10'): 91,
 ('MART1', 'SOX10', 'MITF'): 93,
 ('MART1', 'PRAME', 'MITF'): 170,
 ('MART1', 'PRAME', 'SOX10'): 178
 }

# Extract sorted keys and values, and clean up labels
combinations = ['+'.join(key) for key in sorted_dict.keys()]
counts = list(sorted_dict.values())


# Create a vertical bar plot
plt.figure(figsize=(10, 100))
plt.barh(combinations, counts, height=0.8)
plt.xlabel('Counts')
plt.ylabel('Combinations')
plt.title('Combinations Ordered by Counts (Vertical)')
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/figures_raw/combinations_plot_smaller populations.pdf',  format='pdf')


adata.obs['phenotype'].value_counts()
adata = sm.hl.classify(adata, pos=['pMLC2'],classify_label='passed_classify', phenotype='phenotype', subclassify_phenotype=['Tissue T','CD8 T'], collapse_failed=False, label='classify')
adata.obs['classify'].value_counts()


# Show the plot
plt.show()


# un sorted

# Extract keys and values
combinations = [str(key) for key in data.keys()]
counts = list(data.values())

# Create a vertical bar plot
plt.figure(figsize=(10, 8))
plt.barh(combinations, counts)
plt.xlabel('Counts')
plt.ylabel('Combinations')
plt.title('Combinations and Their Counts (Vertical)')
plt.tight_layout()

# Show the plot
plt.show()


# color mapping 
import matplotlib.pyplot as plt
import numpy as np

data = np.full((63, 6), 0.7)


# Create a custom heatmap using Seaborn
fig, ax = plt.subplots()
sns.heatmap(data, cmap='gray',cbar=False, linewidths=1, linecolor='white', ax=ax)

# Optionally, you can add labels to the rows and columns
ax.set_xlabel("Columns")
ax.set_ylabel("Rows")

# Show the heatmap
plt.show()

