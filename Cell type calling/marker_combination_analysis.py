# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:02:56 2024
@author: aj
Marker combination analysis
"""


import anndata as ad
import pandas as pd
import scimap as sm
import matplotlib.pyplot as plt
import numpy as np

# load data
# load data
MIS = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification5.csv'
IM = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iic-quantificationV3.csv'
MIS = pd.read_csv(MIS)
IM = pd.read_csv(IM)
IM.rename(columns={'GZMB_spots': 'GZMB_SPOTS'}, inplace=True)
IM.rename(columns={'LAG3_spots': 'LAG3_SPOTS'}, inplace=True)
# Find the common columns between the two data frames
common_columns = MIS.columns.intersection(IM.columns)
# Create a list of columns excluding 'LAG3_SPOTS' and 'GZMB_SPOTS'
columns_to_keep = [col for col in common_columns if col not in ['LAG3_SPOTS', 'GZMB_SPOTS']]
columns_to_keep.extend(['LAG3_SPOTS', 'GZMB_SPOTS'])

# fix sox9
MIS['SOX9'] = MIS.apply(lambda row: 0 if row['SOX9 for gating'] == 0 else row['SOX9'], axis=1)

# Subset both data frames based on the common columns
MIS = MIS[columns_to_keep]
IM = IM[columns_to_keep]

MIS.insert(0, 'CellID', MIS.index) # add cellid
IM.insert(0, 'CellID', IM.index) # add cellid

MIS.dropna(inplace=True)
IM.dropna(inplace=True)

# write out CSV
MIS.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME_sept2024.csv', index=False) 
IM.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/INV_PRAME_sept2024.csv', index=False) 


feature_table_path = ['/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME_sept2024.csv',
                      '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/INV_PRAME_sept2024.csv']
adata = sm.pp.mcmicro_to_scimap (feature_table_path, log=False, split='volume')


# load gates
manual_gate = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-gates.csv')

# rescale data
adata = sm.pp.rescale (adata, gate=manual_gate, log=False)

# cell phenotyping
phenotype =  pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/phenotype_workflow_nocd15.csv')
# phenotype cells
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype)
adata.obs['phenotype'].value_counts()
# rename cells
rename= {'CD4 T': ['T cells']}
adata = sm.hl.rename(adata, rename, from_column='phenotype', to_column='phenotype')

bdata = adata[adata.obs['phenotype'] == 'Tumor']
bdata = bdata[bdata.obs['imageid'] == 'MIS_PRAME_sept2024']


# calculate the combinations

from itertools import combinations, chain
import pandas as pd

# Define the list of elements
elements = ['5hmc', 'MART1', 'PRAME', 'SOX9', 'SOX10', 'MITF']

# Generate all possible combinations of elements
all_combinations = chain(*[combinations(elements, i) for i in range(len(elements) + 1)])
# Convert to list
all_combinations_list = list(all_combinations)

# Create a DataFrame (replace this with your actual data)
data = pd.DataFrame(bdata.X, index=bdata.obs.index, columns=bdata.var.index)

# NEW
# Initialize an empty dictionary to store the results
filtered_row_counts = {}

# Initialize an empty list to store the data for the dataframe
dataframe_rows = []

# Iterate over each combination
for pos_comb in all_combinations_list:
    # Define the positive condition (values > 0.5 for selected elements)
    if len(pos_comb) > 0:
        pos_condition = data[list(pos_comb)].apply(lambda x: all(x > 0.5), axis=1)
    else:
        pos_condition = pd.Series([True] * len(data), index=data.index)  # If no positive elements, no restriction on pos_condition
    
    # Define the negative condition (values < 0.5 for remaining elements)
    remaining_elements = [e for e in elements if e not in pos_comb]
    if len(remaining_elements) > 0:
        neg_condition = data[remaining_elements].apply(lambda x: all(x < 0.5), axis=1)
    else:
        neg_condition = pd.Series([True] * len(data), index=data.index)  # If no negative elements, no restriction on neg_condition

    # Combine both conditions (positive AND negative)
    filtered_rows = data[pos_condition & neg_condition]
    
    # Store the count of filtered rows
    filtered_row_counts[pos_comb] = len(filtered_rows)

    # Get the indices of rows that pass the condition
    passed_indices = filtered_rows.index.tolist()

    # If no rows pass the condition, append a 0 count for the combination
    if len(filtered_rows) == 0:
        dataframe_rows.append({
            'Condition': pos_comb,
            'Index': None,  # No index because no row passed
            'Total Passed': 0,
            **{marker: 0 for marker in elements}  # Binarized column for each marker (0)
        })
    else:
        # Store the data for each passed row
        for index in passed_indices:
            # Get the row values for binarized markers (1 if value > 0.5, else 0)
            binarized_values = {marker: int(data.at[index, marker] > 0.5) for marker in elements}
            
            dataframe_rows.append({
                'Condition': pos_comb,
                'Index': index,
                'Total Passed': len(filtered_rows),
                **binarized_values  # Add binarized values for each marker
            })

# Create the dataframe from the list of dictionary records
df_results = pd.DataFrame(dataframe_rows)




# =============================================================================
# # Initialize an empty dictionary to store the results
# filtered_row_counts = {}
# 
# # Initialize an empty list to store the data for the dataframe
# dataframe_rows = []
# 
# # Iterate over each combination
# for pos_comb in all_combinations_list:
#     # Define the positive condition (values > 0.5 for selected elements)
#     if len(pos_comb) > 0:
#         pos_condition = data[list(pos_comb)].apply(lambda x: all(x > 0.5), axis=1)
#     else:
#         pos_condition = pd.Series([True] * len(data), index=data.index)  # If no positive elements, no restriction on pos_condition
#     
#     # Define the negative condition (values < 0.5 for remaining elements)
#     remaining_elements = [e for e in elements if e not in pos_comb]
#     if len(remaining_elements) > 0:
#         neg_condition = data[remaining_elements].apply(lambda x: all(x < 0.5), axis=1)
#     else:
#         neg_condition = pd.Series([True] * len(data), index=data.index)  # If no negative elements, no restriction on neg_condition
# 
#     # Combine both conditions (positive AND negative)
#     filtered_rows = data[pos_condition & neg_condition]
#     
#     # Store the count of filtered rows
#     filtered_row_counts[pos_comb] = len(filtered_rows)
# 
#     # Get the indices of rows that pass the condition
#     passed_indices = filtered_rows.index.tolist()
# 
#     # Store the data for the dataframe
#     for index in passed_indices:
#         dataframe_rows.append({
#             'Condition': pos_comb,
#             'Index': index,
#             'Total Passed': len(filtered_rows)
#         })
# 
# # Create the dataframe from the list of dictionary records
# df_results = pd.DataFrame(dataframe_rows)
# =============================================================================


df_results = df_results.sort_values(by='Total Passed', ascending=True)
# save the results
df_results.to_csv('/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/projects/2023_3D/3-NatureMet Transfer-July2024/Figures/individual panels/tumor combination analysis/index_of_combination_with_binarized_values.csv',index=False)



# Output the filtered_row_counts dictionary
print(filtered_row_counts)


# plot the combinations
def plot_combinations(results, figsize=(10, 6), ylim=None, ytick_interval=None):
    # Prepare the data for plotting
    combinations = ['+'.join(comb) if comb else 'None' for comb in results.keys()]
    counts = list(results.values())
    
    # Create the bar plot with a specified figure size
    plt.figure(figsize=figsize)
    plt.bar(combinations, counts, color='grey', edgecolor='black')

    # Add labels and title
    plt.xlabel('Marker Combinations', fontsize=12)
    plt.ylabel('Number of Cells', fontsize=12)
    plt.title('Number of Cells for Each Marker Combination', fontsize=14)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90)
    
    # Set y-axis limit if provided
    if ylim:
        plt.ylim(ylim)
    
    # Set y-axis ticks interval if provided
    if ytick_interval:
        max_y = ylim[1] if ylim else max(counts)
        plt.yticks(np.arange(0, max_y + ytick_interval, ytick_interval))

    # Display the plot
    plt.tight_layout()
    plt.show()

# Call the function with your results, setting figure size, y-limit, and y-axis tick interval
plot_combinations(results=filtered_row_counts, figsize=(12, 12), ylim=(0, 110), ytick_interval=5)






