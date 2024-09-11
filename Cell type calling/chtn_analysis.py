#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 16:02:55 2024
@author: aj
CHTN data analysis
"""

import scimap as sm
import pandas as pd
import os
import matplotlib.pyplot as plt
import anndata as ad
import seaborn as sns
import numpy as np
import pingouin as pg
import itertools

%matplotlib qt

plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
plt.rcParams['font.family']='sans serif'
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype']=42

# take in raw data and export the clened data
def process_csv(csv_path, output_directory, filter_column='filtermask', filter_value=0, target_column=None):
    """
    Processes a CSV file: adds a 'CellID' column filled with index numbers, capitalizes all column names, 
    drops a specified column if provided, drops rows based on a filter column and NaN values, 
    and saves the processed file.

    :param csv_path: Path to the input CSV file.
    :param output_directory: Directory where the processed file will be saved.
    :param filter_column: Name of the column to filter by. Defaults to 'filtermask'. If None, no filtering is applied.
    :param filter_value: Value to retain in the filter column. Defaults to 0.
    :param target_column: Name of the column to be dropped, or None to skip dropping any column.
    """
    # Load the CSV file
    df = pd.read_csv(csv_path)
    
    # Insert 'CellID' column as the first column, filled with index numbers starting from 1
    df.insert(0, 'CellID', range(1, 1 + len(df)))
    print("CellID column added.")
    
    # Capitalize all column names
    df.columns = [col.upper() for col in df.columns]
    print("All column names have been capitalized.")
    
    # Update filter_column to uppercase to match the updated column names, if it's not None
    if filter_column is not None:
        filter_column = filter_column.upper()
    
    # Drop the specified target_column, if it exists and is not None
    if target_column is not None:
        target_column = target_column.upper()  # Update target_column to uppercase
        if target_column in df.columns:
            df = df.drop(columns=[target_column])
            print(f"Column '{target_column}' has been dropped.")
    
    # Drop rows where the filter column value is not equal to filter_value, if filter_column is specified
    if filter_column is not None and filter_column in df.columns:
        df = df[df[filter_column] == filter_value]
    
    # Drop rows that contain any NaN values
    df = df.dropna()
    
    # Prepare the output file path
    base_name = os.path.basename(csv_path)
    name, ext = os.path.splitext(base_name)
    output_file_path = os.path.join(output_directory, f"{name}_cleaned{ext}")
    
    # Save the processed file
    df.to_csv(output_file_path, index=False)
    print(f"File saved to {output_file_path}")



# usage
csv_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTNa-quantification.csv'
output_directory = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/'
process_csv(csv_path, output_directory, target_column='TIM3.1')


# remove 1st row as MX1 column contained a letter
def remove_first_row_and_overwrite(csv_path):
    # Load the CSV file normally to get column names
    df = pd.read_csv(csv_path)
    # Remove the first row of data, keeping the header intact
    #df = df.iloc[1:].reset_index(drop=True)
    # Overwrite the original CSV file with the modified DataFrame
    df.to_csv(csv_path, index=False)

# Example usage
csv_path_mod = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CHTNa-quantification_cleaned.csv'
remove_first_row_and_overwrite(csv_path_mod)

# check
df = pd.read_csv(csv_path_mod)
df.dtypes
df.shape



# check number of + cells
def count_rows_above_value(df, column_name, value):
    """
    Counts the number of rows in a DataFrame where the value in a specified column is greater than a user-specified value.

    :param df: pandas DataFrame to search.
    :param column_name: The name of the column to check.
    :param value: The value to compare against.
    :return: The count of rows where the column's value is greater than the specified value.
    """
    if column_name not in df.columns:
        print(f"Column '{column_name}' does not exist in the DataFrame.")
        return 0
    
    count = (df[column_name] > value).sum()
    return count

# run
count_rows_above_value(df=df, column_name='MART1', value=600)


##############################################################################
# SCIMAP
##############################################################################

image_path = ['/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CHTNa-quantification_cleaned.csv']
adata = sm.pp.mcmicro_to_scimap (image_path, drop_markers = ["HOECHST"], log=False, split='LAG3SPOTS', CellId='CELLID')
adata

# load the manual gates
manual_gate = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/manual_gates.csv')
adata = sm.pp.rescale (adata, gate=manual_gate, log=False)

# phenotype
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/phenotype_workflow_broad.csv')
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, label="phenotype_broad") 
adata.obs['phenotype_broad'].value_counts()


# rename T cells to CD T
rename= {'CD4 T': ['T cells']}
adata = sm.hl.rename (adata, rename, from_column='phenotype_broad', to_column='phenotype_broad')
    
# save h5ad
adata.write('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.h5ad')

    
# export CSV
data = sm.hl.scimap_to_csv (adata, data_type='scaled')
data.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.csv', index=False)

# check
df1 = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.csv')
count_rows_above_value(df=df1, column_name='MART1', value=0.5)



##############################################################################
# Finer classification of phenotypes for spatial analysis
##############################################################################

# load adata and add 2 new columns with regards to exhaustion and proliferation
adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.h5ad')
# extract data
data = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
meta = adata.obs.copy()
# load the original CSV file and binarize it
binarize_data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.csv', index_col=0)
drop_columns = ['CellID', 'VOLUME (VOXELS)', 'CENTROID_X', 'CENTROID_Y', 'CENTROID_Z', 'FILTERMASK', 'CELLID', 'imageid','phenotype_broad']
marker_columns = ['CD3E', 'MART1', 'PCNA', 'CD8A', 'SOX10', 'FOXP3', 'CD103','CD11C', 'PD1', 'TIM3', 'CD4', 'CD45RO', 'TCF1', 'KI67']
spot_columns = [ 'LAG3SPOTS','GZMBSPOTS', 'PDL1SPOTS', 'MX1SPOTS']
# function to binarize the matrix
def binarize_dataframe(df, drop_columns=None, marker_columns=None, threshold_marker=0.50, spot_columns=None, threshold_spot=1):

    # Drop specified columns
    if drop_columns is not None:
        df = df.drop(columns=drop_columns, errors='ignore')
    
    # Convert marker columns
    if marker_columns is not None:
        for column in marker_columns:
            if column in df.columns:
                df[column] = (df[column] >= threshold_marker).astype(int)
    
    # Convert spot columns
    if spot_columns is not None:
        for column in spot_columns:
            if column in df.columns:
                df[column] = (df[column] > threshold_spot).astype(int)
    
    return df


binarized_data = binarize_dataframe(df = binarize_data, drop_columns = drop_columns, marker_columns = marker_columns, spot_columns = spot_columns)
binarized_data.columns

# Update the binarized_data with combined columns
def update_csv_with_conditions(df, pattern_dict):    
    # Step 3: Iterate over the dictionary to update the DataFrame
    for key, columns in pattern_dict.items():
        # Using any() along axis=1 checks if any of the specified columns have a value of 1
        # The result is then converted to int (True -> 1, False -> 0) and assigned to the new column
        df[key] = df[columns].any(axis=1).astype(int)
    return df


pattern_dict = {
  "TIM3_LAG3": ['TIM3', 'LAG3SPOTS'],
  "PCNA_KI67": ['PCNA', 'KI67'],
}

updated_binarized_data = update_csv_with_conditions(binarized_data, pattern_dict)

# Now add the new columns to data
add_columns = ['GZMBSPOTS','PDL1SPOTS', 'MX1SPOTS', 'TIM3_LAG3', 'PCNA_KI67']
updated_data = pd.concat([data, updated_binarized_data[add_columns]], axis=1)

# creare anndata object
adata = ad.AnnData(updated_binarized_data)
adata.obs = meta

# finer phenotyping
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/phenotype_workflow_fine.csv')
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, label="phenotype_fine") 
adata.obs['phenotype_fine'].value_counts()

# rename T cells to CD T
rename= {'CD4 T': ['T cells']}
adata = sm.hl.rename (adata, rename, from_column='phenotype_fine', to_column='phenotype_fine') 
adata.obs['phenotype_fine'].value_counts()


def replace_brackets_in_column(df, column_name):
    # Ensure the column exists in DataFrame
    if column_name in df.columns:
        # Replace [ and ] with _ in the specified column
        df[column_name] = df[column_name].str.replace('(', '_', regex=False).str.replace(')', '_', regex=False)
    else:
        print(f"Column '{column_name}' not found in DataFrame.")

adata.obs['phenotype_fine_mod'] = adata.obs['phenotype_fine'].values
replace_brackets_in_column(adata.obs, column_name='phenotype_fine_mod')
adata.obs['phenotype_fine_mod'].value_counts()


rename_dict = {'CD103- Tpex': ['TCF1+ CD103-  T _Tex_','TCF1+ CD103-  T _Tp_'], # tex and Tp
               'CD103+ Tpex': ['TCF1+ CD103+  T _Tex_','TCF1+ CD103+  T _Tp_'], # tex and Tp
               'Tmem': ['TCF1- CD103+  T _Tex_','TCF1- CD103+  T _Tp_','TCF1- CD103+  T _activated_','TCF1- CD103+ T _naive_'], # naive, activated, tex and Tp
               'Teff': ['TCF1- CD103-  T _Tp_'], #tp
               'Tex': ['TCF1- CD103-  T _Tex_']} #tex


adata = sm.hl.rename(adata, rename=rename_dict, from_column='phenotype_fine_mod', to_column='phenotype_fine_majorcategories')
adata.obs['phenotype_fine_majorcategories'].value_counts()


# write
adata.write('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.h5ad')

##############################################################################
# Binarize the table for finner sb classifications
##############################################################################

binarize_data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.csv', index_col=0)
binarize_data.columns

# function to binarize the matrix
def binarize_dataframe(df, drop_columns=None, marker_columns=None, threshold_marker=0.50, spot_columns=None, threshold_spot=1):
    """
    Process the DataFrame according to specified operations.

    :param df: pandas DataFrame to process.
    :param drop_columns: List of columns to drop. If None, no columns are dropped.
    :param marker_columns: List of columns to apply marker threshold conversion. If None, skip this operation.
    :param threshold_marker: Threshold value for marker columns. Values >= this threshold are set to 1, else to 0.
    :param spot_columns: List of columns to apply spot threshold conversion. If None, skip this operation.
    :param threshold_spot: Threshold value for spot columns. Values > this threshold are set to 1, values <= to 0.
    """
    
    # Drop specified columns
    if drop_columns is not None:
        df = df.drop(columns=drop_columns, errors='ignore')
    
    # Convert marker columns
    if marker_columns is not None:
        for column in marker_columns:
            if column in df.columns:
                df[column] = (df[column] >= threshold_marker).astype(int)
    
    # Convert spot columns
    if spot_columns is not None:
        for column in spot_columns:
            if column in df.columns:
                df[column] = (df[column] > threshold_spot).astype(int)
    
    return df


# run

drop_columns = ['CellID', 'VOLUME (VOXELS)', 'CENTROID_X',
                'CENTROID_Y', 'CENTROID_Z', 'FILTERMASK', 'CELLID', 'imageid',
                'phenotype_broad']

marker_columns = ['CD3E', 'MART1', 'PCNA', 'CD8A', 'SOX10', 'FOXP3', 'CD103',
                  'CD11C', 'PD1', 'TIM3', 'CD4', 'CD45RO', 'TCF1', 'KI67']

spot_columns = [ 'LAG3SPOTS','GZMBSPOTS', 'PDL1SPOTS', 'MX1SPOTS']



# subset CD8 T
CD8T = binarize_data[binarize_data['phenotype_broad'] == 'CD8 T']
CD8T_binarized_data = binarize_dataframe(df = CD8T, drop_columns = drop_columns, marker_columns = marker_columns, spot_columns = spot_columns)
CD8T_binarized_data.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD8T_binary.csv', index=True)


# subset CD4 T
CD4T = binarize_data[binarize_data['phenotype_broad'] == 'CD4 T']
CD4T_binarized_data = binarize_dataframe(df = CD4T, drop_columns = drop_columns, marker_columns = marker_columns, spot_columns = spot_columns)
CD4T_binarized_data.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD4T_binary.csv', index=True)



##############################################################################
# Binary classifier
##############################################################################



# function to add columns to the binary martix 
def update_csv_with_conditions(csv_path, pattern_dict):
    # Step 2: Read the CSV file into a DataFrame
    df = pd.read_csv(csv_path,index_col=0)
    
    # Step 3: Iterate over the dictionary to update the DataFrame
    for key, columns in pattern_dict.items():
        # Using any() along axis=1 checks if any of the specified columns have a value of 1
        # The result is then converted to int (True -> 1, False -> 0) and assigned to the new column
        df[key] = df[columns].any(axis=1).astype(int)
    
    # Step 4: Overwrite the original CSV file with the modified DataFrame
    df.to_csv(csv_path, index=True)

# Example usage
pattern_dict = {
  "TIM3_LAG3": ['TIM3', 'LAG3SPOTS'],
  "PCNA_KI67": ['PCNA', 'KI67'],
}
csv_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD8T_binary.csv'
update_csv_with_conditions(csv_path, pattern_dict)
csv_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD4T_binary.csv'
update_csv_with_conditions(csv_path, pattern_dict)




# example file
CD8T_binarized_data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD8T_binary.csv',index_col=0)
CD4T_binarized_data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD4T_binary.csv',index_col=0)

# check
count_rows_above_value(df=CD8T_binarized_data, column_name='TCF1', value=0)

# works


# =============================================================================
# def hierarchical_classification_archive (data_path, condition_path, fraction_type='all', return_type='percentage'):
#     # Load datasets
#     data = pd.read_csv(data_path, index_col=0)
#     conditions = pd.read_csv(condition_path)
#     
#     # Prepare a DataFrame to hold classification flags for each cell
#     flags = pd.DataFrame(index=data.index)
#     flags['All'] = True  # All cells start as part of the 'All' category
#     
#     # Enhanced method to track and classify based on explicit parent-child relationships
#     def classify_and_split(base_flag, conditions, level=1, prefix='', parent_flag='All'):
#         nonlocal flags
#         current_level_conditions = conditions[conditions['level'] == level]
#         
#         if current_level_conditions.empty:
#             return
#         
#         for _, condition_row in current_level_conditions.iterrows():
#             for marker, condition in condition_row.dropna().items():
#                 if marker == 'level':  # Skip the 'level' column
#                     continue
#                 
#                 pos_flag = f'{prefix}{marker}+'
#                 neg_flag = f'{prefix}{marker}-'
#                 
#                 flags[pos_flag] = base_flag & (data[marker] == 1)
#                 flags[neg_flag] = base_flag & (data[marker] == 0)
#                 
#                 # Recursively classify and split with updated parent flag
#                 classify_and_split(flags[pos_flag], conditions, level + 1, prefix=pos_flag + ' ', parent_flag=pos_flag)
#                 classify_and_split(flags[neg_flag], conditions, level + 1, prefix=neg_flag + ' ', parent_flag=neg_flag)
#     
#     classify_and_split(flags['All'], conditions, level=1)
#     
#     # Initialize the results Series
#     results = pd.Series(dtype=float if return_type == 'percentage' else int)
#     
#     # Initialize a list to hold classification types and individual indices
#     classification_indices = []
#     
#     for flag in flags.columns[1:]:  # Exclude 'All'
#         current_population = flags[flag].sum()
#         
#         if fraction_type == 'all':
#             base_population = len(data) if return_type == 'percentage' else 1
#         elif fraction_type == 'preceding':
#             parent_flag = flag.rsplit(' ', 2)[0] if ' ' in flag else 'All'
#             base_population = flags[parent_flag].sum() if parent_flag in flags else 1
#         
#         if return_type == 'percentage':
#             results[flag] = (current_population / base_population * 100) if base_population else 0
#         else:
#             results[flag] = current_population
#         
#         # Loop through each index for the current classification and append individually
#         for index in flags[flags[flag]].index:
#             classification_indices.append((flag, index))
#     
#     # Convert the list to a DataFrame
#     classification_df = pd.DataFrame(classification_indices, columns=['Classification', 'Index'])
#     
#     return results, classification_df
# =============================================================================


# Use to derive number of terminal populations too
def hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='percentage'):
    # Load datasets
    data = pd.read_csv(data_path, index_col=0)
    conditions = pd.read_csv(condition_path)
    
    # Prepare a DataFrame to hold classification flags for each cell
    flags = pd.DataFrame(index=data.index)
    flags['All'] = True  # All cells start as part of the 'All' category
    
    def classify_and_split(base_flag, conditions, level=1, prefix='', parent_flag='All'):
        nonlocal flags
        current_level_conditions = conditions[conditions['level'] == level]
        
        if current_level_conditions.empty:
            return
        
        for _, condition_row in current_level_conditions.iterrows():
            for marker, condition in condition_row.dropna().items():
                if marker == 'level':  # Skip the 'level' column
                    continue
                
                pos_flag = f'{prefix}{marker}+'
                neg_flag = f'{prefix}{marker}-'
                
                flags[pos_flag] = base_flag & (data[marker] == 1)
                flags[neg_flag] = base_flag & (data[marker] == 0)
                
                classify_and_split(flags[pos_flag], conditions, level + 1, prefix=pos_flag + ' ', parent_flag=pos_flag)
                classify_and_split(flags[neg_flag], conditions, level + 1, prefix=neg_flag + ' ', parent_flag=neg_flag)
    
    classify_and_split(flags['All'], conditions, level=1)
    
    results = pd.Series(dtype=float if return_type == 'percentage' else int)
    classification_indices = []
    
    for flag in flags.columns[1:]:  # Exclude 'All'
        current_population = flags[flag].sum()
        
        if fraction_type == 'all':
            base_population = len(data) if return_type == 'percentage' else 1
        elif fraction_type == 'preceding':
            parent_flag = flag.rsplit(' ', 2)[0] if ' ' in flag else 'All'
            base_population = flags[parent_flag].sum() if parent_flag in flags else 1
        
        if return_type == 'percentage':
            results[flag] = (current_population / base_population * 100) if base_population else 0
        else:
            results[flag] = current_population
        
        for index in flags[flags[flag]].index:
            classification_indices.append((flag, index))
    
    classification_df = pd.DataFrame(classification_indices, columns=['Classification', 'Index'])
    
    # Determine terminal populations that are also present in results before filtering
    terminal_flags = [flag for flag in flags.columns if not any(flags.columns.str.startswith(flag + ' ')) and flag in results.index]
    terminal_populations = results[terminal_flags]
    terminal_populations_gt_1_percent = terminal_populations[terminal_populations > 1]

    # The number of terminal populations greater than 1%
    num_terminal_populations_gt_1_percent = len(terminal_populations_gt_1_percent)
    
    # Create a DataFrame for terminal populations and their percentages
    terminal_populations_df = pd.DataFrame(terminal_populations_gt_1_percent).reset_index()
    terminal_populations_df.columns = ['Classification', 'Percentage']
    
    return results, classification_df, num_terminal_populations_gt_1_percent, terminal_populations_df


# Usage example
# CD8 T
data_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD8T_binary.csv'
condition_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD8_binary_classifier.csv'
#percentages_abs, classifications_abs = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='absolute')
#percentages_per, classifications_per = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='percentage')
percentages_abs, classifications_abs, terminal_populations, df_ter = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='absolute')
percentages_per, classifications_per, terminal_populations, per_df_ter = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='percentage')
# write out the results
percentages_abs.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_absolute.csv', index=True)
percentages_per.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_percentage.csv', index=True)
per_df_ter.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_terminal_populations.csv', index=True)
classifications_abs.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_cell_index.csv', index=True)

# CD4 T
data_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD4T_binary.csv'
condition_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/CD4_binary_classifier.csv'
#percentages_abs, classifications_abs = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='absolute')
#percentages_per, classifications_per = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='percentage')
percentages_abs, classifications_abs, terminal_populations, df_ter = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='absolute')
percentages_per, classifications_per, terminal_populations, per_df_ter = hierarchical_classification(data_path, condition_path, fraction_type='all', return_type='percentage')
# write out the results
percentages_abs.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_absolute.csv', index=True)
percentages_per.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_percentage.csv', index=True)
per_df_ter.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_terminal_populations.csv', index=True)
classifications_abs.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_cell_index.csv', index=True)



# plot the results
%matplotlib qt

def plot_filtered_data(csv_path, cutoff):
    # Step 1: Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Assuming the first column as the category and the second column as the value
    # Update column names based on actual CSV structure if needed
    category_col = df.columns[0]
    value_col = df.columns[1]
    
    # Step 2: Remove rows where the value is below the cutoff
    filtered_df = df[df[value_col] >= cutoff]
    
    # Sort the DataFrame based on the value column in descending order
    sorted_df = filtered_df.sort_values(by=value_col, ascending=False)
    
    # Step 3: Plot
    plt.figure(figsize=(10, 6))
    plt.bar(sorted_df[category_col], sorted_df[value_col])
    plt.xlabel(category_col)
    plt.ylabel('Percentage')
    plt.xticks(rotation=90)  # Rotate x-axis labels by 90 degrees
    
    # Setting y-axis to show every 5 percent
    y_max = sorted_df[value_col].max()  # Find the maximum value
    y_min = 0  # Assuming you want to start from 0
    plt.yticks(range(int(y_min), int(y_max) + 5, 5))  # Adjust the range and step as necessary
    
    plt.title('Bar Plot of CD8 Subclasses Percentage Above Cutoff')
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels
    plt.show()


# Path to the uploaded CSV file
csv_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_percentage.csv'
csv_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_percentage.csv'

# Execute the function with the provided CSV path and cutoff
plot_filtered_data(csv_path, cutoff=1)



##############################################################################
# Find out what percent of cells are close to tumor cells
##############################################################################

adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.h5ad')

from sklearn.neighbors import BallTree
data = pd.DataFrame({'x': adata.obs['CENTROID_X'], 'y': adata.obs['CENTROID_Y'], 'z': adata.obs['CENTROID_Z'], 'phenotype': adata.obs['phenotype_broad']})

kdt = BallTree(data[['x','y','z']], metric='euclidean') 
ind = kdt.query(data[['x','y','z']], k=6, return_distance=False)
#for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
neighbours_pheno = neighbours.copy()
phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
for i in neighbours_pheno.columns:
            neighbours_pheno[i] = neighbours_pheno[i].map(phenomap, na_action='ignore')
indexmap = dict(zip(list(range(len(ind))), data.index)) # Used for mapping
for i in neighbours.columns:
            neighbours[i] = neighbours[i].map(indexmap, na_action='ignore')         
            
# keep only tumor cells
tumor_cells = neighbours_pheno[neighbours_pheno[0] == 'Tumor']
neighbours = neighbours.loc[tumor_cells.index]
# drop the tumor column
tumor_cells = tumor_cells.drop(columns=[tumor_cells.columns[0]])
neighbours = neighbours.drop(columns=[neighbours.columns[0]])

# indetify a list of unique cells that are close to tumor cells. 
tumor_neighbours = pd.unique(neighbours.values.ravel())#.tolist()

# for each cell type in the classifications_abs find the percent of that are close to a tumor cell
def percent_cells_near_tumor (classifications_abs, tumor_neighbours):
    results = pd.DataFrame(columns=['cell type', 'total_cells', 'proximate_n', 'percent'])
    cell_types = classifications_abs['Classification'].unique()
    tumor_neighbours_set = set(tumor_neighbours)
    for i in cell_types:
        filtered_indices = classifications_abs[classifications_abs['Classification'] == i]['Index']
        index_set = set(filtered_indices)
        common_elements = index_set.intersection(tumor_neighbours_set)
        percent = (len(common_elements) / len(filtered_indices)) * 100
        # add to results
        temp_df = pd.DataFrame([[i, len(filtered_indices),len(common_elements),  percent]], columns=['cell type', 'total_cells', 'proximate_n', 'percent'])
        results = pd.concat([results, temp_df], ignore_index=True)
    return results


# percent of CD8 T cells
cd8_classifications_abs = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_cell_index.csv', index_col=0)
cd4_classifications_abs = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD4_subclasses_cell_index.csv', index_col=0)

tumor_proximate_CD8 = percent_cells_near_tumor (cd8_classifications_abs, tumor_neighbours)
tumor_proximate_CD4 = percent_cells_near_tumor (cd4_classifications_abs, tumor_neighbours)

# write it out
tumor_proximate_CD8.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/tumor_proximate_CD8_subclasses.csv', index=False)
tumor_proximate_CD4.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/tumor_proximate_CD4_subclasses.csv', index=False)


# some testing

def filter_rows_by_markers(df, markers, threshold, operation):
    """
    Filters rows based on the specified operation across given markers.
    
    Parameters:
    - df (pd.DataFrame): The input dataframe.
    - markers (list of str): Column names to apply the threshold on.
    - threshold (float): The threshold value.
    - operation (str): The operation to apply ("and" or "or").
    
    Returns:
    - pd.DataFrame: A filtered dataframe based on the operation.
    """
    if operation == "and":
        # Use all() across axis=1 to find rows where all markers exceed the threshold
        mask = df[markers].apply(lambda x: (x > threshold).all(), axis=1)
    elif operation == "or":
        # Use any() across axis=1 to find rows where any marker exceeds the threshold
        mask = df[markers].apply(lambda x: (x > threshold).any(), axis=1)
    else:
        raise ValueError("Operation must be 'and' or 'or'.")

    return df[mask]


# testing with known combinations
filtered_df = filter_rows_by_markers(df=CD8T_binarized_data, markers=['PCNA','KI67'], threshold=0, operation='or')
len(filtered_df)

filtered_df = filter_rows_by_markers(df=CD8T_binarized_data, markers=['PCNA_KI67'], threshold=0, operation='or')
len(filtered_df)

len(filter_rows_by_markers(df=df1, markers=['CD4'], threshold=0.49999, operation='or'))


##############################################################################
# Investigate the expression profiles of PD1 and LAG3 in CD8 T cells
##############################################################################

adata = ad.read('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.h5ad')
t_cells = adata[adata.obs['phenotype_broad'] == 'CD8 T']
mat = pd.DataFrame(np.log1p(t_cells.raw.X), index=t_cells.obs.index, columns=t_cells.var.index)
mat1 = pd.DataFrame(np.log1p(t_cells.X), index=t_cells.obs.index, columns=t_cells.var.index)



def plot_kde_distributions(df, columns, colors=None, x_min=None, x_max=None):
    """
    Plots the KDE (Kernel Density Estimate) distributions of specified columns in a DataFrame on the same plot,
    using either a provided list of colors or generating a color palette automatically. Allows specifying x-axis limits.

    """
    plt.figure(figsize=(10, 6))
    from scipy.stats import gaussian_kde

    
    # Generate a color palette if no color list is provided
    if colors is None:
        colors = plt.cm.viridis(np.linspace(0, 1, len(columns)))
    elif len(colors) < len(columns):
        print("Warning: Not enough colors provided, some columns will use the same color.")
    
    for i, column in enumerate(columns):
        if column in df.columns:
            data = df[column].dropna()  # Drop NA values for plotting
            if data.empty:
                print(f"Column '{column}' contains only NA values.")
                continue
            kde = gaussian_kde(data)
            x_range = np.linspace(data.min() if x_min is None else x_min, 
                                  data.max() if x_max is None else x_max, 1000)
            plt.plot(x_range, kde(x_range), label=column, color=colors[i % len(colors)])
        else:
            print(f"Column '{column}' not found in DataFrame.")
    
    # Set x-axis limits if specified
    if x_min is not None or x_max is not None:
        plt.xlim(x_min, x_max)
    
    plt.title('KDE of Columns')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(False)
    plt.show()


# plot
plot_kde_distributions(mat, ['TCF1','PD1'], colors=['red','blue'])
plot_kde_distributions(mat1, ['TCF1','PD1'], colors=['red','blue'],x_min=0, x_max=1)



##############################################################################
# Investigate the % of LAG3 vs TIm3 + cells
##############################################################################

def plot_threshold_comparison_save_pdf(df, col1, col2, threshold1=0.5, threshold2=0.5, directory='.', filename='plot.pdf'):
    """
    Calculates and plots the percentages of rows with values greater than specified thresholds in two columns,
    including the percentage of rows meeting the threshold in each column and common to both. Saves the plot as a PDF file.

    Parameters:
    - df: pandas.DataFrame
        The DataFrame containing the data.
    - col1: str
        The name of the first column to analyze.
    - col2: str
        The name of the second column to analyze.
    - threshold1: float, optional
        The threshold value for the first column (default is 0.5).
    - threshold2: float, optional
        The threshold value for the second column (default is 0.5).
    - directory: str, optional
        The directory where the PDF file will be saved (default is current directory).
    - filename: str, optional
        The name of the PDF file (default is 'plot.pdf').

    Returns:
    None
    """
    # Check if columns exist in DataFrame
    if col1 not in df.columns or col2 not in df.columns:
        print(f"One or both columns '{col1}' and '{col2}' not found in DataFrame.")
        return
    
    # Calculate the conditions based on the specified thresholds
    condition_col1 = df[col1] > threshold1
    condition_col2 = df[col2] > threshold2
    
    # Count the number of rows meeting the conditions
    count_col1_only = ((condition_col1) & (~condition_col2)).sum()
    count_col2_only = ((condition_col2) & (~condition_col1)).sum()
    count_both = (condition_col1 & condition_col2).sum()
    
    # Calculate total rows for normalization
    total_rows = len(df)
    
    # Calculate percentages
    percent_col1_only = (count_col1_only / total_rows) * 100
    percent_col2_only = (count_col2_only / total_rows) * 100
    percent_both = (count_both / total_rows) * 100
    
    # Plotting
    categories = [col1, col2, 'Common to Both']
    percentages = [percent_col1_only, percent_col2_only, percent_both]
    
    plt.figure(figsize=(8, 6))
    plt.grid(False)
    plt.bar(categories, percentages, color=['blue', 'green', 'purple'])
    plt.ylabel('Percentage (%)')
    plt.title('Percentage of Rows Exceeding Thresholds')
    plt.ylim(0, max(percentages) + 10)  # Adjust y-axis to fit labels
    for i, v in enumerate(percentages):
        plt.text(i, v + 1, f"{v:.2f}%", ha='center', va='bottom')
    
    # Saving the plot as a PDF
    full_path = f"{directory}/{filename}"
    plt.savefig(full_path, format='pdf')
    
    plt.show()
    print(f"Plot saved as {full_path}")



# manual plot
t_cells = adata[adata.obs['phenotype_broad'] == 'CD8 T']
mat2 = pd.DataFrame(np.log1p(t_cells.X), index=t_cells.obs.index, columns=t_cells.var.index)
mat2 = pd.concat([mat2,t_cells.obs['LAG3SPOTS']], axis=1)
plot_threshold_comparison_save_pdf(mat2, col1='LAG3SPOTS', col2='TIM3', threshold1=1, directory=directory)

# loop
adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_phenotype.h5ad')
bdata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.h5ad')
adata.obs['phenotype_fine_majorcategories'] = bdata.obs['phenotype_fine_majorcategories'].values
# save dir
directory = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/2-Nature revisions- Jan2024/Figures/pdf_analysis_plots/LAG3 vs TIM3'
cell_intrest = ['CD103+ Tpex', 'CD103- Tpex', 'Tex','Tmem', 'Teff']
for i in cell_intrest:
    print(str(i))
    t_cells = adata[adata.obs['phenotype_fine_majorcategories'] == i]
    mat2 = pd.DataFrame(np.log1p(t_cells.X), index=t_cells.obs.index, columns=t_cells.var.index)
    mat2 = pd.concat([mat2,t_cells.obs['LAG3SPOTS']], axis=1)
    plot_threshold_comparison_save_pdf(mat2, col1='LAG3SPOTS', col2='TIM3', threshold1=2, directory=directory, filename = str(i) + '.pdf')
    

# check if CD45RO is present in CD103+ cells
cd103_pos_tpex = adata[adata.obs['phenotype_fine_majorcategories'] == 'CD103+ Tpex']
cd103_pos_tmem = adata[adata.obs['phenotype_fine_majorcategories'] == 'Tmem']

cd103_pos_tpex = pd.DataFrame(cd103_pos_tpex.X, index= cd103_pos_tpex.obs.index, columns=cd103_pos_tpex.var.index)
cd103_pos_tmem = pd.DataFrame(cd103_pos_tmem.X, index= cd103_pos_tmem.obs.index, columns=cd103_pos_tmem.var.index)

count_rows_above_value(df=cd103_pos_tpex, column_name='CD45RO', value=0.50) / len(cd103_pos_tpex) * 100
count_rows_above_value(df=cd103_pos_tmem, column_name='CD45RO', value=0.50) / len(cd103_pos_tmem) * 100


##############################################################################
# Breakdown of LAG3+ spots in a given population
##############################################################################


def plot_histogram_pandas_with_threshold(df, column, threshold, population_of_interest=None):
    """
    Plots a histogram of the distribution of values in a specified column using matplotlib,
    omitting rows with values below the threshold. Bars will have no space between them,
    and all numbers are listed on the x-axis.

    Parameters:
    - df: pandas.DataFrame
        The DataFrame containing the data.
    - column: str
        The name of the column to plot.
    - threshold: float
        The threshold value; rows with column values below this will be omitted.

    Returns:
    None
    """
    import matplotlib.ticker as ticker
    if column in df.columns:
        filtered_df = df[df[column] > threshold]
        plt.figure(figsize=(10, 6))
        # Plot histogram with bars touching each other
        plt.hist(filtered_df[column].dropna(), bins=None, edgecolor='k', alpha=0.7, rwidth=1)  # rwidth set to 1
        
        # Set x-axis to have a tick for every integer value within the range of the data
        plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        
        plt.xlabel(column)
        plt.grid(False)
        plt.ylabel('Frequency')
        if population_of_interest is not None: 
            plt.title(f'Histogram of {population_of_interest}')
        plt.xticks(rotation=90)  # Rotate x-axis labels for better readability if needed
        plt.show()
    else:
        print(f"Column '{column}' not found in DataFrame.")



populaiton = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/CD8_subclasses_cell_index.csv')
population_of_interest = ['TCF1+ PD1+ TIM3_LAG3+']
population_of_interest = ['TCF1- PD1+ TIM3_LAG3+ GZMBSPOTS-']
pop_intrest = populaiton[populaiton['Classification'].isin(population_of_interest)]
mat4 = adata[list(pop_intrest['Index'])].obs

#plot
plot_histogram_pandas_with_threshold(mat4, 'LAG3SPOTS', threshold=0,population_of_interest=population_of_interest)




def plot_histogram_for_category(data, category_column, category_value, value_column, bins=30, threshold=None, folder_path=''):
    """
    Plots a histogram for the distribution of values in a specific column, 
    filtered for a particular category in another column and values that pass a specified threshold, 
    and saves the plot as a PDF in a given folder.
    
    Parameters:
    - data: pandas DataFrame containing the data.
    - category_column: The name of the column to filter the category by.
    - category_value: The specific value of the category to filter for.
    - value_column: The name of the column for which to plot the distribution.
    - bins: The number of bins to use for the histogram. Default is 30.
    - threshold: The minimum value (inclusive) to include in the histogram. If None, all values are included.
    - folder_path: The path to the folder where the PDF file should be saved. If empty, the current directory is used.
    """
    # Filter the dataset for the category of interest
    filtered_data = data[data[category_column] == category_value]
    
    # Further filter the dataset if a threshold is specified
    if threshold is not None:
        filtered_data = filtered_data[filtered_data[value_column] >= threshold]
    
    # Ensure there is data to plot after filtering
    if filtered_data.empty:
        print(f"No data available for {category_value} in {category_column} meeting the threshold {threshold}.")
        return
    
    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(filtered_data[value_column], bins=bins, color='skyblue', edgecolor='black')
    plt.title(f'Distribution of {value_column} for {category_value} in {category_column}')
    plt.xlabel(value_column)
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Ensure folder path is available
    if not os.path.exists(folder_path) and folder_path != '':
        os.makedirs(folder_path)

    # Save the figure as a PDF in the specified folder
    filename = os.path.join(folder_path, f"{category_value.replace(' ', '_')}.pdf")  # Replace spaces with underscores for filename
    plt.savefig(filename, format='pdf')
    
    # Show the plot
    plt.show()
    print(f"Plot saved as {filename}")

# Example usage
# Ensure your data is correctly defined and then use the function like this:
# plot_histogram_for_category(data, 'YourCategoryColumn', 'YourCategoryValue', 'YourValueColumn', bins=30, threshold=YourThreshold, folder_path='YourFolderPath')

adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.h5ad')
data = adata.obs.copy()
# subset to cells of interest
cell_intrest = ['CD103- Tpex','CD103+ Tpex','Tmem','Teff','Tex']
d3 = data[data['phenotype_fine_majorcategories'].isin(cell_intrest)]
d3['phenotype_fine_majorcategories'] = d3['phenotype_fine_majorcategories'].astype('str')

folder_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/2-Nature revisions- Jan2024/Figures/pdf_analysis_plots/LAG3 spot distribution for each cell type'
for i in np.unique(d3['phenotype_fine_majorcategories']):
    plot_histogram_for_category(data=d3, category_column='phenotype_fine_majorcategories', category_value=i, value_column='LAG3SPOTS', bins=30, folder_path=folder_path)

data = adata.obs.copy()
plot_histogram_for_category(data, category_column='imageid', category_value='CHTNa-quantification_cleaned', value_column='LAG3SPOTS', bins=30, folder_path=folder_path, threshold =1)


folder_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/2-Nature revisions- Jan2024/Figures/pdf_analysis_plots/GZMB spot distribution for each cell type'
for i in np.unique(d3['phenotype_fine_majorcategories']):
    plot_histogram_for_category(data=d3, category_column='phenotype_fine_majorcategories', category_value=i, value_column='GZMBSPOTS', bins=30, folder_path=folder_path)


folder_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/2-Nature revisions- Jan2024/Figures/pdf_analysis_plots/MX1 spot distribution for each cell type'
for i in np.unique(d3['phenotype_fine_majorcategories']):
    plot_histogram_for_category(data=d3, category_column='phenotype_fine_majorcategories', category_value=i, value_column='MX1SPOTS', bins=30, folder_path=folder_path)


# t test spots

def t_test_spots(df, value_column, cat_column):
    """
    Performs pairwise t-tests between values in value_column grouped by the unique
    categories in cat_column using Pingouin's t-test.
    
    Parameters:
    - df: DataFrame containing the numeric values and categories.
    - value_column: String, name of the column with numeric data.
    - cat_column: String, name of the categorical column.
    
    Returns:
    - A DataFrame where rows and columns represent categories, and values are p-values
      from the t-tests.
    """
    # Identify all unique categories in the cat_column
    categories = df[cat_column].unique()
    
    # Prepare an empty DataFrame for the p-values
    p_values_df = pd.DataFrame(index=categories, columns=categories, dtype=float)
    
    # Generate all possible pairs of categories
    for category1, category2 in itertools.combinations(categories, 2):
        # Extract values for each category
        values_cat1 = df[df[cat_column] == category1][value_column].dropna()
        values_cat2 = df[df[cat_column] == category2][value_column].dropna()
        
        # Perform t-test between these values using Pingouin
        test_result = pg.ttest(values_cat1, values_cat2, paired=False)
        
        # Store p-value in the DataFrame
        p_val = test_result['p-val'].values[0]
        p_values_df.at[category1, category2] = p_val
        p_values_df.at[category2, category1] = p_val  # Ensure the matrix is symmetric

    # Correctly fill the diagonal with NaNs
    np.fill_diagonal(p_values_df.values, np.nan)
    
    return p_values_df



t_test = t_test_spots (df= adata.obs, value_column='MX1SPOTS', cat_column='phenotype_fine_majorcategories')
t_test.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/subclasses_results/MX1SPOTS_TTEST.csv',decimal='.')


##############################################################################
# spatial analysis of finer categories
##############################################################################

# spatial scatterplot
sm.pl.spatial_scatterPlot (adata, 
                     colorBy = ['phenotype_fine'], 
                     x_coordinate='CENTROID_X', y_coordinate='CENTROID_Y',
                     figsize=(4,4),
                     s=2,
                     customColors=customColors,
                     plotLegend=True,
                     fontsize=3)

plotly (adata,phenotype='phenotype_finer',x='CENTROID_X',y='CENTROID_Y',size=5)
plotly (adata,phenotype='phenotype_fine',x='CENTROID_X',y='CENTROID_Y',size=5)



# simple spatial analysis
adata = sm.tl.spatial_distance (adata, 
                               x_coordinate='CENTROID_X', y_coordinate='CENTROID_Y', 
                               z_coordinate='CENTROID_Z', 
                               phenotype='phenotype_fine_majorcategories', 
                               subset=None, 
                               imageid='imageid', 
                               label='spatial_distance')
# plot
distance_to = ['CD103+ Tpex', 'CD103- Tpex', 'Tex', 'Tmem', 'Teff']
sm.pl.spatial_distance (adata, method='numeric',distance_from='Tumor', distance_to=distance_to, phenotype='phenotype_fine_majorcategories', log=True)
sm.pl.spatial_distance (adata, method='numeric',distance_from='CD103+ Tpex', distance_to=distance_to, phenotype='phenotype_fine_majorcategories', log=True)
sm.pl.spatial_distance (adata, method='numeric',distance_from='CD103- Tpex', distance_to=distance_to, phenotype='phenotype_fine_majorcategories', log=True)


# fine categories
adata = sm.tl.spatial_distance (adata, 
                               x_coordinate='CENTROID_X', y_coordinate='CENTROID_Y', 
                               z_coordinate='CENTROID_Z', 
                               phenotype='phenotype_fine', 
                               subset=None, 
                               imageid='imageid', 
                               label='spatial_distance')




distance_to = ['TCF1- CD103+ T (naive)','TCF1- CD103+  T (activated)','TCF1- CD103+  T (Tp)', 'TCF1- CD103+  T (Tex)']
sm.pl.spatial_distance (adata, method='numeric',distance_from='TCF1+ CD103+  T (Tp)', distance_to=distance_to, phenotype='phenotype_fine', log=True)

distance_to = ['TCF1- CD103- T (naive)','TCF1- CD103-  T (activated)','TCF1- CD103-  T (Tp)', 'TCF1- CD103-  T (Tex)']
sm.pl.spatial_distance (adata, method='numeric',distance_from='TCF1+ CD103-  T (Tp)', distance_to=distance_to, phenotype='phenotype_fine', log=True)

np.unique(adata.obs['phenotype_fine'])

# t test of distances

def compare_categories_by_df2(df1, df2, cat_column):
    """
    Performs pairwise t-tests between columns in df1, where the columns to compare
    are defined by the unique categories in the cat_column of df2.
    
    Parameters:
    - df1: DataFrame where each column contains numeric data for a category.
    - df2: DataFrame that contains the cat_column specifying categories.
    - cat_column: String, name of the categorical column in df2.
    
    Returns:
    - A DataFrame where rows and columns represent categories, and values are p-values
      from the t-tests.
    """
    # Identify all unique categories in the cat_column of df2
    categories = df2[cat_column].unique()
    
    # Prepare an empty DataFrame for the p-values
    p_values_df = pd.DataFrame(index=categories, columns=categories, dtype=float)
    
    # Generate all possible pairs of categories
    for category1, category2 in itertools.combinations(categories, 2):
        # Perform t-test between the corresponding columns in df1
        test_result = pg.ttest(df1[category1].dropna(), df1[category2].dropna(), paired=False)
        
        # Store p-value in the DataFrame
        p_values_df.at[category1, category2] = test_result['p-val'].values[0]
        p_values_df.at[category2, category1] = test_result['p-val'].values[0]  # Ensure the matrix is symmetric

    # Correctly fill the diagonal with NaNs
    np.fill_diagonal(p_values_df.values, np.nan)
    
    return p_values_df


# test
r1 = compare_categories_by_df2(df1=adata.uns['spatial_distance'], df2=adata.obs, cat_column='phenotype_fine')
r1 = compare_categories_by_df2(df1=adata.uns['spatial_distance'], df2=adata.obs, cat_column='phenotype_fine')



adata[adata.obs['phenotype_fine_majorcategories'] == 'Tex'].obs['MX1SPOTS'].median()

# sptial interaction
adata = sm.tl.spatial_interaction(adata, 
                                  method='knn', 
                                  x_coordinate='CENTROID_X', y_coordinate='CENTROID_Y', 
                                  z_coordinate='CENTROID_Z', 
                                  phenotype='phenotype_fine_majorcategories', 
                                  knn=5, 
                                  permutation=1000,
                                  label='spatial_interaction_knn')

# plot
sm.pl.spatial_interaction(adata, 
                          summarize_plot=True, 
                          binary_view=False,
                          spatial_interaction='spatial_interaction_knn',
                          row_cluster=False, linewidths=0.75, linecolor='black')



# look at the distribution of cetains elements between cell types
import itertools
data = adata.obs.copy()
data['PCNA_KI67'] = list(itertools.chain(*adata[:, ["PCNA_KI67"]].X))
 
 
# subset to cells of interest
cell_intrest = ['CD103- Tpex','CD103+ Tpex','Tmem','Teff','Tex'] 
d3 = data[data['phenotype_fine_majorcategories'].isin(cell_intrest)]
d3['phenotype_fine_majorcategories'] = d3['phenotype_fine_majorcategories'].astype('str')

                                    
def plot_grouped_boxenplot_ordered(data, column1, column2):
    """
    Groups the dataframe by column1 and plots boxen plots for the distribution of column2,
    with x-axis ordered by the median of column2 in each group. Additionally, annotates
    the plot with the number of observations in each category of column1, positioned
    at the top of the plot.
    
    Parameters:
    - data: pandas DataFrame containing the data.
    - column1: The name of the column to group the data by.
    - column2: The name of the column for which to plot the distribution.
    """
    # Calculate medians and order the groups
    order = data.groupby(column1)[column2].mean().sort_values().index
    
    # Calculate the number of observations in each group
    counts = data[column1].value_counts().reindex(order)
    
    # Creating the boxenplot with ordered groups
    plt.figure(figsize=(10, 6))
    sns.boxenplot(x=column1, y=column2, data=data, order=order)
    
    # Determine the top of the plot for placing count annotations
    y_max = data[column2].max()
    plt.ylim(0, y_max + y_max * 0.1)  # Adjust y-limit to make space for annotations
    
    # Annotating the plot with the number of observations
    for i, count in enumerate(counts):
        plt.text(i, y_max + y_max * 0.05, f'n={count}', ha='center', va='bottom')
    
    plt.xlabel(column1)
    plt.ylabel(column2)
    plt.title(f'Distribution of {column2} by {column1}')
    plt.xticks(rotation=90)  # Rotates the labels on the x-axis for better readability
    plt.show()


# plot
plot_grouped_boxenplot_ordered (data=d3, column1='phenotype_fine_majorcategories', column2='MX1SPOTS') #PCNA_KI67 #MX1SPOTS

# stacked barplot to show% postive cells for a given column (spots colun)
def plot_percent_stacked_barplot_with_annotations_ordered(data, column1, column2, threshold):
    """
    Creates a percent stacked barplot for groups in column1, classifying each observation
    as 'Positive' or 'Negative' based on a threshold in column2. The barplot shows
    the percentage of 'Positive' and 'Negative' in each group, with annotations indicating
    the number of 'Positive' observations. Groups are ordered by the median value of column2.
    
    Parameters:
    - data: pandas DataFrame containing the data.
    - column1: The name of the column to group the data by.
    - column2: The name of the column to classify as 'Positive' or 'Negative'.
    - threshold: The value used to classify observations in column2.
    """
    # Classify each observation
    data['Classification'] = ['Positive' if x > threshold else 'Negative' for x in data[column2]]
    
    # Calculate medians to order the groups
    median_order = data.groupby(column1)[column2].mean().sort_values().index
    
    # Calculate the percentage of each classification within each group
    classification_counts = data.groupby([column1, 'Classification']).size().unstack(fill_value=0)
    classification_percentages = classification_counts.div(classification_counts.sum(axis=1), axis=0) * 100
    
    # Reorder based on the median values
    classification_percentages = classification_percentages.reindex(median_order)
    classification_counts = classification_counts.reindex(median_order)
    
    # Plotting
    ax = classification_percentages.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='viridis', width=0.8)
    plt.ylabel('Percentage')
    plt.title(f'Percentage of Positive and Negative Classifications in {column1} based on {column2}')
    plt.legend(title='Classification', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation='vertical')  # Set x-axis labels to be vertical
    plt.tight_layout()  # Adjust layout to not cut off labels

    # Annotate with the number of 'Positive' cells
    for i, (index, row) in enumerate(classification_counts.iterrows()):
        pos_count = row.get('Positive', 0)  # Safely get the count, defaulting to 0 if not present
        if pos_count > 0:
            # Annotate above each bar
            ax.annotate(f'n={pos_count}', (i, 100), textcoords="offset points", xytext=(0,10), ha='center')

    plt.show()


# Example usage
plot_percent_stacked_barplot_with_annotations_ordered (data=d3, column1='phenotype_fine_majorcategories', column2='GZMBSPOTS', threshold=1)


# neighbourhood analysis
adata = sm.tl.spatial_count (adata,x_coordinate='CENTROID_X',
                             y_coordinate='CENTROID_Y',
                             z_coordinate='CENTROID_Z',
                             phenotype='phenotype_fine',method='knn',
                             knn=5, imageid='imageid',
                             subset=None,label='spatial_count_knn')

adata = sm.tl.spatial_cluster (adata, k= 10, method = 'kmeans', df_name='spatial_count_knn') # results will be saved under adata.obs['spatial_kmeans']

# view on image
plotly (adata,phenotype='spatial_kmeans',x='CENTROID_X',y='CENTROID_Y',size=5)
plotly (adata,phenotype='phenotype_fine',x='CENTROID_X',y='CENTROID_Y',size=5)
plotly (adata,phenotype='phenotype_fine_majorcategories',x='CENTROID_X',y='CENTROID_Y',size=5)




# stacked barplot
sm.pl.stacked_barplot (adata,x_axis='spatial_kmeans',y_axis='phenotype_fine',
                     method='percent',plot_tool='plotly',
                     color_discrete_sequence=px.colors.qualitative.Alphabet)

# save
adata.write('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.h5ad')



# spatial scattter plot
np.unique(adata.obs['phenotype_fine_majorcategories'])


customColors = { 'Unknown' : '#eff1ed',
                'Treg' : '#eff1ed',
                'TCF1- CD103- T _naive_' : '#eff1ed',
                'TCF1- CD103-  T _activated_' : '#eff1ed',
                'TCF1+ CD103- T _naive_': '#eff1ed',
                'TCF1+ CD103-  T _activated_' : '#eff1ed',
                'TCF1+ CD103+ T _naive_' : '#eff1ed',
                'TCF1+ CD103+  T _activated_' : '#eff1ed',
                'Myeloid cells' : '#eff1ed',
                'Immune' : '#eff1ed',
                'CD4 T' : '#eff1ed',

                'CD103+ Tpex' : '#e05abf',
                'CD103- Tpex' : '#f6aa1c',
                'Teff' : '#bc4749',
                'Tex' : '#708d81',
                'Tmem' : '#6b5e62',
                'Tumor' : '#d6d6d6',
                
    }

spatial_scatterPlot (adata, 
                colorBy = ['phenotype_fine_majorcategories'],
                topLayer = ['CD103+ Tpex','CD103- Tpex','Teff','Tex','Tmem'],
                x_coordinate='CENTROID_X', y_coordinate='CENTROID_Y',
                 figsize=(5,5),
                 s=1,
                 plotLegend=True,
                 fontsize=3,
                 dpi=350,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 outputFileName='Tcells_location.png',
                 outputDir='/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/2-Nature revisions- Jan2024/Figures/pdf_analysis_plots'
                 )


# Identify some example index of Tpex and Tex cells
adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.h5ad')
adata.obs['phenotype_fine'].value_counts()
data = sm.hl.scimap_to_csv (adata, data_type='scaled')
data.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/CHTN_analysis/with_finer_phenotype.csv', index=False)
adata


#uMAP
# remove unknown cells
adata = adata[adata.obs['phenotype_broad'] != 'Unknown']
adata = sm.tl.umap(adata,  use_raw=True, log=True)
sm.pl.umap(adata, color='phenotype_broad', s=0.01)

# cluster with leiden
adata = sm.tl.cluster(adata, method='leiden',use_raw=True, log=True)

sm.pl.umap(adata, color='leiden', s=0.1)





























