#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 17:57:12 2023
@author: UMAP of 3D imaging data
"""


import matplotlib
import scimap as sm
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from PyComplexHeatmap import *

%matplotlib qt
sc.set_figure_params(scanpy=True, fontsize=14,  figsize=(15,10))
matplotlib.rcParams['pdf.fonttype'] = 42 



##############################################################################
# combined analysis
##############################################################################

x = adata[adata.obs['phenotype'] == 'T cells']
x.obs['imageid'].value_counts()


# read or write adata
adata = ad.read("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_INV.h5ad")
adata = ad.read("C:/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_INV_copy.h5ad")

# write
adata.write("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_INV.h5ad")

adata.write(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\MIS_INV.h5ad")

# convert both MIS and INV to equal columns
MIS = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification5.csv'
IM = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iic-quantificationV3.csv'

MIS = 'C:/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification4.csv'
IM = 'C:/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iic-quantificationV2.csv'


MIS = pd.read_csv(MIS)
IM = pd.read_csv(IM)
IM.rename(columns={'GZMB_spots': 'GZMB_SPOTS'}, inplace=True)
IM.rename(columns={'LAG3_spots': 'LAG3_SPOTS'}, inplace=True)
# Find the common columns between the two data frames
common_columns = MIS.columns.intersection(IM.columns)

# fix sox9
MIS['SOX9'] = MIS.apply(lambda row: 0 if row['SOX9 for gating'] == 0 else row['SOX9'], axis=1)

# Create a list of columns excluding 'LAG3_SPOTS' and 'GZMB_SPOTS'
columns_to_keep = [col for col in common_columns if col not in ['LAG3_SPOTS', 'GZMB_SPOTS']]
columns_to_keep.extend(['LAG3_SPOTS', 'GZMB_SPOTS'])

# deop SOX9
#common_columns = common_columns.drop(['SOX9'])
# Subset both data frames based on the common columns
MIS = MIS[columns_to_keep]
IM = IM[columns_to_keep]

MIS.insert(0, 'CellID', MIS.index) # add cellid
IM.insert(0, 'CellID', IM.index) # add cellid

MIS.dropna(inplace=True)
IM.dropna(inplace=True)




# write out CSV
MIS.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME.csv', index=False) 
IM.to_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/INV_PRAME.csv', index=False) 


# create adata object
feature_table_path = ['/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/MIS_PRAME.csv',
                      '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/INV_PRAME.csv']

adata = sm.pp.mcmicro_to_scimap (feature_table_path, log=False, split='volume')

adata.obs['phenotype'].value_counts()

# scale and phenotype the cells
manual_gate = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/manual_gates.csv')
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
adata[adata.obs['phenotype'] == 'T cells'].obs['imageid'].value_counts()
# rename cells
rename= {'CD4 T': ['T cells']}
adata = sm.hl.rename(adata, rename, from_column='phenotype', to_column='phenotype')

# add spots into phenotyping
# Define a function to annotate new values
def annotate_row(row):
    if row['phenotype'] == 'CD8 T' and row['LAG3_SPOTS'] > 0:
        return 'LAG3+ CD8 T'
    elif row['phenotype'] == 'PD1+ T' and row['LAG3_SPOTS'] > 0:
        return 'Exhausted T'
    else:
        return row['phenotype']

df = adata.obs.copy()
adata.obs['phenotype_1'] = df.apply(annotate_row, axis=1)

# percent of CD103+ cells that were PD1/ LAG3 +
adata = sm.hl.classify(adata, pos=['PD1'],  phenotype='phenotype',  classify_label='PD1+ Tissue T', subclassify_phenotype=['Tissue T'], label='classify_tissueT',collapse_failed=False)
adata.obs['phenotype'].value_counts()
adata.obs['classify_tissueT'].value_counts()
# PD1+ Tissue T           7570
# Tissue T                1047
# find out how many of these are LAG3+?
def annotate_row(row):
    if row['classify_tissueT'] == 'Tissue T' and row['LAG3_SPOTS'] > 0:
        return 'LAG3+ Tissue T'
    else:
        return row['classify_tissueT']
    
df = adata.obs.copy()
adata.obs['subclassify_tissueT'] = df.apply(annotate_row, axis=1)
adata.obs['subclassify_tissueT'].value_counts()
# Tissue T                 618
# LAG3+ Tissue T           429
# total tissue T that are PD1 or LAG3 +
((8617-618) / 8617) * 100 # -> 93%

# Check if CD103- cells are PD1+
adata = sm.hl.classify(adata, pos=['PD1'],  phenotype='phenotype',  classify_label='PD1+ T', subclassify_phenotype=['CD8 T'], label='classify_CD8T',collapse_failed=False)
adata.obs['classify_CD8T'].value_counts()
# total CDT that is PD1
(110 / 1423) * 100 # -> 8%

# How many of TPEX are GZMB?
rename= {'Tpex': ['LAG3+ Tissue T', 'PD1+ Tissue T']}
adata = sm.hl.rename (adata, rename, from_column='subclassify_tissueT', to_column='tpex')
adata.obs['tpex'].value_counts()
def annotate_row(row):
    if row['tpex'] == 'Tpex' and row['GZMB_SPOTS'] > 0:
        return 'GZMB+ Tpex'
    else:
        return row['tpex']
df = adata.obs.copy()
adata.obs['GZMB+ Tpex'] = df.apply(annotate_row, axis=1)
adata.obs['GZMB+ Tpex'].value_counts()

(1568 / 7999) * 100


# KI67+ Tpex
adata.obs['tpex'].value_counts()
adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='tpex',  classify_label='Ki67+ Tpex', subclassify_phenotype=['Tpex'], label='classify_Tpex',collapse_failed=False)
adata.obs['classify_Tpex'].value_counts()
adata.obs['classify_tissueT'].value_counts()

inv = adata[adata.obs['imageid'] == 'INV_PRAME']
inv.obs['phenotype'].value_counts()

#KI67+ GZMB + Tpex
adata.obs['GZMB+ Tpex'].value_counts()
adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='GZMB+ Tpex',  classify_label='Ki67+ GZMB+ Tpex', subclassify_phenotype=['GZMB+ Tpex'], label='classify_Tpex_gzmb',collapse_failed=False)
adata.obs['classify_Tpex_gzmb'].value_counts()

#

# no of
# total Tissue T -> 8840
# total Tpex -> 8275
# Tpex + ki67 -> 706
# Tpex + GZMB -> 1612
# TPex + GZMB + ki67 -> 243
# classically exhausted cells
# CD8 T -> 1100
# T cells + PD1 -> 367
# T cells + PD1 + LAG3 -> 138

# identify total number of clasically exhausted cellls
adata.obs['phenotype'].value_counts()
adata = sm.hl.classify(adata, pos=['PD1'],  phenotype='phenotype',  classify_label='CD8+PD1+', subclassify_phenotype=['CD8 T'], label='CD8+PD1+',collapse_failed=False)
adata.obs['CD8+PD1+'].value_counts()

def annotate_row(row):
    if row['CD8+PD1+'] == 'CD8+PD1+' and row['LAG3_SPOTS'] > 0:
        return 'CD8+PD1+LAG3+'
    else:
        return row['CD8+PD1+']
df = adata.obs.copy()
adata.obs['CD8+PD1+LAG3+'] = df.apply(annotate_row, axis=1)
adata.obs['CD8+PD1+LAG3+'].value_counts()


# CD4+ Ki67+ cells


# do the same for CD4 T cells
adata.obs['phenotype'].value_counts()
#pd1+ cd4 t
adata = sm.hl.classify(adata, pos=['PD1'],  phenotype='phenotype',  classify_label='PD1+CD4T', subclassify_phenotype=['CD4 T'], label='classify_CD4T',collapse_failed=False)
adata.obs['classify_CD4T'].value_counts()

adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='classify_CD4T',  classify_label='PD1+KI67+CD4T', subclassify_phenotype=['PD1+CD4T'], label='classify_CD4T_ki67',collapse_failed=False)
adata.obs['classify_CD4T_ki67'].value_counts()

def annotate_row(row):
    if row['classify_CD4T'] == 'PD1+CD4T' and row['GZMB_SPOTS'] > 0:
        return 'CD4T+PD1+GZMB'
    else:
        return row['classify_CD4T']
df = adata.obs.copy()
adata.obs['CD4 T + PD1+ GZMB'] = df.apply(annotate_row, axis=1)
adata.obs['CD4 T + PD1+ GZMB'].value_counts()

adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='CD4 T + PD1+ GZMB',  classify_label='PD1+KI67+CD4T+GZMB', subclassify_phenotype=['CD4T+PD1+GZMB'], label='classify_CD4T_ki67_gzmb',collapse_failed=False)
adata.obs['classify_CD4T_ki67_gzmb'].value_counts()

########################
# do same without PD1
adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='phenotype',  classify_label='ki67+CD4T', subclassify_phenotype=['CD4 T'], label='classify_CD4T_ki67',collapse_failed=False)
adata.obs['classify_CD4T_ki67'].value_counts()

def annotate_row(row):
    if row['phenotype'] == 'CD4 T' and row['GZMB_SPOTS'] > 0:
        return 'CD4T+GZMB'
    else:
        return row['phenotype']
df = adata.obs.copy()
adata.obs['CD4T+GZMB'] = df.apply(annotate_row, axis=1)
adata.obs['CD4T+GZMB'].value_counts()

adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='CD4T+GZMB',  classify_label='KI67+CD4T+GZMB', subclassify_phenotype=['CD4T+GZMB'], label='classify_CD4T_ki67_gzmb',collapse_failed=False)
adata.obs['classify_CD4T_ki67_gzmb'].value_counts()


# CD4+ T cells -> 1285
# CD4 T PD1+ -> 712
# CD4 T + PD1+ KI67 -> 43
# CD4 T + PD1+ GZMB -> 57
# CD4 T + PD1+ GZMB + KI67 -> 3

# CD4 T + KI67 -> 72
# CD4 T + GZMB -> 93
# CD4 T + GZMB + KI67 -> 5


# identify  cells close to melanocytes and then look at proportion of tpex
# How many of TPEX are GZMB?
# percent of CD103+ cells that were PD1/ LAG3 +
bdata = sm.hl.classify(bdata, pos=['PD1'],  phenotype='phenotype',  classify_label='PD1+ Tissue T', subclassify_phenotype=['Tissue T'], label='classify_tissueT',collapse_failed=False)
bdata.obs['phenotype'].value_counts()
bdata.obs['classify_tissueT'].value_counts()
# PD1+ Tissue T           2450
# Tissue T                2739
# find out how many of these are LAG3+?
def annotate_row(row):
    if row['classify_tissueT'] == 'Tissue T' and row['LAG3_SPOTS'] > 0:
        return 'LAG3+ Tissue T'
    else:
        return row['classify_tissueT']
    
df = bdata.obs.copy()
bdata.obs['subclassify_tissueT'] = df.apply(annotate_row, axis=1)
bdata.obs['subclassify_tissueT'].value_counts()
# LAG3+ Tissue T           141


# How many of TPEX are GZMB?
rename= {'Tpex': ['LAG3+ Tissue T', 'PD1+ Tissue T']}
bdata = sm.hl.rename (bdata, rename, from_column='subclassify_tissueT', to_column='tpex')
bdata.obs['tpex'].value_counts()
def annotate_row(row):
    if row['tpex'] == 'Tpex' and row['GZMB_SPOTS'] > 0:
        return 'GZMB+ Tpex'
    else:
        return row['tpex']
df = bdata.obs.copy()
bdata.obs['GZMB+ Tpex'] = df.apply(annotate_row, axis=1)
bdata.obs['GZMB+ Tpex'].value_counts()


# KI67+ Tpex
bdata.obs['tpex'].value_counts()
bdata = sm.hl.classify(bdata, pos=['Ki67'],  phenotype='tpex',  classify_label='Ki67+ Tpex', subclassify_phenotype=['Tpex'], label='classify_Tpex',collapse_failed=False)
bdata.obs['classify_Tpex'].value_counts()

inv = bdata[bdata.obs['imageid'] == 'INV_PRAME']
inv.obs['phenotype'].value_counts()

#KI67+ GZMB + Tpex
bdata = sm.hl.classify(bdata, pos=['Ki67'],  phenotype='GZMB+ Tpex',  classify_label='Ki67+ GZMB+ Tpex', subclassify_phenotype=['GZMB+ Tpex'], label='classify_Tpex_gzmb',collapse_failed=False)
bdata.obs['classify_Tpex_gzmb'].value_counts()

# no of
# total Tissue T -> 2739
# total Tpex -> 2591
# Tpex + ki67 -> 259
# Tpex + GZMB -> 726
# TPex + GZMB + ki67 -> 116



adata.write(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\10-30-2023-OFFICECP.h5ad")


# split MIS and INV
bdata = adata[adata.obs['imageid'] == 'INV']

from sklearn.neighbors import BallTree
data = pd.DataFrame({'x': bdata.obs['X_centroid'], 'y': bdata.obs['Y_centroid'], 'z': bdata.obs['Z_centroid'], 'phenotype': bdata.obs['phenotype']})

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

##
n = pd.DataFrame(tumor_cells.stack(), columns = ["neighbour_phenotype"])
n.index = n.index.get_level_values(0) # Drop the multi index
n = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
n['cell'] = n.index
# do same for indexes
nn = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
nn.index = nn.index.get_level_values(0) # Drop the multi index
nn = nn.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
nn['cell'] = nn.index

# keep only tissue T cell neighbours
tissue_t = n[n['neighbour_phenotype'] == 'Tissue T']
tissue_t_index = nn.loc[tissue_t.index]

# extract all indexes and plot
# Find unique elements in 'neighbour_phenotype' and 'cell' columns
unique_neighbour_phenotype = tissue_t_index['neighbour_phenotype'].unique().tolist()
unique_cell = tissue_t_index['cell'].unique().tolist()
# Combine the unique elements into a single list
combined_unique_elements = list(set(unique_neighbour_phenotype + unique_cell))

# add to data
# Convert the list to a set for faster look-up
combined_unique_elements_set = set(combined_unique_elements)
# Add 'tumor_neighborhood' column based on the condition
data['tumor_neighborhood'] = ['neighbour' if idx in combined_unique_elements_set else 'other' for idx in data.index]
bdata.obs['tumor_neighborhood'] = data['tumor_neighborhood']

# MIS cells (note have to run the above scriopt twice with varying bdata to get this)
mis_cells = combined_unique_elements_set.copy()
inv_cells = combined_unique_elements_set.copy()

# cobine 
neigh_cells = mis_cells | inv_cells
adata.obs['tumor_neighborhood'] = ['neighbour' if idx in neigh_cells else 'other' for idx in adata.obs.index]

# scattter plot
customColors = { 'other' : '#e5e5e5',
                    'neighbour' : '#073b4c',}
sm.pl.spatial_scatterPlot (adata, 
                     colorBy = ['tumor_neighborhood'], 
                     subset = 'INV',
                     figsize=(4,4),
                     customColors=customColors,
                     s=0.5)


# RERUN all calculations with respect to tpex
bdata = adata[adata.obs['tumor_neighborhood'] == 'neighbour']




customColors = { 'Tumor' : '#073b4c',
                'Tissue T' : '#e5e5e5',
                'Unknown' : '#e5e5e5',
                'Dendritic cells' : '#e5e5e5',
                'Macrophage' : '#e5e5e5',
                'endothelial' : '#e5e5e5',
                'CD4 T' : '#e5e5e5',
                'CD8 T' : '#e5e5e5',
                'keratinocytes' : '#e5e5e5',
                'T reg' : '#e5e5e5',
                'CD11B+ CD11C- cells' : '#e5e5e5',
                'B cells' : '#e5e5e5',
                'Langerhan cells' : '#e5e5e5',}
sm.pl.spatial_scatterPlot (adata, 
                     colorBy = ['phenotype'], 
                     subset = 'MIS',
                     figsize=(4,4),
                     customColors=customColors,
                     s=0.5)

n_freq = n.groupby(['cell','neighbour_phenotype']).size().unstack().fillna(0)

# subset the neighbourhoods of intrest
n_freq_subset = n_freq.loc[modified_list]
n_freq_subset_percentage = n_freq_subset.div(n_freq_subset.sum(axis=1), axis=0) * 100

# sort by CD8
sort_column = 'CD8 T'
n_freq_subset_percentage = n_freq_subset_percentage.sort_values(by=sort_column, ascending=True)
n_freq_subset_percentage = pd.concat([n_freq_subset_percentage[sort_column], n_freq_subset_percentage.drop(sort_column, axis=1)], axis=1)

# plot them
# Plot the stacked bar plot
ax = n_freq_subset_percentage.plot(kind='bar', stacked=True, figsize=(20, 10), width=1, grid=False)
ax.legend(title="Category", bbox_to_anchor=(1.04, 1), loc='upper left')
plt.tight_layout()
plt.show()

# remove tumor cells
columns_to_drop = ['Unknown', 'Tumor']
n_freq_subsete_dropped = n_freq_subset.drop(columns=columns_to_drop, inplace=False)
n_freq_subset_percentage_dropped = n_freq_subsete_dropped.div(n_freq_subsete_dropped.sum(axis=1), axis=0) * 100
n_freq_subset_percentage_dropped = n_freq_subset_percentage_dropped.sort_values(by=sort_column, ascending=True)
n_freq_subset_percentage_dropped = pd.concat([n_freq_subset_percentage_dropped[sort_column], n_freq_subset_percentage_dropped.drop(sort_column, axis=1)], axis=1)


# plot them
ax = n_freq_subset_percentage_dropped.plot(kind='bar', stacked=True, figsize=(20, 10), width=1, grid=False)
ax.legend(title="Category", bbox_to_anchor=(1.04, 1), loc='upper left')
plt.tight_layout()
plt.show()




# proportion of KI67+ T cells
rename= {'T': ['Tissue T', 'CD8 T', 'CD4 T', 'T cells', 'T reg']}
adata = sm.hl.rename (adata, rename, from_column='phenotype', to_column='T')
adata = sm.hl.classify(adata, pos=['Ki67'],  phenotype='T',  classify_label='Ki67 T', subclassify_phenotype=['T'], label='classify_T',collapse_failed=False)
adata.obs['classify_T'].value_counts()
adata.obs['T'].value_counts()

bdata = adata[adata.obs['classify_T'] == 'Ki67 T']
len(bdata[bdata.obs['imageid'] == 'INV_PRAME']) / len(bdata[bdata.obs['imageid'] == 'MIS_PRAME'])


# export phenotyping
bdata = adata.copy()
bdata = bdata[bdata.obs['imageid'] == 'MIS']
exportPhenotype = pd.DataFrame({'CellID': bdata.obs['CellID'], 'phenotype': bdata.obs['phenotype']}) 
exportPhenotype.to_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\phenotypes.csv",index=False)
exportPhenotype.to_csv("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/phenotypes_MIS.csv",index=False)


# remove the langer han cells from INV
condition = ~((adata.obs['imageid'] == 'INV') & (adata.obs['phenotype'] == 'Langerhan cells'))
adata = adata[condition]

# plot the proportion of immune cells between the two regions
plotly (adata,phenotype='phenotype',image_id='INV',x='X_centroid',y='Y_centroid',size=8)


y_axis= ['Tissue T', 'CD8 T', 'LAG3+ CD8 T', 'T reg', 'CD4 T',  'B cells', 
         'Langerhan cells','Dendritic cells', 'Macrophage','CD11B+ CD11C- cells',
         'endothelial',
         'keratinocytes', 'Tumor'
         ]
x_axis=['MIS', 'INV']
# plot
sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='phenotype_1',subset_yaxis=y_axis, order_xaxis=x_axis, order_yaxis=y_axis,color=customColors)

sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='phenotype',subset_yaxis=y_axis, order_xaxis=x_axis, order_yaxis=y_axis,color=customColors,method='absolute')

sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='phenotype',method='absolute')

adata[adata.obs['imageid'] == 'INV'].obs['phenotype'].value_counts()


x = sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='phenotype_1',return_data = True)


# remove unknown cells
adata = adata[adata.obs['phenotype'] != 'Unknown']

# scatter plot
customColors = { 'Tissue T' : '#272323',
                'CD8 T' : '#457b9d',
                'T reg' : '#2a9d8f',
                'CD4 T' : '#ce4257',
                'B cells' : '#390099',
                'Langerhan cells' : '#d8e2dc',
                'Dendritic cells' : '#ffb703',    
                'Macrophage' : '#0c5149', 
                'CD11B+ CD11C- cells' : '#90e0ef',
                'endothelial' : '#a3b18a',
                'keratinocytes' : '#99582a',
                'Tumor' : '#faedcd',
                'Unknown': '#e9ecef',
                'LAG3+ CD8 T': '#7209b7'
    }

sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype_1'], 
                 subset = 'MIS',
                 figsize=(8,6),
                 s=4,
                 customColors=customColors,
                 plotLegend=True,
                 fontsize=7,
                 #outputDir=r'C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\figures_raw'
                 outputDir='/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/figures_raw')


# rename for major CD4 and CD8
adata.obs['phenotype'].value_counts()
rename= {'CD8 T': ['Tissue T']}
adata = sm.hl.rename (adata, rename, from_column='phenotype', to_column='major_cd8_phenotype')


# Granzyme expression
df = adata.obs.copy()
adata.obs['GZMB_classification'] = df.apply(lambda row: str(row['major_cd8_phenotype']) + '_GZMB+' if row['GZMB_SPOTS'] > 0 else str(row['major_cd8_phenotype']) + '_GZMB-', axis=1)

bdata = adata[adata.obs['imageid'] == 'MIS']
adata.obs['GZMB_classification'].value_counts()
# CD4 T_GZMB-                   205
# CD4 T_GZMB+                    13
# MIS and INV
# CD4 T_GZMB-                   1128
# CD4 T_GZMB+                     86

# plot
y_axis= ['CD8 T_GZMB+', 'CD4 T_GZMB+',]
x_axis=['MIS', 'INV']
sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='GZMB_classification',subset_yaxis=y_axis, order_xaxis=x_axis, order_yaxis=y_axis)
A
# find the percentage of T_GZMB in all T cells
mis = adata[adata.obs['imageid'] == 'MIS']
l1 = len(mis[mis.obs['phenotype_1'].isin([
                                            'Tissue T', 
                                            'CD8 T', 
                                            'LAG3+ CD8 T'
                                            ])])

l2 = len(mis[mis.obs['phenotype_1'].isin([
                                            'Tissue T', 
                                            #'CD8 T', 
                                            #'LAG3+ CD8 T'
                                            ])])

l2 = len(mis[mis.obs['GZMB_classification'].isin([
                                                  'Tissue T_GZMB+', 
                                                  #'CD8 T_GZMB+', 
                                                  #'LAG3+ CD8 T_GZMB+'
                                                  ])])
l2/l1 * 100

# spatial analysis


x = adata[adata.obs['imageid'] == 'INV']
x[x.obs['phenotype'] == 'CD4 T']
x[x.obs['GZMB_classification'] == 'CD8 T_GZMB+']
x[x.obs['GZMB_classification'] == 'CD4 T_GZMB+']


y = adata.obs[['GZMB_classification','GZMB_SPOTS']]



# perform spatial neighbourhood analysis
adata = spatial_count (adata,   x_coordinate='X_centroid',y_coordinate='Y_centroid',z_coordinate='Z_centroid', method = 'knn', knn=5)

adata = sm.tl.spatial_cluster (adata, df_name='spatial_count', k= 10, method = 'kmeans') # results will be saved under adata.obs['spatial_kmeans']


subset_yaxis = ['B cells', 'CD11B+ CD11C- cells', 'CD4 T', 'CD8 T',
       'Dendritic cells', 'Langerhan cells', 'Macrophage', 'T reg',
       'Tissue T', 'Tumor', 'endothelial', 'keratinocytes']

sm.pl.stacked_barplot (adata,x_axis='spatial_kmeans',y_axis='phenotype', plot_tool='plotly',subset_yaxis=subset_yaxis)
sm.pl.stacked_barplot (adata,x_axis='spatial_kmeans_consolidated',y_axis='phenotype',subset_yaxis=subset_yaxis)

plotly (adata,phenotype='spatial_kmeans_consolidated',image_id='MIS',x='X_centroid',y='Y_centroid',size=7)
plotly (adata,phenotype='phenotype',image_id='MIS',x='X_centroid',y='Y_centroid',size=7)

# consolidate neighbourhoods into larger clusters
rename= {'RCN1': ['1','8'], #tumor
         'RCN2': ['0'], # epidermis
         'RCN3': ['2','6'],
         'RCN4': ['3'], # epidermal tissue T
         'RCN5': ['5','9'], # lymphonets
         'RCN6': ['4'], # perivascular
         'RCN7': ['7'] # dermal macrophages
         }


adata = sm.hl.rename (adata, rename, from_column='spatial_kmeans', 
                                    to_column='spatial_kmeans_consolidated')

sm.pl.spatial_scatterPlot (adata, colorBy = 'spatial_kmeans_consolidated', subset = 'MIS', s=2)

customColors = { 'RCN1' : '#ced4da',
                'RCN2' : '#ffd166',
                'RCN3' : '#0d3b66',
                'RCN4' : '#c1121f',
                'RCN5' : '#219ebc',
                'RCN6' : '#f77f00',
                'RCN7' : '#000000'               
    }

# greying
customColors = { 'RCN1' : '#ced4da',
                'RCN2' : '#ffd166', #
                'RCN3' : '#0d3b66',
                'RCN4' : '#c1121f',
                'RCN5' : '#219ebc',
                'RCN6' : '#f77f00', #
                'RCN7' : '#000000'               
    }

# edede9 grey
sm.pl.spatial_scatterPlot (adata, colorBy = 'spatial_kmeans_consolidated', subset = 'MIS', s=2, customColors=customColors)
sm.pl.spatial_scatterPlot (adata, colorBy = 'spatial_kmeans_consolidated', subset = 'INV', s=2, customColors=customColors)
       
# custom scatter plot
bdata = adata[adata.obs["imageid"] == 'MIS']
bdata = adata[adata.obs["imageid"] == 'INV']
df = bdata.obs[['X_centroid','Y_centroid','phenotype','spatial_kmeans_consolidated']]

# Create a scatter plot for all points, color them grey
plt.figure(figsize=(20, 20))
s=20
# Color mapping dictionary for 'spatial_kmeans'
color_dict = {'RCN5': '#219ebc', 'RCN4': '#c1121f'}
color_dict = {'RCN1': '#bc4749', 'RCN2': '#ffd166', 'RCN6': '#3a5a40'}
# plot
plt.scatter(df['X_centroid'], df['Y_centroid'], c='#e5e5e5', s=s,alpha=0.5)
# Filter DataFrame to only include rows where 'spatial_kmeans' is 9, 5, or 3
filtered_df = df[df['spatial_kmeans_consolidated'].isin(['RCN1', 'RCN2', 'RCN6'])]
#filtered_df = df[df['spatial_kmeans_consolidated'].isin(['RCN5', 'RCN4'])]
# Overlay scatter plot for 'spatial_kmeans' 9, 5, 3 with specified colors
#for k in ['RCN5', 'RCN4']:
for k in ['RCN1', 'RCN2', 'RCN6']:
    subset = filtered_df[filtered_df['spatial_kmeans_consolidated'] == k]
    plt.scatter(subset['X_centroid'], subset['Y_centroid'], c=color_dict[k], label=f'spatial_kmeans_consolidated {k}', s=s)
# Add a solid circle around points that have 'macrophage' in the 'phenotype' column
macrophage_df = df[df['phenotype'] == 'Macrophage']
#plt.scatter(macrophage_df['X_centroid'], macrophage_df['Y_centroid'], s=s, alpha=0.65, facecolors='none', edgecolors='#6c584c', linewidth=1, label='Macrophage')
# Add labels and title
plt.legend(loc='upper left', bbox_to_anchor=(1, 1)); plt.gca().invert_yaxis(); plt.grid(False); plt.xticks([]); plt.yticks([])
plt.savefig(r'C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\figures_raw\spatial_neighbourhood_RCN126_MIS.png')


# color each RCN\
# custom scatter plot
bdata = adata[adata.obs["imageid"] == 'INV']
df = bdata.obs[['X_centroid','Y_centroid','phenotype','spatial_kmeans_consolidated']]

# Create a scatter plot for all points, color them grey
plt.figure(figsize=(20, 20))
s=20
RCN = ['RCN7']
# Color mapping dictionary for 'spatial_kmeans'
color_dict = {RCN[0]: '#252422'}
# plot
plt.scatter(df['X_centroid'], df['Y_centroid'], c='#e5e5e5', s=s,alpha=0.5)
# Filter DataFrame to only include rows where 'spatial_kmeans' is 9, 5, or 3
filtered_df = df[df['spatial_kmeans_consolidated'].isin(RCN)]
# Overlay scatter plot for 'spatial_kmeans' 9, 5, 3 with specified colors
for k in RCN:
    subset = filtered_df[filtered_df['spatial_kmeans_consolidated'] == k]
    plt.scatter(subset['X_centroid'], subset['Y_centroid'], c=color_dict[k], label=f'spatial_kmeans_consolidated {k}', s=s)
# Add labels and title
plt.legend(loc='upper left', bbox_to_anchor=(1, 1)); plt.gca().invert_yaxis(); plt.grid(False); plt.xticks([]); plt.yticks([])
#plt.savefig('C:/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/figures_raw/RCNs/' + RCN[0] + '_INV.pdf')
plt.savefig('C:/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/figures_raw/RCNs/' + RCN[0] + '_INV.png')

# Add a solid circle around points that have 'macrophage' in the 'phenotype' column
macrophage_df = df[df['phenotype'] == 'Macrophage']
plt.scatter(macrophage_df['X_centroid'], macrophage_df['Y_centroid'], s=s, alpha=0.65, facecolors='none', edgecolors='#6c584c', linewidth=1, label='Macrophage')





# 9, 5, 3 - 7

# 3 neighborhoods
# macrophages
# all other cedlls


sm.pl.spatial_scatterPlot (adata, colorBy='spatial_kmeans', subset= 'MIS', s=3)
sm.pl.spatial_scatterPlot (adata, colorBy='phenotype', subset= 'MIS', s=3)

sm.pl.voronoi(adata, color_by='spatial_kmeans', subset= 'MIS')


bdata = adata[adata.obs['phenotype'] == 'Tumor']
bdata = bdata[bdata.obs['imageid'] == 'MIS_PRAME']


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
    



##############################################################################
# Individial analysis
##############################################################################

# read h5ad
adata = ad.read(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iia-quantification3.h5ad")
adata = ad.read("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification3.h5ad")

# read invasive
adata = ad.read(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iic-quantification.h5ad")
adata = ad.read("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iic-quantification.h5ad")


# MIS
feature_table_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-quantification4.csv'
feature_table_path = r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iia-quantification3.csv"

# Invasive
feature_table_path = '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iic-quantificationV2.csv'
feature_table_path = r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iic-quantification.csv"


data = pd.read_csv(feature_table_path)
data.dropna(inplace=True)
# drop sox9 column
data.columns

# remove high expression cells
columns_to_check = ['5hmc', 'MART1', 'HLA-AB','SOX10', 'S100B', 'panCK', 'laminABC', 'S100A', 'CD31',
       'CD206', 'pMLC2', 'CD4', 'CD20', 'CD163', 'IRF1', 'CD3E', 'CD8a',
       'CD11b', 'FOXP3', 'Ki67', 'CD11c', 'SOX9', 'CD103', 'b-actin',
       'podoplanin', 'yH2AX']
data = data[(data[columns_to_check] <= 100).all(axis=1)]

sm.pl.distPlot(adata, 
             layer=None, 
             markers=columns_to_check, 
             plotGrid=True, 
             ncols=5)

reorder = ['5hmc', 'MART1', 'SOX10', 'S100B', 'panCK', 'laminABC', 'S100A', 'CD31',
       'CD206', 'pMLC2', 'CD4', 'CD20', 'CD163', 'IRF1', 'CD3E', 'CD8a',
       'CD11b', 'FOXP3', 'Ki67', 'CD11c', 'SOX9', 'CD103', 'b-actin',
       'podoplanin', 'yH2AX', 'PD1', 'HLA-AB', 'PRAME',
       'SOX9 for gating','MX1_SPOTS', 'yH2AX_SPOTS', 'yH2AX.1', 'LAG3_SPOTS', 'GZMB_SPOTS', 
       'volume', 'X_centroid', 'Y_centroid', 'Z_centroid', '1_eigen',
       '2_eigen', '3-eigen', '1_axislength', '2_axislength', '3_axislength']

data = data[reorder]

# add cell ID
data.insert(0, 'CellID', data.index) # add cellid

# bring HLA-AB to the front
cols = list(data.columns)
cols.insert(3, cols.pop(cols.index('HLA-AB'))) # Move 'HLA-AB' to the 3rd position (index 2)
data = data[cols]

# bring PD1 to the front
cols = list(data.columns)
cols.insert(18, cols.pop(cols.index('PD1'))) # Move 'HLA-AB' to the 3rd position (index 2)
data = data[cols]

# fix the SOX9 gating issue
data['SOX9'] = data.apply(lambda row: 0 if row['SOX9 for gating'] == 0 else row['SOX9'], axis=1)

data.to_csv(feature_table_path, index=False) # replace original data

# load it into scimap
adata = sm.pp.mcmicro_to_scimap (feature_table_path, log=False, split='SOX9 for gating')
adata = sm.pp.mcmicro_to_scimap (feature_table_path, log=False, split='volume')

# write h5ad
# MIS
adata.write(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\MIS_individual.h5ad")
# invasive
adata.write(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iic-quantification.h5ad")

# gate
# MIS gates
manual_gate = pd.read_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iia-gates.csv")
manual_gate = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/F8iia-gates.csv')

# invasive gates
manual_gate = pd.read_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\F8iic-gates.csv")


# rescale
adata = sm.pp.rescale (adata, gate=manual_gate, log=False)

# cell phenotyping
phenotype =  pd.read_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\phenotype_workflow_nocd15.csv")
phenotype =  pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/phenotype_workflow_nocd15.csv')

# phenotype cells
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype)
adata.obs['phenotype'].value_counts()

# rename cells
rename= {'CD4 T': ['T cells']}
adata = sm.hl.rename(adata, rename, from_column='phenotype', to_column='phenotype')

# stacked barplot
sm.pl.stacked_barplot (adata,x_axis='imageid',y_axis='phenotype')

# remove unknown cells
adata = adata[adata.obs['phenotype'] != 'Unknown']

# remove tumor
adata = adata[adata.obs['phenotype'] != 'Tumor']

# remove tumor
adata = adata[adata.obs['phenotype'] != 'keratinocytes']

# morphology
# log volume and eigne
adata.obs['eigen ratio'] = adata.obs['1_eigen'] / adata.obs['3-eigen']
adata.obs['axis length ratio'] = adata.obs['1_axislength'] / adata.obs['3_axislength']

# single cell hetamp
tumor = adata[adata.obs['phenotype'] == 'Tumor']
markers = ['SOX10', 'MART1', 'PRAME', '5hmc', 'S100B']

sc.pp.neighbors(tumor, n_neighbors=20)

sc.tl.leiden(tumor, key_added='clusters', resolution=0.5)

sc.pl.heatmap(tumor, markers, groupby='clusters', swap_axes=True, use_raw=False, cmap='vlag', dendrogram=True, vmin=0, vmax=1)


exportPhenotype = pd.DataFrame({'CellID': tumor.obs['CellID'], 'cluster': tumor.obs['clusters']}) 
exportPhenotype.to_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\phenotypes.csv",index=False)
exportPhenotype.to_csv("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/tumor_clusters.csv",index=False)

%matplotlib qt

# combat
adata = sm.pp.combat (adata,
                batch='imageid',
                layer='raw',
                log=True,
                replaceOriginal=False,
                label='combat')

bdata = ad.AnnData(adata.layers['combat'])
bdata.obs = adata.obs.copy()

sc.tl.leiden(bdata, resolution = 0.1) # Clustering the neighborhood graph
sc.pl.umap(bdata, color=['leiden'], cmap= 'vlag', use_raw=False, s=100) # Plot the UMAP
bdata = bdata[bdata.obs['leiden'] != '5']



# UMAP 
sc.pp.neighbors(bdata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(bdata,n_components=3) # Build a UMAP to visualize the neighbourhood graph
sc.pl.umap(bdata, color=['phenotype'], palette= 'tab20', use_raw=False, size=0.1, projection='3d') # Plot the UMAP

# replace bdata.X
cdata = adata[bdata.obs.index].copy()
bdata.X = cdata.X.copy()
bdata.var.index = adata.var.index.copy()
sc.pl.umap(bdata, color=['MART1','CD4', 'CD8a', 'CD3E',  'CD103','FOXP3','CD31','CD11c','CD206'], cmap= 'vlag', projection='3d', use_raw=False, s=25, vmin=0, vmax=1) # Plot the UMAP
sc.pl.umap(bdata, color=['MART1'], cmap= 'vlag', projection='3d', use_raw=False, s=25, vmin=0, vmax=1) # Plot the UMAP
sc.pl.umap(bdata, color=[ 'CD3E'], cmap= 'vlag', projection='3d', use_raw=False, s=25, vmin=0, vmax=1) # Plot the UMAP
sc.pl.umap(bdata, color=[ 'panCK'], cmap= 'vlag', projection='3d', use_raw=False, s=25, vmin=0, vmax=1) # Plot the UMAP
sc.pl.umap(bdata, color=['CD31'], cmap= 'vlag', projection='3d', use_raw=False, s=25, vmin=0, vmax=1) # Plot the UMAP




# leiden clustering
sc.tl.leiden(adata, resolution = 0.5) # Clustering the neighborhood graph
sc.pl.umap(adata, color=['SOX10','CD31','CD3E', 'CD8a','FOXP3','CD11c', 'CD103','leiden'], cmap= 'vlag', use_raw=False, s=25) # Plot the UMAP
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, fontsize=10)

# export phenotyping
exportPhenotype = pd.DataFrame({'CellID': adata.obs['CellID'], 'phenotype': adata.obs['phenotype']}) 
exportPhenotype = pd.DataFrame({'CellID': adata.obs['CellID'], 'phenotype': adata.obs['phenotype'], 'cluster': adata.obs['leiden']}) 
exportPhenotype.to_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\phenotypes.csv",index=False)
exportPhenotype.to_csv("/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/projects/2023_3D/data/phenotypes.csv",index=False)
# export invasive
exportPhenotype.to_csv(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\data\invasive_phenotypes.csv",index=False)


sc.tl.dendrogram(adata, groupby='phenotype')


# matrix plot
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag", vmin=0, vmax=1, standard_scale='var')
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='leiden', dendrogram=True, use_raw=False, cmap="vlag", vmin=0, vmax=1, standard_scale='var')


sc.pl.umap(adata, color=['SOX10','MART1','CD31','CD20','CD3E', 'CD8a','FOXP3','CD11c', 'CD206', 'CD103'], cmap= 'vlag', use_raw=False, s=10, ncols=2) # Plot the UMAP
sc.pl.umap(adata, color=['volume_log','1_eigen_log','2_eigen_log','3-eigen_log'], cmap= 'vlag', use_raw=False, s=20, ncols=3) # Plot the UMAP
sc.pl.umap(adata, color=['volume_log'], cmap= 'vlag', use_raw=False, s=20, ncols=3) # Plot the UMAP

sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag", vmin=0, vmax=1, standard_scale='obs')


# custom heatmap
d1= pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
d2 = adata.obs[['collagen distance', 'volume','1_eigen', '2_eigen', '3-eigen', 'MX1spots_log', 'granzymeBspots_log', 'LAG3spots_log',
'lysozymeCspots_log', 'catalasespots_log', 'SOX9spots_log']]
d2 = adata.obs[['volume','1_eigen', '2_eigen', '3-eigen']]
combined = pd.concat([d1,d2], axis=1)
combined['eigen1/3'] = combined['1_eigen'] / combined['3-eigen'] 
combined['eigen1/2'] = combined['1_eigen'] / combined['2_eigen']

bdata = ad.AnnData (combined)
bdata.obs = adata.obs[['X_centroid', 'Y_centroid', 'Z_centroid','CellID', 'imageid', 'phenotype']]
sc.pl.matrixplot(bdata, var_names= bdata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag", vmin=0, vmax=1, standard_scale='var')

cdata = adata.copy()
sc.tl.umap(cdata, n_components=3) # Build a UMAP to visualize the neighbourhood graph
sc.pl.umap(cdata, color=['phenotype'], cmap= 'vlag', use_raw=False, s=300, projection='3d') # Plot the UMAP


# selected cells that are something
selected_cells = ["4433","5040","4117","3198","3920","1658","4395","1867","1286","2527","1519","629","1261","910","4761","4898","10729","10329","2481","7079","2884","6755","8222","2868","1223","5636","2562","2434","1617","1012","3669","7154","1084","1961","1838","2111","6233"]

def add_string_to_list_elements(input_list, string_to_add):
    result_list = [string_to_add + element for element in input_list]
    return result_list
# usage:
string_to_append = 'F8iia-quantification3_'
modified_list = add_string_to_list_elements(selected_cells, string_to_append)

from sklearn.neighbors import BallTree
data = pd.DataFrame({'x': adata.obs['X_centroid'], 'y': adata.obs['Y_centroid'], 'z': adata.obs['Z_centroid'], 'phenotype': adata.obs['phenotype']})
kdt = BallTree(data[['x','y','z']], metric='euclidean') 
ind = kdt.query_radius(data[['x','y','z']], r=50, return_distance=False)
for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
neighbours = neighbours.dropna(how='all')

n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
n.index = n.index.get_level_values(0) # Drop the multi index
n = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
n['cell'] = n.index
n_freq = n.groupby(['cell','neighbour_phenotype']).size().unstack().fillna(0)

# subset the neighbourhoods of intrest
n_freq_subset = n_freq.loc[modified_list]
n_freq_subset_percentage = n_freq_subset.div(n_freq_subset.sum(axis=1), axis=0) * 100

# sort by CD8
sort_column = 'CD8 T'
n_freq_subset_percentage = n_freq_subset_percentage.sort_values(by=sort_column, ascending=True)
n_freq_subset_percentage = pd.concat([n_freq_subset_percentage[sort_column], n_freq_subset_percentage.drop(sort_column, axis=1)], axis=1)

# plot them
# Plot the stacked bar plot
ax = n_freq_subset_percentage.plot(kind='bar', stacked=True, figsize=(20, 10), width=1, grid=False)
ax.legend(title="Category", bbox_to_anchor=(1.04, 1), loc='upper left')
plt.tight_layout()
plt.show()

# remove tumor cells
columns_to_drop = ['Unknown', 'Tumor']
n_freq_subsete_dropped = n_freq_subset.drop(columns=columns_to_drop, inplace=False)
n_freq_subset_percentage_dropped = n_freq_subsete_dropped.div(n_freq_subsete_dropped.sum(axis=1), axis=0) * 100
n_freq_subset_percentage_dropped = n_freq_subset_percentage_dropped.sort_values(by=sort_column, ascending=True)
n_freq_subset_percentage_dropped = pd.concat([n_freq_subset_percentage_dropped[sort_column], n_freq_subset_percentage_dropped.drop(sort_column, axis=1)], axis=1)


# plot them
ax = n_freq_subset_percentage_dropped.plot(kind='bar', stacked=True, figsize=(20, 10), width=1, grid=False)
ax.legend(title="Category", bbox_to_anchor=(1.04, 1), loc='upper left')
plt.tight_layout()
plt.show()




# rename clusters
rename= {'T1': ['0'],
         'T2': ['5','7'],
         'T3': ['8','11'],
         'CD103 cells': ['1','2'],
         'c3': ['3'],
         'macrophages': ['4','13'],
         'dendritic': ['6'],
         'endothelial': ['9'],
         'IRF1': ['10'],
         'c12': ['12'],
         'treg': ['14'],
         'proliferating dendritic': ['15'],
         'B cells': ['16']
         }
adata = sm.hl.rename (idata, rename, from_column='leiden', 
                                    to_column='leiden_phenotype')

sc.pl.matrixplot(idata, var_names= adata.var.index, groupby='leiden_phenotype', dendrogram=True, use_raw=False, cmap="vlag", standard_scale='var')
sc.pl.umap(idata, color=['leiden_phenotype'], cmap= 'vlag', use_raw=False, s=25) # Plot the UMAP

# perform UMAP
sc.pp.neighbors(mdata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(mdata) # Build a UMAP to visualize the neighbourhood graph
sc.pl.umap(mdata, color=['phenotype'], cmap= 'vlag', use_raw=False, s=25) # Plot the UMAP
sc.pl.matrixplot(mdata, var_names= adata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag", standard_scale='var')










# 3D scattter plot
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'browser'


def plotlyplot (adata, phenotype='phenotype',size=2):
    df = pd.DataFrame({'x':adata.obs['X_centroid'], 'y':adata.obs['Y_centroid'], 'z':adata.obs['Z_centroid'], 'col': adata.obs[phenotype]})
    fig = px.scatter_3d(df, x='x', y='y', z='z',
                  color='col')
    fig.update_traces(marker=dict(size=size),selector=dict(mode='markers'),hoverlabel = dict(namelength = -1))
    fig.update_yaxes(autorange="reversed", tickformat='g')
    fig.update_xaxes(tickformat='g')
    #fig.update_xaxes(showgrid=False,zeroline=False)
    #fig.update_yaxes(showgrid=False,zeroline=False)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'})
    fig.write_html("/Users/aj/Downloads/3D.html")
    fig.show()

plotlyplot (adata,size=3)
plotlyplot (mdata,size=4)

def plotly (adata,phenotype='phenotype',image_id=None,x='X_centroid',y='Y_centroid',size=2, **kwargs):
    if image_id is not None:
        adata = adata[adata.obs['imageid'] == image_id]    
    data = pd.DataFrame({'x':adata.obs[x], 'y':adata.obs[y],'col': adata.obs[phenotype]})
    data = data.sort_values(by=['col'])
    fig = px.scatter(data, x="x", y="y", color="col", **kwargs)
    fig.update_traces(marker=dict(size=size),selector=dict(mode='markers'),hoverlabel = dict(namelength = -1))
    fig.update_yaxes(autorange="reversed", tickformat='g')
    fig.update_xaxes(tickformat='g')
    #fig.update_xaxes(showgrid=False,zeroline=False)
    #fig.update_yaxes(showgrid=False,zeroline=False)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'})
    return fig


plotly (adata,size=10,phenotype='phenotype_1',image_id='MIS')
plotly (mdata,size=10)



# spatial analysis
# spatial distance
adata = sm.tl.spatial_distance (adata, 
                               x_coordinate='X_centroid', y_coordinate='Y_centroid', 
                               z_coordinate='Z_centroid', 
                               phenotype='phenotype', 
                               subset=None, 
                               imageid='imageid', 
                               label='spatial_distance')

# plot spatial distance
sm.pl.spatial_distance (adata,heatmap_standard_scale=None, log=True, heatmap_cmap='bone')
sm.pl.spatial_distance (adata, method='numeric',distance_from='Tissue T', log=True)
sm.pl.spatial_distance (adata, method='numeric',distance_from='CD8 T', log=True)


# modify collogen distance
c1 = adata.uns['spatial_distance']
c2 = adata.obs['collagen distance']
cobined_distance = pd.concat ([c1,c2], axis=1)
adata.uns['spatial_distance_mod'] = cobined_distance

# plot
adata.obs['phenotype'].unique()
distance_to = ['Tumor', 'B cells', 'keratinocytes', 'Tissue T', 'Macrophage',
       'CD4 T', 'CD8 T', 'Dendritic cells', 'endothelial', 'T cells',
       'Myeloid', 'T reg']
spatial_distance (adata, spatial_distance='spatial_distance_mod', method='numeric',distance_from='collagen distance', distance_to=distance_to, log=True)





adata = sm.pp.mcmicro_to_scimap(feature_table_path= r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\ztempData\exampleSpatialTable.csv")
adata.write(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\ztempData\example.h5ad")



adata = ad.read(r"C:\Users\aj\Downloads\Z147_1_750.h5ad")

adata.obs['phenotype'].value_counts()

































