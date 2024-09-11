# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 17:23:43 2023
@author: aj
Complex heatmap function
"""


from PyComplexHeatmap import *

plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
plt.rcParams['font.family']='sans serif'
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype']=42


# drop unknown cells
adata = adata[adata.obs['phenotype'] != 'Unknown']


groupby = 'phenotype'
morphology = ['eigen ratio', 'axis length ratio']
morphology = ['eigen ratio', 'volume']

data = pd.DataFrame (adata.X, index=adata.obs.index, columns=adata.var.index)
meta = adata.obs.copy()

# summarize based on a fiven column groupby

merged_df = data.merge(meta, left_index=True, right_index=True)
# Group by the 'groupby' column in the merged DataFrame and compute the mean for each group
grouped_avg = merged_df.groupby(groupby).mean()

# expression data
groupby_data = grouped_avg[data.columns]
# morphology data
groupby_morph = grouped_avg[morphology]
scaled_groupby_morph = (groupby_morph - groupby_morph.min()) / (groupby_morph.max() - groupby_morph.min())



plt.figure(figsize=(10, 5))

row_ha = HeatmapAnnotation( EigenRatio=anno_barplot(scaled_groupby_morph['eigen ratio'],
                                                    colors = '#adb5bd',
                                                    legend=False), 
                           axis=0)

ClusterMapPlotter(data=groupby_data,show_rownames=True,show_colnames=True,
                            cmap='vlag', row_cluster=True, col_cluster=True, 
                            standard_scale=1, plot_legend=True,
                            linewidths=0.05, linecolor='black',
                            legend_gap=7,
                            vmin=0, vmax=1,
                            right_annotation=row_ha)

plt.savefig(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\projects\2023_3D\figures_raw\Invasive_heatmap.pdf",bbox_inches='tight')


scaled_groupby_morph.iloc[:,0].nunique()
len(['grey'] * len(scaled_groupby_morph))












#Annotate the rows with average > 0.3
df_rows = df_heatmap.apply(lambda x:x.name if x.sample4 > 0.5 else None,axis=1)
df_rows=df_rows.to_frame(name='Selected')
df_rows['XY']=df_rows.index.to_series().apply(lambda x:'A' if int(x.replace('Fea',''))>=15 else 'B')

row_ha = HeatmapAnnotation(
                           #Scatter=anno_scatterplot(df_heatmap.sample4.apply(lambda x:round(x,2)),height=12,cmap='jet',legend=False),
                           Bar=anno_barplot(df_heatmap.sample4.apply(lambda x:round(x,2)), height=15,cmap='rainbow',legend=False),
                           #selected=anno_label(df_rows,colors='red',relpos=(-0.05,0.4)),
                           #label_kws={'rotation':30,'horizontalalignment':'left','verticalalignment':'bottom'},
                           # axis=0,verbose=0
                           )

col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=10),
                           AB=anno_simple(df.AB,add_text=True),axis=1,
                           CD=anno_simple(df.CD,add_text=True),
                           EF=anno_simple(df.EF,add_text=True,
                                            legend_kws={'frameon':True}),
                           G=anno_boxplot(df_box, cmap='jet',legend=False),
                           verbose=0)

plt.figure(figsize=(5.5, 6.5))
cm = ClusterMapPlotter(data=df_heatmap, top_annotation=col_ha,
                       right_annotation=row_ha,
                       col_cluster=True,row_cluster=True,
                       col_split=df.AB,row_split=2,
                       col_split_gap=0.5,row_split_gap=0.8,
                       label='values',row_dendrogram=True,
                       show_rownames=False,show_colnames=True,
                       tree_kws={'row_cmap': 'Set1'},verbose=0,legend_gap=5,
                       cmap='RdYlBu_r',xticklabels_kws={'labelrotation':-90,'labelcolor':'blue'})

plt.show()



























