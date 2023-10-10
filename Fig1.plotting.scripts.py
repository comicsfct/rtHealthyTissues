# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 01:15:41 2023

@author: paulo
"""

# import more than a few packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# to format plot axis and spines
def BarplotAxisFormat(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.tick_params(direction = 'out', top=False, right = False, bottom = True)
    ax.spines['left'].set_position(('axes',-0.02)) 

# to make tissue names shorter; 
# to give a specific color to each tissue;
from utils import tissues_abbrev_colors
tissue_acronyms = tissues_abbrev_colors.tissue_acronyms
tissue_colors = tissues_abbrev_colors.tissue_colors


# import master table
master_table = pd.read_table('master.table.GTEx.samples.RT.analysis.tab')

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 1C
# plot number of expressed genes of each tissue classified according to our computational approach:
# RT genes (red); NRT genes (blue); and undefined (UND genes, gray).

# reshape data
all_counts_merged_tissue = master_table.groupby(['tissue','RT']).count()['gene_id']\
    .reset_index().pivot(columns = 'RT', index = 'tissue', values = 'gene_id')
    
all_counts_merged_tissue.columns.name = ''
all_counts_merged_tissue['ALL']  = all_counts_merged_tissue.sum(axis = 1)
all_counts_merged_tissue['RT%']  = all_counts_merged_tissue['RT']/all_counts_merged_tissue['ALL'] * 100
all_counts_merged_tissue['UND%'] = all_counts_merged_tissue['UND']/all_counts_merged_tissue['ALL'] * 100

data = all_counts_merged_tissue.reset_index()

fig, ax = plt.subplots(figsize = (4,2), dpi = 150)
#plt.title('Total expressed genes and RT genes per tissue \n', fontsize = 8)

sns.barplot(data = data, y = 'ALL', x = 'tissue', color = 'steelblue', label = 'NRT', 
            edgecolor = 'black', linewidth = .6);

sns.barplot(data = data, y = 'UND', x = 'tissue', color = 'grey', label = 'UND',
            edgecolor = 'black', linewidth = .6);

sns.barplot(data = data, y = 'RT', x = 'tissue', color = 'crimson', label = 'RT',
            edgecolor = 'black', linewidth = .6);

# format legend
plt.legend(frameon = False, fontsize = 6, ncol = 3, loc = 9, bbox_to_anchor = (0.5, 1.05))

# format axis
plt.ylabel('Expressed Genes (#)', fontsize = 7); plt.xlabel('')
ax.set_ylim([0,20000]); ax.set_xlim([-1, data.shape[0]]);

# format text labels
plt.xticks(rotation = 90, fontsize = 6, ha="center");
BarplotAxisFormat(ax); plt.tick_params(axis='y', labelsize = 6)

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_xticklabels()]).map(tissue_acronyms).values
ax.set_xticklabels(labels);

# # add % of DoGs
# dogs_perc = data['RT%'] 
# heigths = data['UND']

# for i in range(data.shape[0]):
#     #ax.text(i, 3000, str(int(dogs_perc[i])) + '\n%', color = 'darkred', ha = 'center', fontsize = 5, weight = 'bold')
#     ax.text(i +.1, heigths[i]+ 800, str(int(dogs_perc[i])) + ' %', color = 'darkred', 
#             ha = 'center', fontsize = 4, weight = 'bold', rotation = 'vertical')
    
plt.show()
 
# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 1D
# Find common genes and dogs across tissues

common_dogs_list = []
common_genes_list = []

for tissue in master_table['tissue'].unique():
    
    tissue_expgenes = master_table[(master_table['tissue'] == tissue)]['gene_id'].tolist()
    tissue_dogs = master_table[(master_table['tissue'] == tissue) & (master_table['RT']=='RT')]['gene_id'].tolist()
    
    common_genes_list.append(tissue_expgenes)                              
    common_dogs_list.append(tissue_dogs)
    
# find intersection ammong all tissues
common_dogs = set(common_dogs_list[0]).intersection(*common_dogs_list)
common_genes = set(common_genes_list[0]).intersection(*common_genes_list)

plt.show()

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 1D
# Find relative % of RT gene type

# import annotation file
annot_bed = pd.read_table('utils/human.gencode.v37.genes.bed', header = None)
annot_bed.columns = ['chrm','start','end','gene_id','gene_name','strand','gene_type']

#merge with annotation file
master_table = master_table.merge(annot_bed)

# create a table containing RT genes only
rt_genes_df = master_table[master_table['RT']=='RT']

# counts per tissue and per gene type
rt_genes_df_genetype_counts = rt_genes_df.value_counts(['tissue', 'gene_type']).unstack().fillna(0)

# clean gene types
from utils import combine_gene_types
rt_genes_df_genetype_counts = combine_gene_types.CombineGeneTypes(rt_genes_df_genetype_counts)

# normalize by total number of counts for each tissue
rt_genes_df_genetype_counts = rt_genes_df_genetype_counts.apply(lambda x: x/rt_genes_df_genetype_counts.sum(axis = 1)) * 100

# show final table
rt_genes_df_genetype_counts

# plotting
fig, ax = plt.subplots(figsize = (4, 2), dpi = 150)
#plt.title('readthrough genes type \n', fontsize = 8)

rt_genes_df_genetype_counts.plot(kind = 'bar', ax = ax, stacked = True, color = plt.cm.Set2(np.arange(10)), 
                                 width = 0.8, edgecolor = 'black', linewidth = 0.6, clip_on = False)

plt.legend(bbox_to_anchor = (1,1), frameon= False, fontsize = 6);

# format axis
plt.ylabel('RT Gene Type (%)', fontsize = 8); plt.xlabel('')
ax.set_ylim([0,100]); ax.set_xlim([-1, rt_genes_df_genetype_counts.shape[0]]);

# format text labels
plt.xticks(rotation = 90, fontsize = 6);
BarplotAxisFormat(ax); plt.tick_params(axis='y', labelsize = 6)

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_xticklabels()]).map(tissue_acronyms).values
ax.set_xticklabels(labels);

plt.show()


# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 1E
# plot heatmap showing % of RT genes in common between pairs of tissues

common_genes_pairwise = pd.DataFrame([])

for n, tissue in enumerate(master_table['tissue'].unique()):
    
    common_rt_genes_pair = []
    
    for i, tis in enumerate(common_dogs_list):

        n_com = len(set(common_dogs_list[n]).intersection(common_dogs_list[i]))
        perc_com = n_com/len(common_dogs_list[n]) * 100
        common_rt_genes_pair.append(perc_com)
        
    common_genes_pairwise = pd.concat([common_genes_pairwise, pd.DataFrame(common_rt_genes_pair)], axis =1 )
    
common_genes_pairwise.columns = master_table['tissue'].unique()
common_genes_pairwise.index = master_table['tissue'].unique()

# PLOT HEATMAP

fig, ax = plt.subplots(figsize = (5,4), dpi = 120)
plt.title('% of shared RT genes', fontsize = 8)

sns.heatmap(common_genes_pairwise, cmap='Blues', #plt.cm.Spectral_r
            linewidth = 0.1, alpha = 0.7, square = True,
            yticklabels=True, xticklabels=True, cbar = True, vmin=0, vmax=100,
            cbar_kws = {"orientation": "vertical", "shrink": 0.8, "pad": 0.04,
                        'ticks': [0,25,50,75,100]});

plt.yticks(fontsize = 5); plt.xticks(fontsize = 5);

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_xticklabels()]).map(tissue_acronyms).values
ax.set_xticklabels(labels);ax.set_yticklabels(labels)

plt.show()