# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:02:20 2023

@author: paulo
"""

# import more than a few packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob,os
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

def BarplotAxisFormat(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.tick_params(direction = 'out', top=False, right = False, bottom = True)
    ax.spines['left'].set_position(('axes',-0.02)) 

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 3A
# Chromatin State enrichment

all_heatmaps = pd.DataFrame([])

states_names = pd.read_table('epigenetics/epigenomic_states_info.txt')


for file in glob.glob('epigenetics/OverlapEnrichment/*txt'):
    
    print(file)
    # import overlap enrichment results
    chrmHMM_result = pd.read_table(file)[:-1]

    # re-shape names
    index = states_names.astype(str)[["MNEMONIC", "STATE NO."]].apply(" ".join, axis = 1)
    chrmHMM_result.index = index
    chrmHMM_result.drop('State (Emission order)', axis = 1, inplace = True)
    chrmHMM_result.drop('Genome %', axis = 1, inplace = True)
        
    # merge with main dataframe
    all_heatmaps = pd.concat([all_heatmaps, chrmHMM_result], axis = 1)

# normalize per column    
all_heatmaps = all_heatmaps.apply(lambda x: x - np.mean(x)/ np.std(x), axis = 0);

# re-organize column order
all_heatmaps = pd.concat([all_heatmaps.filter(like = '_NRT'), 
                          all_heatmaps.filter(like = '_RT'), all_heatmaps.filter(like = 'flank')], axis = 1)

# plot for each region

for i in ['TTSminus', 'TTSplus']:
    
    data = all_heatmaps.filter(like = i)
    
    fig, ax = plt.subplots(figsize = (10, 5), dpi = 150)
    plt.title('Overlap Enrichment -> {}'.format(i), fontsize = 8)
    
    sns.heatmap(data, cmap='BuPu', ax = ax, square = True, cbar=True, vmin = 0,
                linecolor='w', linewidth = 0.01, cbar_kws={'shrink':0.8})
    
    plt.ylabel(''); plt.yticks(fontsize = 6); plt.xticks(fontsize = 6);
    
    # new labels
    labels = data.columns.str.split('_SRX').str[0] + '_' + data.columns.str.split('_').str[-2]
    ax.set_xticklabels(labels);
    
   
# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 3B
# Chromatin State enrichment   

# create expression-matched sub samples

def MatchingSubSample(df_A, df_B, conf):
    '''creates a 'CONF'-matched subsample for df_B with size = len(df_A)'''

    df_A = df_A.sort_values(conf).reset_index(drop=True)
    df_B = df_B.sort_values(conf).reset_index(drop=True)

    # Compute the absolute differences between each element in df_A and df_B
    differences = np.abs(df_A[conf].values[:, None] - df_B[conf].values)

    # Find the index of the minimum difference for each element in df_A
    indices = np.argmin(differences, axis=1)

    # Create a set to store used indices
    used_indices = set()

    # Use the indices to extract the corresponding rows from df_B
    subsampleB = []
    for idx in indices:
        while idx in used_indices:
            
            # If the index is already used or out of bounds, find the next available index
            idx -= 1
            if idx < 0:
                # If all indices are out of bounds, break the loop
                break

        # Check if the index is within the valid range
        if 0 <= idx < len(df_B):
            
            # Add the row to the subsampleB list
            subsampleB.append(df_B.iloc[idx])
            used_indices.add(idx)

    # Convert the subsampleB list to a DataFrame
    subsampleB = pd.DataFrame(subsampleB)

    return subsampleB.reset_index(drop=True)
    
# select folder with output from bedtools intersect (= number of peaks found in each region)
# select tissues of interest (strings that exist within filename in the folder)

chromatin_counts_folder = 'epigenetics/chromatin_marks_counts_imputed/'
tissues = ['psoas', 'colon_sigmoid', 'esophagus', 'lung', 'thymus']

# select folder containing expressed genes per tissue (to create expression-matched subsamples of genes)
expgenes_folder = 'epigenetics/expgenes_beds/'

# create two expression-matched subsamples of genes for RT and NRT with same size
size = 250    # size of each subsample
perm = 1000     # number of permutations per tissue

perm_chromatin_marks_count = []

for tissue in tissues: # for each tissue in the list of interest
    
    for file in os.listdir(expgenes_folder): # find the file containing all expgenes

        if tissue in file:
                           
            tissue_expgenes = pd.read_table(expgenes_folder + file)
            
            # create one dataframe for each group of genes
            expgenes_RT = tissue_expgenes[tissue_expgenes['RT']=='RT']
            expgenes_NRT = tissue_expgenes[tissue_expgenes['RT']=='NRT']
            
            # filter universe of RT genes for DoGs > 4kbp
            dog_len_info = pd.read_table('epigenetics/bed_files_dog_coord/{}_dog_coord.tab'.format(tissue))
            big_dogs = dog_len_info[dog_len_info['length'] > 2500]['gene_id']
            
            expgenes_RT = expgenes_RT[expgenes_RT['gene_id'].isin(big_dogs)]
            
            # create subsamples of each dataframe with same size and same distribution of explevels

            for i in range(0, perm):
                
                print(tissue, 'perm {}/{}'.format(i+1,perm), ' '*20, end = '\r')
                
                expgenes_RT_subsample = expgenes_RT.sample(size)
                expgenes_NRT_subsample = MatchingSubSample(expgenes_RT_subsample, expgenes_NRT, 'fpkm')
                
                # function is written below
#                 SavePermutationSampleDistribution()
                # for each subsample of genes perform the analysis
                
                ### NEW BLOCK ###
                # for each pair of expression matched samples count how many chrm marks where found per region
                # code builds a table for each each file in the folder, separates by region, rt status, and ref epigenome
        
                chromatin_marks_count = []

                for file in os.listdir(chromatin_counts_folder):
                    
                    if tissue == 'psoas': tissue = 'muscle'
                    
                    if tissue in file:

                        epigenome = file.split('-')[0]
                        chrom_mark = file.replace('-','_').split('_')[1] 
                        region = file.split('_')[-1].split('.')[0]
                        rt = file.split('_')[-2]; 
                        
                        # count number of genes with a given chrom mark
                        
                        try:
                                
                            mark_output = pd.read_table(chromatin_counts_folder + file, header = None)

                            if rt == 'RT': # filter for RT genes in subsample
                                mark_output = mark_output[mark_output[3].isin(expgenes_RT_subsample['gene_id'])]

                            elif rt == 'NRT': # filter for NRT genes in subsample
                                mark_output = mark_output[mark_output[3].isin(expgenes_NRT_subsample['gene_id'])]

                            else: # filter for genes in RT subsample
                                mark_output = mark_output[mark_output[3].isin(expgenes_RT_subsample['gene_id'])]

                            # count number of hits for each region
                            mark_count = mark_output.groupby(3).count().shape[0]

                        except:
                            mark_count = 0
                        
                        
                        chromatin_marks_count.append([tissue, epigenome, chrom_mark,region, rt, mark_count])

                chromatin_marks_count = pd.DataFrame(chromatin_marks_count, columns=['tissue','epi','mark','region','type','gene_hits'])

                perm_chromatin_marks_count.append(chromatin_marks_count)

all_perm_combined = pd.DataFrame([])

# for tissue in tissues:

for tissue in tissues:
    
    if tissue == 'psoas': tissue = 'muscle'    
    perm_gene_hits = pd.DataFrame([]) # table to save all gene_hits from each permutation
    tissue_table_comb = pd.DataFrame([]) # table to save all genes_hits combined into a single df

    for table in perm_chromatin_marks_count:

        # find all table for a given tissue
        if table['tissue'].str.contains(tissue).iloc[0] == True:

            # combine all gene_hits arrays into a single one
            perm_gene_hits = pd.concat([perm_gene_hits, table['gene_hits']], axis = 1, ignore_index= True)

            # combine all columns into one
            perm_gene_hits_comb = perm_gene_hits.agg(np.array, axis=1)
            perm_gene_hits_mean = perm_gene_hits.mean(numeric_only=True, axis = 1).round(3)
            perm_gene_hits_std = perm_gene_hits.std(numeric_only=True, axis = 1).round(3)
            
            # merge with main table for the given tissue
            tissue_table_comb = pd.concat([table.drop('gene_hits', axis = 1), perm_gene_hits_comb], axis = 1)
            tissue_table_comb.columns = ['tissue','epi','mark','region','type','gene_hits_perm']

            # merge with main table for the given tissue
            tissue_table_comb = pd.concat([table.drop('gene_hits', axis = 1), 
                                           perm_gene_hits_comb, perm_gene_hits_mean, perm_gene_hits_std], axis = 1)

            tissue_table_comb.columns = ['tissue','epi','mark','region','type','gene_hits_perm', 'gene_hits_mean', 'gene_hits_std']

    all_perm_combined = pd.concat([all_perm_combined, tissue_table_comb], ignore_index=True)
    
# FILTER chromatin marks of interest: H3K4me3 -> H3K27ac -> H3K36me3
MOI = ['H3K4me3', 'H3K79me2', 'H3K36me3', 'H3K4me1', 'H3K27ac', 'H3K4me2', 'H3K9ac','DNase']
all_perm_combined = all_perm_combined[all_perm_combined['mark'].isin(MOI)]

# Show table
all_perm_combined['tissue'].value_counts()        

# barplots for one tissue
tissue = 'muscle'
region = 'TTSplus'

main_table = all_perm_combined

# save source data for the reviewers
source_data_folder = 'source_data/'
main_table.to_csv(source_data_folder + 'source.data.Fig3b-c.csv', index = False)

# filter data for plotting
data = main_table[(main_table['region'] == region) & (main_table['tissue']==tissue)]
# plt.hist(data['gene_hits_perm'].iloc[0])

fig, ax = plt.subplots(figsize = (4, 2), dpi = 150)
plt.title('chromatin marks hits {} {}\n'.format(region, tissue), fontsize = 7)

# order = ['H3K4me1','H3K27ac','H3K36me3', 'H3K9ac', 'H3K4me2', 'DNase']
data = data.sort_values(['type', 'mark'])

sns.boxplot(data = data.explode('gene_hits_perm').reset_index(drop = True), 
            x = 'mark', y = 'gene_hits_perm',  hue = 'type', dodge = True,
            palette = ['steelblue', 'crimson'], order = MOI,  linewidth = 1, fliersize = 0)

#legend
ax.legend(frameon = False, fontsize = 7)

# format axis
plt.ylabel('chromatin marks hits (#)', fontsize = 7); plt.xlabel('')
ax.set_ylim([0, 180]); #ax.set_xlim([-1, data.groupby('tissue').count().shape[0]]);

# format text labels
plt.xticks(rotation = 90, fontsize = 6);
BarplotAxisFormat(ax); plt.tick_params(axis='y', labelsize = 6)

