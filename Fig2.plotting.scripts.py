# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 12:37:47 2023

@author: paulo
"""

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
import glob,os
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

# to format plot axis and spines
def BarplotAxisFormat(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.tick_params(direction = 'out', top=False, right = False, bottom = True)
    ax.spines['left'].set_position(('axes',-0.02)) 
    
def ChrmLexicalOrder(table, delim = 'chr'):
    
    name = table.index.name
    table[name + '_num'] = table.reset_index()[name].str.split(delim).str[1].values
    index_reorder = pd.to_numeric(table[name +'_num'], errors='coerce').sort_values().index
    table_reorder = table.loc[index_reorder].drop(name +'_num', axis=1)
    
    return table_reorder

# to make tissue names shorter; 
# to give a specific color to each tissue;
from utils import tissues_abbrev_colors
tissue_acronyms = tissues_abbrev_colors.tissue_acronyms
tissue_colors = tissues_abbrev_colors.tissue_colors


# import master table
master_table = pd.read_table('master.table.GTEx.samples.RT.analysis.tab')
master_table = master_table[(master_table['RT'] == 'RT')|(master_table['RT'] == 'NRT')]

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2A
# plot expression levels of genebody vs. RT tail

tissue_to_plot = 'muscle-skeletal'

print('Correlation Plot for {}'.format(tissue_to_plot))

fig, ax = plt.subplots(figsize = (3,2), dpi = 150)
plt.title('RT expression vs gene expression \n', fontsize = 8)

data = master_table[master_table['tissue'] == tissue_to_plot]
data = data[data['dogFPKM'] > 0]
data['log2dogFPKM'] = np.log2(data['dogFPKM'])
data['log2geneFPKM'] = np.log2(data['geneFPKM'])

sns.scatterplot(data = data, y = 'log2geneFPKM', x = 'log2dogFPKM', color = '#FF9159',
                s = 8, alpha = 0.6, edgecolor = 'black');

#legend
#ax.legend(fontsize = 4, ncol = 3, labelspacing = 1, frameon = False, markerscale = 0.5)

BarplotAxisFormat(ax)
#ax.set_ylim([0, 12]); ax.set_xlim([0, 15])
plt.yticks(fontsize = 6); plt.xticks(fontsize = 6);
plt.ylabel('RT expression (FPKM)', fontsize = 7); 
plt.xlabel('gene expression (FPKM)', fontsize = 7);

corr, pval = spearmanr(data['log2dogFPKM'], data['log2geneFPKM'])

#ax.text(1,1.3, 'corr = {}'.format(round(corr,3)), fontsize = 7);
print('corr = {:.3}'.format(corr))
print('pval = {:.3}'.format(pval))

# save source data for the reviewers
source_data_folder = 'source_data/'
data.to_csv(source_data_folder + 'source.data.Fig2a.csv', index = False)

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2B
# Compute Correlations for all tissues

data = master_table.copy()
data = data[data['dogFPKM'] > 0]
data['log2dogFPKM'] = np.log2(data['dogFPKM'])
data['log2geneFPKM'] = np.log2(data['geneFPKM'])

corr_plot = []

for tissue in data['tissue'].unique():
    
    tissue_data = data[data['tissue'] == tissue]
    
    corr, pval = spearmanr(tissue_data['log2dogFPKM'], tissue_data['log2geneFPKM'])

    corr_plot.append([tissue, corr, pval])
    
corr_plot = pd.DataFrame(corr_plot, columns=['tissue','corr','pval'])
#corr_plot = corr_plot.set_index('tissue')
corr_plot = corr_plot.iloc[::-1]

fig, ax = plt.subplots(figsize = (2 ,3), dpi = 150)

# prepare data
corr_values = corr_plot['corr'].values
pvals = corr_plot['pval'].values
tissues = corr_plot['tissue'].values

# plot lines and dots (scatter)
ax.scatter(y = tissues, x = corr_values, marker = 'o', s = 18, edgecolor = 'black', lw = 1, zorder=4, c = 'salmon')
ax.hlines(y = tissues, xmin = 0, xmax = corr_values, ls = '--', lw = 1, color = 'black', alpha = 0.5, zorder=1)

plt.yticks(fontsize = 6); plt.xticks(fontsize = 7);
ax.set_xlabel('correlation (spearman)', fontsize = 7)
ax.tick_params(axis = 'y', width = 1, length = 2)
ax.set_xlim([0, 1])

# remove spines 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_yticklabels()]).map(tissue_acronyms).values
ax.set_yticklabels(labels);

plt.show()

# save source data for the reviewers
source_data_folder = 'source_data/'
corr_plot.to_csv(source_data_folder + 'source.data.Fig2b.csv', index = False)

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2C
# Plot % of High Expressed Genes RT vs. NRT

rt_high_expgenes_df = master_table[master_table['exp_group'] == 'high'
                                  ].groupby(['tissue','RT']).count()['gene_id'].reset_index()


data_stacked = rt_high_expgenes_df.pivot(columns = 'RT', values = 'gene_id', index = 'tissue')
data_stacked = data_stacked.apply(lambda x: x/data_stacked.sum(axis = 1) * 100)

fig, ax = plt.subplots(figsize = (4,2), dpi = 150)
#plt.title('% of high expgenes with TRT \n', fontsize = 8)

colors = ['steelblue','crimson']

data_stacked.plot(kind = 'bar', ax = ax, stacked = True, color = colors,
                  width = 0.8, edgecolor = 'black', linewidth = 0.8, clip_on = False)


# format legend
plt.legend(ncol = 2, loc = 10, bbox_to_anchor = (0.5, 1.1), frameon = False, fontsize = 7)

plt.ylabel('high expressed genes (%)', fontsize = 7); plt.xlabel('')
ax.set_ylim([0,100]);

# format text labels
plt.xticks(rotation = 90, fontsize = 6, ha="center");
BarplotAxisFormat(ax); plt.tick_params(axis = 'y', labelsize = 7)

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_xticklabels()]).map(tissue_acronyms).values
ax.set_xticklabels(labels);

# save source data for the reviewers
source_data_folder = 'source_data/'
data_stacked.reset_index().to_csv(source_data_folder + 'source.data.Fig2c.csv', index = False)

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2D
# Plot % of High Expressed Genes RT vs. NRT

master_table_with_length = pd.DataFrame([])

for tissue in master_table['tissue'].unique():
    
    print(tissue, ' '*30, end = '\r')
    #import dog info from each samples for a given tissue
    tissue_folder_dir = 'artdeco-out/' + tissue
    dog_files = glob.glob(tissue_folder_dir + '/dogs-filtered/*out.dogs.fpkm.txt')
    
    # get list of RT genes for a given tissues
    RT_genes_tissue = master_table[(master_table['tissue'] == tissue) & (master_table['RT'] == 'RT')]
    
    all_samples_dogs = pd.DataFrame([])
    
    for sample in dog_files:
        sample_dogs = pd.read_table(sample)
        sample_dogs.columns = ['gene_id','dogLength','-']
        all_samples_dogs = pd.concat([all_samples_dogs, sample_dogs])
    
    # same dog length varies in diff samples; here we consider the longest version
    all_samples_dogs = all_samples_dogs.groupby('gene_id').max().reset_index()
    all_samples_dogs = all_samples_dogs.merge(RT_genes_tissue).drop('-', axis =1)
    
    master_table_with_length = pd.concat([master_table_with_length, all_samples_dogs])

fig, ax = plt.subplots(figsize = (2,3), dpi = 150)
plt.title('readthrough length distribution', fontsize = 7)

sns.boxplot(data = master_table_with_length, y = 'tissue', x = 'dogLength', ax=ax, palette = tissue_colors, orient = 'h',
             width = 0.7, linewidth = 0.8, fliersize = 0.1)

# sns.stripplot(data = data, y = 'tissue', x = 'dogLength', ax = ax,  palette = tissue_colors, orient='h',
#               size = 0.8, jitter = True, dodge = True, alpha = 0.8)

# format text labels
plt.xticks(fontsize = 6); plt.yticks(fontsize = 4);
plt.tick_params(axis = 'y', labelsize = 6)
plt.tick_params(axis = 'x', labelsize = 6)

# replace tissue labels with shorter names
labels = pd.Series([label.get_text() for label in ax.get_yticklabels()]).map(tissue_acronyms).values
ax.set_yticklabels(labels);

# remove spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_position(('axes',-0.02))

# format logarithmic scale
plt.xscale('symlog'); ax.set_xlim([1000, 100000]);
ax.set_xlabel('Log Length (bp)', fontsize = 7);

# save source data for the reviewers
source_data_folder = 'source_data/'
master_table_with_length.to_csv(source_data_folder + 'source.data.Fig2d.csv', index = False)

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2E
# Heatmap representing the enrichment (red) or depletion (blue) of RT gene proportions across different chromosomes for each tissue.

# import all gene coordiantes
all_genes_bed = pd.read_table('utils/human.gencode.v37.genes.bed', header=None)
all_genes_bed.columns = ['chrm','start','end','gene_id','gene_name','strand','gene_type']

# combine with master table
master_table_plus = master_table.merge(all_genes_bed)

# extract RT and NRT entries from the master table
all_tissues_RT = master_table_plus[(master_table_plus['RT'] == 'RT')]
all_tissues_NRT = master_table_plus[(master_table_plus['RT'] == 'NRT')]

# reshape each dataframe to count the number of genes in each chromosome for each tissue
chrm_RT = all_tissues_RT.groupby(['chrm','tissue']).count()['gene_id'].reset_index()
chrm_RT_mod = chrm_RT.pivot(index = 'chrm', columns = 'tissue', values = 'gene_id')
chrm_RT_mod = ChrmLexicalOrder(chrm_RT_mod)

chrm_NRT = all_tissues_NRT.groupby(['chrm','tissue']).count()['gene_id'].reset_index()
chrm_NRT_mod = chrm_NRT.pivot(index = 'chrm', columns = 'tissue', values = 'gene_id').sort_index()
chrm_NRT_mod = ChrmLexicalOrder(chrm_NRT_mod)

# Apply Fisher Test to all combinations of Tissues and Chromosomes
all_fisher_tests = []
all_fisher_tests_oddsratio = []
all_cont_tables = []

for tissue in chrm_RT_mod.columns:
    print(tissue + ' computing fisher for all chromosomes                       ', end = '\r')
    
    for chrm in chrm_RT_mod.index:
        
        # cuild a contingency table
        
        cont_table = pd.DataFrame([[chrm_RT_mod.loc[chrm, tissue], # dogs in a given chromosome
                                    chrm_RT_mod.loc[:, tissue].drop(chrm).sum()], # dogs in all chromosomes except the given chromosome
                                         
                                     [chrm_NRT_mod.loc[chrm, tissue], # expressed genes in given chromosome except dogs
                                      chrm_NRT_mod.loc[:,tissue].drop(chrm).sum()]]) # expressed genes in the remaining chromosomes execept without dogs

        cont_table.columns = ['genes in '+ chrm, 'genes not in '+ chrm]
        cont_table.index = ['TRT_genes','no_TRT_genes']

        # apply fischer test
        oddsratio, pvalue = fisher_exact(cont_table, alternative = 'two-sided')
        
        # append to list
        all_fisher_tests.append([tissue, chrm, pvalue])
        all_fisher_tests_oddsratio.append([tissue, chrm, oddsratio])
        
        # save contigency tables for inspection
        all_cont_tables.append([tissue, chrm, cont_table])

# format dataframes 
all_fisher_tests = pd.DataFrame(all_fisher_tests, columns = ['tissues', 'chrm', 'pvals'])
all_fisher_tests = all_fisher_tests.pivot(index = 'chrm', columns = 'tissues', values = 'pvals')
all_fisher_tests.rename_axis(None, axis=1, inplace=True)

all_fisher_tests_oddsratio = pd.DataFrame(all_fisher_tests_oddsratio, columns = ['tissues', 'chrm', 'oddsratio'])
all_fisher_tests_oddsratio = all_fisher_tests_oddsratio.pivot(index = 'chrm', columns = 'tissues', values = 'oddsratio')
all_fisher_tests_oddsratio.rename_axis(None, axis=1, inplace=True)

all_fisher_tests = ChrmLexicalOrder(all_fisher_tests)

# filter only significative and plot 
all_fisher_tests_significant = all_fisher_tests[all_fisher_tests < 0.1]
all_fisher_tests_adjusted = all_fisher_tests.apply(lambda x: fdrcorrection(x)[1])

# only significant values
all_fisher_tests_significant_adjusted = all_fisher_tests_adjusted[all_fisher_tests_adjusted < 0.05]

fig, ax = plt.subplots(figsize = (5, 4), dpi = 120)
#ax.set_title("chromosome enrichment (fisher's test ; adjusted pvals) \n", fontsize = 8)

# plot enriched - significative p-value where odds ratio is larger than 1
enriched = all_fisher_tests_significant_adjusted[all_fisher_tests_oddsratio > 1]

sns.heatmap(data = enriched, ax = ax, 
            cmap = plt.cm.Reds_r, linewidth = 0.0, alpha = 0.7, square = True,
            yticklabels=True, xticklabels=True, cbar = True, vmin=0, vmax=0.05,
            cbar_kws = {"orientation": "vertical", "shrink": 0.6, "pad": -0.11,
                       'ticks': []});

# plot decreased - significative p-value where odds ratio is smaller than 1
decreased = all_fisher_tests_significant_adjusted[all_fisher_tests_oddsratio < 1]

sns.heatmap(data = decreased, ax = ax, 
            cmap = plt.cm.Blues_r, linewidth = 0.1, alpha = 0.7, square = True,
            yticklabels=True, xticklabels=True, cbar = True, vmin=0, vmax=0.05,
            cbar_kws = {"orientation": "vertical", "shrink": 0.6, "pad": 0.03,
                       'ticks': []});

# ax.figure.axes[0].set_ylabel('enrichment', color='black', size=8, fontweight = 'normal')
# ax.figure.axes[].set_ylabel('sent', color='black', size=8, fontweight = 'normal')

# format cbar text manually
ax.text(enriched.shape[1] + 0.5, 2,'adj.pval', fontsize = 8)
ax.text(enriched.shape[1] + 0.9, 4,'0.05', fontsize = 8)
ax.text(enriched.shape[1] + 0.9, 20.4,'0.00', fontsize = 8)

# format axis
ax.set_xlabel(''); ax.set_ylabel('');
plt.tick_params(axis = 'y', labelsize = 7); plt.xticks(fontsize = 8);
ax.spines['left'].set_position(('axes',-0.03))
#ax.spines['bottom'].set_position(('axes',-0.03))


# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 2F
# Scatterplot RT transcripts percentage in all tissues and density of expressed genes for each chromosome

# compute number of expgenes per chrmosome in each tissue
chrm_expgenes = master_table_plus.groupby(['chrm','tissue']).count()['gene_id'].reset_index()
chrm_expgenes_mod = chrm_expgenes.pivot(index = 'chrm', columns = 'tissue', values = 'gene_id')
chrm_expgenes_mod = ChrmLexicalOrder(chrm_expgenes_mod)

# import info regarding chm sizes
genes_per_chrm = all_genes_bed.groupby('chrm').count()['gene_id'].reset_index()
chrm_len = pd.read_table('utils/human.gencode.v37.chrom.sizes', header = None, names = ['chrm','len'])

# compute chrm density (number of genes divided by chrm lenght in megabases)
chrm_density = genes_per_chrm.merge(chrm_len)
chrm_density['dens'] = chrm_density['gene_id']/(chrm_density['len']* 1e-6) # lenght in megabases
chrm_density = chrm_density.sort_values('dens', ascending = False).reset_index(drop = True)

# compute % of RT genes IN RELATION TO ALL EXPGENES per chrm and average across tissues
rt_percent_chrm = chrm_RT_mod/chrm_expgenes_mod * 100

rt_percent_chrm_mean = rt_percent_chrm.mean(axis = 1).reset_index()

rt_percent_chrm_mean.columns = ['chrm','RT%']
rt_percent_chrm_mean['std'] = rt_percent_chrm.std(axis = 1).reset_index()[0]

# combine data into a single dataframe
chrm_density = rt_percent_chrm_mean.merge(chrm_density).sort_values('dens')

fig, ax = plt.subplots(figsize = (4,3), dpi = 120)
#plt.title('Chrm gene density vs. RT%', fontsize = 8)

sns.scatterplot(data = chrm_density,  y = 'RT%', x = 'dens', color = 'darkcyan',
                s = 40, marker = 'o', alpha = 0.7, edgecolor = 'black');

BarplotAxisFormat(ax)
plt.yticks(fontsize = 8); plt.xticks(fontsize = 8);
plt.ylim([0,30]); plt.xlim([0,60])
plt.ylabel('RT genes (%)', fontsize = 8); plt.xlabel('Gene Density (#genes/chrm length in Mb)', fontsize = 8);

# add some labels
for i, chrm in chrm_density.iterrows():
    if chrm['chrm'] in ['chr16','chr17','chr19','chr4','chrY', 'chr13',]:
        ax.text(chrm['dens']-2.5, chrm['RT%']+2, chrm['chrm'], fontsize = 7, alpha = 0.8)
    
# statistics
corr, pval = spearmanr(chrm_density['RT%'], chrm_density['dens'])