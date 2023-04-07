# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 23:53:49 2023
@author: paulo
"""

""" 
This script takes as input a folder containing one or multiple ARTDeco output folders (one for each tissue in this scenario) 
and uses the two subfolders ("quantification" and "dogs") to combine expression levels from genes and readthrough regions
into a simple dataframe for all the downstream analysis. 

The folder called "dogs-filtered" is the result of our script to filter non-stranded data
and contains all the same files in the "dogs" folder, but exluding all the genes that overlapped with genes 
being expressed in the opposite strand (bedtools intersect was used in this step)
"""

import pandas as pd
import numpy as np
import glob, os
import sys

artdeco_dir = sys.argv[1] # artdeco output folder
#artdeco_dir = 'artdeco-out-example/'

def main():

    # read input folder
    artdeco_dir = sys.argv[1] # artdeco output folder
    #artdeco_dir = 'artdeco-out-example/'
    
    # create list of lists containing all the info
    all_expression_data = GetAllDataInfo(artdeco_dir, fpkm_thres = 1)
    
    # create master table containg info for all downstream analysis
    master_table = create_master_table(all_expression_data, artdeco_dir)
    
    # save table
    master_table.to_csv('master.table.GTEx.samples.RT.analysis.tab', sep = '\t', index = False)
    print("all info compacted into a single dataframe: master.table.GTEx.samples.RT.analysis.tab")
    return master_table
    
def ExpressedGenes(tissue_folder, artdeco_out_dir, fpkm = 1):
    '''returns a list of lists containing sample name [0], gene IDs [1] 
    and corresponding gene fpkm [2], for each sample in the quantification folder.
    
    uses info from the fpkm transcripts table and max_isoform table to find the main isoform and their fpkm'''
    
    tissue_folder_dir = artdeco_out_dir + tissue_folder
    transcripts_fpkm = pd.read_table(tissue_folder_dir + '/quantification/gene.exp.fpkm.txt')
    max_isoforms = pd.read_table(tissue_folder_dir + '/quantification/max_isoform.txt')
    
    gene_isoforms = transcripts_fpkm[transcripts_fpkm['ID'].isin(max_isoforms['Transcript ID'])]
    gene_isoforms.columns = ['Transcript ID', *gene_isoforms.columns[1:]]
    #gene_isoforms.insert(0, 'Gene ID', isoforms['Gene ID'])

    gene_isoforms = gene_isoforms.merge(max_isoforms, on = 'Transcript ID', how = 'inner')
    gene_isoforms = gene_isoforms[['Gene ID','Transcript ID', *gene_isoforms.columns[1:-1]]]

    # builds a list for each sample with gene ids and gene fpkms
    # final result is a list of lists

    gene_ids = []

    for sample in gene_isoforms.iloc[:,3:]:
        gene_fpkm = gene_isoforms[sample][gene_isoforms[sample] >= fpkm]
        gene_id = gene_isoforms.iloc[gene_fpkm.index]['Gene ID']
        gene_ids.append([sample.rsplit('_Aligned', 1)[0], gene_id, gene_fpkm])
        
    # create a list of all expressed genes across samples
    # remove genes that show up in less than 25% of the samples
    
    all_gene_ids = [i for j in [sample[1] for sample in gene_ids] for i in j]
    all_gene_ids = list(np.unique(RemoveElements(all_gene_ids, int(len(gene_ids)/4))))
    
    return gene_ids, all_gene_ids

def GetGenesWithDogs(tissue_folder, artdeco_out_dir, expressed_genes_list = None, fpkm = 0.15):
    
    '''returns a list of lists containig [0] sample name, [1] gene IDs containing dogs
    and [2] the number of dogs per sample. expressed_genes_list contains the Gene ID for all expressed genes'''
    
    tissue_folder_dir = artdeco_out_dir + tissue_folder
    dog_samples = glob.glob(tissue_folder_dir + '/dogs-filtered/*out.dogs.fpkm.txt')
    
    dog_stats = [] # will save all stats for each sample

    for file in dog_samples:

        # read each dog file; get file name and clean
        dog_file = pd.read_table(file)
        sample_name = dog_file.iloc[:,-1].name.rsplit('_', 1)[0]
                
        # use only genes inside expressed_genes_list
        if type(expressed_genes_list) == list:
            dog_file = dog_file[dog_file['ID'].isin(expressed_genes_list)]
        
        #dog_number = dog_file.shape[0]
        dog_genes = dog_file['ID']
        dog_fpkm = dog_file[dog_file.iloc[:,-1] >= fpkm].iloc[:,-1]        
        dog_stats.append([sample_name, dog_genes, dog_fpkm])
    
    # list with all dogs expressed across samples
    
    all_dogs_ids = list(np.unique([i for j in [sample[1] for sample in dog_stats] for i in j]))
    
    return dog_stats, all_dogs_ids
   
def RemoveElements(lst, k):
    '''construct a dictionary mapping value to counts 
    and then use a list comprehension to filter for counts larger than a specified value
    lst: list of elements with repetead elements
    k: elements with a count number inferior to this value will be discarded'''
    
    from collections import Counter
    counted = Counter(lst)
    return [elem for elem in lst if counted[elem] >= k]
        
def GetAllDataInfo(artdeco_dir, fpkm_thres = 0):
    '''uses ExpressedGenes and GetGenesWithDogs functions to import all into into a
    single list of lists. 
    output: a list of lists for each tissue. each sublist contains:
    tissue name, expgenes_stats, expgenes_list, dogs_stats, and dogs_ids
    '''
    
    all_tissues_exp_genes = []

    for tissue in os.listdir(artdeco_dir):
        print(tissue + ' processing info                                              ', end = '\r')

        expgenes_fpkm, expgenes_list_ids = ExpressedGenes(tissue, artdeco_dir, fpkm = fpkm_thres)
        dogs_stats, dogs_list_ids = GetGenesWithDogs(tissue, artdeco_dir, expgenes_list_ids)

        # append to a geral list
        tissue_name = tissue
        all_tissues_exp_genes.append([tissue_name, 
                                      expgenes_fpkm, expgenes_list_ids,
                                      dogs_stats, dogs_list_ids])
    
    return all_tissues_exp_genes

def clean_zeros(table, col):

    new_col = []

    for tissue in table[col]:
    
        new_col.append([elem for elem in tissue if elem !=0])
    
    table[col] = pd.Series(new_col)
    
    return table

def create_master_table(all_expression_data, artdeco_dir):
    ''' 
    takes the info in the list of lists all_expression_data

    Parameters
    ----------
    all_expression_data : list of lists

    Returns
    -------
    dataframe containing all expressed gene IDs, their expression and
    readthough levels (if they exist) and the exp group (high, mid, low)

    '''
        
    all_expression_data_df = pd.DataFrame([])
    
    for tissue_info in all_expression_data[:]:
    
        # print(tissue_info[0].replace('\n','-'), '                     ', end = '\r')
    
        # extract gene expression information per sample
        sample_names = [samples[0] for samples in tissue_info[1]]
        sample_fpkms = [pd.Series(samples[2].values, index = samples[1].values) for samples in tissue_info[1]]
        sample_fpkms_df = pd.concat(sample_fpkms, axis = 1) 
        sample_fpkms_df.columns = sample_names
        
        # compute median across all samples
        sample_fpkms_median = sample_fpkms_df.median(axis = 1).reset_index()
        sample_fpkms_median.columns = ['gene_id','geneFPKM']
        
        # sort and reset the index
        sample_fpkms_median = sample_fpkms_median.sort_values('geneFPKM', ascending = False).dropna()
        
        # filter for only expressed genes in at least 25% of all samples
        sample_fpkms_median = sample_fpkms_median[sample_fpkms_median['gene_id'].isin(tissue_info[2])]
        
        # add expression group information
        # spli data into 3 groups of expression: low, medium and high: attribute a label to each gene accordingly
        high, med, low = np.array_split(sample_fpkms_median, 3)
    
        high.insert(2, 'exp_group', 'high')
        med.insert(2, 'exp_group', 'med')
        low.insert(2, 'exp_group', 'low')
    
        gene_exp_groups = pd.concat([high, med, low])
        
        # extract readthrough expression information (already filtered to only expressed genes in 25% of the samples)

        sample_names = [samples[0] for samples in tissue_info[3]]
        dogs_fpkms = [pd.Series(samples[2].values, index = samples[1].values) for samples in tissue_info[3]]
        dogs_fpkms_df = pd.concat(dogs_fpkms, axis = 1) 
        dogs_fpkms_df.columns = sample_names
     
        # compute median across all samples
        dogs_fpkms_median = dogs_fpkms_df.median(axis = 1).reset_index()
        dogs_fpkms_median.columns = ['gene_id','dogFPKM']
         
        # merge dogs with expression data
        dogs_exp_groups = dogs_fpkms_median.merge(gene_exp_groups, how = 'outer').fillna(0).sort_values('geneFPKM', ascending = False)
        dogs_exp_groups['RT'] = dogs_exp_groups['dogFPKM'].map(lambda x: 'NRT' if x == 0 else 'RT')
        dogs_exp_groups = dogs_exp_groups.reset_index(drop = True)
        
        # use unfiltered data from ARTDeco output to classify dubious cases
        dogs_raw = pd.read_table(artdeco_dir + tissue_info[0] + '/dogs/all_dogs.raw.txt')['ID']
        # get list of dubious cases by subtracting our filtered dogs from the raw data
        all_filtered_dogs = tissue_info[4]
        dogs_dubious = dogs_raw[~dogs_raw.isin(all_filtered_dogs)]
        # make sure dubious cases are expressed in at least 25% of the samples
        dogs_dubious = dogs_dubious[dogs_dubious.isin(tissue_info[2])]
         
        # label dubious genes overlapping NRT label for those genes
        dogs_dubious_label = dogs_exp_groups['gene_id'].isin(dogs_dubious)
        dogs_exp_groups.loc[dogs_dubious_label, 'RT'] = 'UND'
         
        # add label with tissue name
        dogs_exp_groups['tissue'] = tissue_info[0]
        
        # combine into a data frame
        all_expression_data_df = pd.concat([all_expression_data_df, dogs_exp_groups])
        
    return all_expression_data_df

if __name__ == "__main__":
	main()