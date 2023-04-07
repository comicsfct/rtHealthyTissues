# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 02:04:08 2023

@author: paulo
"""

import pandas as pd

def CombineGeneTypes(genetype_counts):     
    
    genetype_counts_new = pd.DataFrame([])
    
    # combine pseudogene types into one category
    pseudo_list = [genetype for genetype in genetype_counts.columns if 'pseudogene' in genetype]
    pseudo_list_sum = genetype_counts[pseudo_list].sum(axis = 1)                                   
    genetype_counts = genetype_counts.drop(pseudo_list, axis = 1)                              
    genetype_counts['pseudogenes'] = pseudo_list_sum
    
    # combine small RNAs into one category
    small_rna_notation = ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA']
    small_rnas_list = [genetype for genetype in genetype_counts.columns if genetype in small_rna_notation]
    small_rnas_list_sum = genetype_counts[small_rnas_list].sum(axis = 1)
    genetype_counts = genetype_counts.drop(small_rnas_list, axis = 1)
    genetype_counts['small_RNAs'] = small_rnas_list_sum

    # combine all other types into one
    others_notation = small_rna_notation + pseudo_list + ['protein_coding', 'lncRNA', 'small_RNAs', 'pseudogenes']
    others_list = [genetypes for genetypes in genetype_counts.columns if genetypes not in others_notation]
    others_list_sum = genetype_counts[others_list].sum(axis = 1) 
    genetype_counts = genetype_counts.drop(others_list, axis = 1)
    genetype_counts['other_genes'] = others_list_sum

    # add to the main dataframe
    genetype_counts_new = pd.concat([genetype_counts_new, genetype_counts], axis = 1)
    
    # re-shape
    genetype_counts_new =  genetype_counts_new[['protein_coding','lncRNA','pseudogenes','small_RNAs','other_genes']]
    
    return genetype_counts_new