# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 23:25:28 2024

@author: paulo
"""

import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import scipy
from statsmodels.stats.multitest import fdrcorrection as fdr

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

# import gene annotation file
all_genes_bed = pd.read_table('utils/human.gencode.v37.genes.bed', header= None)
all_genes_bed.columns = ['chrm','start','end','gene_id','gene_name','strand','gene_type']

# import master table
master_table = pd.read_table('master.table.GTEx.samples.RT.analysis.tab')
master_table_gene_info = master_table.merge(all_genes_bed)
 

# ----- # ----- # ----- # ----- # ----- # ----- # ----- # ----- #
# REPRODUCE FIG 4e
# Explore RT Sequences downstream genes

