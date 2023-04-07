import pandas as pd
import os, glob
import sys

tissue = sys.argv[1] # artdeco folder, e.g "path/to/liver-samples/artdeco_out/"

suffix = ".filtered.intersect"
folder_to_filter = tissue + '/dogs/'
new_folder = tissue + '/dogs-filtered/'

for file in glob.glob('{}*{}'.format(folder_to_filter, suffix)):
    
    print('CREATING {} FILES'.format(os.path.basename(file)[:-len(suffix)-3]))
    
    dogs_filtered_ids = pd.read_table(file, header = None)[3]
    
    # create a folder for the tissue if it does not exist yet
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)
    
    # create a new fpkm.txt file filtered
    dogs_fpkm_file_name = os.path.basename(file)[:-len(suffix)-3] + 'fpkm.txt'
    dogs_fpkm_file = pd.read_table(folder_to_filter + dogs_fpkm_file_name)
    dogs_fpkm_file_filtered = dogs_fpkm_file[dogs_fpkm_file['ID'].isin(dogs_filtered_ids)]
    dogs_fpkm_file_filtered.to_csv(new_folder + dogs_fpkm_file_name, sep = '\t', index = False)

    # create a new raw.txt file filtered
    dogs_raw_file_name = os.path.basename(file)[:-len(suffix)-3] + 'raw.txt'
    dogs_raw_file = pd.read_table(folder_to_filter + dogs_raw_file_name)
    dogs_raw_file_filtered = dogs_raw_file[dogs_raw_file['ID'].isin(dogs_filtered_ids)]
    dogs_raw_file_filtered.to_csv(new_folder + dogs_raw_file_name, sep = '\t', index = False)
    
    # create a new .bed file file filtered
    dogs_bed_file_name = os.path.basename(file)[:-len(suffix)-3] + 'bed'
    dogs_bed_file = pd.read_table(folder_to_filter + dogs_bed_file_name, header = None)
    dogs_bed_file_filtered = dogs_bed_file[dogs_bed_file[3].isin(dogs_filtered_ids)]
    dogs_bed_file_filtered.to_csv(new_folder + dogs_bed_file_name, sep = '\t', index = False, header = None)
    
    os.remove(file)