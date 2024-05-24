# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:31:38 2023

@author: Michael
"""


import scvelo as scv
import glob
import pandas as pd
import argparse


## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i','--input_dir', help='Directory containing input files', required=True)
parser.add_argument('-m','--metadata_dir', help='Directory containing "cellranger_aggr_cell_metadata.tsv" and \
                    "cellranger_aggr_input.csv" files', required=True)

args = parser.parse_args()
in_dir = args.input_dir
meta_dir = args.metadata_dir



### functions

def merge_loom(in_dir, meta_dir):

    barcode_info_path = glob.glob(meta_dir + '/cellranger_aggr_cell_metadata*.tsv')
    barcode_info = pd.read_table(*barcode_info_path, header=0)
    loom_files = glob.glob(in_dir + '/*.loom')
    
    if in_dir + '/merged.loom' in loom_files:
        loom_files.remove(in_dir + '/merged.loom')
    
    # loop over all loom files in in_dir
    count = 1
    for loom in loom_files:
        
        ldata = scv.read(loom, cache=True)
        
        # rename cell barcodes
        sample_name = loom.split('/')[-1].split('.')[0]
        print('Processing: ' + sample_name)
        sample_barcode = barcode_info.loc[barcode_info['sample_id']==sample_name,'barcode'].iloc[0]
        sample_number = sample_barcode.split('-')[-1]
        print('Sample number: ' + sample_number)
        
        barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
        barcodes = [bc[0:len(bc)-1] + '-' + sample_number for bc in barcodes]
        ldata.obs.index = barcodes
        
        # make variable names unique
        ldata.var_names_make_unique()
        
        if count == 1:
            ldata_merge = ldata
        else:
            ldata_merge = ldata_merge.concatenate([ldata])
        
        #print(ldata_merge.obs)
        count += 1
    
    # reformat cell barcodes in final file
    index_part_1 = [bc.split('-')[0] for bc in ldata_merge.obs.index.tolist()]
    index_part_2 = [bc.split('-')[1] for bc in ldata_merge.obs.index.tolist()]
    index_new = [i + '-' + j for i, j in zip(index_part_1, index_part_2)]
    ldata_merge.obs.index = index_new
    print(ldata_merge.obs)
    
    ldata_merge.write_loom(in_dir + '/merged.loom')
    
    return(ldata_merge)
    



### analysis
loom_merged = merge_loom(in_dir=in_dir, meta_dir=meta_dir)



