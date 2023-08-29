#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:12:39 2022

@author: peiyaowu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, os

meta_df_path = sys.argv[1]
wt_rna_path = sys.argv[2]
wt_atac_rep1_path = sys.argv[3]
wt_atac_rep2_path = sys.argv[4]
sample_rna_atac_folder_path = sys.argv[5]
epitensor_result_path = sys.argv[6]
write_path = sys.argv[7]


#%% MODIFY THIS BLOCK
###################################################################################################################################
# input.yml file required by Taiji, go to the website (https://taiji-pipeline.github.io/) for more information

meta_df = pd.read_csv(meta_df_path, sep='\t')


input_dt = 'RNA-seq:\n'
# add WT RNA data
slice_input = (' '*2 + '- id: WT_RNA\n')
slice_input += (' '*4 + 'group: WT\n')
slice_input += (' '*4 + 'replicates:\n')
slice_input += (' '*6 + '- rep: 1\n')
slice_input += (' '*8 + 'files:\n')
slice_input += (' '*10 + '- path: {}\n'.format(wt_rna_path))     
slice_input += (' '*12 + 'tags: [\'GeneQuant\']\n')
input_dt += slice_input
# add sample RNA data
for clsid in meta_df['cluster_id']:
    cls_type, cls_number = clsid.split('-')
    clsid_new = '_'.join([cls_type, cls_number])
    print(clsid)
    if cls_type=='noness':
        slice_input = (' '*2 + '- id: {}_RNA\n'.format(clsid_new))
        slice_input += (' '*4 + 'group: {}\n'.format(clsid_new))
        slice_input += (' '*4 + 'replicates:\n')
        slice_input += (' '*6 + '- rep: 1\n')
        slice_input += (' '*8 + 'files:\n')
        slice_input += (' '*10 + '- path: {}/{}.rna-expr-pscounts.txt\n'.format(sample_rna_atac_folder_path, clsid))      
        slice_input += (' '*12 + 'tags: [\'GeneQuant\']\n')
        input_dt += slice_input
    if cls_type=='ess':
        for rep_id in ['orig', ]:
            clsid_new_new = clsid_new + '_sub_' + rep_id
            slice_input = (' '*2 + '- id: {}_RNA\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'group: {}\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'replicates:\n')
            slice_input += (' '*6 + '- rep: 1\n')
            slice_input += (' '*8 + 'files:\n')
            slice_input += (' '*10 + '- path: {}/{}.rna-expr-pscounts.txt\n'.format(sample_rna_atac_folder_path, clsid))      
            slice_input += (' '*12 + 'tags: [\'GeneQuant\']\n')
            input_dt += slice_input
                
input_dt += 'ATAC-seq:\n'
# add WT atac data
slice_input = (' '*2 + '- id: WT_ATAC\n')
slice_input += (' '*4 + 'group: WT\n')
slice_input += (' '*4 + 'replicates:\n')
slice_input += (' '*6 + '- rep: 1\n')
slice_input += (' '*8 + 'files:\n')
slice_input += (' '*10 + '- path: {}\n'.format(wt_atac_rep1_path))     
slice_input += (' '*12 + 'format: NarrowPeak\n')
slice_input += (' '*6 + '- rep: 2\n')
slice_input += (' '*8 + 'files:\n')
slice_input += (' '*10 + '- path: {}\n'.format(wt_atac_rep2_path))     
slice_input += (' '*12 + 'format: NarrowPeak\n')
input_dt += slice_input

# add sample atac data
for clsid in meta_df['cluster_id']:
    cls_type, cls_number = clsid.split('-')
    clsid_new = '_'.join([cls_type, cls_number])
    print(clsid)
    if cls_type=='noness':
        slice_input = (' '*2 + '- id: {}_ATAC\n'.format(clsid_new))
        slice_input += (' '*4 + 'group: {}\n'.format(clsid_new))
        slice_input += (' '*4 + 'replicates:\n')
        slice_input += (' '*6 + '- rep: 1\n')
        slice_input += (' '*8 + 'files:\n')
        slice_input += (' '*10 + '- path: {}/{}.atac-fragments.bed.gz\n'.format(sample_rna_atac_folder_path, clsid))     
        slice_input += (' '*12 + 'tags: [\'PairedEnd\']\n')
        input_dt += slice_input
    if cls_type=='ess':
        for rep_id in ['orig', ]:
            clsid_new_new = clsid_new + '_sub_' + rep_id
            slice_input = (' '*2 + '- id: {}_ATAC\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'group: {}\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'replicates:\n')
            slice_input += (' '*6 + '- rep: 1\n')
            slice_input += (' '*8 + 'files:\n')
            if rep_id=='orig':
                slice_input += (' '*10 + '- path: {}/{}.atac-fragments.bed.gz\n'.format(sample_rna_atac_folder_path, clsid))       
            slice_input += (' '*12 + 'tags: [\'PairedEnd\']\n')
            input_dt += slice_input

# add WT and sample HiC
input_dt += 'HiC:\n'
for clsid in ['WT'] + list(meta_df['cluster_id']):
    if clsid=='WT':
        clsid_new=clsid
    else:
        cls_type, cls_number = clsid.split('-')
        clsid_new = '_'.join([cls_type, cls_number])
    print(clsid)
    if cls_type!='ess':
        slice_input = (' '*2 + '- id: {}_HiC\n'.format(clsid_new))
        slice_input += (' '*4 + 'group: {}\n'.format(clsid_new))
        slice_input += (' '*4 + 'replicates:\n')
        slice_input += (' '*6 + '- rep: 1\n')
        slice_input += (' '*8 + 'files:\n')
        slice_input += (' '*10 + '- path: {}\n'.format(epitensor_result_path))
        slice_input += (' '*12 + 'tags: [\'ChromosomeLoop\']\n')
        input_dt += slice_input
    else:
        for rep_id in ['orig', ]:
            clsid_new_new = clsid_new + '_sub_' + rep_id
            slice_input = (' '*2 + '- id: {}_HiC\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'group: {}\n'.format(clsid_new_new))
            slice_input += (' '*4 + 'replicates:\n')
            slice_input += (' '*6 + '- rep: 1\n')
            slice_input += (' '*8 + 'files:\n')
            slice_input += (' '*10 + '- path: {}\n'.format(epitensor_result_path))
            slice_input += (' '*12 + 'tags: [\'ChromosomeLoop\']\n')
            input_dt += slice_input
    

f = open(write_path, 'w')
f.write(input_dt)
f.close()



