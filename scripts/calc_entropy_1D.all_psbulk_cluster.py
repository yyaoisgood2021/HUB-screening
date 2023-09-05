#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 21:24:19 2023

this script extracts the the atac peaks/read counts from the psbulk.clusters of the taiji inputs

use bedtools intersect to find counts in the 5kbp genomic bins



@author: peiyaowu
"""

# import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys


from scipy.stats import entropy


#################### atac
dt_folder = sys.argv[1] # overlap files
sv_folder = sys.argv[2] # result to save
os.makedirs(sv_folder, exist_ok=True, )

meta_df_path = sys.argv[3]
overlap_file_base = sys.argv[4]
# res = 5000


meta_df = pd.read_csv(meta_df_path,
    )

chrid_list = (['chr{}'.format(i) for i in range(1,9)] + ['chr{}'.format(i) for i in range(10,22)] + ['chrX']
    + ['before_translocg.chr9', 'before_translocg.chr22', 'after_translocg.chr9', 'after_translocg.chr22'])

sample_type_dict = dict(zip(
    list(meta_df['taiji_id']),
    list(meta_df['km.cluster']))
)



def calc_entropy_psbulk_1(sample_type, dt_folder, chrid_list, overlap_file_base):
    """
    this function return the combined dt 
    normalized thru the entire genome
    """
    # sum_hw = 0 # peak height*peak width
    dt_all = pd.DataFrame()
    for chrid in chrid_list:
        # open overlapped peak dt
        dt = pd.read_csv(os.path.join(
            dt_folder,
            '{}.{}.{}.txt'.format(overlap_file_base, chrid, sample_type)
            ), sep='\t', header=None, 
            dtype={6:'string', }
        )
        dt = dt[dt[4]!=-1].copy().reset_index(drop=True)[[1,6,7]]
        dt[6] = dt[6].astype(float)
        dt['hw'] = dt[6]*dt[7]
        dt = dt.groupby(1).agg({'hw':'sum'})
        dt = dt[dt['hw']>0].copy().reset_index(drop=True)
        dt['c'] = chrid
        dt_all = pd.concat((dt_all, dt.copy()), ignore_index=True, )
    dt_all['hw_norm'] = dt_all['hw'] / np.sum(dt_all['hw'])

    return dt_all


def calc_entropy_psbulk_2(dt_all, chrid_list):
    """
    this function returns a list of entropy values 

    """
    entropy_list = []
    for chrid in chrid_list:
        dt_all_ = dt_all[dt_all['c']==chrid].copy().reset_index(drop=True)
        entropy_list.append(entropy(dt_all_['hw_norm']))

    return entropy_list

# dt_all = calc_entropy_psbulk('noness_10', dt_folder, chrid_list)


#%%
sample_type_list = list(meta_df['taiji_id'])
N_sample = len(sample_type_list)
result_df = pd.DataFrame(data = np.zeros((N_sample, len(chrid_list)+2)), columns = ['psbulk.sample','km.cls']+chrid_list)
#%% load dt
for j, (sample_type,kmclstype) in enumerate(sample_type_dict.items()):
    print(j, (sample_type,kmclstype))
    if sample_type=='WT':
        continue

    result_df.iloc[j, :2] = [sample_type, str(kmclstype)]
    
    dt_all = calc_entropy_psbulk_1(sample_type, dt_folder, chrid_list, overlap_file_base)
    entropy_list = calc_entropy_psbulk_2(dt_all, chrid_list)

    result_df.iloc[j, 2:27] = entropy_list

    # chr_entropy = []
    # intensity_list = []
    # for chrid in chrid_list:
    #     print(chrid)
    #     sample_fname = 'atac_peaks.{}.{}.hg38.bed-{}.txt'.format(sample_type, chrid, chrid)
    #     dt = pd.read_csv(os.path.join(
    #         dt_folder, sample_fname
    #         ), sep='\t', header=None, )
    #     dt_area = get_area_array(dt, peak_intensity=True)
    #     intensity_list.extend(list(dt_area['area']))
    #     chr_entropy.append(array_2_entropy(dt_area))

    # result_df.iloc[j, :21] = chr_entropy
    # result_df.iloc[j, 21] =  array_2_entropy(pd.DataFrame({'area':intensity_list}))

# result_df['sample'] = sample_type_list
# result_df['km.cls'] = [sample_type_dict[i] for i in sample_type_list]

result_df.to_csv(os.path.join(sv_folder, 'atac_entropy.result.intensity.txt'), sep='\t', index=None, )







