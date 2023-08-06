#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:28:43 2023

this script find the variable edges for a given TF



@author: peiyaowu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

from scipy.stats import kruskal

#%%

ext_number = 5 # the edge must be in >= `ext_number` psbulk samples in any one kmcls, else we remove this edge
no_exist_val = 0 # if this edge does not exist in some psbulk sample, we set the perc rank value `no_exist_val`

def split_list_kmcls(lst, meta_df):
    # input is a list with a length = len(meta_df), split this list based on number of each kmcls
    # for expected psbulk_number = [8,7,6,6,7] for 5 km.cls
    
    cls_ids = np.array(meta_df['km.cls']).astype(int)
    cls_id_max = np.max(cls_ids) # Num = cls_id_max + 1
    
    splitted_lst = [] # the results to return
    cummulated_num = 0
    
    for _i in range(cls_id_max+1):
        _n = (cls_ids==_i).sum()
        splitted_lst.append(
            lst[cummulated_num:cummulated_num+_n]
            )
        cummulated_num = cummulated_num + _n
    assert cummulated_num==len(meta_df)
    return splitted_lst


tf_id = int(sys.argv[1])
base_folder = sys.argv[2]

TF_of_interest = pd.read_csv(os.path.join(base_folder, 'Network_info', 'step1.filtered_edges.TF', 'TF_of_interest.txt'), header=None)
TF_of_interest = list(TF_of_interest[0])

if tf_id < len(TF_of_interest):
    TF = TF_of_interest[tf_id]

    meta_df = pd.read_csv(os.path.join(base_folder, 'cls-5.rep-0', 'meta_df.2.txt'),)
    kmcls_dict = dict(zip(meta_df['taiji_id_old'],meta_df['km.cls'],))
    kmcls_psbulk_number_dict = {i:len(meta_df[meta_df['km.cls']==i]) for i in range(5)}

    fname = os.path.join(base_folder, 'Network_info', 'step2.combined_edges', 'combined_edges.{}.txt'.format(TF))

    edge_df = pd.read_csv(fname, index_col=0)
    edge_df['km.cls'] = [kmcls_dict[i] for i in edge_df['psbulk']]

    end_hist = []
    result_hist = []

    for i, (reg, df_) in enumerate(edge_df.groupby('End')):
        # if i>100:
        #     break
        # print(reg, df_)
        # this edge must be in at least one km.cls for >= `ext_number` times (eligible edge requirements). else, it means that this edge is too weak (not exist, do not need to consider)
        l_list = [len(df__) for _, df__ in df_.groupby('km.cls')] # length list
        if max(l_list) < ext_number:
            continue
        # else, this edge is eligible in at least one km.cls, then we should check this edge's change using kruskal-wallis test
        end_hist.append(reg) # append regulatee to history
        edge_result = [] 
        for kmcls in range(5): # 0WT,1common,2common,3noness,4ess
            ext_list = list(df_[df_['km.cls']==kmcls]['perc'])
            edge_result.extend(
                ext_list + [no_exist_val]*(kmcls_psbulk_number_dict[kmcls]-len(ext_list)))
        result_hist.append(edge_result)

    result_df = pd.DataFrame(data = np.array(result_hist)) # column will be 0~33, will use this index to perform kruskal-wallis test
    result_df['end'] = end_hist
    result_df['start'] = [TF]*len(result_df)

    test_p_hist = []
    test_stat_hist = []
    for i in range(len(result_df)):
        s,p = kruskal(
            *split_list_kmcls(list(result_df.iloc[i,:len(meta_df)].values), meta_df)
            )
        test_p_hist.append(p)
        test_stat_hist.append(s)

    result_df['kruskal_p'] = test_p_hist
    result_df['kruskal_stat'] = test_stat_hist

    # result_df = result_df[result_df['kruskal_p']<=0.05].reset_index(drop=True)

    result_df.to_csv(os.path.join(
        base_folder, 'Network_info', 'step3.changed_edges.unadjusted', 'raw.eligible_edges.{}.txt'.format(TF)),
        sep='\t', index=None)










