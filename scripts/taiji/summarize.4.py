#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:28:43 2023

this script find the variable edges for a given TF



@author: peiyaowu
"""


#%% pval correct
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

from scipy.stats import kruskal


def pminf(array):
    x = 1
    pmin_list = []
    N = len(array)
    for index in range(N):
        if array[index] < x:
            pmin_list.insert(index, array[index])
        else:
            pmin_list.insert(index, x)
    return pmin_list
#end function

def cumminf(array):
    cummin = []
    cumulative_min = array[0]
    for p in array:
        if p < cumulative_min:
            cumulative_min = p
        cummin.append(cumulative_min)
    return cummin
#end

def cummaxf(array):
    cummax = []
    cumulative_max = array[0]
    for e in array:
        if e > cumulative_max:
            cumulative_max = e
        cummax.append(cumulative_max)
    return cummax
#end

def order(*args):
    if len(args) > 1:
        if args[1].lower() == 'false':# if ($string1 eq $string2) {
            return sorted(range(len(args[0])), key = lambda k: args[0][k])
        elif list(args[1].lower()) == list('true'):
            return sorted(range(len(args[0])), key = lambda k: args[0][k], reverse = True)
        else:
            print ("%s isn't a recognized parameter" % args[1])
            sys.exit()
    elif len(args) == 1:
        return sorted(range(len(args[0])), key = lambda k: args[0][k])
#end

def p_adjust(*args):
    method = "bh"
    pvalues = args[0]
    if len(args) > 1:
        methods = {"bh", "fdr", "by", "holm", "hommel", "bonferroni", "hochberg"}
        metharg = arg[1].lower()
        if metharg in methods:
            method = metharg
    lp = len(pvalues)
    n = lp
    qvalues = []
    if method == 'hochberg':#already all lower case
        o = order(pvalues, 'TRUE')
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, (index+1)*pvalues[o[index]])
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        ro = order(o)
        qvalues = [pmin[i] for i in ro]
    elif method == 'bh':
        o = order(pvalues, 'TRUE')
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, (n/(n-index))* pvalues[o[index]])
        ro = order(o)
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        qvalues = [pmin[i] for i in ro]
    elif method == 'by':
        q = 0.0
        o = order(pvalues, 'TRUE')
        ro = order(o)
        for index in range(1, n+1):
            q += 1.0 / index;
        cummin_input = []
        for index in range(n):
            cummin_input.insert(index, q * (n/(n-index)) * pvalues[o[index]])
        cummin = cumminf(cummin_input)
        pmin = pminf(cummin)
        qvalues = [pmin[i] for i in ro]
    elif method == 'bonferroni':
        for index in range(n):
            q = pvalues[index] * n
            if (0 <= q) and (q < 1):
                qvalues.insert(index, q)
            elif q >= 1:
                qvalues.insert(index, 1)
            else:
                print ('%g won\'t give a Bonferroni adjusted p' % q)
                sys.exit()
    elif method == 'holm':
        o = order(pvalues)
        cummax_input = []
        for index in range(n):
            cummax_input.insert(index, (n - index) * pvalues[o[index]])
        ro = order(o)
        cummax = cummaxf(cummax_input)
        pmin = pminf(cummax)
        qvalues = [pmin[i] for i in ro]
    elif method == 'hommel':
        i = range(1,n+1)
        o = order(pvalues)
        p = [pvalues[index] for index in o]
        ro = order(o)
        pa = []
        q = []
        smin = n*p[0]
        for index in range(n):
            temp = n*p[index] / (index + 1)
            if temp < smin:
                smin = temp
        for index in range(n):
            pa.insert(index, smin)
            q.insert(index, smin)
        for j in range(n-1,1,-1):
            ij = range(1,n-j+2)
            for x in range(len(ij)):
                ij[x] -= 1
            I2_LENGTH = j - 1
            i2 = []
            for index in range(I2_LENGTH+1):
                i2.insert(index, n - j + 2 + index - 1)
            q1 = j * p[i2[0]] / 2.0
            for index in range(1,I2_LENGTH):
                TEMP_Q1 = j * p[i2[index]] / (2.0 + index)
                if TEMP_Q1 < q1:
                    q1 = TEMP_Q1
            for index in range(n - j + 1):
                q[ij[index]] = min(j * p[ij[index]], q1)
            for index in range(I2_LENGTH):
                q[i2[index]] = q[n-j]
            for index in range(n):
                if pa[index] < q[index]:
                    pa[index] = q[index]
            qvalues = [pa[index] for index in ro]
    else:
        print ("method %s isn't defined." % method)
        sys.exit()
    return qvalues


#%%


ext_number = 5 # the edge must be in >= `ext_number` psbulk samples in any one kmcls, else we remove this edge
no_exist_val = 0 # if this edge does not exist in some psbulk sample, we set the perc rank value `no_exist_val`

base_folder = sys.argv[1]



# km.cls
meta_df = pd.read_csv(os.path.join(base_folder, 'cls-5.rep-0', 'meta_df.2.txt'),)
kmcls_dict = dict(zip(meta_df['taiji_id_old'],meta_df['km.cls'],))
kmcls_psbulk_number_dict = {i:len(meta_df[meta_df['km.cls']==i]) for i in range(5)}

# TF
TF_of_interest = pd.read_csv(os.path.join(base_folder, 'Network_info', 'step1.filtered_edges.TF', 'TF_of_interest.txt'),
    header=None)
TF_of_interest = list(TF_of_interest[0])


# TF = sys.argv[1]

final_result_df = pd.DataFrame()

for TF in TF_of_interest:
    fname = os.path.join(base_folder, 'Network_info', 'step3.changed_edges.unadjusted', 'raw.eligible_edges.{}.txt'.format(TF))
    result_df = pd.read_csv(fname, sep='\t')
    result_df['mean_cls_0'] = np.mean(result_df.iloc[:, :8].values, axis=1)
    result_df['mean_cls_1'] = np.mean(result_df.iloc[:, 8:15].values, axis=1)
    result_df['mean_cls_2'] = np.mean(result_df.iloc[:, 15:21].values, axis=1)
    result_df['mean_cls_3'] = np.mean(result_df.iloc[:, 21:27].values, axis=1)
    result_df['mean_cls_4'] = np.mean(result_df.iloc[:, 27:34].values, axis=1)
    # find the largest difference between cls3 and 4
    result_df['diff_3_4'] = result_df['mean_cls_3'] - result_df['mean_cls_4'] 
    result_df = result_df.sort_values(by=['diff_3_4']).reset_index(drop=True)
    # result_df = result_df[(result_df['diff_3_4']>0.5) | (result_df['diff_3_4']<-0.5)].reset_index(drop=True)
    final_result_df = pd.concat((final_result_df, result_df), ignore_index=True)

p_list = list(final_result_df['kruskal_p'])
q_list = p_adjust(p_list)

final_result_df['qval'] = q_list
final_result_df_new = final_result_df[(final_result_df['diff_3_4']>0.5) | (final_result_df['diff_3_4']<-0.5)].reset_index(drop=True).copy()
final_result_df_new = final_result_df_new[final_result_df_new['qval']<=0.05].reset_index(drop=True).copy()


final_result_df_new.to_csv(os.path.join(
    base_folder, 'cls-5.rep-0', 'edge_results', 'all_significantly_changed_edges.txt'),
    sep='\t', index=None)










