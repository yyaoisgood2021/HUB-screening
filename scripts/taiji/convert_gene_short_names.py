#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 20:30:11 2022

@author: peiyaowu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

#%% load gencode file and make a dict of enc-symbol and gene short names          

gencode_file_path = sys.argv[1]
data_to_transfer_path = sys.argv[2] 
# {prefix}-q30-ReadsPerGene.out.tab
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
data_to_save_path = sys.argv[3]


gencode = pd.read_csv(gencode_file_path, skiprows=5, sep='\t', header=None)
gencode = gencode[[2,8]]
gencode = list(gencode[gencode[2]=='gene'][8])

gdict = {}
for ln in gencode:
    ln_ = ln.split('; ')
    for ln__ in ln_:
        ln___ = ln__.split(' ')[0]
        if ln___== 'gene_id':
            gid = ln__.split('"')[1]
        if ln___== 'gene_name':
            gname = ln__.split('"')[1]
    gdict [gid] = gname
        
# convert WT's gene counts (generated by star) to another txt file      
def convert_g_names(in_file_path, save_file_path, gdict):

    dt = pd.read_csv(in_file_path, sep='\t', header=None, skiprows=4)       

    gid_list = list(dt[0])
    gname_list = [gdict[i] for i in gid_list]

    dt_ = pd.DataFrame({0:gname_list,
        1:dt[1]})

    dt_.to_csv(save_file_path, sep='\t', header=None, index=None)
    return


convert_g_names(data_to_transfer_path, data_to_save_path, gdict)




