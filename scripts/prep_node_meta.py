#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 21:24:19 2023

this script converts the {feature}.{chrid}.bed to node_meta

@author: peiyaowu
"""

import pandas as pd
import os, sys

coord_folder = sys.argv[1]
overlap_folder = sys.argv[2]
result_sv_folder = sys.argv[3]

chrid_list = ['chr{}'.format(i) for i in range(1,9)] + ['chr{}'.format(i) for i in range(10,22)] + ['chrX']
chrid_list = chrid_list + ['before_translocg.chr22', 'before_translocg.chr9', 'after_translocg.chr22', 'after_translocg.chr9']

chip_dt = pd.read_csv('resources/download_chip_dt.txt', sep='\t')
chip_dt['rep'] = chip_dt['rep'].astype(str)
chip_dt = chip_dt.iloc[:-1, ].copy()
chip_dt['name'] = (chip_dt['feature'] + ['.']*len(chip_dt) + 
                    chip_dt['peak_type'] + ['.rep-']*len(chip_dt) +
                    chip_dt['rep'])

for chrid in chrid_list:
    node_f = pd.read_csv(os.path.join(coord_folder, 'eligible_coordinates.{}.hg19.bed'.format(chrid)), sep='\t', header=None)
    node_f[1] = node_f[1].astype(str)
    node_f['str'] = node_f[0] + ['-']*len(node_f) + node_f[1]
    for chip in list(chip_dt['name'])+['ATAC']:
        overlap_dt = pd.read_csv(os.path.join(overlap_folder, '{}.{}.bed'.format(chip, chrid)), sep='\t', header=None, )
        no_over_site_coord = list(overlap_dt[overlap_dt[4]==-1][1].astype(str)) # no overlapping site's coordinate
        node_f[chip] = 1 # set all node as 1 (with peaks)
        node_f.iloc[node_f[1].isin(no_over_site_coord), -1] = 0 # use bool indexer to set 0 : overlap==-1 

    # save node_f
    node_f.to_csv(os.path.join(result_sv_folder, '{}.txt'.format(chrid)),
                sep='\t', index=None)

print('done')





