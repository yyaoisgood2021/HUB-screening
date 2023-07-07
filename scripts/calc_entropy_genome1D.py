#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 21:24:19 2023
this script calculates entropy for each chr based on the input `overlapped data` 

to calculate entropy, you need to: 

1. prepare a `bedgraph file` for your signal (experiment)
format of `bedgraph file`: 
    four-column bed file with chrid, start, end, signal_value (intensity, for read count file, this column is universally 1)
    chr1 10000 10015 20
    chr1 10008 10025 27
    ...
    chr4 10234 10270 15

2. prepare a `coord file` for your genomic loci of interest
format of `coord file`: 
    three-column bed file of equal-sized (resolution) bins with no overlaps, I typically remove centromere and telomere regions
    chr1 10000 15000 
    chr1 15000 20000 
    ...
    chr4 15000 20000 

3. run bedtools intersect to map your signals to each genomic bin with the command:
    bedtools intersect -wao -a `coord file` -b `bedgraph file` > `overlap file save path`
you will get a file of the overlapped results, and it should have 8 columns

4. run this script the command:
    python ./calc_entropy_genome1D.py `overlap file save path` `entropy result save path`



@author: peiyaowu
"""

import numpy as np
import pandas as pd
import os, sys

# import pybedtools
from scipy.stats import entropy

# load dt
dt = pd.read_csv(sys.argv[1], sep='\t', header=None, dtype={6:'string', }, )
dt.columns = dt.columns.astype(str)

# normalize thru the entire genome
dt = dt = dt[dt['4']!=-1].copy().reset_index(drop=True) # chrid, coord_s, intensity, width
dt['6'] = dt['6'].astype(float)
dt['hw'] = dt['6']*dt['7']

dt = dt.groupby(['0', '1']).agg({'hw':'sum'}).reset_index()

dt['hw_norm'] = dt['hw'] / np.sum(dt['hw'])

# calc entropy for each chr
chrid_hist = []
entr_hist = []
print(dt.columns)
for chrid, dt_ in dt.groupby('0'):
    chrid_hist.append(chrid)
    entr_hist.append(
        entropy(dt_['hw_norm'])
        )

result = pd.DataFrame({
    'chr':chrid_hist,
    'entropy':entr_hist
})

print(result)

result.to_csv(sys.argv[2], sep='\t', index=None, )


