#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

ref_file_path = sys.argv[1]
output_sv_path = sys.argv[2]
input_files = list(sys.argv[3:])


dt_all = pd.read_csv(ref_file_path, sep='\t', header=None, )
dt_all.columns = ['g', 'r']
dt_all['r']= 0

for i, f in enumerate(input_files):
    dt_new = pd.read_csv(f, sep='\t', header=None, )
    dt_new.columns = ['g', str(i)]
    dt_all = dt_all.merge(dt_new, on='g', how='left')

dt_all.fillna(0, inplace=True, )
sum_val = dt_all.iloc[:, 1:].sum(axis=1)
dt_all['sum'] = sum_val

dt_all[['g','sum']].to_csv(output_sv_path, sep='\t', header=None, index=None, )


