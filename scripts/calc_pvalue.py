import math, scipy, sys
from scipy.stats import poisson
import numpy as np
import pandas as pd


def pvalue(df, ):
	# df with 5 columns s1, s2, vc_obs, vc_o/e, vc_exp
	result_df = df.copy()
	p_list = []
	for s1, s2, obs_, _, exp_ in result_df.iloc[:].values:
		if int(s1)==int(s2):
			p_list.append(9)
		else:
			p1=poisson.logsf(obs_,exp_,loc=0)
			if np.isinf(p1):
				p_list.append(-1000)
			else:
				p_list.append(round(p1, 4))
	result_df[0] = result_df[0].astype(int)
	result_df[1] = result_df[1].astype(int)
	result_df[2] = np.around(result_df[2], 4) # obs valuues
	result_df[4] = np.around(result_df[4], 4) # exp valuues
	result_df[5] = p_list
	return result_df[[0,1,2,4,5]]


def save_inter_files(input_path, save_path):
	df = pd.read_csv(input_path, sep='\t', header=None)
	df = pvalue(df, )
	df.to_csv(save_path, sep='\t', header=None, index=None)
	return



if __name__=="__main__":
	input_path = sys.argv[1]
	save_path = sys.argv[2]
	save_inter_files(input_path, save_path)
	
		

	
