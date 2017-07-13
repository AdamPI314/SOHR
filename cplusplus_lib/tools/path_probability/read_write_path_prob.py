#!usr/bin/env/ python
import numpy as np
import os

def read_path_prob(filename):
	path_prob=np.loadtxt(filename)
	return path_prob

def read_write_path_prob(filename_in, filename_out):
	path_prob= read_path_prob(filename_in)
	with open(filename_out, 'a') as f_handle:
		#np.savetxt(f_handle, path_prob, fmt='%f')
		np.savetxt(f_handle, path_prob, fmt='%.15e')

#file_dir= "/home/invictus/Documents/eclipse/RxnNetwork_TDDM_alpha_14.6"
#filename_in= os.path.join(file_dir,"output/pathway_prob.csv")
#filename_out=os.path.join(file_dir,"output/pathway_prob_all.csv")
#read_write_path_prob(filename_in, filename_out)
