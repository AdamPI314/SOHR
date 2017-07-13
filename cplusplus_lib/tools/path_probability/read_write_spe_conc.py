#!/usr/bin/env python
import numpy as np
import os

def read_spe_conc(filename):
	spe_conc= np.loadtxt(filename, delimiter=",")
	return spe_conc

def read_write_spe_conc(filename_in, filename_out):
	spe_conc= read_spe_conc(filename_in)
	with open(filename_out, 'a') as f_handle:
		np.savetxt(f_handle, spe_conc, fmt= '%.15e');

#file_dir= "/home/invictus/Documents/eclipse/RxnNetwork_TDDM_alpha_14.7"
#filename_in= os.path.join(file_dir,"output/spe_conc.csv")
#filename_out=os.path.join(file_dir,"output/spe_conc_all.csv")
#read_write_spe_conc(filename_in, filename_out)
