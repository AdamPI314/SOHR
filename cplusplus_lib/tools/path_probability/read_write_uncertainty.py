#!/usr/bin/env python
import numpy as np
import os
def read_uncertainty(filename):
    uncertainty=np.loadtxt(filename,delimiter=",")
    return uncertainty[:,-1]

def read_write_uncertainty(filename_in, filename_out):
    uncertainty=read_uncertainty(filename_in)
    with open(filename_out,'a') as f_handle:
        #np.savetxt(f_handle,uncertainty, fmt='%f')
        np.savetxt(f_handle,uncertainty, fmt='%.15e')

#file_dir= "/home/invictus/Documents/eclipse/RxnNetwork_TDDM_alpha_14.6"
#filename_in= os.path.join(file_dir,"output/uncertainties_random.csv")
#filename_out=os.path.join(file_dir,"output/uncertainties_random_all.csv")
#read_write_uncertainty(filename_in, filename_out)
