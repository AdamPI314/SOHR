#!/usr/bin/env python
import numpy as np
import os
def read_target_time(filename):
    target_time=np.loadtxt(filename)
    return target_time

def read_write_target_time(filename_in, filename_out):
    target_time=read_target_time(filename_in)
    with open(filename_out,'a') as f_handle:
		f_handle.write(str(target_time)+"\n")

#file_dir= "/home/invictus/Documents/eclipse/RxnNetwork_TDDM_alpha_14.6"
#filename_in= os.path.join(file_dir,"output/target_time.csv")
#filename_out=os.path.join(file_dir,"output/target_time_all.csv")
#read_write_target_time(filename_in, filename_out)
