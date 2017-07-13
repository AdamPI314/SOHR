#!/usr/bin/env python
import os
import numpy as np
import sys
import re

if __name__== '__main__':
	if(len(sys.argv)<3):
		print "plz type in the appropriate input and output file name\n"
		exit()	
	
	input_filename= sys.argv[1]
	output_filename= sys.argv[2]

	path_data= np.loadtxt(input_filename, dtype=str)

	re.findall('(S(\d)+R(\d)+S(\d)+R(\d)+)+', 'S12R2424S1R10S2R38S1R10S2R38S1R10S2R38S1R10S2R38S1')
