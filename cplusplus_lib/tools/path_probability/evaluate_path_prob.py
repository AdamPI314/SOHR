#!/usr/bin/env python
import read_write_uncertainty as rwu
import read_write_target_time as rwtt
import read_write_spe_conc as rwsc
import read_write_path_prob as rwpp
import os
import numpy as np
import subprocess as sp

def execute_pathway_generator(num_processor):
	cmd= "make run "+"P="+str(num_processor)
	print cmd
	pid= sp.Popen(cmd, shell= True, stdout= sp.PIPE, stderr= sp.STDOUT)
	pid.wait()


###########################################################
########################## pathway job runner ###########################
###########################################################
def pathway_job_runner():
	#file root directory
	file_dir= os.getcwd()
	#file name
	#uncertainty
	uncertainty_filename_in= os.path.join(file_dir,"output/uncertainties_random.csv")
	uncertainty_filename_out=os.path.join(file_dir,"output/uncertainties_random_all.csv")
	if os.path.isfile(uncertainty_filename_out):
		os.remove(uncertainty_filename_out)

	#target time
	tt_filename_in= os.path.join(file_dir,"output/target_time.csv")
	tt_filename_out= os.path.join(file_dir,"output/target_time_all.csv")
	if os.path.isfile(tt_filename_out):
		os.remove(tt_filename_out)
	
	#species concentration
	sc_filename_in= os.path.join(file_dir,"output/spe_conc.csv")
	sc_filename_out= os.path.join(file_dir,"output/spe_conc_all.csv")
	if os.path.isfile(sc_filename_out):
		os.remove(sc_filename_out)

	# path probability
	pp_filename_in= os.path.join(file_dir,"output/pathway_prob.csv")
	pp_filename_out= os.path.join(file_dir,"output/pathway_prob_all.csv")
	if os.path.isfile(pp_filename_out):
		os.remove(pp_filename_out)

	#num of processor
	num_processor= 8

	#number of sets of k's
	NK= 1000
	for i in xrange(NK):
		#generate pathway
		execute_pathway_generator(num_processor)
		#write uncertainties
		rwu.read_write_uncertainty(uncertainty_filename_in, uncertainty_filename_out)
		#target time
		rwtt.read_write_target_time(tt_filename_in, tt_filename_out)
		#species concentration
		rwsc.read_write_spe_conc(sc_filename_in, sc_filename_out)
		#write probability of pathway
		rwpp.read_write_path_prob(pp_filename_in, pp_filename_out)
