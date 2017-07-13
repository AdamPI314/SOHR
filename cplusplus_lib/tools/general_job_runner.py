import time
import os
import sys

path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'queues_check_exit'))
if not path in sys.path:
	    sys.path.append(path)
path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'path_probability'))
if not path in sys.path:
	    sys.path.append(path)

import job_check_condition as JCC
import evaluate_path_prob as EPP

if __name__== "__main__":
	time_beg= time.time()
	#terminated time if status is off
	criteria_time= 60*60
	#check period
	period= 60
	flag= False

	while((not flag) and ((time.time()-time_beg)<criteria_time)):
		flag= JCC.check_condition("./input/status.log", "on")
		if (flag):
			break
		time.sleep(period)

	#if status is on, run the pathway job
	if(flag):
		EPP.pathway_job_runner()
