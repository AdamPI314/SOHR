#!/usr/bin/env python
import os
import re

def check_string(filename, str):
	#no need to pass arguments to function if you're not using them
	w = str

	#open the file using `with` context manager, it'll automatically close the file for you
	with open(filename) as f_handle:
		found = False
		for line in f_handle:
			#string or substring, any position
			if re.search(".*{0}.*".format(w),line):
				found = True
				break
	return found


def check_condition(filename, str):
	#file found
	if(os.path.isfile(filename)):
		#contains the string
		if check_string(filename, str):
			return True
		#doesn't contain the string
		else:
			return False
	#file not found
	else:
		return False


##	main test
#if __name__== "__main__":
#	if check_condition("USDOI.log", "Libertyon"):
#		print "file contains string!"
#	else:
#		print "file doesn't exit or doesn't contain the string!"
