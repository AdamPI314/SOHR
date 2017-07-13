#/usr/bin/env python
import subprocess as sp
import numpy as np
import re
import signal

def check_host():
	cmd= "echo `hostname`"		
	pid=sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr= sp.PIPE)
	out, err= pid.communicate()

	if "linus" in out:
		host= "linus"
	elif "pople" in out:
		host= "pople"
	pid.wait()

	return host

def check_nodes():
	#cmd= "python sr_check_nodes.py"
	cmd= "/tcghome/sbai/sr_tools/my_bin/check-nodes-timeout.sh -t 300"
	pid=sp.Popen(cmd, shell= True, stdout=sp.PIPE, stderr= sp.PIPE)
	out, err= pid.communicate()
	return out

def parse_and_w2f(filename, check_nodes_out, criteria):
	#with open("check_nodes_out.log", 'w') as f:
	#	f.write(check_nodes_out)
	#f.close()

	#with open("check_nodes_out.log", 'r') as f:
	#	check_nodes_out=f.read()
	#f.close()

	#match match the whole string from ^ to $
	#search return the first match condition
	#findall return all match condition

	#m = re.search("(compute-\d+-\d+|pop\d+):\s+.+load average:\s+(\d+.\d+),\s+(\d+.\d+),\s+(\d+.\d+)",
	#		check_nodes_out)
	m = re.findall("(compute-\d+-\d+|pop\d+):\s+.+load average:\s+(\d+.\d+),\s+(\d+.\d+),\s+(\d+.\d+)",
			check_nodes_out)
	if m:
		#print m.groups()
		avail_nodes= [node for node in m if float(node[-3])<=criteria and float(node[-2])<=criteria and float(node[-1])<=criteria]

	#got to sort node based on nodes[-1]+nodes[-2]+nodes[-3]
	avail_nodes= sorted(avail_nodes, key= lambda x:float(x[-1])+float(x[-2])+float(x[-3]))

	with open(filename, 'w') as f_host:
		f_host.write("# This is a sample host file\n")
		for node in avail_nodes:
			f_host.write(node[0]+":4     # The next 4 procs run on this host, "+node[-3]+" "+node[-2]+" "+node[-1]+"\n")
	f_host.close()


if __name__== "__main__":
	#print check_host()
	print "check nodes (walltime=300 seconds)..."	
	print "might take longer on pople...\n"	
	check_nodes_out=[]
	check_nodes_out= check_nodes()

	criteria=15.0; filename= "hosts"
	print "search nodes with current load less than ", criteria, "..."
	print "write to file ", filename, "...\n"
	parse_and_w2f(filename, check_nodes_out, criteria)

