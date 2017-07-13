import os
import sys
import numpy as np


def path_up_to_N(n):
    S="S0"
    for i in xrange(n):
        S = S+"R"+str(i%2)+"S"+str((i+1)%2)
    S +="\n"
    return S

if __name__ == '__main__':
    N = 1000+1
    print os.getcwd()
    file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "input"))
    print file_dir

    with open(os.path.join(file_dir, "pathway_name.csv"), 'w') as f:
        for j in xrange(N):
            f.write(path_up_to_N(j))

    print os.listdir(file_dir)