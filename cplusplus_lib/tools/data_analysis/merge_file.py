#!/usr/bin/env python
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import sys
# import re
# import numpy as np
import pandas as pd
# import scipy as sp
# import matplotlib.pyplot as plt
import shutil


# <codecell>
def execute(file_dir):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print file_dir

    # <codecell>

    file_dir1 = "/home/invictus/Documents/Summer_2013_Argonne/DATA_WORKDIR/TDDM_var_alpha_Mar_25_13_36_07_2014_PP"
    file_dir2 = "/home/invictus/Documents/Summer_2013_Argonne/DATA_WORKDIR/TDDM_var_alpha_Mar_25_10_15_51_2014_PP"
    # file_dir1= "D:\Documents\Documents_to_research\Summer_2013_Argonne\TDDM_var_alpha_Mar_25_13_36_07_2014_PP"
    # file_dir2= "D:\Documents\Documents_to_research\Summer_2013_Argonne\TDDM_var_alpha_Mar_25_10_15_51_2014_PP"

    # <markdowncell>

    # # copy from file_dir2 to file_dir
    # * check whether file_dir contains the files gonna to be transfered

    # <codecell>

    for filename in os.listdir(os.path.join(file_dir, "output")):
        if not filename.find("_1") == -1:
            os.remove(os.path.join(file_dir, "output", filename))

    for filename in os.listdir(os.path.join(file_dir1, "output")):
        shutil.copyfile(os.path.join(file_dir1, "output", filename), \
                        os.path.join(file_dir, "output", filename.split(".")[0] + "_1." + filename.split(".")[1]))

    # <markdowncell>

    # # copy from file_dir2 to file_dir
    # * check whether file_dir contains the files gonna to be transfered

    # <codecell>


    for filename in os.listdir(os.path.join(file_dir, "output")):
        if not filename.find("_2") == -1:
            os.remove(os.path.join(file_dir, "output", filename))

    for filename in os.listdir(os.path.join(file_dir2, "output")):
        shutil.copyfile(os.path.join(file_dir2, "output", filename), \
                        os.path.join(file_dir, "output", filename.split(".")[0] + "_2." + filename.split(".")[1]))

    # <markdowncell>

    # * delete file whoes name contains "_3" and "_4"

    # <codecell>

    for filename in os.listdir(os.path.join(file_dir, "output")):
        if (not filename.find("_3") == -1) or (not filename.find("_4") == -1):
            os.remove(os.path.join(file_dir, "output", filename))

    # <markdowncell>

    # # re-index "_2" file if needed
    # * filename of files need to be re-index

    # <codecell>

    # source files 1
    filename_src1 = []
    for filename in os.listdir(os.path.join(file_dir, "output")):
        if (not filename.find("_1") == -1) and (not filename.find("pathway") == -1) and (
        not filename.find("csv") == -1):
            filename_src1.append(filename)

    file_t = pd.read_csv(os.path.join(file_dir, "output", filename_src1[0]), delimiter="\t",
                         names=["pathway", "number"])
    mask = file_t.pathway.str.contains("#")
    MaxIndex = file_t[mask].tail(1).pathway.str.replace("#", "").tolist()[0]
    MaxIndex = int(MaxIndex) + 1
    print MaxIndex
    print filename_src1

    # source files 2
    filename_src2 = map(lambda x: x.replace("1", "2"), filename_src1)
    print filename_src2

    # source files 3
    filename_src3 = map(lambda x: x.replace("1", "3"), filename_src1)
    print filename_src3

    # write index-updated files 2 into files 3
    for i in xrange(len(filename_src2)):
        file_t = file_t = pd.read_csv(os.path.join(file_dir, "output", filename_src2[i]), \
                                      delimiter="\t", names=["pathway", "number"])
        mask = file_t.pathway.str.contains("#")
        file_t.loc[mask, "pathway"] = file_t.loc[mask, "pathway"].map(
            lambda x: "#" + str(int(x.split("#")[1]) + MaxIndex))
        file_t.to_csv(os.path.join(file_dir, "output", filename_src3[i]), encoding='utf-8', index=False, sep="\t",
                      header=False)

    # <markdowncell>

    # # other src files that have no need to re-index
    # * no re-index
    # * file name

    # <codecell>

    # other src file 1
    filename_src_other_1 = []
    for filename in os.listdir(os.path.join(file_dir, "output")):
        if (not filename.find("_1") == -1) and (filename.find("csv") == -1 or filename.find("pathway") == -1):
            filename_src_other_1.append(filename)
    print filename_src_other_1

    # other src file 2
    filename_src_other_2 = map(lambda x: x.replace("_1", "_2"), filename_src_other_1)
    print filename_src_other_2

    # <markdowncell>

    # * merge file

    # <codecell>

    def merge_files(in_file1, in_file2, out_file):
        with open(out_file, "w") as outfile:
            with open(in_file1, "r") as infile1:
                for line1 in infile1:
                    outfile.write(line1)
            with open(in_file2, "r") as infile2:
                for line2 in infile2:
                    outfile.write(line2)

    for i in xrange(len(filename_src1)):
        filename_t = filename_src1[i].replace("_1", "_4")
        print filename_src1[i], filename_src3[i], filename_t
        merge_files(os.path.join(file_dir, "output", filename_src1[i]), \
                    os.path.join(file_dir, "output", filename_src3[i]), \
                    os.path.join(file_dir, "output", filename_t))
    for i in xrange(len(filename_src_other_1)):
        filename_t = filename_src_other_1[i].replace("_1", "_4")
        print filename_src_other_1[i], filename_src_other_2[i], filename_t
        merge_files(os.path.join(file_dir, "output", filename_src_other_1[i]), \
                    os.path.join(file_dir, "output", filename_src_other_2[i]), \
                    os.path.join(file_dir, "output", filename_t))

    # <markdowncell>

    # # delete files whose names contain "_1.", "_2." or "_3."

    # <codecell>

    for filename in os.listdir(os.path.join(file_dir, "output")):
        if (not filename.find("_1.") == -1) or (not filename.find("_2.") == -1) or (not filename.find("_3.") == -1):
            os.remove(os.path.join(file_dir, "output", filename))

    # <markdowncell>

    # # rename files whose names contain "_4"
    # * "_4"

    # <codecell>

    for filename in os.listdir(os.path.join(file_dir, "output")):
        if not filename.find("_4") == -1:
            os.rename(os.path.join(file_dir, "output", filename),
                      os.path.join(file_dir, "output", filename.replace("_4", "")))


if __name__ == "__main__":
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    execute(file_dir)
