# -*- coding: utf-8 -*-
import numpy as np
import os

import my_utility as mu
import fit_least_square_regression as flsr

import conc_error as ce


def least_square_regression(file_dir, species_type="H2O", filename_conc_in="./", tau=0.8):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    uncertainty = mu.uncertainty_c(file_dir)
    # print "file_dir:\t", file_dir
    print "uncertainty data shape:\t", np.shape(uncertainty.data)

    # # Import target data-target time
    # * $\tau(k_1, k_2, ..., k_{n_k})=\tau_0+\sum_{i=1}^{n_k}{A_i(k_i)}+\sum_{i<j}^{n_k}{\sum_{j=1}^{n_k}{B_{ij}(k_i,k_j)}}$

    target_time = mu.target_time_c(file_dir)
    print "target data shape:\t", np.shape(target_time.data)
    print "mean of target data:\t", np.mean(target_time.data)

    # # Import conc_error data 
    # * $C(k_1, k_2, ..., k_{n_k})=C_0+\sum_{i=1}^{n_k}{A_i(k_i)}+\sum_{i<j}^{n_k}{\sum_{j=1}^{n_k}{B_{ij}(k_i,k_j)}}$

    #    conc_error_all= ce.conc_error_c.read_conc_error(file_dir, "conc_error.csv")
    spe_conc_all = ce.conc_error_c.read_spe_conc(file_dir, "spe_conc_all.csv")
    # print np.shape(spe_conc_all)

    # # Target function, a specific species, at a specific time
    # * name_list=['H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2']
    #     * name_list=[0->'H2O', 1->'H2', 2->'H2O2', 3->'H', 4->'OH', 5->'HO2']
    # * definition of concentration
    #     * initial concentration H2->2.0/3, O2->1.0/3, H->4.0/3
    #     * relative concentration= $\frac{concentration}{total ~number ~of ~H ~atoms \times 3.0/4}$
    name_index = {"H2O": 0, "H2": 1, "H2O2": 2, "H": 3, "OH": 4, "HO2": 5}
    ith_species = int(name_index[species_type])
    # lsode_z= ce.conc_error_c.get_lsode_conc_ith_species(conc_error_all, ith_species, period= 6)
    lsode_z = ce.conc_error_c.get_lsode_conc_ith_species_v2(spe_conc_all, ith_species, period=6)
    print "lsode conc data shape:\t", np.shape(lsode_z)

    #    index_2nd_all= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23]
    index_2nd_all = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19]

    # tau=str(0.8);
    lsode_or_pathway = "lsode";
    spe_index = 0;
    time_index = None
    data_all = uncertainty.data[:, index_2nd_all];
    print np.shape(data_all)

    #    #########################################################################################################
    #    #fit 1D and 2D simulataneously for 19 Ks
    print "fitting..."

    # fit_1D_2D_all= flsr.fit_1D_2D_all_c(data_all, lsode_z, 5, 2)
    fit_info = np.loadtxt(os.path.join(file_dir, "tools/data_analysis", "fit_config.inp"), dtype=int);
    # print fit_info
    N_variable = fit_info[0];
    Nth_order_1st = fit_info[1];
    Nth_order_2nd = fit_info[2]

    fit_1D_2D_all = flsr.fit_1D_2D_all_c(data_all, lsode_z, Nth_order_1st, Nth_order_2nd)
    print "fit DONE! extract data to zero_order, first_order and second_order..."
    zero_order_t, first_order_t, second_order_t = fit_1D_2D_all.split_1D_coef_array()

    print "writting coef to file..."
    # write coef to file
    curr_dir_coef = os.path.join(file_dir, "output", "spe_conc_coef")
    if not os.path.exists(curr_dir_coef):
        os.makedirs(curr_dir_coef)
    # filename_t= os.path.join(curr_dir_coef, species_type+"_"+lsode_or_pathway+"_"+tau+"_conc_coef.csv")
    filename_t = os.path.join(curr_dir_coef, filename_conc_in)
    # delete file if exist
    if os.path.isfile(filename_t):
        os.remove(filename_t)

    with open(filename_t, "w") as f:
        f.write("#" + species_type + " " + str(spe_index) + " " + str(tau) + " " + str(time_index) + "\n")
        f.write("#" + "var" + "\n")
        f.write(str(np.var(lsode_z)))
        f.write("\n#" + "zero" + "\n")
        zero_order_t.tofile(f, sep=", ", format="%s")
        f.write("\n#" + "first" + "\n")
        first_order_t.tofile(f, sep=", ", format="%s")
        f.write("\n#" + "second" + "\n")
        second_order_t.tofile(f, sep=", ", format="%s")
        f.write("\n")

        # pathway_z= ce.conc_error_c.get_pathway_conc_ith_species(conc_error_all, 0, period= 6)
        # print "pathway conc data shape:\t", np.shape(pathway_z)
        #
        # species_type="H2O"; tau=str(0.8); lsode_or_pathway= "pathway"; spe_index=0; time_index=None
        # data_all= uncertainty.data[:, index_2nd_all];

#    #########################################################################################################
#    #fit 1D and 2D simulataneously for 23 Ks
#    print "fit..."
#    fit_1D_2D_all= flsr.fit_1D_2D_all_c(data_all, pathway_z, 5, 2)
#    print "fit DONE! extract data to zero_order, first_order and second_order..."
#    zero_order_t, first_order_t, second_order_t= fit_1D_2D_all.split_1D_coef_array()
#    
#    print "write coef to file..."
#    # write coef to file
#    curr_dir_coef= os.path.join(file_dir, "output", "conc_error_coef")
#    if not os.path.exists(curr_dir_coef):
#        os.makedirs(curr_dir_coef)
#    with open(os.path.join(curr_dir_coef, species_type+"_"+lsode_or_pathway+"_"+tau+"_conc_coef.csv"), "w") as f:
#        f.write("#"+species_type+" "+str(spe_index)+" "+tau+" "+str(time_index)+"\n")
#        f.write("#"+"var"+"\n")
#        f.write(str(np.var(pathway_z)))
#        f.write("\n#"+"zero"+"\n")
#        zero_order_t.tofile(f, sep=", ", format="%s")
#        f.write("\n#"+"first"+"\n")
#        first_order_t.tofile(f, sep=", ", format="%s")
#        f.write("\n#"+"second"+"\n")
#        second_order_t.tofile(f, sep=", ", format="%s")
#        f.write("\n")
#
