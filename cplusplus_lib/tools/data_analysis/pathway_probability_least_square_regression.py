# -*- coding: utf-8 -*-

import numpy as np
import os

import my_utility as mu
import fit_least_square_regression as flsr
import pathway_probability as pp


# # Import uncertainty data
# * $k_i=2\times\frac{k_i-k_i^{MIN}}{k_i^{MAX}-k_i^{MIN}}-1$
# * in the range of (-1, 1)

def least_square_regression(file_dir, species_type, filename_in, top_N=10, N_pathway=50):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print file_dir
    uncertainty = mu.uncertainty_c(file_dir)

    # # Pathway Probability and Species concentration conversion factor
    # * Initial concentration $H_2$= 2/3, $O_2$= 1/3
    # * $\frac{\#H_2O}{\#H} \equiv$ $\frac{conc~of~H_2O}{conc~of~H_2 \times 2} \equiv$ $\frac{conc~of~H_2O}{2/3 \times 2}
    # \equiv$ $\frac{\#H_2O}{\#tracjectories} \equiv$ $\frac{\#trajectories~terminate~ with~ H_2O/2}{\#tracjectories}
    # \equiv$ $\frac{Probability ~of~Pathway~terminated~with~H_2O}{2}$
    # # $conc~of~H_2O \equiv$ $Probability ~of~Pathway~terminated~with~H_2O$ $\times \frac{2}{3}$

    name_factor = {'H2O': 2, 'H2': 2, 'H2O2': 2, 'H': 1, 'OH': 1, 'HO2': 1}

    # # Import pathway prob data-definition
    # * $P(k_1, k_2, ..., k_{n_k})=P_0+\sum_{i=1}^{n_k}{A_i(k_i)}+\sum_{i<j}^{n_k}{\sum_{j=1}^{n_k}{B_{ij}(k_i,k_j)}}$

    #    file_name_pathway_prob= os.path.join(file_dir, "output", "pathway_"+species_type+".csv")
    #
    #    pathway_frame_object= pp.pathway_prob_c.read_pathway_as_pandas_frame_object(file_name_pathway_prob)
    #    pathway_frame_object_list= pp.pathway_prob_c.convert_pandas_frame_object_to_list(pathway_frame_object)
    #
    #    all_pathway_name, all_pathway_prob= pp.pathway_prob_c.get_name_prob_of_all_pathway(pathway_frame_object)
    filename_pathway_prob = os.path.join(file_dir, "output", "pathway_prob_all.csv")
    pathway_prob_topN = np.loadtxt(filename_pathway_prob)
    # pathway_prob_topN= np.reshape(pathway_prob_topN, [-1, 50])
    pathway_prob_topN = np.reshape(pathway_prob_topN, [-1, N_pathway])
    print np.shape(pathway_prob_topN)

    # # Pathway PP CDF
    # * CDF

    #    curr_dir_pp= os.path.join(file_dir, "output", "dist", species_type)
    #    if not os.path.exists(curr_dir_pp):
    #        os.makedirs(curr_dir_pp)
    #    # Just use the first set of K
    #	pp.pathway_prob_c.plot_cdf(pathway_prob_topN[0,:], species_type, curr_dir_pp, 50)
    #
    #    # * Write pathway probablity data into file
    #    pathway_prob_topN= [pp.pathway_prob_c.get_list_of_prob_of_a_pathway(pathway_frame_object_list, all_pathway_name, i) \
    #                        for i in xrange(top_N)]
    #    print "pathway probability topN data shape:\t", np.shape(pathway_prob_topN)
    ####################################################################################################################
    # convert pathway probability to corresponding concentration
    ####################################################################################################################
    for i in xrange(np.shape(pathway_prob_topN)[0]):
        for j in xrange(np.shape(pathway_prob_topN)[1]):
            pathway_prob_topN[i][j] *= (4.0 / 3 / name_factor[species_type])
            ############################################################################################################
            #
            #    curr_dir_t= os.path.join(file_dir, "output", "coef")
            #    if not os.path.exists(curr_dir_t):
            #        os.makedirs(curr_dir_t)
            #    filename_t= os.path.join(curr_dir_t, species_type+"_pathway_probability_topN.csv")
            #    np.savetxt(filename_t, pathway_prob_topN, delimiter=",", fmt= "%.18e")
            #    pathway_prob_topN_t= np.loadtxt(filename_t, delimiter=",", dtype=float)

    # # Least square regression
    # * least square regression
    # * fit all variable at the same time to $2^{nd}$ order
    # * save to file

    # got to eliminate the duplicated reactions
    # index_2nd_all= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23]
    index_2nd_all = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19]

    # write coef to file
    curr_dir_coef = os.path.join(file_dir, "output", "pathway_prob_coef")
    if not os.path.exists(curr_dir_coef):
        os.makedirs(curr_dir_coef)
        # filename_t= os.path.join(curr_dir_coef, species_type+"_coef.csv")
    filename_t = os.path.join(curr_dir_coef, filename_in)
    # delete file if exist
    if os.path.isfile(filename_t):
        os.remove(filename_t)

    fit_info = np.loadtxt(os.path.join(file_dir, "tools/data_analysis", "fit_config.inp"), dtype=int);
    # print fit_info
    N_variable = fit_info[0]
    Nth_order_1st = fit_info[1]
    Nth_order_2nd = fit_info[2]
    # topN pathways
    # top_N = 3
    for i in xrange(top_N):
        ################################################################################################################
        # parse data, calculate first order sensitivity index
        # which pathway
        Nth_pathway = i
        print "fitting the Nth pathway:", i
        prob_of_a_pathway = pathway_prob_topN[:, Nth_pathway]
        # data set x
        data_all = uncertainty.data[:, index_2nd_all];
        z = prob_of_a_pathway

        # fit 1D and 2D simulataneously for 23 Ks
        fit_1D_2D_all = flsr.fit_1D_2D_all_c(data_all, z, Nth_order_1st, Nth_order_2nd)
        zero_order_t, first_order_t, second_order_t = fit_1D_2D_all.split_1D_coef_array()

        with open(filename_t, "a") as f:
            f.write("#" + str(Nth_pathway))
            f.write("\n#" + "var" + "\n")
            f.write(str(np.var(z)))
            f.write("\n#" + "zero" + "\n")
            zero_order_t.tofile(f, sep=", ", format="%s")
            f.write("\n#" + "first" + "\n")
            first_order_t.tofile(f, sep=", ", format="%s")
            f.write("\n#" + "second" + "\n")
            second_order_t.tofile(f, sep=", ", format="%s")
            f.write("\n")
