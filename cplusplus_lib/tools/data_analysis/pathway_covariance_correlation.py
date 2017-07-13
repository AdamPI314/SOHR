# -*- coding: utf-8 -*-

import numpy as np
import os
import math

import pathway_probability as pp

import fit_least_square_regression as flsr
import variance_correlation_SI as vcs
import parse_regression_coef as prc

import sensitivity_plot as splot


# # import conc coef and var from lsode calculation

def calculate_plot_covariance_correlation(file_dir, species_type="H2O", filename_conc_in="./"):
    # * PATHWAY name, prob

    file_name_PATHWAY_prob = os.path.join(file_dir, "output", "pathway_" + species_type + ".csv")

    # PATHWAY_frame_object= pp.pathway_prob_c.read_pathway_as_pandas_frame_object(file_name_PATHWAY_prob)

    # PATHWAY_prob is normalized
    # PATHWAY_name, PATHWAY_prob= pp.pathway_prob_c.get_name_prob_of_all_pathway(PATHWAY_frame_object)

    PATHWAY_name = np.loadtxt(os.path.join(file_dir, "output", "pathway_name.csv"), dtype=np.str)
    PATHWAY_prob_t = np.loadtxt(os.path.join(file_dir, "output", "pathway_prob_all.csv"))
    PATHWAY_prob_t = np.reshape(PATHWAY_prob_t, [-1, 50]);
    # calculate the average of pathway probability, average all sets of K's
    PATHWAY_prob = [np.average(PATHWAY_prob_t[:, i]) for i in range(np.shape(PATHWAY_prob_t)[1])]
    print np.shape(PATHWAY_prob)

    # * LSODE

    # file_name_LSODE= os.path.join(file_dir, "output", "conc_error_coef", species_type+"_lsode_0.8_conc_coef.csv")
    file_name_LSODE = os.path.join(file_dir, "output", "spe_conc_coef", filename_conc_in)
    var_coef_frame_object_LSODE = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_LSODE)
    var_coef_list_LSODE = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object_LSODE)

    var_LSODE, coef_LSODE = prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_LSODE, 0)

    # # Import pathway least square regression coef
    # * zeroth order, first order and second order
    # * variance of pathway probability

    # file_name_PATHWAY= os.path.join(file_dir, "output", "coef", species_type+"_coef.csv")
    file_name_PATHWAY = os.path.join(file_dir, "output", "pathway_prob_coef", "H2O_coef.csv")
    var_coef_frame_object_PATHWAY = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_PATHWAY)
    var_coef_list_PATHWAY = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(
        var_coef_frame_object_PATHWAY)

    var_list_PATHWAY = map(lambda x: x[0],
                           [prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_PATHWAY, i) \
                            for i in xrange(np.shape(var_coef_list_PATHWAY)[0])])
    coef_list_PATHWAY = map(lambda x: x[1],
                            [prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_PATHWAY, i) \
                             for i in xrange(np.shape(var_coef_list_PATHWAY)[0])])
    topN = np.shape(var_coef_list_PATHWAY)[0]

    # N_variable=23; Nth_order_1st=5; Nth_order_2nd=2
    # N_variable=19; Nth_order_1st=4; Nth_order_2nd=2
    fit_info = np.loadtxt(os.path.join(file_dir, "tools/data_analysis", "fit_config.inp"), dtype=int);
    # print fit_info
    N_variable = fit_info[0];
    Nth_order_1st = fit_info[1];
    Nth_order_2nd = fit_info[2]

    coef_topN_PATHWAY = [flsr.fit_1D_2D_all_c.split_1D_coef_array_static( \
        coef_list_PATHWAY[i], N_variable, Nth_order_1st, Nth_order_2nd) \
                         for i in xrange(len(coef_list_PATHWAY))]

    first_order_topN_PATHWAY = map(lambda x: x[1], coef_topN_PATHWAY)
    second_order_topN_PATHWAY = map(lambda x: x[2], coef_topN_PATHWAY)

    # # Import pathway concentration

    # curr_dir_t= os.path.join(file_dir, "output", "coef")
    # if not os.path.exists(curr_dir_t):
    #    os.makedirs(curr_dir_t)
    # filename_t= os.path.join(curr_dir_t, species_type+"_pathway_probability_topN.csv")
    # pathway_prob_topN_t= np.loadtxt(filename_t, delimiter=",", dtype=float)
    PATHWAY_prob = np.loadtxt(os.path.join(file_dir, "output", "pathway_prob_all.csv"))
    pathway_prob_topN_t = np.reshape(PATHWAY_prob, [-1, 50]);
    pathway_prob_topN_t = np.transpose(pathway_prob_topN_t)

    print "pathway probability topN data shape:\t", np.shape(pathway_prob_topN_t)

    # * Calculate pathway covariance
    # * Plot pathway covariance

    PATHWAY_covariance = np.cov(pathway_prob_topN_t)
    print "pathway covariance shape:\t", np.shape(PATHWAY_covariance)
    PATHWAY_covariance_t = []
    for i in xrange(topN):
        for j in xrange(i + 1, topN):
            PATHWAY_covariance_t.append(PATHWAY_covariance[i, j])
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_covariance_all_" + species_type + ".png")
    splot.bar3D_2D_PATHWAY_covariance_v2(PATHWAY_covariance_t, topN, topN, filename, text_in="Total")

    # * Calculate pathway correlation
    # * Plot pathway correlation coefficient

    PATHWAY_correlation = np.corrcoef(pathway_prob_topN_t)
    print "pathway probability correlation data shape:\t", np.shape(PATHWAY_correlation)

    PATHWAY_correlation_t = []
    for i in xrange(topN):
        for j in xrange(i + 1, topN):
            PATHWAY_correlation_t.append(PATHWAY_correlation[i, j])
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_correlation_all_" + species_type + ".png")
    splot.bar3D_2D_PATHWAY_correlation_v2(PATHWAY_correlation_t, topN, topN, filename, text_in="Total")

    # # SI from pathway correlation
    # * pathway correlation
    # * $\langle P_{i}P_{j}\rangle$

    # * Calculate $0^{th}$ order covariance
    # * Plot $0^{th}$ order covariance

    zeroth_order_covariance_t = []
    for i in xrange(topN):
        for j in xrange(i + 1, topN):
            zeroth_order_covariance_t.append(0.0)
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_0th_covariance_" + species_type + ".png")
    splot.bar3D_2D_PATHWAY_covariance_v2(zeroth_order_covariance_t, topN, topN, filename)

    # * Calculate $0^{th}$ order pathway correlation coefficient
    # * Plot $0^{th}$ order pathway correlation coefficient

    zeroth_order_correlation_t = []
    for i in xrange(topN):
        for j in xrange(i + 1, topN):
            zeroth_order_correlation_t.append(0.0)
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_0th_correlation_" + species_type + ".png")
    splot.bar3D_2D_PATHWAY_correlation_v2(zeroth_order_correlation_t, topN, topN, filename)

    # * Calculate $1^{st}$ order covariance
    # * write $1^{st}$ order co-variance into file "H2O_covariance_1st.csv", delimited by ","

    PATHWAY_pair_ind = []
    for i in xrange(topN):
        for j in xrange(i + 1, topN):
            PATHWAY_pair_ind.append((i, j))

    first_order_covariance = []
    for rxn_ind in xrange(N_variable):
        first_order_covariance_t = []
        for i in xrange(topN):
            for j in xrange(i + 1, topN):
                first_order_covariance_t.append(vcs.correlation_1st_c(first_order_topN_PATHWAY[i][rxn_ind], \
                                                                      first_order_topN_PATHWAY[j][rxn_ind], 1.0).data)
        first_order_covariance.append(first_order_covariance_t)

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename1 = os.path.join(curr_dir, species_type + "_1st_covariance.csv")
    np.savetxt(filename1, first_order_covariance, delimiter=",", fmt="%.18e")
    first_order_covariance_t = np.loadtxt(filename1, delimiter=",", dtype=float)

    print "first order covariance data shape:\t", np.shape(first_order_covariance_t)
    print "sum of first order covariance data, has shape:\t", sum(
        [sum(first_order_covariance_t[i]) for i in xrange(N_variable)])

    # * Plot of $1^{st}$ order covariance

    ####################################################################################################################
    # Plot
    ####################################################################################################################
    """
    * $1^{st}$ order pathway covariance plot
    """
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    # The nth pair of reactions
    Nth_item = 10
    filename1 = os.path.join(curr_dir, "PATHWAY_1st_covariance_" + species_type + "_" + str(Nth_item) + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(first_order_covariance_t[Nth_item], topN, topN, filename1)
    splot.bar3D_2D_PATHWAY_covariance_v2(first_order_covariance_t[Nth_item], topN, topN, filename1, \
                                         text_in="reaction: " + str(Nth_item))

    # * Sum $1^{st}$ order covariance
    # * $\sum _{i=1}^{25}P_m(k_i)P_n(k_i)$ for every pathway pair $(m,n)$
    #     * Where i stands for the ith reaction
    # * Plot sum of $1^{st}$ order covariance

    first_order_covariance_25 = [sum(first_order_covariance_t.T[i]) for i in
                                 xrange(np.shape(first_order_covariance_t.T)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_1st_covariance_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(first_order_covariance_25, topN, topN, filename)
    splot.bar3D_2D_PATHWAY_covariance_v2(first_order_covariance_25, topN, topN, filename, \
                                         text_in="reaction: " + "ALL")

    # * Calculate $1^{st}$ order pathway correlation
    # * Plot $1^{st}$ order pathway correlation

    first_order_correlation_t = np.ones_like(first_order_covariance_t)

    for i in xrange(np.shape(first_order_covariance_t)[0]):
        for j in xrange(np.shape(first_order_covariance_t)[1]):
            first_order_correlation_t[i][j] = first_order_covariance_t[i][j] / (math.sqrt( \
                var_list_PATHWAY[PATHWAY_pair_ind[j][0]] * var_list_PATHWAY[PATHWAY_pair_ind[j][1]]))

    ####################################################################################################################
    # Plot
    ####################################################################################################################
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    # The nth pair of reactions
    Nth_item = 10
    filename1_2 = os.path.join(curr_dir, "PATHWAY_1st_correlation_" + species_type + "_" + str(Nth_item) + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(first_order_correlation_t[Nth_item], topN, topN, filename1)
    splot.bar3D_2D_PATHWAY_correlation_v2(first_order_correlation_t[Nth_item], topN, topN, filename1_2, \
                                          text_in="reaction: " + str(Nth_item))

    # * Sum $1^{st}$ order correlation
    # * $\sum _{i=1}^{25}P_m(k_i)P_n(k_i)$ for every pathway pair $(m,n)$
    #     * Where i stands for the ith reaction
    # * Plot sum of $1^{st}$ order correlation

    first_order_correlation_25 = [sum(first_order_correlation_t.T[i]) for i in
                                  xrange(np.shape(first_order_correlation_t.T)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_1st_correlation_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_correlation_v1(first_order_correlation_25, topN, topN, filename)
    splot.bar3D_2D_PATHWAY_correlation_v2(first_order_correlation_25, topN, topN, filename, \
                                          text_in="reaction: " + "ALL")

    # * Calculate $2^{nd}$ order covariance
    # * write $2^{nd}$ order co-variance into file "H2O_covariance_2nd.csv", delimited by ","

    second_order_covariance = []
    N_pairs = np.shape(second_order_topN_PATHWAY[0])[0]
    print "N_pairs:\t", N_pairs
    for rxn_ind in xrange(N_pairs):
        index_2nd_covariance = []
        second_order_covariance_t = []
        for i in xrange(topN):
            for j in xrange(i + 1, topN):
                index_2nd_covariance.append((i, j))
                second_order_covariance_t.append(vcs.correlation_2nd_c(second_order_topN_PATHWAY[i][rxn_ind], \
                                                                       second_order_topN_PATHWAY[j][rxn_ind], 1.0).data)
        second_order_covariance.append(second_order_covariance_t)

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename2 = os.path.join(curr_dir, species_type + "_2nd_covariance.csv")
    np.savetxt(filename2, second_order_covariance, delimiter=",", fmt="%.18e")
    second_order_covariance_t = np.loadtxt(filename2, delimiter=",", dtype=float)

    print "second order covariance data shape:\t", np.shape(second_order_covariance_t)
    print "sum of second order covariance data, has shape:\t", sum(
        [sum(second_order_covariance_t[i]) for i in xrange(np.shape(second_order_covariance_t)[0])])

    # * Plot of $2^{nd}$ order covariance

    ####################################################################################################################
    # Plot
    ####################################################################################################################
    """
    * $2^{nd}$ order pathway covariance plot
    """
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    # The nth pair of reactions
    Nth_item = 10
    filename2 = os.path.join(curr_dir, "PATHWAY_2nd_covariance_" + species_type + "_" + str(Nth_item) + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(second_order_covariance_t[Nth_item], topN, topN, filename2)
    splot.bar3D_2D_PATHWAY_covariance_v2(second_order_covariance_t[Nth_item], topN, topN, filename2, \
                                         text_in="reaction: " + str(Nth_item))

    # * Sum $2^{nd}$ order covariance
    # * $\sum _{i=1}^{253}P_m(permutation_i)P_n(permutation_i)$ for every pathway pair $(m,n)$
    #     * where permutation i stands for ith reaction pair
    # * Plot sum of $2^{nd}$ order covariance

    second_order_covariance_253 = [sum(second_order_covariance_t.T[i]) for i in
                                   xrange(np.shape(second_order_covariance_t.T)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_2nd_covariance_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(second_order_covariance_253, topN, topN, filename, \
    #                                        text_in= "reaction: "+"ALL")
    splot.bar3D_2D_PATHWAY_covariance_v2(second_order_covariance_253, topN, topN, filename, \
                                         text_in="reaction: " + "ALL")

    # * Calculate $2^{nd}$ order pathway correlation
    # * Plot $2^{nd}$ order pathway correlation

    second_order_correlation_t = np.ones_like(second_order_covariance_t)
    for i in xrange(np.shape(second_order_covariance_t)[0]):
        for j in xrange(np.shape(second_order_covariance_t)[1]):
            second_order_correlation_t[i][j] = second_order_covariance_t[i][j] / (math.sqrt( \
                var_list_PATHWAY[PATHWAY_pair_ind[j][0]] * var_list_PATHWAY[PATHWAY_pair_ind[j][1]]))
    ####################################################################################################################
    # Plot
    ####################################################################################################################
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    # The nth pair of reactions
    Nth_item = 10
    filename2_2 = os.path.join(curr_dir, "PATHWAY_2nd_correlation_" + species_type + "_" + str(Nth_item) + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(second_order_correlation_t[Nth_item], topN, topN, filename1)
    splot.bar3D_2D_PATHWAY_correlation_v2(second_order_correlation_t[Nth_item], topN, topN, filename2_2, \
                                          text_in="reaction: " + str(Nth_item))

    # * Sum $2^{nd}$ order correlation
    # * $\sum _{i=1}^{253}P_m(permutation_i)P_n(permutation_i)$ for every pathway pair $(m,n)$
    #     * where permutation i stands for ith reaction pair
    # * Plot sum of $2^{nd}$ order correlation

    second_order_correlation_253 = [sum(second_order_correlation_t.T[i]) for i in
                                    xrange(np.shape(second_order_correlation_t.T)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_2nd_correlation_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_correlation_v1(second_order_correlation_253, topN, topN, filename)
    splot.bar3D_2D_PATHWAY_correlation_v2(second_order_correlation_253, topN, topN, filename, \
                                          text_in="reaction: " + "ALL")

    # * Sum $1^{st}$ and $2^{nd}$ order covariance
    # * $\sum _{i=1^{st}}^{2^{nd}}P_m(i)P_n(i)$ for every pathway pair $(m,n)$
    #     * where permutation i stands for $1^{st}$ or $2^{nd}$ term
    # * Plot sum of $1^{st}$ and $2^{nd}$ order covariance

    first_second_order_covariance_25_253 = [second_order_covariance_253[i] + first_order_covariance_25[i] \
                                            for i in xrange(np.shape(first_order_covariance_25)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_1st_2nd_covariance_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_covariance_v1(first_second_order_covariance_25_253, topN, topN, filename)
    splot.bar3D_2D_PATHWAY_covariance_v2(first_second_order_covariance_25_253, topN, topN, filename, \
                                         text_in="reaction: " + "ALL" + ", $1^{st}$ and $2^{nd}$ order")

    # * Sum $1^{st}$ and $2^{nd}$ order correlation
    # * $\sum _{i=1^{st}}^{2^{nd}}P_m(i)P_n(i)$ for every pathway pair $(m,n)$
    #     * where permutation i stands for $1^{st}$ or $2^{nd}$ term
    # * Plot sum of $1^{st}$ and $2^{nd}$ order correlation

    first_second_order_correlation_25_253 = [second_order_correlation_253[i] + first_order_correlation_25[i] \
                                             for i in xrange(np.shape(first_order_correlation_25)[0])]
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    filename = os.path.join(curr_dir, "PATHWAY_1st_2nd_correlation_" + species_type + ".png")
    # splot.bar3D_2D_PATHWAY_correlation_v1(first_second_order_correlation_25_253, topN, topN, filename)
    splot.bar3D_2D_PATHWAY_correlation_v2(first_second_order_correlation_25_253, topN, topN, filename,
                                          text_in="reaction: " + "ALL" + ", $1^{st}$ and $2^{nd}$ order")

    # # Calculate sensitivity index from pathway covariance
    # * $1_{st}$ order
    # * $2_{nd}$ order

    # * LSODE-Global sensitivity, SI=$\frac{cov(pathway)}{var(LSODE~conc)}$
