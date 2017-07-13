# -*- coding: utf-8 -*-

import numpy as np
import os
# import pathway_probability as pp

import fit_least_square_regression as flsr
import variance_correlation_SI as vcs
import parse_regression_coef as prc

import sensitivity_plot as splot


# # import conc coef and var from lsode calculation

def calculate_plot_pathway_individual_SI(file_dir, species_type, filename_LSODE_in, N_pathway=50):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))

    #    # * PATHWAY name, prob
    #
    #    file_name_PATHWAY_prob= os.path.join(file_dir, "output", "pathway_"+species_type+".csv")
    #
    #    PATHWAY_frame_object= pp.pathway_prob_c.read_pathway_as_pandas_frame_object(file_name_PATHWAY_prob)
    #
    # PATHWAY_prob is normalized
    # PATHWAY_name, PATHWAY_prob= pp.pathway_prob_c.get_name_prob_of_all_pathway(PATHWAY_frame_object)
    PATHWAY_name = np.loadtxt(os.path.join(file_dir, "output", "pathway_name.csv"), dtype=np.str)
    PATHWAY_prob_t = np.loadtxt(os.path.join(file_dir, "output", "pathway_prob_all.csv"))
    PATHWAY_prob_t = np.reshape(PATHWAY_prob_t, [-1, N_pathway]);
    # calculate the average of pathway probability, average all sets of K's and normalize
    PATHWAY_prob = [np.average(PATHWAY_prob_t[:, i]) for i in range(np.shape(PATHWAY_prob_t)[1])]
    print np.shape(PATHWAY_prob)
    PATHWAY_prob_total = np.sum(PATHWAY_prob)
    PATHWAY_prob = [x / PATHWAY_prob_total for x in PATHWAY_prob]

    # * LSODE

    # file_name_LSODE= os.path.join(file_dir, "output", "spe_conc_coef", "H2O_lsode_0.8_conc_coef.csv")
    file_name_LSODE = os.path.join(file_dir, "output", "spe_conc_coef", filename_LSODE_in)
    var_coef_frame_object_LSODE = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_LSODE)
    var_coef_list_LSODE = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object_LSODE)

    var_LSODE, coef_LSODE = prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_LSODE, 0)

    # # Import pathway least square regression coef
    # * zeroth order, first order and second order
    # * variance of pathway probability

    # file_name_PATHWAY= os.path.join(file_dir, "output", "coef", "H2O_coef.csv")
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

    # N_variable=23; Nth_order_1st=5; Nth_order_2nd=2
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

    # # Calculate sensitivity index from pathway probability
    # * $1^{st}$ order
    # * $2^{nd}$ order

    # * self-Global sensitivity, SI=$\frac{var(pathway)}{var(pathway)}$
    # * Write to file

    topN = np.shape(var_coef_list_PATHWAY)[0];
    Variance_1st_all_PATHWAY = [None] * topN
    Variance_2nd_all_PATHWAY = [None] * topN

    for Nth_PATHWAY in xrange(topN):
        # for Nth_PATHWAY in xrange(0, 0+1):
        ################################################################################################################
        # calculate sensitivity index, 1st order and 2nd order
        ################################################################################################################

        SI_1st_object_all_PATHWAY = [vcs.SI_1st_c(1.0, first_order_topN_PATHWAY[Nth_PATHWAY][i], \
                                                  np.shape(first_order_topN_PATHWAY[Nth_PATHWAY])[1]) \
                                     for i in xrange(np.shape(first_order_topN_PATHWAY[Nth_PATHWAY])[0])]
        Variance_1st_all_PATHWAY[Nth_PATHWAY] = map(lambda x: x.data, SI_1st_object_all_PATHWAY)

        SI_2nd_object_all_PATHWAY = [vcs.SI_2nd_c(1.0, second_order_topN_PATHWAY[Nth_PATHWAY][i]) \
                                     for i in xrange(np.shape(second_order_topN_PATHWAY[Nth_PATHWAY])[0])]
        Variance_2nd_all_PATHWAY[Nth_PATHWAY] = map(lambda x: x.data, SI_2nd_object_all_PATHWAY)

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename1_1 = os.path.join(curr_dir, species_type + "_PATHWAY_variance_1st" + ".csv")
    print "PATHWAY_variance_1st filename:\t", filename1_1
    np.savetxt(filename1_1, Variance_1st_all_PATHWAY, delimiter=",", fmt="%.18e")
    filename2_2 = os.path.join(curr_dir, species_type + "_PATHWAY_variance_2nd" + ".csv")
    print "PATHWAY_variance_2nd filename:\t", filename2_2
    np.savetxt(filename2_2, Variance_2nd_all_PATHWAY, delimiter=",", fmt="%.18e")

    Variance_1st_all_PATHWAY_t = np.loadtxt(filename1_1, delimiter=",", dtype=float)
    Variance_2nd_all_PATHWAY_t = np.loadtxt(filename2_2, delimiter=",", dtype=float)
    SI_1st_all_PATHWAY_SELF_t = Variance_1st_all_PATHWAY_t
    SI_2nd_all_PATHWAY_SELF_t = Variance_2nd_all_PATHWAY_t

    for i in xrange(np.shape(SI_1st_all_PATHWAY_SELF_t)[0]):
        SI_1st_all_PATHWAY_SELF_t[i] = map(lambda x: x / var_list_PATHWAY[i], Variance_1st_all_PATHWAY_t[i])
        SI_2nd_all_PATHWAY_SELF_t[i] = map(lambda x: x / var_list_PATHWAY[i], Variance_2nd_all_PATHWAY_t[i])

    for Nth_PATHWAY in xrange(topN):
        # for Nth_PATHWAY in xrange(0, 0+1):
        ################################################################################################################
        # Plot
        ################################################################################################################
        """
        * $1^{st}$ order sensitivity index plot
        """
        curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI", "pathway_" + str(Nth_PATHWAY))
        if not os.path.exists(curr_dir):
            os.makedirs(curr_dir)
        filename1 = os.path.join(curr_dir, "SI_1st_PATHWAY_" + species_type + ".png")
        splot.bar_1D_PATHWAY_SI(SI_1st_all_PATHWAY_SELF_t[Nth_PATHWAY], PATHWAY_name[Nth_PATHWAY],
                                PATHWAY_prob[Nth_PATHWAY], \
                                filename1)

        """
        * $2^{nd}$ order sensitivty index plot
        """
        filename2 = os.path.join(curr_dir, "SI_2nd_PATHWAY_" + species_type + ".png")
        splot.bar_2D_PATHWAY_SI(SI_2nd_all_PATHWAY_SELF_t[Nth_PATHWAY], N_variable, N_variable, filename2, \
                                text_in="PATHWAY: #" + str(Nth_PATHWAY) + ", SELF")

        """
        * $2^{nd}$ order sensitivity index plot-bar3d
        """
        filename3 = os.path.join(curr_dir, "SI_2nd_PATHWAY_" + species_type + "_3D.png")
        splot.bar3D_2D_PATHWAY_SI(SI_2nd_all_PATHWAY_SELF_t[Nth_PATHWAY], N_variable, N_variable, filename3, \
                                  text_in="PATHWAY: #" + str(Nth_PATHWAY) + ", SELF")

    # * 2D plot, plot all $1^{st}$ PATHWAYs SELF SI in one figure, Version 1

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_SELF_PATHWAY_SI_1st_all_v1" + ".png")
    print "SELF_PATHWAY_SI_1st_all_v1 filename:\t", filename
    splot.bar3D_2D_PATHWAY_SI_v2(SI_1st_all_PATHWAY_SELF_t.T, filename, text_in="PATHWAY: #" + "ALL" + ", SELF")

    # * D plot, plot all $1^{st}$ PATHWAYs SELF SI in one figure, Version 2

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_SELF_PATHWAY_SI_1st_all_v2" + ".png")
    print "SELF_PATHWAY_SI_1st_all_v2 filename:\t", filename
    splot.bar_1D_PATHWAY_SI_v2(SI_1st_all_PATHWAY_SELF_t, filename, set_ybound_1=False)

    # * 2D plot, plot all $2^{nd}$ PATHWAYs SELF SI in one figure

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_SELF_PATHWAY_SI_2nd_all" + ".png")

    splot.bar3D_2D_PATHWAY_SI_v3(SI_2nd_all_PATHWAY_SELF_t, N_variable, N_variable, filename, \
                                 text_in="PATHWAY: #" + "ALL" + ", SELF")

    # * LSODE-Global sensitivity, SI=$\frac{var(pathway)}{var(LSODE~conc)}$

    topN = np.shape(var_coef_list_PATHWAY)[0];
    Variance_1st_all_PATHWAY_t = np.loadtxt(filename1_1, delimiter=",", dtype=float)
    Variance_2nd_all_PATHWAY_t = np.loadtxt(filename2_2, delimiter=",", dtype=float)
    SI_1st_all_PATHWAY_LSODE_t = Variance_1st_all_PATHWAY_t
    SI_2nd_all_PATHWAY_LSODE_t = Variance_2nd_all_PATHWAY_t

    for i in xrange(np.shape(SI_1st_all_PATHWAY_LSODE_t)[0]):
        SI_1st_all_PATHWAY_LSODE_t[i] = map(lambda x: x / var_LSODE, Variance_1st_all_PATHWAY_t[i])
        SI_2nd_all_PATHWAY_LSODE_t[i] = map(lambda x: x / var_LSODE, Variance_2nd_all_PATHWAY_t[i])

    filename_3 = os.path.join(curr_dir, species_type + "_PATHWAY_individual_SI_1st" + ".csv")
    np.savetxt(filename_3, SI_1st_all_PATHWAY_LSODE_t, delimiter=",", fmt="%.18e")
    filename_4 = os.path.join(curr_dir, species_type + "_PATHWAY_individual_SI_2nd" + ".csv")
    np.savetxt(filename_4, SI_2nd_all_PATHWAY_LSODE_t, delimiter=",", fmt="%.18e")

    for Nth_PATHWAY in xrange(topN):
        # for Nth_PATHWAY in xrange(0, 0+1):
        ################################################################################################################
        # Plot
        ################################################################################################################
        """
        * $1^{st}$ order sensitivity index plot
        """
        curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI", "pathway_" + str(Nth_PATHWAY))
        if not os.path.exists(curr_dir):
            os.makedirs(curr_dir)
        filename1 = os.path.join(curr_dir, "SI_1st_PATHWAY_LSODE_" + species_type + ".png")
        splot.bar_1D_PATHWAY_SI_v3(SI_1st_all_PATHWAY_LSODE_t[Nth_PATHWAY], PATHWAY_name[Nth_PATHWAY],
                                   PATHWAY_prob[Nth_PATHWAY], \
                                   filename1, set_ybound_1=False)

        """
        * $2^{nd}$ order sensitivty index plot
        """
        filename2 = os.path.join(curr_dir, "SI_2nd_PATHWAY_LSODE_" + species_type + ".png")
        splot.bar_2D_PATHWAY_SI(SI_2nd_all_PATHWAY_LSODE_t[Nth_PATHWAY], N_variable, N_variable, filename2, \
                                text_in="PATHWAY: #" + str(Nth_PATHWAY) + ", LSODE")

        """
        * $2^{nd}$ order sensitivity index plot-bar3d
        """
        filename3 = os.path.join(curr_dir, "SI_2nd_PATHWAY_LSODE_" + species_type + "_3D.png")
        splot.bar3D_2D_PATHWAY_SI(SI_2nd_all_PATHWAY_LSODE_t[Nth_PATHWAY], N_variable, N_variable, filename3, \
                                  text_in="PATHWAY: #" + str(Nth_PATHWAY) + ", LSODE")

    # * 2D plot, plot all $1^{st}$ PATHWAYs LSODE SI in one figure

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_LSODE_PATHWAY_SI_1st_all_v1" + ".png")
    print "LSODE_PATHWAY_SI_1st_all_v1 filename:\t", filename
    splot.bar3D_2D_PATHWAY_SI_v2(SI_1st_all_PATHWAY_LSODE_t.T, filename, text_in="PATHWAY: #" + "ALL" + ", LSODE")

    # * 1D plot, plot all $1^{st}$ PATHWAYs LSODE SI in one figure

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_LSODE_PATHWAY_SI_1st_all_v2" + ".png")
    print "LSODE_PATHWAY_SI_1st_all_v2 filename:\t", filename
    splot.bar_1D_PATHWAY_SI_v2(SI_1st_all_PATHWAY_LSODE_t, filename, set_ybound_1=False)

    # * 2D plot, plot all $2^{nd}$ PATHWAYs LSODE SI in one figure

    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename = os.path.join(curr_dir, species_type + "_LSODE_PATHWAY_SI_2nd_all" + ".png")

    splot.bar3D_2D_PATHWAY_SI_v3(SI_2nd_all_PATHWAY_LSODE_t, N_variable, N_variable, filename, \
                                 text_in="PATHWAY: #" + "ALL" + ", LSODE")
