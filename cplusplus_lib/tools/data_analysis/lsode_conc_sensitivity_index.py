# -*- coding: utf-8 -*-

import numpy as np
import os

import my_utility as mu
import reaction as rxn

import fit_least_square_regression as flsr
import variance_correlation_SI as vcs
import parse_regression_coef as prc

import sensitivity_plot as splot


# # import conc coef and var from lsode calculation

def calculate_plot_SI(file_dir, conc_coef_filename, species_type="H2O"):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))

    # * LSODE   

    file_name_LSODE = os.path.join(file_dir, "output", "spe_conc_coef", conc_coef_filename)
    var_coef_frame_object_LSODE = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_LSODE)
    var_coef_list_LSODE = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object_LSODE)

    var_LSODE, coef_LSODE = prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_LSODE, 0)

    print "LSODE var:\t", var_LSODE

    # # Calculate sensitivity index from lsode conc
    # * $1_{st}$ order
    # * $2_{nd}$ order

    # N_variable=23; Nth_order_1st=5; Nth_order_2nd=2
    # N_variable=19; Nth_order_1st=4; Nth_order_2nd=2
    fit_info = np.loadtxt(os.path.join(file_dir, "tools/data_analysis", "fit_config.inp"), dtype=int);
    # print fit_info
    N_variable = fit_info[0];
    Nth_order_1st = fit_info[1];
    Nth_order_2nd = fit_info[2]

    zero_order_LSODE, first_order_LSODE, second_order_LSODE = \
        flsr.fit_1D_2D_all_c.split_1D_coef_array_static(coef_LSODE, N_variable, Nth_order_1st, Nth_order_2nd)

    SI_1st_object_all_LSODE = [vcs.SI_1st_c(var_LSODE, first_order_LSODE[i], np.shape(first_order_LSODE)[1]) \
                               for i in xrange(np.shape(first_order_LSODE)[0])]
    SI_1st_all_LSODE = map(lambda x: x.data, SI_1st_object_all_LSODE)
    print "sum of 1st order SI:\t", sum(SI_1st_all_LSODE)
    SI_1st_all_topN_ele_LSODE, SI_1st_all_topN_ind_LSODE = \
        mu.my_utility_c.topN_element_and_index(SI_1st_all_LSODE, N_variable)

    ####################################################################################################################
    SI_2nd_object_all_LSODE = [vcs.SI_2nd_c(var_LSODE, second_order_LSODE[i]) \
                               for i in xrange(np.shape(second_order_LSODE)[0])]
    SI_2nd_all_LSODE = map(lambda x: x.data, SI_2nd_object_all_LSODE)
    print "sum of 2nd order SI:\t", np.sum(SI_2nd_all_LSODE)

    topN_ith_pair_LSODE = 5
    SI_2nd_all_topN_ele_LSODE, SI_2nd_all_topN_ind_LSODE = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_LSODE, topN_ith_pair_LSODE)

    print "topN 1st order SI index:\t", map(lambda x: rxn.reaction_c.mapped_rxn_ind[x], SI_1st_all_topN_ind_LSODE)

    rxn_pair_ind = SI_2nd_all_topN_ind_LSODE
    print "topN 2nd order SI index:\t", map(
        lambda x: map(lambda y: rxn.reaction_c.mapped_rxn_ind[y], rxn.reaction_c.index_pair[x]), rxn_pair_ind)
    print "topN 2nd order SI:\t", SI_2nd_all_topN_ele_LSODE[0:topN_ith_pair_LSODE]

    curr_dir = os.path.join(file_dir, "output", species_type + "_conc_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)

    """
    * $1^{st}$ order sensitivity index plot
    """
    filename1 = os.path.join(curr_dir, "SI_1st_conc_" + species_type + ".png")
    splot.bar_1D_LSODE_SI(SI_1st_all_LSODE, filename1)

    """
    * $2^{nd}$ order sensitivty index plot
    """
    filename2 = os.path.join(curr_dir, "SI_2nd_conc_" + species_type + ".png")
    splot.bar_2D_LSODE_SI(SI_2nd_all_LSODE, N_variable, N_variable, filename2)

    """
    * $2^{nd}$ order sensitivity index plot-bar3d
    """
    filename3 = os.path.join(curr_dir, "SI_2nd_conc_" + species_type + "_3D.png")
    splot.bar3D_2D_LSODE_SI(SI_2nd_all_LSODE, N_variable, N_variable, filename3)
