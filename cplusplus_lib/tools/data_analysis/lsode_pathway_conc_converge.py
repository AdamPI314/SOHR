# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import os

import my_utility as mu
import parse_regression_coef as prc
import evaluate_concentration as ec
import conc_error as ce

from matplotlib import cm


# # Import uncertainty data
# * $k_i=2\times\frac{k_i-k_i^{MIN}}{k_i^{MAX}-k_i^{MIN}}-1$
# * in the range of (-1, 1)

def plot_lsode_pathway_conc(file_dir, species_type="H2O"):
    # file_dir= os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    uncertainty = mu.uncertainty_c(file_dir)

    # # import conc coef and var from lsode calculation
    # * LSODE

    file_name_LSODE = os.path.join(file_dir, "output", "spe_conc_coef", "H2O_lsode_0.8_conc_coef.csv")
    var_coef_frame_object_LSODE = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_LSODE)
    var_coef_list_LSODE = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object_LSODE)

    var_LSODE, coef_all_LSODE = prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_LSODE, 0)

    # # Import least square regression coef
    # * zeroth order, first order and second order
    # * variance of pathway probability

    file_name = os.path.join(file_dir, "output", "coef", "H2O_coef.csv")
    var_coef_frame_object = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name)
    var_coef_list = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object)

    var_list = map(lambda x: x[0], [prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list, i) \
                                    for i in xrange(np.shape(var_coef_list)[0])])
    coef_list = map(lambda x: x[1], [prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list, i) \
                                     for i in xrange(np.shape(var_coef_list)[0])])

    print "pathway var:\t", var_list
    print "coef data shape:\t", np.shape(coef_list)

    # # Calculate relative concentration from pathway probability
    # * concentration

    N_variable = 23;
    # Nth_order_1st=5; Nth_order_2nd=2
    topN = np.shape(var_coef_list)[0];
    N_pathway = 10
    coef_N = ec.total_coef_topN_pathway(coef_list, N_pathway)
    print "sum of topN coef data, has shape:\t", np.shape(coef_N)

    # pre-factor, pre-coefficient
    # convert before
    # coef_N= [2.0/3.0* coef_t for coef_t in coef_N]

    # # Calculate concentration from LSODE directly
    # * LSODE

    uncertainties_K0 = mu.uncertainty_c.read_nominal_Ks(file_dir)
    LSODE_t = ec.evaluate_conc_at_K_all(uncertainties_K0, coef_all_LSODE, N_variable=23, Nth_order_1st=5,
                                        Nth_order_2nd=2)
    PATHWAY_t = ec.evaluate_conc_at_K_all(uncertainties_K0, coef_N, N_variable=23, Nth_order_1st=5, Nth_order_2nd=2)

    print "LSODE:\t", LSODE_t
    print "PATHWAY:\t", PATHWAY_t
    print "relative error:\t", (PATHWAY_t - LSODE_t) / LSODE_t

    data_all = [np.random.rand() * 2 - 1.0 for i in xrange(N_variable)]
    LSODE_t = ec.evaluate_conc_at_K_all(data_all, coef_all_LSODE, N_variable=23, Nth_order_1st=5, Nth_order_2nd=2)
    PATHWAY_t = ec.evaluate_conc_at_K_all(data_all, coef_N, N_variable=23, Nth_order_1st=5, Nth_order_2nd=2)

    print "LSODE:\t", LSODE_t
    print "PATHWAY:\t", PATHWAY_t
    print "relative error:\t", (PATHWAY_t - LSODE_t) / LSODE_t

    # # Import conc_error data 
    # * $C(k_1, k_2, ..., k_{n_k})=C_0+\sum_{i=1}^{n_k}{A_i(k_i)}+\sum_{i<j}^{n_k}{\sum_{j=1}^{n_k}{B_{ij}(k_i,k_j)}}$

    conc_error_all = ce.conc_error_c.read_conc_error(file_dir, "conc_error.csv")

    # # Target function, a specific species, at a specific time
    # * name_list=['H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2']
    #     * name_list=[0->'H2O', 1->'H2', 2->'H2O2', 3->'H', 4->'OH', 5->'HO2']
    # * definition of concentration
    #     * initial concentration H2->2.0/3, O2->1.0/3, H->4.0/3
    #     * relative concentration= $\frac{concentration}{total ~number ~of ~H ~atoms \times 3.0/4}$

    lsode_z = ce.conc_error_c.get_lsode_conc_ith_species(conc_error_all, 0, period=6)

    pathway_z = ce.conc_error_c.get_pathway_conc_ith_species(conc_error_all, 0, period=6)

    index_2nd_all = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23]
    uncertainty_data_all = uncertainty.data[:, index_2nd_all];

    LSODE_fit_z = [ec.evaluate_conc_at_K_all(uncertainty_data_all[i], coef_all_LSODE, N_variable=23, Nth_order_1st=5,
                                             Nth_order_2nd=2) \
                   for i in xrange(np.shape(uncertainty_data_all)[0])]

    topNset = 50  # how many data points
    fig, ax = plt.subplots()

    colors = cm.rainbow(np.linspace(0, 1, topN + 3))

    ax.plot(lsode_z[0:topNset], label="lsode", color=colors[-1])
    ax.plot(LSODE_fit_z[0:topNset], label="lsode regression", color=colors[-2])
    ax.plot(pathway_z[0:topNset], label="all pathways", color=colors[-3])

    for N_pathway in xrange(topN):
        coef_N = ec.total_coef_topN_pathway(coef_list, N_pathway)
        # pre-factor, pre-coefficient
        # convert before
        # coef_N= [2.0/3.0* coef_t for coef_t in coef_N]
        PATHWAY_fit_z = [
            ec.evaluate_conc_at_K_all(uncertainty_data_all[i], coef_N, N_variable=23, Nth_order_1st=5, Nth_order_2nd=2) \
            for i in xrange(topNset)]

        ax.plot(PATHWAY_fit_z[0:topNset], label="top " + str(N_pathway + 1) + " pathways", color=colors[N_pathway])
        ax.legend(loc=2, bbox_to_anchor=(1.0, 1.0), prop={'size': 5})

    ax.set_title("[$H_2O$] from LSODE and PATHWAY Analysis")
    ax.set_xlabel("Random Points in RateConstant HyperCube")
    ax.set_ylabel("[$H_2O$]")
    fig.subplots_adjust(right=0.85)
    curr_dir = os.path.join(file_dir, "output", "conc_converge")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    # fig.savefig(os.path.join(curr_dir, "H2O_top_"+str(topNset)+".jpg"), dpi=600, bbox_inches= "tight")
    fig.savefig(os.path.join(curr_dir, "H2O_top_" + str(topNset) + ".jpg"), dpi=600)
