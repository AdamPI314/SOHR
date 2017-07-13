import numpy as np
# import re
import sys
# import scipy as sp
import time
import os
import parse_regression_coef as prc

if __name__ == "__main__":
    start_time = time.time()

    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # tau, species_type, topN pathways
    tau = 0.8;
    species_type = "H2O";
    topN = 50;
    # LSODE coef filename
    filename_conc = species_type + "_lsode_" + str(tau) + "_conc_coef.csv"

    # file_name_LSODE= os.path.join(file_dir, "output", "spe_conc_coef", "H2O_lsode_0.8_conc_coef.csv")
    file_name_LSODE = os.path.join(file_dir, "output", "spe_conc_coef", filename_conc)
    var_coef_frame_object_LSODE = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_LSODE)
    var_coef_list_LSODE = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(var_coef_frame_object_LSODE)

    var_LSODE, coef_LSODE = prc.parse_regression_coef_c.var_zero_first_second_in_a_list(var_coef_list_LSODE, 0)

    #################################################################################################
    ################### Pathway individual sensitivity index ########################################
    #################################################################################################
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_SI")
    if not os.path.exists(curr_dir):
        os.makedirs(curr_dir)
    filename1_1 = os.path.join(curr_dir, species_type + "_PATHWAY_variance_1st" + ".csv")
    print "PATHWAY_variance_1st filename:\t", filename1_1
    filename2_2 = os.path.join(curr_dir, species_type + "_PATHWAY_variance_2nd" + ".csv")
    print "PATHWAY_variance_2nd filename:\t", filename2_2

    file_name_PATHWAY = os.path.join(file_dir, "output", "pathway_prob_coef", "H2O_coef.csv")
    var_coef_frame_object_PATHWAY = prc.parse_regression_coef_c.read_pathway_as_pandas_frame_object(file_name_PATHWAY)
    var_coef_list_PATHWAY = prc.parse_regression_coef_c.convert_pandas_frame_object_to_list(
        var_coef_frame_object_PATHWAY)
    topN = np.shape(var_coef_list_PATHWAY)[0];

    Variance_1st_all_PATHWAY_t = np.loadtxt(filename1_1, delimiter=",", dtype=float)
    Variance_2nd_all_PATHWAY_t = np.loadtxt(filename2_2, delimiter=",", dtype=float)
    SI_1st_all_PATHWAY_LSODE_t = Variance_1st_all_PATHWAY_t
    SI_2nd_all_PATHWAY_LSODE_t = Variance_2nd_all_PATHWAY_t

    for i in xrange(np.shape(SI_1st_all_PATHWAY_LSODE_t)[0]):
        SI_1st_all_PATHWAY_LSODE_t[i] = map(lambda x: x / var_LSODE, Variance_1st_all_PATHWAY_t[i])
        SI_2nd_all_PATHWAY_LSODE_t[i] = map(lambda x: x / var_LSODE, Variance_2nd_all_PATHWAY_t[i])

    # print np.shape(SI_1st_all_PATHWAY_LSODE_t)
    # print SI_1st_all_PATHWAY_LSODE_t
    # print np.shape(SI_2nd_all_PATHWAY_LSODE_t)
    filename_3 = os.path.join(curr_dir, species_type + "_PATHWAY_individual_SI_1st" + ".csv")
    np.savetxt(filename_3, SI_1st_all_PATHWAY_LSODE_t, delimiter=",", fmt="%.18e")
    filename_4 = os.path.join(curr_dir, species_type + "_PATHWAY_individual_SI_2nd" + ".csv")
    np.savetxt(filename_4, SI_2nd_all_PATHWAY_LSODE_t, delimiter=",", fmt="%.18e")

    #################################################################################################
    ################### Pathway correlation sensitivity index #######################################
    #################################################################################################
    curr_dir = os.path.join(file_dir, "output", species_type + "_pathway_correlation_SI")
    filename1 = os.path.join(curr_dir, species_type + "_1st_covariance.csv")

    first_order_covariance_t = np.loadtxt(filename1, delimiter=",", dtype=float)
    first_order_covariance_SI_t = first_order_covariance_t / var_LSODE
    print "first order covariance SI data shape:\t", np.shape(first_order_covariance_SI_t)
    filename_5 = os.path.join(curr_dir, species_type + "_correlation_SI_1st" + ".csv")
    np.savetxt(filename_5, first_order_covariance_SI_t, delimiter=",", fmt="%.18e")

    filename2 = os.path.join(curr_dir, species_type + "_2nd_covariance.csv")
    second_order_covariance_t = np.loadtxt(filename2, delimiter=",", dtype=float)
    second_order_covariance_SI_t = second_order_covariance_t / var_LSODE
    print "second order covariance SI data shape:\t", np.shape(second_order_covariance_SI_t)
    filename_6 = os.path.join(curr_dir, species_type + "_correlation_SI_2nd" + ".csv")
    np.savetxt(filename_6, second_order_covariance_SI_t, delimiter=",", fmt="%.18e")

    elapsed_time = time.time() - start_time
    print elapsed_time, "sec"
