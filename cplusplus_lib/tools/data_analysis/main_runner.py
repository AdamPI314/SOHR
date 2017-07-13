#!/usr/bin/env python
import argparse
import os
import sys
import time

import lsode_pathway_conc_least_square_regression as lpclsr
import lsode_conc_sensitivity_index as lcsi
import pathway_probability_least_square_regression as pplsr
# import lsode_pathway_conc_converge as lpcc
import pathway_individual_sensitivity_index as pisi
import pathway_covariance_correlation as pcc
import pathway_covariance_sensitivity_index as pcsi

import pathway_probability as pp
import parse_spe_reaction_info as psri


def get_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', help="value of argument")

    return parser.parse_args(argv)


if __name__ == "__main__":
    start_time = time.time()

    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # '''
    #     1. initialization
    # '''
    # # tau
    # tau = 0.8
    # species_type = "H2O"
    # # LSODE coef filename
    # filename_conc = species_type + "_lsode_" + str(tau) + "_conc_coef.csv"
    # # topN pathway
    # topN = 50
    # # total number of pathway from simulation
    # Npathway = 50
    #
    # '''
    #     2. least square regression of concentration data, save to file
    # '''
    # filename_spe_conc = os.path.join(file_dir, "output/spe_conc_coef", filename_conc)
    # print filename_spe_conc
    # if not os.path.isfile(filename_spe_conc):
    #     lpclsr.least_square_regression(file_dir, species_type, filename_conc, tau)
    #     print "least square regression of concentration data, save to file. DONE!"
    # else:
    #     print filename_spe_conc, "exists already!"
    #
    # '''
    #     3. calculate and plot LSODE concentration sensitivity index
    # '''
    # lcsi.calculate_plot_SI(file_dir, filename_conc, species_type)
    # print "calculate and plot LSODE concentration sensitivity index. DONE!"
    #
    # '''
    #     4. least square regression of pathway probability data, rescale pathway
    #         probability first, save to file
    # '''
    # # pathway prob coef filename
    # filename_pp = species_type + "_coef.csv"
    # print "least square regression of pathway probability data. DONE!"
    # filename_spe_pp = os.path.join(file_dir, "output/pathway_prob_coef", filename_pp)
    # print filename_spe_pp
    # if not os.path.isfile(filename_spe_pp):
    #     pplsr.least_square_regression(file_dir, species_type, filename_pp, top_N=topN, N_pathway=Npathway)
    #     print "least square regression of pathway probability data, save to file. DONE!"
    # else:
    #     print filename_spe_pp, "exists already!"
    #
    # elapsed_time = time.time() - start_time
    # print "elapsed time:\t", elapsed_time / 60, " minutes"
    #
    # '''
    #     5. calculate pathway individual sensitivity index
    # '''
    # pisi.calculate_plot_pathway_individual_SI(file_dir, species_type, filename_conc, Npathway)
    # print "calculate pathway individual sensitivity index. DONE!"
    #
    # '''
    #     6. calculate pathway covariance and correlation
    # '''
    # pcc.calculate_plot_covariance_correlation(file_dir, species_type, filename_conc)
    # print "calculate pathway covariance and correlation. DONE!"
    #
    # '''
    #     7. calculate pathway covariance sensitivity index and plot
    # '''
    # pcsi.calculate_plot_pathway_covariance_SI(file_dir, species_type, filename_conc)
    # print "calculate pathway covariance sensitivity index and plot. DONE!"

    '''
        8. pathway filter
    '''
    spe_ind_name_dict, spe_name_ind_dict = psri.parse_spe_info(
        os.path.join(file_dir, "output", "species_labelling.csv"))

    pp.pathway_prob_c.read_pathname_s2f_path_endswith(os.path.join(file_dir, "output", "pathway_stat_backup.csv"),
                                                      os.path.join(file_dir, "output", "pathname_endswith_OH.csv"),
                                                      spe_name_ind_dict['OH'])

    psri.read_pathname_convert_2_real_spe_rxn_v2(os.path.join(file_dir, "output", "species_labelling.csv"),
                                                 os.path.join(file_dir, "output", "reaction_labelling.csv"),
                                                 os.path.join(file_dir, "output", "pathname_endswith_OH.csv"),
                                                 os.path.join(file_dir, "output", "real_pathname_endswith_OH.csv"),
                                                 10)
    elapsed_time = time.time() - start_time
    print "elapsed time:\t", elapsed_time / 60, " minutes"
