import numpy as np
import os


class conc_error_c:
    """
    read data, parse data
    conc from lsode
    conc from all related pathway
    error between them
    """
    # static variables
    name_list = ['H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2']

    @staticmethod
    def read_conc_error(file_dir, name="conc_error.csv"):
        conc_error_filename_in = os.path.join(file_dir, "output", name)
        conc_error = np.loadtxt(conc_error_filename_in)
        return conc_error

    @staticmethod
    def get_lsode_conc_ith_species(conc_error_all, ith, period=6):
        lsode_conc_ith_data = conc_error_all[:, 0]
        lsode_conc_ith_data = np.reshape(lsode_conc_ith_data, (-1, period))
        return lsode_conc_ith_data[:, ith]

    @staticmethod
    def get_pathway_conc_ith_species(conc_error_all, ith, period=6):
        pathway_conc_ith_data = conc_error_all[:, 1]
        pathway_conc_ith_data = np.reshape(pathway_conc_ith_data, (-1, period))
        return pathway_conc_ith_data[:, ith]

    @staticmethod
    def get_relative_error_ith_species(conc_error_all, ith, period=6):
        relative_error_ith_data = conc_error_all[:, 2]
        relative_error_ith_data = np.reshape(relative_error_ith_data, (-1, period))
        return relative_error_ith_data[:, ith]

    @staticmethod
    def read_spe_conc(file_dir, name="spe_conc_all.csv"):
        spe_conc_filename_in = os.path.join(file_dir, "output", name)
        spe_conc = np.loadtxt(spe_conc_filename_in)
        return spe_conc

    @staticmethod
    # have 8 species, exclude O2 and O
    def get_lsode_conc_ith_species_v2(spe_conc_all, ith, period=6):
        lsode_conc_ith_data = np.reshape(spe_conc_all, (-1, period + 2))
        lsode_conc_ith_data = lsode_conc_ith_data[:, 1:7]
        return lsode_conc_ith_data[:, ith]
