import os
import numpy as np


class conc_speciation_c:
    """
    speciation of concentration, read data, parse data
    """

    @staticmethod
    def read_target_conc(file_dir, name="conc_all.csv"):
        target_conc_filename_in = os.path.join(file_dir, "output", name)
        target_conc = np.loadtxt(target_conc_filename_in)
        #     target_time=target_time-np.mean(target_time)
        return target_conc

    @staticmethod
    def get_conc_ith_species(target_conc_all, ith, period=37):
        ith_conc_data = target_conc_all[:, ith]
        ith_conc_data = np.reshape(ith_conc_data, (-1, period))
        return ith_conc_data

    @staticmethod
    def get_conc_ith_species_jth_time(target_conc_all, ith, jth, period=37):
        ith_conc_data = target_conc_all[:, ith]
        ith_conc_data = np.reshape(ith_conc_data, (-1, period))
        ith_jth_conc_data = ith_conc_data[:, jth]
        return ith_jth_conc_data
