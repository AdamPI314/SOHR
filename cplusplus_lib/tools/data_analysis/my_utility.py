import os
import numpy as np


########################################################################################################################
class my_utility_c:
    # return top N element and their index

    # static variables
    # marker_all= ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p',\
    #          '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
    marker_all = ['.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', '|',
                  '_', 'D', 'd']

    @staticmethod
    def topN_element_and_index(in_list, topN=3):
        #     by using dict, will delete duplicated items
        #     dict_t= {in_list[i]:i for i in xrange(len(in_list))}
        #     dict_t2= sorted(dict_t.items(), key=lambda x: x[0], reverse=True)
        dict_t = [(in_list[i], i) for i in xrange(len(in_list))]
        dict_t2 = sorted(dict_t, reverse=True)
        ele = [x[0] for x in dict_t2[0:topN]];
        ind = [x[1] for x in dict_t2[0:topN]]
        return ele, ind


########################################################################################################################
class uncertainty_c:
    """
    read in uncertainty and normalize them in the range of (-1, 1)
    """

    @staticmethod
    def read_uncertainty(file_dir):
        ura_filename_in = os.path.join(file_dir, "output/uncertainties_random_all.csv")
        ura_filename_in_const = os.path.join(file_dir, "input/uncertainties.inp")
        uncertainties = np.loadtxt(ura_filename_in)
        uncertainties_const = np.loadtxt(ura_filename_in_const, delimiter=" ")[:, -1]
        # modification due to duplicated reaction
        # uncertainties_const[24]=uncertainties_const[23]

        # uncertainties_new=np.reshape(uncertainties, (-1, 25))
        uncertainties_new = np.reshape(uncertainties, (-1, 21))

        # data modification, logendre polynomial is orthonormal in the range of (-1, 1)
        for i in xrange(len(uncertainties_const)):
            uncertainties_new[:, i] = (uncertainties_new[:, i] - 1 / uncertainties_const[i]) \
                                      / (uncertainties_const[i] - 1 / uncertainties_const[i])
            uncertainties_new[:, i] = uncertainties_new[:, i] * 2 - 1
        return uncertainties_new

    """
    read uncertainty const, reshape, re-normalized
    """

    @staticmethod
    def read_nominal_Ks(file_dir):
        ura_filename_in_const = os.path.join(file_dir, "input/uncertainties.inp")
        uncertainties_const = np.loadtxt(ura_filename_in_const, delimiter=" ")[:, -1]
        # modification due to duplicated reaction
        # uncertainties_const[24]=uncertainties_const[23]

        # index_2nd_all= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23]
        index_2nd_all = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19]

        uncertainties_denominator = uncertainties_const - np.divide(1.0, uncertainties_const)
        uncertainties_numerator = 1.0 - np.divide(1.0, uncertainties_const)
        uncertainties_K0 = np.divide(uncertainties_numerator, uncertainties_denominator)
        uncertainties_K0 = uncertainties_K0 * 2 - 1

        uncertainties_K0 = uncertainties_K0[index_2nd_all]
        return uncertainties_K0

    def __init__(self, file_dir):
        self.data = self.read_uncertainty(file_dir)


class target_time_c:
    """
    read in target time
    """

    @staticmethod
    def read_target_time(file_dir):
        target_time_filename_in = os.path.join(file_dir, "output/target_time_all.csv")
        target_time = np.loadtxt(target_time_filename_in)
        #     target_time=target_time-np.mean(target_time)
        return target_time

    def __init__(self, file_dir):
        self.data = self.read_target_time(file_dir)
