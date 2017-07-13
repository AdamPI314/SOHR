import numpy as np
import fit_least_square_regression as flsr

from operator import add


def total_coef_topN_pathway(coef_list, topN):
    coef_t = [0] * len(coef_list[0])
    for i in xrange(0, topN):
        coef_t = map(add, coef_list[i], coef_t)
    return coef_t


def evaluate_conc_at_K_all(data_all, coef_N, N_variable=23, Nth_order_1st=5, Nth_order_2nd=2):
    zero_order_t, first_order_t, second_order_t = flsr.fit_1D_2D_all_c.split_1D_coef_array_static(coef_N, \
                                                                                                  N_variable,
                                                                                                  Nth_order_1st,
                                                                                                  Nth_order_2nd)
    # return value
    # 0th_order
    value_t = zero_order_t
    # 1st_order
    for i in xrange(N_variable):
        # insert dummy value, exclude the constant terms
        coef_t = np.insert(first_order_t[i], 0, 0, axis=0)
        value_t += np.polynomial.legendre.legval(data_all[i], coef_t)
    # 2nd_order
    order_index = 0  # label the 23 matrice
    for i in xrange(0, N_variable):
        for j in xrange(i + 1, N_variable):
            coef_2nd_order_matrix = np.reshape(second_order_t[order_index], (Nth_order_2nd, Nth_order_2nd))
            coef_2nd_order_matrix = np.insert(np.insert(coef_2nd_order_matrix, 0, 0, axis=0), 0, 0, axis=1)
            value_t += np.polynomial.legendre.legval2d(data_all[i], data_all[j], coef_2nd_order_matrix)
            order_index += 1

    return value_t
