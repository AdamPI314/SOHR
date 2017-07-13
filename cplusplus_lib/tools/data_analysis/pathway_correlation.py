from scipy.integrate import quad
from scipy.integrate import dblquad
import numpy as np


########################################################################################################################
class correlation_1st_c:
    """
    calculate 1st order correlation of pathway, maybe
    only the terms wrt. the same rate constant K will contribute
    """
    xmin = -1
    xmax = 1

    @staticmethod
    def func_leg(x, coef):
        coef_t = np.insert(coef, 0, 0, axis=0)
        return np.polynomial.legendre.legval(x, coef_t)

    def density_uniform(self, x, a, b):
        return 1.0 / (b - a)

    def func_1st(self, x, coef1, coef2):
        return self.func_leg(x, coef1) * self.func_leg(x, coef2) * self.density_uniform(x, self.xmin, self.xmax)

    def correlation_1st(self, coef1, coef2, TotalVariance=1.0):
        I0 = quad(self.func_1st, self.xmin, self.xmax, args=(coef1, coef2))
        return I0[0] / TotalVariance

    def __init__(self, coef1, coef2, TotalVariance=1.0):
        self.data = self.correlation_1st(coef1, coef2, TotalVariance)


########################################################################################################################
class correlation_2nd_c:
    """
    calculate 2nd order correlation of pathway, maybe
    only the terms wrt. the same rate constant pair (i,j) will contribute
    """
    # static variables
    Nth_order = 2
    xmin = -1
    xmax = 1

    # evaluate at grid point
    def my_legendre_polynomial_2D_val_x_y(self, x, y, *coef1D_in):
        # input parameter- pointer got to be pointer to 1D array
        # reshape to a 1D array and a 2D array
        coef2D = np.reshape(coef1D_in, (self.Nth_order, self.Nth_order))
        # exclude 0th and 1st order coef
        coef2D = np.insert(np.insert(coef2D, 0, 0, axis=0), 0, 0, axis=1)
        # value
        return np.polynomial.legendre.legval2d(x, y, coef2D)

    def density_uniform(self, x, a, b):
        return 1.0 / (b - a)

    def func_2nd(self, x, y, coef1, coef2):
        return self.my_legendre_polynomial_2D_val_x_y(x, y, coef1) * self.my_legendre_polynomial_2D_val_x_y(x, y,
                                                                                                            coef2) * \
               self.density_uniform(x, self.xmin, self.xmax) * self.density_uniform(y, self.xmin, self.xmax)

    def SI_2nd(self, xmin, xmax, coef1, coef2, TotalVariance=1.0):
        I0 = dblquad(self.func_2nd, xmin, xmax, lambda x: xmin, lambda x: xmax, args=(coef1, coef2))
        return I0[0] / TotalVariance

    def __init__(self, my_coef2D1, my_coef2D2, TotalVariance=1.0):
        self.Nth_order = int(np.sqrt(len(my_coef2D1)))
        self.data = self.SI_2nd(self.xmin, self.xmax, my_coef2D1, my_coef2D2, TotalVariance)
