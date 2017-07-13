import os
import numpy as np
# from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! for plot 3D
import mpl_toolkits.mplot3d.axes3d  # register 3d projection
import matplotlib.pylab as plt


########################################################################################################################
class plot_1D_c:
    """
    plot 1D
    save to file
    """
    xmin = -1
    xmax = 1

    @staticmethod
    def func_leg(x, *coef):
        # exclude the 0th order coef, which is a constant
        coef_t = np.insert(coef, 0, 0, axis=0)
        return np.polynomial.legendre.legval(x, coef_t)

    def __init__(self, data_x, target_time, coef, zero_order, index, file_dir, prefix="",
                 file_name="target_vs_K_1D.jpg"):
        file_name1 = prefix + file_name.split(".")[0] + "_" + str(index) + "." + file_name.split(".")[-1]
        fig, axes = plt.subplots(1, 1, sharex=True, sharey=True)
        x = data_x;
        y = target_time
        axes.plot(x, y, 'o', lw=2, color='#FF7F00')
        axes.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        axes.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        x2 = np.linspace(self.xmin, self.xmax, 50);
        y2 = self.func_leg(x2, *coef)
        axes.plot(x2, y2)
        y3 = y2 + zero_order
        axes.plot(x2, y3)
        axes.set_xlabel("K" + str(index))
        axes.set_ylabel("target")

        fig.savefig(os.path.join(file_dir, file_name1), dpi=600)


########################################################################################################################
class plot_2D_c:
    """
    plot 2D
    save to file
    """

    # The definition below is only for 2D legendre value evaluation
    def my_legendre_polynomial_2D_val(self, data_x_y, *coef2D_in):
        # input parameter- pointer got to be pointer to 1D array
        # reshape to a 2D array
        coef2D = np.reshape(coef2D_in, (self.Nth_order, self.Nth_order))
        # exclude 0th and 1st order coef
        coef2D = np.insert(np.insert(coef2D, 0, 0, axis=0), 0, 0, axis=1)
        # value
        return np.polynomial.legendre.legval2d(data_x_y[0], data_x_y[1], coef2D)

    def __init__(self, data_x_y, target_time, my_coef2D, index, file_dir, prefix="", file_name="target_vs_K_2D.jpg"):
        self.Nth_order = int(np.sqrt(len(my_coef2D)))
        x = np.linspace(-1, 1)
        y = np.linspace(-1, 1)
        xg, yg = np.meshgrid(x, y)

        fig, ax = plt.subplots()
        #         plot1= ax.imshow(self.my_legendre_polynomial_2D_val([xg, yg], *my_coef2D), extent=[-1,1,1,-1], origin="upper")
        plot1 = ax.imshow(self.my_legendre_polynomial_2D_val([xg, yg], *my_coef2D), extent=[-1, 1, 1, -1],
                          origin="upper")

        cb = fig.colorbar(plot1, ax=ax)
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2, 2))
        cb.update_ticks()

        ax.set_xlabel("K" + str(index[1]))
        ax.set_ylabel("K" + str(index[0]))
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

        file_name1 = prefix + file_name.split(".")[0] + "_" + str(index[0]) + "_" + str(index[1]) + "_contour." + \
                     file_name.split(".")[-1]
        fig.savefig(os.path.join(file_dir, file_name1), dpi=600)

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        plt.hold(True)

        ax2.scatter(data_x_y[:, 0], data_x_y[:, 1], target_time, c='#d3ffce', marker='o', s=1.0, alpha=0.8)

        ax2.plot_surface(xg, yg, self.my_legendre_polynomial_2D_val([xg, yg], *my_coef2D), rstride=1, cstride=1,
                         cmap=plt.cm.jet, linewidth=0, antialiased=False, alpha=0.8)
        ax2.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax2.ticklabel_format(axis='z', style='sci', scilimits=(-2, 2))

        ax2.set_xlabel("K" + str(index[0]))
        ax2.set_ylabel("K" + str(index[1]))
        ax2.set_zlabel("target")
        ax2.set_title('Fit function')
        ax2.view_init(37.5, -30)

        # ax2.view_init(20,70)
        # ax2.view_init(100,-30)

        fig2.subplots_adjust(left=0.01, right=0.90, top=0.95)
        file_name2 = prefix + file_name.split(".")[0] + "_" + str(index[0]) + "_" + str(index[1]) + "." + \
                     file_name.split(".")[-1]
        fig2.savefig(os.path.join(file_dir, file_name2), dpi=600)


########################################################################################################################
class plot_1D_2D_c:
    """
    plot 1D and 2D together
    save to file
    """
    xmin = -1
    xmax = 1

    # 1D
    def func_leg_1D(self, data_x_y, *zero_first_in):
        # exclude the 0th order coef, which is a constant
        # zero_t= zero_first_in[0]
        coef0 = zero_first_in[1:1 + self.Nth_order_1st];
        coef1 = zero_first_in[1 + self.Nth_order_1st:]
        coef_t0 = np.insert(coef0, 0, 0, axis=0);
        coef_t1 = np.insert(coef1, 0, 0, axis=0)
        return np.polynomial.legendre.legval(data_x_y[0], coef_t0) + np.polynomial.legendre.legval(data_x_y[1], coef_t1)
        # return np.polynomial.legendre.legval(data_x_y[0], coef_t0)+np.polynomial.legendre.legval(data_x_y[1], coef_t1)+zero_t

    # The definition below is only for 2D legendre value evaluation
    def my_legendre_polynomial_2D_val(self, data_x_y, *coef2D_in):
        # input parameter- pointer got to be pointer to 1D array
        # reshape to a 2D array
        coef2D = np.reshape(coef2D_in, (self.Nth_order, self.Nth_order))
        # exclude 0th and 1st order coef
        coef2D = np.insert(np.insert(coef2D, 0, 0, axis=0), 0, 0, axis=1)
        # value
        return np.polynomial.legendre.legval2d(data_x_y[0], data_x_y[1], coef2D)

    # The definition below is only for 1D and 2D legendre value evaluation
    def my_legendre_polynomial_1D_2D_val(self, data_x_y, *zero_first_second_in):
        # input parameter- pointer got to be pointer to 1D array
        zero_t = zero_first_second_in[0]
        coef0 = zero_first_second_in[1:1 + self.Nth_order_1st]
        coef1 = zero_first_second_in[1 + self.Nth_order_1st:self.Nth_order_1st + 1 + self.Nth_order_1st]
        coef_t0 = np.insert(coef0, 0, 0, axis=0);
        coef_t1 = np.insert(coef1, 0, 0, axis=0)

        coef2D_in = zero_first_second_in[self.Nth_order_1st + 1 + self.Nth_order_1st:]
        # reshape to a 2D array
        coef2D = np.reshape(coef2D_in, (self.Nth_order_2nd, self.Nth_order_2nd))
        # exclude 0th and 1st order coef
        coef2D = np.insert(np.insert(coef2D, 0, 0, axis=0), 0, 0, axis=1)
        # value
        return np.polynomial.legendre.legval(data_x_y[0], coef_t0) + np.polynomial.legendre.legval(data_x_y[1],
                                                                                                   coef_t1) + zero_t + \
               np.polynomial.legendre.legval2d(data_x_y[0], data_x_y[1], coef2D)

    def __init__(self, data_x_y, target_time, zero_order, my_coef1D_1, my_coef1D_2, my_coef2D, \
                 index, file_dir, prefix="", file_name="target_vs_K_1D_2D.jpg"):
        zero_first = [zero_order] + list(my_coef1D_1) + list(my_coef1D_2)
        #         zero_first_second= np.concatenate([[zero_order], my_coef1D_1, my_coef1D_2, my_coef2D])
        zero_first_second = [zero_order] + list(my_coef1D_1) + list(my_coef1D_2) + list(my_coef2D)

        self.Nth_order_1st = len(my_coef1D_1)
        self.Nth_order_2nd = int(np.sqrt(len(my_coef2D)))

        self.Nth_order = int(np.sqrt(len(my_coef2D)))
        x = np.linspace(self.xmin, self.xmax)
        y = np.linspace(self.xmin, self.xmax)
        xg, yg = np.meshgrid(x, y)

        fig, ax = plt.subplots()
        #         plot1= ax.imshow(self.my_legendre_polynomial_2D_val([xg, yg], *my_coef2D), extent=[-1,1,1,-1], origin="upper")
        plot1 = ax.imshow(self.func_leg_1D([xg, yg], *zero_first), extent=[-1, 1, 1, -1], origin="upper")
        #         plot1= ax.imshow(self.my_legendre_polynomial_1D_2D_val([xg, yg], *zero_first_second), extent=[-1,1,1,-1], origin="upper")

        cb = fig.colorbar(plot1, ax=ax)
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2, 2))
        cb.update_ticks()

        ax.set_xlabel("K" + str(index[0]))
        ax.set_ylabel("K" + str(index[1]))
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

        file_name1 = prefix + file_name.split(".")[0] + "_" + str(index[0]) + "_" + str(index[1]) + "_contour." + \
                     file_name.split(".")[-1]
        fig.savefig(os.path.join(file_dir, file_name1), dpi=600)

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')

        ax2.scatter(data_x_y[:, 0], data_x_y[:, 1], target_time, c='#d3ffce', marker='o', s=1.0, alpha=0.8)

        ax2.plot_surface(xg, yg, self.my_legendre_polynomial_2D_val([xg, yg], *my_coef2D), rstride=1, cstride=1,
                         cmap=plt.cm.jet, linewidth=0, antialiased=False, alpha=0.8)
        ax2.plot_surface(xg, yg, self.func_leg_1D([xg, yg], *zero_first), rstride=1, cstride=1,
                         cmap=plt.cm.jet, linewidth=0, alpha=0.8)
        ax2.plot_surface(xg, yg, self.my_legendre_polynomial_1D_2D_val([xg, yg], *zero_first_second), rstride=1,
                         cstride=1,
                         cmap=plt.cm.jet, linewidth=0, alpha=0.8)

        ax2.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax2.ticklabel_format(axis='z', style='sci', scilimits=(-2, 2))

        ax2.set_xlabel("K" + str(index[0]))
        ax2.set_ylabel("K" + str(index[1]))
        ax2.set_zlabel("target")
        ax2.set_title('Fit function')
        #         ax2.view_init(37.5,-30)
        ax2.view_init(10.5, -15)

        #         ax2.view_init(20,70)
        #         ax2.view_init(100,-30)

        fig2.subplots_adjust(left=0.01, right=0.90, top=0.95)
        file_name2 = prefix + file_name.split(".")[0] + "_" + str(index[0]) + "_" + str(index[1]) + "." + \
                     file_name.split(".")[-1]
        fig2.savefig(os.path.join(file_dir, file_name2), dpi=600)
