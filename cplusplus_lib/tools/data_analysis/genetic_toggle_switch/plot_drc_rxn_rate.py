import os
# import sys
import matplotlib.pylab as plt
import numpy as np
from matplotlib.lines import Line2D

if __name__ == '__main__':
    # file_dir= "/home/shirong/Documents/summer_2015_dalian/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.12"
    file_dir = "D:\\fall_2015_boulder\SRC_EXE_DIR\RxnNetwork_TDDM_alpha_19.0\RxnNetwork_TDDM_alpha_19.0"
    print file_dir

    # # import data
    # * import data

    # time from SOHR solution
    time_SOHR = np.loadtxt(os.path.join(file_dir, "output", "time_SOHR.csv"))
    print np.shape(time_SOHR)

    # rxn_rate_0
    rxn_rate_0 = np.loadtxt(os.path.join(file_dir, "output", "reaction_rate_SOHR_time_0.csv"))
    print np.shape(rxn_rate_0)

    # drc 0 from SOHR solution
    drc_SOHR_0 = np.loadtxt(os.path.join(file_dir, "output", "drc_SOHR_time_0.csv"))
    print np.shape(drc_SOHR_0)

    # int_drc 0 from SOHR solution
    int_drc_SOHR_0 = np.loadtxt(os.path.join(file_dir, "output", "int_drc_SOHR_time_0.csv"))
    print np.shape(int_drc_SOHR_0)

    # # rxn_rate_1
    # rxn_rate_1 = np.loadtxt(os.path.join(file_dir, "output", "reaction_rate_SOHR_time_1.csv"))
    # print np.shape(rxn_rate_1)

    # # drc 1 from SOHR solution
    # drc_SOHR_1 = np.loadtxt(os.path.join(file_dir, "output", "drc_SOHR_time_1.csv"))
    # print np.shape(drc_SOHR_1)
    #
    # # int_drc 1 from SOHR solution
    # int_drc_SOHR_1 = np.loadtxt(os.path.join(file_dir, "output", "int_drc_SOHR_time_1.csv"))
    # print np.shape(int_drc_SOHR_1)

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                markers.append(m)
        except TypeError:
            pass
    styles = markers + [
        r'$\lambda$',
        r'$\bowtie$',
        r'$\circlearrowleft$',
        r'$\clubsuit$',
        r'$\checkmark$']

    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    linestyles = ['_', '-', '--', ':']

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    counter = 0

    # ax.plot(time_SOHR, rxn_rate_0[:,0], '--', color='r', linewidth= 1.0, label= 'rxn_rate_0 1')
    # ax.plot(time_SOHR, rxn_rate_0[:,1], '--', color='b', linewidth= 1.0, label= 'rxn_rate_0 2')
    # ax.plot(time_SOHR, rxn_rate_0[:,2], '--', color='c', linewidth= 1.0, label= 'rxn_rate_0 3')

    ax.plot(time_SOHR, drc_SOHR_0[:,0], linestyles[counter], markersize=10, color=colors[counter], label= 'drc A 0'); counter+=1
    ax.plot(time_SOHR, drc_SOHR_0[:,1], linestyles[counter], markersize=10, color=colors[counter], label= 'drc X 0'); counter+=1
    ax.plot(time_SOHR, drc_SOHR_0[:,2], linestyles[counter], markersize=10, color=colors[counter], label= 'drc Y 0'); counter+=1
    ax.plot(time_SOHR, drc_SOHR_0[:,3], linestyles[counter], markersize=10, color=colors[counter], label= 'drc B 0'); counter+=1

    # ax.plot(time_SOHR, int_drc_SOHR_0[:, 0], '--', color='r', linewidth=1.0, label='int drc A 0')
    # ax.plot(time_SOHR, int_drc_SOHR_0[:, 1], '--', color='b', linewidth=1.0, label='int drc X 0')
    # ax.plot(time_SOHR, int_drc_SOHR_0[:, 2], '--', color='y', linewidth=1.0, label='int drc Y 0')
    # ax.plot(time_SOHR, int_drc_SOHR_0[:, 3], '--', color='c', linewidth=1.0, label='int drc B 0')

    # ax.plot(time_SOHR, rxn_rate_1[:,0], '-', color='r', linewidth= 1.0, label= 'rxn_rate_1 1')
    # ax.plot(time_SOHR, rxn_rate_1[:,1], '-', color='b', linewidth= 1.0, label= 'rxn_rate_1 2')
    # ax.plot(time_SOHR, rxn_rate_1[:,2], '-', color='y', linewidth= 1.0, label= 'rxn_rate_1 3')

    # ax.plot(time_SOHR, drc_SOHR_1[:,0], '-', color='r', linewidth= 1.0, label= 'drc A 1')
    # ax.plot(time_SOHR, drc_SOHR_1[:,1], '-', color='b', linewidth= 1.0, label= 'drc X 1')
    # ax.plot(time_SOHR, drc_SOHR_1[:,2], '-', color='y', linewidth= 1.0, label= 'drc Y 1')
    # ax.plot(time_SOHR, drc_SOHR_1[:,3], '-', color='c', linewidth= 1.0, label= 'drc B 1')
    #
    # ax.plot(time_SOHR, int_drc_SOHR_1[:, 0], '-', color='r', linewidth=1.0, label='int drc A 1')
    # ax.plot(time_SOHR, int_drc_SOHR_1[:, 1], '-', color='b', linewidth=1.0, label='int drc X 1')
    # ax.plot(time_SOHR, int_drc_SOHR_1[:, 2], '-', color='y', linewidth=1.0, label='int drc Y 1')
    # ax.plot(time_SOHR, int_drc_SOHR_1[:, 3], '-', color='c', linewidth=1.0, label='int drc B 1')

    ax.set_ybound(0, 2.5)
    ax.legend(loc=2, prop={'size': 8.0})
    # plt.show()

    fig.savefig(os.path.join(file_dir, "output", "SOHR_drc_rxn_rate.jpg"), dpi=500)
    plt.close()
