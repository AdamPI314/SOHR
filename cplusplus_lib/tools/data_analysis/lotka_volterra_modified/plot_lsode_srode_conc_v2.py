import os
# import sys
import matplotlib.pylab as plt
import numpy as np
from matplotlib.lines import Line2D

if __name__ == '__main__':
    # file_dir= "/home/shirong/Documents/summer_2015_dalian/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.12"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.15"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.16"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.17"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.24"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.25"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_18.0"
    # file_dir= "/home/shirong/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_18.4"
    # file_dir = "D:\\fall_2015_boulder\SRC_EXE_DIR\RxnNetwork_TDDM_alpha_19.0\RxnNetwork_TDDM_alpha_19.0"
    file_dir = "/home/invictus/Documents/spring_2016_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_20.0.0/RxnNetwork_TDDM_alpha_19.0"
    print file_dir

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
    # linestyles = ['_', '-', '--', ':']
    linestyles = Line2D.lineStyles.keys()
    print linestyles
    print markers

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    # # import data
    # * import data

    # time
    time = np.loadtxt(os.path.join(file_dir, "output", "time.csv"), delimiter=",")
    print np.shape(time)

    # time from srode solution
    time_srode = np.loadtxt(os.path.join(file_dir, "output", "time_srode.csv"), delimiter=",")
    print np.shape(time_srode)

    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration.csv"), delimiter=",")
    print np.shape(concentration)

    # concentration from srode solution
    concentration_srode = np.loadtxt(os.path.join(file_dir, "output", "concentration_srode.csv"), delimiter=",", dtype=float)
    print np.shape(concentration_srode)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    ax.plot(time, concentration[:, 0], '--', color='r', linewidth=1.0, label='lsode A')
    ax.plot(time, concentration[:, 1], '--', color='b', linewidth=1.0, label='lsode X')
    ax.plot(time, concentration[:, 2], '--', color='g', linewidth=1.0, label='lsode Y')
    ax.plot(time, concentration[:, 3], '--', color='c', linewidth=1.0, label='lsode B')

    counter = 0
    ax.plot(time_srode, concentration_srode[:, 0], linestyle='-', color=colors[counter % len(colors)], label='drc A 0'); counter+=1
    ax.plot(time_srode, concentration_srode[:, 1], linestyle='-', color=colors[counter % len(colors)], label='drc X 0'); counter+=1
    ax.plot(time_srode, concentration_srode[:, 2], linestyle='-', color=colors[counter % len(colors)], label='drc Y 0'); counter+=1
    ax.plot(time_srode, concentration_srode[:, 3], linestyle='-', color=colors[counter % len(colors)], label='drc Z 0'); counter+=1

    ax.set_xlim([time_srode[0], time_srode[-1]])

    ax.legend(loc=0, prop={'size': 10.0})
    plt.xlabel("time", size=20)
    plt.ylabel("concentration", size=20)
    plt.title("N=1", size=20)
    # plt.show()
    ax.set_ybound(-1.0, 4.0)

    fig.savefig(os.path.join(file_dir, "output", "lsode_srode_concentration.jpg"), dpi=500)
    plt.close()
