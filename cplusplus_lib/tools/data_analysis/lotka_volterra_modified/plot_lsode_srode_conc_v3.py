import os
import sys
import matplotlib.pylab as plt
import numpy as np
from matplotlib.lines import Line2D
import argparse
import json

if __name__ == '__main__':
    # argparser = argparse.ArgumentParser()
    # argparser.add_argument("--n", help = "Number of iterations", type=int, default=1, required=False)

    # args= argparser.parse_args()

    # file_dir= "/home/shirong/Documents/summer_2015_dalian/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.12"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.15"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.16"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.17"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.24"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_17.25"
    # file_dir= "/home/invictus/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_18.0"
    # file_dir= "/home/shirong/Documents/fall_2015_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_18.4"
    # file_dir = "D:\\fall_2015_boulder\SRC_EXE_DIR\RxnNetwork_TDDM_alpha_19.0\RxnNetwork_TDDM_alpha_19.0"
    # file_dir = "/home/invictus/Documents/spring_2016_boulder/SRC_EXE_DIR/RxnNetwork_TDDM_alpha_20.0.0/RxnNetwork_TDDM_alpha_19.0"
    # file_dir = "D:\VS_workspace\RxnNetwork_TDDM_alpha_20.1.0\RxnNetwork_TDDM_alpha_19.0"
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    print file_dir

    with open(os.path.join(file_dir, "input", "setting.json"), "r") as settingFileHandler:
        inputData = json.load(settingFileHandler)
    print "iterationNumber", type(inputData["SOHR_init"]["iterationNumber"])

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
    time = np.loadtxt(os.path.join(file_dir, "output", "time_dlsode_M.csv"), delimiter=",")
    print np.shape(time)

    # time from srode solution
    time_srode = np.loadtxt(os.path.join(file_dir, "output", "time_srode_0.csv"), delimiter=",")
    print np.shape(time_srode)

    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration_dlsode_M.csv"), delimiter=",")
    print np.shape(concentration)

    # concentration from srode solution
    concentration_srode = np.loadtxt(os.path.join(file_dir, "output", "concentration_srode_0.csv"), delimiter=",",
                                     dtype=float)
    print np.shape(concentration_srode)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    counter = 0
    ax.plot(time, concentration[:, 0], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode A');
    counter += 1
    ax.plot(time, concentration[:, 1], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode X');
    counter += 1
    ax.plot(time, concentration[:, 2], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode Y');
    counter += 1
    ax.plot(time, concentration[:, 3], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode B');
    counter += 1

    counter = 0
    ax.plot(time_srode, concentration_srode[:, 0], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
            label='srode A');
    counter += 1
    ax.plot(time_srode, concentration_srode[:, 1], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
            label='srode X');
    counter += 1
    ax.plot(time_srode, concentration_srode[:, 2], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
            label='srode Y');
    counter += 1
    ax.plot(time_srode, concentration_srode[:, 3], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
            label='srode B');
    counter += 1

    ax.set_xlim([time_srode[0], time_srode[-1]])

    ax.legend(loc=0, prop={'size': 10.0})
    plt.xlabel("time", size=30)
    plt.ylabel("concentration", size=25)
    # plt.title("N="+str(args.n), size=25)
    plt.title("N=" + str(inputData["SOHR_init"]["iterationNumber"]), size=25)

    # plt.show()
    ax.set_ybound(-1.0, 4.0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    # fig.savefig(os.path.join(file_dir, "output", "lsode_srode_concentration_"+str(args.n)+".jpg"), dpi=500)
    fig.savefig(os.path.join(file_dir, "output",
                             "lsode_srode_concentration_" + str(inputData["SOHR_init"]["iterationNumber"]) + ".jpg"),
                dpi=500)
    plt.close()
