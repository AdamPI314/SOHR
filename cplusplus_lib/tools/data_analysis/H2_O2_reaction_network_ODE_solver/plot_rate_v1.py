import os
import sys

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import numpy as np
from matplotlib.lines import Line2D
import argparse
import json

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", help="which iteration",
                           type=int, default=0, required=False)
    # which species
    argparser.add_argument("-r", help="which reaction",
                           type=int, default=0, required=False)

    args = argparser.parse_args()

    file_dir = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))

    print(file_dir)

    iteration = 0
    if args.i == 0:
        with open(os.path.join(file_dir, "input", "setting.json"), "r") as settingFileHandler:
            inputData = json.load(settingFileHandler)
        iteration = inputData["SOHR_init"]["iterationNumber"]
    else:
        iteration = args.i

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
    linestyles = Line2D.lineStyles.keys()

    # dlsode time
    time_dlsode = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_fraction.csv"), delimiter=",")
    # srode time
    time_srode = np.loadtxt(os.path.join(
        file_dir, "output", "time_srode_fraction_" + str(iteration) + ".csv"), delimiter=",")
    concentration_dlsode = np.loadtxt(os.path.join(
        file_dir, "output", "reaction_rate_dlsode_fraction.csv"), delimiter=",", dtype=float)
    concentration_srode = np.loadtxt(
        os.path.join(file_dir, "output", "reaction_rate_srode_fraction_" + str(iteration) + ".csv"),
        delimiter=",",
        dtype=float)

    labels = ["R"+str(i) for i in range(np.shape(concentration_dlsode)[1])]

    s_str = ""
    for i in range(len(labels)):
        s_str += "(" + str(i) + "," + labels[i] + ")"
    print(s_str)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)
    ax.plot(time_dlsode, concentration_dlsode[:, args.r], '--', color=colors[args.r % len(colors)], linewidth=2.5,
            label='dlsode ' + labels[args.r])
    ax.plot(time_srode, concentration_srode[:, args.r], linestyle='-', color=colors[args.r % len(colors)],
            linewidth=2.5,
            label='srode ' + labels[args.r])

    ax.legend(loc=0, prop={'size': 15.0})
    plt.xlabel("time", size=30)
    plt.ylabel("reaction rate", size=25)
    plt.title("ITERATION #" + str(iteration), size=25)
    ax.grid()
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))

    # plt.show()

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_rotation(45)

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", "rate_dlsode_srode_" + str(iteration) + "_" +
                             labels[args.r] + ".jpg"),
                dpi=500)
    plt.close()
