import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import numpy as np
from matplotlib.lines import Line2D
import json

if __name__ == '__main__':
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

    # time
    time = np.loadtxt(os.path.join(file_dir, "output", "time_dlsode_M.csv"), delimiter=",")
    print np.shape(time)
    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration_dlsode_M.csv"), delimiter=",")
    print np.shape(concentration)

    # time from SOHR solution
    time_SOHR = np.loadtxt(os.path.join(file_dir, "output", "time_SOHR_M_all.csv"), delimiter=",")
    print np.shape(time_SOHR)
    # concentration from SOHR solution
    concentration_SOHR = np.loadtxt(os.path.join(file_dir, "output", "concentration_SOHR_M_all.csv"), delimiter=",",
                                     dtype=float)
    print np.shape(concentration_SOHR)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)


    N1 = 1 
    N2 = 1
    counter = 0
    # ax.plot(time[::N1], concentration[::N1, 0], 'o', color=colors[counter % len(colors)], linewidth=1.5, label='EXACT A')
    # counter += 1
    ax.plot(time[::N1], concentration[::N1, 1], '+', color=colors[counter % len(colors)], linewidth=1.5, label='EXACT X')
    counter += 1
    ax.plot(time[::N1], concentration[::N1, 2], 'o', color=colors[counter % len(colors)], linewidth=1.5, label='EXACT Y')
    counter += 1
    # ax.plot(time[::N1], concentration[::N1, 3], 'o', color=colors[counter % len(colors)], linewidth=1.5, label='EXACT B')
    # counter += 1

    counter = 0
    # ax.plot(time_SOHR[::N2], concentration_SOHR[::N2, 0], '--', color=colors[counter % len(colors)], linewidth=1.5, label='SOHR A')
    # counter += 1
    ax.plot(time_SOHR[::N2], concentration_SOHR[::N2, 1], '--', color=colors[counter % len(colors)], linewidth=1.5, label='SOHR X')
    counter += 1
    ax.plot(time_SOHR[::N2], concentration_SOHR[::N2, 2], '--', color=colors[counter % len(colors)], linewidth=1.5, label='SOHR Y')
    counter += 1
    # ax.plot(time_SOHR[::N2], concentration_SOHR[::N2, 3], '--', color=colors[counter % len(colors)], linewidth=1.5, label='SOHR B')
    # counter += 1

    # diff = concentration_SOHR[::N2, 1]+concentration_SOHR[::N2, 2]+concentration_SOHR[::N2, 3]-concentration[::N1, 1]-concentration[::N1, 2]-concentration[::N1, 3]
    # print diff

    # ax.plot(time_SOHR[::N2], diff, '--', color=colors[counter % len(colors)], linewidth=1.5, label='DIFFERENCE')
    # counter += 1

    ax.set_xbound(0, time_SOHR[-1])
    # ax.set_xbound(0, 20)
    # ax.set_ybound(-0.5, 1.05 * np.max(concentration_SOHR[::N2, -1]))
    ax.grid('on')

    ax.legend(loc=0, prop={'size': 10.0})
    plt.xlabel("time", size=30)
    plt.ylabel("concentration", size=25)
    # plt.title("N=" + str(inputData["SOHR_init"]["iterationNumber"]), size=25)
    plt.title("Concentrations of intermediates X and Y")

    # plt.show()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    fig.savefig(os.path.join(file_dir, "output", "LSODE_SOHR_X_Y.jpg"), dpi=500)
    plt.close()
