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
    iterationNumber = inputData["SOHR_init"]["iterationNumber"]

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                markers.append(m)
        except TypeError:
            pass
    styles = markers + [r'$\lambda$',
        r'$\bowtie$',
        r'$\circlearrowleft$',
        r'$\clubsuit$',
        r'$\checkmark$']

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
    # SOHR concentration over iterations
    concentration_SOHR = [None] * iterationNumber
    for itr in range(0, iterationNumber):
        concentration_SOHR[itr] = np.loadtxt(os.path.join(file_dir, "output", "concentration_SOHR_M_all_" + str(itr + 1) + ".csv"),
            delimiter=",", dtype=float)
    print np.shape(concentration_SOHR)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    N1 = 1
    N2 = 1
    N3 = 1 

    # colors = ('r', 'k', 'b', 'c', 'g', 'y', 'm')
    Ncolors = iterationNumber / N3 + 2 
    cmap = plt.get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, Ncolors)]

    counter = 0
    for itr in range(0, iterationNumber, N3):
        ax.plot(time_SOHR[::N2], concentration_SOHR[itr][::N2, 1], '--', color=colors[counter % len(colors)],
                linewidth=1.5,
                label='SOHR X' + ' n=' + str(itr + 1))
        counter += 1

    ax.plot(time[::N1], concentration[::N1, 1], '*', color=colors[counter % len(colors)], linewidth=1.5,
            markersize=7, label='EXACT X')

    ax.set_xbound(0, time_SOHR[-1])
    ax.grid('on')

    leg = ax.legend(loc=0, prop={'size': 10.0})
    leg.get_frame().set_alpha(0.5)

    plt.xlabel("time", size=30)
    plt.ylabel("concentration", size=25)
    plt.title("Concentration of intermediate X")

    # plt.show()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    fig.savefig(os.path.join(file_dir, "output", "LSODE_SOHR_X_iterations.jpg"), dpi=500)
    plt.close()
