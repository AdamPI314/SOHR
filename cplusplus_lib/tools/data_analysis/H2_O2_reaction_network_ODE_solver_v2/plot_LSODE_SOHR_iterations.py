""" Plot iterations. """
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import numpy as np
from matplotlib.lines import Line2D
import json


def plot_SOHR_LSODE_iteration(speIndex):
    """
        polt SOHR LSODE iterations results
    """
    file_dir = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
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

    # time
    time = np.loadtxt(os.path.join(file_dir, "output",
                                   "time_dlsode_fraction.csv"), delimiter=",")
    # concentration
    concentration = np.loadtxt(os.path.join(
        file_dir, "output", "concentration_dlsode_fraction.csv"), delimiter=",")

    # time from SOHR solution
    time_SOHR = np.loadtxt(os.path.join(
        file_dir, "output", "time_SOHR_fraction_all.csv"), delimiter=",")
    # anchor time points
    timeN1 = inputData["SOHR_init"]["timeN1"]
    tn = len(time_SOHR) / timeN1
    anchor_time_point_index = [(i + 1) * timeN1 for i in range(tn)]

    # SOHR concentration over iterations
    concentration_SOHR = [None] * iterationNumber
    for itr in range(0, iterationNumber):
        concentration_SOHR[itr] = np.loadtxt(os.path.join(file_dir, "output", "concentration_SOHR_fraction_all_" + str(itr + 1) + ".csv"),
                                             delimiter=",", dtype=float)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    step1 = 10
    step2 = 1
    step3 = 1
    # [0-->O2, 1-->H2O, 2-->H2, 3-->H2O2, 4-->H, 5-->OH, 6-->HO2, 7-->O]
    speName = ["O2", "H2O", "H2", "H2O2", "H", "OH", "HO2", "O"]
    # speIndex = 0

    # colors = ('r', 'k', 'b', 'c', 'g', 'y', 'm')
    Ncolors = iterationNumber / step3 + 2
    cmap = plt.get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, Ncolors)]

    # Initial guess and EXACT
    lengendHandlers = [None] * (iterationNumber + 2)
    initialGuessIndex = 0
    iterationCounter = 1
    for itr in range(0, iterationNumber, step3):
        startIndex = 0
        for endIndex in anchor_time_point_index:
            if itr == 0:
                iHandler, = ax.plot([time_SOHR[startIndex], time_SOHR[endIndex]], [concentration_SOHR[itr][startIndex, speIndex], concentration_SOHR[itr][startIndex, speIndex]] ,
                    color=colors[initialGuessIndex % len(colors)],
                    linewidth=1.5,
                    label ='SOHR initial guess')

            pHandler, = ax.plot(time_SOHR[startIndex:endIndex+1:step2], concentration_SOHR[itr][startIndex:endIndex+1:step2, speIndex], '--',
                color=colors[iterationCounter % len(colors)],
                linewidth=1.5,
                label ='SOHR' + ' n=' + str(itr + 1))

            # add legends, only for the first interval
            if endIndex == anchor_time_point_index[0]:
                lengendHandlers[itr+1] = pHandler
                if itr == 0:
                    lengendHandlers[0] = iHandler


            startIndex = endIndex + 1
        iterationCounter += 1

    exactHandler, = ax.plot(time[::step1], concentration[::step1, speIndex], '*', color=colors[iterationCounter % len(colors)], linewidth=1.5,
            markersize=5, label='EXACT')
    lengendHandlers[-1] = exactHandler

    ax.set_xbound(0, time_SOHR[-1])
    ax.grid('on')

    leg = ax.legend(handles = lengendHandlers, loc=0, prop={'size': 10.0})
    leg.get_frame().set_alpha(0.5)

    plt.xlabel("time", size=30)
    plt.ylabel("concentration", size=25)
    plt.title("Concentration of "+ speName[speIndex])

    # plt.show()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    fig.savefig(os.path.join(file_dir, "output",
                             "LSODE_SOHR_" + speName[speIndex] + "_iterations.jpg"), dpi=500)
    plt.close()


if __name__ == "__main__":
    for i in range(8):
        plot_SOHR_LSODE_iteration(i)
    print "Job Finished."
