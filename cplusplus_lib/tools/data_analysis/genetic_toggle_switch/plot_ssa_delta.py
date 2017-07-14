import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import numpy as np
from matplotlib.lines import Line2D
import json
from matplotlib import ticker

def plot_ssa(file_dir, _fileCounter):
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    print(file_dir)

    with open(os.path.join(file_dir, "input", "setting.json"), "r") as settingFileHandler:
        inputData = json.load(settingFileHandler)
    print("iterationNumber", inputData["SOHR_init"]["iterationNumber"])

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

    # linestyles = ['_', '-', '--', ':']
    linestyles = Line2D.lineStyles.keys()

    # time
    time = np.loadtxt(os.path.join(file_dir, "output", "time_ssa_number.csv"), delimiter=",")
    print(np.shape(time))

    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration_ssa_number.csv"), delimiter=",")
    print(np.shape(concentration))

    # colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    Nspecies = np.shape(concentration)[1]
    cmap = plt.get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, Nspecies)]

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    step1 = 1

    data_A = concentration[::step1, 0] + concentration[::step1, 2] + 2*concentration[::step1, 4] + 2*concentration[::step1, 6]
    data_B = concentration[::step1, 0] + concentration[::step1, 3] + 2*concentration[::step1, 5] + 2*concentration[::step1, 7]
    data = data_A - data_B


    N1 = -10
    counter = 0
    ax.plot(time[N1:]-time[N1], data[N1:], '-', color=colors[counter % len(colors)], linewidth=1.0, label='$N_A-N_B$')

    ax.set_xbound(0, time[-1]-time[N1])
    ax.grid('on')

    ax.legend(loc=0, prop={'size': 10.0})
    plt.xlabel("time", size=30)
    plt.ylabel("$\Delta$", size=25)

    # ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))

    # plt.show()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        # tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        # tick.label.set_rotation(45)

    fig.tight_layout()

    fig.savefig(os.path.join(file_dir, "output", "ssa_delta_" + str(_fileCounter+1) + ".jpg"), dpi=500)
    plt.close()
    print("Job Finished.")
