import os
import sys
import matplotlib.pylab as plt
import numpy as np
from matplotlib.lines import Line2D
import json

if __name__ == '__main__':
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    print(file_dir)

    with open(os.path.join(file_dir, "input", "setting.json"), "r") as settingFileHandler:
        inputData = json.load(settingFileHandler)
    print("iterationNumber:\t", inputData["SOHR_init"]["iterationNumber"])

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

    # time
    time = np.loadtxt(os.path.join(file_dir, "output", "time_dlsode_M.csv"), delimiter=",")

    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration_dlsode_M.csv"), delimiter=",")

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    counter = 0
    step1 = 1
    ax.plot(time[::step1], concentration[::step1, 0], '-', color=colors[counter % len(colors)], linewidth=1.0, label='X0')
    counter += 1
    ax.plot(time[::step1], concentration[::step1, 1], '-', color=colors[counter % len(colors)], linewidth=1.0, label='X1')
    counter += 1
    ax.plot(time[::step1], concentration[::step1, 2], '-', color=colors[counter % len(colors)], linewidth=1.0, label='X2')
    counter += 1

    ax.set_xbound(0, time[-1])
    ax.grid('on')

    ax.legend(loc=0, prop={'size': 10.0})

    # plt.show()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    fig.savefig(os.path.join(file_dir, "output", "ssa_concentration.jpg"), dpi=500)
    plt.close()
    print("Job Finished.")
