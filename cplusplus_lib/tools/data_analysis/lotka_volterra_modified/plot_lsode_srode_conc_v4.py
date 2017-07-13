import os
import sys
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

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    # # import data
    # time
    time = np.loadtxt(os.path.join(file_dir, "output", "time_dlsode_M.csv"), delimiter=",")
    print np.shape(time)

    # # time from srode solution
    # time_srode = np.loadtxt(os.path.join(file_dir, "output", "time_srode_0.csv"), delimiter=",")
    # print np.shape(time_srode)

    # concentration
    concentration = np.loadtxt(os.path.join(file_dir, "output", "concentration_dlsode_M.csv"), delimiter=",")
    print np.shape(concentration)

    # # concentration from srode solution
    # concentration_srode = np.loadtxt(os.path.join(file_dir, "output", "concentration_srode_0.csv"), delimiter=",",
    #                                  dtype=float)
    # print np.shape(concentration_srode)

    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True)

    counter = 0
    # ax.plot(time, concentration[:, 0], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode A1');
    # counter += 1
    # ax.plot(time, concentration[:, 1], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode A2');
    # counter += 1
    # ax.plot(time, concentration[:, 2], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode A3');
    counter += 1
    ax.plot(time, concentration[:, 3], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode X1');
    counter += 1
    ax.plot(time, concentration[:, 4], '--', color=colors[counter % len(colors)], linewidth=1.5, label='lsode X2');
    counter += 1

    # counter = 0
    # ax.plot(time_srode, concentration_srode[:, 0], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
    #         label='srode A1');
    # counter += 1
    # ax.plot(time_srode, concentration_srode[:, 1], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
    #         label='srode A2');
    # counter += 1
    # ax.plot(time_srode, concentration_srode[:, 2], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
    #         label='srode A3');
    # counter += 1
    # ax.plot(time_srode, concentration_srode[:, 3], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
    #         label='srode X1');
    # counter += 1
    # ax.plot(time_srode, concentration_srode[:, 4], linestyle='-', color=colors[counter % len(colors)], linewidth=1.5,
    #         label='srode X2');
    # counter += 1
    #
    # ax.set_xlim([time_srode[0], time_srode[-1]])

    ax.legend(loc=0, prop={'size': 10.0})
    ax.set_xbound(0.0, 30.0)
    ax.grid('on')

    plt.xlabel("time", size=30)
    plt.ylabel("concentration", size=25)
    # plt.title("N="+str(args.n), size=25)
    # plt.title("N=" + str(inputData["SOHR_init"]["iterationNumber"]), size=25)
    plt.title("Intermediate X1 and X2", size=25)


    # plt.show()
    # ax.set_ybound(-1.0, 4.0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
        tick.label.set_rotation(45)

    fig.tight_layout()

    # fig.savefig(os.path.join(file_dir, "output", "lsode_srode_concentration_"+str(args.n)+".jpg"), dpi=500)
    fig.savefig(os.path.join(file_dir, "output",
                             "lsode_concentration.jpg"),
                dpi=500)
    plt.close()
