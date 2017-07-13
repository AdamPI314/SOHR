import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import re

import my_utility as mu
import reaction as rxn

from mpl_toolkits.mplot3d import Axes3D

from matplotlib import ticker
from matplotlib import cm
import math

"""
sensitivity index plot, 1D, LSODE
"""


def bar_1D_LSODE_SI(SI_1st_all_LSODE, filename):
    ####################################################################################################################
    # figure object
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
    reaction_obj = rxn.reaction_c()

    N = len(SI_1st_all_LSODE)
    x = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars: can also be len(x) sequence

    topN_SI = N;
    SI_ele, SI_ind = mu.my_utility_c.topN_element_and_index(SI_1st_all_LSODE, topN_SI)
    colors = colors_t = matplotlib.cm.rainbow(np.linspace(0, 1, N))

    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, 10))  # top 10
    colors_t2 = colors_t2[::-1]

    for i in xrange(np.shape(colors_t)[0]):
        if i < np.shape(colors_t2)[0]:
            colors_t[i] = colors_t2[i]
        else:
            colors_t[i] = colors_t2[-1]

    for i in xrange(np.shape(colors)[0]):
        colors[SI_ind[i]] = colors_t[i]
    bar1 = ax.bar(x, SI_1st_all_LSODE, width, color=colors, bottom=0, align="center")

    ax.set_xticks(x)
    ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(N)]
    ax.set_xticklabels(ticks, rotation=45, fontsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax.margins(0.01)
    # ax.set_xbound(0, 23)
    ax.set_ybound(0, 1.0)

    ax.set_xlabel("rxn")
    ax.set_ylabel('$1^{st}$' + "SI")
    ax.set_title('$1^{st}$' + "SI" + " of LSODE conc wrt K's")
    ax.text(1, 1.0 - 0.1, "SUM($1^{st}~$SI): " + '%1.4g' % (np.sum(SI_1st_all_LSODE)))

    # text of reaction involved
    matched_rxn = SI_ind[0:5]
    rx = 1;
    ry = 1.0 - 0.2;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, ry - (i + 1) * step_t, '$R_{' + str(rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]) + "}$" + \
                ': ' + reaction_obj.rxns_latex[rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]], color='#1B9E77',
                alpha=0.8)

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%.3g' % (height),
                    ha='center', va='bottom', fontsize=4)

    autolabel(bar1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    #   os.makedirs(curr_dir)
    fig.savefig(filename, dpi=1000, bbox_inches='tight')
    plt.close(fig)


"""
sensitivity index plot, 2D, LSODE
"""


def bar_2D_LSODE_SI(SI_2nd_all_LSODE, N_variable1, N_variable2, filename2):
    """
    Problem:
    I spent last few days trying to find a way to remove tiny margins from axes in a 3D plot. \
    I tried ax.margins(0) and ax.autoscale_view('tight') and other approaches, but these small margins are still there. \
    In particular, I don't like that the bar histograms are elevated, i.e., \
    their bottom is not at the zero level -- see example image.
    Solution:
    There is not property or method that can modify this margins. You need to patch the source code.
    http://stackoverflow.com/questions/16488182/removing-axes-margins-in-3d-plot
    """
    ###patch start###
    from mpl_toolkits.mplot3d.axis3d import Axis
    if not hasattr(Axis, "_get_coord_info_old"):
        def _get_coord_info_new(self, renderer):
            mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
            mins += deltas / 4
            maxs -= deltas / 4
            return mins, maxs, centers, deltas, tc, highs

        Axis._get_coord_info_old = Axis._get_coord_info
        Axis._get_coord_info = _get_coord_info_new
    ###patch end###

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_LSODE)):
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][0]][rxn.reaction_c.index_pair[i][1]] = SI_2nd_all_LSODE[i]
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][1]][rxn.reaction_c.index_pair[i][0]] = SI_2nd_all_LSODE[i]

    # color Matrix element   
    topN_ith_pair_LSODE = 5
    SI_2nd_all_topN_ele_LSODE, SI_2nd_all_topN_ind_LSODE = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_LSODE, topN_ith_pair_LSODE)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_LSODE))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_LSODE)):
        pair_ind = SI_2nd_all_topN_ind_LSODE[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][0]][rxn.reaction_c.index_pair[pair_ind][1]] = colors_t2[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][1]][rxn.reaction_c.index_pair[pair_ind][0]] = colors_t2[i]

    ## plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # iterate over row
    for i in xrange(np.shape(Matrix_2nd_SI)[0]):
        xs = np.arange(np.shape(Matrix_2nd_SI)[0])
        # plot column
        ys = Matrix_2nd_SI[:, i]

        # You can provide either a single color or an array. To demonstrate this,
        cs = [c for c in Matrix_2nd_color[:, i]]
        # ax.bar(xs, ys, zs=i, zdir='y', color=cs, alpha=0.8, align= "center")
        ax.bar(xs, ys, zs=i, zdir='y', color=cs, alpha=0.8)

    ax.set_xlabel('R')
    ax.set_ylabel('R')
    ax.set_zlabel('$2^{nd}$ SI')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(-2, 2))

    x = y = np.arange(np.shape(Matrix_2nd_SI)[0])  # the x locations for the groups
    x = x + 0.5;
    ax.set_xticks(x)
    x_ticks = y_ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    ax.set_xticklabels(x_ticks, rotation=45, fontsize=5)

    y = y;
    ax.set_yticks(y)
    ax.set_yticklabels(y_ticks, rotation=-15, fontsize=5)
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_LSODE[0] * 0.8,
            "SUM($2^{nd}~$SI): " + '%1.4g' % (np.sum(Matrix_2nd_SI) / 2.0))

    ## To change the color of the grid lines
    # for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    #    axis._axinfo['grid']['color']  = 0.7, 1.0, 0.7, 1.0
    # # Set X & Y lim
    # ax.set_ylim3d(0, np.shape(Matrix_2nd_SI)[0])
    # ax.set_xlim3d(0, np.shape(Matrix_2nd_SI)[1])
    ax.margins(0.0)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename2= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+".png")
    fig.savefig(filename2, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, LSODE
"""


def bar3D_2D_LSODE_SI(SI_2nd_all_LSODE, N_variable1, N_variable2, filename3):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_LSODE)):
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][0]][rxn.reaction_c.index_pair[i][1]] = SI_2nd_all_LSODE[i]
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][1]][rxn.reaction_c.index_pair[i][0]] = SI_2nd_all_LSODE[i]

    # color Matrix element   
    topN_ith_pair_LSODE = 5
    SI_2nd_all_topN_ele_LSODE, SI_2nd_all_topN_ind_LSODE = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_LSODE, topN_ith_pair_LSODE)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_LSODE))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_LSODE)):
        pair_ind = SI_2nd_all_topN_ind_LSODE[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][0]][rxn.reaction_c.index_pair[pair_ind][1]] = colors_t2[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][1]][rxn.reaction_c.index_pair[pair_ind][0]] = colors_t2[i]

    column_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]
    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()

    colors = Matrix_2nd_color.flatten()

    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.5, edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.w_zaxis.set_major_formatter(formatter)
    # ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) # doesn't work


    ax.set_xlabel('R')
    ax.set_ylabel('R')
    ax.set_zlabel('$2^{nd}$ order SI')
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_LSODE[0] * 0.8,
            "SUM($2^{nd}~$SI): " + '%1.4g' % (np.sum(Matrix_2nd_SI) / 2.0))
    ax.margins(0.0)

    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    # plt.show()

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
converted pathway name to regex format
"""


def PATHWAY_name_to_regex(pathway_name):
    matched_RS = re.findall("[R|S]", pathway_name)
    matched_N = re.findall("\d+[*]?", pathway_name)
    str_t = "$"
    for i in xrange(len(matched_RS)):
        if "*" in matched_N[i]:
            str_t += matched_RS[i] + "_{" + matched_N[i][0:-1] + "}" + "^{" + matched_N[i][-1] + "}"
        else:
            str_t += matched_RS[i] + "_{" + matched_N[i] + "}"
    str_t += "$"
    return str_t


"""
sensitivity index plot, 1D; include pathway name in the figure
"""


def bar_1D_PATHWAY_SI(SI_1st_all_LSODE, pathway_name, pathway_prob, filename, set_ybound_1=True):
    #############################################################################################################################
    # figure object
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
    reaction_obj = rxn.reaction_c()

    N = len(SI_1st_all_LSODE)
    x = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars: can also be len(x) sequence

    topN_SI = N;
    SI_ele, SI_ind = mu.my_utility_c.topN_element_and_index(SI_1st_all_LSODE, topN_SI)
    colors = colors_t = matplotlib.cm.rainbow(np.linspace(0, 1, N))

    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, 10))  # top 10
    colors_t2 = colors_t2[::-1]

    for i in xrange(np.shape(colors_t)[0]):
        if i < np.shape(colors_t2)[0]:
            colors_t[i] = colors_t2[i]
        else:
            colors_t[i] = colors_t2[-1]

    for i in xrange(np.shape(colors)[0]):
        colors[SI_ind[i]] = colors_t[i]
    bar1 = ax.bar(x, SI_1st_all_LSODE, width, color=colors, bottom=0, align="center", alpha=0.8)

    ax.set_xticks(x)
    ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(N)]
    ax.set_xticklabels(ticks, rotation=45, fontsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax.margins(0.01)
    # ax.set_xbound(0, 23)

    ax.set_xlabel("rxn")
    ax.set_ylabel('$1^{st}$' + "SI")
    ax.set_title('$1^{st}$' + "SI" + " of pathway prob wrt K's")

    if set_ybound_1 == True:
        ax.set_ybound(0, 1.0)
    else:
        ax.set_ybound(0, 1.5 * max(SI_1st_all_LSODE))
    ybound = ax.get_ybound()[1]
    pathway_name_t = PATHWAY_name_to_regex(reaction_obj.re_label_index(pathway_name))
    ax.text(1, (1.0 - 0.08) * ybound, "PATHWAY:", color="b", alpha=0.8)
    ax.text(1, (1.0 - 0.14) * ybound, pathway_name_t, color="b", alpha=0.8)
    ax.text(1, (1.0 - 0.22) * ybound, "PATHWAY prob: " + '%1.4g' % (pathway_prob), alpha=0.8)
    ax.text(1, (1.0 - 0.28) * ybound, "SUM($1^{st}~$SI): " + '%1.4g' % (np.sum(SI_1st_all_LSODE)), alpha=0.8)

    # text of reaction involved
    matched_rxn = re.findall("(?<=R)\d+", reaction_obj.re_label_index(pathway_name))
    rx = 1;
    ry = 1.0 - 0.36;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, (ry - (i + 1) * step_t) * ybound, '$R_{' + str(matched_rxn[i]) + "}$" + \
                ': ' + reaction_obj.rxns_latex[int(matched_rxn[i])], color='#F15A22', alpha=0.8)

    # text of reaction with topN SI
    matched_rxn = SI_ind[0:5]
    rx = 12;
    ry = 1.0 - 0.36;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, (ry - (i + 1) * step_t) * ybound,
                '$R_{' + str(rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]) + "}$" + \
                ': ' + reaction_obj.rxns_latex[rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]], color='#1B9E77',
                alpha=0.8)

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%1.3g' % (height),
                    ha='center', va='bottom', fontsize=4)

    autolabel(bar1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    #   os.makedirs(curr_dir)
    fig.savefig(filename, dpi=1000, bbox_inches='tight')
    plt.close(fig)


"""
sensitivity index plot, 1D; include pathway name in the figure, v2, all pathway together, sum up
"""


# def bar_1D_PATHWAY_SI_v2(SI_1st_all_LSODE, pathway_name, pathway_prob, filename, set_ybound_1= True):
def bar_1D_PATHWAY_SI_v2(SI_1st_all_LSODE, filename, set_ybound_1=True):
    ####################################################################################################################
    # figure object
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
    reaction_obj = rxn.reaction_c()

    lx = np.shape(SI_1st_all_LSODE)[0]
    ly = np.shape(SI_1st_all_LSODE)[1]

    bottoms = np.zeros(ly)
    sum_height = np.zeros(ly)

    colors_t = matplotlib.cm.rainbow(np.linspace(0, 1, lx))  # top 10
    colors_t = colors_t[::-1]

    for i in xrange(lx):
        x = np.arange(ly)  # the x locations for the groups
        width = 0.35  # the width of the bars: can also be len(x) sequence

        topN_SI = ly;
        SI_ele, SI_ind = mu.my_utility_c.topN_element_and_index(SI_1st_all_LSODE[i], topN_SI)
        colors = [colors_t[i] for color_i in xrange(ly)]

        #         for i in xrange(np.shape(colors_t)[0]):
        #             if i < np.shape(colors_t2)[0]:
        #                 colors_t[i]= colors_t2[i]
        #             else:
        #                 colors_t[i]=colors_t2[-1]

        #         for i in xrange(np.shape(colors)[0]):
        #             colors[SI_ind[i]]= colors_t[i]

        sum_height = np.add(sum_height, SI_1st_all_LSODE[i])
        bar1 = ax.bar(x, SI_1st_all_LSODE[i], width, color=colors, bottom=bottoms, align="center", alpha=0.8)
        bottoms = np.add(bottoms, SI_1st_all_LSODE[i])

    ax.set_xticks(x)
    ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(ly)]
    ax.set_xticklabels(ticks, rotation=45, fontsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax.margins(0.01)
    # ax.set_xbound(0, 23)

    ax.set_xlabel("rxn")
    ax.set_ylabel('$1^{st}$' + "SI")
    ax.set_title('$1^{st}$' + "SI" + " of pathway prob wrt K's")

    if set_ybound_1 == True:
        ax.set_ybound(0, 1.0)
    else:
        ax.set_ybound(0, 1.5 * max(sum_height))
    ybound = ax.get_ybound()[1]
    #     pathway_name_t= PATHWAY_name_to_regex(reaction_obj.re_label_index(pathway_name))
    #     ax.text(1, (1.0-0.08)*ybound ,"PATHWAY:", color= "b", alpha= 0.8)
    #     ax.text(1, (1.0-0.14)*ybound ,pathway_name_t, color= "b", alpha= 0.8)
    #     ax.text(1, (1.0-0.22)*ybound ,"PATHWAY prob: "+'%1.4g'%(pathway_prob), alpha= 0.8)
    ax.text(1, (1.0 - 0.15) * ybound, "SUM($1^{st}~$SI): " + '%1.4g' % (np.sum(SI_1st_all_LSODE)), alpha=0.8)

    #     #text of reaction involved
    #     matched_rxn= re.findall("(?<=R)\d+", reaction_obj.re_label_index(pathway_name))
    #     rx= 1; ry= 1.0-0.36; step_t= 0.05
    #     for i in xrange(len(matched_rxn)):
    #         ax.text(rx, (ry-(i+1)*step_t)*ybound, '$R_{'+str(matched_rxn[i])+"}$"+\
    #                 ': '+reaction_obj.rxns_latex[int(matched_rxn[i])], color='#F15A22', alpha=0.8)

    # text of reaction with topN SI
    matched_rxn = SI_ind[0:5]
    rx = 12;
    ry = 1.0 - 0.1;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, (ry - (i + 1) * step_t) * ybound,
                '$R_{' + str(rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]) + "}$" + \
                ': ' + reaction_obj.rxns_latex[rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]], color='#1B9E77',
                alpha=0.8)

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height() + rect.get_y()
            ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%1.3g' % (height),
                    ha='center', va='bottom', fontsize=4)

    autolabel(bar1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    #   os.makedirs(curr_dir)
    fig.savefig(filename, dpi=1000, bbox_inches='tight')
    plt.close(fig)


"""
sensitivity index plot, 1D; include pathway name in the figure
"""


def bar_1D_PATHWAY_SI_v3(SI_1st_all_LSODE, pathway_name, pathway_prob, filename, set_ybound_1=True):
    ####################################################################################################################
    # figure object
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
    reaction_obj = rxn.reaction_c()

    N = len(SI_1st_all_LSODE)
    x = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars: can also be len(x) sequence

    topN_SI = N;
    SI_ele, SI_ind = mu.my_utility_c.topN_element_and_index(SI_1st_all_LSODE, topN_SI)
    colors = colors_t = matplotlib.cm.rainbow(np.linspace(0, 1, N))

    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, 10))  # top 10
    colors_t2 = colors_t2[::-1]

    for i in xrange(np.shape(colors_t)[0]):
        if i < np.shape(colors_t2)[0]:
            colors_t[i] = colors_t2[i]
        else:
            colors_t[i] = colors_t2[-1]

    for i in xrange(np.shape(colors)[0]):
        colors[SI_ind[i]] = colors_t[i]
    bar1 = ax.bar(x, SI_1st_all_LSODE, width, color=colors, bottom=0, align="center", alpha=0.8)

    ax.set_xticks(x)
    ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(N)]
    ax.set_xticklabels(ticks, rotation=45, fontsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax.margins(0.01)
    # ax.set_xbound(0, 23)

    ax.set_xlabel("reaction index")
    ax.set_ylabel('$1^{st}$' + "order SI")
    ax.set_title('$1^{st}$' + "order SI" + " of pathway prob wrt. K's")

    if set_ybound_1 == True:
        ax.set_ybound(0, 1.0)
    else:
        ax.set_ybound(0, 1.5 * max(SI_1st_all_LSODE))
    ybound = ax.get_ybound()[1]
    pathway_name_t = PATHWAY_name_to_regex(reaction_obj.re_label_index(pathway_name))
    ax.text(0.5, (1.0 - 0.08) * ybound, "PATHWAY:", color="b", alpha=0.8)
    ax.text(0.5, (1.0 - 0.14) * ybound, pathway_name_t, color="b", alpha=0.8)
    ax.text(0.5, (1.0 - 0.22) * ybound, "PATHWAY prob: " + '%1.4g' % (pathway_prob), alpha=0.8)
    ax.text(0.5, (1.0 - 0.28) * ybound, "SUM($1^{st}~$order SI): " + '%1.4g' % (np.sum(SI_1st_all_LSODE)), alpha=0.8)

    # text of reaction involved
    matched_rxn_with_star = re.findall("(?<=R)\d+\*?", reaction_obj.re_label_index(pathway_name))
    matched_rxn = re.findall("(?<=R)\d+", reaction_obj.re_label_index(pathway_name))

    rx = 0.5;
    ry = 1.0 - 0.36;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, (ry - (i + 1) * step_t) * ybound, '$R_{' + str(matched_rxn[i]) + "}" + "^{" \
                + str("" if len(matched_rxn_with_star[i]) == len(matched_rxn[i]) else "*") + "}$" \
                + ': ' + reaction_obj.rxns_latex[int(matched_rxn[i])], color='#F15A22', alpha=0.8)

    # text of reaction with topN SI
    matched_rxn = SI_ind[0:5]
    rx = 10.5;
    ry = 1.0 - 0.36;
    step_t = 0.05
    for i in xrange(len(matched_rxn)):
        ax.text(rx, (ry - (i + 1) * step_t) * ybound,
                '$R_{' + str(rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]) + "}$" + \
                ': ' + reaction_obj.rxns_latex[rxn.reaction_c.mapped_rxn_ind[matched_rxn[i]]], color='#1B9E77',
                alpha=0.8)

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%1.3g' % (height),
                    ha='center', va='bottom', fontsize=4)

    autolabel(bar1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    #     plt.show()
    fig.savefig(filename, dpi=1000, bbox_inches='tight')
    plt.close(fig)


"""
sensitivity index plot, 2D, redefine colors
"""


def bar_2D_PATHWAY_SI(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename2, text_in=""):
    """
    Problem:
    I spent last few days trying to find a way to remove tiny margins from axes in a 3D plot. \
    I tried ax.margins(0) and ax.autoscale_view('tight') and other approaches, but these small margins are still there. \
    In particular, I don't like that the bar histograms are elevated, i.e., \
    their bottom is not at the zero level -- see example image.
    Solution:
    There is not property or method that can modify this margins. You need to patch the source code.
    http://stackoverflow.com/questions/16488182/removing-axes-margins-in-3d-plot
    """
    ###patch start###
    from mpl_toolkits.mplot3d.axis3d import Axis
    if not hasattr(Axis, "_get_coord_info_old"):
        def _get_coord_info_new(self, renderer):
            mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
            mins += deltas / 4
            maxs -= deltas / 4
            return mins, maxs, centers, deltas, tc, highs

        Axis._get_coord_info_old = Axis._get_coord_info
        Axis._get_coord_info = _get_coord_info_new
    ###patch end###

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        # Matrix_2nd_SI[rxn.reaction_c.index_pair[i][0]][rxn.reaction_c.index_pair[i][1]]= SI_2nd_all_PATHWAY[i]
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][1]][rxn.reaction_c.index_pair[i][0]] = SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        pair_ind = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][0]][rxn.reaction_c.index_pair[pair_ind][1]] = colors_t2[i]
        # Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][0]][rxn.reaction_c.index_pair[pair_ind][1]]= [1.0, 1.0, 1.0, 0.0]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][1]][rxn.reaction_c.index_pair[pair_ind][0]] = colors_t2[i]

    ## plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # iterate over row
    for i in xrange(np.shape(Matrix_2nd_SI)[0]):
        xs = np.arange(np.shape(Matrix_2nd_SI)[0])
        # plot column
        ys = Matrix_2nd_SI[:, i]

        # You can provide either a single color or an array. To demonstrate this,
        cs = [c for c in Matrix_2nd_color[:, i]]
        # ax.bar(xs, ys, zs=i, zdir='y', color=cs, alpha=0.8, align= "center")
        ax.bar(xs, ys, zs=i, zdir='y', color=cs, alpha=0.8)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('$2^{nd}$ SI')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(-2, 2))

    x = y = np.arange(np.shape(Matrix_2nd_SI)[0])  # the x locations for the groups
    x = x + 0.5;
    ax.set_xticks(x)
    x_ticks = y_ticks = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    ax.set_xticklabels(x_ticks, rotation=45, fontsize=5)

    y = y;
    ax.set_yticks(y)
    ax.set_yticklabels(y_ticks, rotation=-15, fontsize=5)
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.0,
            "SUM($2^{nd}~$SI): " + '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.1, text_in)

    ## To change the color of the grid lines
    # for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    #    axis._axinfo['grid']['color']  = 0.7, 1.0, 0.7, 1.0
    # # Set X & Y lim
    # ax.set_ylim3d(0, np.shape(Matrix_2nd_SI)[0])
    # ax.set_xlim3d(0, np.shape(Matrix_2nd_SI)[1])
    ax.margins(0.0)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename2= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+".png")
    fig.savefig(filename2, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors
"""


def bar3D_2D_PATHWAY_SI(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[rxn.reaction_c.index_pair[i][0]][rxn.reaction_c.index_pair[i][1]] = SI_2nd_all_PATHWAY[i]
    #         Matrix_2nd_SI[rxn.reaction_c.index_pair[i][1]][rxn.reaction_c.index_pair[i][0]]= SI_2nd_all_PATHWAY[i]

    # color Matrix element   
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        pair_ind = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][0]][rxn.reaction_c.index_pair[pair_ind][1]] = colors_t2[i]
        Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][1]][rxn.reaction_c.index_pair[pair_ind][0]] = colors_t2[i]
        # Matrix_2nd_color[rxn.reaction_c.index_pair[pair_ind][1]][rxn.reaction_c.index_pair[pair_ind][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]
    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()

    colors = Matrix_2nd_color.flatten()

    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.5, edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.w_zaxis.set_major_formatter(formatter)
    # ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) # doesn't work


    ax.set_xlabel('R')
    ax.set_ylabel('R')
    ax.set_zlabel('$2^{nd}$ order SI')
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.0,
            "SUM($2^{nd}~$SI): " + '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.1, text_in)
    ax.margins(0.0)

    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    # plt.show()

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, plot all 1st order sensitivity index together
"""


def bar3D_2D_PATHWAY_SI_v2(Matrix_2nd_SI, filename3, text_in=""):
    N_variable1 = np.shape(Matrix_2nd_SI)[0]
    N_variable2 = np.shape(Matrix_2nd_SI)[1]

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))

    # color Matrix element   
    # topN_ith_pair_PATHWAY= 5

    SI_2nd_all_PATHWAY = Matrix_2nd_SI.T.flatten()
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([N_variable1, N_variable2], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        Matrix_2nd_color[SI_2nd_all_topN_ind_PATHWAY[i] / N_variable2][SI_2nd_all_topN_ind_PATHWAY[i] % N_variable2] = \
        colors_t2[i]

    column_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(N_variable1)]
    # PATHWAY_SUM= [sum(Matrix_2nd_SI.T[i]) for i in xrange(np.shape(Matrix_2nd_SI.T)[0])]
    # row_names = [str(i)+': '+'%1.1e'%(PATHWAY_SUM[i]) for i in xrange(N_variable2)]
    row_names = [str(i) for i in xrange(N_variable2)]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = N_variable1  # Work out matrix dimensions
    ly = N_variable2
    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.T.flatten()

    colors = Matrix_2nd_color.flatten()

    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.5, edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, edgecolor=[1.0, 1.0, 1.0, 0.0])
    #     ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.5)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.w_zaxis.set_major_formatter(formatter)
    # ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) # doesn't work


    ax.set_xlabel('R')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('$1^{st}$ order SI')
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.0,
            "SUM($1^{st}~$SI): " + '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, SI_2nd_all_topN_ele_PATHWAY[0] * 1.1, text_in)
    ax.margins(0.0)

    ticksx = np.arange(0.0, N_variable1, 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, N_variable2, 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    # plt.show()

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, plot 2nd order all together, sum topN pathway for each bin
"""


def bar3D_2D_PATHWAY_SI_v3(SI_2nd_all_PATHWAY, lx, ly, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))

    topN = len(SI_2nd_all_PATHWAY)

    fig = plt.figure()
    ax = Axes3D(fig)

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    #     zpos = np.zeros(lx*ly)
    dz = np.zeros(lx * ly)

    max_hist = 0.0
    sum_topN_pathway = 0.0

    for ith_P in xrange(topN):
        # SI Matrix element
        Matrix_2nd_SI = np.zeros(lx * ly).reshape((lx, ly))
        # symmetric matrix
        for i in xrange(len(SI_2nd_all_PATHWAY[ith_P])):
            Matrix_2nd_SI[rxn.reaction_c.index_pair[i][0]][rxn.reaction_c.index_pair[i][1]] = SI_2nd_all_PATHWAY[ith_P][
                i]
            Matrix_2nd_SI[rxn.reaction_c.index_pair[i][1]][rxn.reaction_c.index_pair[i][0]] = SI_2nd_all_PATHWAY[ith_P][
                i]

        # color Matrix element   
        topN_ith_pair_PATHWAY = 2  # just need the highest one
        # topN_ith_pair_PATHWAY= len(SI_2nd_all_PATHWAY[ith_P])
        SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
            mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY[ith_P], topN_ith_pair_PATHWAY)

        max_hist += SI_2nd_all_topN_ele_PATHWAY[0]
        sum_topN_pathway += np.sum(Matrix_2nd_SI)

        # color, just topN kinds of colors
        colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN))  # top 10
        colors_t2 = colors_t2[::-1]

        Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
        Matrix_2nd_color.fill(colors_t2[ith_P]);

        if ith_P == 0:
            zpos = np.zeros(lx * ly)
        elif ith_P > 0:
            zpos += dz

        dx = 0.8 * np.ones_like(zpos)
        # dy = dx.copy()
        dy = 0.05 * np.ones_like(zpos)
        dz = Matrix_2nd_SI.flatten()

        colors = Matrix_2nd_color.flatten()

        # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.5, edgecolor= [1.0, 1.0, 1.0, 0.0])
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, edgecolor=[1.0, 1.0, 1.0, 0.0])

    column_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(rxn.reaction_c.mapped_rxn_ind[i]) for i in xrange(np.shape(Matrix_2nd_SI)[0])]

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.w_zaxis.set_major_formatter(formatter)
    # ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) # doesn't work


    ax.set_xlabel('R')
    ax.set_ylabel('R')
    ax.set_zlabel('$2^{nd}$ order SI')
    ax.text(1, 1.0 - 0.1, max_hist * 1.0, "SUM($2^{nd}~$SI): " + '%1.4g' % (sum_topN_pathway))
    ax.text(1, 1.0 - 0.1, max_hist * 1.1, text_in)
    ax.margins(0.0)

    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    # curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    # if not os.path.exists(curr_dir):
    # os.makedirs(curr_dir)
    # filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, pathway covariance, version 1
"""


def bar3D_2D_PATHWAY_covariance_v1(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    #     formatter = ticker.ScalarFormatter(useMathText=True)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    pair_ind = []
    for i in xrange(N_variable1):
        for j in xrange(i + 1, N_variable2):
            pair_ind.append((i, j))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[pair_ind[i][0]][pair_ind[i][1]] = SI_2nd_all_PATHWAY[i]
        Matrix_2nd_SI[pair_ind[i][1]][pair_ind[i][0]] = SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        i_t = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[pair_ind[i_t][0]][pair_ind[i_t][1]] = colors_t2[i]
        Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]] = colors_t2[i]
    #         Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[1])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]

    xpos_t = np.arange(-1, lx + 1, 1)  # Set up a mesh of positions
    # xpos_t = np.arange(0,lx+2,1)    # Set up a mesh of positions
    ypos_t = np.arange(0, ly + 2, 1)
    xpos_t, ypos_t = np.meshgrid(xpos_t + 0.25, ypos_t + 0.25)  # for plot surface

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()
    colors = Matrix_2nd_color.flatten()

    offset = min(dz)
    zpos = [abs(offset) - abs(x) if x < 0 else abs(offset) for x in dz]
    #     dz= [abs(offset)+abs(x) if x<0 else abs(offset)+abs(x) for x in dz]
    dz = [abs(x) for x in dz]
    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.6, zsort='average', edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, zsort='average', edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.plot_surface(xpos_t, ypos_t, [[abs(offset)] * np.shape(xpos_t)[0]] * np.shape(ypos_t)[1], \
                    alpha=0.2, shade=True, cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.1, antialiased=False)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    #     ax.w_zaxis.set_major_formatter(formatter)
    #     ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) #doesn't work
    #     ax.w_zaxis.set_major_locator(ticker.LinearLocator(10))
    #     zticks_name= [ax.get_zlim()[1]-abs(offset) for i in xrange(10)]
    #     zticks_name= np.linspace(ax.get_zlim()[0]-abs(offset), ax.get_zlim()[1]-abs(offset), 10)
    #     ax.w_zaxis.set_ticklabels(zticks_name)

    #     formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
    order_of_magnitude = math.floor(math.log10(ax.get_zlim()[1]))
    #    print order_of_magnitude
    formatter = ticker.FuncFormatter(lambda x, pos: ("%.1f") % ((x + offset) / (10 ** order_of_magnitude)))
    #     formatter._set_offset(offset)
    #     formatter.set_useOffset(-offset)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    ax.w_zaxis.set_major_formatter(formatter)

    ax.set_xlabel('PATHWAY')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('covariance' + "$/10^{" + str(int(order_of_magnitude)) + "}$")
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0] + abs(offset)) * 1.0, "SUM(covariance): " + \
            '%1.2e' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0] + abs(offset)) * 1.1, text_in)
    ax.margins(0.0)
    #     ax.view_init(5.0, -27.5)


    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    #    plt.show()

    #     curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    #     if not os.path.exists(curr_dir):
    #         os.makedirs(curr_dir)
    #     filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, pathway covariance, version 2
"""


def bar3D_2D_PATHWAY_covariance_v2(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    pair_ind = []
    for i in xrange(N_variable1):
        for j in xrange(i + 1, N_variable2):
            pair_ind.append((i, j))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[pair_ind[i][0]][pair_ind[i][1]] = SI_2nd_all_PATHWAY[i]
        # Matrix_2nd_SI[pair_ind[i][1]][pair_ind[i][0]]= SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        i_t = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[pair_ind[i_t][0]][pair_ind[i_t][1]] = colors_t2[i]
        Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]] = colors_t2[i]
    #         Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[1])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]

    xpos_t = np.arange(-1, lx + 1, 1)  # Set up a mesh of positions
    # xpos_t = np.arange(0,lx+2,1)    # Set up a mesh of positions
    ypos_t = np.arange(0, ly + 2, 1)
    xpos_t, ypos_t = np.meshgrid(xpos_t + 0.25, ypos_t + 0.25)  # for plot surface

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()
    colors = Matrix_2nd_color.flatten()

    #     offset= min(dz)
    #     zpos= [abs(offset)-abs(x) if x<0 else abs(offset) for x in dz]
    # #     dz= [abs(offset)+abs(x) if x<0 else abs(offset)+abs(x) for x in dz]
    #     dz= [abs(x) for x in dz]
    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.6, zsort='average', edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, zsort='average', edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.plot_surface(xpos_t, ypos_t, [[0.0] * np.shape(xpos_t)[0]] * np.shape(ypos_t)[1], \
                    alpha=0.2, shade=True, cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.1, antialiased=False)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    #     ax.w_zaxis.set_major_formatter(formatter)
    #     ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) #doesn't work
    #     ax.w_zaxis.set_major_locator(ticker.LinearLocator(10))
    #     zticks_name= [ax.get_zlim()[1]-abs(offset) for i in xrange(10)]
    #     zticks_name= np.linspace(ax.get_zlim()[0]-abs(offset), ax.get_zlim()[1]-abs(offset), 10)
    #     ax.w_zaxis.set_ticklabels(zticks_name)

    #     formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
    #     order_of_magnitude= math.floor(math.log10(ax.get_zlim()[1]))
    #     formatter=ticker.FuncFormatter(lambda x, pos: ("%.1f")%((x+offset)/(10**order_of_magnitude)))

    #     formatter._set_offset(offset)
    #     formatter.set_useOffset(-offset)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    ax.w_zaxis.set_major_formatter(formatter)

    ax.set_xlabel('PATHWAY')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('covariance')
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.0, "SUM(covariance): " + \
            '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.1, text_in)
    ax.margins(0.0)
    ax.set_zbound(min(dz), max(dz))
    #     ax.view_init(5.0, -27.5)


    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    #    plt.show()

    #     curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    #     if not os.path.exists(curr_dir):
    #         os.makedirs(curr_dir)
    #     filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, pathway correlation, version 2
"""


def bar3D_2D_PATHWAY_correlation_v2(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    pair_ind = []
    for i in xrange(N_variable1):
        for j in xrange(i + 1, N_variable2):
            pair_ind.append((i, j))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[pair_ind[i][0]][pair_ind[i][1]] = SI_2nd_all_PATHWAY[i]
        # Matrix_2nd_SI[pair_ind[i][1]][pair_ind[i][0]]= SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        i_t = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[pair_ind[i_t][0]][pair_ind[i_t][1]] = colors_t2[i]
        Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]] = colors_t2[i]
    #         Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[1])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]

    xpos_t = np.arange(-1, lx + 1, 1)  # Set up a mesh of positions
    # xpos_t = np.arange(0,lx+2,1)    # Set up a mesh of positions
    ypos_t = np.arange(0, ly + 2, 1)
    xpos_t, ypos_t = np.meshgrid(xpos_t + 0.25, ypos_t + 0.25)  # for plot surface

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()
    colors = Matrix_2nd_color.flatten()

    #     offset= min(dz)
    #     zpos= [abs(offset)-abs(x) if x<0 else abs(offset) for x in dz]
    # #     dz= [abs(offset)+abs(x) if x<0 else abs(offset)+abs(x) for x in dz]
    #     dz= [abs(x) for x in dz]
    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.6, zsort='average', edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, zsort='average', edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.plot_surface(xpos_t, ypos_t, [[0.0] * np.shape(xpos_t)[0]] * np.shape(ypos_t)[1], \
                    alpha=0.2, shade=True, cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.1, antialiased=False)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    #     ax.w_zaxis.set_major_formatter(formatter)
    #     ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) #doesn't work
    #     ax.w_zaxis.set_major_locator(ticker.LinearLocator(10))
    #     zticks_name= [ax.get_zlim()[1]-abs(offset) for i in xrange(10)]
    #     zticks_name= np.linspace(ax.get_zlim()[0]-abs(offset), ax.get_zlim()[1]-abs(offset), 10)
    #     ax.w_zaxis.set_ticklabels(zticks_name)

    #     formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
    #     order_of_magnitude= math.floor(math.log10(ax.get_zlim()[1]))
    #     formatter=ticker.FuncFormatter(lambda x, pos: ("%.1f")%((x+offset)/(10**order_of_magnitude)))

    #     formatter._set_offset(offset)
    #     formatter.set_useOffset(-offset)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    ax.w_zaxis.set_major_formatter(formatter)

    ax.set_xlabel('PATHWAY')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('correlation')
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.0, "SUM(correlation): " + \
            '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.1, text_in)
    ax.margins(0.0)
    ax.set_zbound(min(dz), max(dz))
    #     ax.view_init(5.0, -27.5)


    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    #    plt.show()

    #     curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    #     if not os.path.exists(curr_dir):
    #         os.makedirs(curr_dir)
    #     filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, SI from pathway covariance
"""


def bar3D_2D_PATHWAY_covariance_SI(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    pair_ind = []
    for i in xrange(N_variable1):
        for j in xrange(i + 1, N_variable2):
            pair_ind.append((i, j))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[pair_ind[i][0]][pair_ind[i][1]] = SI_2nd_all_PATHWAY[i]
        # Matrix_2nd_SI[pair_ind[i][1]][pair_ind[i][0]]= SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        i_t = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[pair_ind[i_t][0]][pair_ind[i_t][1]] = colors_t2[i]
        Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]] = colors_t2[i]
    #         Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[1])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]

    xpos_t = np.arange(-1, lx + 1, 1)  # Set up a mesh of positions
    # xpos_t = np.arange(0,lx+2,1)    # Set up a mesh of positions
    ypos_t = np.arange(0, ly + 2, 1)
    xpos_t, ypos_t = np.meshgrid(xpos_t + 0.25, ypos_t + 0.25)  # for plot surface

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()
    colors = Matrix_2nd_color.flatten()

    #     offset= min(dz)
    #     zpos= [abs(offset)-abs(x) if x<0 else abs(offset) for x in dz]
    # #     dz= [abs(offset)+abs(x) if x<0 else abs(offset)+abs(x) for x in dz]
    #     dz= [abs(x) for x in dz]
    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.6, zsort='average', edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, zsort='average', edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.plot_surface(xpos_t, ypos_t, [[0.0] * np.shape(xpos_t)[0]] * np.shape(ypos_t)[1], \
                    alpha=0.2, shade=True, cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.1, antialiased=False)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    #     ax.w_zaxis.set_major_formatter(formatter)
    #     ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) #doesn't work
    #     ax.w_zaxis.set_major_locator(ticker.LinearLocator(10))
    #     zticks_name= [ax.get_zlim()[1]-abs(offset) for i in xrange(10)]
    #     zticks_name= np.linspace(ax.get_zlim()[0]-abs(offset), ax.get_zlim()[1]-abs(offset), 10)
    #     ax.w_zaxis.set_ticklabels(zticks_name)

    #     formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
    #     order_of_magnitude= math.floor(math.log10(ax.get_zlim()[1]))
    #     formatter=ticker.FuncFormatter(lambda x, pos: ("%.1f")%((x+offset)/(10**order_of_magnitude)))

    #     formatter._set_offset(offset)
    #     formatter.set_useOffset(-offset)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    ax.w_zaxis.set_major_formatter(formatter)

    ax.set_xlabel('PATHWAY')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('covariance SI')
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.0, "SUM(covariance SI): " + \
            '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.1, text_in)
    ax.margins(0.0)
    ax.set_zbound(min(dz), max(dz))
    #     ax.view_init(5.0, -27.5)


    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    #    plt.show()

    #     curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    #     if not os.path.exists(curr_dir):
    #         os.makedirs(curr_dir)
    #     filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)


"""
sensitivity index plot, 2D, bar3D, redefine colors, SI from pathway covariance
Symmetric
"""


def bar3D_2D_PATHWAY_covariance_SI_v2(SI_2nd_all_PATHWAY, N_variable1, N_variable2, filename3, text_in=""):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    pair_ind = []
    for i in xrange(N_variable1):
        for j in xrange(i + 1, N_variable2):
            pair_ind.append((i, j))

    # SI Matrix element
    Matrix_2nd_SI = np.zeros(N_variable1 * N_variable2).reshape((N_variable1, N_variable2))
    # symmetric matrix
    for i in xrange(len(SI_2nd_all_PATHWAY)):
        Matrix_2nd_SI[pair_ind[i][0]][pair_ind[i][1]] = SI_2nd_all_PATHWAY[i]
        Matrix_2nd_SI[pair_ind[i][1]][pair_ind[i][0]] = SI_2nd_all_PATHWAY[i]

        # color Matrix element
    #     topN_ith_pair_PATHWAY= 5
    topN_ith_pair_PATHWAY = len(SI_2nd_all_PATHWAY)
    SI_2nd_all_topN_ele_PATHWAY, SI_2nd_all_topN_ind_PATHWAY = \
        mu.my_utility_c.topN_element_and_index(SI_2nd_all_PATHWAY, topN_ith_pair_PATHWAY)
    # color, just topN kinds of colors
    colors_t2 = matplotlib.cm.rainbow(np.linspace(0, 1, topN_ith_pair_PATHWAY))  # top 10
    colors_t2 = colors_t2[::-1]

    Matrix_2nd_color = np.empty([np.shape(Matrix_2nd_SI)[0], np.shape(Matrix_2nd_SI)[1]], dtype=np.ndarray)
    Matrix_2nd_color.fill(colors_t2[-1]);
    # Modify color matrix, symmetric matrix
    for i in xrange(len(SI_2nd_all_topN_ind_PATHWAY)):
        i_t = SI_2nd_all_topN_ind_PATHWAY[i]
        Matrix_2nd_color[pair_ind[i_t][0]][pair_ind[i_t][1]] = colors_t2[i]
        Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]] = colors_t2[i]
    #         Matrix_2nd_color[pair_ind[i_t][1]][pair_ind[i_t][0]]= [1.0, 1.0, 1.0, 0.0]

    column_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[0])]
    row_names = [str(i) for i in xrange(np.shape(Matrix_2nd_SI)[1])]

    fig = plt.figure()
    ax = Axes3D(fig)

    lx = np.shape(Matrix_2nd_SI)[0]  # Work out matrix dimensions
    ly = np.shape(Matrix_2nd_SI)[1]

    xpos_t = np.arange(-1, lx + 1, 1)  # Set up a mesh of positions
    # xpos_t = np.arange(0,lx+2,1)    # Set up a mesh of positions
    ypos_t = np.arange(0, ly + 2, 1)
    xpos_t, ypos_t = np.meshgrid(xpos_t + 0.25, ypos_t + 0.25)  # for plot surface

    xpos = np.arange(0, lx, 1)  # Set up a mesh of positions
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    xpos = xpos.flatten()  # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    dx = 0.8 * np.ones_like(zpos)
    # dy = dx.copy()
    dy = 0.05 * np.ones_like(zpos)
    dz = Matrix_2nd_SI.flatten()
    colors = Matrix_2nd_color.flatten()

    #     offset= min(dz)
    #     zpos= [abs(offset)-abs(x) if x<0 else abs(offset) for x in dz]
    # #     dz= [abs(offset)+abs(x) if x<0 else abs(offset)+abs(x) for x in dz]
    #     dz= [abs(x) for x in dz]
    # ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=colors, alpha= 0.6, zsort='average', edgecolor= [1.0, 1.0, 1.0, 0.0])
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.6, zsort='average', edgecolor=[1.0, 1.0, 1.0, 0.0])

    ax.plot_surface(xpos_t, ypos_t, [[0.0] * np.shape(xpos_t)[0]] * np.shape(ypos_t)[1], \
                    alpha=0.2, shade=True, cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.1, antialiased=False)

    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    #     ax.w_zaxis.set_major_formatter(formatter)
    #     ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2)) #doesn't work
    #     ax.w_zaxis.set_major_locator(ticker.LinearLocator(10))
    #     zticks_name= [ax.get_zlim()[1]-abs(offset) for i in xrange(10)]
    #     zticks_name= np.linspace(ax.get_zlim()[0]-abs(offset), ax.get_zlim()[1]-abs(offset), 10)
    #     ax.w_zaxis.set_ticklabels(zticks_name)

    #     formatter = ticker.ScalarFormatter(useMathText=True, useOffset=True)
    #     order_of_magnitude= math.floor(math.log10(ax.get_zlim()[1]))
    #     formatter=ticker.FuncFormatter(lambda x, pos: ("%.1f")%((x+offset)/(10**order_of_magnitude)))

    #     formatter._set_offset(offset)
    #     formatter.set_useOffset(-offset)
    #     formatter.set_scientific(True)
    #     formatter.set_powerlimits((-2,2))
    ax.w_zaxis.set_major_formatter(formatter)

    ax.set_xlabel('PATHWAY')
    ax.set_ylabel('PATHWAY')
    ax.set_zlabel('covariance SI')
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.0, "SUM(covariance SI): " + \
            '%1.4g' % (np.sum(Matrix_2nd_SI)))
    ax.text(1, 1.0 - 0.1, (SI_2nd_all_topN_ele_PATHWAY[0]) * 1.1, text_in)
    ax.margins(0.0)
    ax.set_zbound(min(dz), max(dz))
    #     ax.view_init(5.0, -27.5)


    ticksx = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.xticks(ticksx, column_names, fontsize=6, rotation=45)

    ticksy = np.arange(0.0, np.shape(Matrix_2nd_SI)[0], 1) + 0.5
    plt.yticks(ticksy, row_names, fontsize=6, rotation=-15)

    #    plt.show()

    #     curr_dir= os.path.join(file_dir, "output", species_type+"_conc_SI")
    #     if not os.path.exists(curr_dir):
    #         os.makedirs(curr_dir)
    #     filename3= os.path.join(curr_dir, "SI_2nd_conc_"+species_type+"_3D.png")
    fig.savefig(filename3, dpi=1000, bbox_inches="tight")
    plt.close(fig)
