import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Ellipse
import os


class pathway_prob_c:
    # read pathway data as pandas frame object
    @staticmethod
    def read_pathway_as_pandas_frame_object(filename):
        return pd.read_csv(filename, delimiter=",", names=["pathway", "number"])

    # convert pandas frame object to list
    @staticmethod
    def convert_pandas_frame_object_to_list(pathway_data):
        # parse the last symbol "#" and get the max index, #K's =maxIndex+1
        mask = pathway_data.pathway.str.contains("#")
        MaxIndex = pathway_data[mask].tail(1).pathway.str.replace("#", "").tolist()[0]
        MaxIndex = int(MaxIndex)

        # convert Object to list, Store the pathway data in a list indexed with their own index
        pathway = []
        for i in xrange(MaxIndex):
            mask_split1 = pathway_data.index >= pathway_data[pathway_data.pathway == ("#" + str(i))].index
            mask_split2 = pathway_data.index < pathway_data[pathway_data.pathway == ("#" + str(i + 1))].index
            pathway.append(pathway_data[(mask_split1) & (mask_split2)])

        # the last element
        mask_split2 = pathway_data.index >= pathway_data[pathway_data.pathway == ("#" + str(MaxIndex))].index
        pathway.append(pathway_data[mask_split2])

        return pathway

    # get name and normalized average prob of all pathway from pandas frame object
    @staticmethod
    def get_name_prob_of_all_pathway(pathway_data):
        mask = pathway_data.pathway.str.contains("#")
        pathway_t = pathway_data[~mask]
        pathway_name = pathway_t.pathway.tolist()
        pathway_number = pathway_t.number.tolist()

        # pathway dictionary, no duplicated pathway
        pathway_dict = dict()

        for i in xrange(len(pathway_name)):
            if not pathway_dict.has_key(pathway_name[i]):
                pathway_dict[pathway_name[i]] = pathway_number[i]
            else:
                pathway_dict[pathway_name[i]] += pathway_number[i]

        # sort pathway_dict, return sorted tuple pair
        pathway_tuple_sorted = sorted(pathway_dict.items(), key=lambda x: x[1], reverse=True)
        # get pathway list
        all_pathway_name = [x[0] for x in pathway_tuple_sorted]
        all_pathway_number = [x[1] for x in pathway_tuple_sorted]
        all_pathway_number = all_pathway_number / np.sum(all_pathway_number)

        return all_pathway_name, all_pathway_number

    # get a list of prob of Nth pathway wrt. to all sets of K's
    @staticmethod
    def get_list_of_prob_of_a_pathway(pathway, all_pathway_name, Nth_pathway):
        prob_of_a_pathway = np.zeros(len(pathway))
        for i in xrange(len(pathway)):
            # this K's contains current pathway 
            if pathway[i].apply(lambda x: all_pathway_name[Nth_pathway] in x.values, axis=1).any():
                prob_of_a_pathway[i] = pathway[i][pathway[i].pathway == all_pathway_name[Nth_pathway]].number.tolist()[
                    0]

        return prob_of_a_pathway

    # filters
    # get name and normalized average prob of all pathway from pandas frame object
    @staticmethod
    def get_name_prob_of_pathways_contains_R21(pathway_data):
        mask1 = pathway_data.pathway.str.contains("#")
        # reaction R21->R40, R21*->R41 in pathway network labelling
        mask20 = pathway_data.pathway.str.contains("R40")
        mask21 = pathway_data.pathway.str.contains("R41")
        pathway_t = pathway_data[(~mask1) & (mask20 | mask21)]
        pathway_name = pathway_t.pathway.tolist()
        pathway_number = pathway_t.number.tolist()

        # pathway dictionary, no duplicated pathway
        pathway_dict = dict()

        for i in xrange(len(pathway_name)):
            if not pathway_dict.has_key(pathway_name[i]):
                pathway_dict[pathway_name[i]] = pathway_number[i]
            else:
                pathway_dict[pathway_name[i]] += pathway_number[i]

        # sort pathway_dict, return sorted tuple pair
        pathway_tuple_sorted = sorted(pathway_dict.items(), key=lambda x: x[1], reverse=True)
        # get pathway list
        # comprehensive list
        all_pathway_name = [x[0] for x in pathway_tuple_sorted]
        all_pathway_number = [x[1] for x in pathway_tuple_sorted]
        all_pathway_number = all_pathway_number / np.sum(all_pathway_number)

        return all_pathway_name, all_pathway_number

        #  plot pathway probability cdf

    @staticmethod
    def plot_cdf(all_pathway_prob_list, species_type, file_dir, topN=20):
        # new allocation, apply for new memory
        all_pathway_prob_cdf = list(all_pathway_prob_list)
        for i in xrange(1, len(all_pathway_prob_cdf)):
            all_pathway_prob_cdf[i] = all_pathway_prob_cdf[i - 1] + all_pathway_prob_cdf[i]
        topN_CDF = topN
        all_pathway_prob_cdf = all_pathway_prob_cdf[0:topN_CDF]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(np.arange(topN_CDF), all_pathway_prob_cdf[0:topN_CDF], marker='o', color='#D95F02', label="cdf",
                 markersize=1.0)
        ax1.plot(np.arange(topN_CDF), all_pathway_prob_list[0:topN_CDF], marker='o', color='#7bbfea', label="pdf",
                 markersize=1.0)
        criteria = 0.90;
        ax1.hlines(criteria, 0, topN_CDF, color='#1d953f', label=str(criteria))
        #    criteria= 0.99; ax1.hlines(criteria, 0, topN_CDF, color= 'r', label=str(criteria))

        ax1.legend(bbox_to_anchor=(1.0, 1.0), loc=2, borderaxespad=0.)

        ax1.set_title(species_type)
        ax1.set_xlabel("pathway index")
        ax1.set_ylabel("pathway prob")

        # inset plot
        ax2 = plt.axes([.3, .3, .5, .5], axisbg='y')
        topN_CDF_t = 15
        ax2.plot(np.arange(topN_CDF_t), all_pathway_prob_cdf[0:topN_CDF_t], marker='o', color='#D95F02', label="cdf",
                 markersize=1.0)
        ax2.plot(np.arange(topN_CDF_t), all_pathway_prob_list[0:topN_CDF_t], marker='o', color='#7bbfea', label="pdf",
                 markersize=1.0)
        criteria = 0.90
        ax2.hlines(criteria, 0, topN_CDF_t, color='#1d953f', label=str(criteria))
        #    criteria= 0.95; ax1.hlines(criteria, 0, topN_CDF_t, color= '#1d953f', label=str(criteria))

        # arrow, annotate
        src_pos, dst_pos = (0, 0.5), (50, 0.6)
        #     src_patch = plt.Rectangle(src_pos, 5, .10, color='r', clip_on=False)
        # src_patch = Ellipse(src_pos, width= 60, height= 1.2, angle= 0, color='#FFFFCC', edgecolor= 'r', clip_on=False)
        src_patch = Ellipse(src_pos, width=60, height=1.2, angle=0, facecolor='#FFFFFF', edgecolor='#FF6633',
                            clip_on=False)
        ax1.add_patch(src_patch)
        # dst_patch = Ellipse(dst_pos, width= 23, height= 0.1, angle= 0, color='#66FFFF', clip_on=False)
        # ax1.add_patch(dst_patch)
        arrowprops = dict(
            arrowstyle='fancy',
            connectionstyle='arc3,rad=-0.8',
            facecolor='r',
            #         patchA=dst_patch,
            patchA=src_patch,
            shrinkA=1,
            shrinkB=1)
        #     ant = ax1.annotate('haha', xy=src_pos, xytext= dst_pos, arrowprops=arrowprops)
        ax1.annotate('', xy=dst_pos, xytext=src_pos, arrowprops=arrowprops)

        fig.subplots_adjust(right=0.83)

        # curr_dir= os.path.join(file_dir, "output", "dist", species_type)
        # if not os.path.exists(curr_dir):
        #    os.makedirs(curr_dir)
        fig.savefig(os.path.join(file_dir, species_type + "_pdf_cdf.jpg"), dpi=600)

    @staticmethod
    def read_pathname_s2f_path_endswith(filename1, filename2, species_index):
        p_pfo = pathway_prob_c.read_pathway_as_pandas_frame_object(filename1)

        mask1 = p_pfo.pathway.str.endswith("S" + str(species_index))
        p_pfo[mask1].to_csv(filename2, header=False, index=False, sep=',', columns=['pathway', 'number'])
