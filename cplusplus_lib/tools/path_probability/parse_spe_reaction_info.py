import pandas as pd
import numpy as np
import re
#import os
#import sys

# parse species
def parse_spe_info(filename):
    """
    parse species info from file= "os.path.join(file_dir, "output", "species_labelling.dat")"
    """
    line_content= np.genfromtxt(filename, dtype=str, delimiter='\n')
    
    matched_str= [re.findall("(\d+)\t-->\t([\w|\-|(|)]+)", line)[0] for line in line_content]
    matched_str_reverse=[(x[1], x[0]) for x in matched_str]
    
    spe_ind_name_dict= dict(matched_str)
    spe_name_ind_dict= dict(matched_str_reverse)
    
    return spe_ind_name_dict, spe_name_ind_dict
    
    
# parse reactions
def parse_rxn_and_its_index(filename):
    """
    parse reaction info from file= "os.path.join(file_dir, "output", "reaction_labelling.dat")"
    """
    #load data
    # line_content= [x.rstrip('\n') for x in open(filename)]
    line_content= np.genfromtxt(filename, dtype=str, delimiter='\n')
    #use regex to parse reaction itself
    matched_rxn_str= [re.findall("(\d+):\s+([\w|+|=|>|<|(|)|\-|,]+)", line)[0] for line in line_content \
                  if len(re.findall("(\d+):\s+([\w|+|=|>|<|(|)|\-|,]+)", line))!=0]
    # reaction index
    matched_index_str= [re.findall("[=|<|>]+\s+(\d+)\s*(\d+)?", line)[0] for line in line_content \
                    if len(re.findall("[=|<|>]+\s+(\d+)\s*(\d+)?", line))!=0]
    #map the old and new reaction index, basically create two dictionaries
    old_new_index_dict= dict([(matched_rxn_str[i][0], matched_index_str[i]) for i in xrange(len(matched_rxn_str))])
    #map the new old reaction index
    new_old_index_dict= dict()
    for i in xrange(len(matched_rxn_str)):
        new_old_index_dict.update({matched_index_str[i][0]:str(matched_rxn_str[i][0])})
        if(matched_index_str[i][1]!=''):
            new_old_index_dict.update({matched_index_str[i][1]:str(matched_rxn_str[i][0])+'*'})
    #reaction's new index and reaction itself
    #parse reaction arrow first
    matched_rxn_arrow= [re.findall("([=|>|<]+)", ind_rxn[1]) for ind_rxn in matched_rxn_str]
    #split reactant and product
    reactant_product= [re.findall("([\w|+|(|)|\-|,]+)[=|>|<]+([\w|+|(|)|\-|,]+)", ind_rxn[1]) \
                   for ind_rxn in matched_rxn_str]
    #map reaction new reaction label and the exact reaction
    new_ind_rxn_dict=dict()
    for i in xrange(len(matched_index_str)):
        if matched_rxn_arrow[i][0]=='=' or matched_rxn_arrow[i][0]=='<=>':
            new_ind_rxn_dict.update({matched_index_str[i][0]:reactant_product[i][0][0]+'=>'+reactant_product[i][0][1]})
            new_ind_rxn_dict.update({matched_index_str[i][1]:reactant_product[i][0][1]+'=>'+reactant_product[i][0][0]})
        elif matched_rxn_arrow[i][0]=='=>':
            new_ind_rxn_dict.update({matched_index_str[i][0]:reactant_product[i][0][0]+'=>'+reactant_product[i][0][1]})
        elif matched_rxn_arrow[i][0]=='<=':
            new_ind_rxn_dict.update({matched_index_str[i][1]:reactant_product[i][0][1]+'=>'+reactant_product[i][0][0]})
            
    return old_new_index_dict, new_old_index_dict, new_ind_rxn_dict
    

def PATH_to_real_spe_rxn(spe_ind_name_dict, new_ind_rxn_dict, pathway_name):
    """
    converted path to their real species name and reaction format instead of index
    """
    matched_S= re.findall("S(\d+)", pathway_name)
    matched_R= re.findall("R(\d+)", pathway_name)
    #print matched_S
    #print matched_R
    # always starts from species
    str_t= '['+spe_ind_name_dict[matched_S[0]]+'] '
    for i in xrange(len(matched_R)):
        str_t+=new_ind_rxn_dict[matched_R[i]]
        str_t+="-->"
        str_t+='['+spe_ind_name_dict[matched_S[i+1]]+'] '
        
    return str_t
    
def read_pathname_convert_2_real_spe_rxn(filename_spe, filename_rxn, filename_p, filename_p_out, topN= 50):
    """
    read species and reaction info, convert path info into real species and reaction instead of index and write to file
    """
    #load path data
    path_data= pd.read_csv(filename_p, names=['path', 'prob'])
    total_prob= sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    # renormalize
    path_data['prob']/=total_prob
    
    #load spe and reaction info
    spe_ind_name_dict, spe_name_ind_dict= parse_spe_info(filename_spe)
    old_new_index_dict, new_old_index_dict, new_ind_rxn_dict= parse_rxn_and_its_index(filename_rxn)
    
    #convert species reaction index to real species and reactions
    path_data['path']= path_data['path'].apply(lambda x:PATH_to_real_spe_rxn(spe_ind_name_dict, new_ind_rxn_dict, x))
    
    #write to file
    path_data[0:topN].to_csv(filename_p_out, header= False, index= False, sep= ';', columns=['path', 'prob'])
