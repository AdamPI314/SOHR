
"""
fast transitions
"""
from collections import OrderedDict


def get_fast_rxn_chattering_spe():
    """
    return fast reaction and chattering species
    better make them in the same order
    """
    fast_transitions = [
        {}
    ]

    fast_reaction_list = []
    chattering_species_list = []
    for _, r_s in enumerate(fast_transitions):
        print(r_s)
        fast_reaction_list.append(
            tuple([str(r_s['rxn'][0]), int(r_s['rxn'][1])]))
        chattering_species_list.append(
            tuple([str(r_s['spe'][0]), int(r_s['spe'][1])]))

    fast_reaction = OrderedDict(fast_reaction_list)
    chattering_species = OrderedDict(chattering_species_list)
    # print(fast_reaction, chattering_species)
    return fast_reaction, chattering_species


if __name__ == '__main__':
    get_fast_rxn_chattering_spe()
