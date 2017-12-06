
"""
fast transitions
"""
from collections import OrderedDict


def get_fast_rxn_trapped_spe():
    """
    return fast reaction and trapped species
    better make them in the same order
    """
    fast_transitions = [
        {}
    ]

    fast_reaction_list = []
    trapped_species_list = []
    for _, r_s in enumerate(fast_transitions):
        print(r_s)
        fast_reaction_list.append(
            tuple([str(r_s['rxn'][0]), int(r_s['rxn'][1])]))
        trapped_species_list.append(
            tuple([str(r_s['spe'][0]), int(r_s['spe'][1])]))

    fast_reaction = OrderedDict(fast_reaction_list)
    trapped_species = OrderedDict(trapped_species_list)
    # print(fast_reaction, trapped_species)
    return fast_reaction, trapped_species


if __name__ == '__main__':
    get_fast_rxn_trapped_spe()
