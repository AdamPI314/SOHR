"""
local settings, include settings from C++ code and fast transitions
"""
from collections import OrderedDict


def get_local_settings():
    """
    local settings for C++ codes
    """
    setting = {
        # end time
        "end_t": 1.0e08,
        # critical time, after which print out more data points
        "critical_t": 1.0e08,
        # reference time, to a combustion system, this is gonna be the ignition delay time
        # for Propane, time when temperature=1800K
        "max_tau": 38.8,
        # exact time = tau*max_tau
        "tau": 1.0,
        # species oriented, if true, pick pathways ending with top_n species,
        #  if False, just top n pathway
        "spe_oriented": False,
        # condense species path, no reactions
        "species_path": False,
        # atom followed
        "atom_f": "X",
        "init_s": 0,
        # end species index, either None, or [] or [14, 15]
        "end_s_idx": [],
        # top n path
        "top_n_p": 10,
        # top n path for gephi to generate coordinates
        "top_n_p_gephi": 10,
        # top n species
        "top_n_s": 10,
        # number of trajectory used to generate pathway list running mc simulation
        "mc_n_traj": 10000000,
        # path integral number of trajectory
        "pi_n_traj": 10000,
        # number of time points when prepare path integral time points
        "pi_n_time": 10,
        # tag, M or fraction
        "tag": "ssa_number"
    }
    return setting


def get_fast_rxn_trapped_spe():
    """
    return fast reaction and trapped species
    better make them in the same order
    """
    fast_transitions = [
        {
            # 0	1	2A=>A2
            # reactants	2	A	products	4	A2
            # 2	3	A2=>2A
            # reactants	4	A2	products	2	A
            "rxn": [0, 2],
            "spe": [2, 4]
        },
        {
            # 4	5	O+A2=>OA2
            # reactants	4	A2	1	O	products	6	OA2
            # 6	7	OA2 = >O + A2
            # reactants	6	OA2	products	4	A2	1	O
            "rxn": [4, 6],
            "spe": [4, 6]
        },
        {
            # 1	2	2B = >B2
            # reactants	3	B	products	5	B2
            # 3	4	B2 = >2B
            # reactants	5	B2	products	3	B
            "rxn": [1, 3],
            "spe": [3, 5]
        },
        {
            # 5	6	O + B2= >OB2
            # reactants	5	B2	1	O	products	7	OB2
            # 7	8	OB2 = >O + B2
            # reactants	7	OB2	products	5	B2	1	O
            "rxn": [5, 7],
            "spe": [5, 7]
        }
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
    print(fast_reaction, trapped_species)
    return fast_reaction, trapped_species


if __name__ == '__main__':
    get_fast_rxn_trapped_spe()
