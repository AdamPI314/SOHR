"""
local settings, include settings from C++ code and fast transitions
"""
from collections import OrderedDict


def get_local_settings():
    """
    local settings for C++ codes
    """
    setting = {
        "system": {
            "condition": "cv",
            "initializer": "dlsode"
        },
        "network": {
            "merge_chatterings": "yes"
        },
        "propagator": {
            "primary_type": "from_file",
            "type": "dlsode",
            "sub_type": "time_propagator_cv_s2m_pgt",
            "convert_molar_concentration_to_mole_fraction": "no",
            "normalize_initial_concentration": "yes"
        },
        # trajectory max time, used to solve referene trajectory
        "traj_max_t": 4.375545109e-3,
        # trajectory critical time, after which print out more data points
        "traj_critical_t": 4.375545109e-3,
        # reference time, to a combustion system, this is gonna be the ignition delay time
        # for Propane, time when temperature=1800K
        "tau": 4.375545109e-3,
        # time at which using MC to generate pathway list, time=mc_t*tau
        "mc_t": 0.8,
        # beginning time, for pathway or for trajectory, exact time = begin_t*tau
        "begin_t": 0.0,
        # end time, for pathway or for trajectory, exact time = end_t*tau
        # here 0.25718313951098054 is actually 0.2 seconds
        "end_t": 0.8,
        # "end_t": 0.9,
        # species oriented, if true, pick pathways ending with top_n species,
        #  if False, just top n pathway
        "spe_oriented": False,
        # condense species path, no reactions
        "species_path": False,
        # atom followed
        "atom_f": "H",
        "init_s": 2,
        # terminal species for file ./setting.json, either None, or [] or [14, 15]
        "terminal_spe": [],
        # end species index, either None, or [] or [14, 15]
        "end_s_idx": [],
        # top n path
        "top_n_p": 500,
        # top n path for gephi to generate coordinates
        "top_n_p_gephi": 500,
        # top n species
        "top_n_s": 10,
        # number of trajectory used to generate pathway list running mc simulation
        "mc_n_traj": 1000000,
        # path integral number of trajectory
        "pi_n_traj": 10000,
        # number of time points when prepare path integral time points
        "pi_n_time": 50,
        # tag, M or fraction
        "tag": "M"
    }
    return setting


def get_chattering_species(atom_followed="C"):
    """
    return chattering species
    the chatteing reaction infomation is just for reference, will not use it
    as long as the paired chattering species is provided, should be fine
    better make them in the same order
    """
    fast_transitions = [{}]

    trapped_species_list = []
    for _, r_s in enumerate(fast_transitions):
        print(r_s)
        if 'spe' not in r_s:
            continue
        if atom_followed not in r_s['spe']:
            continue
        if len(r_s['spe'][atom_followed]) != 2:
            continue

        trapped_species_list.append(
            [int(r_s['spe'][atom_followed][0]), int(r_s['spe'][atom_followed][1])])
    print(trapped_species_list)

    chattering_species = {}
    for idx, val in enumerate(trapped_species_list):
        print(idx, val)
        chattering_species.update({str(idx + 1): val})

    chattering_species = OrderedDict(chattering_species)
    print(chattering_species)
    return chattering_species


if __name__ == '__main__':
    get_chattering_species("HA2")
