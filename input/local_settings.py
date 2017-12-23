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
        "traj_max_t": 0.85,
        # trajectory critical time, after which print out more data points
        "traj_critical_t": 0.751999999880706205,
        # reference time, to a combustion system, this is gonna be the ignition delay time
        # for Propane, time when temperature=1800K
        "tau": 0.85,
        # beginning time, for pathway or for trajectory, exact time = begin_t*max_tau
        "begin_t": 0.0,
        # end time, for pathway or for trajectory, exact time = end_t*max_tau
        "end_t": 1.0,
        # species oriented, if true, pick pathways ending with top_n species,
        #  if False, just top n pathway
        "spe_oriented": True,
        # condense species path, no reactions
        "species_path": False,
        # atom followed
        "atom_f": "C",
        "init_s": 62,
        # end species index, either None, or [] or [14, 15]
        "end_s_idx": [14],
        # top n path
        "top_n_p": 10,
        # top n path for gephi to generate coordinates
        "top_n_p_gephi": 10,
        # top n species
        "top_n_s": 10,
        # number of trajectory used to generate pathway list running mc simulation
        "mc_n_traj": 1000000,
        # path integral number of trajectory
        "pi_n_traj": 10000,
        # number of time points when prepare path integral time points
        "pi_n_time": 1,
        # tag, M or fraction
        "tag": "M"
    }
    return setting


def get_fast_rxn_trapped_spe():
    """
    return fast reaction and trapped species
    better make them in the same order
    """
    fast_transitions = [
        {
            # 1068    549     O2+npropyl=npropyloo
            # reactants       9       O2      60      npropyl products        78      npropyloo
            # 1069    -549    O2+npropyl=npropyloo
            "rxn": [1068, 1069],
            "spe": [60, 78]
        },
        # 1096    565     O2+ipropyl=ipropyloo
        # reactants       9       O2      61      ipropyl products        80      ipropyloo
        # 1097    -565    O2+ipropyl=ipropyloo
        {
            "rxn": [1096, 1097],
            "spe": [61, 80]
        },
        # 132     69      CH3+O2(+M)=CH3OO(+M)
        # reactants       25      CH3     9       O2      products        27      CH3OO
        # 133     -69     CH3+O2(+M)=CH3OO(+M)
        {
            "rxn": [132, 133],
            "spe": [25, 27]
        },
        # 1116    575     O2+QOOH_1=well_1
        # reactants       9       O2      87      QOOH_1  products        90      well_1
        # 1117    -575    O2+QOOH_1=well_1
        {
            "rxn": [1116, 1117],
            "spe": [87, 90]
        },
        # 348     180     C2H5+O2=CH3CH2OO
        # reactants       39      C2H5    9       O2      products        50      CH3CH2OO
        # 349     -180    C2H5+O2=CH3CH2OO
        {
            "rxn": [348, 349],
            "spe": [39, 50]
        },
        # 1080    556     npropyloo=QOOH_1        557     npropyloo=QOOH_1
        # reactants       78      npropyloo       products        87      QOOH_1
        # 1081    -556    npropyloo=QOOH_1        -557    npropyloo=QOOH_1
        {
            "rxn": [1080, 1081],
            "spe": [78, 87]
        },
        # 586     300     O2C2H4OH=CH2CH2OH+O2
        # reactants       85      O2C2H4OH        products        54      CH2CH2OH        9       O2
        # 587     -300    O2C2H4OH=CH2CH2OH+O2
        {
            "rxn": [586, 587],
            "spe": [85, 54]
        },
        # 1042    536     allyloxy=vinoxylmethyl
        # reactants       72      allyloxy        products        108     vinoxylmethyl
        # 1043    -536    allyloxy=vinoxylmethyl
        {
            "rxn": [1042, 1043],
            "spe": [72, 108]
        },
        # 434     224     acetyl+O2=acetylperoxy
        # reactants       9       O2      45      acetyl  products        47      acetylperoxy
        # 435     -224    acetyl+O2=acetylperoxy
        {
            "rxn": [434, 435],
            "spe": [45, 47]
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
    # print(fast_reaction, trapped_species)
    return fast_reaction, trapped_species


if __name__ == '__main__':
    get_fast_rxn_trapped_spe()
