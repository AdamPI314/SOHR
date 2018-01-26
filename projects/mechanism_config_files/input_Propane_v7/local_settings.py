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
        "traj_max_t": 0.779074999626780951,
        # trajectory critical time, after which print out more data points
        "traj_critical_t": 0.751999999880706205,
        # reference time, to a combustion system, this is gonna be the ignition delay time
        # for Propane, time when temperature=1800K
        "tau": 0.777655955130997,
        # time at which using MC to generate pathway list, time=mc_t*tau
        "mc_t": 0.25718313951098054,
        # beginning time, for pathway or for trajectory, exact time = begin_t*tau
        "begin_t": 0.0,
        # end time, for pathway or for trajectory, exact time = end_t*tau
        # here 0.25718313951098054 is actually 0.2 seconds
        "end_t": 1.001,
        # "end_t": 0.9,
        # species oriented, if true, pick pathways ending with top_n species,
        #  if False, just top n pathway
        "spe_oriented": False,
        # condense species path, no reactions
        "species_path": False,
        # atom followed
        "atom_f": "HA3",
        "init_s": 60,
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
    fast_transitions = [
        # 1068    549     O2+npropyl=npropyloo
        # reactants       9       O2      60      npropyl products        78      npropyloo
        # 1069    -549    O2+npropyl=npropyloo
        {
            "rxn": [1068, 1069],
            "spe": {
                "H": [60, 78],
                "O": [78, 9],
                "C": [60, 78],
                "HA1": [60, 78],
                "HA2": [60, 78],
                "HA3": [60, 78],
                "HA4": [60, 78]
            }
        },

        # 1096    565     O2+ipropyl=ipropyloo
        # reactants       9       O2      61      ipropyl products        80      ipropyloo
        # 1097    -565    O2+ipropyl=ipropyloo
        {
            "rxn": [1096, 1097],
            "spe": {
                "H": [61, 80],
                "O": [80, 9],
                "C": [61, 80],
                "HA1": [61, 80],
                "HA2": [61, 80],
                "HA3": [61, 80],
                "HA4": [61, 80]
            }
        },

        # 1116    575     O2+QOOH_1=well_1
        # reactants       9       O2      87      QOOH_1  products        90      well_1
        # 1117    -575    O2+QOOH_1=well_1
        {
            "rxn": [1116, 1117],
            "spe": {
                "H": [87, 90],
                "O": [90, 9],
                "C": [87, 90],
                "HA1": [87, 90],
                "HA2": [87, 90],
                "HA3": [87, 90],
                "HA4": [87, 90]
            }
        },

        # 1080    556     npropyloo=QOOH_1        557     npropyloo=QOOH_1
        # reactants       78      npropyloo       products        87      QOOH_1
        # 1081    -556    npropyloo=QOOH_1        -557    npropyloo=QOOH_1
        {
            "rxn": [1080, 1081],
            "spe": {
                "H": [78, 87],
                "C": [78, 87],
                "HA1": [78, 87],
                "HA2": [78, 87],
                "HA3": [78, 87],
                "HA4": [78, 87]
            }
        },

        # 1124	579	O2 + QOOH_2 = well_2
        # reactants	9	O2	88	QOOH_2	products	91	well_2
        # net_reactants	9	O2	88	QOOH_2	net_products	91	well_2
        {
            "rxn": [1124, 1125],
            "spe": {
                "H": [88, 91],
                "C": [88, 91],
                "HA1": [88, 91],
                "HA2": [88, 91],
                "HA3": [88, 91],
                "HA4": [88, 91]
            }
        },

        # 1146	590	O2 + QOOH_3 = well_3
        # reactants	9	O2	89	QOOH_3	products	92	well_3
        # net_reactants	9	O2	89	QOOH_3	net_products	92	well_3
        {
            "rxn": [1146, 1147],
            "spe": {
                "H": [89, 92],
                "C": [89, 92],
                "HA1": [89, 92],
                "HA2": [89, 92],
                "HA3": [89, 92],
                "HA4": [89, 92]
            }
        },

        # 1214	624	prod_1=frag_1+OH
        # reactants	94	prod_1	products	10	OH	101	frag_1
        # net_reactants	94	prod_1	net_products	10	OH	101	frag_1
        {
            "rxn": [1214, 1215],
            "spe": {
                "H": [94, 101],
                "C": [94, 101],
                "HA1": [94, 101],
                "HA2": [94, 101],
                "HA3": [94, 101],
                "HA4": [94, 101]
            }
        },

        # 1042    536     allyloxy=vinoxylmethyl
        # reactants       72      allyloxy        products        108     vinoxylmethyl
        # 1043    -536    allyloxy=vinoxylmethyl
        {
            "rxn": [1042, 1043],
            "spe": {
                "H": [72, 108],
                "C": [72, 108],
                "O": [72, 108],
                "HA1": [72, 108],
                "HA2": [72, 108],
                "HA3": [72, 108],
                "HA4": [72, 108]
            }
        },

        # 348     180     C2H5+O2=CH3CH2OO
        # reactants       39      C2H5    9       O2      products        50      CH3CH2OO
        # 349     -180    C2H5+O2=CH3CH2OO
        {
            "rxn": [348, 349],
            "spe": {
                "H": [39, 50],
                "C": [39, 50],
                "O": [9, 50],
                "HA1": [39, 50],
                "HA2": [39, 50],
                "HA3": [39, 50],
                "HA4": [39, 50]
            }
        },

        # 132     69      CH3+O2(+M)=CH3OO(+M)
        # reactants       25      CH3     9       O2      products        27      CH3OO
        # 133     -69     CH3+O2(+M)=CH3OO(+M)
        {
            "rxn": [132, 133],
            "spe": {
                "H": [25, 27],
                "C": [25, 27],
                "O": [25, 27],
                "HA1": [25, 27],
                "HA2": [25, 27],
                "HA3": [25, 27],
                "HA4": [25, 27]
            }
        }

        # # 586     300     O2C2H4OH=CH2CH2OH+O2
        # # reactants       85      O2C2H4OH        products        54      CH2CH2OH        9       O2
        # # 587     -300    O2C2H4OH=CH2CH2OH+O2
        # {
        #     "rxn": [586, 587],
        #     "spe": {
        #         "C": [85, 54]
        #     }
        # },

        # # 434     224     acetyl+O2=acetylperoxy
        # # reactants       9       O2      45      acetyl  products        47      acetylperoxy
        # # 435     -224    acetyl+O2=acetylperoxy
        # {
        #     "rxn": [434, 435],
        #     "spe": {
        #         "C": [45, 47]
        #     }
        # }

    ]

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
