
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
        }
        # # 1080    556     npropyloo=QOOH_1        557     npropyloo=QOOH_1
        # # reactants       78      npropyloo       products        87      QOOH_1
        # # 1081    -556    npropyloo=QOOH_1        -557    npropyloo=QOOH_1
        # {
        #     "rxn": [1080, 1081],
        #     "spe": [78, 87]
        # }
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
