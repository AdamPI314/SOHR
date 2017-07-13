import os
import sys
import json


def read_configuration(fn):
    with open(fn, 'r') as fp1:
        s_d = json.load(fp1)
    return s_d


def write_configuration(sd, fn):
    with open(fn, 'w') as fp2:
        json.dump(sd, fp2, sort_keys=True, indent=4)

def read_iteration_Number(fn1):
    # fn1 = os.path.join(file_dir, "input", "setting.json")
    setting_d = read_configuration(fn1)
    return setting_d['SOHR_init']['iterationNumber']


if __name__ == '__main__':
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    # print(file_dir)

    fn1 = os.path.join(file_dir, "input", "setting.json")
    # print(fn1)
    setting_d = read_configuration(fn1)

    print setting_d['chem_init']['species_index_concentration']
    print setting_d['chem_init']['init_temperature']
    print setting_d['chem_init']['pressure_atm']

    setting_d['chem_init']['species_index_concentration']['0'] = 0.33
    setting_d['chem_init']['init_temperature'] = 1000
    setting_d['chem_init']['pressure_atm'] = 5

    fn2 = os.path.join(file_dir, "input", "test.json")
    write_configuration(setting_d, fn2)
