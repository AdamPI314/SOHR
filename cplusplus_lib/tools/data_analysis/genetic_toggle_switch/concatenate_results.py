import numpy as np
import os


def concatenate_time(file_dir, counter=0):
    time = np.loadtxt(os.path.join(file_dir, 'output', 'time_ssa_number.csv'))
    if counter == 0:
        np.savetxt(os.path.join(file_dir, 'output', 'time_ssa_number_all.csv'), time)
    else:
        with open(os.path.join(file_dir, 'output', 'time_ssa_number_all.csv'), 'ab') as fp:
            np.savetxt(fp, time)


def concatenate_concentration(file_dir, counter=0):
    concentration = np.loadtxt(
        os.path.join(file_dir, 'output', 'concentration_ssa_number.csv'), delimiter=',')
    if counter == 0:
        np.savetxt(os.path.join(file_dir, 'output', 'concentration_ssa_number_all.csv'), concentration,
                   delimiter=',')
    else:
        with open(os.path.join(file_dir, 'output', 'concentration_ssa_number_all.csv'), 'ab') as fp:
            np.savetxt(fp, concentration, delimiter=',')


def concatenate_temperature(file_dir, counter=3):
    temperature = np.loadtxt(os.path.join(file_dir, 'output', 'temperature_ssa_number.csv'))
    if counter == 0:
        np.savetxt(os.path.join(file_dir, 'output', 'temperature_ssa_number_all.csv'), temperature)
    else:
        with open(os.path.join(file_dir, 'output', 'temperature_ssa_number_all.csv'), 'ab') as fp:
            np.savetxt(fp, temperature)


def concatenate_pressure(file_dir, counter=3):
    pressure = np.loadtxt(os.path.join(file_dir, 'output', 'pressure_ssa_number.csv'))
    if counter == 0:
        np.savetxt(os.path.join(file_dir, 'output', 'pressure_ssa_number_all.csv'), pressure)
    else:
        with open(os.path.join(file_dir, 'output', 'pressure_ssa_number_all.csv'), 'ab') as fp:
            np.savetxt(fp, pressure)

