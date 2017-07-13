import os
import sys
import numpy as np
import read_write_configuration as rwc


def concatenate_time(file_dir, p_iter_N=3):
    setting_d = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    time_iteration_number = setting_d['SOHR_init']['timeIterationNumber']
    if time_iteration_number == 1:
        time = np.loadtxt(os.path.join(file_dir, 'output', 'time_SOHR_fraction_' + str(p_iter_N) + '.csv'))
        np.savetxt(os.path.join(file_dir, 'output', 'time_SOHR_fraction_all.csv'), time)
    else:
        time_o = np.loadtxt(os.path.join(file_dir, 'output', 'time_SOHR_fraction_all.csv'))
        time_n = np.loadtxt(os.path.join(file_dir, 'output', 'time_SOHR_fraction_' + str(p_iter_N) + '.csv'))
        time_n += time_o[-1]

        with open(os.path.join(file_dir, 'output', 'time_SOHR_fraction_all.csv'), 'ab') as fp:
            np.savetxt(fp, time_n[1::])


def concatenate_concentration(file_dir, p_iter_N=3):
    setting_d = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    time_iteration_number = setting_d['SOHR_init']['timeIterationNumber']
    for p_iter_i in range(1, p_iter_N + 1):
        concentration = np.loadtxt(os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_' + str(p_iter_i) + '.csv'),
                delimiter=',')
        if time_iteration_number == 1:
            np.savetxt(os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_all_' + str(p_iter_i) + '.csv'), concentration,
                    delimiter=',')
            # last iteration, save one more copy as all
            if p_iter_i == p_iter_N: 
                np.savetxt(os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_all.csv'), concentration,
                        delimiter=',')
        else:
            with open(os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_all_' + str(p_iter_i) + '.csv'), 'ab') as fp:
                np.savetxt(fp, concentration[1::, :], delimiter=',')
            # last iteration, save one more copy as all
            if p_iter_i == p_iter_N: 
                with open(os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_all.csv'), 'ab') as fp:
                    np.savetxt(fp, concentration[1::, :], delimiter=',')



def concatenate_temperature(file_dir, p_iter_N=3):
    setting_d = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    time_iteration_number = setting_d['SOHR_init']['timeIterationNumber']
    if time_iteration_number == 1:
        temperature = np.loadtxt(os.path.join(file_dir, 'output', 'temperature_SOHR_fraction_' + str(p_iter_N) + '.csv'))
        np.savetxt(os.path.join(file_dir, 'output', 'temperature_SOHR_fraction_all.csv'), temperature)
    else:
        temperature = np.loadtxt(os.path.join(file_dir, 'output', 'temperature_SOHR_fraction_' + str(p_iter_N) + '.csv'))
        with open(os.path.join(file_dir, 'output', 'temperature_SOHR_fraction_all.csv'), 'ab') as fp:
            np.savetxt(fp, temperature[1::])


def concatenate_pressure(file_dir, p_iter_N=3):
    setting_d = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    time_iteration_number = setting_d['SOHR_init']['timeIterationNumber']
    if time_iteration_number == 1:
        pressure = np.loadtxt(os.path.join(file_dir, 'output', 'pressure_SOHR_fraction_' + str(p_iter_N) + '.csv'),
                delimiter=',')
        np.savetxt(os.path.join(file_dir, 'output', 'pressure_SOHR_fraction_all.csv'), pressure, delimiter=',')
    else:
        pressure = np.loadtxt(os.path.join(file_dir, 'output', 'pressure_SOHR_fraction_' + str(p_iter_N) + '.csv'),
                delimiter=',')
        with open(os.path.join(file_dir, 'output', 'pressure_SOHR_fraction_all.csv'), 'ab') as fp:
            np.savetxt(fp, pressure[1::, :], delimiter=',')
