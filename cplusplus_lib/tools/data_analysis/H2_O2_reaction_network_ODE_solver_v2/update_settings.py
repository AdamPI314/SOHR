import os
import sys
import read_write_configuration as rwc
import read_last_line as rll


def update_dlsode_settings(file_dir, time=None):
    # there will always be a current setting
    curr_settings = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    settings_0 = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting_0.json'))
    settings_0['SOHR_init']['timeIterationNumber'] = curr_settings['SOHR_init']['timeIterationNumber']

    if time is not None:
        time_n = time
    else:
        time_n = rll.read_time(os.path.join(file_dir, 'output', 'time_SOHR_fraction_all.csv'))

    settings_0['time']['critical_time'] = time_n
    settings_0['time']['max_time'] = time_n
    settings_0['time']['path_end_time'] = time_n

    settings_0['job']['job_type'] = "solve_ODEs_for_concentration_using_LSODE"
    rwc.write_configuration(settings_0, os.path.join(file_dir, 'input', 'setting.json'))


# time iteration number t_iter_N
def update_settings(file_dir, t_iter_N=None, dt=None):
    # there will always be a current setting
    curr_settings = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting.json'))
    if t_iter_N == 0:
        tin = 0
    else:
        tin = curr_settings['SOHR_init']['timeIterationNumber']

    # the zero time interval
    if tin == 0:
        settings_0 = rwc.read_configuration(os.path.join(file_dir, 'input', 'setting_0.json'))
        curr_settings = settings_0
        # other time intervals
    else:
        concentration = rll.read_concentration(
                os.path.join(file_dir, 'output', 'concentration_SOHR_fraction_all.csv'))
        curr_settings['chem_init']['species_index_concentration'] = concentration

        temperature = rll.read_temperature(
                os.path.join(file_dir, 'output', 'temperature_SOHR_fraction_all.csv'))
        curr_settings['chem_init']['init_temperature'] = temperature

        pressure = rll.read_pressure(
                os.path.join(file_dir, 'output', 'pressure_SOHR_fraction_all.csv'))
        curr_settings['chem_init']['pressure_atm'] = pressure

    curr_settings['SOHR_init']['timeIterationNumber'] += 1

    if dt is not None:
        curr_settings['time']['critical_time'] = dt
        curr_settings['time']['max_time'] = dt
        curr_settings['time']['path_end_time'] = dt

    rwc.write_configuration(curr_settings, os.path.join(file_dir, 'input', 'setting.json'))


if __name__ == '__main__':
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    update_settings(file_dir)
