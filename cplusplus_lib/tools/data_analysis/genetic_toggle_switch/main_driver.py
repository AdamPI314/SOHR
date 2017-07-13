import os
import sys
import update_settings as us
import job_drivers
import concatenate_results as cr
import read_write_configuration as rwc
import time

if __name__ == '__main__':
    time_i = time.time()

    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))

    time_iter_counter = 0 
    time_iter_N = 3 
    dt = 10.0 

    while time_iter_counter < time_iter_N:
        # update settings
        us.update_settings(file_dir, time_iter_counter, dt=dt)

        # run jobs
        job_drivers.delete_SOHR_temp_files(file_dir)
        try:
            job_drivers.path_ode_run(file_dir)
        except RuntimeError:
            print("RuntimeError")
            exit()

        p_iter_N = rwc.read_iteration_Number(os.path.join(file_dir, "input", "setting.json"))
        # organize files. concatenate into one file
        cr.concatenate_time(file_dir, p_iter_N=p_iter_N)
        cr.concatenate_concentration(file_dir, p_iter_N=p_iter_N)
        cr.concatenate_temperature(file_dir, p_iter_N=p_iter_N)
        cr.concatenate_pressure(file_dir, p_iter_N=p_iter_N)

        time_iter_counter += 1

    # run dlosde
    job_drivers.run_dlsode(file_dir)

    # make figures
    job_drivers.make_figures(file_dir)

    # send email
    job_drivers.send_email(file_dir)

    time_e = time.time()
    print("running time:\t" + str("{:.2f}".format((time_e - time_i)/3600.0)) + " hours\n")
