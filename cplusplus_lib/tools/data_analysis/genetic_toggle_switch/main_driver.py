import os
import sys
import job_drivers
import concatenate_results as cr
import time

if __name__ == '__main__':
    time_i = time.time()

    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))

    counter = 0 
    totalN = 25 

    while counter < totalN:
        # run jobs
        try:
            job_drivers.path_make_run(file_dir)
        except RuntimeError:
            print("RuntimeError")
            exit()

        # make a figure
        job_drivers.make_a_figure(file_dir, counter)

        # organize files. concatenate into one file
        cr.concatenate_time(file_dir, counter)
        cr.concatenate_concentration(file_dir, counter)
        cr.concatenate_temperature(file_dir, counter)
        cr.concatenate_pressure(file_dir, counter)

        counter += 1

    # send email
    job_drivers.send_email(file_dir)

    time_e = time.time()
    print("running time:\t" + str("{:.2f}".format((time_e - time_i)/3600.0)) + " hours\n")
