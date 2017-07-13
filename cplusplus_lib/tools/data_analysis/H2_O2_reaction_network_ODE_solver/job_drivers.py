import subprocess
import os
import update_settings as us


def delete_non_dlsode_files(file_dir):
    os.chdir(file_dir)
    cmd = ["find", "./output", "-type", "f", "!", "-name", "*dlsode*", "-delete"]

    # Open/Create the output file
    outFile = open(os.path.join(file_dir, 'output', 'output_all.txt'), 'a+')
    errorFile = open(os.path.join(file_dir, 'output', 'error_all.txt'), 'a+')

    try:
        result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFile)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        outFile.write(out)

    outFile.close()
    errorFile.close()


def delete_srode_temp_files(file_dir):
    os.chdir(file_dir)
    cmd = ["find", "./output", "-type", "f", "!", "-name", "*all*", "-delete"]

    # Open/Create the output file
    outFile = open(os.path.join(file_dir, 'output', 'output_all.txt'), 'a+')
    errorFile = open(os.path.join(file_dir, 'output', 'error_all.txt'), 'a+')

    try:
        result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFile)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        outFile.write(out)

    outFile.close()
    errorFile.close()


def path_ode_run(file_dir):
    os.chdir(file_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    outFile = open(os.path.join(file_dir, 'output', 'output_all.txt'), 'a+')
    errorFile = open(os.path.join(file_dir, 'output', 'error_all.txt'), 'a+')

    try:
        result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFile)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        outFile.write(out)

    outFile.close()
    errorFile.close()


def run_dlsode(file_dir):
    os.chdir(file_dir)
    us.update_dlsode_settings(file_dir, time=None)
    path_ode_run(file_dir)


# http://stackoverflow.com/questions/3000724/running-matlab-in-the-background
def make_a_figure(file_dir, ind):
    os.chdir(file_dir)
    matlab_script_dir = os.path.join(file_dir, "tools/data_analysis/H2_O2_reaction_network_ODE_solver")
    matlab_script_filename = "plot_concentration_v2"
    matlab_cmd = "cd " + matlab_script_dir + "; " + matlab_script_filename + "(" + str(
            ind) + ")" + "; cd ../../..; exit;"
    cmd = ["nohup", "matlab", "-nosplash", "-nodisplay", "-r", matlab_cmd]

    # Open/Create the output file
    outFile = open(os.path.join(file_dir, 'output', 'output_all.txt'), 'a+')
    errorFile = open(os.path.join(file_dir, 'output', 'error_all.txt'), 'a+')

    try:
        result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFile)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        outFile.write(out)

    outFile.close()
    errorFile.close()


def make_figures(file_dir):
    for i in range(1, 9):
        make_a_figure(file_dir, i)


def send_email(file_dir):
    os.chdir(file_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
            "-u", "RUNNING JOB", "-m", "JOB FINISHED."+"\n"+file_dir,
            "-a", "./output/O2.png", "-a", "./output/H2O.png", "-a", "./output/H2.png", "-a", "./output/H2O2.png",
            "-a", "./output/H.png", "-a", "./output/OH.png", "-a", "./output/HO2.png", "-a", "./output/O.png"]

    # Open/Create the output file
    outFile = open(os.path.join(file_dir, 'output', 'output_all.txt'), 'a+')
    errorFile = open(os.path.join(file_dir, 'output', 'error_all.txt'), 'a+')

    try:
        result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFile)
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        outFile.write(out)

    outFile.close()
    errorFile.close()
