import subprocess
import os
from plot_ssa_delta import plot_ssa


def delete_temp_files(file_dir):
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


def path_make_run(file_dir):
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

def make_a_figure(file_dir, counter):
    os.chdir(file_dir)
    plot_ssa(file_dir, counter)

def send_email(file_dir):
    os.chdir(file_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
           "-u", "RUNNING JOB", "-m", "JOB FINISHED." + "\n" + file_dir]

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
