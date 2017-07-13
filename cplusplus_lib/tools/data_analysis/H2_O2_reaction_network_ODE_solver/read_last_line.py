from collections import deque
import csv
import os
import sys


def get_last_row(csv_filename):
    with open(csv_filename, 'r') as f:
        try:
            last_row = deque(csv.reader(f), 1)[0]
        except IndexError:  # empty file
            last_row = None
        return last_row


def read_time(fn0):
    time = get_last_row(fn0)[0]
    return time


def read_concentration(fn1):
    lr1 = get_last_row(fn1)
    lr1_d = dict()
    for i in range(len(lr1)):
        lr1_d[unicode(i)] = lr1[i]
    return lr1_d


def read_temperature(fn2):
    lr2 = get_last_row(fn2)[0]
    return lr2


def read_pressure(fn3):
    lr3 = get_last_row(fn3)[0]
    return lr3


if __name__ == '__main__':
    file_dir = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir))
    fn1 = os.path.join(file_dir, 'output', 'concentration_srode_fraction_3.csv')

    lr1 = get_last_row(fn1)
    lr1_d = dict()
    for i in range(len(lr1)):
        lr1_d[unicode(i)] = lr1[i]
    print(lr1_d)

    fn2 = os.path.join(file_dir, 'output', 'temperature_srode_fraction_3.csv')
    lr2 = get_last_row(fn2)[0]
    print(lr2)

    fn3 = os.path.join(file_dir, 'output', 'pressure_srode_fraction_3.csv')
    lr3 = get_last_row(fn3)[0]
    print(lr3)
