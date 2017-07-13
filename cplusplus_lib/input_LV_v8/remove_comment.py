#!/usr/bin/env python
import os
import sys
import argparse

if __name__== '__main__':
    parser= argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input filename")
    parser.add_argument("-o", "--output", help="output filename")
    args= parser.parse_args()

    print args.input
    print args.output

    with open(args.input) as f:
        with open(args.output, "w") as f2:
            for line in f:
                if not "//" in line:
                    f2.write(line)
