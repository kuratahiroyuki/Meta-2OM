#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--infile1', help='file')
    parser.add_argument('-i2', '--infile2', help='file')
    parser.add_argument('-o1', '--outfile1', help='file')
    parser.add_argument('-o2', '--outfile2', help='file')

    infile1 = parser.parse_args().infile1
    infile2 = parser.parse_args().infile2

    outfile1 = parser.parse_args().outfile1
    outfile2 = parser.parse_args().outfile2

    df1 = pd.read_csv(infile1, sep='\t', header=None ) #CV
    df2 = pd.read_csv(infile2, sep='\t', header=None) #test

    df1.iloc[ (df1[1] == -1), 1 ] = 0
    df2.iloc[ (df2[1] == -1), 1 ] = 0 

    df1.to_csv(outfile1, header=None, index=None)
    df2.to_csv(outfile2, header=None, index=None)


