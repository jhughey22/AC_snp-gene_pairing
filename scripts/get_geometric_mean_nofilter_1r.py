#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Thu Mar  7 14:16:23 2019

@author: jordanhughey

This script gets geometric mean histone and dnase RPM counts
"""

import argparse
parser = argparse.ArgumentParser(description='convert read counts to reads per million')
parser.add_argument('--rpm_file', action='store', dest='infile', default=None)
parser.add_argument('--out_file', action='store', dest='outf', default=None)
args = parser.parse_args()

import csv
#import sys
import math

infile = args.infile
#infile = sys.argv[1]

out = args.outf
outfile = open(out, "wb")
out = csv.writer(outfile, delimiter="\t")

with open(infile) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            if line_list[0] == "Geneid":
                out_list = line_list[:4] + ["H3K27ac_comb_rpm","dnase_comb_rpm","geometric_mean"]
                out.writerow(out_list)
            else:
                gm = math.sqrt(float(line_list[4])*float(line_list[5]))
                out_list = line_list + [gm]
                out.writerow(out_list)
                
outfile.close()
