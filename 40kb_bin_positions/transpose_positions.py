#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Tue Jan  8 12:22:24 2019

@author: jordanhughey
"""

#This script will be put in a file to run over many files
#Turn positions into row instead of column
import csv
import sys
#infile = "./bin_positions/chr22_bin_positions_only.txt"
infile = sys.argv[1]

with open(infile) as f:
    lines = f.read().splitlines()
    
print len(lines)
print lines[1]
lines.insert(0,0)
print len(lines)
print lines[1]

outfile = open(infile[:-8] + "colheader.txt", "wb")
out = csv.writer(outfile, delimiter="\t")
out.writerow(lines)
outfile.close()
