#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Tue Mar  5 15:35:34 2019

@author: jordanhughey

This script converts read counts data to reads per million
"""

import csv
#import sys

import argparse
parser = argparse.ArgumentParser(description='convert read counts to reads per million')
parser.add_argument('--counts_file', action='store', dest='infile', default=None)
parser.add_argument('--rep1_total', action='store', dest='rep1_tot', default=None)
parser.add_argument('--rep2_total', action='store', dest='rep2_tot', default=None)
parser.add_argument('--out_file', action='store', dest='outfile', default=None)
args = parser.parse_args()


#infile="K562_H3K27ac_counts.txt"
infile = args.infile
outf = args.outfile
#rep1_tot = 2818742
rep1_tot = args.rep1_tot
#rep2_tot = 1223781
rep2_tot = args.rep2_tot

comb_tot = int(rep1_tot) + int(rep2_tot)
#Number of total reads divided by 1,000,000
#This is the per million scaling factor
scale_fac = comb_tot/float(1000000)
print "scale factor = %s"%scale_fac

#out = infile[:-4] + "_rpm.txt"
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")

#This block loops through every line of histone counts file converts counts to rpm
with open(infile) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            if line_list[0].startswith("#"):
                continue
		#this block adds to header row
            elif line_list[0] == "Geneid":
                out_list = line_list[:5] + ["combined_rpm"] + line_list[6:]
                out.writerow(out_list)
		#this conditional converts to reads per million
            else:
                rpm = (int(line_list[6]) + int(line_list[7])) / scale_fac
                out_list = line_list[:5] + [rpm] + line_list[6:]
                out.writerow(out_list)
                
outfile.close()
            
