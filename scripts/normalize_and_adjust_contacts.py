#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Tue Mar 19 12:52:47 2019

@author: jordanhughey

This script normalizes Hi-C contact matrices using 
normalization vectors and adds expected counts to quantify contact score
"""

import csv
import argparse
parser = argparse.ArgumentParser(description='normalizes and adds expected pseduocount to hic contact matrices')
parser.add_argument('--raw_file', action='store', dest='rawmat', default=None)
parser.add_argument('--norm_vec_file', action='store', dest='normvec', default=None)
parser.add_argument('--expected_vec', action='store', dest='expvec', default=None)
parser.add_argument('--out_file', action='store', dest='outf', default=None)
parser.add_argument('--resolution', action='store', dest='res', default=None)
args = parser.parse_args()

#rawmat = "chr22_5kb.RAWobserved"
rawmat = args.rawmat
#normvec = "chr22_5kb.KRnorm"
normvec = args.normvec
#expvec = "chr22_5kb.KRexpected"
expvec = args.expvec
#res = 5000
res = int(args.res)
chrom = rawmat[rawmat.index('E30/')+4:rawmat.index('_5kb')]

#initialize norm factor dict
fac_dict = {}

#this block loops through the normvec file and stores factors in a dictionary
#based of index
with open(normvec) as f:
        i=1
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            if line_list[0] == "NaN":
                fac_dict[i] = 1
                i+=1
            else:
                fac_dict[i] = float(line_list[0])
                i+=1
#test length of normalization factor dictionary
print len(fac_dict)

#initialize exp factor dict
exp_dict = {}
#does the same as above storing expected values in dictionary
#based off index
with open(expvec) as f:
        i=1
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            exp_dict[i] = float(line_list[0])
            i+=1

print len(exp_dict)

#opens file to write too
outf = args.outf
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")
#header for writing file
header = ["Chrom","locus1","locus2","raw_observed","dist","i","j","ifac","jfac","normalized_score","expected","contact_score"]
out.writerow(header)

#1mb expected value
exp_mb_index = (1000000/res) + 1

with open(rawmat) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            #gets index for i and j
            i = (int(line_list[0])/res) + 1
            j = (int(line_list[1])/res) + 1
            #gets factor from dict using index
            ifac = fac_dict[i]
            jfac = fac_dict[j]
            #norm score from normalization factor conversion
            norm_score = float(line_list[2])/(ifac*jfac)
            #get index for expected; basically distance/res
            d = ((int(line_list[1]) - int(line_list[0]))/res) + 1
            dist = abs(int(line_list[1]) - int(line_list[0]))
            #gets exp from dict using index
            if d <= exp_mb_index:
                ex = exp_dict[exp_mb_index]
            else:
                ex = exp_dict[d]
            contact = min(100,norm_score+ex)
            out_list = [chrom] + line_list + [dist,i,j,ifac,jfac,norm_score,ex,contact]
            out.writerow(out_list)

#finally close writing file
outfile.close()

                
