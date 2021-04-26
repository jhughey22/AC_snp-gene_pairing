#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Tue Mar 26 15:25:07 2019

@author: jordanhughey

This script puts hic normalized matrix into sparse triplet format
"""
#set or parse inputs
import csv
import argparse
parser = argparse.ArgumentParser(description='processes contact matrix')
parser.add_argument('--matrix', action='store', dest='matrix', default=None)
parser.add_argument('--out_file', action='store', dest='outf', default=None)
args = parser.parse_args()

#matrix = "IMR90.nor.chr22.mat.final"
matrix = args.matrix
#outf = "IMR90_chr22_hic_mat"
outf = args.outf
chrom = matrix[matrix.index('chr'):]

#This block puts matrix into lower triangle format with diag
trip_name = outf + "_sparse_triplet.txt"
trip_file = open(trip_name, "wb")
out = csv.writer(trip_file, delimiter="\t")

#intialize index dict
dex_dict = {}

#Prints each loci pair and the contact frequency in the matrix
#Stops at the diagonal to only print lower triangle
i=1
with open(matrix) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            if line_list[0] == "0":
                for item in line_list:
                    dex_dict[line_list.index(item)] = item
            else:
                for item in range(1,i):
                    if float(line_list[item]) > 0:
                        #continue
                        out_list = [chrom, line_list[0],dex_dict[item],line_list[item]]
                        out.writerow(out_list)
                    else:
                        continue
            i+=1

trip_file.close()

print '%s  %s'%(chrom, len(dex_dict))
#print dex_dict[1]
#print dex_dict[1283]

