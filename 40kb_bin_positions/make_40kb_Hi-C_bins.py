#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:18:41 2020

@author: jordanhughey

This script makes bin positions for 10kb Hi-C

"""
import csv

chrom_sizes = 'hg19.chrom.sizes.txt'
res = 40000

size_dict = {}

with open(chrom_sizes) as f:
    inlines = f.readlines()
    #loop through inlines
    for line in inlines:
        #split each line into list based off tab delimeter
        line_list = line.strip().split("\t")
        size_dict[line_list[0]] = int(line_list[1])
        
print size_dict

"""for i in range(1,23):
    outf = 'chr' + str(i) + '_bin_positions.txt'
    outfile = open(outf, "wb")
    out = csv.writer(outfile, delimiter="\t")
    
    chrom = 'chr' + str(i)
    print chrom,size_dict[chrom]
    c_size = size_dict[chrom]
    bins = range(1, c_size, res)
    for j in range(len(bins)):
        if j < len(bins) - 1:
            out_list = [chrom, bins[j], bins[j]+res]
            out.writerow(out_list)
        elif j == len(bins) - 1:
            out_list = [chrom, bins[j], c_size]
            out.writerow(out_list)
        
    outfile.close()"""

'''for i in range(1,23):
    outf = 'chr' + str(i) + '_bin_positions_colheader.txt'
    outfile = open(outf, "wb")
    out = csv.writer(outfile, delimiter="\t")
    
    chrom = 'chr' + str(i)
    print chrom,size_dict[chrom]
    c_size = size_dict[chrom]
    bins = range(1, c_size, res)
    bins_final = [0] + bins
    out.writerow(bins_final)
        
    outfile.close()
'''

for i in ['X']:
    outf = 'chr' + str(i) + '_bin_positions.txt'
    outfile = open(outf, "wb")
    out = csv.writer(outfile, delimiter="\t")
    
    chrom = 'chr' + str(i)
    print chrom,size_dict[chrom]
    c_size = size_dict[chrom]
    bins = range(1, c_size, res)
    for j in range(len(bins)):
        if j < len(bins) - 1:
            out_list = [chrom, bins[j], bins[j]+res]
            out.writerow(out_list)
        elif j == len(bins) - 1:
            out_list = [chrom, bins[j], c_size]
            out.writerow(out_list)
        
    outfile.close()

for i in ['X']:
    outf = 'chr' + str(i) + '_bin_positions_colheader.txt'
    outfile = open(outf, "wb")
    out = csv.writer(outfile, delimiter="\t")
    
    chrom = 'chr' + str(i)
    print chrom,size_dict[chrom]
    c_size = size_dict[chrom]
    bins = range(1, c_size, res)
    bins_final = [0] + bins
    out.writerow(bins_final)
        
    outfile.close()



