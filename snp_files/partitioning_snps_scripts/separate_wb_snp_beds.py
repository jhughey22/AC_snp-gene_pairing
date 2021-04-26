#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 08:41:27 2020

@author: jordanhughey

This script takes SCZ snps breaks them into individual chrom beds

"""

import csv
import sys

#hmag='Adult_brain.genes.annot'
hmag = sys.argv[1]
#snp_annot='GTEx_V7_snp_annot_final_allchrs.bed'
#now dbsnp file
#chro = id_file[id_file.index(".chr")+1:]
out_file = sys.argv[2]
chro = out_file[out_file.index("_chr")+1:out_file.index(".bed")]
chrom = out_file[out_file.index("_chr")+4:out_file.index(".bed")]
print "Processing %s"%chro
print "Processing %s"%chrom


outf = out_file
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")

with open(hmag) as f:
    inlines = f.readlines()
    #loop through inlines
    for line in inlines:
        #split each line into list based off tab delimeter
        line_list = line.strip().split("\t")
        if line_list[0] == chrom:
            chr_final = 'chr' + line_list[0]
            #out_list = [chr_final, line_list[1], int(line_list[1])+1, line_list[6]]
            out_list = [chr_final, line_list[1], int(line_list[1])+1, line_list[2]]
            out.writerow(out_list)
            
        
outfile.close()
