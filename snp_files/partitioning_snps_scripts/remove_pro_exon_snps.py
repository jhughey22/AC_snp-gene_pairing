#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 08:41:27 2020

@author: jordanhughey

This script takes hmag prefrontal cortex gene-snp target pairs and makes

a bed file with each snp position and its target gene as 5th field

"""

import csv
import sys

#hmag='Adult_brain.genes.annot'
hmag = sys.argv[2]
#snp_annot='GTEx_V7_snp_annot_final_allchrs.bed'
#now dbsnp file
id_file = sys.argv[1]
#id_file = 'pro_exon_rsid_list.txt'
#chro = id_file[id_file.index(".chr")+1:]
out_file = sys.argv[3]

#print "Processing %s"%chro


id_dict = {}
with open(id_file) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            rsid = line_list[0]
            if rsid not in id_dict:
            	id_dict[rsid] = 1
print len(id_dict) 
        
        

#out_file = 'Adult_brain.genes.short.ens.snps.bed'

outf = out_file
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")

with open(hmag) as f:
    inlines = f.readlines()
    #loop through inlines
    for line in inlines:
        #split each line into list based off tab delimeter
        line_list = line.strip().split("\t")
        rsid = line_list[3]
        if rsid in id_dict:
        	continue
        else:
        	out.writerow(line_list)         
        
outfile.close()
