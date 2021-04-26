#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#jhughey
#Liu Group
"""
Created on Mon Apr  1 13:58:05 2019

@author: jordanhughey

This script matches snps with 1st position in contact matrix

2nd position will be another script for target gene
"""

import csv
import argparse
parser = argparse.ArgumentParser(description='matches enhancer to contact matrix based off index')
parser.add_argument('--contact_ref_file', action='store', dest='contact_ref', default=None)
parser.add_argument('--enhancer_file', action='store', dest='enh_file', default=None)
parser.add_argument('--out_file', action='store', dest='outf', default=None)
parser.add_argument('--resolution', action='store', dest='res', default=None)
parser.add_argument('--chr', action='store', dest='chrom', default=None)
args = parser.parse_args()

#contact_ref = "IMR90_hic_mat_chr22_ref.txt"
contact_ref = args.contact_ref
#enh_file = "IMR90_ABC_gm.txt"
enh_file = args.enh_file
res = int(args.res)
chrom = 'chr' + args.chrom


#initialize ref dict
#index is key and contact score is value
#one key can have multiple contact values (list of lists form)
ref_dict = {}

#dictionary of all non-zero contact loci is created to match with snps of interest
with open(contact_ref) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            #if line_list[0] == 'locus1':
                #continue
            #i = ((int(line_list[1]) - 1) / res) + 1
            i = (int(line_list[1]) / res) + 1
            #j = ((int(line_list[2]) - 1) / res) + 1
            j = (int(line_list[2]) / res) + 1
            if i not in ref_dict:
                ref_dict[i] = [line_list]
            else:
                ref_dict[i].append(line_list)
            if j not in ref_dict:
                ref_dict[j] = [[line_list[0], line_list[2], line_list[1], line_list[3], line_list[4]]]
            else:
                ref_dict[j].append([line_list[0], line_list[2], line_list[1], line_list[3], line_list[4]])

print len(ref_dict)
#print ref_dict[404]

outf = args.outf + '_enh_match.' + chrom
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")

#This block goes through enhancer(regulatory snp) file line by line
#based off if snp overlaps with Hi-C locus
#contact locus is added with contact score
with open(enh_file) as f:
        inlines = f.readlines()
        for line in inlines:
            line_list = line.strip().split("\t")
            """if line_list[0] == 'Geneid':
                out_list = line_list[:4] + [line_list[6]] + ["locus1","locus2","contact_score"]
                out.writerow(out_list)"""
            if line_list[0] == chrom:
		#middle of enhancer
                snp_pos = int(line_list[1])
		#index based off resolution
                ind = (snp_pos/res) + 1
                if ind in ref_dict:
                    for item in ref_dict[ind]:
                        out_list = line_list + item
                        out.writerow(out_list)
                    
outfile.close()
