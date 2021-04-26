#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:15:16 2020

@author: jordanhughey

This script takes signficant activity values positions and outputs 

bins or snp-gene pairs that overlap them
"""

import csv
import sys

#add resolution as input
activity=sys.argv[1]
contacts=sys.argv[2]
out_file = sys.argv[3]
res = int(sys.argv[4])

#This keeps the bin for the activity significant region
activity_dict = {}
with open(activity) as f:
    inlines = f.readlines()
    #loop through inlines
    for line in inlines:
        #split each line into list based off tab delimeter
        line_list = line.strip().split("\t")
        if line_list[0] == 'Geneid' or line_list[0] == 'Chr':
            continue
        pos = '%s_%s'%(line_list[2],line_list[3])
        if pos not in activity_dict:
            activity_dict[pos] = [float(line_list[4]), float(line_list[5])]
            
print len(activity_dict)


outf = out_file
outfile = open(outf, "wb")
out = csv.writer(outfile, delimiter="\t")

#We match the first locus of the Hi-C data with the activity signficant loci
#If they match then we output that row as it is contact and activity significant
with open(contacts) as f:
    inlines = f.readlines()
    #loop through inlines
    for line in inlines:
        #split each line into list based off tab delimeter
        line_list = line.strip().split("\t")
        pos1 = int(line_list[1]) + 1
        #pos2 = int(line_list[1]) + 40000
        pos2 = int(line_list[1]) + res
        pos = '%i_%i'%(pos1,pos2)
        if pos in activity_dict:
            #if activity_dict[pos][0] >= h_avg and activity_dict[pos][1] >= d_avg:
            #if activity_dict[pos][0] > 0.77 and activity_dict[pos][1] > 0.77:
            out_list = [line_list[9], line_list[13]]
            out.writerow(out_list)

outfile.close()

