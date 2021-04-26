#!/bin/bash
#Jordan Hughey
#Liu Lab

#This script parses 10kb bin positions file for just the bin positions to add to Hi-C matrix

for i in `seq 1 22`;
do
	cat chr${i}_bin_positions.txt | cut -f 2 > chr${i}_bin_positions_only.txt
done
