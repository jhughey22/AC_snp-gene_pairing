#!/bin/bash
#Jordan Hughey
#Liu Lab

#This script gets the activity quantification (RPM) for H3K27ac bam files
#The peaks SAF file have already been defined in another script

usage() { echo "Usage: bash $0 [-a <rep1 bam file>] [-p <dnase_peaks_file>] [-c <cell type>] [-o <outfile.txt>]" 1>&2; exit 1; }

while getopts "a:b:p:c:o:" opt; do
	case "$opt" in
	a)
		a=${OPTARG}
		;;
	c)
		c=${OPTARG}
		;;
	p)
		p=${OPTARG}
		;;
	o)
		o=${OPTARG}
		;;
	*)
		usage
		;;
	esac
done


rep1=${a}
peaks=${p}
cell=${c}
outfile=${o%.txt*}

#tool that counts reads based of feature positions
featureCounts -a ${peaks} -F SAF -o ${outfile}_H3K27ac_counts.txt ${rep1} > ${outfile}_activity_histone_processing.log 

#Get the totals for reads for rep1 and rep2 (input to python script below)
count_r1=`cat ${outfile}_H3K27ac_counts.txt.summary | grep -v 'Status' | awk 'total=total+$2 {print total}' | tail -1`
echo "rep1 reads = ${count_r1}" >> ${outfile}_activity_histone_processing.log


#Python script converts read counts from above to reads per million
python2 ./scripts/convert_rpm_1r.py --counts_file ${outfile}_H3K27ac_counts.txt --rep1_total ${count_r1} --out_file ${outfile}_histone_rpm.txt >> ${outfile}_activity_histone_processing.log

#place intermediate files in folder
mv ${outfile}_H3K27ac_counts.txt ${outfile}_H3K27ac_counts.txt.summary ${outfile}_activity_histone_processing.log activity_scripts/intermediate
mv ${outfile}_histone_rpm.txt activity_scripts/rpm_files
