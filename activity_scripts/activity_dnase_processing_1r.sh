#!/bin/bash
#Jordan Hughey
#Liu Lab

#This script gets the activity quantification (RPM) for DNase bam files
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
featureCounts -a ${peaks} -F SAF -o ${outfile}_dnase_counts.txt ${rep1} > ${outfile}_activity_dnase_processing.log

#Get the totals for reads for rep1 and rep2 (input to python script below)
count_r1=`cat ${outfile}_dnase_counts.txt.summary | grep -v 'Status' | awk 'total=total+$2 {print total}' | tail -1`
echo "rep1 reads = ${count_r1}" >> ${outfile}_activity_dnase_processing.log

#Python script converts read counts from above to reads per million
python2 ./scripts/convert_rpm_1r.py --counts_file ${outfile}_dnase_counts.txt --rep1_total ${count_r1} --out_file ${outfile}_dnase_rpm.txt >> ${outfile}_activity_dnase_processing.log

#place intermediate files in folder
mv ${outfile}_dnase_counts.txt ${outfile}_dnase_counts.txt.summary ${outfile}_activity_dnase_processing.log activity_scripts/intermediate
mv ${outfile}_dnase_rpm.txt activity_scripts/rpm_files
