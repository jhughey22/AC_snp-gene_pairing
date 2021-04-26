#!/usr/bin/env bash
#Jordan Hughey
#Liu Lab

#This script gets the geometric mean (activity quantity) between dnase and histone rpms

usage() { echo "Usage: bash $0 [-a <H3K27ac_rpm_file>] [-b <dnase_rpm_file>] [ -c cell_type ] [ -o outfile ]" 1>&2; exit 1; }

while getopts "a:b:c:o:" opt; do
	case "$opt" in
	a)
		a=${OPTARG}
		;;
	b)
		b=${OPTARG}
		;;
	c)
		c=${OPTARG}
		;;
	o)
		o=${OPTARG}
		;;
	*)
		usage
		;;
	esac
done


histone_rpm=${a}
dnase_rpm=${b}
cell=${c}
outfile=${o%.txt*}

#cuts dnase rpm to add to histone file
cat ${dnase_rpm} | cut -f 6 > ${outfile}_dnase_rpm_for_paste.txt
paste -d '\t' ${histone_rpm} ${outfile}_dnase_rpm_for_paste.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$6,$9}' > ${outfile}_rpm_combined.txt

#now get geometric mena with python script
python get_geometric_mean_nofilter_2r.py --rpm_file ${outfile}_rpm_combined.txt --out_file ${outfile}_gm.txt

mv ${outfile}_dnase_rpm_for_paste.txt ${outfile}_rpm_combined.txt intermediate
