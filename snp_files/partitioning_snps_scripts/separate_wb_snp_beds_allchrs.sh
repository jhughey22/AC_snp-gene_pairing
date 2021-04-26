#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#This script reformats each dbsnp vcf file by chromosome
#infile=$1

wd=`pwd`

for i in `seq 1 22` X
do
	qsub -d ${wd} -N wb_beds_${i} -v snp_annot="GTEx_V7_snp_annot_final_allchrs.txt",out="GTEx_V7_snps_chr${i}.bed" separate_wb_snp_beds.pbs
done
