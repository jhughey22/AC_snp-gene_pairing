#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#This script reformats each dbsnp vcf file by chromosome
#infile=$1

wd=`pwd`

for i in `seq 1 22` X
do
	qsub -d ${wd} -N remove_pe_${i} -v rsid="pro_exon_rsid_list.txt",snp_annot="./full_snps_beds/GTEx_V7_snps_chr${i}.sorted.bed",out="GTEx_V7_snps_intergenic_chr${i}.sorted.bed" remove_pro_exon_snps.pbs
	#python /gpfs/group/dxl46/default/private/jordan/vcf_by_chr/add_rsid.py /gpfs/group/dxl46/default/private/jordan/vcf_by_chr/00-All_112315.vcf_nocomments.firstfields.chr${i} ${infile}
done
