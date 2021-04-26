#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#This script reformats each dbsnp vcf file by chromosome
#infile=$1

for i in `seq 1 22` X
do
	cat ./full_snps_beds/GTEx_V7_snps_chr${i}.bed | sort -k1,1 -k2,2n > ./full_snps_beds/GTEx_V7_snps_chr${i}.sorted.bed
	bedtools intersect -a ./full_snps_beds/GTEx_V7_snps_chr${i}.sorted.bed -b gencode.v26lift37.annotation.promoter.exon.bed -wa -wb -sorted > GTEx_V7_snps_chr${i}_pro_exon_intersect.txt
	cat GTEx_V7_snps_chr${i}_pro_exon_intersect.txt | cut -f 4 | sort | uniq | wc -l >> pro_exon_intersect.log
	cat GTEx_V7_snps_chr${i}_pro_exon_intersect.txt | cut -f 4 | sort | uniq >> pro_exon_rsid_list.txt
done
