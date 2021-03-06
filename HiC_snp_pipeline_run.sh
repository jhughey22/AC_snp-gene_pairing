#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#First blocks to get opts from command line
ARGUMENT_LIST=(
    "resolution"
    "resolution_short"
    "tissue"
    "window"
    "cutoff"
)

# read arguments
opts=$(getopt \
    --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
    --name "$(basename "$0")" \
    --options "" \
    -- "$@"
)

eval set --$opts

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resolution)
            resolution=$2
            shift 2
            ;;

        --resolution_short)
            resolution_short=$2
            shift 2
            ;;

        --tissue)
            tissue=$2
            shift 2
            ;;

        --window)
            window=$2
            shift 2
            ;;

        --cutoff)
            cutoff=$2
            shift 2
            ;;

        *)
            break
            ;;
    esac
done

#need to add input for get opts
res=${resolution}
res_short=${resolution_short}
tiss=${tissue}
window=${window}
cutoff=${cutoff}

#make directories to store intermediate results
mkdir intermediate_files
mkdir intermediate_files/log_files


#Script 1 is to add hic positions to matrix per chromosome
mkdir intermediate_files/final_mat_files
hic=${tiss}
#hic=SX1

for i in `seq 20 22` 
do
        echo "starting chr${i}" >> ./intermediate_files/log_files/${hic}_add_pos.log
        hic_file=${hic}.nor.chr${i}.mat
        #out=${tiss}.40kb.mat
	out=${tiss}.${res_short}.mat
        #paste the index col to the matrix
        paste -d ' ' ${res_short}_bin_positions/chr${i}_bin_positions_only.txt ./contact_maps/${hic_file} | tr ' ' '\t' > ./intermediate_files/final_mat_files/${out}.intermediate.chr${i}
        #prints first index row to final mat file
        cat ${res_short}_bin_positions/chr${i}_bin_positions_colheader.txt > ./intermediate_files/final_mat_files/${out}.final.chr${i}
        #prints the above ind col and matrix to final file
        cat ./intermediate_files/final_mat_files/${out}.intermediate.chr${i} >> ./intermediate_files/final_mat_files/${out}.final.chr${i}
        echo "Finished with chr${i}" >> ./intermediate_files/log_files/${hic}_add_pos.log
done


#Convert to sparse triplet format
outdir=intermediate_files/${tiss}_sparse_triplet_format
mkdir ${outdir}
matrix=${tiss}.${res_short}.mat.final.chr
out=${tiss}_hic_mat_chr

parallel python2 ./scripts/convert_sparse_triplet_format.py --matrix ./intermediate_files/final_mat_files/${matrix}{} --out_file ${outdir}/${out}{} >> ./intermediate_files/log_files/${out}_all_convert_triplet.log ::: {20..22}


#Limit sparse triplet format to 1mb dist contacts

mkdir intermediate_files/${tiss}_sparse_triplet_format_${window}
#spt=${tiss}_sparse_triplet_pre/chr
spt=intermediate_files/${tiss}_sparse_triplet_format/${tiss}_hic_mat_chr
out=intermediate_files/${tiss}_sparse_triplet_format_${window}/${tiss}_${res_short}_sparse_triplet_${window}_chr

doawk () {
	chr=$1
	spt=$2
	out=$3
	cutoff=$4
	cat ${spt}${chr}_sparse_triplet.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$2-$3,$4}' | awk -F '\t' -v OFS='\t' -v cutoff="${cutoff}" '$4 <= cutoff {print}' > ${out}${chr}.txt
}
export -f doawk
parallel doawk {} ${spt} ${out} ${cutoff} ::: {20..22}

#Snp to bin matching for both SCZ snps and genotypes for PANTS
#need to make sure snps have been preprocessed for this

outdir=intermediate_files/${tiss}_snp_match_enhancers_${window}
mkdir ${outdir}
out_file=${outdir}/${tiss}_${window}_new_ABC
enh_file=./snp_files/intergenic_snps_beds/GTEx_V7_snps_intergenic_chr
contact_ref=intermediate_files/${tiss}_sparse_triplet_format_${window}/${tiss}_${res_short}_sparse_triplet_${window}_chr


parallel python2 ./scripts/match_enhancer_new_ABC.py --contact_ref_file ${contact_ref}{}.txt --enhancer_file ${enh_file}{}.sorted.bed --out_file ${out_file} --resolution ${res} --chr {} ">>" ./intermediate_files/log_files/${tiss}_${window}_new_ABC_chr{}.log ::: {20..22} 

#Do the same with SZC snps as in above block

outdir=intermediate_files/SCZ_snp_match_enhancers_${window}
mkdir ${outdir}
out_file=${outdir}/SCZ2_${window}_new_ABC
enh_file=./snp_files/scz_intergenic_snps_beds/scz2_snps_intergenic_chr
contact_ref=intermediate_files/${tiss}_sparse_triplet_format_${window}/${tiss}_${res_short}_sparse_triplet_${window}_chr

parallel python2 ./scripts/match_enhancer_new_ABC.py --contact_ref_file ${contact_ref}{}.txt --enhancer_file ${enh_file}{}.bed --out_file ${out_file} --resolution ${res} --chr {} ">>" ./intermediate_files/log_files/SCZ_${window}_new_ABC_chr{}.log ::: {20..22}


#Clean Hi-c bins to run distribution fitting

mkdir intermediate_files/${tiss}_SCZ2_match_enhancers_${window}
hmag=intermediate_files/${tiss}_snp_match_enhancers_${window}/${tiss}_${window}_new_ABC_enh_match
#echo $hmag
scz=intermediate_files/SCZ_snp_match_enhancers_${window}/SCZ2_${window}_new_ABC_enh_match
hmag_out=intermediate_files/${tiss}_SCZ2_match_enhancers_${window}/${tiss}_enh_match_${window}
scz_out=intermediate_files/${tiss}_SCZ2_match_enhancers_${window}/SCZ_enh_match_${window}

for i in `seq 20 22`
do
	echo "Chrom     hic_locus1      hic_locus2      dist    contact" > ${hmag_out}.chr${i}
	echo "Chrom     hic_locus1      hic_locus2      dist    contact" > ${scz_out}.chr${i}
done

parallel "cat ${hmag}.chr{} | cut -f 1,6- | sort | uniq >> ${hmag_out}.chr{}" ::: {20..22}
parallel "cat ${scz}.chr{} | cut -f 1,7- | sort | uniq >> ${scz_out}.chr{}" ::: {20..22}

#Fit Hi-C contact distributions and call significant contacts

mkdir intermediate_files/${tiss}_SCZ2_match_enhancers_significant_${window}
hmag=intermediate_files/${tiss}_SCZ2_match_enhancers_${window}/${tiss}_enh_match_${window}
scz=intermediate_files/${tiss}_SCZ2_match_enhancers_${window}/SCZ_enh_match_${window}
out=intermediate_files/${tiss}_SCZ2_match_enhancers_significant_${window}/${tiss}_SCZ2_enh_match_significant_${window}

parallel Rscript ./scripts/fit_contact_dist_tfinal.R ${hmag}.chr{} ${scz}.chr{} ${out}.chr{} ::: {20..22}


#make significant contacts into bed format

mkdir intermediate_files/${tiss}_SCZ2_match_enhancers_beds_${window}
spt=intermediate_files/${tiss}_SCZ2_match_enhancers_significant_${window}/${tiss}_SCZ2_enh_match_significant_${window}
out=intermediate_files/${tiss}_SCZ2_match_enhancers_beds_${window}/${tiss}_SCZ2_enh_match_significant_${window}

doawk () {
        chr=$1
        spt=$2
        out=$3
	res=$4
        cat ${spt}.chr${chr} | grep -v 'Chrom' | awk -F '\t' -v OFS='\t' -v res="${res}" 'res_final=res-1 {print $1,$2,$2+res_final,$0}' | sort -k1,1 -k2,2n > ${out}.chr${chr}.bed
}
export -f doawk

parallel doawk {} ${spt} ${out} ${res} ::: {20..22}


#run bedtools intersect to match snp with hi-c bin

mkdir intermediate_files/${tiss}_snp_intersect_${window}
sig=intermediate_files/${tiss}_SCZ2_match_enhancers_beds_${window}/${tiss}_SCZ2_enh_match_significant_${window}
input=./snp_files/intergenic_snps_beds/GTEx_V7_snps_intergenic
out=intermediate_files/${tiss}_snp_intersect_${window}/${tiss}_sig_contacts_intersect_snp_${window}

parallel "bedtools intersect -a ${sig}.chr{}.bed -b ${input}_chr{}.sorted.bed -wa -wb -sorted | cut -f 4-12,14-15 > ${out}.chr{}" ::: {20..22}


#make snp intersect into bed format

mkdir intermediate_files/${tiss}_snp_intersect_beds_${window}
spt=intermediate_files/${tiss}_snp_intersect_${window}/${tiss}_sig_contacts_intersect_snp_${window}
out=intermediate_files/${tiss}_snp_intersect_beds_${window}/${tiss}_sig_contacts_intersect_snp_${window}

doawk () {
        chr=$1
        spt=$2
        out=$3
	res=$4
        cat ${spt}.chr${chr} | awk -F '\t' -v OFS='\t' -v res="${res}" 'res_final=res-1 {print $1,$3,$3+res_final,$0}' | sort -k1,1 -k2,2n > ${out}.chr${chr}.bed
	#cat ${spt}.chr${chr} | grep -v 'Chrom' | awk -F '\t' -v OFS='\t' '{print $1,$2,$2+39999,$0}' | sort -k1,1 -k2,2n > ${out}.chr${chr}.bed
}
export -f doawk

parallel doawk {} ${spt} ${out} ${res} ::: {20..22}


#run bedtools intersect to match gene with hi-c bin

mkdir intermediate_files/${tiss}_gene_intersect_${window}
sig=intermediate_files/${tiss}_snp_intersect_beds_${window}/${tiss}_sig_contacts_intersect_snp_${window}
out=intermediate_files/${tiss}_gene_intersect_${window}/${tiss}_sig_contacts_intersect_snp_gene_${window}

#parallel "bedtools intersect -a ${sig}.chr{}.bed -b ${input}_chr{}.sorted.bed -wa -wb -sorted | cut -f 4-12,14-15 > ${out}.chr{}" ::: {1..22}
parallel "bedtools intersect -a ${sig}.chr{}.bed -b ./snp_files/gencode.v26lift37.annotation.promoter.exon.bed -wa -wb -sorted | cut -f 4- > ${out}.chr{}" ::: {20..22}




