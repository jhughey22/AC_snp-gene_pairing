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
    "dnase"
    "histone"
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

        --dnase)
            window=$2
            shift 2
            ;;

        --histone)
            cutoff=$2
            shift 2
            ;;

        *)
            break
            ;;
    esac
done


res=${resolution}
res_short=${resolution_short}
tiss=${tissue}
window=${window}
cutoff=${cutoff}
dnase=${dnase}
histone=${histone}

#Create saf format from gene intersect hic
mkdir intermediate_files/${tiss}_saf_files_${window}
sig=intermediate_files/${tiss}_gene_intersect_${window}/${tiss}_sig_contacts_intersect_snp_gene_${window}.chr
out=intermediate_files/${tiss}_saf_files_${window}/${tiss}_snp_gene_pairs_bins_${window}.chr

for i in `seq 20 22`
do
	echo "GeneID    Chr     Start   End     Strand" > ${out}${i}.saf
done

doawk () {
        chr=$1
        sig=$2
        out=$3
	res=$4
	cat ${sig}${chr} | awk -F '\t' -v OFS='\t' -v res="${res}" '{print $1,$2+1,$2+res,"."}' | sort | uniq | awk -F '\t' -v OFS='\t' 'id=id+1 {print id,$0}' >> ${out}${chr}.saf
        #cat ${spt}${chr}_5kb_GM12878_ref.txt | grep -v 'Chrom' | cut -f 1,2,3,5,10 | awk -F '\t' -v OFS='\t' -v cutoff="${cutoff}" '$4 <= cutoff {print}' > ${out}${chr}.txt
}
export -f doawk
parallel doawk {} ${sig} ${out} ${res} ::: {20..22}


#get H3K27ac and DNase counts
#mkdir activity_scripts
mkdir activity_scripts/${tiss}_gm_${window}
mkdir activity_scripts/rpm_files
mkdir activity_scripts/intermediate

saf=intermediate_files/${tiss}_saf_files_${window}/${tiss}_snp_gene_pairs_bins_${window}.chr
out=${tiss}_activity_${window}_chr
#hrep1=Spleen_ENCFF407QBM_H3K27ac.bam
#drep1=Spleen_ENCFF039NMI_DNase.bam
hrep1=${histone}
drep1=${dnase}

#first script gets rpms from histone bam files for saf positions
parallel bash ./activity_scripts/activity_histone_processing_1r.sh -a ${hrep1} -p ${saf}{}.saf -c ${tiss} -o ${out}{}.txt ::: {20..22}

#second script gets rpms from dnase bam files for saf positions
parallel bash ./activity_scripts/activity_dnase_processing_1r.sh -a ${drep1} -p ${saf}{}.saf -c ${tiss} -o ${out}{}.txt ::: {20..22}

#final script for activity gets geometric mean between histone and dnase rpm
#this is activity value
parallel bash ./activity_scripts/get_geometric_mean_nofilter_1r.sh -a activity_scripts/rpm_files/${out}{}_histone_rpm.txt -b activity_scripts/rpm_files/${out}{}_dnase_rpm.txt -c ${tiss} -o ${out}{}.txt -w ${window} ::: {20..22}


#Create saf format for scz background
mkdir intermediate_files/SCZ_saf_files_${window}
sig=intermediate_files/${tiss}_SCZ2_match_enhancers_${window}/SCZ_enh_match_${window}.chr
out=intermediate_files/SCZ_saf_files_${window}/SCZ_enh_match_${window}.chr

for i in `seq 20 22`
do
        echo "GeneID    Chr     Start   End     Strand" > ${out}${i}.saf
done

doawk () {
        chr=$1
        sig=$2
        out=$3
	res=$4
	cat ${sig}${chr} | awk -F '\t' -v OFS='\t' -v res="${res}" '{print $1,$2+1,$2+res,"."}' | grep -v "Chrom" | sort | uniq | awk -F '\t' -v OFS='\t' 'id=id+1 {print id,$0}' >> ${out}${chr}.saf
}
export -f doawk
parallel doawk {} ${sig} ${out} ${res} ::: {20..22}




#get H3K27ac and DNase counts for background

saf=intermediate_files/SCZ_saf_files_${window}/SCZ_enh_match_${window}.chr
out=SCZ_activity_${window}_chr

#first script gets rpms from histone bam files for saf positions
parallel bash ./activity_scripts/activity_histone_processing_1r.sh -a ${hrep1} -p ${saf}{}.saf -c ${tiss} -o ${out}{}.txt ::: {20..22}

#second script gets rpms from dnase bam files for saf positions
parallel bash ./activity_scripts/activity_dnase_processing_1r.sh -a ${drep1} -p ${saf}{}.saf -c ${tiss} -o ${out}{}.txt ::: {20..22}

#final script for activity gets geometric mean between histone and dnase rpm
#this is activity value
parallel bash ./activity_scripts/get_geometric_mean_nofilter_1r.sh -a activity_scripts/rpm_files/${out}{}_histone_rpm.txt -b activity_scripts/rpm_files/${out}{}_dnase_rpm.txt -c ${tiss} -o ${out}{}.txt -w ${window} ::: {20..22}


#now run fit distribution for each chrom
mkdir intermediate_files/${tiss}_activity_significant_${window}
mkdir intermediate_files/${tiss}_activity_filter_significant_${window}
hmag=./activity_scripts/${tiss}_gm_${window}/${tiss}_activity_${window}_chr
scz=./activity_scripts/${tiss}_gm_${window}/SCZ_activity_${window}_chr
out=intermediate_files/${tiss}_activity_significant_${window}/${tiss}_SCZ2_activity_significant_${window}.chr

parallel Rscript ./scripts/fit_activity_dist_tfinal.R ${hmag}{}_gm.txt ${scz}{}_gm.txt ${out}{} ::: {20..22}


#merge significant activity loci with snp-gene significant loci
act=intermediate_files/${tiss}_activity_significant_${window}/${tiss}_SCZ2_activity_significant_${window}.chr
out=intermediate_files/${tiss}_activity_filter_significant_${window}/${tiss}_new_act_filt_snp-genes_${window}.chr
con=intermediate_files/${tiss}_gene_intersect_${window}/${tiss}_sig_contacts_intersect_snp_gene_${window}.chr
parallel python2 ./scripts/activity_filtering_significant_matching.py ${act}{} ${con}{} ${out}{} ${res} ::: {20..22}

#Combine snp-gene pair lists
out=${tiss}_snp-gene_pairs_sig_act_filt_${window}.txt

for i in `seq 20 22`
do
        echo "Start chr${i}"
        f=./intermediate_files/${tiss}_activity_filter_significant_${window}/${tiss}_new_act_filt_snp-genes_${window}.chr${i}
        cat ${f} | sort | uniq >> ${out}
        f2=./snp_files/promoter_exon_snps_beds/GTEx_V7_snps_chr${i}_pro_exon_intersect.txt
        cat ${f2} | cut -f 4,8 | sort | uniq >> ${out}
        echo "Finished chr${i}"
done






