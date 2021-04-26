#!/usr/bin/env Rscript
#Jordan Hughey
#Liu Group

#This script fits background distributions for dnase and H3K27ac rpms
#uses distributions to call significant loci

suppressMessages(library(fitdistrplus))
suppressMessages(library(qvalue))
"%&%" <- function(a,b) paste(a, b, sep = "")

args<-commandArgs(TRUE)

#read in hi-c file
#hmag_file <- '/gpfs/group/dxl46/default/private/jordan/ABC_dir/AC_hub_test/activity_scripts/SX_gm_1mb/SX_activity_1mb_chr13_gm.txt'
hmag_file <- args[1]
#scz_file <- '/gpfs/group/dxl46/default/private/jordan/ABC_dir/AC_hub_test/activity_scripts/SX_gm_1mb/SCZ_activity_1mb_chr13_gm.txt'
scz_file <- args[2]


hic_hmag <- read.table(hmag_file, header = TRUE, row.names = 'Geneid', stringsAsFactors = FALSE)
hic_scz <- read.table(scz_file, header = TRUE, stringsAsFactors = FALSE)


out_file <- args[3]
cat(out_file %&% ' \n')
#out_file <- '/gpfs/group/dxl46/default/private/jordan/ABC_dir/AC_hub_test/SX_SCZ2_activity_significant_1mb.chr13'
header <- c(colnames(hic_hmag), 'his_pvalue', 'his_qvalue', 'dnase_pvalue', 'dnase_qvalue')

write(header, file = out_file, ncolumns = 11, sep = "\t")

#Filter out top 5% to not bias distribution
top_histone_95 <- quantile(hic_scz$H3K27ac_comb_rpm, 0.95)
hic_filt <- subset(hic_scz, hic_scz$H3K27ac_comb_rpm < top_histone_95 & hic_scz$H3K27ac_comb_rpm > 0)

#fit weibull distribution of all contacts of certain distance
fw <- fitdist(hic_filt$H3K27ac_comb_rpm, "weibull")

#will use these parameters to get pvalues
w_shape <- fw$estimate['shape']
w_scale <- fw$estimate['scale']

hic_hmag$his_pvalue <- pweibull(hic_hmag$H3K27ac_comb_rpm, shape = w_shape, scale = w_scale, lower.tail = FALSE)
print(nrow(subset(hic_hmag, hic_hmag$his_pvalue <= 0.05)))
print(nrow(subset(hic_hmag, hic_hmag$his_pvalue > 0.05)))
pvals <- hic_hmag$his_pvalue
#qvals <- qvalue(p = pvals)
qvals <- tryCatch(qvalue(p = pvals), error = function(cond) {message('Error'); message(geterrmessage()); list()})
if (length(qvals) > 0) {
  hic_hmag$his_qvalue <- qvals$qvalues
  print(nrow(subset(hic_hmag, hic_hmag$his_qvalue <= 0.05)))
  print(nrow(subset(hic_hmag, hic_hmag$his_qvalue > 0.05)))
} else {
  qvals <- qvalue(p = pvals, pi0 = 1)
  hic_hmag$his_qvalue <- qvals$qvalues
  print(nrow(subset(hic_hmag, hic_hmag$his_qvalue <= 0.05)))
  print(nrow(subset(hic_hmag, hic_hmag$his_qvalue > 0.05)))
}


############################################################################################
#now repeat with DNase counts

top_dnase_95 <- quantile(hic_scz$dnase_comb_rpm, 0.95)
hic_filt_dnase <- subset(hic_scz, hic_scz$dnase_comb_rpm < top_dnase_95 & hic_scz$dnase_comb_rpm > 0)

#fit weibull distribution of all contacts of certain distance
fw_dnase <- fitdist(hic_filt_dnase$dnase_comb_rpm, "weibull")

#will use these parameters to get pvalues
w_dnase_shape <- fw_dnase$estimate['shape']
w_dnase_scale <- fw_dnase$estimate['scale']

hic_hmag$dnase_pvalue <- pweibull(hic_hmag$dnase_comb_rpm, shape = w_dnase_shape, scale = w_dnase_scale, lower.tail = FALSE)
print(nrow(subset(hic_hmag, hic_hmag$dnase_pvalue <= 0.05)))
print(nrow(subset(hic_hmag, hic_hmag$dnase_pvalue > 0.05)))
pvals <- hic_hmag$dnase_pvalue

#qvals <- qvalue(p = pvals)
qvals <- tryCatch(qvalue(p = pvals), error = function(cond) {message('Error'); message(geterrmessage()); list()})
if (length(qvals) > 0) {
	hic_hmag$dnase_qvalue <- qvals$qvalues
	print(nrow(subset(hic_hmag, hic_hmag$dnase_qvalue <= 0.05)))
	print(nrow(subset(hic_hmag, hic_hmag$dnase_qvalue > 0.05)))
        hic_sig <- subset(hic_hmag, hic_hmag$his_qvalue <= 0.05 & hic_hmag$dnase_qvalue <= 0.05)
	write.table(hic_sig, file = out_file, append = TRUE, sep = "\t", row.names=FALSE, 
	            col.names=FALSE, quote=F)
} else {
  qvals <- qvalue(p = pvals, pi0 = 1)
  hic_hmag$dnase_qvalue <- qvals$qvalues
  print(nrow(subset(hic_hmag, hic_hmag$dnase_qvalue <= 0.05)))
  print(nrow(subset(hic_hmag, hic_hmag$dnase_qvalue > 0.05)))
  hic_sig <- subset(hic_hmag, hic_hmag$his_qvalue <= 0.05 & hic_hmag$dnase_qvalue <= 0.05)
  write.table(hic_sig, file = out_file, append = TRUE, sep = "\t", row.names=FALSE, 
              col.names=FALSE, quote=F)
}






