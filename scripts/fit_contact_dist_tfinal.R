#!/usr/bin/env Rscript
#Jordan Hughey
#Liu Group

#This script fits distributions on background Hi-C frequencies 
#and uses distributions to call significant contacts

suppressMessages(library(fitdistrplus))
suppressMessages(library(qvalue))
"%&%" <- function(a,b) paste(a, b, sep = "")

args<-commandArgs(TRUE)

#read in hi-c file
hmag_file <- args[1]
scz_file <- args[2]
hic_hmag <- read.table(hmag_file, header = TRUE, stringsAsFactors = FALSE)
hic_scz <- read.table(scz_file, header = TRUE, stringsAsFactors = FALSE)

#distribution fitting done at each distance
distances <- sort(unique(hic_scz$dist))

out_file <- args[3]
header <- c('Chrom', 'hic_locus1', 'hic_locus2', 'dist', 'contact', 'c_pvalue', 'c_qvalue')
write(header, file = out_file, ncolumns = 7, sep = "\t")

#loop to fit distribution at each distance
#futher distances should have weaker signals
i <- 1
for (i in 1:length(distances)) {
  #print which part of loop were on
  cat(i, "/", length(distances), "\n")
  dis <- distances[i]
  cat(dis, '\n')
  #subset hic dataframe to only contacts of same distance
  hic_hmag_dist <- subset(hic_hmag, hic_hmag$dist == dis)
  hic_scz_dist <- subset(hic_scz, hic_scz$dist == dis)

  #normalize hic scores
  #hic_dist$c_norm <- as.vector(scale(hic_dist$contact))

  #now subset based off lower 95 percentile
  top_95 <- quantile(hic_scz_dist$contact, 0.95)
  hic_filt <- subset(hic_scz_dist, hic_scz_dist$contact < top_95 & hic_scz_dist$contact > 0)
  
  #fit weibull distribution of all contacts of certain distance
  fw <- fitdist(hic_filt$contact, "weibull")
  
  #will use these parameters to get pvalues
  w_shape <- fw$estimate['shape']
  w_scale <- fw$estimate['scale']
  
  hic_hmag_dist$c_pvalue <- pweibull(hic_hmag_dist$contact, shape = w_shape, scale = w_scale, lower.tail = FALSE)
  print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_pvalue <= 0.05)))
  print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_pvalue > 0.05)))
  pvals <- hic_hmag_dist$c_pvalue
  #qvals <- qvalue(p = pvals)
  qvals <- tryCatch(qvalue(p = pvals), error = function(cond) {message('tryCatch - Warning'); message(geterrmessage()); list()})
  if (length(qvals) > 0) {
    hic_hmag_dist$c_qvalue <- qvals$qvalues
    print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_qvalue <= 0.01)))
    print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_qvalue > 0.01)))
    hic_sig <- subset(hic_hmag_dist, hic_hmag_dist$c_qvalue <= 0.05)
    write.table(hic_sig, file = out_file, append = TRUE, sep = "\t", row.names=FALSE, 
                col.names=FALSE, quote=F)
  } else {
    qvals <- qvalue(p = pvals, pi0 = 1)
    hic_hmag_dist$c_qvalue <- qvals$qvalues
    print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_qvalue <= 0.01)))
    print(nrow(subset(hic_hmag_dist, hic_hmag_dist$c_qvalue > 0.01)))
    hic_sig <- subset(hic_hmag_dist, hic_hmag_dist$c_qvalue <= 0.05)
    write.table(hic_sig, file = out_file, append = TRUE, sep = "\t", row.names=FALSE, 
                col.names=FALSE, quote=F)
  }
}
