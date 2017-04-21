setwd("/Users/merrittpolomsky/Documents/bds-files/chapter-09-working-with-range-data")
rm(list=ls())

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicRanges")

library(BiocInstaller)
biocLite("GenomicFeatures")
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
biocLite("IRanges")

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

library(IRanges)

library(rtracklayer)

dbsnp137 <- import('mm10_snp137_chr1_trunc.bed.gz')

collasped_exons <- reduce(exons(txdb), ignore.strand=TRUE)
chr1_collasped_exons <- collasped_exons[seqnames(collasped_exons) == "chr1"]


summary(width(dbsnp137))


dbsnp137_resized <- dbsnp137
zw_i <- width(dbsnp137_resized) == 0
dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)


hits <- findOverlaps(dbsnp137_resized, chr1_collasped_exons, ignore.strand=TRUE)

length(unique(queryHits(hits)))
length(unique(queryHits(hits)))/length(dbsnp137_resized)

var_counts <- countOverlaps(chr1_collasped_exons, dbsnp137_resized, ignore.strand=TRUE)

chr1_collasped_exons$num_vars <- var_counts

export.bed(chr1_collasped_exons, "chr1_collasped_exons.bed")
