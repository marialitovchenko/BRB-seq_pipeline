#!/usr/bin/env Rscript

# FILE: barcodesPerLibrary.R --------------------------------------------------
#
# USAGE: Rscript barcodesPerLibrary.R path_to_sample_tab path_to_barcodes_tab
#
# DESCRIPTION: creates barcode files for the current sequencing runs. Replaces
#              plate coordinates with actual sample names.
#
# ARGUMENTS:  1) path to table with samples 2) path table with barcodes
#
# REQUIREMENTS: data.table
# BUGS: --
# NOTES: --
# AUTHOR:  Maria Litovchenko, maria.litovchenko@gmail.com
# COMPANY:  Alithea Genomics, Lausanne, Switzerland
# VERSION:  1
# CREATED:  05.07.2020
# REVISION: 05.07.2020
library(data.table)

sampleTabPath <- commandArgs(trailingOnly = T)[1]
barcodeTabPath <- commandArgs(trailingOnly = T)[2]

# read in sample table
sampleTab <- fread(sampleTabPath, header = T, stringsAsFactors = F)
sampleTab <- sampleTab[, c('pos', 'RunID', 'LibraryID', 'SampleID', 
                           'SampleName')]
setkey(sampleTab, pos)

# read in barcode table
barcodeTab <- fread(barcodeTabPath, header = T, stringsAsFactors = F)
colnames(barcodeTab) <- c('pos', 'B1')
setkey(barcodeTab, pos)

# assign barcodes to samples 
samplesBarcodes <- merge(sampleTab, barcodeTab, all.x = T)
setnames(samplesBarcodes, 'SampleName', 'Name')
samplesBarcodes[, ID := paste(RunID, LibraryID, SampleID, sep = '_')]
samplesBarcodes <- samplesBarcodes[order(RunID, LibraryID, SampleID, Name)]

for (seqRunLib in unique(samplesBarcodes$ID)) {
  write.table(data.frame(Name = character(), B1 = character()), 
              paste0(seqRunLib, '.txt'), append = F, 
              quote = F, sep = '\t', col.names = T, row.names = F)
}

for (i in 1:nrow(samplesBarcodes)){
  write.table(samplesBarcodes[, c('Name', 'B1')][i, ], 
              paste0(samplesBarcodes$ID[i], '.txt'), append = T, quote = F,
              sep = '\t', col.names = F, row.names = F)
}
