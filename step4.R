# FILE:  --------------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION:
#
# OPTIONS:  none
# REQUIREMENTS:  ggplot2, data.table
# BUGS: --
# NOTES: 
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  29.05.2020
# REVISION: 29.05.2020

# Libraries, colors, plotting themes ------------------------------------------
library(data.table)
library(DT)
library(ggpubr)
library(ggplot2)
library(shiny)

# plotting theme
mashaGgplot2Theme <- list(
  theme_classic(base_size = 18) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.text.x = element_text(colour = 'black', size = 8),
          axis.line.y = element_line(colour = 'black', size = 0.5,
                                     linetype ='solid'),
          axis.text.y = element_text(colour = 'black', size = 12),
          panel.grid.minor = element_line(colour = "grey", size = 0.5,
                                          linetype = 2),
          strip.background = element_blank())
)

# Inputs ----------------------------------------------------------------------
countTabPath <- '../../NXT0562_NF/countTables/readsCombined_NXT0562_nxid13251.csv'
infoTabPath <- '../../MEC1Mice_BRBseq_samples_design.csv'

# names of the QC rows in the count table
qcNames <- c("__alignment_not_unique", "__ambiguous", "__no_feature",
             "__not_aligned", "__too_low_aQual")
# coverage cutoff.
minCoverage <- 2 * 10^6

# read in count table produced by BRB-seq tools
countTab <- fread(file = countTabPath, header = T, stringsAsFactors = F)
countTab <- countTab[, colnames(countTab) != 'undetermined', with = F]
# read in table with information about samples
infoTab <- fread(file = infoTabPath, header = T)

infoTab[, ID := apply(infoTab, 1, 
                      function(x) paste(x['Color'], x['Ctr or Axin2'], 
                                        x['Experiment batch'],
                                        x['Mouse'],
                                        sep = '_'))]
setkey(infoTab, `Plate position BRBseq`)
# restrict count table only to samples from info tab
countTab <- countTab[, c("Gene_id","Gene_name", 
                         infoTab$`Plate position BRBseq`), with = F]

metaData <- copy(infoTab)

# Remove samples with low coverage ------------------------------------
# I put cut off on the amount of USABLE reads
usableReads <- colSums(countTab[!Gene_id %in% qcNames][, -2:-1])
# samples to keep
samplesToKeep <- names(usableReads)[which(usableReads > minCoverage)]
# inform, that some samples didn't pass the cutoff
if (!all(samplesToKeep %in% colnames(countTab))){
  notPass <- setdiff(colnames(countTab), samplesToKeep)
  msg <- paste('Removed', paste(notPass, collapse = ', '), 'samples because',
               'they did not have enough of usable reads')
  message(msg)
  countTab <- countTab[, c("Gene_id","Gene_name", samplesToKeep)]
}

# Remove not-expressed genes & check # of expressed genes ---------------------
notExpr <- rowSums(countTab[, -2:-1]) == 0
# in case of at least one gene not being expressed - we remove it and inform
if (any(notExpr)) {
  msg <- paste0('Removed ', sum(notExpr), ' not expressed genes. This is ',
                round(100 * sum(notExpr) / nrow(countTab), 2), 
                '% of all genes')
  message(msg)
  countTab <- countTab[!notExpr, ]
}

# plot number of expressed genes per sample
numbGenesInSamples <- apply(countTab[!Gene_id %in% qcNames, -2:-1], 2,
                            function(x) sum(x != 0))
names(numbGenesInSamples) <- infoTab[names(numbGenesInSamples)]$ID
numbGenesInSamples <- data.table(ID = names(numbGenesInSamples),
                                 numbGenesInSamples = numbGenesInSamples)
setkey(metaData, ID)
setkey(numbGenesInSamples, ID)
metaData <- merge(metaData, numbGenesInSamples, all.x = T)
ggplot(metaData, aes(x = ID, y = numbGenesInSamples, fill = Color)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Number of genes") + 
  mashaGgplot2Theme + facet_grid(. ~ `Experiment batch`, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

