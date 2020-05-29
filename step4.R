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

# Inputs ----------------------------------------------------------------------
countTabPath <- '../../NXT0562_NF/countTables/readsCombined_NXT0562_nxid13251.csv'
infoTabPath <- '../../MEC1Mice_BRBseq_samples_design.csv'

# AAAA ------
countTab <- fread(countTabPath, header = T, stringsAsFactors = F)
countTab <- countTab[, colnames(countTab) != 'undetermined', with = F]

infoTab <- fread(infoTabPath, header = T)

sortedSamples <- colnames(countTab)
sortedSamples <- sortedSamples[!sortedSamples %in% c('Gene_id', 'Gene_name')]

# Number of expressed genes in each sample ------------------------------------
numbExprGenes <- apply(countTab[, sortedSamples, with = F], 2, 
                       function(x) sum(x != 0))
