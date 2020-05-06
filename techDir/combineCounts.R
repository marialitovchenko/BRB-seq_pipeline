#!/usr/bin/env Rscript

# FILE: combineCounts.R -------------------------------------------------------
#
# USAGE: Called from the nextflow pipeline in order to merge counts from 
#        different demultiplex mapped files into one.
#        Takes 2 arguments: 
#        1) what to merge: UMI or reads
#        2) path to the table containing info about samples
#        3) output prefix
# Rscript --vanilla combineCounts.R reads path/to/count/folder this_are_reads
#
# DESCRIPTION: Merging counts from different demultiplex mapped files into one
#              with use of data.table. Names of the samples are inferred from 
#              the paths. All count files in the given folder are used.
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  17.04.2017
# REVISION: 17.04.2017

library(data.table)

# Functions -------------------------------------------------------------------
#' extractCountsColumn
#' Checks in the count file which column has the counts, checks if it's the 
#' same as given sample name, if not - uses column with the given sample name.
#' Informs if ALL columns in the file has 0 or if > 1 column in the file has 
#' counts.
#' @param filePath path to the counts file
#' @param sampleName string, sample name
#' @return data.table with column names Gene_id and name of the sample
#' @note I don't need to know which sample I'm considering as the only sample
#'       with non-zero counts will be the one which was submitted to countDGE
#'       from BRB-seq tools. However, if all samples failed and none of the 
#'       reads were assigned to the genes, I have a problem. So still request a
#'       sample name as an input.
extractCountsColumn <- function(filePath, sampleName) {
  oneSampleCounts <- fread(filePath, header = T, stringsAsFactors = F)
  setnames(oneSampleCounts, 'Unknown_Barcode', 'undetermined')
  maxCounts <- as.integer(apply(oneSampleCounts, 2, max))
  if (max(maxCounts, na.rm = T) != 0) {
    # here it's != 0 and not compared to maximum, in order to catch if only
    # one sample has assigned counts
    sampleIndex <- which(maxCounts != 0)
    if (length(sampleIndex) == 1) {
      guessName <- colnames(oneSampleCounts)[sampleIndex]
      if (guessName == sampleName) {
        result <- oneSampleCounts[, c(1:2, sampleIndex), with = F]
        setkey(result, Gene_id, Gene_name)
      } else {
        msg <- paste0("Given sample name (", sampleName, ") doesn't",
                      " correspond to the sample with the counts in ",
                      "the table (", colnames(oneSampleCounts)[sampleIndex],
                      "). Assigning given sample name (", sampleName, ")",
                      " and taking counts from corresponding column.")
        warning(msg)
        result <- oneSampleCounts[, c('Gene_id', 'Gene_name', sampleName)]
        setkey(result, Gene_id, Gene_name)
      }
    } else {
      msg <- paste0("Multiple samples have counts in this table, although",
                    sampleName, "is assigned as main sample. ",
                    "Assigning given sample name (", sampleName, ")",
                    " and taking counts from corresponding column.")
      warning(msg)
      result <- oneSampleCounts[, c('Gene_id', 'Gene_name', sampleName)]
      setkey(result, Gene_id, Gene_name)
    }
  } else {
    msg <- paste0("No reads falling within the genes was found for any",
                  "sample.")
    warning(msg)
    result <- oneSampleCounts[, c('Gene_id', 'Gene_name', sampleName)]
    setkey(result, Gene_id, Gene_name)
  }
  result
}

#' mergeCountsInto1Tab
#' Merges columns from individual count files into 1 count table using 
#' extractCountsColumn function
#' @param countsPaths vector of strings, paths to count files
#' @return data.table with columns Gene_id, Gene_name + as many columns as 
#'         files
mergeCountsInto1Tab <- function(countsPaths) {
  # extract names of the samples from the paths
  sampleNames <- sapply(countsPaths, function(x) gsub('[.].*', '', x))
  sampleNames <- sapply(sampleNames, function(x) gsub('.*[/]', '', x))
  
  countCols <- lapply(1:length(sampleNames), 
                      function(x) extractCountsColumn(countsPaths[[x]],
                                                      sampleNames[x]))
  countTab <- countCols[[1]]
  setkey(countTab, Gene_id, Gene_name)
  for (i in 2:length(countCols)) {
    countTab <- merge(countTab, countCols[[i]], all = T)
    setkey(countTab, Gene_id, Gene_name)
  }
  setcolorder(countTab, sort(colnames(countTab)))
  countTab
}

# Perform merging -------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
infoTabPath <- args[1]
mode <- args[2]
outputPref <- args[3]

# check that all arguments were submitted
if (is.na(infoTabPath)) {
  msg <- paste('[ERROR] in combineCounts.R:',
               'A path to the info table has to be submitted!')
  stop(msg)
} else {
  if (is.na(mode)) {
    msg <- paste('[ERROR] in combineCounts.R:',
                 'A mode of merging needs to be submitted!')
    stop(msg)
  } else {
    if (!mode %in% c('UMI', 'reads')) {
      msg <- paste('[ERROR] in combineCounts.R:',
                   'A mode should be one of "reads" or "UMI"')
      stop(msg)
    } else {
      if (is.na(outputPref)) {
        outputPref <- ''
        msg <- paste('[WARNING] in combineCounts.R:',
                     "Output prefix wasn't submitted,",
                     "name clash is possible!")
        warning(msg)
      }
    }
  }
}

# read-in info table containing info about runs, libraries, genomes, counts
infoTab <- fread(infoTabPath, header = F, stringsAsFactors = F)
if (nrow(infoTab) == 0) {
  msg <- paste('[ERROR] in combineCounts.R:',
               'An empty info table was submitted!')
  stop(msg)
}
infoNames <- c('Run', 'Library', 'Sample', 'Specie', 'Genome', 
               'Trimmed_Index_Fq', 'Trimmed_Data_Fq', 'Demultiplex_Fq',
               'Mapped', 'Map_stats', 'Counts_UMI', 'Counts_Reads')
setnames(infoTab, colnames(infoTab), infoNames)

# perform the merge
for (run in unique(infoTab$Run)) {
  for (lib in unique(infoTab[Run == run]$Library)) {
    sampsInLibRun <- infoTab[Run == run & Library == lib]
    mergedTab <- switch(mode, 
                        "reads" = mergeCountsInto1Tab(sampsInLibRun$Counts_UMI),
                        "UMI" = mergeCountsInto1Tab(sampsInLibRun$Counts_UMI))
    write.table(mergedTab, paste0(outputPref, '_', run, '_', lib, '.csv'),
                append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  }
}
