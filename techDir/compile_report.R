#!/usr/bin/env Rscript

# FILE: compile_report.R ------------------------------------------------------
#
# USAGE: Rscript compile_report.R [ARGS]
#
# DESCRIPTION: prepares input files needed for Generate_UserReport.Rmd and runs
#              report generation
#
# ARGUMENTS:  1) path to Rmarkdown file 2) user name 3) pi name 4) path to the 
#             folder with results from nextflow pipeline 5) Run ID 
#             6) Library ID 7) Sample ID 8) Specie 9) Genome. 
#             ATTENTION! None of the arguments should have spaces inside them
# REQUIREMENTS:  ggplot2, data.table, rmarkdown
# BUGS: --
# NOTES: --
# AUTHOR:  Maria Litovchenko, maria.litovchenko@gmail.com
# COMPANY:  Alithea Genomics, Lausanne, Switzerland
# VERSION:  1
# CREATED:  03.07.2020
# REVISION: 03.07.2020

# Unzip / get paths to files for markdown -------------------------------------
cmdArgs <- commandArgs(trailingOnly = T)
# exctract path to markdown script straight away
markdownScript <- cmdArgs[1]
cmdArgs <- cmdArgs[-1]

names(cmdArgs) <- c('User', 'PI', 'ResultFolder', 'RunID', 'LibraryID', 
                    'SampleID', 'Specie', 'Genome', 'submissionTab', 
                    'outputDirPath')

subFolders <- paste(cmdArgs['RunID'], cmdArgs['LibraryID'], 
                    cmdArgs['SampleID'], sep = '/') 

# get fastqc text report for R1 and R2: we need to unzip them first
fastqcDir <- paste(cmdArgs['ResultFolder'], 'fastQC', subFolders, sep = '/')
fastqcR1zip <- list.files(fastqcDir, pattern = '_R1_.*zip$', full.names = T) 
fastqcR2zip <- list.files(fastqcDir, pattern = '_R2_.*zip$', full.names = T) 
system(paste('unzip -d', fastqcDir, '-u', fastqcR1zip))
system(paste('unzip -d', fastqcDir, '-u', fastqcR2zip))
fastqcR1file <- list.files(gsub('.zip', '/', fastqcR1zip), 
                           pattern = 'fastqc_data.txt', full.names = T)
fastqcR2file <- list.files(gsub('.zip', '/', fastqcR2zip), 
                           pattern = 'fastqc_data.txt', full.names = T)

# get trimming report from trim_galore and fastqc for trimmed R2
trimDir <-  paste(cmdArgs['ResultFolder'], 'trimmed', subFolders, sep = '/')
trimRepR2 <- list.files(trimDir, pattern = '_trimming_report.txt$', 
                        full.names = T)
trimFastqcR2 <- list.files(trimDir, pattern = 'val_2_fastqc.zip$', 
                           full.names = T)
system(paste('unzip -d', trimDir, '-u', trimFastqcR2))
trimFastqcR2 <- list.files(gsub('.zip', '/', trimFastqcR2), 
                           pattern = 'fastqc_data.txt', full.names = T)  

# demultiplexing statistics
demultiplStatsFile <-  paste(cmdArgs['ResultFolder'], 'demultiplexed', 
                             subFolders, 'stats.txt', sep = '/')
# mapping statistics
mapStatsFile <- paste(cmdArgs['ResultFolder'], 'mapStatsTab.csv', sep = '/')
# count table
countTabFile <- paste0(cmdArgs['ResultFolder'], '/countTables/',
                       cmdArgs['RunID'], '_', cmdArgs['LibraryID'], '_',
                       cmdArgs['SampleID'], '_', cmdArgs['Genome'],
                       '_readsCombined.csv')

# Run Rmarkdown script to generate user report --------------------------------
inputArgs <- c(cmdArgs[c(-10, -3)], fastqcR1file, fastqcR2file, trimRepR2, 
               trimFastqcR2, demultiplStatsFile, mapStatsFile, countTabFile)
names(inputArgs) <- c(names(cmdArgs)[c(-10, -3)], 'fastqcR1', 'fastqcR2', 
                      'trimInfoR2', 'fastqcTrimR2', 'demultStats', 'mapStats',
                      'countTab')
resultHTML <- paste0(cmdArgs['RunID'], '_', cmdArgs['LibraryID'], '_',
                    cmdArgs['SampleID'], '_', cmdArgs['Genome'], '.html')
# for beginning, each html will be saved in subfolder: RunID/LibraryID/SampleID
# of the result folder. Otherwise, since Rmd first creates md file with the 
# same name as markdown (Generate_UserReport.md) then upon a parallel execution
# of same markdown file for different input parameters, clushes could occur and
# not all reports will be created. Subfolders solve this problem.
resultHTML <- paste(cmdArgs['outputDirPath'], subFolders, cmdArgs['Genome'], 
                    resultHTML, sep = '/')

# recursively create directories
system(paste('mkdir -p', dirname(resultHTML)))

system(paste0("R -e \"rmarkdown::render(\'", markdownScript,
              "\', output_file=\'",resultHTML, "\')\" --args ", 
              paste(inputArgs, collapse = ' ')))

# Move the tables to the output folder ----------------------------------------
currScriptDir <- dirname(markdownScript)
tabsToMove <- paste0(cmdArgs['RunID'], '_', cmdArgs['LibraryID'], '_',
                     cmdArgs['SampleID'], '*.csv')
tabsToMove <- paste0(currScriptDir, '/', tabsToMove)
system(paste('mv', tabsToMove, cmdArgs['outputDirPath']))

# Move html to the output folder, remove tmp subfolders -----------------------
system(paste('mv', resultHTML, cmdArgs['outputDirPath']))


