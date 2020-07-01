#!/usr/bin/env Rscript

cmdArgs <- c('Test User', 'Test PI', 'theResult', 'NXT0570', 'nxid13448', 'BRB_AM_100',
                 'danio_rerio', 'GRCz11.100_GFP')
names(cmdArgs) <- c('User', 'PI', 'ResultFolder', 'RunID', 'LibraryID', 
                    'SampleID', 'Specie', 'Genome')

# get fastqc text report for R1 and R2
fastqcDir <- paste(cmdArgs['ResultFolder'], 'fastQC', cmdArgs['LibraryID'], 
                   cmdArgs['SampleID'], sep = '/')
fastqcR1zip <- list.files(fastqcDir, pattern = '_R2_.*zip$', full.names = T) 
fastqcR2zip <- list.files(fastqcDir, pattern = '_R2_.*zip$', full.names = T) 
system(paste('unzip -d', fastqcDir, '-u', fastqcR1zip))
system(paste('unzip -d', fastqcDir, '-u', fastqcR2zip))
fastqcR1file <- list.files(gsub('.zip', '/', fastqcR1zip), 
                           pattern = 'fastqc_data.txt', full.names = T)
fastqcR2file <- list.files(gsub('.zip', '/', fastqcR2zip), 
                           pattern = 'fastqc_data.txt', full.names = T)

# get trimming report from trim_galore and fastqc for trimmed R2
trimDir <-  paste(cmdArgs['ResultFolder'], 'trimmed', cmdArgs['LibraryID'],
                   cmdArgs['SampleID'], sep = '/')
trimRepR2 <- list.files(trimDir, pattern = '_trimming_report.txt$', 
                        full.names = T)
trimFastqcR2 <- list.files(trimDir, pattern = 'val_2_fastqc.zip$', 
                           full.names = T)
system(paste('unzip -d', trimDir, '-u', trimFastqcR2))
trimFastqcR2 <- list.files(gsub('.zip', '/', trimFastqcR2), 
                           pattern = 'fastqc_data.txt', full.names = T)  

# demultiplexing statistics
demultiplStatsFile <-  paste(cmdArgs['ResultFolder'], 'demultiplexed', 
                             cmdArgs['LibraryID'], cmdArgs['SampleID'], 
                             'stats.txt', sep = '/')
# mapping statistics
mapStatsFile <- paste(cmdArgs['ResultFolder'], 'mapStatsTab.csv', sep = '/')
# count table
countTabFile <- paste0(cmdArgs['ResultFolder'], '/countTables/',
                       cmdArgs['RunID'], '_', cmdArgs['LibraryID'], '_',
                       cmdArgs['Genome'], '_', cmdArgs['SampleID'], 
                       '_readsCombined.csv')

print(fastqcR1file)
print(fastqcR2file)
print(trimFastqcR2)
print(demultiplStatsFile)
print(mapStatsFile)
print(countTabFile)

inputArgs <- c('Test User', 'Test PI', 'NXT0570', 'nxid13448', 'BRB_AM_100',
               'danio_rerio', 'GRCz11.100_GFP', 
               'test_input/fastQC/nxid12916/BRBseq_v3_plate_1_S25/BRBseq_v3_plate_1_S25_R1_001_fastqc/fastqc_data.txt',
               'test_input/fastQC/nxid12916/BRBseq_v3_plate_1_S25/BRBseq_v3_plate_1_S25_R1_001_fastqc/fastqc_data.txt',
               'test_input/BRBseq_v3_plate_1_S25_R2_001.fastq.gz_trimming_report.txt',
               'test_input/BRBseq_v3_plate_1_S25_val_1_fastqc/fastqc_data.txt',
               'test_input/stats.txt', 'test_input/mapStatsTab.csv',
               'test_input/NXT0570_nxid13444_GRCz11.100_GFP_BRB_AM_50_readsCombined.csv')
names(inputArgs) <- c('User', 'PI', 'RunID', 'LibraryID', 'SampleID', 'Specie',
                      'Genome', 'fastqcR1', 'fastqcR2', 'trimInfoR2', 'fastqcTrimR2',
                      'demultStats', 'mapStats', 'countTab')

# do some sort of processing/error checking
#  e.g. you could investigate the optparse/getopt packages which
#   allow for much more sophisticated arguments e.g. named ones.
rmarkdown::render('test.rmd')
