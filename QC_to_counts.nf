#!/usr/bin/env nextflow

def helpMessage() {

    log.info """ 
    B R B - s e q   N E X T F L O W   P I P E L I N E    
    ================================================================================
    Welcome to the Nextflow BRB-seq analysis command line pipeline!

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow forTest.nf --inputTab table.csv --FQdir fqDir --genomeDir allGenomes --outputDir theResult

    or

    nextflow forTest.nf --inputTab table.csv

    Mandatory arguments:
      --inputTab        Path to the table containing information about input 
                        data. The table should have following columns: RunID,
                        (i.e. NXT0540), LibraryID (i.e. nxid12916), SampleID
                        (i.e. BRBseq_v3_plate_1_S25), Specie (i.e. Hsapiens),
                        Genome (i.e. hg38). Specie and Genome indicate to which
                        genome version of which specie sample should be aligned
                        to.  

                        If no FQdir is provided (see below), the system will 
                        assume that input fastq files are located in 
                        [current dir]/RunID/LibraryID. 

                        If no genomeDir is provided (see below), the system 
                        will that STAR indexed genome is located in
                        [current dir]/Specie/Genome

                        If no outputDir is provided (see below), the system 
                        will output files in the current directory

    
    Optional arguments:
    This arguments are not going to be needed with use of graphical user
    interface
      --FQdir           Path to the directory containing folders (one per run) 
                        with fastq files
      --genomeDir       Path to the directory containing STAR indexed genomes
      --outputDir       Path to the output directory
      """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/* ----------------------------------------------------------------------------
* Input handling
*----------------------------------------------------------------------------*/
// path to the input table with samples
sampleTabPath = file(params.inputTab)

//user directory: directory which contains all the fastq files
params.FQdir = file('.')
userDir = file(params.FQdir)
// genomes directory: directory with all compiled STAR indexed genomes
params.genomeDir = file('.')
genomePath = file(params.genomeDir)
// output folder
params.outputDir = file('.')
outputDir = file(params.outputDir) 

// also should be in tech dir
brbseqTools="/home/litovche/bin/BRBseqTools.1.5.jar"
numbOfProc = 4

rInputTab="/home/litovche/Desktop/BRB-seq_pipeline/rInputTab.csv"
mapStatsTab = outputDir + "/mapStatsTab.csv"


/* ----------------------------------------------------------------------------
* LOG: inform user about all the inputs
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabInfoCh = Channel.fromPath( sampleTabPath )
sampleTabInfoCh
    .splitCsv(header: true, sep:'\t')
    .into { logFqFiles; logGenomesFiles }

logFqFiles
    .map{ row ->  userDir.toString() + '/' + row.RunID.toString() + '/' + row.LibraryID.toString() }
    .unique()
    .toList()
    .set{logFqFiles}

logGenomesFiles
    .map{ row ->  genomePath.toString() + '/' + row.Specie.toString() + '/' + row.Genome.toString() }
    .unique()
    .toList()
    .set{logGenomesFiles}

log.info """\
            B R B - s e q   N E X T F L O W   P I P E L I N E    
         ================================================================================
         Submitted input table          : ${sampleTabPath}
         Expect to find fastq-s in      : ${logFqFiles.toString().replaceAll(/DataflowVariable.value../, '').replaceAll(/..$/, '')}
         Expect to find genomes in      : ${logGenomesFiles.toString().replaceAll(/DataflowVariable.value../, '').replaceAll(/..$/, '')}
         BRBseq tools in                : ${brbseqTools}
         Output folder                  : ${outputDir}

         Upon completion, following folders are going to be created:
         ${outputDir}/trimmed
         ${outputDir}/demultiplexed
         ${outputDir}/mapped
         ${outputDir}/mapStats
         ${outputDir}/counts
         ${outputDir}/countTables

         Mapping statistics could be found in: ${mapStatsTab}

         ================================================================================

         L E T' S   G O ! ! !    
         """
         .stripIndent()

/* ----------------------------------------------------------------------------
* Read input table
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabCh = Channel.fromPath( sampleTabPath )
sampleTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.RunID, row.LibraryID, row.SampleID, row.Specie, 
                       row.Genome) }
    .set{ sampleTab }

/* ----------------------------------------------------------------------------
* Trim reads by quality and adapterss with trimgalore
*----------------------------------------------------------------------------*/
process trimReads {
    publishDir "${outputDir}/trimmed", pattern: '*_val_*.fq.gz' 

    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome from sampleTab

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, 
          "${SampleID}_val_1.fq.gz", 
          "${SampleID}_val_2.fq.gz" into trimmedFiles

    shell:
    '''
    # full paths for R1 and R2
    R1path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R1code} | \
             grep "!{params.fastqExtens}")
    R2path=$(find !{userDir}'/'!{RunID} -type f | grep !{LibraryID} | \
             grep !{SampleID} | grep !{params.R2code} | \
             grep "!{params.fastqExtens}")
    # perform trimming with trim galore
    trim_galore --paired $R1path $R2path --basename !{SampleID} \
                !{params.trimGalore_allParams}
    '''
}

/* ----------------------------------------------------------------------------
* !!! IMPORTANT NOTE: !!!
* With use of BRB-seq tools, it is not nessecary to demultiplex and map files
* one by one in order to get to the count table. Count table can be obtained
* directly from maped trimmed file. However, it's not possible then to derive
* percentage of unmapped reads, multiple mapping percentage, etc, from it. This
* is why we still need to map the individual demultiplexed files. But it 
* doesn't make sense to do two mapping runs: for just trimmed bam and for the
* demultiplexed one, it will take twice much time. For example, just trimming +
* mapping + counting takes 14 minutes for the test file, and this is without
* demultiplexing and mapping individual files. On the other hand, trimming +
* demultiplexing + mapping + counting takes 26 minutes in total. This is due to
* the higher degree of parallelization. So, I will demultiplex and
* map and count and then constract one count table per library.
*----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
* Demultiplex reads
*----------------------------------------------------------------------------*/
process demultiplex {
    publishDir "${outputDir}/demultiplexed/${LibraryID}/${SampleID}", pattern: '*.fastq.gz'

    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, 
          trimmedR2 from trimmedFiles

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, 
          trimmedR2, path('*.fastq.gz') into demultiplexBundle

    shell:
    '''
    java -jar !{brbseqTools} Demultiplex -r1 !{trimmedR1} -r2 !{trimmedR2} \
                               -c !{params.barcodefile} -o "." \
                               !{params.brbseqTools_commonParams}
    '''
}

// fork into #(of demultiplexed file) channels, preserving metadata
demultiplexBundle
    .flatMap { item ->
        RunID = item[0];
        LibraryID = item[1];
        SampleID = item[2];
        Specie = item[3];
        Genome = item[4];
        trimmedR1 = item[5];
        trimmedR2 = item[6];
        files  = item[7];
        files.collect { onefile -> return [ RunID, LibraryID, SampleID, Specie,
                        Genome, trimmedR1, trimmedR2, onefile ] }
    }
    .set { demultiplexFiles }

/* ----------------------------------------------------------------------------
* Map demultiplexed reads to reference genome with STAR
*----------------------------------------------------------------------------*/
process mapWithStar {
    publishDir "${outputDir}/mapped/${LibraryID}/${SampleID}", 
                pattern: '*.sortedByCoord.out.bam'
    publishDir "${outputDir}mapStats/${LibraryID}/${SampleID}", 
                pattern: '*_Log.final.out'

    // STAR is hungry for memory, so I give more; tries 3 times, gives us
    // afterwards
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, trimmedR2,
          demultiplexfq from demultiplexFiles

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, trimmedR2,
          demultiplexfq, path('*.sortedByCoord.out.bam'),
          path('*_Log.final.out') into mappedBundle

    shell:
    '''
    mapPrefName=`basename !{demultiplexfq} | sed 's/[.].*//g'`
    mapPrefName=$mapPrefName"_"
    STAR --runMode alignReads --readFilesIn !{demultiplexfq} \
         --genomeDir !{genomePath} --outFileNamePrefix $mapPrefName \
         !{params.star_allParams}
    '''
}

/* ----------------------------------------------------------------------------
* I will split the channel here, one will go to the aggregation of mapping
* statistics, and the other one - into counting reads
*----------------------------------------------------------------------------*/
mappedBundle.into{ mappedForStats; mappedForCounts }

/* ----------------------------------------------------------------------------
* Aggregats mapping statistics
* I'll create a string with all sample info and will append to it mapping stats
* taken from STAR
*----------------------------------------------------------------------------*/
process aggregateMapStats {
    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, trimmedR2,
          demultiplexfq, mappedBam, mappedLog from mappedForStats

    output:
    stdout mappingStatsAggr
 
    shell:
    '''
    # initial sample info
    statsAggr=(!{RunID})
    statsAggr+=(!{LibraryID})
    statsAggr+=(!{SampleID})
    statsAggr+=(!{Specie})
    statsAggr+=(!{Genome})
    statsAggr+=(!{demultiplexfq})

    # append mapping stats info
    statsAggr+=(`grep "Number of input reads" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Uniquely mapped reads number" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads mapped to multiple loci" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads mapped to too many loci" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: too many mismatches" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: too short" !{mappedLog} | sed 's/.*|//'`)
    statsAggr+=(`grep "Number of reads unmapped: other" !{mappedLog} | sed 's/.*|//'`)
    echo "${statsAggr[@]}"
    '''
}

mappingStatsAggr
    .collectFile(name: mapStatsTab, newLine: false)

/* ----------------------------------------------------------------------------
* Count reads in demultiplexed trimmed bams
*----------------------------------------------------------------------------*/
process countReads {
    publishDir "${outputDir}/counts/${LibraryID}/${SampleID}", pattern: '{*.detailed.txt}'
    errorStrategy 'ignore'

    // Hungry for memory, so I give more, tries 3 times, gives us afterwards
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, trimmedR2,
          demultiplexfq, mappedBam, mappedLog from mappedForCounts

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, trimmedR1, trimmedR2,
          demultiplexfq, mappedBam, mappedLog, path('*.umis.detailed.txt'), 
          path('*.reads.detailed.txt') into countedBundle

    shell:
    '''
    gtfPath=`find !{genomePath} | grep .gtf$`
    java -jar -Xmx2g !{brbseqTools} CreateDGEMatrix -f !{trimmedR1} \
         -b !{mappedBam} -c !{params.barcodefile} -o "." \
         -gtf $gtfPath !{params.brbseqTools_commonParams}
   
    samplName=`basename !{mappedBam} | sed 's/_Aligned.sortedByCoord.out.bam/.count/g'`
    mv output.dge.reads.detailed.txt $samplName".dge.reads.detailed.txt"
    mv output.dge.reads.txt $samplName".dge.reads.txt"
    mv output.dge.umis.detailed.txt $samplName".dge.umis.detailed.txt"
    mv output.dge.umis.txt $samplName".dge.umis.txt"
    '''
}

/* ----------------------------------------------------------------------------
* Output paths to files in the table for R
*----------------------------------------------------------------------------*/
countedBundle
    .flatMap { item ->
        item[0].toString() + ' ' + item[1].toString() + ' ' +
        item[2].toString() + ' ' + item[3].toString() + ' ' +
        item[4].toString() + ' ' + item[5].toString() + ' ' +
        item[6].toString() + ' ' + item[7].toString() + ' ' +
        item[8].toString() + ' ' + item[9].toString() + ' ' +
        item[10].toString() + ' ' + item[11].toString()
    }
    .collectFile(name: rInputTab, newLine: true)
    .set{fileForR}

/* ----------------------------------------------------------------------------
* Merge count tables per sample into 1 count table
*----------------------------------------------------------------------------*/
process mergeCounts {
    publishDir "${outputDir}/countTables", pattern: '{*Combined.csv}'

    input:
        path inputForR from fileForR

    output:
        path('*Combined.csv') into countTables

    shell:
    '''
    Rscript --vanilla /home/litovche/Desktop/BRB-seq_pipeline/combineCounts.R !{inputForR} reads readsCombined
    Rscript --vanilla /home/litovche/Desktop/BRB-seq_pipeline/combineCounts.R !{inputForR} UMI umiCombined
    '''
}
