#!/usr/bin/env nextflow

// path to the input table with samples
sampleTabPath = 'test_input/sampleTable.csv'
numbOfProc = 4

fastqExtens=".fastq\$\\|.fq\$\\|.fq.gz\$\\|.fastq.gz\$"

R1code = "_R1_"
R2code = "_R2_"

brbseqTools="/home/litovche/bin/BRBseqTools.1.5.jar"
barcodefile="test_input/barcodes_v3.txt"
umiLen=10

gtfPath="/home/litovche/Documents/RefGen/chr21human/hg38.refGene.gtf"

rInputTab="/home/litovche/Desktop/BRB-seq_pipeline/rInputTab.txt"


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
    publishDir "trimmed", pattern: '*_val_*.fq.gz' 

    input:
    tuple RunID, LibraryID, SampleID, Specie, Genome from sampleTab

    output:
    tuple RunID, LibraryID, SampleID, Specie, Genome, 
          "${SampleID}_val_1.fq.gz", 
          "${SampleID}_val_2.fq.gz" into trimmedFiles

    shell:
    '''
    # full paths for R1 and R2
    R1path=$(find "../../../""!{RunID}" -type f | grep "!{LibraryID}" | \
             grep "!{SampleID}" | grep "!{R1code}" | grep "!{fastqExtens}")
    R2path=$(find "../../../""!{RunID}" -type f | grep "!{LibraryID}" | \
             grep "!{SampleID}" | grep "!{R2code}" | grep "!{fastqExtens}")
    # perform trimming with trim galore
    trim_galore -q 20 --length 20 --paired $R1path $R2path --fastqc \
                --basename "!{SampleID}"
    '''
}

/* ----------------------------------------------------------------------------
* !!! IMPORTANT NOTE: !!!
* With use of BRB-seq tools, it is not nessecary to demultiplex and map files
* one by one in order to get to the count table. Count table can be obtained
* firectly from maped trimmed file. However, it's not possible then to derive
* percentage of unmapped reads, multiple mapping percentage, etc, from it. This
* is why we still need to map the individual demultiplexed files. 
*
* Therefore, here I split channel into 2, there first will go without 
* deduplication directly to mapping and read counting, and second channel will
* go through demultiplexing and mapping in order to get mapping statistics
*----------------------------------------------------------------------------*/
trimmedFiles
     .into{ trimmedForCounts; trimmedForMapStats }

/* ----------------------------------------------------------------------------
* !!! PROCESSING 1: FROM TRIMMED INTO COUNTS, no deduplication
*----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
* A) Map trimmed reads to reference genome with STAR
*----------------------------------------------------------------------------*/
process mapTrimmedWithStar {
    publishDir "trimmedMapped", 
               saveAs: { LibraryID + "_" + SampleID + "_" + Genome + 
                        "_Aligned.sortedByCoord.out.bam" },
               pattern: '{*.sortedByCoord.out.bam}'
    cpus 4

    // STAR is hungry for memory, so I give more; tries 3 times, gives us
    // afterwards
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2) from trimmedForCounts

    output:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2),
        path('*.sortedByCoord.out.bam') into trimmedMappedBundle

    shell:
    '''
    STAR --runMode alignReads --runThreadN 4 \
                    --genomeDir /home/litovche/Documents/RefGen/chr21human/ \
                    --outFilterMultimapNmax 1 \
                    --readFilesCommand zcat \
                    --outSAMtype BAM SortedByCoordinate \
                    --readFilesIn "!{trimmedR2}"
    '''
}

/* ----------------------------------------------------------------------------
* B) Count reads in mapped trimmed bams
*----------------------------------------------------------------------------*/
process countReads {
    publishDir "counts/${LibraryID}/${SampleID}",
               pattern: '{*.dge.umis.detailed.txt, *.dge.reads.detailed.txt}'

    // Hungry for memory, so I give more, tries 3 times, gives us afterwards
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2), 
          path(mappedBam) from trimmedMappedBundle

    output:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2), path(mappedBam),
          path('*.dge.umis.detailed.txt'),
          path('*.dge.reads.detailed.txt') into countedBundle

    shell:
    '''
    java -jar -Xmx2g "!{brbseqTools}" CreateDGEMatrix -f "!{trimmedR1}" \
         -b "!{mappedBam}" -c "../../../""!{barcodefile}" \
         -o "." -gtf "!{gtfPath}" -p BU -UMI "!{umiLen}"
    '''
}

/* ----------------------------------------------------------------------------
* !!! PROCESSING 2: DEMULTIPLEXING TRIMMED, MAPPING, GETTING MAP STATS
*----------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------
* Demultiplex reads
*----------------------------------------------------------------------------*/
process demultiplex {
    publishDir "demultiplexed/${LibraryID}/${SampleID}", pattern: '*.fastq.gz'

    input:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2) from trimmedForMapStats

    output:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2),
          path('*.fastq.gz') into demultiplexBundle

    shell:
    '''
    java -jar "!{brbseqTools}" Demultiplex \
                               -r1 "!{trimmedR1}" -r2 "!{trimmedR2}" \
                               -c "../../../""!{barcodefile}" \
                               -p BU -UMI "!{umiLen}" -o "."

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
    publishDir "mapped/${LibraryID}/${SampleID}",
               pattern: '{*.sortedByCoord.out.bam, *._Log.final.out}'

    // STAR is hungry for memory, so I give more; tries 3 times, gives us
    // afterwards
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2),
          path(demultiplexfq) from demultiplexFiles

    output:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2), path(demultiplexfq),
          path('*.sortedByCoord.out.bam'),
          path('*_Log.final.out') into mappedBundle

    shell:
    '''
    mapPrefName=`basename "!{demultiplexfq}" | sed 's/[.].*//g'`
    mapPrefName="!{LibraryID}""_""!{SampleID}""_"$mapPrefName"_"
    STAR --runMode alignReads --runThreadN 1 \
                    --genomeDir /home/litovche/Documents/RefGen/chr21human/ \
                    --outFilterMultimapNmax 1 \
                    --readFilesCommand zcat \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFileNamePrefix $mapPrefName \
                    --readFilesIn "!{demultiplexfq}"
    '''
}
