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

/* ----------------------------------------------------------------------------
* Read inputs
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabCh = Channel.fromPath( sampleTabPath )
sampleTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.RunID, row.LibraryID, row.SampleID, row.Specie, 
                       row.Genome) }
    .set{ sampleTab }

/* ----------------------------------------------------------------------------
* Trim reads by quality and adapterss
*----------------------------------------------------------------------------*/
process trimReads {
    input:
    set RunID, LibraryID, SampleID, Specie, Genome from sampleTab

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
* Demultiplex reads
*----------------------------------------------------------------------------*/
process demultiplex {
    input:
    tuple val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
          path(trimmedR1), path(trimmedR2) from trimmedFiles

    output:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
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
* Map reads to reference genome with STAR
*----------------------------------------------------------------------------*/
process mapWithStar {
    // STAR is hungry for memory, so I give more
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    // tries 3 times, gives us afterwards
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2), 
        path(demultiplexfq) from demultiplexFiles

    output:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2), path(demultiplexfq), 
        path('*.sortedByCoord.out.bam'), 
        path('Log.final.out') into mappedBundle

    shell:
    '''
    STAR --runMode alignReads --runThreadN 1 \
                    --genomeDir /home/litovche/Documents/RefGen/chr21human/ \
                    --outFilterMultimapNmax 1 \
                    --readFilesCommand zcat \
                    --outSAMtype BAM SortedByCoordinate \
                    --readFilesIn "!{demultiplexfq}"
    '''
}

/* ----------------------------------------------------------------------------
* Count reads
*----------------------------------------------------------------------------*/
process countReads {
    // hungry for memory, so I give more
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    // tries 3 times, gives us afterwards
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2), path(demultiplexfq), path(mappedBam),
        path(mappedLog) from mappedBundle

    output:
    set val(RunID), val(LibraryID), val(SampleID), val(Specie), val(Genome),
        path(trimmedR1), path(trimmedR2), path(demultiplexfq), path(mappedBam),
        path(mappedLog), path('*.dge.umis.detailed.txt'),
        path('*.dge.reads.detailed.txt') into countedBundle

    shell:
    '''
    java -jar -Xmx2g "!{brbseqTools}" CreateDGEMatrix -f "!{trimmedR1}" \
         -b "!{mappedBam}" -c "../../../""!{barcodefile}" \
         -o "." -gtf "!{gtfPath}" -p BU -UMI "!{umiLen}"
    '''
}


countedBundle
    .flatMap { item ->
        forPrint = item[7].toString() + ' ' + item[1].toString();
    }
    .collectFile(name: 'sample.AGA2.txt', newLine: true)
    .set{fileForR}
