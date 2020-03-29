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
    path '*.fastq.gz' into demultiplexBundle 

    shell:
    '''
    java -jar "!{brbseqTools}" Demultiplex \
                               -r1 "!{trimmedR1}" -r2 "!{trimmedR2}" \
                               -c "../../../""!{barcodefile}" \
                               -p BU -UMI "!{umiLen}" -o "."
 
    '''
}

process mapWithStar {
    input: val(x) from demultiplexBundle.flatMap()

    output: stdout resultB

    """
    echo $x
    """
}
