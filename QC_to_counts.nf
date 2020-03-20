#!/usr/bin/env nextflow

// path to the input table with samples
sampleTabPath = 'test_input/sampleTable.csv'
numbOfProc = 4
R1code = "_R1_"
R2code = "_R2_"
//fastqExtens=".fastq\$\|.fq\$\|.fq.gz\$\|.fastq.gz\$"
fastqExtens=".fastq\$\\|.fq\$"

/* ----------------------------------------------------------------------------
* Read inputs
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
sampleTabCh = Channel.fromPath( sampleTabPath )
sampleTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.RunID, row.LibraryID, row.SampleID, row.Specie, 
                       row.Genome) }
    .set { sampleTab }

/* ----------------------------------------------------------------------------
* Read inputs
*----------------------------------------------------------------------------*/
process trimReads {
    input:
    set RunID, LibraryID, SampleID, Specie, Genome from sampleTab

    shell:
    '''
    R1path=$(find "!{RunID}" -type f | grep "!{LibraryID}" | \
             grep "!{SampleID}" | grep "!{R1code}" | grep "!{fastqExtens}")
    echo !R1path
    '''
}
