#!/usr/bin/env nextflow

// path to the input table with samples
sampleTabPath = 'test_input/sampleTable.csv'


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
process sayHello {

    """
    echo 'Hello world!' > aga.txt
    """

}


process trimReads {
    input:
    set RunID, LibraryID, SampleID, Specie, Genome from sampleTab

    shell:
    '''
    R1path=$(find "!{RunID}" -type f | grep "!{LibraryID}" | grep "!{SampleID}")
    echo !R1path
    '''
}
