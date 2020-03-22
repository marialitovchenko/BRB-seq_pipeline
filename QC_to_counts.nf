#!/usr/bin/env nextflow


// path to the input table with samples
sampleTabPath = 'test_input/sampleTable.csv'
numbOfProc = 4
R1code = "_R1_"
R2code = "_R2_"
outputDir = '.'
trimmedDir=outputDir + "/trimmed"
fastqExtens=".fastq\$\\|.fq\$\\|.fq.gz\$\\|.fastq.gz\$"

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
* Trim reads by quality and adapterss
*----------------------------------------------------------------------------*/
process trimReads {
    input:
    set RunID, LibraryID, SampleID, Specie, Genome from sampleTab

    shell:
    '''
    R1path=$(find "../../../""!{RunID}" -type f | grep "!{LibraryID}" | \
             grep "!{SampleID}" | grep "!{R1code}" | grep "!{fastqExtens}")
    R2path=$(find "../../../""!{RunID}" -type f | grep "!{LibraryID}" | \
             grep "!{SampleID}" | grep "!{R2code}" | grep "!{fastqExtens}")

    echo "["$( date +%Y-%m-%d,%H-%M-%S )"]: Started trimming of "$R1path" & " \
         $R2path 
    #trim_galore -q 20 --length 20 --paired $R1path $R2path \
                -o !trimmedDir --fastqc
    echo "["$( date +%Y-%m-%d,%H-%M-%S )"]: Finished trimming of "$R1path" & "\
         $R2path
    '''
}
