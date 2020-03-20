#!/usr/bin/env nextflow

// path to the input table with samples
sampleTabPath = 'test_input/sampleTable.csv'
// create channel which reads from the input table with samples
sampleTabCh = Channel.fromPath( sampleTabPath )
sampleTabCh
    .splitCsv(header: true, sep:'\t')
    .subscribe {row ->
       print "${row.RunID} - ${row.LibraryID} - ${row.SampleID} - "
       println "${row.Specie} - ${row.Genome}"	
    }

