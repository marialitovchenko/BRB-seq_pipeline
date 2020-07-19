env {
  NXF_OPTS="-Xms512M -Xmx2G"
}

executor {
    name = 'local'
    // The maximum number of CPUs made available by the underlying system
    // (only used by the local executor).
    cpus = 6
    // The maximum amount of memory made available by the underlying system
    //(only used by the local executor).
    memory = '20 GB'
}

process {
    //all workflow processes use 2 cpus if not otherwise specified in the
    //workflow script
    cpus = 2
    memory = { 2.GB * task.attempt }
    time = { 48.h * task.attempt }

    maxRetries = 3

    // Process-specific resource requirements
    withLabel: low_memory {
        memory = { 1.GB * task.attempt }
    }
    withLabel: mid_memory {
        memory = { 2.GB * task.attempt }
        time = { 60.h * task.attempt }
    }
    withLabel: high_memory {
        cpus = 3
        memory = { 2.GB * task.attempt }
        time = { 72.h * task.attempt }
    }

    withName: trimReads {
        time = { 72.h * task.attempt }
    }
    withName: mapWithStar {
        errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    }
    withName: countReads {
        errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    }
}

params {
    // general parameters
    R1code = "_R1_"
    R2code = "_R2_"
    fastqExtens = ".fastq\$\\|.fq\$\\|.fq.gz\$\\|.fastq.gz\$"
    
 
    // cutadapt parameters.
    // parameters, like trimGalore_gzip
    cutadapt_quality = 20
    cutadapt_adapter = 'AAGCAGTGGTATCAACGCAGAGTAC' // this is BRB-seq adapter
    cutadapt_allParams = '-B ' + params.cutadapt_adapter + 
                         ' -q ' + params.cutadapt_quality + ' ' +
                         '--cores=' + process.cpus

    // STAR parameters, not all, but main ones
    star_readFilesCommand = 'zcat'
    star_outSAMtype = 'BAM'
    star_sortBy = 'SortedByCoordinate'
    star_outFilterScoreMinOverLread = 0.66
    star_outFilterMatchNminOverLread = 0.66
    star_outFilterMultimapNmax = 1
    // assemble STAR parameters
    star_allParams = '--outFilterMultimapNmax ' + params.star_outFilterMultimapNmax + ' ' +
                     '--readFilesCommand ' + params.star_readFilesCommand + ' ' +
                     '--outSAMtype ' + params.star_outSAMtype + ' ' + params.star_sortBy + ' ' +
                     '--outFilterScoreMinOverLread ' + params.star_outFilterScoreMinOverLread + ' ' +
                     '--outFilterMatchNminOverLread ' + params.star_outFilterMatchNminOverLread

}
