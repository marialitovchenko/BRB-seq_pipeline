#!/usr/bin/env nextflow

def helpMessage() {

    log.info """
    -      \033[41m B R B - s e q   N E X T F L O W  I N D E X  G E N O M E v1.0\033[0m-
    ================================================================================
    Welcome to the Nextflow BRB-seq genome indexing pipeline. This pipeline can be 
    used in two modes: 1) creating and indexing a "usual" = from ENSEML genome and
    2) creating and indexing a custom genome, i.e. then you have plasmids included
    into the genome or specie-mixed genome. For the mode 1 you don't need to 
    download anything, the download will be done automatically. For the mode 2 you
    need to provide fasta file and gtf file. 

    Usage:
    The \033[1;91mtypical\033[0m command for running the pipeline is as follows:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m folderWithGenomes 

    or

    to put the pipeline into \033[1;91mbackground\033[0m mode:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m folderWithGenomes \\
                        \033[1;91m-bg\033[0m \\
                        \033[1;91m-N\033[0m your.email@gmail.com

    \033[1;91mMandatory\033[0m arguments:
      \033[1;91m--inputTab\033[0m        Path to the table containing information about input
                        data.
                        In \033[93m mode 1 \033[0m (creating and indexing a "usual" = from ENSEML)
                        table should have just one column: GenomeCode (i.e. hs19,
                        hs38, ... , etc; check proper code for your specie of
                        interest on
                        https://hgdownload.soe.ucsc.edu/downloads.html).
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome) table
                        should have 3 columns: GenomeCode (desired name, i.e.
                        MixedHsMM), Fasta (name of the fasta file for your genome,
                        i.e. "MixedHsMM.fa", it must be located in the folder you
                        give as genomeDir argument(folderWithGenomes))) and GTF
                        (name of the GTF file for your genome, i.e. "MixedHsMM.gtf",
                        it must be located in the folder you give as genomeDir
                        argument, (folderWithGenomes))
      \033[1;91m--genomeDir\033[0m       Directory
                        In \033[93m mode 1 \033[0m (creating and indexing a "usual" = from ENSEML):
                        just a directory name. It will be created.
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome):
                        a directory containing Fasta and GTF files of your custom genome.

    \033[1;91mOptional\033[0m arguments:
    This arguments are not going to be needed with use of graphical user
    interface
      \033[1;91m--help\033[0m            Displays this message
      \033[1;91m-bg\033[0m               Puts execution of the pipeline into background mode
      \033[1;91m-N\033[0m                email adress in order to get notified upon pipeline complition.
                        Do not use epfl email address, because emails can't pass firewall.
                        Use gmail.
      \033[1;91m-resume\033[0m           Resumes execution of the pipeline from the moment it
                        was interrupted
      """.stripIndent()
}

// Show help message
params.help = ''
if (params.help) {
    helpMessage()
    exit 0
}

/* ----------------------------------------------------------------------------
* Input handling
*----------------------------------------------------------------------------*/
// path to the input table with specie name, genome version
genomeTabPath = file(params.inputTab)

// genomes directory (basically output folder): directory where indexed genome
// will be put to
params.genomeDir = file('.')
genomePath = file(params.genomeDir)

// common inputs: path to UCSC genomes
goldenPath="https://hgdownload.soe.ucsc.edu/goldenPath/"
// path to picard jar
picardJar=/home/litovche/bin/picard.jar

threadsNumb=6

/* ----------------------------------------------------------------------------
* Read input table
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
genomeTabCh = Channel.fromPath( genomeTabPath )
genomeTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.Specie, row.GenomeVersion, row.genomeFasta, row.genomeGTF) }
    .set{ genomeTab }

/* ----------------------------------------------------------------------------
* Download the genome if user wants new specie/version which we don't have
*----------------------------------------------------------------------------*/
//process downloadGenome {
//    publishDir genomeDir, pattern: '*.fa'
//
//    input:
//    tuple goldenPath, genomeVersion
//
//    output:
//    tuple '*.fa', '*.gtf' into genomeFiles

//    shell:
//    '''
//    # Download genome from UCSC
//    genomeFasta="!{genomeVersion}"".fa.gz"
//    wget "!{goldenPath}""!{genomeVersion}""/chromosomes/chr21.fa.gz"
//    gunzip chr21.fa.gz
//
//    # Download gene annotation file
//    wget "!{goldenPath}""!{genomeVersion}""/bigZips/genes/""!{genomeVersion}"".refGene.gtf.gz"
//    gunzip "!{genomeVersion}"".refGene.gtf.gz"
//    '''
//}

/* ----------------------------------------------------------------------------
* Index genome for use with STAR
*----------------------------------------------------------------------------*/
process indexGenomeForSTAR {
    publishDir genomeDir

    input:
    tuple Specie, GenomeVersion, genomeFasta, genomeGTF from genomeFiles
    //tuple genomeFasta, genomeGTF from genomeFiles

    shell:
    '''
    # index with STAR
    STAR --runMode genomeGenerate --runThreadN "!{threadsNumb}" --genomeDir "!{genomeDir}" \
    --genomeFastaFiles "!{genomeFasta}" --sjdbGTFfile "!{genomeGTF}"
    # created files: chrLength.txt, chrStart.txt, exonGeTrInfo.tab,
    # genomeParameters.txt, SA, sjdbList.fromGTF.out.tab,
    # chrNameLength.txt, exonInfo.tab, geneInfo.tab, Log.out, SAindex,
    # sjdbList.out.tab, chrName.txt, Genome, sjdbInfo.txt, transcriptInfo.tab

    # index with samtools
    samtools faidx "!{genomeFasta}"
    # created files :  *.fai

    # index picard
    dictFile=$(echo "!{genomeFasta}" | sed "s/fa$/dict/g")
    java -jar "!{picardJar}" CreateSequenceDictionary \
              REFERENCE="!{genomeFasta}" OUTPUT=$dictFile
    # created files: *.dict
    '''
}
