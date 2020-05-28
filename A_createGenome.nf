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
                        table should have just two columns: Specie (i.e. homo_sapience),
                        and GenomeCode (i.e. hs19, hs38, ... , etc; check proper code 
                        for your specie of interest on
                        https://hgdownload.soe.ucsc.edu/downloads.html).
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome) table
                        should have 4 columns: Specie (i.e. homo_sapience), GenomeCode 
                        (desired name, i.e. MixedHsMM), Fasta (name of the fasta file for
                        your genome, i.e. "MixedHsMM.fa", it must be located in the 
                        genomeDir/specie_name/GenomeCode 
                        (i.e. folderWithGenomes/homo_sapience/myCustomGenome) 
                        and GTF (name of the GTF file for your genome, i.e. "MixedHsMM.gtf",
                        it must be located in the genomeDir/specie_name
                        (i.e. folderWithGenomes/homo_sapience/myCustomGenome)
                        \033[93m You can mix two modes in one table, i.e. \033[0m
                        Specie  GenomeCode  Fasta   GTF
                        homo_sapience   hg19
                        homo_sapience   myCustomGenome    customGenome.fa    customGenome.gtf
      \033[1;91m--genomeDir\033[0m       Directory
                        In \033[93m mode 1 \033[0m (creating and indexing a "usual" = from ENSEML):
                        just a directory name. It will be created.
                        In \033[93m mode 2 \033[0m (creating and indexing a custom genome):
                        a directory containing specie_name directory containing Fasta and 
                        GTF files of your custom genome (i.e. folderWithGenomes/homo_sapience) 

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
// path to the input table (with columns in dependence with mode )
genomeTabPath = file(params.inputTab)

// genomes directory (basically output folder): directory where indexed genome
// will be put to
params.genomeDir = file('.')
genomePath = file(params.genomeDir)

// common inputs: path to UCSC genomes
goldenPath="https://hgdownload.soe.ucsc.edu/goldenPath/"
// path to picard jar
picardJar="/home/litovche/bin/picard.jar"

threadsNumb=6

/* ----------------------------------------------------------------------------
* Read input table
*----------------------------------------------------------------------------*/
// create channel which reads from the input table with samples
genomeTabCh = Channel.fromPath( genomeTabPath )
genomeTabCh
    .splitCsv(header: true, sep:'\t')
    .map{ row -> tuple(row.Specie, row.GenomeCode, row.Fasta, row.GTF) }
    .into{ genomeTab_download; genomeTab_custom }

genomeTab_download
    .filter{ it[2] == null }
    .filter{ it[3] == null}
    .set{genomeTab_download_flt}

genomeTab_custom
    .filter{ it[2] != null }
    .filter{ it[3] != null}
    .set{genomeTab_custom_flt}

/* ----------------------------------------------------------------------------
* Download the genome if user wants new specie/version which we don't have
*----------------------------------------------------------------------------*/
process downloadGenome {
    publishDir "${genomePath}/${Specie}/${GenomeCode}",
                mode: 'copy', overwrite: true

    input:
    tuple Specie, GenomeCode, Fasta, GTF from genomeTab_download_flt

    output:
    tuple Specie, GenomeCode, "*.fa", 
          "*.gtf" optional true into genomes_ensembl

    shell:
    '''
    # if both are null: this is the case to download from ENSEMBL
    if [ !{Fasta} = "null" ] && [ !{GTF} = "null" ]; then
        # Download genome from UCSC
        wget !{goldenPath}!{GenomeCode}/bigZips/!{GenomeCode}.fa.gz
        gunzip $genomeFasta

        # Download gene annotation file
        wget !{goldenPath}!{GenomeCode}/bigZips/genes/!{GenomeCode}.refGene.gtf.gz
        gunzip !{GenomeCode}".refGene.gtf.gz"
    fi
    '''
}

genomeTab_custom_flt
    .map { item ->
        GenomeCode = item[0];
        Fasta = genomePath + "/" + item[1];
        GTF = genomePath + "/" + item[2];
        return [ GenomeCode, Fasta, GTF ] 
    }
    .mix(genomes_ensembl)
    .set{genomesToIndex}

