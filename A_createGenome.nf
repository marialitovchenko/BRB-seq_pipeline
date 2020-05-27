#!/usr/bin/env nextflow

def helpMessage() {

    log.info """
    -      \033[41m B R B - s e q   N E X T F L O W  I N D E X  G E N O M E v1.0\033[0m-
    ================================================================================
    Welcome to the Nextflow BRB-seq genome indexing pipeline

    Usage:
    The \033[1;91mtypical\033[0m command for running the pipeline is as follows:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m allGenomes 

    or

    to put the pipeline into \033[1;91mbackground\033[0m mode:
    nextflow forTest.nf \033[1;91m--inputTab\033[0m table.csv \\
                        \033[1;91m--genomeDir\033[0m allGenomes \\
                        \033[1;91m-bg\033[0m \\
                        \033[1;91m-N\033[0m your.email@gmail.com

    \033[1;91mMandatory\033[0m arguments:
      \033[1;91m--inputTab\033[0m        Path to the table containing information about input
                        data. The table should have following columns: RunID,
                        (i.e. NXT0540), LibraryID (i.e. nxid12916), SampleID
                        (i.e. BRBseq_v3_plate_1_S25), Specie (i.e. Hsapiens),
                        Genome (i.e. hg38). Specie and Genome indicate to which
                        genome version of which specie sample should be aligned
                        to.
      \033[1;91m--genomeDir\033[0m       Path to the directory containing all your genome
                        versions for all your species. For example, a valid
                        genome directory TestGenomeDir would contain two
                        folders names mus_musculus and homo_sapiens.
                        Consequently, homo_sapiens folder would contain
                        GRCh37.75 and GRCh38.99, and mus_musculus would contain
                        GRCm38.68 and GRCm38.98. \033[93m Please use then homo_sapiens
                        or mus_musculus in a Specie column of your input table,
                        and use GRCh37.75/GRCh38.99/GRCm38.68/GRCm38.98 in a
                        Genome column.\033[0m

    \033[1;91mOptional\033[0m arguments:
    This arguments are not going to be needed with use of graphical user
    interface
      \033[1;91m--help\033[0m            Displays this message
      \033[1;91m-bg\033[0m               Puts execution of the pipeline into background mode
      \033[1;91m-N\033[0m                email adress in order to get notified upon pipeline complition.Do not use epfl email address, because emails can't pass firewall. Use gmail.
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
