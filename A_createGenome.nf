#!/usr/bin/env nextflow

/* ----------------------------------------------------------------------------
* Inputs
*----------------------------------------------------------------------------*/
// genome version, it contains specie inside
genomeVersion="hg38"
genomeDir="/home/litovche/Documents/RefGen/chr21human/"

// common inputs: path to UCSC genomes
goldenPath="https://hgdownload.soe.ucsc.edu/goldenPath/"
// path to picard jar
picardJar=/home/litovche/bin/picard.jar

threadsNumb=12

/* ----------------------------------------------------------------------------
* Index genome for use with STAR
*----------------------------------------------------------------------------*/
process downloadGenome {
    publishDir genomeDir, pattern: '*.fa'

    input:
    tuple goldenPath, genomeVersion

    output:
    tuple '*.fa', '*.gtf' into genomeFiles

    shell:
    '''
    # Download genome from UCSC
    genomeFasta="!{genomeVersion}"".fa.gz"
    wget "!{goldenPath}""!{genomeVersion}""/chromosomes/chr21.fa.gz"
    gunzip chr21.fa.gz

    # Download gene annotation file
    wget "!{goldenPath}""!{genomeVersion}""/bigZips/genes/""!{genomeVersion}"".refGene.gtf.gz"
    gunzip "!{genomeVersion}"".refGene.gtf.gz"
    '''
}

process indexGenomeForSTAR {
    publishDir genomeDir

    input:
    tuple genomeFasta, genomeGTF from genomeFiles

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
