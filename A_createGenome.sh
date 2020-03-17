#!/bin/bash
###############################################################################
# Script    : A_createGenome.sh
# Description    :
# Args          :
# Example       : ./A_createGenome.sh
# Author     : Maria Litovchenko
# Email         : maria.litovchenko@epfl.ch
###############################################################################

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# genome version, it contains specie inside
genomeVersion=$1
genomeDir=$2

genomeVersion="hg38"
genomeDir="/home/litovche/Documents/RefGen/chr21human/"

# common inputs: path to UCSC genomes
goldenPath="https://hgdownload.soe.ucsc.edu/goldenPath/"
# path to picard jar
picardJar=/home/litovche/bin/picard.jar

# go to the future genome directory
cd $genomeDir

# -----------------------------------------------------------------------------
# Index genome for use with STAR
# -----------------------------------------------------------------------------
# 1) Download genome from UCSC
    genomeFasta=$genomeVersion".fa.gz"
	wget $goldenPath$genomeVersion"/chromosomes/chr21.fa.gz"
    gunzip chr21.fa.gz
    genomeFasta=$genomeDir"chr21.fa"

# 2) Download gene annotation file
     wget $goldenPath$genomeVersion"/bigZips/genes/"$genomeVersion".refGene.gtf.gz"
     gunzip $genomeVersion".refGene.gtf.gz"
     genomeGTF=$genomeDir$genomeVersion".refGene.gtf"

# 3) Index with STAR
	STAR --runMode genomeGenerate --runThreadN 12 --genomeDir $genomeDir \
         --genomeFastaFiles $genomeFasta --sjdbGTFfile $genomeGTF
# created files: chrLength.txt, chrStart.txt, exonGeTrInfo.tab, 
# genomeParameters.txt, SA, sjdbList.fromGTF.out.tab,
# chrNameLength.txt, exonInfo.tab, geneInfo.tab, Log.out, SAindex, 
# sjdbList.out.tab, chrName.txt, Genome, sjdbInfo.txt, transcriptInfo.tab

# 4) Index with samtools
	samtools faidx $genomeFasta
# created files :  *.fai

# 5) Index with Picard
    dictFile=$(echo $genomeFasta | sed "s/fa$/dict/g")
	java -jar $picardJar CreateSequenceDictionary \
    		REFERENCE=$genomeFasta \
    		OUTPUT=$dictFile
# created files: *.dict

exit 0;
