#!/bin/bash
###############################################################################
# Script	: QC_to_counts.sh                                                                                       
# Description	: Performs QC, alignment, demultiplexing and reads falling 
#		  within genes counting for BRB-seq samples sequenced on 
#		  Illumina machine.
# Args          : A tab-separated file with columns RunID, LibraryID, SampleID,
# 		  Specie, Genome; all other columns are ignored. Should contain
#                 header.                                                                                         
# Author 	: Maria Litovchenko                             
# Email         : maria.litovchenko@epfl.ch                                      
###############################################################################

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# source the library of functions
source 0_functions.sh

# path to the input table
sampleTab=$1

# if user didn't submit output directory, inform about it
if (( $# == 1 )); 
then
   msg=": No output directory was submitted, using current directory as output"
   echo $( currentTime ) $msg
   outputDir="."
else
   outputDir=$2
   mkdir $outputDir
   msg=": Created directory "$outputDir
   echo $( currentTime ) $msg
fi

numbOfProc=4

# -----------------------------------------------------------------------------
# Read inputs 
# -----------------------------------------------------------------------------
# create arrays which will hold information about input files
RUNS=(); LIBRARIES=(); SAMPLES=(); SPECIES=(); GENOMES=();
while read runid library sampleid specie genome; do 
	RUNS+=($runid)
	LIBRARIES+=($library)
	SAMPLES+=($sampleid)
	SPECIES+=($specie)
	GENOMES+=($genome)
done < $sampleTab

echo $( currentTime )  ": Read file" $sampleTab
echo "	Number of submitted runs:	" $( numbUniqItems "${RUNS[@]}" )
echo "	Number of submitted libraries:	" $( numbUniqItems "${LIBRARIES[@]}" )
echo "	Number of submitted samples:	" $( numbUniqItems "${SAMPLES[@]}" )
echo "	Number of submitted species:	" $( numbUniqItems "${SPECIES[@]}" )
echo "	Number of submitted genomes:	" $( numbUniqItems "${GENOMES[@]}" )

# Note: I do not run fastqc here because it's going to be done after trimming
# anyway

# -----------------------------------------------------------------------------
# Trim reads 
# -----------------------------------------------------------------------------
echo $( currentTime )  ": Started trimming reads"

libCount=${#LIBRARIES[@]}
for (( i=1; i<$libCount; i+=$numbOfProc )); do
   # bash arrays are 0-indexed
   pids=""
   endInd=$(($numbOfProc - 1))

   for j in $(seq 0 $endInd); do
      position=$(( $i + $j ))

      # Get srr, dgrp and instrument
      currSRR=${SRRs[$position]}
      currDGRP=${DGRPs[$position]}
      currInstrument=${Intrument[$position]}

      echo $refGen $currSRR $currDGRP $currInstrument $fastqDir $trimmedDir $bwaDir $deduplDir

      #trimFastq  1>$currSRR".trimMapDedupl.out" 2>$currSRR".trimMapDedupl.err" &
      pids="$pids $!"
   done
   waitall $pids
done

trimFastq 

STAR --runMode alignReads --twopassMode Basic --outSAMmapqUnique 60 --runThreadN 8 --genomeDir $genomeDir/STAR_Index --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $bamDir/ --readFilesIn $inputDir/${fastqName}_R2_001.fastq.gz

# Demultiplex and generate output count/UMI matrix
# INPUT: R1.fastq and barcodes
# java -jar /software/BRBseqTools.1.5.jar CreateDGEMatrix -f $inputDir/${fastqName}_R1_001.fastq.gz -b $bamDir/Aligned.out.bam -o $resultDir -c $barcodefile -gtf $inputGTF -p BU???? -UMI 10

exit 0;
