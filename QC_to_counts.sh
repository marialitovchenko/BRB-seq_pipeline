#!/bin/bash
###############################################################################
# Script	: QC_to_counts.sh
# Description	: Performs QC, alignment, demultiplexing and reads falling 
#		  within genes counting for BRB-seq samples sequenced on
#		  Illumina machine.
# Args          : 1) A tab-separated file with columns RunID, LibraryID, SampleID,
# 		  Specie, Genome; all other columns are ignored. Should contain
#                 header.
#                 2) Number of processes to run simulteniously
# Example       : ./QC_to_counts.sh 
# Author 	: Maria Litovchenko
# Email         : maria.litovchenko@epfl.ch
# TO DO		: put fool defense: if not all cells in the input table are 
#                 filled, if there is no such file in the directory, etc
# find: ftsopen: No such file or directory
###############################################################################

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# source the library of functions
source 0_functions.sh
# source the config file to get paths to all software etc
source 1_config.sh

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
R1code="_R1_"
R2code="_R2_"
trimmedDir=$outputDir"/trimmed"
createIfNotExist $trimmedDir

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

sampleCount=${#SAMPLES[@]}
for (( i=1; i<$sampleCount; i+=$numbOfProc )); do 
   pids="" # processes IDs
   endInd=$(($numbOfProc - 1))

   # lunch numbOfProc "jobs" at the same time
   for j in $(seq 0 $endInd); do
     position=$(( $i + $j ))
     
     # get fastq files corresponding to the sample
     currRun=${RUNS[$position]}
     currLib=${LIBRARIES[$position]}
     currSample=${SAMPLES[$position]}

     R1path=$(find "$currRun" -type f | grep "$currLib" | grep "$currSample" | grep "$R1code")
     R2path=$(find "$currRun" -type f | grep "$currLib" | grep "$currSample" | grep "$R2code")

     # perform trimming
     trimFastq $trimmedDir $R1path $R2path & 
     pids="$pids $!"
   done
   waitall $pids
done


#STAR --runMode alignReads --twopassMode Basic --outSAMmapqUnique 60 --runThreadN 8 --genomeDir $genomeDir/STAR_Index --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $bamDir/ --readFilesIn $inputDir/${fastqName}_R2_001.fastq.gz

# Demultiplex and generate output count/UMI matrix
# INPUT: R1.fastq and barcodes
# java -jar /software/BRBseqTools.1.5.jar CreateDGEMatrix -f $inputDir/${fastqName}_R1_001.fastq.gz -b $bamDir/Aligned.out.bam -o $resultDir -c $barcodefile -gtf $inputGTF -p BU???? -UMI 10

exit 0;
