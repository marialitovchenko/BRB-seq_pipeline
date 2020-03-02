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

# Functions -------------------------------------------------------------------

#===  FUNCTION  ===============================================================
# NAME		:  numbUniqItems
# DESCRIPTION	: returns number of unique items in array
# PARAMETER  1	: array
# RETURN	: a number
#==============================================================================
function numbUniqItems () {
  arr=("$@")
  numb=$(echo "${arr[@]}" | tr ' ' '\n' | sort | uniq | wc -l)
  numb="$(($numb-1))"
  echo $numb
}

# Inputs ----------------------------------------------------------------------
# path to the input table
sampleTab=$1

# Read inputs -----------------------------------------------------------------
# create arrays which will hold information about input files
RUNS=(); LIBRARIES=(); SAMPLES=(); SPECIES=(); GENOMES=();
while read runid library sampleid specie genome; do 
	RUNS+=($runid)
	LIBRARIES+=($library)
	SAMPLES+=($sampleid)
	SPECIES+=($specie)
	GENOMES+=($genome)
done < $sampleTab

echo "Read file" $sampleTab
echo "	Number of submitted runs:	" $( numbUniqItems "${RUNS[@]}" )
echo "	Number of submitted libraries:	" $( numbUniqItems "${LIBRARIES[@]}" )
echo "	Number of submitted samples:	" $( numbUniqItems "${SAMPLES[@]}" )
echo "	Number of submitted species:	" $( numbUniqItems "${SPECIES[@]}" )
echo "	Number of submitted genomes:	" $( numbUniqItems "${GENOMES[@]}" )

