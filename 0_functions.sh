#!/bin/bash
###############################################################################
# Script	: 0_functions.sh                                                                                       
# Description	: 
# Args          :                                                                                        
# Author 	: Maria Litovchenko                             
# Email         : maria.litovchenko@epfl.ch                                      
###############################################################################

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

