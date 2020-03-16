#!/bin/bash
###############################################################################
# Script	: 0_functions.sh                                                                                       
# Description	: 
# Args          :                                                                                        
# Author 	: Maria Litovchenko                             
# Email         : maria.litovchenko@epfl.ch                                      
###############################################################################

#===  FUNCTION  ===============================================================
# NAME          : currentTime
# DESCRIPTION   : displays current time
# RETURN        : a string
#==============================================================================
function currentTime () {
  echo "["$( date +%Y-%m-%d,%H-%M-%S )"]"
}

#===  FUNCTION  ===============================================================
# NAME		:  numbUniqItems
# DESCRIPTION	: returns number of unique items in array
# PARAMETER  1	: array
# RETURN	: a number
#==============================================================================
function numbUniqItems () {
  local arr=("$@")
  local numb=$(echo "${arr[@]}" | tr ' ' '\n' | sort | uniq | wc -l)
  local numb="$(($numb-1))"
  echo $numb
}

#===  FUNCTION  ===============================================================
# NAME:  trimFastq
# DESCRIPTION: trims reads with trim-galore based on quality
# PARAMETER  1: path to folder to put trimmed files
# PARAMETER  2: path to fastq file - read 1
# PARAMETER  3: path to fastq file - read 2
#==============================================================================
trimFastq() {
   if (( $# == 2 )); then
      echo "Started trimming of "$2" at "currentTime
      trim_galore -q 20 $2 -o $1 --fastqc
      echo "Finished trimming of "$2" at "currentTime
   fi

   if (( $# == 3 )); then
      echo "Started trimming of "$2" and "$3" at "currentTime
      trim_galore -q 20 --paired $2 $3 -o $1 --fastqc
      echo "Finished trimming of "$2" and "$3" at "currentTime
   fi
}


#===  FUNCTION  ===============================================================
# NAME:  waitall
# DESCRIPTION: waits to run the next batch of processes till current batch is
#	       finished
#==============================================================================
waitall() { # PID...
  local errors=0
  while :; do
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
        set -- "$@" "$pid"
      elif wait "$pid"; then
        echo "COMPLETED: $pid exited with zero exit status." >&2;
      else
        echo "ERROR: $pid exited with non-zero exit status." >&2;
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    sleep ${WAITALL_DELAY:-1}
   done
  ((errors == 0))
}

