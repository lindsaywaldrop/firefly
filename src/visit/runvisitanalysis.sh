#!/bin/bash

WD=${1:?Please provide a top-level working directory}
species=${2:?Please provide a species name}
startrun=${3:?Provide a start run number}
endrun=${4:?Provide an end run number}

for i in `seq ${startrun} ${endrun}`; do
  if [ -d "$WD"/results/visit/${species}/sim${i}/ ]; then
    rm "$WD"/results/visit/${species}/sim${i}/*.curve
  else
    mkdir -p "$WD"/results/visit/${species}/sim${i}/
  fi
done


# Script for running several VisIt python scripts

# Check to see if all files are present first!! 

/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s Shear-lineout.py \
"$WD"/results/ibamr/${species} "$WD"/results/visit/${species} $startrun $(($endrun+1))

/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s Across-lineout.py \
"$WD"/results/ibamr/${species} "$WD"/results/visit/${species} $startrun $(($endrun+1))