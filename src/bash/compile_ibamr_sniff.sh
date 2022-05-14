#!/bin/bash

#Shell script for compiling the main2d program in IBAMR. To be run on Keck Cluster.
WD=${1:?Please provide a path for the main directory}
Species=${2:?Please provide a species name}

cd "$WD"
echo "Setting up directories..."
echo $WD
echo "for "
echo $Species
mkdir bin/
mkdir results/
mkdir results/ibamr
mkdir results/ibamr/runs/
mkdir results/ibamr/log-files/
mkdir data/input2d-files/
mkdir data/parameters-files/
mkdir data/csv-files/
mkdir data/csv-files/${Species}/
mkdir data/ibamr-files/
mkdir data/ibamr-files/${Species}/

cd "$WD"/bin
echo "Compiling main2d..."
cp "$WD"/src/ibamr/* .
make main2d
rm *.o stamp-2d *.C Makefile *.h

cd "$WD"/src/bash/

echo "Complete"
