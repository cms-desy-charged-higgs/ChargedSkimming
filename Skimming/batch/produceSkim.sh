#!/bin/bash
source $CHDIR/ChargedAnalysis/setenv.sh CMSSW
cd $TMP

INFILE=$1
ISDATA=$2
CHANNEL=$3
XSEC=$4
OUT=$5
JOB=$6

##Copy file from DAS
xrdcp $1 nanoFile.root

##Do the skimming
NanoSkim nanoFile.root $ISDATA "$CHANNEL" $XSEC ${OUT}_${JOB}.root
rm nanoFile.root

##Move output to base dir
mkdir -p $CHDIR/Skim2/$OUT/output/
mv ${OUT}_${JOB}.root $CHDIR/Skim2/$OUT/output
