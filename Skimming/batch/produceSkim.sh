#!/bin/bash

##Source relevant stuff
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
export SCRAM_ARCH=slc6_amd64_gcc700

##Set CMSSW
eval `scramv1 project CMSSW CMSSW_10_4_0`
cd CMSSW_10_4_0/src/
eval `scramv1 runtime -sh`

##Move ChargeHiggs framework into CMSSW src dir
mv ../../ChargedHiggs/ ./
mv ../../x509 ./
export X509_USER_PROXY=$CMSSW_BASE/src/x509

scram b

##Copy file from DAS
xrdcp $1 nanoFile.root

##Do the skimming
nanoskimmer.py --filename nanoFile.root --out-name $2 --channel ${@:3}
rm nanoFile.root

##Move output to base dir
mv *.root ../../
