#!/bin/bash

source $CHDIR/ChargedAnalysis/setenv.sh CMSSW

hadd -f "$@"
