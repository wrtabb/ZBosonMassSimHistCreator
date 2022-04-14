#!/bin/bash

echo "Start processing at " $(date)
pwd
ls -l

# Set up software and environment
echo "Set up environment"
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW CMSSW_10_6_1
cd CMSSW_10_6_1/src
eval `scramv1 runtime -sh`
cd -

#directory to save output data
mkdir output_data

# Run Processing code 
echo "Run the tree analysis script"
echo "ProcessFiles.py $1 $2 $3 $4"
python3 ProcessFiles.py $1 $2 $3 $4

echo "Ending processing at " $(date)

