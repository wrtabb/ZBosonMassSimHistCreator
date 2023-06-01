#!/bin/bash

# two options:
# $1 - the directory where the root files are located. By default, this directory is called 'output_data' because it is output by the condor jobs
# $2 - the directory where you want to store the combined files
input_dir=$1
output_dir=$2

hadd $2/unfolding_histograms.root $1/*.root 
