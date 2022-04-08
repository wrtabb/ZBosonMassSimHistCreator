#!/bin/bash

nbins=$1
low_bin=$2
high_bin=$3

root -l << EOF
    .L produceHistograms.C
    produceHistograms($nbins,$low_bin,$high_bin)
EOF

