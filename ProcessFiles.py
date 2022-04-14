#!/usr/bin/env python
  
import sys
import subprocess

if (len(sys.argv)!=5) :
    print("ERROR: this script expects exactly 4 arguments")
    sys.exit(0)

fileName    = sys.argv[1]
nbins       = sys.argv[2]
lowbin      = sys.argv[3]
highbin     = sys.argv[4]

print("Input file: "+fileName)
print("number of bins: "+nbins)
print("lower histogram boundary: "+lowbin)
print("upper histogram boundary: "+highbin)

scriptSpecs = 'produceHistograms.C+("'+str(fileName)+'",'+str(nbins)+','+str(lowbin)+','+str(highbin)+')'
print(scriptSpecs)
rootCommand = ['root']
rootCommand.append('-l')
rootCommand.append('-b')
rootCommand.append('-q')
rootCommand.append(scriptSpecs)

(out,err) = subprocess.Popen(rootCommand,stdout=subprocess.PIPE).communicate()
print(out)
print("Any errors?")
print(err)
