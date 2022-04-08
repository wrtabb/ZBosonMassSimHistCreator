# ZBosonMassSimHistCreator

This code is not written very efficiently, but I tried to at least make the most important parts easy to understand. It takes a while to run though, maybe 90-120 minutes. 

This code needs the ROOT data analysis framework to run. I wrote it using CMSSW 10.6.1.

The data for producing the histograms is located in the storage space accessible from t3.unl.edu. If you want to run this code and cannot access the data, let me know and I'll help you figure it out.

If you want to produce the histograms with different binning, simply run the ProduceHistograms.sh bash script. It takes three arguments in this order: number_of_bins, low_bin, high_bin:
```
./ProduceHistograms.sh number_of_bins low_bin high_bin
```

Then this will run the analysis code to make the histograms with the number of bins and range provided in the options. It will take a while to run.

To make some plots to see the distributions, run 'makePlots.C' using ROOT:
```
root -l makePlots.C
```
