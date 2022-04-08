# ZBosonMassSimHistCreator

This code is not written very efficiently, but I tried to at least make the most important parts easy to understand. It takes a while to run though, maybe 90-120 minutes. 

This code needs the ROOT data analysis framework to run. I wrote it using CMSSW 10.6.1.

The data for producing the histograms is located in the storage space accessible from t3.unl.edu. If you want to run this code and cannot access the data, let me know and I'll help you figure it out.

If you want to produce the histograms with different binning, simply run the ProduceHistograms.sh bash script. It takes three arguments in this order: number_of_bins, low_bin, high_bin. For example:
```
./ProduceHistograms.sh 150 50 200
```
This will make histograms with 150 bins starting at 50 GeV and up to 200 GeV.

Then this will run the analysis code to make the histograms with the number of bins and range provided in the options. It will take a while to run.

To make some plots to see the distributions, run MakePlots.sh with the subdirectory name of the root file as an option. For example:
```
./MakePlots.sh bins_150_low_50_high_200
```
Because the 150 bin samples are located in 'output_data/bins_150_low_50_high_200/'
