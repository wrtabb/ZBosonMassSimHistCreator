# ZBosonMassSimHistCreator

I was needing to remake the plots again with different binning, and the length of time it was taking added to the fact that it had to run locally on my machine was frustrating. So I updated the package to instead run using Condor. 

One must have a certificate to run the code, but it takes far less time now to run. 

Previously, the code would run locally. It would load all ~900 root files and combine the trees into one TChain and then run sequentially over all events. Now a Condor job is submitted for each of the files and they are all run on a distributed network with many files being processed in parallel. This makes a tremendous difference in run time.

To run the package, just do:
```
./Submit.sh nbins lowbin highbin
```

When the code is finished, it outputs the results into root files in the 'output_data' directory. These files can be quickly combined using the CombineFiles.sh script.

Then the 'unfolding_histograms.root' file can be moved to a directory and the ~900 individual files can be deleted.

Where nbins is the number of bins, low bin is the value of the lower edge of the histogram and high bin is the value of the upper edge of the histogram.

I quickly made three different scenarios:
```
nbins = 150
lowbin = 50
highbin = 200
```
```
nbins = 50
lowbin = 20
highbin = 500
```
```
nbins = 50
lowbin = 50
highbin = 150
```

These are located in the subdirectories within 'output_data'

The matrix made with a much larger range looks more diagonal, so may be easier to solve. Also having wider bins helps with this as well. As can be seen from the case with 150 bins, we have a large blob at the z-boson mass peak with many off-diagonal elements. This may make it more difficult to solve.

If you have any issues running the code, feel free to ask any questions you have of me. I am also perfectly happy to run the code for you again to make distributions fitting whatever parameters you have in mind.
