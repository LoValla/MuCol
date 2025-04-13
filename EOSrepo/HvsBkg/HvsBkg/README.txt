Folder that contains everything that is needed to analyse different variables of the H->tautau, Z->tautau and mumu->tautaununu events.
The Histograms.root files contain all the details of the events of the different categories (signal and different backgrounds).
One need to run the TauPairsAnalysis.cxx file, followed by TauPairsPlot.cxx to plot many different variables of the tau pairs for the different categories.
One can also train and test a BDT for classification, running the TMVAClassificationTot.C file.
Finally, by running fit_TMVATot.C, a combined histogram fit for bkg+signal can be performed, to extract an estimate of the signal events and the associated statistical uncertainty.