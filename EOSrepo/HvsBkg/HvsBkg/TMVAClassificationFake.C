

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Tools.h"

int TMVAClassificationFake(TString myMethodList = "") {

  TString fake_label = "MuMu";
  TString decay0 = "had";
  TString decay1 = "had";

  ROOT::EnableImplicitMT();
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string, int> Use;

  // Boosted Decision Trees
  Use["BDT"] = 1;  // uses Adaptive Boost
  Use["BDTG"] = 0; // uses Gradient Boost
  Use["BDTB"] = 0; // uses Bagging
  Use["BDTD"] = 0; // decorrelation + Adaptive Boost
  Use["BDTF"] = 0; // allow usage of fisher discriminant for node splitting
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end();
         it++)
      it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');
    for (UInt_t i = 0; i < mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod
                  << "\" not known in TMVA under this name. Choose among the "
                     "following:"
                  << std::endl;
        for (std::map<std::string, int>::iterator it = Use.begin();
             it != Use.end(); it++)
          std::cout << it->first << " ";
        std::cout << std::endl;
        return 1;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // Here the preparation phase begins

  // Read training and test data
  // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
  TFile *input(0);
  TString fname = "./Hvs" + fake_label + "_Histograms.root";
  input = TFile::Open(fname); // check if file in local directory exists

  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- TMVAClassification       : Using input file: "
            << input->GetName() << std::endl;

  // Register the training and test trees

  // Register the training and test trees
  TTree *signalTree = (TTree *)input->Get("Htree_" + decay0 + "+" + decay1);
  TTree *ttmmbkgTree = (TTree *)input->Get("Ztree_" + decay0 + "+" + decay1);

  Int_t signalMCEntries =
      static_cast<TH1F *>(
          input->Get("H_MC_histograms")->FindObject("h_mcMInv"))
          ->GetEntries();
  Int_t ttmmMCEntries =
      static_cast<TH1F *>(
          input->Get("Z_MC_histograms")->FindObject("h_mcMInv"))
          ->GetEntries();

  Int_t signalRecoEntries = signalTree->GetEntries();
  Int_t ttmmRecoEntries = ttmmbkgTree->GetEntries();

  Double_t intlumi = 1.E3;    // in fb-1
  Double_t signalXsec = 31.5; // in fb
  Double_t ttmmbkgXsec = 185.3;

  Double_t signalEff = Double_t(signalRecoEntries) / Double_t(signalMCEntries);
  Double_t ttmmbkgEff = Double_t(ttmmRecoEntries) / Double_t(ttmmMCEntries);

  std::cout << "Signal efficiency: " << signalEff << '\n';
  std::cout << "mumu->tautaumumu efficiency: " << ttmmbkgEff << '\n';

  Double_t n_signalExp = intlumi * signalXsec * signalEff;
  Double_t n_ttmmbkgExp = intlumi * ttmmbkgXsec * ttmmbkgEff;

  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight = Double_t(n_signalExp) / Double_t(signalRecoEntries);
  Double_t ttmmbkgWeight = Double_t(n_ttmmbkgExp) / Double_t(ttmmRecoEntries);

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName("TMVAC" + fake_label + ".root");
  TFile *outputFile = TFile::Open(outfileName, "RECREATE");

  TMVA::Factory *factory =
      new TMVA::Factory("TMVAClassification", outputFile,
                        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;"
                        "P;G,D:AnalysisType=Classification");

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  dataloader->AddVariable("recMInv", "#tau#tau invariant Mass", "GeV", 'F');
  dataloader->AddVariable("recDeltaR", "#tau#tau #DeltaR", "", 'F');
  dataloader->AddVariable("recPtPair", "#tau#tau system pt", "GeV/c", 'F');
  dataloader->AddVariable("recEtaPair", "#tau#tau system eta", "", 'F');
  dataloader->AddVariable("recDeltaPhi", "#tau#tau #Delta#Phi", "", 'F');
  dataloader->AddVariable("recDeltaEta", "#tau#tau #Delta#eta", "", 'F');
  dataloader->AddVariable("recSepAngle", "#tau#tau Separation Angle", "", 'F');
  dataloader->AddVariable("recPtSingle_0", "1st #tau pt", "GeV/c", 'F');
  dataloader->AddVariable("recPtSingle_1", "2nd #tau pt", "GeV/c", 'F');
  //dataloader->AddVariable("recE_0", "1st #tau E", "GeV", 'F');
  //dataloader->AddVariable("recE_1", "2nd #tau E", "GeV", 'F');
  //dataloader->AddVariable("recEtaSingle_0", "1st #tau #eta", "", 'F');
  //dataloader->AddVariable("recEtaSingle_1", "2nd #tau #eta", "", 'F');
  //dataloader->AddVariable("recDM_0", "Reco Decay Mode 1st tau", "", 'I');
  //dataloader->AddVariable("recDM_1", "Reco Decay Mode 2nd tau", "", 'I');
  dataloader->AddVariable("recChAsymm_0", "1st #tau Charge Asymmetry", "", 'F');
  dataloader->AddVariable("recChAsymm_1", "2nd #tau Charge Asymmetry", "", 'F');
  dataloader->AddVariable("recD0overSigma_0", "1st #tau d_{0}/#sigma_{d_{0}}",
                          "", 'F');
  dataloader->AddVariable("recD0overSigma_1", "2nd #tau d_{0}/#sigma_{d_{0}}",
                          "", 'F');
  //dataloader->AddVariable("recNPfos_0", "Number pfos 1st tau", "", 'I');
  //dataloader->AddVariable("recNPfos_1", "Number pfos 2nd tau", "", 'I');

  // You can add so-called "Spectator variables", which are not used in the MVA
  // training, but will appear in the final "TestTree" produced by TMVA. This
  // TestTree will contain the input variables, the response values of all
  // trained MVAs, and the spectator variables

  // dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
  // dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

  // You can add an arbitrary number of signal or bkgTree trees
  dataloader->AddSignalTree(signalTree, signalWeight);
  dataloader->AddBackgroundTree(ttmmbkgTree, ttmmbkgWeight);

  // To give different trees for training and testing, do as follows:
  //
  //     dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight,
  //     "Training" ); dataloader->AddSignalTree( signalTestTree,
  //     signalTestWeight,  "Test" );

  // Use the following code instead of the above two or four lines to add signal
  // and bkgTree training and test events "by hand" NOTE that in this case
  // one should not give expressions (such as "var1+var2") in the input
  //      variable definition, but simply compute the expression before adding
  //      the event
  // ```cpp
  // // --- begin ----------------------------------------------------------
  // std::vector<Double_t> vars( 4 ); // vector has size of number of input
  // variables Float_t  treevars[4], weight;
  //
  // // Signal
  // for (UInt_t ivar=0; ivar<4; ivar++) signalTree->SetBranchAddress( Form(
  // "var%i", ivar+1 ), &(treevars[ivar]) ); for (UInt_t i=0;
  // i<signalTree->GetEntries(); i++) {
  //    signalTree->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is
  //    testing
  //    // note that the weight can also be event-wise
  //    if (i < signalTree->GetEntries()/2.0)
  //    dataloader->AddSignalTrainingEvent( vars, signalWeight ); else
  //    dataloader->AddSignalTestEvent    ( vars, signalWeight );
  // }
  //
  // // bkgTree (has event weights)
  // bkgTree->SetBranchAddress( "weight", &weight );
  // for (UInt_t ivar=0; ivar<4; ivar++) bkgTree->SetBranchAddress( Form(
  // "var%i", ivar+1 ), &(treevars[ivar]) ); for (UInt_t i=0;
  // i<bkgTree->GetEntries(); i++) {
  //    bkgTree->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is
  //    testing
  //    // note that the weight can also be event-wise
  //    if (i < bkgTree->GetEntries()/2)
  //    dataloader->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
  //    else                                dataloader->AddBackgroundTestEvent
  //    ( vars, backgroundWeight*weight );
  // }
  // // --- end ------------------------------------------------------------
  // ```
  // End of tree registration

  // Set individual event weights (the variables must exist in the original
  // TTree)
  // -  for signal    : `dataloader->SetSignalWeightExpression
  // ("weight1*weight2");`
  // -  for bkgTree:
  // `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
  // dataloader->SetBackgroundWeightExpression("weight");

  // Apply additional cuts on the signal and bkgTree samples (can be
  // different)
  TCut mycuts =
      ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  // Tell the dataloader how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used
  // for training, and the other half for testing:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  //
  // To also specify the number of testing events, use:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut,
  //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V"
  //         );
  dataloader->PrepareTrainingAndTestTree(
      mycuts, mycutb, "SplitMode=Random:NormMode=NumEvents:!V");

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod(
        dataloader, TMVA::Types::kBDT, "BDTG",
        "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:"
        "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

  if (Use["BDT"]) // Adaptive Boost
    factory->BookMethod(
        dataloader, TMVA::Types::kBDT, "BDT",
        "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:"
        "AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:"
        "SeparationType=GiniIndex:nCuts=20");

  if (Use["BDTB"]) // Bagging
    factory->BookMethod(
        dataloader, TMVA::Types::kBDT, "BDTB",
        "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod(
        dataloader, TMVA::Types::kBDT, "BDTD",
        "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:"
        "SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

  if (Use["BDTF"]) // Allow Using Fisher discriminant in node splitting for
                   // (strong) linearly correlated variables
    factory->BookMethod(
        dataloader, TMVA::Types::kBDT, "BDTF",
        "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType="
        "AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

  // For an example of the category classifier usage, see:
  // TMVAClassificationCategory
  //
  // --------------------------------------------------------------------------------------------------
  //  Now you can optimize the setting (configuration) of the MVAs using the set
  //  of training events
  // STILL EXPERIMENTAL and only implemented for BDT's !
  //
  //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
  //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
  //
  // --------------------------------------------------------------------------------------------------

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
  delete dataloader;
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch())
    TMVA::TMVAGui(outfileName);

  return 0;
}

int main(int argc, char **argv) {
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i = 1; i < argc; i++) {
    TString regMethod(argv[i]);
    if (regMethod == "-b" || regMethod == "--batch")
      continue;
    if (!methodList.IsNull())
      methodList += TString(",");
    methodList += regMethod;
  }
  return TMVAClassificationFake(methodList);
}
