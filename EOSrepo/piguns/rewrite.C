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

void rewrite() {
  TString filename_in = "lctuple_gun.root";
  TString filename_out = "output_lctuple.root";

  TFile *file_in = new TFile(filename_in, "READ");
  TFile *file_out = new TFile(filename_out, "RECREATE");

  TDirectory *dir = (TDirectory *)file_in->Get("FilterDL_VXDE");

  TTree *tree_in = static_cast<TTree *>(dir->Get("MyLCTuple"));

  TTree *tree_out = tree_in->CloneTree();

  file_in->Close();
  file_in->Delete();

  tree_out->Write();
  file_out->Close();
  file_out->Delete();
}