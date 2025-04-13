#include <fstream>
#include <iostream>

#include "Math/Vector4D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

const int nRecDecayModes = 8;

// histogram declaration

TH1I *h_mcDecayMode;
TH2I *h_DecayMode_McvsRec;
TH2I *h_DecayMode_RecvsMc;

void TauDMAnalysis(const TString filename = "output_EvalTauFinder.root") {

  // histogram initialization
  h_DecayMode_McvsRec = new TH2I(
      "h_DecayMode_McvsRec", ";MC #tau decay mode; Reco #tau decay mode", 9,
      -0.5, 8.5, nRecDecayModes, -0.5, nRecDecayModes - 0.5);
  h_DecayMode_RecvsMc = new TH2I(
      "h_DecayMode_RecvsMc", ";Reco #tau decay mode;MC #tau decay mode",
      nRecDecayModes, -0.5, nRecDecayModes - 0.5, 9, -0.5, 8.5);
  h_mcDecayMode = new TH1I(
      "h_mcDecayMode", "MC #tau decay mode; Decay Mode; Entries", 9, -0.5, 8.5);
  //  --- Open the ntuple file and get the tree

  TFile *input_file = new TFile(filename, "update");

  if (!input_file->IsOpen())
    throw std::invalid_argument("filename not valid");

  TTree *tree = (TTree *)input_file->Get("evtree");
  // TTree *match_tree = (TTree *)input_file->Get("taumatch");
  // TTree *miss_tree = (TTree *)input_file->Get("mcmiss");
  // TTree *fake_tree = (TTree *)input_file->Get("fake");

  // --- MC particles
  int ntau_mc;
  int ntau_rec;

  int evID;
  int runID;

  // int *evt_id = new int[10000];
  int *mcDecayMode = new int[10000000];
  int *tau_charge = new int[10000000];
  int *recNPfos = new int[10000000];
  int *recNQTracks = new int[10000000];
  int *pfosPdg = new int[50000000];

  tree->SetBranchAddress("nTausMC", &ntau_mc);
  tree->SetBranchAddress("nTausRec", &ntau_rec);
  tree->SetBranchAddress("EvID", &evID);
  tree->SetBranchAddress("RunID", &runID);
  tree->SetBranchAddress("mcDecayMode", mcDecayMode);
  tree->SetBranchAddress("recNPfos", recNPfos);
  tree->SetBranchAddress("recNQTracks", recNQTracks);
  tree->SetBranchAddress("charge", tau_charge);
  tree->SetBranchAddress("pfosPdg", pfosPdg);

  const long int nEntries = tree->GetEntries();

  TString plotsfolder = "plots";

  for (int ientry = 0; ientry < nEntries; ++ientry) {
    if (ientry % 10000 == 0)
      std::cout << ientry << " / " << nEntries << '\n';

    tree->GetEntry(ientry);

    if (ntau_mc != 1)
      continue;

    if (mcDecayMode[0] == 4)
      mcDecayMode[0] = 3;
    if (mcDecayMode[0] == 6)
      mcDecayMode[0] = 4;
    if (mcDecayMode[0] == 7)
      mcDecayMode[0] = 5;
    if (mcDecayMode[0] == 8)
      mcDecayMode[0] = 6;

    h_mcDecayMode->Fill(mcDecayMode[0]);

    if (ntau_rec != 1) {
      h_DecayMode_McvsRec->Fill(mcDecayMode[0], 6);
      h_DecayMode_RecvsMc->Fill(6, mcDecayMode[0]);
      continue;
    }

    if (abs(tau_charge[0]) != 1) {
      h_DecayMode_McvsRec->Fill(mcDecayMode[0], 6);
      h_DecayMode_RecvsMc->Fill(6, mcDecayMode[0]);
      continue;
    }

    


    // if (abs(tau_charge[0]) != 1) continue;

    int n_el = 0;
    int n_mu = 0;
    int n_pipl = 0;
    int n_pimi = 0;
    int n_gamma = 0;

    int recDecayMode = 5;

    for (int j = 0; j != recNPfos[0]; j++) {
      if (abs(pfosPdg[j]) == 11)
        n_el++;
      else if (abs(pfosPdg[j]) == 13)
        n_mu++;
      else if (pfosPdg[j] == 22)
        n_gamma++;
      else if (pfosPdg[j] == 211 || pfosPdg[j] == 321)
        n_pipl++;
      else if (pfosPdg[j] == -211 || pfosPdg[j] == -321)
        n_pimi++;
    }

    if ((n_pipl + n_pimi) == 1 && recNPfos[0] == 1)
      recDecayMode = 0;
    else if ((n_pipl + n_pimi) == 1 && (recNPfos[0] - recNQTracks[0]) == 1)
      recDecayMode = 1;
    else if ((n_pipl + n_pimi) == 1 && (recNPfos[0] - recNQTracks[0]) == 2)
      recDecayMode = 2;
    else if ((n_pipl + n_pimi) == 1 && (recNPfos[0] - recNQTracks[0]) >= 3)
      recDecayMode = 3;
    else if (n_pipl + n_pimi == 3)
      recDecayMode = 4;
    else if (n_el == 1 && recNPfos[0] == 1)
      recDecayMode = 5;
    else if (n_mu == 1 && recNPfos[0] == 1)
      recDecayMode = 5;
    else
      recDecayMode = 5;

    h_DecayMode_McvsRec->Fill(mcDecayMode[0], recDecayMode);
    h_DecayMode_RecvsMc->Fill(recDecayMode, mcDecayMode[0]);
  }

  for (int col = 1; col != 9 + 1; col++) {
    int tot_col_events = h_mcDecayMode->GetBinContent(col);
    std::cout << "col" << col << ": " << tot_col_events << '\n';
    if (tot_col_events == 0)
      continue;

    for (int row = 1; row != nRecDecayModes + 1; row++) {
      int bin_events = h_DecayMode_McvsRec->GetBinContent(col, row);

      std::cout << "float(bin_events) / float(tot_col_events) = "
                << float(bin_events) / float(tot_col_events) << '\n';
      h_DecayMode_McvsRec->SetBinContent(
          col, row, float(bin_events) / float(tot_col_events) * 1000);
    }
  }

  {
    const char *rec_labels[nRecDecayModes] = {
        "h^{#pm}",  "h^{#pm} + N", "h^{#pm} + 2N", "h^{#pm} + #geq3N",
        "3h^{#pm}", "Others",  "None",   "#mu^{#pm}"};
    const char *mc_labels[9] = {"h^{#pm}",
                                "h^{#pm} + #pi^{0}",
                                "h^{#pm} + 2#pi^{0}",
                                "3h^{#pm}",
                                "e^{#pm}",
                                "#mu^{#pm}",
                                "h^{#pm} + 3#pi^{0}",
                                "Others",
                                "3h^{#pm} + #pi^{0}",
                                };

    auto c_mcDecayMode =
        new TCanvas("c_mcDecayMode", "c_mcDecayMode", 900, 800);
    h_mcDecayMode->Draw();
    // h_mcDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
    for (int i = 1; i != 9 + 1; i++)
      h_mcDecayMode->GetXaxis()->SetBinLabel(i, mc_labels[i - 1]);
    h_mcDecayMode->GetXaxis()->SetLabelSize(0.045);
    c_mcDecayMode->SaveAs("plots/c_mcDecayMode.png");
    c_mcDecayMode->Close();

    auto c_DecayMode_McvsRec =
        new TCanvas("c_DecayMode_McvsRec", "c_DecayMode_McvsRec", 1000, 800);
    gStyle->SetOptStat(0);
    h_DecayMode_McvsRec->Draw("COL TEXT");
    // h_DecayMode_McvsRec->GetXaxis()->SetBit(TAxis::kLabelsHori);
    h_DecayMode_McvsRec->GetYaxis()->SetRange(1, nRecDecayModes - 1);
    h_DecayMode_McvsRec->GetYaxis()->SetTitle("Reconstructed decay mode");
    h_DecayMode_McvsRec->GetXaxis()->SetRange(1, 4);
    for (unsigned int i = 1; i <= 9; i++)
      h_DecayMode_McvsRec->GetXaxis()->SetBinLabel(i, mc_labels[i - 1]);
    for (unsigned int j = 1; j <= nRecDecayModes; j++)
      h_DecayMode_McvsRec->GetYaxis()->SetBinLabel(j, rec_labels[j - 1]);
    h_DecayMode_McvsRec->GetXaxis()->SetLabelSize(0.045);
    h_DecayMode_McvsRec->GetYaxis()->SetLabelSize(0.045);
    h_DecayMode_McvsRec->GetXaxis()->SetTitleOffset(1.15);
    c_DecayMode_McvsRec->SaveAs("plots/c_DecayMode_McvsRec.png");
    c_DecayMode_McvsRec->Close();
  }

  delete[] tau_charge;
  delete[] mcDecayMode;
  delete[] recNPfos;
  delete[] recNQTracks;
  delete[] pfosPdg;

  h_DecayMode_McvsRec->Delete();
  h_DecayMode_RecvsMc->Delete();
  h_mcDecayMode->Delete();

  tree->Delete();
  input_file->Close();
  input_file->Delete();
}