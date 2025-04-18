#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Math/Vector4D.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

/*Code for decay modes:
0 = 1-prong
1 = 3-prong
2 = e
3 = mu
4 = other
*/

// ---- Histograms booking

TH1F *h_recD0Single, *h_recSigmaD0Single, *h_recD0overSigma;

int eta_classifier(float eta) {
  if (0. <= fabs(eta) && fabs(eta) <= 0.5)
    return 0;
  else if (0.5 <= fabs(eta) && fabs(eta) <= 1.)
    return 1;
  else if (1. <= fabs(eta) && fabs(eta) <= 1.5)
    return 2;
  // else if (1.5 <= fabs(eta) && fabs(eta) <= 2.)
  // return 3;
  else
    return 3;
}

void HistoRecoFiller(const TString &filename,
                     const TString &decay = "had",
                     const TString &boson = "Z") {
  bool corrections = true;

  float pt_cut = 20.;

  double D0min = 1.E-5;
  double D0max = 5.;
  int nlogbins_D0 = 50;
  double logbinedges_D0[nlogbins_D0 + 1];

  for (int i = 0; i <= nlogbins_D0; i++) {
    logbinedges_D0[i] = pow(
        10, TMath::Log10(D0min) + (TMath::Log10(D0max) - TMath::Log10(D0min)) /
                                      double(nlogbins_D0) * double(i));
  }

  double SigmaD0min = 1.E-6;
  double SigmaD0max = 1.E-3;
  int nlogbins_SigmaD0 = 50;
  double logbinedges_SigmaD0[nlogbins_SigmaD0 + 1];

  for (int i = 0; i <= nlogbins_SigmaD0; i++) {
    logbinedges_SigmaD0[i] =
        pow(10, TMath::Log10(SigmaD0min) +
                    (TMath::Log10(SigmaD0max) - TMath::Log10(SigmaD0min)) /
                        double(nlogbins_SigmaD0) * double(i));
  }

  double D0overSigmamin = 1.E-2;
  double D0overSigmamax = 4.E+2;
  int nlogbins_D0overSigma = 50;
  double logbinedges_D0overSigma[nlogbins_D0overSigma + 1];

  for (int i = 0; i <= nlogbins_D0overSigma; i++) {
    logbinedges_D0overSigma[i] = pow(
        10, TMath::Log10(D0overSigmamin) +
                (TMath::Log10(D0overSigmamax) - TMath::Log10(D0overSigmamin)) /
                    double(nlogbins_D0overSigma) * double(i));
  }

  // --- Histograms initialization
    h_recD0Single = new TH1F(
        "h_recD0Single",
        "Reco D0 of the single taus; d_{0} [mm]",
        nlogbins_D0, logbinedges_D0);
    h_recSigmaD0Single =
        new TH1F("h_recSigmaD0Single",
                 "Reco error on D0 of the single "
                         "taus; #sigma_{d_{0}} [mm]",
                 nlogbins_SigmaD0, logbinedges_SigmaD0);
    h_recD0overSigma =
        new TH1F("h_recD0overSigma",
                 "d0 Reco Significance of the "
                         "taus; d_{0}/#sigma_{d_{0}}",
                 nlogbins_D0overSigma, logbinedges_D0overSigma);

  // Open and read energy correction file
  TFile *ECorr_file = new TFile(
      "/eos/user/l/lvalla/MuColl/tauguns/histo_corrections.root",
      "READ");

  if (!ECorr_file->IsOpen())
    throw std::invalid_argument("ECorr filename not valid");

  TList *ECorrHistos_e = (TList *)ECorr_file->Get("e_corrections");
  TList *ECorrHistos_mu = (TList *)ECorr_file->Get("mu_corrections");
  TList *ECorrHistos_1pi = (TList *)ECorr_file->Get("1pi_corrections");
  TList *ECorrHistos_1pi0N =
      (TList *)ECorr_file->Get("1pi0N_corrections");
  TList *ECorrHistos_1piN = (TList *)ECorr_file->Get("1piN_corrections");
  /*TList *ECorrHistos_1pi2N = (TList *)ECorr_file->Get("1pi2N_corrections");
  TList *ECorrHistos_1pi3N = (TList *)ECorr_file->Get("1pi3N_corrections");*/
  TList *ECorrHistos_3pi = (TList *)ECorr_file->Get("3pi_corrections");
  TList *ECorrHistos_other = (TList *)ECorr_file->Get("other_corrections");

  TH1F *histos_e_array[4] = {
      static_cast<TH1F *>(ECorrHistos_e->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_e->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_e->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_e->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_mu_array[4] = {
      static_cast<TH1F *>(ECorrHistos_mu->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_mu->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_mu->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_mu->FindObject("h_ECorr_Eta3"))};
  /*TH1F *histos_1pi_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1pi->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1pi->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1pi->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1pi->FindObject("h_ECorr_Eta3"))};*/
  TH1F *histos_1pi0N_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1pi0N->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1pi0N->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1pi0N->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1pi0N->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_1piN_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1piN->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1piN->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1piN->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1piN->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_3pi_array[4] = {
      static_cast<TH1F *>(ECorrHistos_3pi->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_3pi->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_3pi->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_3pi->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_other_array[4] = {
      static_cast<TH1F *>(ECorrHistos_other->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_other->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_other->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_other->FindObject("h_ECorr_Eta3"))};

  // Open evaltaufinder output file
  TFile *input_file = new TFile(filename, "READ");

  if (!input_file->IsOpen())
    throw std::invalid_argument("filename not valid");

  TTree *inTree = (TTree *)input_file->Get("evtree");

  // --- MC particles
  int ntau_rec, evID, runID;

  int *mcDecayMode = new int[100000];
  int *charge = new int[100000];
  float *recPt = new float[100000];
  float *recEta = new float[100000];
  float *recD0 = new float[100000];
  float *recSigmaD0 = new float[100000];
  int *recNPfos = new int[500000];
  int *recNQTracks = new int[500000];
  int *pfosPdg = new int[500000];
  int *matched = new int[100000];

  inTree->SetBranchAddress("nTausRec", &ntau_rec);
  inTree->SetBranchAddress("EvID", &evID);
  inTree->SetBranchAddress("RunID", &runID);
  inTree->SetBranchAddress("mcDecayMode", mcDecayMode);
  inTree->SetBranchAddress("charge", charge);
  inTree->SetBranchAddress("recNPfos", recNPfos);
  inTree->SetBranchAddress("recNQTracks", recNQTracks);
  inTree->SetBranchAddress("recPt", recPt);
  inTree->SetBranchAddress("recD0", recD0);
  inTree->SetBranchAddress("recSigmaD0", recSigmaD0);
  inTree->SetBranchAddress("pfosPdg", pfosPdg);
  inTree->SetBranchAddress("Matched", matched);

  float recPtSingle[2] = {0., 0.}, recEtaSingle[2] = {0., 0.},
        recD0overSigma[2] = {0., 0.};
  int recDM[2] = {-1, -1};

  int recDM_0 = -1, recDM_1 = -1, recNPfos_0 = 0, recNPfos_1 = 0;

  const long int nEntries = inTree->GetEntries();

  std::cout << boson << " type entries: " << nEntries << '\n';

  for (int ientry = 0; ientry < nEntries; ++ientry) {

    inTree->GetEntry(ientry);

    float ECorr[3] = {1., 1.};
    int n_el[4] = {0, 0};
    int n_mu[4] = {0, 0};
    int n_pipl[4] = {0, 0};
    int n_pimi[4] = {0, 0};
    int n_gamma[4] = {0, 0};
    int n_neutrals[4] = {0, 0};
    int eta_range[4] = {0, 0};

    TH1F *histo_corr = nullptr;
    float bin_correction;
    TString recDecayMode[2];

    for (int t = 0, pfos_counter = 0; t != ntau_rec; ++t) {

      n_neutrals[t] = recNPfos[t] - recNQTracks[t];
      eta_range[t] = eta_classifier(recEta[t]);

      for (int j = 0; j != recNPfos[t]; j++) {
        if (abs(pfosPdg[j + pfos_counter]) == 11)
          n_el[t]++;
        else if (abs(pfosPdg[j + pfos_counter]) == 13)
          n_mu[t]++;
        else if (pfosPdg[j + pfos_counter] == 22)
          n_gamma[t]++;
        else if (pfosPdg[j + pfos_counter] == 211 ||
                 pfosPdg[j + pfos_counter] == 321)
          n_pipl[t]++;
        else if (pfosPdg[j + pfos_counter] == -211 ||
                 pfosPdg[j + pfos_counter] == -321)
          n_pimi[t]++;
      }
      pfos_counter += 20;

      /*if ((n_pimi[t] + n_pipl[t] == 1) && (recNQTracks[t] == 1)) {
        histo_corr = histos_1pi_array[eta_range[t]];
        recDecayMode[t] = "1pi";
      }*/
      if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] == 0)) {
        histo_corr = histos_1pi0N_array[eta_range[t]];
        recDecayMode[t] = "1pi0N";
      } else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] >= 1)) {
        histo_corr = histos_1piN_array[eta_range[t]];
        recDecayMode[t] = "1piN";
      } /*else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] == 2)) {
        histo_corr = histos_1pi2N_array[eta_range[t]];
        recDecayMode[t] = "1pi2N";
      } else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] >= 3)) {
        histo_corr = histos_1pi3N_array[eta_range[t]];
        recDecayMode[t] = "1pi3N";
      }*/
      else if ((n_pimi[t] + n_pipl[t] == 3) && (recNQTracks[t] == 3)) {
        histo_corr = histos_3pi_array[eta_range[t]];
        recDecayMode[t] = "3pi";
      } else if (n_el[t] == 1 && (recNQTracks[t] == 1)) {
        histo_corr = histos_e_array[eta_range[t]];
        recDecayMode[t] = "e";
      } else if (n_mu[t] == 1 && (recNQTracks[t] == 1)) {
        histo_corr = histos_mu_array[eta_range[t]];
        recDecayMode[t] = "mu";
      } else {
        histo_corr = histos_other_array[eta_range[t]];
        recDecayMode[t] = "other";
      }

      bin_correction = histo_corr->GetBinContent(histo_corr->FindBin(recPt[t]));
      bin_correction > 0. ? ECorr[t] = bin_correction : ECorr[t] = 1.;

      if(recDecayMode[t] == "1pi0N" || recDecayMode[t] == "1piN") recDecayMode[t] = "1pi";

      histo_corr = nullptr;

    if (corrections == false){
      ECorr[t] = 1.;
    }

    float recPtCorr = recPt[t] * ECorr[t];

    if (abs(charge[t]) != 1)
      continue;

    if (recPtCorr[t] < pt_cut)
      continue;

    if (recEta[t] > 2.)
      continue;
    
    if(decay == "had"){
      if(recDecayMode[t] != "1pi" && recDecayMode[t] != "3pi")
        continue;
    }

    h_recD0Single->Fill(std::fabs(recD0[t]));
    h_recSigmaD0Single->Fill(std::fabs(recSigmaD0[t]));
    recD0overSigma[t] = std::fabs(recD0[t] / recSigmaD0[t]);
    h_recD0overSigma->Fill(recD0overSigma[t]);

      pfos_counter += 20;
    }
  }
  inTree->Delete();
  input_file->Close();
  input_file->Delete();

  gStyle->SetOptStat(111111);
  h_recD0Single->Scale(1. / h_recD0Single->Integral());
  h_recSigmaD0Single->Scale(1. / h_recSigmaD0Single->Integral());
  h_recD0overSigma->Scale(1. / h_recD0overSigma->Integral());

  {
    delete[] recPt;
    delete[] recEta;
    delete[] recD0;
    delete[] recSigmaD0;
    h_recD0Single->Delete();
    h_recSigmaD0Single->Delete();
    h_recD0overSigma->Delete();
  }
}

void PlotD0Histo(TH1* h1, TH1* h2){
  TCanvas *c1 = new TCanvas("c1", "d0", 800, 600);

  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);

  // Customize histograms
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);

  display_overflow(h1);
  display_overflow(h2);

  // Draw histograms on the same canvas
  h1->Draw("HISTO");
  h2->Draw("HISTO SAME");

  h1->Scale(1. / h1->GetEntries());
  h2->Scale(1. / h2->GetEntries());

  h1->SetMaximum(h2->GetMaximum() * 1.1);

  h1->SetTitle(";d_{0}/#sigma_{d_{0}}; a.u.");

  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetXaxis()->SetLabelSize(0.045);
  h1->GetXaxis()->SetLabelOffset(0.012);
  h1->GetYaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetLabelSize(0.045);
  h1->GetYaxis()->SetLabelOffset(0.0125);

  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);
  c1->SetRightMargin(0.1);

  // Add legend
  TLegend *legend1 = new TLegend(0.62, 0.73, 0.89, 0.87);
  legend1->AddEntry(h1, "true #tau_{h}", "l");
  legend1->AddEntry(h2, "fake #tau_{h} (jets)", "l");
  legend1->SetBorderSize(0);
  legend1->SetTextSize(0.05);
  legend1->Draw();

  // Add TPaveText
  TPaveText *pt1 = new TPaveText(30, 0.056, 120, 0.078);
  pt1->SetFillColor(0);
  pt1->SetMargin(0);
  pt1->SetFillStyle(0);
  pt1->SetBorderSize(0);
  pt1->SetTextAlign(11);
  pt1->SetTextSize(0.045);
  pt1->AddText("Muon Collider");
  pt1->AddText("#bf{#it{Simulation}}");
  pt1->AddText("#bf{#sqrt{s} = 3 TeV}");
  //pt1->AddText("#bf{H #rightarrow #tau_{h}#tau_{h} sample}");
  // pt1->AddText("  ");
  //  pt1->AddText("#bf{}");
  pt1->Draw();

  // Save the canvas as an image
  c1->SaveAs("d0analysis.pdf");
}

void D0Analysis(
    const TString filenameH = "/eos/user/l/lvalla/MuColl/H2tautau10k/"
                              "output_EvalTauFinder_new.root",
    const TString filenameFake = "/eos/user/l/lvalla/MuColl/Z2jets/"
                                 "output_EvalTauFinder.root") {

  TH1* h_truetau = D0HistoFiller(filenameH, "had", "H");
  TH1* h_faketau = D0HistoFiller(filenameFake, "had", "Z");

  PlotD0Histo(h_truetau, h_faketau);
}