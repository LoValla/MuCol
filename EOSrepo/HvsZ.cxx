#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "Math/Vector4D.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// ---- Histograms booking
std::vector<TH1F *> Z_histo_maker() {
  // --- Z Histograms initialization
  TH1F *hZ_recMInv = new TH1F(
      "hZ_recMInv", "Invariant mass of tau pairs; M_{inv} [GeV]; Events", 50,
      0., 150.);
  TH1F *hZ_mcMInv = new TH1F(
      "hZ_mcMInv", "MC Invariant mass of tau pairs; M_{inv} [GeV]; Events", 100,
      50., 130.);

  TH1F *hZ_recPtPair =
      new TH1F("hZ_recPtPair", "p_{T} of reco tau pairs; p_{T} [GeV]; Events",
               70, 0., 500.);
  TH1F *hZ_mcPtPair =
      new TH1F("hZ_mcPtPair", "p_{T} of MC tau pairs; p_{T} [GeV]; Events", 70,
               0., 500.);

  TH1F *hZ_recEtaPair = new TH1F(
      "hZ_recEtaPair", "#eta of reco tau pairs; #eta; Events", 40, -3, 3);
  TH1F *hZ_mcEtaPair = new TH1F(
      "hZ_mcEtaPair", "#eta of MC tau pairs; #eta; Events", 40, -3., 3.);

  TH1F *hZ_recPtSingle =
      new TH1F("hZ_recPtSingle",
               "p_{T} of Single reco taus; p_{T} [GeV]; Events", 80, 0., 500.);
  TH1F *hZ_mcPtSingle =
      new TH1F("hZ_mcPtSingle", "p_{T} of Single MC taus; p_{T} [GeV]; Events",
               80, 0., 500.);

  TH1F *hZ_recDeltaPhi = new TH1F(
      "hZ_recDeltaPhi", "#Delta#Phi between the two taus; |#Delta#Phi|; Events",
      40, 0., 180.);
  TH1F *hZ_mcDeltaPhi = new TH1F(
      "hZ_mcDeltaPhi",
      "#Delta#Phi between the two MC taus; |#Delta#Phi|; Events", 40, 0., 180.);

  TH1F *hZ_recDeltaR =
      new TH1F("hZ_recDeltaR", "#DeltaR between the two taus; #DeltaR; Events",
               50, 0., 7.5);
  TH1F *hZ_mcDeltaR =
      new TH1F("hZ_mcDeltaR",
               "#DeltaR between the two MC taus; #DeltaR; Events", 50, 0., 7.5);

  TH1F *hZ_recSepAngle = new TH1F(
      "hZ_recSepAngle", "Angle between the two rec taus; #alpha [rad]; Events",
      50, 0., 3.1416);
  TH1F *hZ_mcSepAngle = new TH1F(
      "hZ_mcSepAngle", "Angle between the two MC taus; #alpha [rad]; Events",
      50, 0., 3.1416);

  return {hZ_recMInv,     hZ_mcMInv,     hZ_recPtPair,   hZ_mcPtPair,
          hZ_recEtaPair,  hZ_mcEtaPair,  hZ_recPtSingle, hZ_mcPtSingle,
          hZ_recDeltaPhi, hZ_mcDeltaPhi, hZ_recDeltaR,   hZ_mcDeltaR,
          hZ_recSepAngle, hZ_mcSepAngle};
}

std::vector<TH1F *> H_histo_maker() {
  // --- H Histograms initialization
  TH1F *hH_recMInv = new TH1F(
      "hH_recMInv", "Invariant mass of tau pairs; M_{inv} [GeV]; Events", 50,
      0., 150.);
  TH1F *hH_mcMInv = new TH1F(
      "hH_mcMInv", "MC Invariant mass of tau pairs; M_{inv} [GeV]; Events", 100,
      50., 130.);

  TH1F *hH_recPtPair =
      new TH1F("hH_recPtPair", "p_{T} of reco tau pairs; p_{T} [GeV]; Events",
               70, 0., 500.);
  TH1F *hH_mcPtPair =
      new TH1F("hH_mcPtPair", "p_{T} of MC tau pairs; p_{T} [GeV]; Events", 50,
               0., 500.);
  TH1F *hH_recEtaPair = new TH1F(
      "hH_recEtaPair", "#eta of reco tau pairs; #eta; Events", 30, -3, 3);
  TH1F *hH_mcEtaPair =
      new TH1F("hH_mcEtaPair", "#eta of MC tau pairs; #eta; Events", 30, -3, 3);

  TH1F *hH_recPtSingle =
      new TH1F("hH_recPtSingle",
               "p_{T} of Single reco taus; p_{T} [GeV]; Events", 80, 0., 500.);
  TH1F *hH_mcPtSingle =
      new TH1F("hH_mcPtSingle", "p_{T} of Single MC taus; p_{T} [GeV]; Events",
               80, 0., 500.);

  TH1F *hH_recDeltaPhi = new TH1F(
      "hH_recDeltaPhi", "#Delta#Phi between the two taus; |#Delta#Phi|; Events",
      40, 0., 180.);
  TH1F *hH_mcDeltaPhi = new TH1F(
      "hH_mcDeltaPhi",
      "#Delta#Phi between the two MC taus; |#Delta#Phi|; Events", 40, 0., 180.);

  TH1F *hH_recDeltaR =
      new TH1F("hH_recDeltaR", "#DeltaR between the two taus; #DeltaR; Events",
               50, 0., 7.5);
  TH1F *hH_mcDeltaR =
      new TH1F("hH_mcDeltaR",
               "#DeltaR between the two MC taus; #DeltaR; Events", 50, 0., 7.5);

  TH1F *hH_recSepAngle = new TH1F(
      "hH_recSepAngle", "Angle between the two rec taus; #alpha [rad]; Events",
      50, 0., 3.1416);
  TH1F *hH_mcSepAngle = new TH1F(
      "hH_mcSepAngle", "Angle between the two MC taus; #alpha [rad]; Events",
      50, 0., 3.1416);

  return {hH_recMInv,     hH_mcMInv,     hH_recPtPair,   hH_mcPtPair,
          hH_recEtaPair,  hH_mcEtaPair,  hH_recPtSingle, hH_mcPtSingle,
          hH_recDeltaPhi, hH_mcDeltaPhi, hH_recDeltaR,   hH_mcDeltaR,
          hH_recSepAngle, hH_mcSepAngle};
}

int eta_classifier(float eta) {
  if (0. <= fabs(eta) && fabs(eta) <= 0.5)
    return 0;
  else if (0.5 <= fabs(eta) && fabs(eta) <= 1.)
    return 1;
  else if (1. <= fabs(eta) && fabs(eta) <= 1.5)
    return 2;
  //else if (1.5 <= fabs(eta) && fabs(eta) <= 2.)
    //return 3;
  else
    return 3;
}

void Histo_filler(std::vector<TH1F *> histos_vec, const TString &filename) {

  // Open and read energy correction file
  TFile *ECorr_file = new TFile(
      "/eos/user/l/lvalla/MuColl/tauguns/histo_corrections.root",
      "READ");

  ECorr_file->ls();

  if (!ECorr_file->IsOpen())
    throw std::invalid_argument("ECorr filename not valid");

  TList *ECorrHistos_e = (TList *)ECorr_file->Get("e_corrections");
  TList *ECorrHistos_mu = (TList *)ECorr_file->Get("mu_corrections");
  TList *ECorrHistos_1prong0NP =
      (TList *)ECorr_file->Get("1prong0NP_corrections");
  TList *ECorrHistos_1prong1NP =
      (TList *)ECorr_file->Get("1prong1NP_corrections");
  TList *ECorrHistos_1prong2NP =
      (TList *)ECorr_file->Get("1prong2NP_corrections");
  TList *ECorrHistos_1prong3NP =
      (TList *)ECorr_file->Get("1prong3NP_corrections");
  TList *ECorrHistos_3prong = (TList *)ECorr_file->Get("3prong_corrections");
  TList *ECorrHistos_others = (TList *)ECorr_file->Get("others_corrections");

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
  TH1F *histos_1prong0NP_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1prong0NP->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1prong0NP->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1prong0NP->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1prong0NP->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_1prong1NP_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1prong1NP->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1prong1NP->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1prong1NP->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1prong1NP->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_1prong2NP_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1prong2NP->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1prong2NP->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1prong2NP->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1prong2NP->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_1prong3NP_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1prong3NP->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1prong3NP->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1prong3NP->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1prong3NP->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_3prong_array[4] = {
      static_cast<TH1F *>(ECorrHistos_3prong->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_3prong->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_3prong->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_3prong->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_others_array[4] = {
      static_cast<TH1F *>(ECorrHistos_others->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_others->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_others->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_others->FindObject("h_ECorr_Eta3"))};

  TH1F *&h_recMInv = histos_vec[0];
  TH1F *&h_mcMInv = histos_vec[1];
  TH1F *&h_recPtPair = histos_vec[2];
  TH1F *&h_mcPtPair = histos_vec[3];
  TH1F *&h_recEtaPair = histos_vec[4];
  TH1F *&h_mcEtaPair = histos_vec[5];
  TH1F *&h_recPtSingle = histos_vec[6];
  TH1F *&h_mcPtSingle = histos_vec[7];
  TH1F *&h_recDeltaPhi = histos_vec[8];
  TH1F *&h_mcDeltaPhi = histos_vec[9];
  TH1F *&h_recDeltaR = histos_vec[10];
  TH1F *&h_mcDeltaR = histos_vec[11];
  TH1F *&h_recSepAngle = histos_vec[12];
  TH1F *&h_mcSepAngle = histos_vec[13];

  TFile *input_file = new TFile(filename, "READ");

  if (!input_file->IsOpen())
    throw std::invalid_argument("filename not valid");

  TTree *tree = (TTree *)input_file->Get("evtree");

  // --- MC particles
  int ntau_rec, evID, runID;

  int *mcDecayMode = new int[100000];

  float *mcE = new float[100000];
  float *mcPx = new float[100000];
  float *mcPy = new float[100000];
  float *mcPz = new float[100000];
  float *mcPt = new float[100000];
  float *mcEta = new float[100000];
  float *mcPhi = new float[100000];

  float *mcE_vis = new float[100000];
  float *mcPx_vis = new float[100000];
  float *mcPy_vis = new float[100000];
  float *mcPz_vis = new float[100000];
  float *mcPt_vis = new float[100000];
  float *mcEta_vis = new float[100000];
  float *mcPhi_vis = new float[100000];
  float *recE = new float[100000];
  float *recPx = new float[100000];
  float *recPy = new float[100000];
  float *recPz = new float[100000];
  float *recPt = new float[100000];
  float *recEta = new float[100000];
  float *recPhi = new float[100000];
  int *recNPfos = new int[500000];
  int *recNQTracks = new int[500000];
  int *pfosPdg = new int[500000];
  int *matched = new int[100000];

  tree->SetBranchAddress("nTausRec", &ntau_rec);
  tree->SetBranchAddress("EvID", &evID);
  tree->SetBranchAddress("RunID", &runID);
  tree->SetBranchAddress("mcE", mcE);
  tree->SetBranchAddress("mcPx", mcPx);
  tree->SetBranchAddress("mcPy", mcPy);
  tree->SetBranchAddress("mcPz", mcPz);
  tree->SetBranchAddress("mcPt", mcPt);
  tree->SetBranchAddress("mcEta", mcEta);
  tree->SetBranchAddress("mcPhi", mcPhi);
  tree->SetBranchAddress("mcDecayMode", mcDecayMode);
  tree->SetBranchAddress("mcE_vis", mcE_vis);
  tree->SetBranchAddress("mcPx_vis", mcPx_vis);
  tree->SetBranchAddress("mcPy_vis", mcPy_vis);
  tree->SetBranchAddress("mcPz_vis", mcPz_vis);
  tree->SetBranchAddress("mcPt_vis", mcPz_vis);
  tree->SetBranchAddress("mcEta_vis", mcEta_vis);
  tree->SetBranchAddress("mcPhi_vis", mcPhi_vis);
  tree->SetBranchAddress("recE", recE);
  tree->SetBranchAddress("recPx", recPx);
  tree->SetBranchAddress("recPy", recPy);
  tree->SetBranchAddress("recPz", recPz);
  tree->SetBranchAddress("recPt", recPt);
  tree->SetBranchAddress("recEta", recEta);
  tree->SetBranchAddress("recPhi", recPhi);
  tree->SetBranchAddress("recNPfos", recNPfos);
  tree->SetBranchAddress("recNQTracks", recNQTracks);
  tree->SetBranchAddress("pfosPdg", pfosPdg);
  tree->SetBranchAddress("Matched", matched);

  const long int nEntries = tree->GetEntries();

  float mcPxPair = 0., mcPyPair = 0., mcPzPair = 0., mcPtPair = 0., mcMInv = 0.,
        mcDeltaPhi = 0., mcDeltaR = 0., mcSepAngle = 0., mcEtaPair = 0.;
  float recPxPair = 0., recPyPair = 0., recPzPair = 0., recPtPair = 0.,
        recMInv = 0., recDeltaPhi = 0., recDeltaR = 0., recSepAngle = 0.,
        recEtaPair = 0.;

  for (int ientry = 0; ientry < nEntries; ++ientry) {

    tree->GetEntry(ientry);

    /*if (mcDecayMode[0] == -1 || mcDecayMode[1] == -1) {
      std::cout << "Run: " << runID << " - "
                << "Ev: " << evID << '\n';
    }*/

    h_mcPtSingle->Fill(mcPt[0]);
    h_mcPtSingle->Fill(mcPt[1]);

    // std::cout << "mcPx[0] = " << mcPx[0] << '\n';
    // std::cout << "mcPx[1] = " << mcPx[1] << '\n';

    mcPxPair = mcPx[0] + mcPx[1];
    mcPyPair = mcPy[0] + mcPy[1];
    mcPzPair = mcPz[0] + mcPz[1];

    mcPtPair = std::sqrt((mcPx[0] + mcPx[1]) * (mcPx[0] + mcPx[1]) +
                         (mcPy[0] + mcPy[1]) * (mcPy[0] + mcPy[1]));
    mcEtaPair =
        std::atanh(mcPzPair / sqrt(mcPxPair * mcPxPair + mcPyPair * mcPyPair +
                                   mcPzPair * mcPzPair));
    mcMInv = std::sqrt((mcE[0] + mcE[1]) * (mcE[0] + mcE[1]) -
                       ((mcPx[0] + mcPx[1]) * (mcPx[0] + mcPx[1]) +
                        (mcPy[0] + mcPy[1]) * (mcPy[0] + mcPy[1]) +
                        (mcPz[0] + mcPz[1]) * (mcPz[0] + mcPz[1])));

    mcDeltaPhi = std::min(std::abs(mcPhi[0] - mcPhi[1]),
                          float(2 * M_PI - std::abs(mcPhi[0] - mcPhi[1])));
    mcDeltaPhi *= 180. / M_PI;
    mcDeltaR = std::sqrt((mcEta[0] - mcEta[1]) * (mcEta[0] - mcEta[1]) +
                         (mcPhi[0] - mcPhi[1]) * (mcPhi[0] - mcPhi[1]));
    mcSepAngle =
        acos((mcPx[0] * mcPx[1] + mcPy[0] * mcPy[1] + mcPz[0] * mcPz[1]) /
             (sqrt(mcPx[0] * mcPx[0] + mcPy[0] * mcPy[0] + mcPz[0] * mcPz[0]) *
              sqrt(mcPx[1] * mcPx[1] + mcPy[1] * mcPy[1] + mcPz[1] * mcPz[1])));
    // SepAngle *= 180. / M_PI;

    h_mcPtPair->Fill(mcPtPair);
    h_mcMInv->Fill(mcMInv);
    h_mcDeltaPhi->Fill(mcDeltaPhi);
    h_mcDeltaR->Fill(mcDeltaR);
    h_mcSepAngle->Fill(mcSepAngle);
    h_mcEtaPair->Fill(mcEtaPair);

    if (ntau_rec != 2)
      continue;
    if (matched[0] == 0 || matched[1] == 0)
      continue;

    float ECorr[2] = {1., 1.};
    int n_el[2] = {0, 0};
    int n_mu[2] = {0, 0};
    int n_pipl[2] = {0, 0};
    int n_pimi[2] = {0, 0};
    int n_gamma[2] = {0, 0};
    int n_neutrals[2] = {recNPfos[0] - recNQTracks[0],
                         recNPfos[1] - recNQTracks[1]};
    int eta_range[2] = {eta_classifier(recEta[0]), eta_classifier(recEta[1])};

    int pfos_counter = 0;

    TH1F *histo_corr = nullptr;
    float bin_correction;

    //std::cout << "pfosPdg[19] = " << pfosPdg[19] << '\n';
    //std::cout << "pfosPdg[20] = " << pfosPdg[20] << '\n';

    for (int t = 0; t != ntau_rec; ++t) {
      for (int j = 0; j != recNPfos[t]; j++) {
        if (abs(pfosPdg[j + pfos_counter]) == 11)
          n_el[t]++;
        else if (abs(pfosPdg[j + pfos_counter]) == 13)
          n_mu[t]++;
        else if (pfosPdg[j + pfos_counter] == 22)
          n_gamma[t]++;
        else if (pfosPdg[j + pfos_counter] == 211)
          n_pipl[t]++;
        else if (pfosPdg[j + pfos_counter] == -211)
          n_pimi[t]++;
      }
      pfos_counter += 20;

      if (n_el[t] == 1 && (n_pimi[t] + n_pipl[t] == 0))
        histo_corr = histos_e_array[eta_range[t]];
      else if (n_mu[t] == 1 && (n_pimi[t] + n_pipl[t] == 0))
        histo_corr = histos_mu_array[eta_range[t]];
      else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] == 0))
        histo_corr = histos_1prong0NP_array[eta_range[t]];
      else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] == 1))
        histo_corr = histos_1prong1NP_array[eta_range[t]];
      else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] == 2))
        histo_corr = histos_1prong2NP_array[eta_range[t]];
      else if ((n_pimi[t] + n_pipl[t] == 1) && (n_neutrals[t] >= 3))
        histo_corr = histos_1prong3NP_array[eta_range[t]];
      else if (n_pimi[t] + n_pipl[t] == 3)
        histo_corr = histos_3prong_array[eta_range[t]];
      else
        histo_corr = histos_others_array[eta_range[t]];

      bin_correction = histo_corr->GetBinContent(histo_corr->FindBin(recE[t]));
      bin_correction > 0.01 ? ECorr[t] = bin_correction : ECorr[t] = 1.;

      histo_corr = nullptr;
    }

    //ECorr[0] = 1.;
    //ECorr[1] = 1.;

    float recPxCorr[2] = {recPx[0] * ECorr[0], recPx[1] * ECorr[1]};
    float recPyCorr[2] = {recPy[0] * ECorr[0], recPy[1] * ECorr[1]};
    float recPzCorr[2] = {recPz[0] * ECorr[0], recPz[1] * ECorr[1]};
    float recPtCorr[2] = {recPt[0] * ECorr[0], recPt[1] * ECorr[1]};
    float recECorr[2] = {recE[0] * ECorr[0], recE[1] * ECorr[1]};

    h_recPtSingle->Fill(recPtCorr[0]);
    h_recPtSingle->Fill(recPtCorr[1]);

    recPxPair = recPxCorr[0] + recPxCorr[1];
    recPyPair = recPyCorr[0] + recPyCorr[1];
    recPzPair = recPzCorr[0] + recPzCorr[1];

    recPtPair = std::sqrt(
        (recPxCorr[0] + recPxCorr[1]) * (recPxCorr[0] + recPxCorr[1]) +
        (recPyCorr[0] + recPyCorr[1]) * (recPyCorr[0] + recPyCorr[1]));
    recEtaPair = std::atanh(recPzPair /
                            sqrt(recPxPair * recPxPair + recPyPair * recPyPair +
                                 recPzPair * recPzPair));
    recMInv = std::sqrt(
        (recECorr[0] + recECorr[1]) * (recECorr[0] + recECorr[1]) -
        ((recPxCorr[0] + recPxCorr[1]) * (recPxCorr[0] + recPxCorr[1]) +
         (recPyCorr[0] + recPyCorr[1]) * (recPyCorr[0] + recPyCorr[1]) +
         (recPzCorr[0] + recPzCorr[1]) * (recPzCorr[0] + recPzCorr[1])));
    recDeltaPhi = std::min(std::abs(recPhi[0] - recPhi[1]),
                           float(2 * M_PI - std::abs(recPhi[0] - recPhi[1])));
    recDeltaPhi *= 180. / M_PI;
    recDeltaR = std::sqrt((recEta[0] - recEta[1]) * (recEta[0] - recEta[1]) +
                          (recPhi[0] - recPhi[1]) * (recPhi[0] - recPhi[1]));
    recSepAngle =
        acos((recPxCorr[0] * recPxCorr[1] + recPyCorr[0] * recPyCorr[1] +
              recPzCorr[0] * recPzCorr[1]) /
             (sqrt(recPxCorr[0] * recPxCorr[0] + recPyCorr[0] * recPyCorr[0] +
                   recPzCorr[0] * recPzCorr[0]) *
              sqrt(recPxCorr[1] * recPxCorr[1] + recPyCorr[1] * recPyCorr[1] +
                   recPzCorr[1] * recPzCorr[1])));

    h_recPtPair->Fill(recPtPair);
    h_recMInv->Fill(recMInv);
    h_recDeltaPhi->Fill(recDeltaPhi);
    h_recDeltaR->Fill(recDeltaR);
    h_recSepAngle->Fill(recSepAngle);
    h_recEtaPair->Fill(recEtaPair);
  }

  delete[] mcE;
  delete[] mcPx;
  delete[] mcPy;
  delete[] mcPz;
  delete[] mcPt;
  delete[] mcEta;
  delete[] mcPhi;
  delete[] mcE_vis;
  delete[] mcPx_vis;
  delete[] mcPy_vis;
  delete[] mcPz_vis;
  delete[] mcPt_vis;
  delete[] mcEta_vis;
  delete[] mcPhi_vis;
  delete[] recE;
  delete[] recPx;
  delete[] recPy;
  delete[] recPz;
  delete[] recPt;
  delete[] recEta;
  delete[] recPhi;
  delete[] recNPfos;
  delete[] recNQTracks;
  delete[] pfosPdg;
  delete[] matched;

  tree->Delete();
  input_file->Close();
  input_file->Delete();
}

void delete_histograms(std::vector<TH1F *> &histos_vec) {
  for (unsigned int i = 0; i != histos_vec.size(); ++i) {
    histos_vec[i]->Delete();
  }
}

void Plots_drawer(std::vector<TH1F *> Hhistos_vec,
                  std::vector<TH1F *> Zhistos_vec) {

  auto mystyle = new TStyle("mystyle", "mystyle Style");

  mystyle->SetLegendBorderSize(0);
  mystyle->SetOptStat(111111);
  mystyle->SetLegendFillColor(-1);
  mystyle->SetLegendTextSize(0.03);
  mystyle->SetTitleBorderSize(0);
  mystyle->SetHistFillStyle(3002);

  auto c_recMInv = new TCanvas("c_recMInv", "c_recMInv", 800, 800);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  Hhistos_vec[0]->Scale(1 / Hhistos_vec[0]->Integral());
  Zhistos_vec[0]->Scale(1 / Zhistos_vec[0]->Integral());
  Hhistos_vec[0]->SetLineColor(kRed);
  Hhistos_vec[0]->SetFillColor(kRed);
  Hhistos_vec[0]->SetFillStyle(3002);
  Zhistos_vec[0]->SetLineColor(kBlue);
  Zhistos_vec[0]->SetFillColor(kBlue);
  Zhistos_vec[0]->SetFillStyle(3003);
  Hhistos_vec[0]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[0]->GetMaximum() * 1.2);
  Hhistos_vec[0]->Draw("HISTO");
  Zhistos_vec[0]->Draw("HISTO, SAME");
  auto leg0 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg0->AddEntry(Hhistos_vec[0], "H#rightarrow#tau#tau", "f");
  leg0->AddEntry(Zhistos_vec[0], "Z#rightarrow#tau#tau", "f");
  leg0->Draw();
  gPad->Modified();
  c_recMInv->SaveAs("./plots/c_recMInv.png");
  c_recMInv->Close();

  auto c_mcMInv = new TCanvas("c_mcMInv", "c_mcMInv", 800, 800);
  Hhistos_vec[1]->Scale(1 / Hhistos_vec[1]->Integral());
  Zhistos_vec[1]->Scale(1 / Zhistos_vec[1]->Integral());
  Hhistos_vec[1]->SetLineColor(kRed);
  Hhistos_vec[1]->SetFillColor(kRed);
  Hhistos_vec[1]->SetFillStyle(3002);
  Zhistos_vec[1]->SetLineColor(kBlue);
  Zhistos_vec[1]->SetFillColor(kBlue);
  Zhistos_vec[1]->SetFillStyle(3003);
  Hhistos_vec[1]->GetXaxis()->SetRangeUser(70., 130.);
  Hhistos_vec[1]->Draw("HISTO");
  Zhistos_vec[1]->Draw("HISTO, SAME");
  auto leg1 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg1->AddEntry(Hhistos_vec[1], "H#rightarrow#tau#tau", "f");
  leg1->AddEntry(Zhistos_vec[1], "Z#rightarrow#tau#tau", "f");
  leg1->Draw();
  gPad->Modified();
  c_mcMInv->SaveAs("./plots/c_mcMInv.png");
  c_mcMInv->Close();

  auto c_recPtPair = new TCanvas("c_recPtPair", "c_recPtPair", 800, 800);
  Hhistos_vec[2]->Scale(1 / Hhistos_vec[2]->Integral());
  Zhistos_vec[2]->Scale(1 / Zhistos_vec[2]->Integral());
  Hhistos_vec[2]->SetLineColor(kRed);
  Hhistos_vec[2]->SetFillColor(kRed);
  Hhistos_vec[2]->SetFillStyle(3002);
  Zhistos_vec[2]->SetLineColor(kBlue);
  Zhistos_vec[2]->SetFillColor(kBlue);
  Zhistos_vec[2]->SetFillStyle(3003);
  Hhistos_vec[2]->GetXaxis()->SetRangeUser(0., 250.);
  Hhistos_vec[2]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[2]->GetMaximum() * 1.2);
  Hhistos_vec[2]->Draw("HISTO");
  Zhistos_vec[2]->Draw("HISTO, SAME");
  auto leg2 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg2->AddEntry(Hhistos_vec[2], "H#rightarrow#tau#tau", "f");
  leg2->AddEntry(Zhistos_vec[2], "Z#rightarrow#tau#tau", "f");
  leg2->Draw();
  gPad->Modified();
  c_recPtPair->SaveAs("./plots/c_recPtPair.png");
  c_recPtPair->Close();

  auto c_mcPtPair = new TCanvas("c_mcPtPair", "c_mcPtPair", 800, 800);
  Hhistos_vec[3]->Scale(1 / Hhistos_vec[3]->Integral());
  Zhistos_vec[3]->Scale(1 / Zhistos_vec[3]->Integral());
  Hhistos_vec[3]->SetLineColor(kRed);
  Hhistos_vec[3]->SetFillColor(kRed);
  Hhistos_vec[3]->SetFillStyle(3002);
  Zhistos_vec[3]->SetLineColor(kBlue);
  Zhistos_vec[3]->SetFillColor(kBlue);
  Zhistos_vec[3]->SetFillStyle(3003);
  Hhistos_vec[3]->GetXaxis()->SetRangeUser(0., 250.);
  Hhistos_vec[3]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[3]->GetMaximum() * 1.2);
  Hhistos_vec[3]->Draw("HISTO");
  Zhistos_vec[3]->Draw("HISTO, SAME");
  auto leg3 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg3->AddEntry(Hhistos_vec[3], "H#rightarrow#tau#tau", "f");
  leg3->AddEntry(Zhistos_vec[3], "Z#rightarrow#tau#tau", "f");
  leg3->Draw();
  gPad->Modified();
  c_mcPtPair->SaveAs("./plots/c_mcPtPair.png");
  c_mcPtPair->Close();

  auto c_recEtaPair = new TCanvas("c_recEtaPair", "c_recEtaPair", 800, 800);
  Hhistos_vec[4]->Scale(1 / Hhistos_vec[4]->Integral());
  Zhistos_vec[4]->Scale(1 / Zhistos_vec[4]->Integral());
  Hhistos_vec[4]->SetLineColor(kRed);
  Hhistos_vec[4]->SetFillColor(kRed);
  Hhistos_vec[4]->SetFillStyle(3002);
  Zhistos_vec[4]->SetLineColor(kBlue);
  Zhistos_vec[4]->SetFillColor(kBlue);
  Zhistos_vec[4]->SetFillStyle(3003);
  Hhistos_vec[4]->SetMinimum(0.);
  Hhistos_vec[4]->GetYaxis()->SetRangeUser(0.,
                                           Hhistos_vec[4]->GetMaximum() * 1.2);
  Hhistos_vec[4]->Draw("HISTO");
  Zhistos_vec[4]->Draw("HISTO, SAME");
  auto leg4 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg4->AddEntry(Hhistos_vec[4], "H#rightarrow#tau#tau", "f");
  leg4->AddEntry(Zhistos_vec[4], "Z#rightarrow#tau#tau", "f");
  leg4->Draw();
  gPad->Modified();
  c_recEtaPair->SaveAs("./plots/c_recEtaPair.png");
  c_recEtaPair->Close();

  auto c_mcEtaPair = new TCanvas("c_mcEtaPair", "c_mcEtaPair", 800, 800);
  Hhistos_vec[5]->Scale(1 / Hhistos_vec[5]->Integral());
  Zhistos_vec[5]->Scale(1 / Zhistos_vec[5]->Integral());
  Hhistos_vec[5]->SetLineColor(kRed);
  Hhistos_vec[5]->SetFillColor(kRed);
  Hhistos_vec[5]->SetFillStyle(3002);
  Zhistos_vec[5]->SetLineColor(kBlue);
  Zhistos_vec[5]->SetFillColor(kBlue);
  Zhistos_vec[5]->SetFillStyle(3003);
  Hhistos_vec[5]->SetMinimum(0.);
  Hhistos_vec[5]->GetYaxis()->SetRangeUser(0.,
                                           Hhistos_vec[5]->GetMaximum() * 1.2);
  Hhistos_vec[5]->Draw("HISTO");
  Zhistos_vec[5]->Draw("HISTO, SAME");
  auto leg5 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg5->AddEntry(Hhistos_vec[5], "H#rightarrow#tau#tau", "f");
  leg5->AddEntry(Zhistos_vec[5], "Z#rightarrow#tau#tau", "f");
  leg5->Draw();
  gPad->Modified();
  c_mcEtaPair->SaveAs("./plots/c_mcEtaPair.png");
  c_mcEtaPair->Close();

  auto c_recPtSingle = new TCanvas("c_recPtSingle", "c_recPtSingle", 800, 800);
  Hhistos_vec[6]->Scale(1 / Hhistos_vec[6]->Integral());
  Zhistos_vec[6]->Scale(1 / Zhistos_vec[6]->Integral());
  Hhistos_vec[6]->SetLineColor(kRed);
  Hhistos_vec[6]->SetFillColor(kRed);
  Hhistos_vec[6]->SetFillStyle(3002);
  Zhistos_vec[6]->SetLineColor(kBlue);
  Zhistos_vec[6]->SetFillColor(kBlue);
  Zhistos_vec[6]->SetFillStyle(3003);
  Hhistos_vec[6]->GetXaxis()->SetRangeUser(0., 250.);
  Hhistos_vec[6]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[6]->GetMaximum() * 1.2);
  Hhistos_vec[6]->Draw("HISTO");
  Zhistos_vec[6]->Draw("HISTO, SAME");
  auto leg6 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg6->AddEntry(Hhistos_vec[6], "H#rightarrow#tau#tau", "f");
  leg6->AddEntry(Zhistos_vec[6], "Z#rightarrow#tau#tau", "f");
  leg6->Draw();
  gPad->Modified();
  c_recPtSingle->SaveAs("./plots/c_recPtSingle.png");
  c_recPtSingle->Close();

  auto c_mcPtSingle = new TCanvas("c_mcPtSingle", "c_mcPtSingle", 800, 800);
  Hhistos_vec[7]->Scale(1 / Hhistos_vec[7]->Integral());
  Zhistos_vec[7]->Scale(1 / Zhistos_vec[7]->Integral());
  Hhistos_vec[7]->SetLineColor(kRed);
  Hhistos_vec[7]->SetFillColor(kRed);
  Hhistos_vec[7]->SetFillStyle(3002);
  Zhistos_vec[7]->SetLineColor(kBlue);
  Zhistos_vec[7]->SetFillColor(kBlue);
  Zhistos_vec[7]->SetFillStyle(3003);
  Hhistos_vec[7]->GetXaxis()->SetRangeUser(0., 250.);
  Hhistos_vec[7]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[7]->GetMaximum() * 1.2);
  Hhistos_vec[7]->Draw("HISTO");
  Zhistos_vec[7]->Draw("HISTO, SAME");
  auto leg7 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg7->AddEntry(Hhistos_vec[7], "H#rightarrow#tau#tau", "f");
  leg7->AddEntry(Zhistos_vec[7], "Z#rightarrow#tau#tau", "f");
  leg7->Draw();
  gPad->Modified();
  c_mcPtSingle->SaveAs("./plots/c_mcPtSingle.png");
  c_mcPtSingle->Close();

  auto c_recDeltaPhi = new TCanvas("c_recDeltaPhi", "c_recDeltaPhi", 800, 800);
  Hhistos_vec[8]->Scale(1 / Hhistos_vec[8]->Integral());
  Zhistos_vec[8]->Scale(1 / Zhistos_vec[8]->Integral());
  Hhistos_vec[8]->SetLineColor(kRed);
  Hhistos_vec[8]->SetFillColor(kRed);
  Hhistos_vec[8]->SetFillStyle(3002);
  Zhistos_vec[8]->SetLineColor(kBlue);
  Zhistos_vec[8]->SetFillColor(kBlue);
  Zhistos_vec[8]->SetFillStyle(3003);
  Hhistos_vec[8]->SetMinimum(0.);
  Hhistos_vec[8]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[8]->GetMaximum() * 1.2);
  Hhistos_vec[8]->Draw("HISTO");
  Zhistos_vec[8]->Draw("HISTO, SAME");
  auto leg8 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg8->AddEntry(Hhistos_vec[8], "H#rightarrow#tau#tau", "f");
  leg8->AddEntry(Zhistos_vec[8], "Z#rightarrow#tau#tau", "f");
  leg8->Draw();
  gPad->Modified();
  c_recDeltaPhi->Update();
  c_recDeltaPhi->SaveAs("./plots/c_recDeltaPhi.png");
  c_recDeltaPhi->Close();

  auto c_mcDeltaPhi = new TCanvas("c_mcDeltaPhi", "c_mcDeltaPhi", 800, 800);
  Hhistos_vec[9]->Scale(1 / Hhistos_vec[9]->Integral());
  Zhistos_vec[9]->Scale(1 / Zhistos_vec[9]->Integral());
  Hhistos_vec[9]->SetLineColor(kRed);
  Hhistos_vec[9]->SetFillColor(kRed);
  Hhistos_vec[9]->SetFillStyle(3002);
  Zhistos_vec[9]->SetLineColor(kBlue);
  Zhistos_vec[9]->SetFillColor(kBlue);
  Zhistos_vec[9]->SetFillStyle(3003);
  Hhistos_vec[9]->SetMinimum(0.);
  Hhistos_vec[9]->GetYaxis()->SetRangeUser(0.,
                                           Zhistos_vec[9]->GetMaximum() * 1.2);
  Hhistos_vec[9]->Draw("HISTO");
  Zhistos_vec[9]->Draw("HISTO, SAME");
  auto leg9 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg9->AddEntry(Hhistos_vec[9], "H#rightarrow#tau#tau", "f");
  leg9->AddEntry(Zhistos_vec[9], "Z#rightarrow#tau#tau", "f");
  leg9->Draw();
  gPad->Modified();
  c_mcDeltaPhi->SaveAs("./plots/c_mcDeltaPhi.png");
  c_mcDeltaPhi->Close();

  auto c_recDeltaR = new TCanvas("c_recDeltaR", "c_recDeltaR", 800, 800);
  Hhistos_vec[10]->Scale(1 / Hhistos_vec[10]->Integral());
  Zhistos_vec[10]->Scale(1 / Zhistos_vec[10]->Integral());
  Hhistos_vec[10]->SetLineColor(kRed);
  Hhistos_vec[10]->SetFillColor(kRed);
  Hhistos_vec[10]->SetFillStyle(3002);
  Zhistos_vec[10]->SetLineColor(kBlue);
  Zhistos_vec[10]->SetFillColor(kBlue);
  Zhistos_vec[10]->SetFillStyle(3003);
  Hhistos_vec[10]->GetYaxis()->SetRangeUser(0., Zhistos_vec[10]->GetMaximum() *
                                                    1.2);
  Hhistos_vec[10]->Draw("HISTO");
  Zhistos_vec[10]->Draw("HISTO, SAME");
  auto leg10 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg10->AddEntry(Hhistos_vec[10], "H#rightarrow#tau#tau", "f");
  leg10->AddEntry(Zhistos_vec[10], "Z#rightarrow#tau#tau", "f");
  leg10->Draw();
  gPad->Modified();
  c_recDeltaR->Update();
  c_recDeltaR->SaveAs("./plots/c_recDeltaR.png");
  c_recDeltaR->Close();

  auto c_mcDeltaR = new TCanvas("c_mcDeltaR", "c_mcDeltaR", 800, 800);
  Hhistos_vec[11]->Scale(1 / Hhistos_vec[11]->Integral());
  Zhistos_vec[11]->Scale(1 / Zhistos_vec[11]->Integral());
  Hhistos_vec[11]->SetLineColor(kRed);
  Hhistos_vec[11]->SetFillColor(kRed);
  Hhistos_vec[11]->SetFillStyle(3002);
  Zhistos_vec[11]->SetLineColor(kBlue);
  Zhistos_vec[11]->SetFillColor(kBlue);
  Zhistos_vec[11]->SetFillStyle(3003);
  Hhistos_vec[11]->GetYaxis()->SetRangeUser(0., Zhistos_vec[11]->GetMaximum() *
                                                    1.2);
  Hhistos_vec[11]->Draw("HISTO");
  Zhistos_vec[11]->Draw("HISTO, SAME");
  auto leg11 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg11->AddEntry(Hhistos_vec[11], "H#rightarrow#tau#tau", "f");
  leg11->AddEntry(Zhistos_vec[11], "Z#rightarrow#tau#tau", "f");
  leg11->Draw();
  gPad->Modified();
  c_mcDeltaR->SaveAs("./plots/c_mcDeltaR.png");
  c_mcDeltaR->Close();

  auto c_recSepAngle = new TCanvas("c_recSepAngle", "c_recSepAngle", 800, 800);
  Hhistos_vec[12]->Scale(1 / Hhistos_vec[12]->Integral());
  Zhistos_vec[12]->Scale(1 / Zhistos_vec[12]->Integral());
  Hhistos_vec[12]->SetLineColor(kRed);
  Hhistos_vec[12]->SetFillColor(kRed);
  Hhistos_vec[12]->SetFillStyle(3002);
  Zhistos_vec[12]->SetLineColor(kBlue);
  Zhistos_vec[12]->SetFillColor(kBlue);
  Zhistos_vec[12]->SetFillStyle(3003);
  Hhistos_vec[12]->GetYaxis()->SetRangeUser(0., Zhistos_vec[12]->GetMaximum() *
                                                    1.2);
  Hhistos_vec[12]->Draw("HISTO");
  Zhistos_vec[12]->Draw("HISTO, SAME");
  auto leg12 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg12->AddEntry(Hhistos_vec[12], "H#rightarrow#tau#tau", "f");
  leg12->AddEntry(Zhistos_vec[12], "Z#rightarrow#tau#tau", "f");
  leg12->Draw();
  gPad->Modified();
  c_recSepAngle->SaveAs("./plots/c_recSepAngle.png");
  c_recSepAngle->Close();

  auto c_mcSepAngle = new TCanvas("c_mcSepAngle", "c_mcSepAngle", 800, 800);
  Hhistos_vec[13]->Scale(1 / Hhistos_vec[13]->Integral());
  Zhistos_vec[13]->Scale(1 / Zhistos_vec[13]->Integral());
  Hhistos_vec[13]->SetLineColor(kRed);
  Hhistos_vec[13]->SetFillColor(kRed);
  Hhistos_vec[13]->SetFillStyle(3002);
  Zhistos_vec[13]->SetLineColor(kBlue);
  Zhistos_vec[13]->SetFillColor(kBlue);
  Zhistos_vec[13]->SetFillStyle(3003);
  Hhistos_vec[13]->GetYaxis()->SetRangeUser(0., Zhistos_vec[13]->GetMaximum() *
                                                    1.2);

  Hhistos_vec[13]->Draw("HISTO");
  Zhistos_vec[13]->Draw("HISTO, SAME");
  auto leg13 = new TLegend(0.73, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg13->AddEntry(Hhistos_vec[13], "H#rightarrow#tau#tau", "f");
  leg13->AddEntry(Zhistos_vec[13], "Z#rightarrow#tau#tau", "f");
  leg13->Draw();
  gPad->Modified();
  c_mcSepAngle->SaveAs("./plots/c_mcSepAngle.png");
  c_mcSepAngle->Close();
}

void HvsZ(const TString filenameH = "/eos/user/l/lvalla/MuColl/H2tautau10k/"
                                    "Outputs/output_EvalTauFinder.root",
          const TString filenameZ = "/eos/user/l/lvalla/MuColl/Z2tautau10k/"
                                    "Outputs/output_EvalTauFinder.root") {

  auto H_histograms = H_histo_maker();
  auto Z_histograms = Z_histo_maker();

  Histo_filler(H_histograms, filenameH);
  Histo_filler(Z_histograms, filenameZ);

  Plots_drawer(H_histograms, Z_histograms);

  delete_histograms(H_histograms);
  delete_histograms(Z_histograms);
}