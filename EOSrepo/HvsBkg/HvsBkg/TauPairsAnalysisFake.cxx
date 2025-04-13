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

TH1F *h_recMInv, *h_mcMInv, *h_recPtPair, *h_mcPtPair, *h_recPtSingle,
    *h_mcPtSingle, *h_recDeltaPhi, *h_mcDeltaPhi, *h_mcDeltaEta, *h_recDeltaEta,
    *h_recDeltaR, *h_mcDeltaR, *h_recSepAngle, *h_mcSepAngle, *h_recEtaPair,
    *h_mcEtaPair, *h_recD0Single, *h_recSigmaD0Single, *h_recChAsymm,
    *h_recD0overSigma;

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

void HistoMCFiller(const TString &filename, const TString &rootfilename,
                   const TString &boson) {
  // --- Histograms initialization
  {
    h_mcMInv = new TH1F(
        "h_mcMInv",
        boson + "#rightarrow#tau#tau MC Invariant mass ; M_{inv} [GeV]", 40,
        80., 100.);
    h_mcPtPair = new TH1F(
        "h_mcPtPair",
        boson + "#rightarrow#tau#tau MC p_{T} of the system; p_{T} [GeV]", 50,
        0., 300.);
    h_mcPtSingle = new TH1F(
        "h_mcPtSingle",
        boson + "#rightarrow#tau#tau MC p_{T} of single taus; p_{T} [GeV]", 50,
        0., 500.);
    h_mcDeltaPhi = new TH1F(
        "h_mcDeltaPhi",
        boson +
            "#rightarrow#tau#tau MC #Delta#Phi between the taus; |#Delta#Phi|",
        40, 0., 180.);
    h_mcDeltaEta = new TH1F(
        "h_mcDeltaEta",
        boson +
            "#rightarrow#tau#tau MC #Delta#eta between the taus; |#Delta#eta|",
        40, 0., 4.);
    h_mcDeltaR = new TH1F(
        "h_mcDeltaR",
        boson + "#rightarrow#tau#tau MC #DeltaR between the taus; #DeltaR", 50,
        0., 4.5);
    h_mcSepAngle = new TH1F(
        "h_mcSepAngle",
        boson + "#rightarrow#tau#tau MC angle between the taus; #alpha [rad]",
        50, 0., 3.1416);
    h_mcEtaPair =
        new TH1F("h_mcEtaPair",
                 boson + "#rightarrow#tau#tau MC #eta of the system; #eta", 50,
                 -2.5, 2.5);
  }

  // Open evaltaufinder output file
  TFile *input_file = new TFile(filename, "READ");

  if (!input_file->IsOpen())
    throw std::invalid_argument("filename not valid");

  TTree *inTree = (TTree *)input_file->Get("evtree");

  int evID, runID;
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

  inTree->SetBranchAddress("EvID", &evID);
  inTree->SetBranchAddress("RunID", &runID);
  inTree->SetBranchAddress("mcE", mcE);
  inTree->SetBranchAddress("mcPx", mcPx);
  inTree->SetBranchAddress("mcPy", mcPy);
  inTree->SetBranchAddress("mcPz", mcPz);
  inTree->SetBranchAddress("mcPt", mcPt);
  inTree->SetBranchAddress("mcEta", mcEta);
  inTree->SetBranchAddress("mcPhi", mcPhi);
  inTree->SetBranchAddress("mcDecayMode", mcDecayMode);
  inTree->SetBranchAddress("mcE_vis", mcE_vis);
  inTree->SetBranchAddress("mcPx_vis", mcPx_vis);
  inTree->SetBranchAddress("mcPy_vis", mcPy_vis);
  inTree->SetBranchAddress("mcPz_vis", mcPz_vis);
  inTree->SetBranchAddress("mcPt_vis", mcPz_vis);
  inTree->SetBranchAddress("mcEta_vis", mcEta_vis);
  inTree->SetBranchAddress("mcPhi_vis", mcPhi_vis);

  const long int nEntries = inTree->GetEntries();

  float mcPxPair = 0., mcPyPair = 0., mcPzPair = 0., mcPtPair = 0., mcMInv = 0.,
        mcDeltaPhi = 0., mcDeltaR = 0., mcDeltaEta = 0., mcSepAngle = 0.,
        mcEtaPair = 0.;

  for (int ientry = 0; ientry < nEntries; ++ientry) {

    inTree->GetEntry(ientry);

    /*if (mcDecayMode[0] == 8 || mcDecayMode[1] == 8) {
      std::cout << "Run: " << runID << " - "
                << "Ev: " << evID << '\n';
    }*/

    h_mcPtSingle->Fill(mcPt[0]);
    h_mcPtSingle->Fill(mcPt[1]);

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

    mcDeltaPhi = std::min(std::fabs(mcPhi[0] - mcPhi[1]),
                          float(2 * M_PI - std::fabs(mcPhi[0] - mcPhi[1])));

    mcDeltaEta = std::fabs(mcEta[0] - mcEta[1]);
    mcDeltaR = std::sqrt(mcDeltaEta * mcDeltaEta + mcDeltaPhi * mcDeltaPhi);
    mcSepAngle =
        acos((mcPx[0] * mcPx[1] + mcPy[0] * mcPy[1] + mcPz[0] * mcPz[1]) /
             (sqrt(mcPx[0] * mcPx[0] + mcPy[0] * mcPy[0] + mcPz[0] * mcPz[0]) *
              sqrt(mcPx[1] * mcPx[1] + mcPy[1] * mcPy[1] + mcPz[1] * mcPz[1])));
    // SepAngle *= 180. / M_PI;

    h_mcPtPair->Fill(mcPtPair);
    h_mcMInv->Fill(mcMInv);
    h_mcDeltaPhi->Fill(mcDeltaPhi * 180. / M_PI);
    h_mcDeltaEta->Fill(mcDeltaEta);
    h_mcDeltaR->Fill(mcDeltaR);
    h_mcSepAngle->Fill(mcSepAngle);
    h_mcEtaPair->Fill(mcEtaPair);
  }

  inTree->Delete();
  input_file->Close();
  input_file->Delete();

  gStyle->SetOptStat(111111);
  h_mcPtSingle->Scale(1. / h_mcPtSingle->Integral());
  h_mcPtPair->Scale(1. / h_mcPtPair->Integral());
  h_mcMInv->Scale(1. / h_mcMInv->Integral());
  h_mcDeltaPhi->Scale(1. / h_mcDeltaPhi->Integral());
  h_mcDeltaR->Scale(1. / h_mcDeltaR->Integral());
  h_mcSepAngle->Scale(1. / h_mcSepAngle->Integral());
  h_mcEtaPair->Scale(1. / h_mcEtaPair->Integral());

  TFile *rootfile = new TFile(rootfilename, "UPDATE");
  TList *histo_list = new TList();
  histo_list->Add(h_mcPtSingle);
  histo_list->Add(h_mcPtPair);
  histo_list->Add(h_mcMInv);
  histo_list->Add(h_mcDeltaPhi);
  histo_list->Add(h_mcDeltaR);
  histo_list->Add(h_mcSepAngle);
  histo_list->Add(h_mcEtaPair);
  histo_list->Write(boson + "_MC_histograms", TObject::kSingleKey);
  rootfile->Close();

  {
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
    h_mcPtPair->Delete();
    h_mcPtSingle->Delete();
    h_mcMInv->Delete();
    h_mcDeltaPhi->Delete();
    h_mcDeltaEta->Delete();
    h_mcDeltaR->Delete();
    h_mcSepAngle->Delete();
    h_mcEtaPair->Delete();
  }
}

void HistoRecoFiller(const TString &filename, const TString &rootfilename,
                     const TString &recdecay0 = "all",
                     const TString &recdecay1 = "all",
                     const TString &boson = "B",
                     const bool corrections = true) {
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
  {
    h_recMInv = new TH1F("h_recMInv",
                         boson + "#rightarrow#tau#tau Reco Invariant mass of "
                                 "the system; m_{vis} (#tau#tau) [GeV/c^{2}]",
                         40, 0., 280.);
    h_recPtPair = new TH1F("h_recPtPair",
                           boson + "#rightarrow#tau#tau Reco p_{T} of the "
                                   "system; p_{T} (#tau#tau) [GeV/c]",
                           25, 0., 300.);
    h_recPtSingle = new TH1F(
        "h_recPtSingle",
        boson + "#rightarrow#tau#tau Reco p_{T} of single taus; p_{T} [GeV/c]",
        30, 0., 300.);
    h_recDeltaPhi =
        new TH1F("h_recDeltaPhi",
                 boson + "#rightarrow#tau#tau Reco #Delta#Phi "
                         "between the taus; |#Delta#Phi (#tau#tau)|",
                 30, 0., 180.);
    h_recDeltaEta = new TH1F("h_recDeltaEta",
                             boson + "#rightarrow#tau#tau Reco #Delta#eta "
                                     "between the taus; |#Delta#eta (#tau#tau)|",
                             40, 0., 4.);
    h_recDeltaR = new TH1F("h_recDeltaR",
                           boson + "#rightarrow#tau#tau Reco #DeltaR between "
                                   "the taus; #DeltaR (#tau#tau)",
                           30, 0., 4.5);
    h_recSepAngle = new TH1F("h_recSepAngle",
                             boson + "#rightarrow#tau#tau Reco angle between "
                                     "the taus; #alpha_{#tau#tau} [rad]",
                             50, 0., 3.1416);
    h_recEtaPair = new TH1F(
        "h_recEtaPair",
        boson + "#rightarrow#tau#tau MC #eta of the system; #eta_{#tau#tau}",
        50, -2.5, 2.5);
    h_recChAsymm = new TH1F("h_recChAsymm",
                            boson + "#rightarrow#tau#tau Reco charge asymmetry "
                                    "of the leading track; #Upsilon",
                            25, -1., 1.);
    h_recD0Single = new TH1F(
        "h_recD0Single",
        boson + "#rightarrow#tau#tau Reco D0 of the single taus; d_{0} [mm]",
        nlogbins_D0, logbinedges_D0);
    h_recSigmaD0Single =
        new TH1F("h_recSigmaD0Single",
                 boson + "#rightarrow#tau#tau Reco error on D0 of the single "
                         "taus; #sigma_{d_{0}} [mm]",
                 nlogbins_SigmaD0, logbinedges_SigmaD0);
    h_recD0overSigma =
        new TH1F("h_recD0overSigma",
                 boson + "#rightarrow#tau#tau Reco Significance of the "
                         "taus; d_{0}/#sigma_{d_{0}}",
                 nlogbins_D0overSigma, logbinedges_D0overSigma);
  }

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
  /*
  TH1F *histos_1pi2N_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1pi2N->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1pi2N->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1pi2N->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1pi2N->FindObject("h_ECorr_Eta3"))};
  TH1F *histos_1pi3N_array[4] = {
      static_cast<TH1F *>(ECorrHistos_1pi3N->FindObject("h_ECorr_Eta0")),
      static_cast<TH1F *>(ECorrHistos_1pi3N->FindObject("h_ECorr_Eta1")),
      static_cast<TH1F *>(ECorrHistos_1pi3N->FindObject("h_ECorr_Eta2")),
      static_cast<TH1F *>(ECorrHistos_1pi3N->FindObject("h_ECorr_Eta3"))};
  */
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
  float *recE = new float[100000];
  float *recPx = new float[100000];
  float *recPy = new float[100000];
  float *recPz = new float[100000];
  float *recPt = new float[100000];
  float *recEta = new float[100000];
  float *recPhi = new float[100000];
  float *recD0 = new float[100000];
  float *recSigmaD0 = new float[100000];
  int *recNPfos = new int[500000];
  int *recNQTracks = new int[500000];
  float *pfosPx = new float[500000];
  float *pfosPy = new float[500000];
  float *pfosPz = new float[500000];
  float *pfosPt = new float[500000];
  int *pfosPdg = new int[500000];
  int *matched = new int[100000];

  inTree->SetBranchAddress("nTausRec", &ntau_rec);
  inTree->SetBranchAddress("EvID", &evID);
  inTree->SetBranchAddress("RunID", &runID);
  inTree->SetBranchAddress("mcDecayMode", mcDecayMode);
  inTree->SetBranchAddress("charge", charge);
  inTree->SetBranchAddress("recE", recE);
  inTree->SetBranchAddress("recPx", recPx);
  inTree->SetBranchAddress("recPy", recPy);
  inTree->SetBranchAddress("recPz", recPz);
  inTree->SetBranchAddress("recPt", recPt);
  inTree->SetBranchAddress("recEta", recEta);
  inTree->SetBranchAddress("recPhi", recPhi);
  inTree->SetBranchAddress("recNPfos", recNPfos);
  inTree->SetBranchAddress("recNQTracks", recNQTracks);
  inTree->SetBranchAddress("recD0", recD0);
  inTree->SetBranchAddress("recSigmaD0", recSigmaD0);
  inTree->SetBranchAddress("pfosPx", pfosPx);
  inTree->SetBranchAddress("pfosPy", pfosPy);
  inTree->SetBranchAddress("pfosPz", pfosPz);
  inTree->SetBranchAddress("pfosPt", pfosPt);
  inTree->SetBranchAddress("pfosPdg", pfosPdg);
  inTree->SetBranchAddress("Matched", matched);

  float recPxPair = 0., recPyPair = 0., recPzPair = 0., recPtPair = 0.,
        recMInv = 0., recDeltaPhi = 0., recDeltaEta = 0., recDeltaR = 0.,
        recSepAngle = 0., recEtaPair = 0., recChAsymm[2] = {0., 0.},
        recPtSingle[2] = {0., 0.}, recEtaSingle[2] = {0., 0.},
        recD0overSigma[2] = {0., 0.};
  int recDM[2] = {-1, -1};

  float recChAsymm_0 = 0., recChAsymm_1 = 0., recPtSingle_0 = 0.,
        recPtSingle_1 = 0., recD0overSigma_0 = 0., recD0overSigma_1 = 0.,
        recEtaSingle_0 = 0., recEtaSingle_1 = 0., recE_0 = 0., recE_1 = 0.;

  int recDM_0 = -1, recDM_1 = -1, recNPfos_0 = 0, recNPfos_1 = 0;

  TFile *rootfile = new TFile(rootfilename, "UPDATE");
  TTree *outTree = new TTree("outTree", "output Tree");
  outTree->Branch("recPtPair", &recPtPair, "recPtPair/F");
  outTree->Branch("recEtaPair", &recEtaPair, "recEtaPair/F");
  outTree->Branch("recMInv", &recMInv, "recMInv/F");
  outTree->Branch("recDeltaPhi", &recDeltaPhi, "recDeltaPhi/F");
  outTree->Branch("recDeltaEta", &recDeltaEta, "recDeltaEta/F");
  outTree->Branch("recDeltaR", &recDeltaR, "recDeltaR/F");
  outTree->Branch("recSepAngle", &recSepAngle, "recSepAngle/F");
  outTree->Branch("recEtaSingle_0", &recEtaSingle_0, "recEtaSingle_0/F");
  outTree->Branch("recEtaSingle_1", &recEtaSingle_1, "recEtaSingle_1/F");
  outTree->Branch("recChAsymm_0", &recChAsymm_0, "recChAsymm_0/F");
  outTree->Branch("recChAsymm_1", &recChAsymm_1, "recChAsymm_1/F");
  outTree->Branch("recPtSingle_0", &recPtSingle_0, "recPtSingle_0/F");
  outTree->Branch("recPtSingle_1", &recPtSingle_1, "recPtSingle_1/F");
  outTree->Branch("recE_0", &recE_0, "recE_0/F");
  outTree->Branch("recE_1", &recE_1, "recE_1/F");
  outTree->Branch("recD0overSigma_0", &recD0overSigma_0, "recD0overSigma_0/F");
  outTree->Branch("recD0overSigma_1", &recD0overSigma_1, "recD0overSigma_1/F");
  outTree->Branch("recDM_0", &recDM_0, "recDM_0/I");
  outTree->Branch("recDM_1", &recDM_1, "recDM_1/I");
  outTree->Branch("recNPfos_0", &recNPfos_0, "recNPfos_0/I");
  outTree->Branch("recNPfos_1", &recNPfos_1, "recNPfos_1/I");

  const long int nEntries = inTree->GetEntries();

  std::cout << boson << " type entries: " << nEntries << '\n';

  for (int ientry = 0; ientry < nEntries; ++ientry) {

    inTree->GetEntry(ientry);

    if (ntau_rec != 2)
      continue;

    // if (matched[0] == 0 || matched[1] == 0) continue;
    if (!((charge[0] == 1 && charge[1] == -1) ||
          (charge[0] == -1 && charge[1] == 1)))
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

    TH1F *histo_corr = nullptr;
    float bin_correction;
    TString recDecayMode[2];

    for (int t = 0, pfos_counter = 0; t != ntau_rec; ++t) {
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

      bin_correction = histo_corr->GetBinContent(histo_corr->FindBin(recE[t]));
      bin_correction > 0. ? ECorr[t] = bin_correction : ECorr[t] = 1.;

      if(recDecayMode[t] == "1pi0N" || recDecayMode[t] == "1piN") recDecayMode[t] = "1pi";

      histo_corr = nullptr;
    }

    if (!(recdecay0 == "all" && recdecay1 == "all")) {
      if (recdecay0 == "all") {
        if (recdecay1 == "had") {
          if (!(recDecayMode[1] == "1pi" || recDecayMode[0] == "1pi" ||
                recDecayMode[1] == "3pi" || recDecayMode[0] == "3pi"))
            continue;
        } else if (recdecay1 == "hade") {
          if (!(recDecayMode[1] == "1pi" || recDecayMode[0] == "1pi" ||
                recDecayMode[1] == "3pi" || recDecayMode[0] == "3pi" ||
                recDecayMode[1] == "e" || recDecayMode[0] == "e"))
            continue;
        } else if (!(recDecayMode[1] == recdecay1 ||
                     recDecayMode[0] == recdecay1))
          continue;
      } else if (recdecay1 == "all") {
        if (recdecay0 == "had") {
          if (!(recDecayMode[1] == "1pi" || recDecayMode[0] == "1pi" ||
                recDecayMode[1] == "3pi" || recDecayMode[0] == "3pi"))
            continue;
        } else if (recdecay0 == "hade") {
          if (!(recDecayMode[1] == "1pi" || recDecayMode[0] == "1pi" ||
                recDecayMode[1] == "3pi" || recDecayMode[0] == "3pi" ||
                recDecayMode[1] == "e" || recDecayMode[0] == "e"))
            continue;
        } else if (!(recDecayMode[1] == recdecay0 ||
                     recDecayMode[0] == recdecay0))
          continue;
      } else if (recdecay0 == "hade" && recdecay1 == "hade") {
        if (!((recDecayMode[1] == "1pi" || recDecayMode[1] == "3pi" ||
               recDecayMode[1] == "e") &&
              (recDecayMode[0] == "1pi" || recDecayMode[0] == "3pi" ||
               recDecayMode[0] == "e")))
          continue;
      } else if (recdecay0 == "had" && recdecay1 == "had") {
        if (!((recDecayMode[1] == "1pi" || recDecayMode[1] == "3pi") &&
              (recDecayMode[0] == "1pi" || recDecayMode[0] == "3pi")))
          continue;
      } else if (recdecay0 == "had") {
        if (!(((recDecayMode[1] == "1pi" || recDecayMode[1] == "3pi") &&
               recDecayMode[0] == recdecay1) ||
              ((recDecayMode[0] == "1pi" || recDecayMode[0] == "3pi") &&
               recDecayMode[1] == recdecay1)))
          continue;
      } else if (recdecay1 == "had") {
        if (!(((recDecayMode[1] == "1pi" || recDecayMode[1] == "3pi") &&
               recDecayMode[0] == recdecay0) ||
              ((recDecayMode[0] == "1pi" || recDecayMode[0] == "3pi") &&
               recDecayMode[1] == recdecay0)))
          continue;
      } else if (!((recDecayMode[0] == recdecay0 &&
                    recDecayMode[1] == recdecay1) ||
                   (recDecayMode[0] == recdecay1 &&
                    recDecayMode[1] == recdecay0)))
        continue;
    }

    for (int i = 0; i != 2; ++i) {
      if (recDecayMode[i] == "1pi")
        recDM[i] = 0;
      else if (recDecayMode[i] == "3pi")
        recDM[i] = 1;
      else if (recDecayMode[i] == "e")
        recDM[i] = 2;
      else if (recDecayMode[i] == "mu")
        recDM[i] = 3;
      else if (recDecayMode[i] == "other")
        recDM[i] = 4;
    }

    if (corrections == false) {
      ECorr[0] = 1.;
      ECorr[1] = 1.;
    }

    float recPxCorr[2] = {recPx[0] * ECorr[0], recPx[1] * ECorr[1]};
    float recPyCorr[2] = {recPy[0] * ECorr[0], recPy[1] * ECorr[1]};
    float recPzCorr[2] = {recPz[0] * ECorr[0], recPz[1] * ECorr[1]};
    float recPtCorr[2] = {recPt[0] * ECorr[0], recPt[1] * ECorr[1]};
    float recECorr[2] = {recE[0] * ECorr[0], recE[1] * ECorr[1]};

    if (recPtCorr[0] < pt_cut || recPtCorr[1] < pt_cut)
      continue;

    recPtSingle[0] = recPtCorr[0];
    recPtSingle[1] = recPtCorr[1];
    h_recPtSingle->Fill(recPtCorr[0]);
    h_recPtSingle->Fill(recPtCorr[1]);

    recPxPair = recPxCorr[0] + recPxCorr[1];
    recPyPair = recPyCorr[0] + recPyCorr[1];
    recPzPair = recPzCorr[0] + recPzCorr[1];

    recPtPair = std::sqrt(
        (recPxCorr[0] + recPxCorr[1]) * (recPxCorr[0] + recPxCorr[1]) +
        (recPyCorr[0] + recPyCorr[1]) * (recPyCorr[0] + recPyCorr[1]));
    recEtaPair = fabs(std::atanh(recPzPair / sqrt(recPxPair * recPxPair +
                                                  recPyPair * recPyPair +
                                 recPzPair * recPzPair)));

    /*recEtaSingle[0] =
        fabs(std::atanh(recPzCorr[0] / sqrt(recPxCorr[0] * recPxCorr[0] +
                                            recPyCorr[0] * recPyCorr[0] +
                                            recPzCorr[0] * recPzCorr[0])));*/
    recEtaSingle[0] = recEta[0];
    recEtaSingle[1] = recEta[1];

    recMInv = std::sqrt(
        (recECorr[0] + recECorr[1]) * (recECorr[0] + recECorr[1]) -
        ((recPxCorr[0] + recPxCorr[1]) * (recPxCorr[0] + recPxCorr[1]) +
         (recPyCorr[0] + recPyCorr[1]) * (recPyCorr[0] + recPyCorr[1]) +
         (recPzCorr[0] + recPzCorr[1]) * (recPzCorr[0] + recPzCorr[1])));
    recDeltaPhi = std::min(std::fabs(recPhi[0] - recPhi[1]),
                           float(2 * M_PI - std::fabs(recPhi[0] - recPhi[1])));

    recDeltaEta = std::fabs(recEta[0] - recEta[1]);
    recDeltaR =
        std::sqrt(recDeltaEta * recDeltaEta + recDeltaPhi * recDeltaPhi);
    recSepAngle =
        acos((recPxCorr[0] * recPxCorr[1] + recPyCorr[0] * recPyCorr[1] +
              recPzCorr[0] * recPzCorr[1]) /
             (sqrt(recPxCorr[0] * recPxCorr[0] + recPyCorr[0] * recPyCorr[0] +
                   recPzCorr[0] * recPzCorr[0]) *
              sqrt(recPxCorr[1] * recPxCorr[1] + recPyCorr[1] * recPyCorr[1] +
                   recPzCorr[1] * recPzCorr[1])));

    h_recPtPair->Fill(recPtPair);
    h_recMInv->Fill(recMInv);
    h_recDeltaPhi->Fill(recDeltaPhi * 180. / M_PI);
    h_recDeltaEta->Fill(recDeltaEta);
    h_recDeltaR->Fill(recDeltaR);
    h_recSepAngle->Fill(recSepAngle);
    h_recEtaPair->Fill(fabs(recEtaPair));

    h_recD0Single->Fill(std::fabs(recD0[0]));
    h_recD0Single->Fill(std::fabs(recD0[1]));
    h_recSigmaD0Single->Fill(std::fabs(recSigmaD0[0]));
    h_recSigmaD0Single->Fill(std::fabs(recSigmaD0[1]));
    recD0overSigma[0] = std::fabs(recD0[0] / recSigmaD0[0]);
    recD0overSigma[1] = std::fabs(recD0[1] / recSigmaD0[1]);
    h_recD0overSigma->Fill(recD0overSigma[0]);
    h_recD0overSigma->Fill(recD0overSigma[1]);

    float ptPfo = 0.;
    float ptLeadPfo = 0.;

    for (int t = 0, pfos_counter = 0; t != ntau_rec; ++t) {
      ptLeadPfo = 0.;
      for (int j = 0; j != recNPfos[t]; j++) {
        ptPfo = std::sqrt(pfosPx[j + pfos_counter] * pfosPx[j + pfos_counter] +
                          pfosPy[j + pfos_counter] * pfosPy[j + pfos_counter]);
        if (ptPfo > ptLeadPfo)
          ptLeadPfo = ptPfo;
      }
      recChAsymm[t] = (2. * ptLeadPfo / recPt[t]) - 1.;
      h_recChAsymm->Fill(recChAsymm[t]);

      pfos_counter += 20;
    }

    int idx_maxPt = 0, idx_minPt = 1;
    if (recPt[1] > recPt[0]) {
      idx_maxPt = 1;
      idx_minPt = 0;
    }

    recPtSingle_0 = recPtSingle[idx_maxPt];
    recPtSingle_1 = recPtSingle[idx_minPt];
    recDM_0 = recDM[idx_maxPt];
    recDM_1 = recDM[idx_minPt];
    recE_0 = recE[idx_maxPt];
    recE_1 = recE[idx_minPt];
    recEtaSingle_0 = recEtaSingle[idx_maxPt];
    recEtaSingle_1 = recEtaSingle[idx_minPt];
    recD0overSigma_0 = recD0overSigma[idx_maxPt];
    recD0overSigma_1 = recD0overSigma[idx_minPt];
    recChAsymm_0 = recChAsymm[idx_maxPt];
    recChAsymm_1 = recChAsymm[idx_minPt];
    recNPfos_0 = recNPfos[idx_maxPt];
    recNPfos_1 = recNPfos[idx_minPt];

    outTree->Fill();
  }
  inTree->Delete();
  input_file->Close();
  input_file->Delete();

  outTree->Write(boson + "tree_" + recdecay0 + "+" + recdecay1,
                 TObject::kSingleKey);
  // outTree->Delete();

  gStyle->SetOptStat(111111);
  h_recPtSingle->Scale(1. / h_recPtSingle->Integral());
  h_recPtPair->Scale(1. / h_recPtPair->Integral());
  h_recMInv->Scale(1. / h_recMInv->Integral());
  h_recDeltaPhi->Scale(1. / h_recDeltaPhi->Integral());
  h_recDeltaEta->Scale(1. / h_recDeltaEta->Integral());
  h_recDeltaR->Scale(1. / h_recDeltaR->Integral());
  h_recSepAngle->Scale(1. / h_recSepAngle->Integral());
  h_recEtaPair->Scale(1. / h_recEtaPair->Integral());
  h_recD0Single->Scale(1. / h_recD0Single->Integral());
  h_recSigmaD0Single->Scale(1. / h_recSigmaD0Single->Integral());
  h_recChAsymm->Scale(1. / h_recChAsymm->Integral());
  h_recD0overSigma->Scale(1. / h_recD0overSigma->Integral());

  TList *histo_list = new TList();
  histo_list->Add(h_recPtSingle);
  histo_list->Add(h_recPtPair);
  histo_list->Add(h_recMInv);
  histo_list->Add(h_recDeltaPhi);
  histo_list->Add(h_recDeltaEta);
  histo_list->Add(h_recDeltaR);
  histo_list->Add(h_recSepAngle);
  histo_list->Add(h_recEtaPair);
  histo_list->Add(h_recD0Single);
  histo_list->Add(h_recSigmaD0Single);
  histo_list->Add(h_recChAsymm);
  histo_list->Add(h_recD0overSigma);

  if (corrections == true) {
    histo_list->Write(boson + "_" + recdecay0 + "+" + recdecay1 + "_histograms",
                      TObject::kSingleKey);
  } else {
    histo_list->Write(boson + "_" + recdecay0 + "+" + recdecay1 +
                          "_histograms_uncorr",
                      TObject::kSingleKey);
  }
  rootfile->Close();

  {
    delete[] recE;
    delete[] recPx;
    delete[] recPy;
    delete[] recPz;
    delete[] recPt;
    delete[] recEta;
    delete[] recPhi;
    delete[] recD0;
    delete[] recSigmaD0;
    delete[] pfosPx;
    delete[] pfosPy;
    delete[] pfosPz;

    h_recPtPair->Delete();
    h_recPtSingle->Delete();
    h_recMInv->Delete();
    h_recDeltaPhi->Delete();
    h_recDeltaEta->Delete();
    h_recDeltaR->Delete();
    h_recSepAngle->Delete();
    h_recEtaPair->Delete();
    h_recD0Single->Delete();
    h_recSigmaD0Single->Delete();
    h_recChAsymm->Delete();
    h_recD0overSigma->Delete();
  }
}

void TauPairsAnalysisFake(
    const TString filenameH = "/eos/user/l/lvalla/MuColl/H2tautau10k/"
                              "output_EvalTauFinder_new.root",
    const TString filenameFake = "/eos/user/l/lvalla/MuColl/Z2jets/"
                                 "output_EvalTauFinder.root",
    const TString rootfilename =
        "/eos/user/l/lvalla/MuColl/HvsZ/Hvsjets_Histograms.root") {
  TFile *rootfile = new TFile(rootfilename, "RECREATE");
  rootfile->Close();

  HistoMCFiller(filenameH, rootfilename, "H");
  HistoMCFiller(filenameFake, rootfilename, "Z");

  bool corrections = true;
  {
    HistoRecoFiller(filenameH, rootfilename, "all", "all", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "mu", "mu", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "e", "e", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "e", "mu", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "e", "1pi", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "1pi", "1pi", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "1pi", "3pi", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "3pi", "3pi", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "had", "had", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "had", "e", "H", corrections);
    HistoRecoFiller(filenameH, rootfilename, "hade", "hade", "H", corrections);
    /*HistoRecoFiller(filenameH, rootfilename, "e", "1pi0N", "H",
    corrections); HistoRecoFiller(filenameH, rootfilename, "e", "1pi1N", "H",
    corrections); HistoRecoFiller(filenameH, rootfilename, "e", "1pi2N", "H",
    corrections); HistoRecoFiller(filenameH, rootfilename, "e", "1pi3N", "H",
    corrections);*/

    HistoRecoFiller(filenameFake, rootfilename, "all", "all", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "mu", "mu", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "e", "e", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "e", "mu", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "e", "1pi", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "1pi", "1pi", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "1pi", "3pi", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "3pi", "3pi", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "had", "had", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "had", "e", "Z", corrections);
    HistoRecoFiller(filenameFake, rootfilename, "hade", "hade", "Z",
                    corrections);
    /*HistoRecoFiller(filenameFake, rootfilename, "e", "1pi0N", "Z",
    corrections); HistoRecoFiller(filenameFake, rootfilename, "e", "1pi1N",
    "Z", corrections); HistoRecoFiller(filenameFake, rootfilename, "e",
    "1pi2N", "Z", corrections); HistoRecoFiller(filenameFake, rootfilename,
    "e", "1pi3N", "Z", corrections);*/
  }
}