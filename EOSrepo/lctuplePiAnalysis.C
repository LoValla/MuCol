#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// --- Histogram declaration

TH1F *h_EvsEta0_mcAll, *h_EvsEta1_mcAll, *h_EvsEta2_mcAll, *h_EvsEta3_mcAll,
    *h_EvsEta4_mcAll, *h_EvsEta0_Reco, *h_EvsEta1_Reco, *h_EvsEta2_Reco,
    *h_EvsEta3_Reco, *h_EvsEta4_Reco, *h_EvsEta0_Trk, *h_EvsEta1_Trk,
    *h_EvsEta2_Trk, *h_EvsEta3_Trk, *h_EvsEta4_Trk;

TH1F *h_EtavsE0_mcAll, *h_EtavsE1_mcAll, *h_EtavsE2_mcAll, *h_EtavsE3_mcAll,
    *h_EtavsE4_mcAll, *h_EtavsE5_mcAll, *h_EtavsE6_mcAll, *h_EtavsE0_Reco,
    *h_EtavsE1_Reco, *h_EtavsE2_Reco, *h_EtavsE3_Reco, *h_EtavsE4_Reco,
    *h_EtavsE5_Reco, *h_EtavsE6_Reco, *h_EtavsE0_Trk, *h_EtavsE1_Trk,
    *h_EtavsE2_Trk, *h_EtavsE3_Trk, *h_EtavsE4_Trk, *h_EtavsE5_Trk,
    *h_EtavsE6_Trk;

TEfficiency *eff_EvsEta0, *eff_EvsEta1, *eff_EvsEta2, *eff_EvsEta3,
    *eff_EvsEta4, *eff_EtavsE0, *eff_EtavsE1, *eff_EtavsE2, *eff_EtavsE3,
    *eff_EtavsE4, *eff_EtavsE5, *eff_EtavsE6;

TEfficiency *eff_trk_EvsEta0, *eff_trk_EvsEta1, *eff_trk_EvsEta2,
    *eff_trk_EvsEta3, *eff_trk_EvsEta4, *eff_trk_EtavsE0, *eff_trk_EtavsE1,
    *eff_trk_EtavsE2, *eff_trk_EtavsE3, *eff_trk_EtavsE4, *eff_trk_EtavsE5,
    *eff_trk_EtavsE6;

TGraph *g_pt_pions, *g_pt_neutrons;

// ===========================================================================

void lctuplePiAnalysis(const TString filename = lctuple_piguns.root") {
  const TString plotsfolder = "plots"; // plots folder name

  double etamax = 2.5;
  double etamin = 0.;
  int nbins_eta = 100;
  double Emax = 300.;
  double Emin = 1.;
  int nbins_E = 20;
  int nlogbins_E = 16;
  // double logbin_edges[nlogbins_E + 1];

  // float binsemiwidth_E = 0.5 * (Emax - Emin) / nbins_E;
  /*
  for (int i = 0; i <= nlogbins_E; i++) {
    logbin_edges[i] =
        pow(10, TMath::Log10(Emin) + (TMath::Log10(Emax) - TMath::Log10(Emin)) /
                                         double(nlogbins_E) * double(i));
    std::cout << logbin_edges[i] << '\n';
  }*/

  double logbin_edges[17] = {5.,   10.,  20.,  30.,  40.,  50.,
                             60.,  80.,  100., 120., 140., 160.,
                             180., 200., 250., 300., 500.};
  /*double logbin_edges[52] = {
      5.,   10.,  15.,  20.,  25.,  30.,  35.,  40.,  45.,  50.,  55.,
      60.,  65.,  70.,  75.,  80.,  85.,  90.,  95.,  100., 105., 110.,
      115., 120., 125., 130., 135., 140., 145., 150., 155., 160., 165.,
      170., 175., 180., 185., 190., 195., 200., 210., 220., 230., 240.,
      250., 260., 270., 280., 290., 300., 400., 500.};*/

  // Histogram initialization
  {
    h_EvsEta0_mcAll =
        new TH1F("h_EvsEta0_mcAll",
                 "Pion E events (All, 0 < |#eta| < 0.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta1_mcAll =
        new TH1F("h_EvsEta1_mcAll",
                 "Pion events (All, 0.5 < |#eta| < 1); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta2_mcAll =
        new TH1F("h_EvsEta2_mcAll",
                 "Pion events (All, 1 < |#eta| < 1.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta3_mcAll =
        new TH1F("h_EvsEta3_mcAll",
                 "Pion events (All, 1.5 < |#eta| < 2); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta4_mcAll =
        new TH1F("h_EvsEta4_mcAll",
                 "Pion events (All, 2 < |#eta| < 2.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);

    h_EvsEta0_Reco =
        new TH1F("h_EvsEta0_Reco",
                 "Pion events (Hits, 0 < |#eta| < 0.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta1_Reco =
        new TH1F("h_EvsEta1_Reco",
                 "Pion events (Hits, 0.5 < |#eta| < 1); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta2_Reco =
        new TH1F("h_EvsEta2_Reco",
                 "Pion events (Hits, 1 < |#eta| < 1.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta3_Reco =
        new TH1F("h_EvsEta3_Reco",
                 "Pion events (Hits, 1.5 < |#eta| < 2); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta4_Reco =
        new TH1F("h_EvsEta4_Reco",
                 "Pion events (Hits, 2 < |#eta| < 2.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta0_Trk =
        new TH1F("h_EvsEta0_Trk",
                 "Pion events (Hits, 0 < |#eta| < 0.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta1_Trk =
        new TH1F("h_EvsEta1_Trk",
                 "Pion events (Hits, 0.5 < |#eta| < 1); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta2_Trk =
        new TH1F("h_EvsEta2_Trk",
                 "Pion events (Hits, 1 < |#eta| < 1.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta3_Trk =
        new TH1F("h_EvsEta3_Trk",
                 "Pion events (Hits, 1.5 < |#eta| < 2); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);
    h_EvsEta4_Trk =
        new TH1F("h_EvsEta4_Trk",
                 "Pion events (Hits, 2 < |#eta| < 2.5); E_{gen} [GeV]; events",
                 nlogbins_E, logbin_edges);

    h_EtavsE0_mcAll =
        new TH1F("h_EtavsE0_mcAll", "Pion events (All, E0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_mcAll =
        new TH1F("h_EtavsE1_mcAll", "Pion events (All, E1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_mcAll =
        new TH1F("h_EtavsE2_mcAll", "Pion events (All, E2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_mcAll =
        new TH1F("h_EtavsE3_mcAll", "Pion events (All, E3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_mcAll =
        new TH1F("h_EtavsE4_mcAll", "Pion events (All, E4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_mcAll =
        new TH1F("h_EtavsE5_mcAll", "Pion events (All, E5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE0_Reco =
        new TH1F("h_EtavsE0_Reco", "Pion events (Hits, E0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_Reco =
        new TH1F("h_EtavsE1_Reco", "Pion events (Hits, E1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_Reco =
        new TH1F("h_EtavsE2_Reco", "Pion events (Hits, E2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_Reco =
        new TH1F("h_EtavsE3_Reco", "Pion events (Hits, E3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_Reco =
        new TH1F("h_EtavsE4_Reco", "Pion events (Hits, E4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_Reco =
        new TH1F("h_EtavsE5_Reco", "Pion events (Hits, E5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE0_Trk =
        new TH1F("h_EtavsE0_Trk", "Pion events (Hits, E0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_Trk =
        new TH1F("h_EtavsE1_Trk", "Pion events (Hits, E1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_Trk =
        new TH1F("h_EtavsE2_Trk", "Pion events (Hits, E2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_Trk =
        new TH1F("h_EtavsE3_Trk", "Pion events (Hits, E3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_Trk =
        new TH1F("h_EtavsE4_Trk", "Pion events (Hits, E4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_Trk =
        new TH1F("h_EtavsE5_Trk", "Pion events (Hits, E5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    g_pt_pions = new TGraph();
    g_pt_neutrons = new TGraph();
  }

  //  --- Open the ntuple file and get the tree

  TFile *input_file = new TFile(filename.Data(), "read");

  TTree *myLCTuple = (TTree *)input_file->Get("MyLCTuple");

  // --- Loop over the ntuple entries

  // --- MC particles
  int n_mcp;
  int *mcPdg_ = new int[1500000];
  int *mcGenCode = new int[1500000];
  float *mcPx = new float[1500000];
  float *mcPy = new float[1500000];
  float *mcPz = new float[1500000];
  float *mcEne_ = new float[1500000];
  float *mcQ = new float[1500000];

  myLCTuple->SetBranchAddress("nmcp", &n_mcp);
  myLCTuple->SetBranchAddress("mcpdg", mcPdg_);
  myLCTuple->SetBranchAddress("mcgst", mcGenCode);
  myLCTuple->SetBranchAddress("mcmox", mcPx);
  myLCTuple->SetBranchAddress("mcmoy", mcPy);
  myLCTuple->SetBranchAddress("mcmoz", mcPz);
  myLCTuple->SetBranchAddress("mcene", mcEne_);
  myLCTuple->SetBranchAddress("mccha", mcQ);

  // additional MC variables
  float *mcPt = new float[10000];
  float *mcPmag = new float[10000];
  float *mcEne = new float[10000];
  float *mcEta = new float[10000];
  float *mcEtaAbs = new float[10000];
  float *mcPhi = new float[10000];
  float *mcPdg = new float[10000];

  // --- RECO particles
  int n_reco;
  int *reco_pdg = new int[10000];
  float *reco_px = new float[10000];
  float *reco_py = new float[10000];
  float *reco_pz = new float[10000];
  float *reco_ene = new float[10000];
  float *reco_q = new float[10000];

  myLCTuple->SetBranchAddress("nrec", &n_reco);
  myLCTuple->SetBranchAddress("rctyp", reco_pdg);
  myLCTuple->SetBranchAddress("rcmox", reco_px);
  myLCTuple->SetBranchAddress("rcmoy", reco_py);
  myLCTuple->SetBranchAddress("rcmoz", reco_pz);
  myLCTuple->SetBranchAddress("rcene", reco_ene);
  myLCTuple->SetBranchAddress("rccha", reco_q);

  // additional reco variables
  float *reco_pmag = new float[10000];
  float *reco_eta = new float[10000];
  float *reco_phi = new float[10000];

  // track variables
  int ts_idx, n_trst;
  float *ts_phi = new float[10000];
  float *ts_tnl = new float[10000];
  float *ts_ome = new float[40000];
  float *ts_cov = new float[150000];
  myLCTuple->SetBranchAddress("ntrst", &n_trst);
  myLCTuple->SetBranchAddress("trsip", &ts_idx);
  myLCTuple->SetBranchAddress("tsphi", ts_phi);
  myLCTuple->SetBranchAddress("tstnl", ts_tnl);
  myLCTuple->SetBranchAddress("tsome", ts_ome);
  myLCTuple->SetBranchAddress("tscov", ts_cov);

  float *ts_theta = new float[10000];
  float *ts_eta = new float[10000];
  float *ts_pt = new float[10000];
  float *ts_ptunc = new float[10000];
  float *ts_p = new float[10000];
  float *ts_punc = new float[10000];

  int n_mcp_stable = 0;

  const long int nEntries = myLCTuple->GetEntries();
  for (int ientry = 0; ientry < nEntries; ++ientry) {
    myLCTuple->GetEntry(ientry);

    n_mcp_stable = 0;

    // --- loop over the Monte Carlo particles
    for (int i = 0; i < n_mcp; ++i) {
      // --- keep only the stable particles and save the information
      if (mcGenCode[i] != 1)
        continue;

      mcPt[n_mcp_stable] = TMath::Sqrt(mcPx[i] * mcPx[i] + mcPy[i] * mcPy[i]);
      mcPmag[n_mcp_stable] = TMath::Sqrt(mcPx[i] * mcPx[i] + mcPy[i] * mcPy[i] +
                                         mcPz[i] * mcPz[i]);
      mcEne[n_mcp_stable] = mcEne_[i];
      mcEta[n_mcp_stable] =
          0.5 * TMath::Log((mcPmag[i] + mcPz[i]) / (mcPmag[i] - mcPz[i]));
      mcEtaAbs[n_mcp_stable] = std::fabs(mcEta[n_mcp_stable]);
      mcPhi[n_mcp_stable] = TMath::ATan2(mcPy[i], mcPx[i]);
      mcPdg[n_mcp_stable] = mcPdg_[i];

      if (mcEtaAbs[n_mcp_stable] >= 0. && mcEtaAbs[n_mcp_stable] <= 0.5)
        h_EvsEta0_mcAll->Fill(mcEne[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 0.5 && mcEtaAbs[n_mcp_stable] <= 1.)
        h_EvsEta1_mcAll->Fill(mcEne[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 1. && mcEtaAbs[n_mcp_stable] <= 1.5)
        h_EvsEta2_mcAll->Fill(mcEne[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 1.5 && mcEtaAbs[n_mcp_stable] <= 2.)
        h_EvsEta3_mcAll->Fill(mcEne[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 2. && mcEtaAbs[n_mcp_stable] <= 2.5)
        h_EvsEta4_mcAll->Fill(mcEne[n_mcp_stable]);

      if (mcEne[n_mcp_stable] >= 0. && mcEne[n_mcp_stable] < 50.)
        h_EtavsE0_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcEne[n_mcp_stable] >= 50. && mcEne[n_mcp_stable] < 100.)
        h_EtavsE1_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcEne[n_mcp_stable] >= 100. && mcEne[n_mcp_stable] < 150.)
        h_EtavsE2_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcEne[n_mcp_stable] >= 150. && mcEne[n_mcp_stable] < 200.)
        h_EtavsE3_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcEne[n_mcp_stable] >= 200. && mcEne[n_mcp_stable] <= 300.)
        h_EtavsE4_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcEne[n_mcp_stable] >= 300.)
        h_EtavsE5_mcAll->Fill(mcEtaAbs[n_mcp_stable]);

      ++n_mcp_stable;
    } // i loop

    // --- loop over the reconstructed particles

    for (int j = 0; j < n_reco; ++j) {

      // check if matched
      if (std::abs(reco_pdg[j]) != 211) {
        continue;
      }

      reco_pmag[j] =
          TMath::Sqrt(reco_px[j] * reco_px[j] + reco_py[j] * reco_py[j] +
                      reco_pz[j] * reco_pz[j]);
      reco_eta[j] = 0.5 * TMath::Log((reco_pmag[j] + reco_pz[j]) /
                                     (reco_pmag[j] - reco_pz[j]));
      reco_phi[j] = TMath::ATan2(reco_py[j], reco_px[j]);

      for (int k = 0; k < n_mcp_stable; ++k) {
        if (reco_pdg[j] == mcPdg[k]) {
          float deltaPhi =
              std::min(std::fabs(reco_phi[j] - mcPhi[k]),
                       float(2 * M_PI - std::fabs(reco_phi[j] - mcPhi[k])));

          if ( // std::fabs(reco_ene[j] - mcEne[k]) < 0.1 * mcEne[k]
              true) {
            if (TMath::Sqrt(pow((reco_eta[j] - mcEta[k]), 2) +
                            pow(deltaPhi, 2)) < 0.01) {

              if (mcEtaAbs[k] >= 0. && mcEtaAbs[k] <= 0.5)
                h_EvsEta0_Reco->Fill(mcEne[k]);
              else if (mcEtaAbs[k] >= 0.5 && mcEtaAbs[k] <= 1.)
                h_EvsEta1_Reco->Fill(mcEne[k]);
              else if (mcEtaAbs[k] >= 1. && mcEtaAbs[k] <= 1.5)
                h_EvsEta2_Reco->Fill(mcEne[k]);
              else if (mcEtaAbs[k] >= 1.5 && mcEtaAbs[k] <= 2.)
                h_EvsEta3_Reco->Fill(mcEne[k]);
              else if (mcEtaAbs[k] >= 2. && mcEtaAbs[k] <= 2.5)
                h_EvsEta4_Reco->Fill(mcEne[k]);

              if (mcEne[k] >= 0. && mcEne[k] < 50.)
                h_EtavsE0_Reco->Fill(mcEtaAbs[k]);
              else if (mcEne[k] >= 50. && mcEne[k] < 100.)
                h_EtavsE1_Reco->Fill(mcEtaAbs[k]);
              else if (mcEne[k] >= 100. && mcEne[k] < 150.)
                h_EtavsE2_Reco->Fill(mcEtaAbs[k]);
              else if (mcEne[k] >= 150. && mcEne[k] < 200.)
                h_EtavsE3_Reco->Fill(mcEtaAbs[k]);
              else if (mcEne[k] >= 200. && mcEne[k] <= 300.)
                h_EtavsE4_Reco->Fill(mcEtaAbs[k]);
              else if (mcEne[k] >= 300.)
                h_EtavsE5_Reco->Fill(mcEtaAbs[k]);

              break;
            }
          }
        }
      }
    }

    /*std::cout << '\n';
    std::cout << "n_trks = " << n_trst / 4;
    std::cout << '\n';*/

    for (int t = ts_idx; t < n_trst; t += 4) {

      ts_theta[t] = TMath::ATan(1. / ts_tnl[t]);
      if (ts_theta[t] < 0.) {
        ts_theta[t] += M_PI;
      }
      ts_eta[t] = -TMath::Log(TMath::Tan(ts_theta[t] / 2.));
      ts_pt[t] =
          3.57 * 1.6022e-19 / (std::fabs(ts_ome[t]) * 1.e3 * 5.344286e-19);
      ts_ptunc[t] = std::fabs(sqrt(ts_cov[t * 15 + 5]) / ts_ome[t]);
      ts_p[t] = std::fabs(ts_pt[t] / sin(ts_theta[t]));
      ts_punc[t] =
          std::sqrt(pow(std::cos(ts_theta[t]) * std::sqrt(ts_cov[t * 15 + 14]) /
                            (1 + pow(ts_tnl[t], -2)) / std::sin(ts_theta[t]),
                        2) +
                    pow(ts_ptunc[t], 2));

      for (int k = 0; k < n_mcp_stable; ++k) {

        float deltaPhi =
            std::min(std::fabs(ts_phi[t] - mcPhi[k]),
                     float(2 * M_PI - std::fabs(ts_phi[t] - mcPhi[k])));
        if ((ts_p[t] - mcPmag[k]) < 0.1 * mcPmag[k]) {
          if (TMath::Sqrt(pow((ts_eta[t] - mcEta[k]), 2) + pow(deltaPhi, 2)) <
              0.01) {

            if (reco_pdg[0] == 211) {
              // h_pt_pions->Fill(ts_pt[k]);
              // h_pt_unc_pions->Fill(ts_ptunc[k]);

              g_pt_pions->AddPoint(ts_pt[k], ts_ptunc[k]);

            } else if (reco_pdg[0] == 2112) {
              // h_pt_neutrons->Fill(ts_pt[k]);
              // h_pt_unc_neutrons->Fill(ts_ptunc[k]);

              g_pt_neutrons->AddPoint(ts_pt[k], ts_ptunc[k]);
            }

            if (mcEtaAbs[k] >= 0. && mcEtaAbs[k] <= 0.5)
              h_EvsEta0_Trk->Fill(mcEne[k]);
            else if (mcEtaAbs[k] >= 0.5 && mcEtaAbs[k] <= 1.)
              h_EvsEta1_Trk->Fill(mcEne[k]);
            else if (mcEtaAbs[k] >= 1. && mcEtaAbs[k] <= 1.5)
              h_EvsEta2_Trk->Fill(mcEne[k]);
            else if (mcEtaAbs[k] >= 1.5 && mcEtaAbs[k] <= 2.)
              h_EvsEta3_Trk->Fill(mcEne[k]);
            else if (mcEtaAbs[k] >= 2. && mcEtaAbs[k] <= 2.5)
              h_EvsEta4_Trk->Fill(mcEne[k]);

            if (mcEne[k] >= 0. && mcEne[k] < 50.)
              h_EtavsE0_Trk->Fill(mcEtaAbs[k]);
            else if (mcEne[k] >= 50. && mcEne[k] < 100.)
              h_EtavsE1_Trk->Fill(mcEtaAbs[k]);
            else if (mcEne[k] >= 100. && mcEne[k] < 150.)
              h_EtavsE2_Trk->Fill(mcEtaAbs[k]);
            else if (mcEne[k] >= 150. && mcEne[k] < 200.)
              h_EtavsE3_Trk->Fill(mcEtaAbs[k]);
            else if (mcEne[k] >= 200. && mcEne[k] <= 300.)
              h_EtavsE4_Trk->Fill(mcEtaAbs[k]);
            else if (mcEne[k] >= 300.)
              h_EtavsE5_Trk->Fill(mcEtaAbs[k]);

            break;
          }
        }
      } // k loop
    }   // t loop
  }     // ientry loop

  // draw the histograms
  {

    if (TEfficiency::CheckConsistency(*h_EvsEta0_Reco, *h_EvsEta0_mcAll)) {
      auto c_EvsEta0_eff =
          new TCanvas("c_EvsEta0_eff", "c_EvsEta0_eff", 900, 800);
      eff_EvsEta0 = new TEfficiency(*h_EvsEta0_Reco, *h_EvsEta0_mcAll);
      // eff_EvsEta0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta0->SetConfidenceLevel(0.68);
      eff_EvsEta0->Draw("APL");
      gPad->SetLogx();
      eff_EvsEta0->SetTitle("Pion ID efficiency (0 < |#eta_{gen}| < 0.5); "
                            "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta0->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta0_eff->SaveAs(plotsfolder + "/c_EvsEta0_eff.png");
      c_EvsEta0_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta1_Reco, *h_EvsEta1_mcAll)) {
      auto c_EvsEta1_eff =
          new TCanvas("c_EvsEta1_eff", "c_EvsEta1_eff", 900, 800);
      eff_EvsEta1 = new TEfficiency(*h_EvsEta1_Reco, *h_EvsEta1_mcAll);
      // eff_EvsEta1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta1->SetConfidenceLevel(0.68);
      eff_EvsEta1->Draw("APL");
      eff_EvsEta1->SetTitle("Pion ID efficiency (0.5 < |#eta_{gen}| < 1); "
                            "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta1->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      gPad->SetLogx();
      c_EvsEta1_eff->SaveAs(plotsfolder + "/c_EvsEta1_eff.png");
      c_EvsEta1_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta2_Reco, *h_EvsEta2_mcAll)) {
      auto c_EvsEta2_eff =
          new TCanvas("c_EvsEta2_eff", "c_EvsEta2_eff", 900, 800);
      eff_EvsEta2 = new TEfficiency(*h_EvsEta2_Reco, *h_EvsEta2_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta2->SetConfidenceLevel(0.68);
      eff_EvsEta2->Draw("APL");
      gPad->SetLogx();
      eff_EvsEta2->SetTitle("Pion ID efficiency (1 < |#eta_{gen}| < 1.5); "
                            "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta2->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta2_eff->SaveAs(plotsfolder + "/c_EvsEta2_eff.png");
      c_EvsEta2_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta3_Reco, *h_EvsEta3_mcAll)) {
      auto c_EvsEta3_eff =
          new TCanvas("c_EvsEta3_eff", "c_EvsEta3_eff", 900, 800);
      eff_EvsEta3 = new TEfficiency(*h_EvsEta3_Reco, *h_EvsEta3_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta3->SetConfidenceLevel(0.68);
      eff_EvsEta3->Draw("APL");
      gPad->SetLogx();
      eff_EvsEta3->SetTitle("Pion ID efficiency (1.5 < |#eta_{gen}| < 2); "
                            "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta3->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta3_eff->SaveAs(plotsfolder + "/c_EvsEta3_eff.png");
      c_EvsEta3_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta4_Reco, *h_EvsEta4_mcAll)) {
      auto c_EvsEta4_eff =
          new TCanvas("c_EvsEta4_eff", "c_EvsEta4_eff", 900, 800);
      eff_EvsEta4 = new TEfficiency(*h_EvsEta4_Reco, *h_EvsEta4_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta4->SetConfidenceLevel(0.68);
      eff_EvsEta4->Draw("APL");
      gPad->SetLogx();
      eff_EvsEta4->SetTitle("Pion ID efficiency (2 < |#eta_{gen}| < 2.5); "
                            "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta4->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta4_eff->SaveAs(plotsfolder + "/c_EvsEta4_eff.png");
      c_EvsEta4_eff->Close();
    }

    for (int b = 0; b != h_EtavsE0_Reco->GetNbinsX() + 2; ++b) {
      std::cout << "h_EtavsE0_Reco (bin " << b
                << ") = " << h_EtavsE0_Reco->GetBinContent(b) << '\n';
      std::cout << "h_EtavsE0_mcAll (bin " << b
                << ") = " << h_EtavsE0_mcAll->GetBinContent(b) << '\n';
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE0_Reco, *h_EtavsE0_mcAll)) {
      auto c_EtavsE0_eff =
          new TCanvas("c_EtavsE0_eff", "c_EtavsE0_eff", 900, 800);
      eff_EtavsE0 = new TEfficiency(*h_EtavsE0_Reco, *h_EtavsE0_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE0->SetConfidenceLevel(0.68);
      eff_EtavsE0->Draw("APL");
      eff_EtavsE0->SetTitle("Pion ID efficiency (0 < E_{gen} < 50 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE0->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE0_eff);
      gPad->Update();
      c_EtavsE0_eff->SaveAs(plotsfolder + "/c_EtavsE0_eff.png");
      c_EtavsE0_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE1_Reco, *h_EtavsE1_mcAll)) {
      auto c_EtavsE1_eff =
          new TCanvas("c_EtavsE1_eff", "c_EtavsE1_eff", 900, 800);
      eff_EtavsE1 = new TEfficiency(*h_EtavsE1_Reco, *h_EtavsE1_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE1->SetConfidenceLevel(0.68);
      eff_EtavsE1->Draw("APL");
      eff_EtavsE1->SetTitle("Pion ID efficiency (50 < E_{gen} < 100 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE1->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE1_eff);
      gPad->Update();
      c_EtavsE1_eff->SaveAs(plotsfolder + "/c_EtavsE1_eff.png");
      c_EtavsE1_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE2_Reco, *h_EtavsE2_mcAll)) {
      auto c_EtavsE2_eff =
          new TCanvas("c_EtavsE2_eff", "c_EtavsE2_eff", 900, 800);
      eff_EtavsE2 = new TEfficiency(*h_EtavsE2_Reco, *h_EtavsE2_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE2->SetConfidenceLevel(0.68);
      eff_EtavsE2->Draw("APL");
      eff_EtavsE2->SetTitle("Pion ID efficiency (100 < E_{gen} < 150 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE2->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE2_eff);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE2_eff->SaveAs(plotsfolder + "/c_EtavsE2_eff.png");
      c_EtavsE2_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE3_Reco, *h_EtavsE3_mcAll)) {
      auto c_EtavsE3_eff =
          new TCanvas("c_EtavsE3_eff", "c_EtavsE3_eff", 900, 800);
      eff_EtavsE3 = new TEfficiency(*h_EtavsE3_Reco, *h_EtavsE3_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE3->SetConfidenceLevel(0.68);
      eff_EtavsE3->Draw("APL");
      eff_EtavsE3->SetTitle("Pion ID efficiency (150 < E_{gen} < 200 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE3->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE3_eff);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE3_eff->SaveAs(plotsfolder + "/c_EtavsE3_eff.png");
      c_EtavsE3_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE4_Reco, *h_EtavsE4_mcAll)) {
      auto c_EtavsE4_eff =
          new TCanvas("c_EtavsE4_eff", "c_EtavsE4_eff", 900, 800);
      eff_EtavsE4 = new TEfficiency(*h_EtavsE4_Reco, *h_EtavsE4_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE4->SetConfidenceLevel(0.68);
      eff_EtavsE4->Draw("APL");
      eff_EtavsE4->SetTitle("Pion ID efficiency (200 < E_{gen} < 300 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE4->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE4_eff);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE4_eff->SaveAs(plotsfolder + "/c_EtavsE4_eff.png");
      c_EtavsE4_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE5_Reco, *h_EtavsE5_mcAll)) {
      auto c_EtavsE5_eff =
          new TCanvas("c_EtavsE5_eff", "c_EtavsE5_eff", 900, 800);
      eff_EtavsE5 = new TEfficiency(*h_EtavsE5_Reco, *h_EtavsE5_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE5->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE5->SetConfidenceLevel(0.68);
      eff_EtavsE5->Draw("APL");
      eff_EtavsE5->SetTitle("Pion ID efficiency (300 < E_{gen} < 500 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_EtavsE5->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE5_eff);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE5_eff->SaveAs(plotsfolder + "/c_EtavsE5_eff.png");
      c_EtavsE5_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta0_Trk, *h_EvsEta0_mcAll)) {
      auto c_EvsEta0_eff_trk =
          new TCanvas("c_EvsEta0_eff_trk", "c_EvsEta0_eff_trk", 900, 800);
      eff_trk_EvsEta0 = new TEfficiency(*h_EvsEta0_Trk, *h_EvsEta0_mcAll);
      // eff_trk_EvsEta0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EvsEta0->SetConfidenceLevel(0.68);
      eff_trk_EvsEta0->Draw("APL");
      gPad->SetLogx();
      eff_trk_EvsEta0->SetTitle("Pion ID efficiency (0 < |#eta_{gen}| < 0.5); "
                                "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EvsEta0->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta0_eff_trk->SaveAs(plotsfolder + "/c_EvsEta0_eff_trk.png");
      c_EvsEta0_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta1_Trk, *h_EvsEta1_mcAll)) {
      auto c_EvsEta1_eff_trk =
          new TCanvas("c_EvsEta1_eff_trk", "c_EvsEta1_eff_trk", 900, 800);
      eff_trk_EvsEta1 = new TEfficiency(*h_EvsEta1_Trk, *h_EvsEta1_mcAll);
      // eff_trk_EvsEta1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EvsEta1->SetConfidenceLevel(0.68);
      eff_trk_EvsEta1->Draw("APL");
      eff_trk_EvsEta1->SetTitle("Pion ID efficiency (0.5 < |#eta_{gen}| < 1); "
                                "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EvsEta1->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      gPad->SetLogx();
      c_EvsEta1_eff_trk->SaveAs(plotsfolder + "/c_EvsEta1_eff_trk.png");
      c_EvsEta1_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta2_Trk, *h_EvsEta2_mcAll)) {
      auto c_EvsEta2_eff_trk =
          new TCanvas("c_EvsEta2_eff_trk", "c_EvsEta2_eff_trk", 900, 800);
      eff_trk_EvsEta2 = new TEfficiency(*h_EvsEta2_Trk, *h_EvsEta2_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EvsEta2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EvsEta2->SetConfidenceLevel(0.68);
      eff_trk_EvsEta2->Draw("APL");
      gPad->SetLogx();
      eff_trk_EvsEta2->SetTitle("Pion ID efficiency (1 < |#eta_{gen}| < 1.5); "
                                "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EvsEta2->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta2_eff_trk->SaveAs(plotsfolder + "/c_EvsEta2_eff_trk.png");
      c_EvsEta2_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta3_Trk, *h_EvsEta3_mcAll)) {
      auto c_EvsEta3_eff_trk =
          new TCanvas("c_EvsEta3_eff_trk", "c_EvsEta3_eff_trk", 900, 800);
      eff_trk_EvsEta3 = new TEfficiency(*h_EvsEta3_Trk, *h_EvsEta3_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EvsEta3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EvsEta3->SetConfidenceLevel(0.68);
      eff_trk_EvsEta3->Draw("APL");
      gPad->SetLogx();
      eff_trk_EvsEta3->SetTitle("Pion ID efficiency (1.5 < |#eta_{gen}| < 2); "
                                "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EvsEta3->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta3_eff_trk->SaveAs(plotsfolder + "/c_EvsEta3_eff_trk.png");
      c_EvsEta3_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta4_Trk, *h_EvsEta4_mcAll)) {
      auto c_EvsEta4_eff_trk =
          new TCanvas("c_EvsEta4_eff_trk", "c_EvsEta4_eff_trk", 900, 800);
      eff_trk_EvsEta4 = new TEfficiency(*h_EvsEta4_Trk, *h_EvsEta4_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EvsEta4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EvsEta4->SetConfidenceLevel(0.68);
      eff_trk_EvsEta4->Draw("APL");
      gPad->SetLogx();
      eff_trk_EvsEta4->SetTitle("Pion ID efficiency (2 < |#eta_{gen}| < 2.5); "
                                "E_{gen} [GeV]; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EvsEta4->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta4_eff_trk->SaveAs(plotsfolder + "/c_EvsEta4_eff_trk.png");
      c_EvsEta4_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE0_Trk, *h_EtavsE0_mcAll)) {
      auto c_EtavsE0_eff_trk =
          new TCanvas("c_EtavsE0_eff_trk", "c_EtavsE0_eff_trk", 900, 800);
      eff_trk_EtavsE0 = new TEfficiency(*h_EtavsE0_Trk, *h_EtavsE0_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE0->SetConfidenceLevel(0.68);
      eff_trk_EtavsE0->Draw("APL");
      eff_trk_EtavsE0->SetTitle("Pion ID efficiency (0 < E_{gen} < 50 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE0->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE0_eff_trk);
      gPad->Update();
      c_EtavsE0_eff_trk->SaveAs(plotsfolder + "/c_EtavsE0_eff_trk.png");
      c_EtavsE0_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE1_Trk, *h_EtavsE1_mcAll)) {
      auto c_EtavsE1_eff_trk =
          new TCanvas("c_EtavsE1_eff_trk", "c_EtavsE1_eff_trk", 900, 800);
      eff_trk_EtavsE1 = new TEfficiency(*h_EtavsE1_Trk, *h_EtavsE1_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE1->SetConfidenceLevel(0.68);
      eff_trk_EtavsE1->Draw("APL");
      eff_trk_EtavsE1->SetTitle("Pion ID efficiency (50 < E_{gen} < 100 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE1->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE1_eff_trk);
      gPad->Update();
      c_EtavsE1_eff_trk->SaveAs(plotsfolder + "/c_EtavsE1_eff_trk.png");
      c_EtavsE1_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE2_Trk, *h_EtavsE2_mcAll)) {
      auto c_EtavsE2_eff_trk =
          new TCanvas("c_EtavsE2_eff_trk", "c_EtavsE2_eff_trk", 900, 800);
      eff_trk_EtavsE2 = new TEfficiency(*h_EtavsE2_Trk, *h_EtavsE2_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE2->SetConfidenceLevel(0.68);
      eff_trk_EtavsE2->Draw("APL");
      eff_trk_EtavsE2->SetTitle("Pion ID efficiency (100 < E_{gen} < 150 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE2->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE2_eff_trk);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE2_eff_trk->SaveAs(plotsfolder + "/c_EtavsE2_eff_trk.png");
      c_EtavsE2_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE3_Trk, *h_EtavsE3_mcAll)) {
      auto c_EtavsE3_eff_trk =
          new TCanvas("c_EtavsE3_eff_trk", "c_EtavsE3_eff_trk", 900, 800);
      eff_trk_EtavsE3 = new TEfficiency(*h_EtavsE3_Trk, *h_EtavsE3_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE3->SetConfidenceLevel(0.68);
      eff_trk_EtavsE3->Draw("APL");
      eff_trk_EtavsE3->SetTitle("Pion ID efficiency (150 < E_{gen} < 200 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE3->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE3_eff_trk);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE3_eff_trk->SaveAs(plotsfolder + "/c_EtavsE3_eff_trk.png");
      c_EtavsE3_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE4_Trk, *h_EtavsE4_mcAll)) {
      auto c_EtavsE4_eff_trk =
          new TCanvas("c_EtavsE4_eff_trk", "c_EtavsE4_eff_trk", 900, 800);
      eff_trk_EtavsE4 = new TEfficiency(*h_EtavsE4_Trk, *h_EtavsE4_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE4->SetConfidenceLevel(0.68);
      eff_trk_EtavsE4->Draw("APL");
      eff_trk_EtavsE4->SetTitle("Pion ID efficiency (200 < E_{gen} < 300 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE4->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE4_eff_trk);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE4_eff_trk->SaveAs(plotsfolder + "/c_EtavsE4_eff_trk.png");
      c_EtavsE4_eff_trk->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE5_Trk, *h_EtavsE5_mcAll)) {
      auto c_EtavsE5_eff_trk =
          new TCanvas("c_EtavsE5_eff_trk", "c_EtavsE5_eff_trk", 900, 800);
      eff_trk_EtavsE5 = new TEfficiency(*h_EtavsE5_Trk, *h_EtavsE5_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE5->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE5->SetConfidenceLevel(0.68);
      eff_trk_EtavsE5->Draw("APL");
      eff_trk_EtavsE5->SetTitle("Pion ID efficiency (300 < E_{gen} < 500 GeV); "
                                "|#eta_{gen}|; Efficiency");
      gPad->Update();
      auto graph = eff_trk_EtavsE5->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      // plot_boundaries(c_EtavsE5_eff_trk);
      gPad->Update();
      graph->SetMaximum(1.);
      c_EtavsE5_eff_trk->SaveAs(plotsfolder + "/c_EtavsE5_eff_trk.png");
      c_EtavsE5_eff_trk->Close();
    }

    TCanvas *c_ptvssigma =
        new TCanvas("c_mom_vs_unc_ts", "c_mom_unc_ts", 1100, 500);
    TGaxis::SetMaxDigits(4);
    c_ptvssigma->Divide(2, 1);
    c_ptvssigma->cd(1);
    gStyle->SetOptStat(111111);
    g_pt_pions->SetTitle(
        "p_{T} tracks (pions);p_{T} [GeV];#sigma_{p_{T}}/^{}p_{T}");
    g_pt_pions->GetYaxis()->SetTitleOffset(0.7);
    g_pt_pions->SetMarkerStyle(26);
    g_pt_pions->SetMarkerSize(0.3);
    g_pt_pions->SetMarkerColor(kBlue);
    gPad->SetLogy();
    // g_pt_pions->GetXaxis()->SetLimits();
    g_pt_pions->SetMaximum(0.5);
    g_pt_pions->Draw("APL");
    c_ptvssigma->cd(2);
    gStyle->SetOptStat(111111);
    g_pt_neutrons->SetTitle(
        "p_{T} tracks (no pions);p_{T} [GeV];#sigma_{p_{T}}/^{}p_{T}");
    g_pt_neutrons->GetYaxis()->SetTitleOffset(0.7);
    g_pt_neutrons->SetMarkerStyle(26);
    g_pt_neutrons->SetMarkerSize(0.3);
    g_pt_neutrons->SetMarkerColor(kRed);
    gPad->SetLogy();
    // g_pt_neutrons->GetXaxis()->SetLimits();
    g_pt_neutrons->SetMaximum(0.5);
    g_pt_neutrons->Draw("APL");
    c_ptvssigma->SaveAs(plotsfolder + "/c_scatter.pdf");
    c_ptvssigma->Close();
  }

  // Clean up the heap
  {
    delete[] mcPdg;
    delete[] mcGenCode;
    delete[] mcPx;
    delete[] mcPy;
    delete[] mcPz;
    delete[] mcQ;

    delete[] reco_pdg;
    delete[] reco_px;
    delete[] reco_py;
    delete[] reco_pz;
    delete[] reco_ene;
    delete[] reco_q;

    delete[] mcPt;
    delete[] mcPmag;
    delete[] mcEta;
    delete[] mcPhi;
    delete[] mcEne;
    delete[] mcPdg;
    delete[] reco_eta;
    delete[] reco_phi;

    h_EvsEta0_mcAll->Delete();
    h_EvsEta1_mcAll->Delete();
    h_EvsEta2_mcAll->Delete();
    h_EvsEta3_mcAll->Delete();
    h_EvsEta4_mcAll->Delete();
    h_EvsEta0_Reco->Delete();
    h_EvsEta1_Reco->Delete();
    h_EvsEta2_Reco->Delete();
    h_EvsEta3_Reco->Delete();
    h_EvsEta4_Reco->Delete();
    h_EvsEta0_Trk->Delete();
    h_EvsEta1_Trk->Delete();
    h_EvsEta2_Trk->Delete();
    h_EvsEta3_Trk->Delete();
    h_EvsEta4_Trk->Delete();
    h_EtavsE0_mcAll->Delete();
    h_EtavsE1_mcAll->Delete();
    h_EtavsE2_mcAll->Delete();
    h_EtavsE3_mcAll->Delete();
    h_EtavsE4_mcAll->Delete();
    h_EtavsE5_mcAll->Delete();
    h_EtavsE6_mcAll->Delete();
    h_EtavsE0_Reco->Delete();
    h_EtavsE1_Reco->Delete();
    h_EtavsE2_Reco->Delete();
    h_EtavsE3_Reco->Delete();
    h_EtavsE4_Reco->Delete();
    h_EtavsE5_Reco->Delete();
    h_EtavsE0_Trk->Delete();
    h_EtavsE1_Trk->Delete();
    h_EtavsE2_Trk->Delete();
    h_EtavsE3_Trk->Delete();
    h_EtavsE4_Trk->Delete();
    h_EtavsE5_Trk->Delete();
    eff_EvsEta0->Delete();
    eff_EvsEta1->Delete();
    eff_EvsEta2->Delete();
    eff_EvsEta3->Delete();
    eff_EvsEta4->Delete();
    eff_EtavsE0->Delete();
    eff_EtavsE1->Delete();
    eff_EtavsE2->Delete();
    eff_EtavsE3->Delete();
    eff_EtavsE4->Delete();
    eff_EtavsE5->Delete();
    eff_trk_EvsEta0->Delete();
    eff_trk_EvsEta1->Delete();
    eff_trk_EvsEta2->Delete();
    eff_trk_EvsEta3->Delete();
    eff_trk_EvsEta4->Delete();
    eff_trk_EtavsE0->Delete();
    eff_trk_EtavsE1->Delete();
    eff_trk_EtavsE2->Delete();
    eff_trk_EtavsE3->Delete();
    eff_trk_EtavsE4->Delete();
    eff_trk_EtavsE5->Delete();
    g_pt_pions->Delete();
    g_pt_neutrons->Delete();
  }

  // --- Close the ntuple file
  input_file->Close();
  input_file->Delete();
}
