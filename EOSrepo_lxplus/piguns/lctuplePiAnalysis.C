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

void plot_boundaries(TCanvas *canvas = nullptr) {
  TLine *line_06 = new TLine(0.6, 0., 0.6, 1.);
  TLine *line_11 = new TLine(1.1, 0., 1.1, 1.);
  TLine *line_2 = new TLine(0.8, 0., 0.8, 1.);
  TLine *line_21 = new TLine(2.1, 0., 2.1, 1.);
  line_06->SetLineColor(kRed);
  line_11->SetLineColor(kRed);
  line_2->SetLineColor(kRed);
  line_21->SetLineColor(kRed);
  line_06->SetLineStyle(1);
  line_11->SetLineStyle(1);
  line_2->SetLineStyle(1);
  line_21->SetLineStyle(1);
  line_06->SetLineWidth(3);
  line_11->SetLineWidth(3);
  line_2->SetLineWidth(3);
  line_21->SetLineWidth(3);
  line_06->Draw();
  line_11->Draw();
  line_2->Draw();
  line_21->Draw();
}

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

TH1F *h_deltaR_tsmc, *h_deltaR_tscl, *h_deltaR_tscl_pi, *h_deltaR_tscl_neu,
    *h_deltaR_tscl_others;
TH1I *h_ncl, *h_nts, *h_mcpstable;
TH1I *h_recopdg;

TH1F *h_mc_eta;
TH1F *h_mc_pt;

// ===========================================================================

void lctuplePiAnalysis(const TString filename = "output_lctuple.root") {
  const TString plotsfolder = "plots"; // plots or plotsTracks

  double etamax = 2.5;
  double etamin = 0.;
  int nbins_eta = 20;
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
    h_EvsEta0_mcAll = new TH1F(
        "h_EvsEta0_mcAll",
        "Pion events (All, 0 < |#eta| < 0.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta1_mcAll = new TH1F(
        "h_EvsEta1_mcAll",
        "Pion events (All, 0.5 < |#eta| < 1); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta2_mcAll = new TH1F(
        "h_EvsEta2_mcAll",
        "Pion events (All, 1 < |#eta| < 1.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta3_mcAll = new TH1F(
        "h_EvsEta3_mcAll",
        "Pion events (All, 1.5 < |#eta| < 2); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta4_mcAll = new TH1F(
        "h_EvsEta4_mcAll",
        "Pion events (All, 2 < |#eta| < 2.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);

    h_EvsEta0_Reco = new TH1F(
        "h_EvsEta0_Reco",
        "Pion events (Hits, 0 < |#eta| < 0.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta1_Reco = new TH1F(
        "h_EvsEta1_Reco",
        "Pion events (Hits, 0.5 < |#eta| < 1); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta2_Reco = new TH1F(
        "h_EvsEta2_Reco",
        "Pion events (Hits, 1 < |#eta| < 1.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta3_Reco = new TH1F(
        "h_EvsEta3_Reco",
        "Pion events (Hits, 1.5 < |#eta| < 2); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta4_Reco = new TH1F(
        "h_EvsEta4_Reco",
        "Pion events (Hits, 2 < |#eta| < 2.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta0_Trk = new TH1F(
        "h_EvsEta0_Trk",
        "Pion events (Hits, 0 < |#eta| < 0.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta1_Trk = new TH1F(
        "h_EvsEta1_Trk",
        "Pion events (Hits, 0.5 < |#eta| < 1); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta2_Trk = new TH1F(
        "h_EvsEta2_Trk",
        "Pion events (Hits, 1 < |#eta| < 1.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta3_Trk = new TH1F(
        "h_EvsEta3_Trk",
        "Pion events (Hits, 1.5 < |#eta| < 2); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);
    h_EvsEta4_Trk = new TH1F(
        "h_EvsEta4_Trk",
        "Pion events (Hits, 2 < |#eta| < 2.5); p_{T, gen} [GeV/c]; events",
        nlogbins_E, logbin_edges);

    h_EtavsE0_mcAll =
        new TH1F("h_EtavsE0_mcAll", "Pion events (All, p_{T}0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_mcAll =
        new TH1F("h_EtavsE1_mcAll", "Pion events (All, p_{T}1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_mcAll =
        new TH1F("h_EtavsE2_mcAll", "Pion events (All, p_{T}2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_mcAll =
        new TH1F("h_EtavsE3_mcAll", "Pion events (All, p_{T}3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_mcAll =
        new TH1F("h_EtavsE4_mcAll", "Pion events (All, p_{T}4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_mcAll =
        new TH1F("h_EtavsE5_mcAll", "Pion events (All, p_{T}5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE0_Reco =
        new TH1F("h_EtavsE0_Reco", "Pion events (Hits, p_{T}0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_Reco =
        new TH1F("h_EtavsE1_Reco", "Pion events (Hits, p_{T}1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_Reco =
        new TH1F("h_EtavsE2_Reco", "Pion events (Hits, p_{T}2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_Reco =
        new TH1F("h_EtavsE3_Reco", "Pion events (Hits, p_{T}3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_Reco =
        new TH1F("h_EtavsE4_Reco", "Pion events (Hits, p_{T}4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_Reco =
        new TH1F("h_EtavsE5_Reco", "Pion events (Hits, p_{T}5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE0_Trk =
        new TH1F("h_EtavsE0_Trk", "Pion events (Hits, p_{T}0); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_Trk =
        new TH1F("h_EtavsE1_Trk", "Pion events (Hits, p_{T}1); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_Trk =
        new TH1F("h_EtavsE2_Trk", "Pion events (Hits, p_{T}2); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_Trk =
        new TH1F("h_EtavsE3_Trk", "Pion events (Hits, p_{T}3); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_Trk =
        new TH1F("h_EtavsE4_Trk", "Pion events (Hits, p_{T}4); |#eta|; events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_Trk =
        new TH1F("h_EtavsE5_Trk", "Pion events (Hits, p_{T}5); |#eta|; events",
                 nbins_eta, etamin, etamax);
    g_pt_pions = new TGraph();
    g_pt_neutrons = new TGraph();

    h_mc_eta = new TH1F("h_mc_eta", "h_mc_eta; |#eta|; events", nbins_eta,
                        etamin, etamax);
    h_mc_pt = new TH1F("h_mc_pt", "h_mc_pt; |pt|; events", 30, 0., 300.);

    int nbins_deltaR = 40;
    /*    float min_deltaR = 1.e-6;
    float max_deltaR = 0.3;
    float deltaR_edges[31];

    for (int i = 0; i <= nbins_deltaR; i++) {
      deltaR_edges[i] =
          pow(10, TMath::Log10(min_deltaR) +
                      (TMath::Log10(max_deltaR) - TMath::Log10(min_deltaR)) /
                          double(nbins_deltaR) * double(i));
    }*/

    float max_deltaR = 0.3;

    h_deltaR_tsmc =
        new TH1F("h_deltaR_tsmc",
                 "#DeltaR between ts (IP) and MC particle; #DeltaR; events",
                 nbins_deltaR, 0., max_deltaR);
    h_deltaR_tscl =
        new TH1F("h_deltaR_tscl",
                 "#DeltaR between ts (cal) and Cluster; #DeltaR; events",
                 nbins_deltaR, 0., max_deltaR);
    h_deltaR_tscl_pi = new TH1F("h_deltaR_tscl_pi",
                                "#DeltaR between ts (cal) and Cluster - only "
                                "#pi^{#pm} reconstructed; #DeltaR; events",
                                nbins_deltaR, 0., max_deltaR);
    h_deltaR_tscl_neu = new TH1F("h_deltaR_tscl_neu",
                                 "#DeltaR between ts (cal) and Cluster - only "
                                 "n reconstructed; #DeltaR; events",
                                 nbins_deltaR, 0., max_deltaR);
    h_deltaR_tscl_others = new TH1F(
        "h_deltaR_tscl_others",
        "#DeltaR between ts (cal) and Cluster - others; #DeltaR; events",
        nbins_deltaR, 0., max_deltaR);
    h_ncl = new TH1I("h_ncl", "Number of Clusters; Number of clusters; events",
                     9, -0.5, 8.5);
    h_nts = new TH1I("h_nts", "Number of Tracks; Number of tracks; events", 5,
                     -0.5, 4.5);
    h_mcpstable = new TH1I(
        "h_mcpstable", "Number of stable MC particles; Number of prts; events",
        9, -0.5, 10.5);
    h_recopdg =
        new TH1I("h_recopdg", "Pdg of reconstructed particles; ; events", 11,
                 -0.5, 10.5);
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
  float *ts_phi = new float[100000];
  float *ts_tnl = new float[100000];
  float *ts_ome = new float[400000];
  float *ts_cov = new float[1500000];
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

  // Pandora Cluster
  int n_cl;
  float *cl_x = new float[100000];
  float *cl_y = new float[100000];
  float *cl_z = new float[100000];
  float *cl_ene = new float[100000];
  float *cl_phi = new float[100000];
  float *cl_eta = new float[100000];
  myLCTuple->SetBranchAddress("nclu", &n_cl);
  myLCTuple->SetBranchAddress("clpox", cl_x);
  myLCTuple->SetBranchAddress("clpoy", cl_y);
  myLCTuple->SetBranchAddress("clpoz", cl_z);
  myLCTuple->SetBranchAddress("clene", cl_ene);

  int n_mcp_stable = 0;

  const long int nEntries = myLCTuple->GetEntries();
  for (int ientry = 0; ientry < nEntries; ++ientry) {
    int ok_track = 0;

    myLCTuple->GetEntry(ientry);

    h_ncl->Fill(n_cl);
    h_nts->Fill(n_trst / 4);

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
        h_EvsEta0_mcAll->Fill(mcPt[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 0.5 && mcEtaAbs[n_mcp_stable] <= 1.)
        h_EvsEta1_mcAll->Fill(mcPt[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 1. && mcEtaAbs[n_mcp_stable] <= 1.5)
        h_EvsEta2_mcAll->Fill(mcPt[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 1.5 && mcEtaAbs[n_mcp_stable] <= 2.)
        h_EvsEta3_mcAll->Fill(mcPt[n_mcp_stable]);
      else if (mcEtaAbs[n_mcp_stable] >= 2. && mcEtaAbs[n_mcp_stable] <= 2.5)
        h_EvsEta4_mcAll->Fill(mcPt[n_mcp_stable]);

      if (mcPt[n_mcp_stable] >= 0. && mcPt[n_mcp_stable] < 50.)
        h_EtavsE0_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcPt[n_mcp_stable] >= 50. && mcPt[n_mcp_stable] < 100.)
        h_EtavsE1_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcPt[n_mcp_stable] >= 100. && mcPt[n_mcp_stable] < 150.)
        h_EtavsE2_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcPt[n_mcp_stable] >= 150. && mcPt[n_mcp_stable] < 200.)
        h_EtavsE3_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcPt[n_mcp_stable] >= 200. && mcPt[n_mcp_stable] <= 300.)
        h_EtavsE4_mcAll->Fill(mcEtaAbs[n_mcp_stable]);
      else if (mcPt[n_mcp_stable] >= 300.)
        h_EtavsE5_mcAll->Fill(mcEtaAbs[n_mcp_stable]);

      h_mc_pt->Fill(mcPt[n_mcp_stable]);
      h_mc_eta->Fill(mcEta[n_mcp_stable]);

      ++n_mcp_stable;
    } // i loop

    h_mcpstable->Fill(n_mcp_stable);

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

          if ( // std::fabs(reco_pt[j] - mcPt[k]) < 0.1 * mcPt[k]
              true) {
            if (TMath::Sqrt(pow((reco_eta[j] - mcEta[k]), 2) +
                            pow(deltaPhi, 2)) < 0.05) {

              if (mcEtaAbs[k] >= 0. && mcEtaAbs[k] <= 0.5)
                h_EvsEta0_Reco->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 0.5 && mcEtaAbs[k] <= 1.)
                h_EvsEta1_Reco->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 1. && mcEtaAbs[k] <= 1.5)
                h_EvsEta2_Reco->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 1.5 && mcEtaAbs[k] <= 2.)
                h_EvsEta3_Reco->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 2. && mcEtaAbs[k] <= 2.5)
                h_EvsEta4_Reco->Fill(mcPt[k]);

              if (mcPt[k] >= 0. && mcPt[k] < 50.)
                h_EtavsE0_Reco->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 50. && mcPt[k] < 100.)
                h_EtavsE1_Reco->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 100. && mcPt[k] < 150.)
                h_EtavsE2_Reco->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 150. && mcPt[k] < 200.)
                h_EtavsE3_Reco->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 200. && mcPt[k] <= 300.)
                h_EtavsE4_Reco->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 300.)
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

    // k = 0 IP, k = 1 first tracker hit, k = 2 last tracker hit, k = 3 entrance
    // ecal
    float deltaRmin_tsmc = 1000.;
    float deltaRmin_tscl = 1000.;
    float association_done = false;

    for (int t = ts_idx; t < n_trst; ++t) {

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
      if (t % 4 == 0) {

        for (int k = 0; k < n_mcp_stable; ++k) {
          if (association_done == true)
            break;

          float deltaPhi =
              std::min(std::fabs(ts_phi[t] - mcPhi[k]),
                       float(2 * M_PI - std::fabs(ts_phi[t] - mcPhi[k])));

          float deltaR =
              TMath::Sqrt(pow((ts_eta[t] - mcEta[k]), 2) + pow(deltaPhi, 2));

          if (deltaR < deltaRmin_tsmc) {
            deltaRmin_tsmc = deltaR;
          }
          if (deltaR < 0.5) {
            if (true
                //(ts_p[t] - mcPmag[k]) < 0.1 * mcPmag[k]
            ) {

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
                h_EvsEta0_Trk->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 0.5 && mcEtaAbs[k] <= 1.)
                h_EvsEta1_Trk->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 1. && mcEtaAbs[k] <= 1.5)
                h_EvsEta2_Trk->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 1.5 && mcEtaAbs[k] <= 2.)
                h_EvsEta3_Trk->Fill(mcPt[k]);
              else if (mcEtaAbs[k] >= 2. && mcEtaAbs[k] <= 2.5)
                h_EvsEta4_Trk->Fill(mcPt[k]);

              if (mcPt[k] >= 0. && mcPt[k] < 50.)
                h_EtavsE0_Trk->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 50. && mcPt[k] < 100.)
                h_EtavsE1_Trk->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 100. && mcPt[k] < 150.)
                h_EtavsE2_Trk->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 150. && mcPt[k] < 200.)
                h_EtavsE3_Trk->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 200. && mcPt[k] <= 300.)
                h_EtavsE4_Trk->Fill(mcEtaAbs[k]);
              else if (mcPt[k] >= 300.)
                h_EtavsE5_Trk->Fill(mcEtaAbs[k]);

              ok_track++;
              association_done = true;
              break;
            }
          }
        } // k loop
      }

      if (n_trst != 4)
        continue;

      // compare track and cluster
      if (t % 4 == 3 && n_cl > 0) {

        for (int cl = 0; cl < n_cl; ++cl) {

          cl_phi[cl] = TMath::ATan2(cl_y[cl], cl_x[cl]);
          cl_eta[cl] =
              0.5 *
              TMath::Log((std::sqrt(cl_x[cl] * cl_x[cl] + cl_y[cl] * cl_y[cl] +
                                    cl_z[cl] * cl_z[cl]) +
                          cl_z[cl]) /
                         (std::sqrt(cl_x[cl] * cl_x[cl] + cl_y[cl] * cl_y[cl] +
                                    cl_z[cl] * cl_z[cl]) -
                          cl_z[cl]));

          float deltaPhi =
              std::min(std::fabs(ts_phi[t] - cl_phi[cl]),
                       float(2 * M_PI - std::fabs(ts_phi[t] - cl_phi[cl])));

          float deltaR =
              TMath::Sqrt(pow((ts_eta[t] - cl_eta[cl]), 2) + pow(deltaPhi, 2));

          if (deltaR < deltaRmin_tscl)
            deltaRmin_tscl = deltaR;
        }
      }

    } // t loop

    if (n_mcp_stable < ok_track) {
      std::cout << "Event " << ientry << ": "
                << "n_mcp_stable = " << n_mcp_stable
                << ", ok_track = " << ok_track << '\n';
    }

    if (n_trst == 0)
      continue;

    if (deltaRmin_tsmc < 500.) {
      h_deltaR_tsmc->Fill(deltaRmin_tsmc);
    }

    if (deltaRmin_tscl < 500.) {
      h_deltaR_tscl->Fill(deltaRmin_tscl);

      if (n_reco == 1) {
        if (std::abs(reco_pdg[0]) == 211) {
          h_recopdg->Fill(0);
          h_deltaR_tscl_pi->Fill(deltaRmin_tscl);
        } else if (std::abs(reco_pdg[0]) == 2112) {
          h_recopdg->Fill(3);
          h_deltaR_tscl_neu->Fill(deltaRmin_tscl);
        } else if (std::abs(reco_pdg[0]) == 11) {
          h_recopdg->Fill(4);
          h_deltaR_tscl_others->Fill(deltaRmin_tscl);
        } else if (reco_pdg[0] == 22) {
          h_recopdg->Fill(5);
          h_deltaR_tscl_others->Fill(deltaRmin_tscl);
        } else {
          h_recopdg->Fill(8);
          h_deltaR_tscl_others->Fill(deltaRmin_tscl);
        }
      } else if (n_reco == 2) {
        if ((std::abs(reco_pdg[0]) == 211 && reco_pdg[1] == 22) ||
            (reco_pdg[0] == 22 && std::abs(reco_pdg[1]) == 211)) {
          h_recopdg->Fill(1);
          h_deltaR_tscl_pi->Fill(deltaRmin_tscl);
        } else if ((reco_pdg[0] == 2112 && reco_pdg[1] == 22) ||
                   (reco_pdg[0] == 22 && reco_pdg[1] == 2112)) {
          h_recopdg->Fill(6);
          h_deltaR_tscl_neu->Fill(deltaRmin_tscl);
        } else if ((std::abs(reco_pdg[0]) == 211 && reco_pdg[1] == 22) ||
                   (reco_pdg[0] == 22 && std::abs(reco_pdg[1]) == 211)) {
          h_recopdg->Fill(9);
          h_deltaR_tscl_pi->Fill(deltaRmin_tscl);
        } else if (reco_pdg[0] == 2112 && reco_pdg[1] == 2112) {
          h_recopdg->Fill(10);
          h_deltaR_tscl_neu->Fill(deltaRmin_tscl);
        } else {
          h_recopdg->Fill(8);
          h_deltaR_tscl_others->Fill(deltaRmin_tscl);
        }
      } else if (n_reco > 2) {
        if (std::abs(reco_pdg[0]) == 211) {
          h_recopdg->Fill(2);
          h_deltaR_tscl_pi->Fill(deltaRmin_tscl);
        } else if (std::abs(reco_pdg[0]) == 11) {
          h_recopdg->Fill(7);
          h_deltaR_tscl_others->Fill(deltaRmin_tscl);
        } else {
          h_recopdg->Fill(8);
          h_deltaR_tscl_neu->Fill(deltaRmin_tscl);
        }

        /*
        std::cout << '\n';
        for (int j = 0; j < n_reco; ++j) {
          std::cout << reco_pdg[j] << "  +  ";
        }
        std::cout << '\n';*/
      }
    }
  } // ientry loop

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
                            "p_{T, gen} [GeV/c]; Efficiency");
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
                            "p_{T, gen} [GeV/c]; Efficiency");
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
                            "p_{T, gen} [GeV/c]; Efficiency");
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
                            "p_{T, gen} [GeV/c]; Efficiency");
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
                            "p_{T, gen} [GeV/c]; Efficiency");
      gPad->Update();
      auto graph = eff_EvsEta4->GetPaintedGraph();
      graph->SetMinimum(0.);
      graph->SetMaximum(1.);
      gPad->Update();
      c_EvsEta4_eff->SaveAs(plotsfolder + "/c_EvsEta4_eff.png");
      c_EvsEta4_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE0_Reco, *h_EtavsE0_mcAll)) {
      auto c_EtavsE0_eff =
          new TCanvas("c_EtavsE0_eff", "c_EtavsE0_eff", 900, 800);
      eff_EtavsE0 = new TEfficiency(*h_EtavsE0_Reco, *h_EtavsE0_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE0->SetConfidenceLevel(0.68);
      eff_EtavsE0->Draw("APL");
      eff_EtavsE0->SetTitle("Pion ID efficiency (0 < p_{T, gen} < 50 GeV); "
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
      eff_EtavsE1->SetTitle("Pion ID efficiency (50 < p_{T, gen} < 100 GeV); "
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
      eff_EtavsE2->SetTitle("Pion ID efficiency (100 < p_{T, gen} < 150 GeV); "
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
      eff_EtavsE3->SetTitle("Pion ID efficiency (150 < p_{T, gen} < 200 GeV); "
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
      eff_EtavsE4->SetTitle("Pion ID efficiency (200 < p_{T, gen} < 300 GeV); "
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
      eff_EtavsE5->SetTitle("Pion ID efficiency (300 < p_{T, gen} < 500 GeV); "
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
                                "p_{T, gen} [GeV/c]; Efficiency");
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
      eff_trk_EvsEta1->SetTitle("Track efficiency (0.5 < |#eta_{gen}| < 1); "
                                "p_{T, gen} [GeV/c]; Efficiency");
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
      eff_trk_EvsEta2->SetTitle("Track efficiency (1 < |#eta_{gen}| < 1.5); "
                                "p_{T, gen} [GeV/c]; Efficiency");
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
      eff_trk_EvsEta3->SetTitle("Track efficiency (1.5 < |#eta_{gen}| < 2); "
                                "p_{T, gen} [GeV/c]; Efficiency");
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
      eff_trk_EvsEta4->SetTitle("Track efficiency (2 < |#eta_{gen}| < 2.5); "
                                "p_{T, gen} [GeV/c]; Efficiency");
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
      eff_trk_EtavsE0->SetTitle("Track efficiency (0 < p_{T, gen} < 50 GeV); "
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
      eff_trk_EtavsE1->SetTitle("Track efficiency (50 < p_{T, gen} < 100 GeV); "
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
    /*
    for (int b = 0; b != h_EtavsE1_Trk->GetNbinsX() + 2; ++b) {
      std::cout << "h_EtavsE1_Trk (bin " << b
                << ") = " << h_EtavsE0_Reco->GetBinContent(b) << '\n';
      std::cout << "h_EtavsE1_mcAll (bin " << b
                << ") = " << h_EtavsE1_mcAll->GetBinContent(b) << '\n';
    }*/

    if (TEfficiency::CheckConsistency(*h_EtavsE2_Trk, *h_EtavsE2_mcAll)) {
      auto c_EtavsE2_eff_trk =
          new TCanvas("c_EtavsE2_eff_trk", "c_EtavsE2_eff_trk", 900, 800);
      eff_trk_EtavsE2 = new TEfficiency(*h_EtavsE2_Trk, *h_EtavsE2_mcAll);
      gStyle->SetOptStat(111);
      // eff_trk_EtavsE2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_EtavsE2->SetConfidenceLevel(0.68);
      eff_trk_EtavsE2->Draw("APL");
      eff_trk_EtavsE2->SetTitle(
          "Track efficiency (100 < p_{T, gen} < 150 GeV); "
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
      eff_trk_EtavsE3->SetTitle(
          "Track efficiency (150 < p_{T, gen} < 200 GeV); "
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
      eff_trk_EtavsE4->SetTitle(
          "Track efficiency (200 < p_{T, gen} < 300 GeV); "
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
      eff_trk_EtavsE5->SetTitle(
          "Track efficiency (300 < p_{T, gen} < 500 GeV); "
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
    g_pt_pions->GetXaxis()->SetLimits(0., 350.);
    g_pt_pions->SetMaximum(0.5);
    g_pt_pions->Draw("AP");
    c_ptvssigma->cd(2);
    gStyle->SetOptStat(111111);
    g_pt_neutrons->SetTitle(
        "p_{T} tracks (no pions);p_{T} [GeV];#sigma_{p_{T}}/^{}p_{T}");
    g_pt_neutrons->GetYaxis()->SetTitleOffset(0.7);
    g_pt_neutrons->SetMarkerStyle(26);
    g_pt_neutrons->SetMarkerSize(0.3);
    g_pt_neutrons->SetMarkerColor(kRed);
    gPad->SetLogy();
    g_pt_neutrons->GetXaxis()->SetLimits(0., 350.);
    g_pt_neutrons->SetMaximum(0.5);
    g_pt_neutrons->Draw("AP");
    c_ptvssigma->SaveAs(plotsfolder + "/c_scatter.pdf");
    c_ptvssigma->Close();

    TCanvas *c_deltaR_tsmc =
        new TCanvas("c_deltaR_tsmc", "c_deltaR_tsmc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    // c_deltaR_tsmc->SetLogx();
    h_deltaR_tsmc->Draw();
    c_deltaR_tsmc->SaveAs("./plots/c_deltaR_tsmc.png");
    c_deltaR_tsmc->Close();

    TCanvas *c_deltaR_tscl =
        new TCanvas("c_deltaR_tscl", "c_deltaR_tscl", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    // c_deltaR_tscl->SetLogx();
    h_deltaR_tscl->GetXaxis()->SetRangeUser(1.e-4, 0.3);
    h_deltaR_tscl->Draw();
    c_deltaR_tscl->SaveAs("./plots/c_deltaR_tscl.png");
    c_deltaR_tscl->Close();

    TCanvas *c_deltaR_tscl_pi =
        new TCanvas("c_deltaR_tscl_pi", "c_deltaR_tscl_pi", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    // c_deltaR_tscl_pi->SetLogx();
    h_deltaR_tscl_pi->GetXaxis()->SetRangeUser(1.e-4, 0.3);
    h_deltaR_tscl_pi->Draw();
    c_deltaR_tscl_pi->SaveAs("./plots/c_deltaR_tscl_pi.png");
    c_deltaR_tscl_pi->Close();

    TCanvas *c_deltaR_tscl_neu =
        new TCanvas("c_deltaR_tscl_neu", "c_deltaR_tscl_neu", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    // c_deltaR_tscl_neu->SetLogx();
    h_deltaR_tscl_neu->GetXaxis()->SetRangeUser(1.e-4, 0.3);
    h_deltaR_tscl_neu->Draw();
    c_deltaR_tscl_neu->SaveAs("./plots/c_deltaR_tscl_neu.png");
    c_deltaR_tscl_neu->Close();

    TCanvas *c_deltaR_tscl_others =
        new TCanvas("c_deltaR_tscl_others", "c_deltaR_tscl_others", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    // c_deltaR_tscl_others->SetLogx();
    h_deltaR_tscl_others->GetXaxis()->SetRangeUser(1.e-4, 0.3);
    h_deltaR_tscl_others->Draw();
    c_deltaR_tscl_others->SaveAs("./plots/c_deltaR_tscl_others.png");
    c_deltaR_tscl_others->Close();

    TCanvas *c_ncl = new TCanvas("c_ncl", "c_ncl", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_ncl->Draw();
    c_ncl->SaveAs("./plots/c_ncl.png");
    c_ncl->Close();

    TCanvas *c_nts = new TCanvas("c_nts", "c_nts", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_nts->Draw();
    c_nts->SaveAs("./plots/c_nts.png");
    c_nts->Close();
    TCanvas *c_mcpstable = new TCanvas("c_mcpstable", "c_mcpstable", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_mcpstable->Draw();
    c_mcpstable->SaveAs("./plots/c_mcpstable.png");
    c_mcpstable->Close();

    const char *reco_labels[11] = {"#pi^{#pm}",
                                   "#pi^{#pm} + #gamma",
                                   "#pi^{#pm} + >1 prts ",
                                   "n",
                                   "e^{#pm}",
                                   "#gamma",
                                   "n + #gamma",
                                   "e^{#pm} + #gamma",
                                   "Others (ns+#gammas)",
                                   "#pi^{#pm} + n",
                                   "2n"};

    TCanvas *c_recopdg = new TCanvas("c_recopdg", "c_recopdg", 800, 800);
    gStyle->SetOptStat(111111);
    h_recopdg->Draw();
    for (int i = 1; i != h_recopdg->GetNbinsX() + 1; i++)
      h_recopdg->GetXaxis()->SetBinLabel(i, reco_labels[i - 1]);
    c_recopdg->SaveAs("./plots/c_recopdg.png");
    c_recopdg->Close();

    TCanvas *c_mc_pt = new TCanvas("c_mc_pt", "c_mc_pt", 800, 800);
    gStyle->SetOptStat(111111);
    h_mc_pt->Draw();
    c_mc_pt->SaveAs("./plots/c_mc_pt.png");
    c_mc_pt->Close();

    TCanvas *c_mc_eta = new TCanvas("c_mc_eta", "c_mc_eta", 800, 800);
    gStyle->SetOptStat(111111);
    h_mc_eta->Draw();
    c_mc_eta->SaveAs("./plots/c_mc_eta.png");
    c_mc_eta->Close();
  }

  // Create a new .root file
  TFile *outputFile = new TFile("../efficiencies.root", "RECREATE");

  // Create a directory (folder) named "piguns" inside the root file
  TDirectory *dir = outputFile->mkdir("piguns");
  dir->cd();

  // Save the TEfficiency objects to the file
  eff_EvsEta0->Write("eff_EvsEta0");
  eff_EvsEta1->Write("eff_EvsEta1");
  eff_EvsEta2->Write("eff_EvsEta2");
  eff_EvsEta3->Write("eff_EvsEta3");
  eff_EvsEta4->Write("eff_EvsEta4");
  eff_EtavsE0->Write("eff_EtavsE0");
  eff_EtavsE1->Write("eff_EtavsE1");
  eff_EtavsE2->Write("eff_EtavsE2");
  eff_EtavsE3->Write("eff_EtavsE3");
  eff_EtavsE4->Write("eff_EtavsE4");
  eff_EtavsE5->Write("eff_EtavsE5");
  // Close the file
  outputFile->Close();
  // Clean up
  outputFile->Delete();

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

    delete[] cl_ene;
    delete[] cl_phi;
    delete[] cl_x;
    delete[] cl_y;
    delete[] cl_z;
    delete[] cl_eta;

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
    h_deltaR_tsmc->Delete();
    h_deltaR_tscl->Delete();
    h_deltaR_tscl_pi->Delete();
    h_deltaR_tscl_neu->Delete();
    h_deltaR_tscl_others->Delete();
    h_ncl->Delete();
    h_nts->Delete();
    h_recopdg->Delete();
  }

  // --- Close the ntuple file
  input_file->Close();
  input_file->Delete();
}
