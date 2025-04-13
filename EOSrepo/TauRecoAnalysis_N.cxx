//analysis and plots of the reconstructed taus, with neutral particles reconstructed too

#include <fstream>
#include <iostream>

#include "Math/Vector4D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// histogram declaration

TH1I *h_mcDecayMode, *h_mymcDecayMode, *h_recDecayMode, *h_myrecDecayMode;

TH1F *h_EvsEta0_mcAll, *h_EvsEta1_mcAll, *h_EvsEta2_mcAll, *h_EvsEta3_mcAll,
    *h_EvsEta4_mcAll, *h_EvsEta0_mcHits, *h_EvsEta1_mcHits, *h_EvsEta2_mcHits,
    *h_EvsEta3_mcHits, *h_EvsEta4_mcHits;

TH1F *h_EtavsE0_mcAll, *h_EtavsE1_mcAll, *h_EtavsE2_mcAll, *h_EtavsE3_mcAll,
    *h_EtavsE4_mcAll, *h_EtavsE5_mcAll, *h_EtavsE6_mcAll, *h_EtavsE0_mcHits,
    *h_EtavsE1_mcHits, *h_EtavsE2_mcHits, *h_EtavsE3_mcHits, *h_EtavsE4_mcHits,
    *h_EtavsE5_mcHits, *h_EtavsE6_mcHits;

TH2F *h_E2d_Eta0, *h_E2d_Eta1, *h_E2d_Eta2, *h_E2d_Eta3, *h_E2d_Eta4,
    *h_EvsEta_mcAll, *h_EvsEta_mcHits, *h_EvsEvis;
TH1F *h_ECorr_Eta0, *h_ECorr_Eta1, *h_ECorr_Eta2, *h_ECorr_Eta3, *h_ECorr_Eta4;
TProfile *EProfile_Eta0, *EProfile_Eta1, *EProfile_Eta2, *EProfile_Eta3,
    *EProfile_Eta4;

TEfficiency *eff_EvsEta, *eff_EvsEta0, *eff_EvsEta1, *eff_EvsEta2, *eff_EvsEta3,
    *eff_EvsEta4, *eff_EtavsE0, *eff_EtavsE1, *eff_EtavsE2, *eff_EtavsE3,
    *eff_EtavsE4, *eff_EtavsE5, *eff_EtavsE6;

/*void plot_boundaries(const TString &decay = "all", TCanvas *canvas = nullptr)
{

  // if (decay == "1pi" || decay == "1pi0N" || decay == "1pi1N" || decay ==
  // "1pi2N" || decay == "1pi3N" || decay == "3pi" || decay == "had"|| decay ==
  // "other") {}
  TLine *end_outertrackerbarrel = new TLine(0.763, 0., 0.763, 1.);
  TLine *begin_outertrackerendcap = new TLine(1.217, 0., 1.217, 1.);
  TLine *end_outertrackerendcap = new TLine(2.0, 0., 2.0, 1.);
  TLine *end_innertrackerbarrel = new TLine(1.062, 0., 1.062, 1.);
  TLine *begin_innertrackerendcap = new TLine(2.134, 0., 2.134, 1.);
  TLine *begin_vertexendcap = new TLine(1.549, 0., 1.549, 1.);
  // TLine *end_vertexbarrel = new TLine();
  TLine *nozzle = new TLine(2.436, 0., 2.436, 1.);
  end_outertrackerbarrel->SetLineColor(kCyan);
  begin_outertrackerendcap->SetLineColor(kCyan);
  end_outertrackerendcap->SetLineColor(kCyan);
  end_innertrackerbarrel->SetLineColor(kGreen);
  begin_innertrackerendcap->SetLineColor(kGreen);
  begin_vertexendcap->SetLineColor(kViolet);
  nozzle->SetLineColor(kRed);
  end_outertrackerbarrel->SetLineStyle(9);
  begin_outertrackerendcap->SetLineStyle(2);
  end_outertrackerendcap->SetLineStyle(3);
  end_innertrackerbarrel->SetLineStyle(9);
  begin_innertrackerendcap->SetLineStyle(2);
  begin_vertexendcap->SetLineStyle(2);
  nozzle->SetLineStyle(1);
  end_outertrackerbarrel->SetLineWidth(3);
  begin_outertrackerendcap->SetLineWidth(3);
  end_outertrackerendcap->SetLineWidth(3);
  end_innertrackerbarrel->SetLineWidth(3);
  begin_innertrackerendcap->SetLineWidth(3);
  begin_vertexendcap->SetLineWidth(3);
  nozzle->SetLineWidth(3);
  end_outertrackerbarrel->Draw();
  begin_outertrackerendcap->Draw();
  end_outertrackerendcap->Draw();
  end_innertrackerbarrel->Draw();
  //begin_innertrackerendcap->Draw();
  begin_vertexendcap->Draw();
  nozzle->Draw();
}*/

void plot_boundaries(const TString &decay = "1piN", TCanvas *canvas = nullptr) {
  /*
  // if (decay == "1pi" || decay == "1pi0N" || decay == "1pi1N" || decay ==
  // "1pi2N" || decay == "1pi3N" || decay == "3pi" || decay == "had"|| decay ==
  // "other") {}
  TLine *line_06 = new TLine(0.6, 0., 0.6, 1.);
  TLine *line_11 = new TLine(1.1, 0., 1.1, 1.);
  TLine *line_2 = new TLine(2., 0., 2., 1.);
  line_06->SetLineColor(kRed);
  line_11->SetLineColor(kRed);
  line_2->SetLineColor(kRed);
  line_06->SetLineStyle(1);
  line_11->SetLineStyle(1);
  line_2->SetLineStyle(1);
  line_06->SetLineWidth(3);
  line_11->SetLineWidth(3);
  line_2->SetLineWidth(3);
  line_06->Draw();
  line_11->Draw();
  line_2->Draw();*/
}

void Histo_filler(const TString &decay = "N",
                  const TString &filename = "output_EvalTauFinder.root") {

  float etamax = 2.5;
  float etamin = 0.;
  int nbins_eta = 20;
  double Emax = 1000.;
  double Emin = 5.;
  int nbins_E = 50;
  int nlogbins_E = 13;
  float binsemiwidth_E = 0.5 * (Emax - Emin) / nbins_E;
  /*double logbin_edges[nlogbins_E + 1];

  for (int i = 0; i <= nlogbins_E; i++) {
    logbin_edges_E[i] =
        pow(10, TMath::Log10(Emin) + (TMath::Log10(Emax) - TMath::Log10(Emin)) /
                                         double(nlogbins_E) * double(i));
    // std::cout << logbin_edges[i] << '\n';
  }*/

  /*double logbin_edges_E[18] = {5.,   10.,  20.,  30.,  40.,  60.,
                               80.,  100., 120., 140., 160., 180.,
                               200., 250., 300., 500., 800., 1000.};*/
  double logbin_edges_E[14] = {5.,   10.,  20.,  30.,  40.,  60.,  80.,
                               100., 120., 150., 200., 300., 500., 800.};

  // histogram initialization
  {
    h_mcDecayMode = new TH1I("h_mcDecayMode", "MC #tau decay mode; Decay Mode",
                             9, -0.5, 8.5);
    h_mymcDecayMode = new TH1I("h_mymcDecayMode",
                               "MC #tau decay mode; Decay Mode", 5, -0.5, 4.5);

    h_recDecayMode = new TH1I("h_recDecayMode",
                              "Rec #tau decay mode; Decay Mode", 8, -0.5, 7.5);

    h_myrecDecayMode = new TH1I(
        "h_myrecDecayMode", "Rec #tau decay mode; Decay Mode", 5, -0.5, 4.5);
    h_E2d_Eta0 =
        new TH2F("h_E2d_Eta0",
                 "Tau reco energy vs Visible gen energy (0 < |#eta| < 0.5); "
                 "p_{T, reco}^{vis} [GeV/c]; p_{T, gen}^{vis} [GeV/c]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta1 =
        new TH2F("h_E2d_Eta1",
                 "Tau reco energy vs Visible gen energy (0.5 < |#eta| < 1); "
                 "p_{T, reco}^{vis} [GeV/c]; p_{T, gen}^{vis} [GeV/c]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta2 =
        new TH2F("h_E2d_Eta2",
                 "Tau reco energy vs Visible gen energy (1. < |#eta| < 1.5); "
                 "p_{T, reco}^{vis} [GeV/c]; p_{T, gen}^{vis} [GeV/c]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta3 =
        new TH2F("h_E2d_Eta3",
                 "Tau reco energy vs Visible gen energy (1.5 < |#eta| < 2); "
                 "p_{T, reco}^{vis} [GeV/c]; p_{T, gen}^{vis} [GeV/c]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta4 =
        new TH2F("h_E2d_Eta4",
                 "Tau reco energy vs Visible gen energy (2 < |#eta|); "
                 "p_{T, reco}^{vis} [GeV/c]; p_{T, gen}^{vis} [GeV/c]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);

    h_EvsEta0_mcAll = new TH1F("h_EvsEta0_mcAll",
                               "Tau E events (All, 0 < |#eta| < 0.5); p_{T, "
                               "gen}^{vis} [GeV/c]; Events",
                               nlogbins_E, logbin_edges_E);
    h_EvsEta1_mcAll = new TH1F(
        "h_EvsEta1_mcAll",
        "Tau events (All, 0.5 < |#eta| < 1); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta2_mcAll = new TH1F(
        "h_EvsEta2_mcAll",
        "Tau events (All, 1 < |#eta| < 1.5); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta3_mcAll = new TH1F(
        "h_EvsEta3_mcAll",
        "Tau events (All, 1.5 < |#eta| < 2); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta4_mcAll = new TH1F(
        "h_EvsEta4_mcAll",
        "Tau events (All, |#eta| > 2); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);

    h_EvsEta0_mcHits = new TH1F(
        "h_EvsEta0_mcHits",
        "Tau events (Hits, 0 < |#eta| < 0.5); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta1_mcHits = new TH1F(
        "h_EvsEta1_mcHits",
        "Tau events (Hits, 0.5 < |#eta| < 1); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta2_mcHits = new TH1F(
        "h_EvsEta2_mcHits",
        "Tau events (Hits, 1 < |#eta| < 1.5); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta3_mcHits = new TH1F(
        "h_EvsEta3_mcHits",
        "Tau events (Hits, 1.5 < |#eta| < 2); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta4_mcHits = new TH1F(
        "h_EvsEta4_mcHits",
        "Tau events (Hits, |#eta| > 2); p_{T, gen}^{vis} [GeV/c]; Events",
        nlogbins_E, logbin_edges_E);

    h_EtavsE0_mcAll =
        new TH1F("h_EtavsE0_mcAll", "Tau events (All, E0); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_mcAll =
        new TH1F("h_EtavsE1_mcAll", "Tau events (All, E1); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_mcAll =
        new TH1F("h_EtavsE2_mcAll", "Tau events (All, E2); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_mcAll =
        new TH1F("h_EtavsE3_mcAll", "Tau events (All, E3); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_mcAll =
        new TH1F("h_EtavsE4_mcAll", "Tau events (All, E4); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_mcAll =
        new TH1F("h_EtavsE5_mcAll", "Tau events (All, E5); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE6_mcAll =
        new TH1F("h_EtavsE6_mcAll", "Tau events (All, E6); |#eta|; Events",
                 nbins_eta, etamin, etamax);

    h_EtavsE0_mcHits =
        new TH1F("h_EtavsE0_mcHits", "Tau events (Hits, E0); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE1_mcHits =
        new TH1F("h_EtavsE1_mcHits", "Tau events (Hits, E1); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE2_mcHits =
        new TH1F("h_EtavsE2_mcHits", "Tau events (Hits, E2); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE3_mcHits =
        new TH1F("h_EtavsE3_mcHits", "Tau events (Hits, E3); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE4_mcHits =
        new TH1F("h_EtavsE4_mcHits", "Tau events (Hits, E4); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE5_mcHits =
        new TH1F("h_EtavsE5_mcHits", "Tau events (Hits, E5); |#eta|; Events",
                 nbins_eta, etamin, etamax);
    h_EtavsE6_mcHits =
        new TH1F("h_EtavsE6_mcHits", "Tau events (Hits, E6); |#eta|; Events",
                 nbins_eta, etamin, etamax);

    h_EvsEta_mcAll =
        new TH2F("h_EvsEta_mcAll",
                 "Tau visible E vs #eta (MC, all particles); "
                 "p_{T, gen}^{vis} [GeV/c]; |#eta_{gen}^{vis}|",
                 nlogbins_E, logbin_edges_E, nbins_eta, etamin, etamax);
    h_EvsEta_mcHits =
        new TH2F("h_EvsEta_mcHits",
                 "Tau visible E vs #eta (MC, reco particles); "
                 "p_{T, gen}^{vis} [GeV/c]; |#eta_{gen}^{vis}|",
                 nlogbins_E, logbin_edges_E, nbins_eta, etamin, etamax);

    h_EvsEvis =
        new TH2F("h_EvsEvis",
                 "#tau E generated vs E visible (MC); "
                 "E_{gen}^{#tau} [GeV]; p_{T, gen}^{vis} [GeV/c];",
                 nbins_E + 1, Emin - binsemiwidth_E, Emax + binsemiwidth_E,
                 nbins_E + 1, Emin - binsemiwidth_E, Emax + binsemiwidth_E);

    h_ECorr_Eta0 = new TH1F("h_ECorr_Eta0",
                            "energy corrections (0 < |#eta| < 0.5); "
                            "p_{T, reco}^{vis}(#tau) [GeV/c]; #frac{p_{T, "
                            "gen}^{vis}}{p_{T, reco}^{vis}}",
                            nlogbins_E, logbin_edges_E);
    h_ECorr_Eta1 = new TH1F("h_ECorr_Eta1",
                            "energy corrections (0.5 < |#eta| < 1); "
                            "p_{T, reco}^{vis}(#tau) [GeV/c]; #frac{p_{T, "
                            "gen}^{vis}}{p_{T, reco}^{vis}}",
                            nlogbins_E, logbin_edges_E);
    h_ECorr_Eta2 = new TH1F("h_ECorr_Eta2",
                            "energy corrections (1 < |#eta| < 1.5); "
                            "p_{T, reco}^{vis}(#tau) [GeV/c]; #frac{p_{T, "
                            "gen}^{vis}}{p_{T, reco}^{vis}}",
                            nlogbins_E, logbin_edges_E);
    h_ECorr_Eta3 = new TH1F("h_ECorr_Eta3",
                            "energy corrections (1.5 < |#eta| < 2.); "
                            "p_{T, reco}^{vis}(#tau) [GeV/c]; #frac{p_{T, "
                            "gen}^{vis}}{p_{T, reco}^{vis}}",
                            nlogbins_E, logbin_edges_E);
    h_ECorr_Eta4 = new TH1F("h_ECorr_Eta4",
                            "energy corrections (2. < |#eta|); "
                            "p_{T, reco}^{vis}(#tau) [GeV/c]; #frac{p_{T, "
                            "gen}^{vis}}{p_{T, reco}^{vis}}",
                            nlogbins_E, logbin_edges_E);
  }

  //  --- Open the ntuple file and get the tree
  TFile *input_file = new TFile(filename, "READ");

  if (!input_file->IsOpen())
    throw std::invalid_argument("filename not valid");

  TTree *tree = (TTree *)input_file->Get("evtree;1");

  // --- MC particles
  int ntau_mc;
  int ntau_rec;

  int evID;
  int runID;

  // int *evt_id = new int[10000];
  float *mcPtvis = new float[10000000];
  float *recE = new float[10000000];
  float *mcPt = new float[10000000];
  float *mcEvis = new float[10000000];
  float *recPt = new float[10000000];
  float *mcEta = new float[10000000];
  float *recEta = new float[10000000];
  float *mcPhi = new float[10000000];
  float *recPhi = new float[10000000];
  int *mcDecayMode = new int[10000000];
  int *tau_charge = new int[10000000];
  int *recNPfos = new int[10000000];
  int *recNQTracks = new int[10000000];
  int *pfosPdg = new int[50000000];
  float *pfosPt = new float[50000000];
  float *pfosDeltaR = new float[50000000];

  tree->SetBranchAddress("nTausMC", &ntau_mc);
  tree->SetBranchAddress("nTausRec", &ntau_rec);
  tree->SetBranchAddress("EvID", &evID);
  tree->SetBranchAddress("RunID", &runID);
  tree->SetBranchAddress("mcDecayMode", mcDecayMode);
  tree->SetBranchAddress("mcPt", mcPt);
  tree->SetBranchAddress("mcPt_vis", mcPtvis);
  tree->SetBranchAddress("mcE_vis", mcEvis);
  tree->SetBranchAddress("recE", recE);
  tree->SetBranchAddress("recPt", recPt);
  tree->SetBranchAddress("mcEta_vis", mcEta);
  tree->SetBranchAddress("recEta", recEta);
  tree->SetBranchAddress("mcPhi_vis", mcPhi);
  tree->SetBranchAddress("recPhi", recPhi);
  tree->SetBranchAddress("recNPfos", recNPfos);
  tree->SetBranchAddress("recNQTracks", recNQTracks);
  tree->SetBranchAddress("charge", tau_charge);
  tree->SetBranchAddress("pfosPdg", pfosPdg);
  tree->SetBranchAddress("pfosPt", pfosPt);
  tree->SetBranchAddress("pfosDeltaR", pfosDeltaR);

  const long int nEntries = tree->GetEntries();

  for (int ientry = 0; ientry < nEntries; ++ientry) {

    tree->GetEntry(ientry);

    int mymcDecayMode[ntau_mc];

    for (unsigned int i = 0; i < ntau_mc; ++i) {

      if (mcDecayMode[i] < 1 || mcDecayMode[i] > 3)
        continue;

      mcEta[i] = std::fabs(mcEta[i]);

      // fix eta, draw histos for E, p
      if (mcEta[i] >= 0. && mcEta[i] <= 0.5)
        h_EvsEta0_mcAll->Fill(mcPtvis[i]);
      else if (mcEta[i] >= 0.5 && mcEta[i] <= 1.)
        h_EvsEta1_mcAll->Fill(mcPtvis[i]);
      else if (mcEta[i] >= 1. && mcEta[i] <= 1.5)
        h_EvsEta2_mcAll->Fill(mcPtvis[i]);
      else if (mcEta[i] >= 1.5 && mcEta[i] <= 2.)
        h_EvsEta3_mcAll->Fill(mcPtvis[i]);
      else if (mcEta[i] >= 2.)
        h_EvsEta4_mcAll->Fill(mcPtvis[i]);

      // fix E (compute efficiencies)
      if (mcPtvis[i] >= 0. && mcPtvis[i] < 50.)
        h_EtavsE0_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 50. && mcPtvis[i] < 100.)
        h_EtavsE1_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 100. && mcPtvis[i] < 150.)
        h_EtavsE2_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 150. && mcPtvis[i] < 200.)
        h_EtavsE3_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 200. && mcPtvis[i] < 300.)
        h_EtavsE4_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 300. && mcPtvis[i] < 500.)
        h_EtavsE5_mcAll->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 500.)
        h_EtavsE6_mcAll->Fill(mcEta[i]);

      h_EvsEta_mcAll->Fill(mcPtvis[i], mcEta[i]);
      h_EvsEvis->Fill(mcPt[i], mcPtvis[i]);

      h_mcDecayMode->Fill(mcDecayMode[i]);
      h_mymcDecayMode->Fill(mymcDecayMode[i]);
    }
  }
  for (int ientry = 0; ientry < nEntries; ++ientry) {

    tree->GetEntry(ientry);

    if (ntau_mc != 1 || ntau_rec != 1)
      continue;

    if (abs(tau_charge[0]) != 1)
      continue;

    float deltaPhi =
        std::min(std::fabs(recPhi[0] - mcPhi[0]),
                 float(2 * M_PI - std::fabs(recPhi[0] - mcPhi[0])));

    if (TMath::Sqrt(pow((recEta[0] - mcEta[0]), 2) + pow(deltaPhi, 2)) >= 0.05)
      continue;

    /*if (std::fabs(recPt[0] - mcPtvis[0]) > 0.1 * mcPtvis[0])
      continue;*/

    int recDecayMode = -1;
    int myrecDecayMode = -1;
    int mymcDecayMode[ntau_mc];

    int n_el = 0;
    int n_mu = 0;
    int n_pipl = 0;
    int n_pimi = 0;
    int n_gamma = 0;
    int n_neutrals = recNPfos[0] - recNQTracks[0];

    for (int j = 0; j != recNPfos[0]; j++) {
      switch (abs(pfosPdg[j])) {
      case 11:
        n_el++;
        break;
      case 13:
        n_mu++;
        break;
      case 22:
        n_gamma++;
        break;
      case 211:
      case 321:
        n_pipl++;
        break;
      case -211:
      case -321:
        n_pimi++;
        break;
      }
    }

    int n_pi = n_pipl + n_pimi;

    if (n_pi == 1 && recNQTracks[0] == 1) {
      recDecayMode = n_neutrals;
    } else if (n_pi == 3 && recNQTracks[0] == 3) {
      recDecayMode = 4;
    } else if (n_el == 1 && n_pi == 0) {
      recDecayMode = 5;
    } else if (n_mu == 1 && n_pi == 0) {
      recDecayMode = 6;
    } else {
      recDecayMode = 7;
    }

    h_recDecayMode->Fill(recDecayMode);
    h_myrecDecayMode->Fill(myrecDecayMode);

    if (recDecayMode < 1 || recDecayMode > 3)
      continue;

    for (unsigned int i = 0; i < ntau_mc; ++i) {

      if (mcDecayMode[i] < 1 || mcDecayMode[i] > 3)
        continue;

      mcEta[i] = std::fabs(mcEta[i]);

      // fix eta, draw histos for E, p
      if (mcEta[i] >= 0. && mcEta[i] <= 0.5) {
        h_E2d_Eta0->Fill(recPt[i], mcPtvis[i]);
        h_EvsEta0_mcHits->Fill(mcPtvis[i]);
      } else if (mcEta[i] >= 0.5 && mcEta[i] <= 1.) {
        h_E2d_Eta1->Fill(recPt[i], mcPtvis[i]);

        h_EvsEta1_mcHits->Fill(mcPtvis[i]);
      } else if (mcEta[i] >= 1. && mcEta[i] <= 1.5) {
        h_E2d_Eta2->Fill(recPt[i], mcPtvis[i]);

        h_EvsEta2_mcHits->Fill(mcPtvis[i]);
      } else if (mcEta[i] >= 1.5 && mcEta[i] <= 2.) {
        h_E2d_Eta3->Fill(recPt[i], mcPtvis[i]);

        h_EvsEta3_mcHits->Fill(mcPtvis[i]);
      } else if (mcEta[i] >= 2.) {
        h_E2d_Eta4->Fill(recPt[i], mcPtvis[i]);

        h_EvsEta4_mcHits->Fill(mcPtvis[i]);
      }

      // fix E (compute efficiencies)
      if (mcPtvis[i] >= 0. && mcPtvis[i] < 50.)
        h_EtavsE0_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 50. && mcPtvis[i] < 100.)
        h_EtavsE1_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 100. && mcPtvis[i] < 150.)
        h_EtavsE2_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 150. && mcPtvis[i] < 200.)
        h_EtavsE3_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 200. && mcPtvis[i] < 300.)
        h_EtavsE4_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 300. && mcPtvis[i] < 500.)
        h_EtavsE5_mcHits->Fill(mcEta[i]);
      else if (mcPtvis[i] >= 500.)
        h_EtavsE6_mcHits->Fill(mcEta[i]);

      h_EvsEta_mcHits->Fill(mcPtvis[i], mcEta[i]);
    }
  }
  tree->Delete();
  input_file->Close();
  input_file->Delete();

  const char *dm_mclabels[9] = {"#pi^{#pm}",
                                "#pi^{#pm} + #pi^{0}",
                                "#pi^{#pm} + 2#pi^{0}",
                                "#pi^{#pm} + 3#pi^{0}",
                                "3#pi^{#pm}",
                                "3#pi^{#pm} + #pi^{0}",
                                "e^{#pm}",
                                "#mu^{#pm}",
                                "Other"};

  const char *dm_reclabels[9] = {
      "#pi^{#pm}", "#pi^{#pm} + 1N", "#pi^{#pm} + 2N", "#pi^{#pm} + 3N",
      "3-prong",   "e^{#pm}",        "#mu^{#pm}",      "Other"};

  const char *dm_mylabels[5] = {"1-prong", "3-prong", "e^{#pm}", "#mu^{#pm}",
                                "Other"};

  auto c_recDecayMode =
      new TCanvas("c_recDecayMode", "c_recDecayMode", 800, 800);
  // gStyle->SetOptStat(111111);
  // h_recDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
  h_recDecayMode->Draw("HISTO");
  // h_recDecayMode->Scale(1./h_recDecayMode->GetEntries());
  gPad->Update();
  for (int i = 1; i <= 8; i++)
    h_recDecayMode->GetXaxis()->SetBinLabel(i, dm_reclabels[i - 1]);
  h_recDecayMode->GetXaxis()->SetLabelSize(0.045);
  h_recDecayMode->GetXaxis()->SetTitleOffset(1.4);
  h_recDecayMode->GetXaxis()->CenterTitle(true);
  c_recDecayMode->SaveAs("plots_reco_" + decay + "/c_recDecayMode.png");
  c_recDecayMode->Close();

  auto c_myrecDecayMode =
      new TCanvas("c_myrecDecayMode", "c_myrecDecayMode", 800, 800);
  h_myrecDecayMode->Draw("HISTO");
  // h_myrecDecayMode->Scale(1. / h_myrecDecayMode->GetEntries());
  //  h_myrecDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
  for (int i = 1; i <= 5; i++)
    h_myrecDecayMode->GetXaxis()->SetBinLabel(i, dm_mylabels[i - 1]);
  h_myrecDecayMode->GetXaxis()->SetLabelSize(0.045);
  h_myrecDecayMode->GetXaxis()->SetTitleOffset(1.4);
  h_myrecDecayMode->GetXaxis()->CenterTitle(true);
  c_myrecDecayMode->SaveAs("plots_reco_" + decay + "/c_myrecDecayMode.png");
  c_myrecDecayMode->Close();

  auto c_mcDecayMode = new TCanvas("c_mcDecayMode", "c_mcDecayMode", 800, 800);
  h_mcDecayMode->Draw("HISTO");
  // h_mcDecayMode->Scale(1. / h_mcDecayMode->GetEntries());
  //  h_mcDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
  for (int i = 1; i <= 9; i++)
    h_mcDecayMode->GetXaxis()->SetBinLabel(i, dm_mclabels[i - 1]);
  h_mcDecayMode->GetXaxis()->SetLabelSize(0.045);
  h_mcDecayMode->GetXaxis()->SetTitleOffset(1.4);
  h_mcDecayMode->GetXaxis()->CenterTitle(true);
  c_mcDecayMode->SaveAs("plots_reco_" + decay + "/c_mcDecayMode.png");
  c_mcDecayMode->Close();

  auto c_mymcDecayMode =
      new TCanvas("c_mymcDecayMode", "c_mymcDecayMode", 800, 800);
  h_mymcDecayMode->Draw("HISTO");
  // h_mymcDecayMode->Scale(1. / h_mymcDecayMode->GetEntries());
  //  h_mymcDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
  for (int i = 1; i <= 5; i++)
    h_mymcDecayMode->GetXaxis()->SetBinLabel(i, dm_mylabels[i - 1]);
  h_mymcDecayMode->GetXaxis()->SetLabelSize(0.045);
  h_mymcDecayMode->GetXaxis()->SetTitleOffset(1.4);
  h_mymcDecayMode->GetXaxis()->CenterTitle(true);
  c_mymcDecayMode->SaveAs("plots_reco_" + decay + "/c_mymcDecayMode.png");
  c_mymcDecayMode->Close();

  auto c_E2d_Eta0 = new TCanvas("c_E2d_Eta0", "c_E2d_Eta0", 800, 800);
  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  h_E2d_Eta0->Draw("COLZ");
  h_E2d_Eta0->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_E2d_Eta0->SetRightMargin(0.15);
  c_E2d_Eta0->SaveAs("plots_reco_" + decay + "/c_E2d_Eta0.png");
  c_E2d_Eta0->Close();

  auto c_E2d_Eta1 = new TCanvas("c_E2d_Eta1", "c_E2d_Eta1", 800, 800);
  h_E2d_Eta1->Draw("COLZ");
  h_E2d_Eta1->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_E2d_Eta1->SetRightMargin(0.15);
  c_E2d_Eta1->SaveAs("plots_reco_" + decay + "/c_E2d_Eta1.png");
  c_E2d_Eta1->Close();

  auto c_E2d_Eta2 = new TCanvas("c_E2d_Eta2", "c_E2d_Eta2", 800, 800);
  h_E2d_Eta2->Draw("COLZ");
  h_E2d_Eta2->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_E2d_Eta2->SetRightMargin(0.15);
  c_E2d_Eta2->SaveAs("plots_reco_" + decay + "/c_E2d_Eta2.png");
  c_E2d_Eta2->Close();

  auto c_E2d_Eta3 = new TCanvas("c_E2d_Eta3", "c_E2d_Eta3", 800, 800);
  h_E2d_Eta3->Draw("COLZ");
  h_E2d_Eta3->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_E2d_Eta3->SetRightMargin(0.15);
  c_E2d_Eta3->SaveAs("plots_reco_" + decay + "/c_E2d_Eta3.png");
  c_E2d_Eta3->Close();

  auto c_E2d_Eta4 = new TCanvas("c_E2d_Eta4", "c_E2d_Eta4", 800, 800);
  h_E2d_Eta4->Draw("COLZ");
  h_E2d_Eta4->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_E2d_Eta4->SetRightMargin(0.15);
  c_E2d_Eta4->SaveAs("plots_reco_" + decay + "/c_E2d_Eta4.png");
  c_E2d_Eta4->Close();

  auto c_EvsEta_mcHits =
      new TCanvas("c_EvsEta_mcHits", "c_EvsEta_mcHits", 800, 800);
  h_EvsEta_mcHits->Draw("COLZ");
  h_EvsEta_mcHits->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  c_E2d_Eta4->SetRightMargin(0.15);
  c_EvsEta_mcHits->SaveAs("plots_reco_" + decay + "/c_EvsEta_mcHits.png");
  c_EvsEta_mcHits->Close();

  auto c_EvsEta_mcAll =
      new TCanvas("c_EvsEta_mcAll", "c_EvsEta_mcAll", 800, 800);
  h_EvsEta_mcAll->Draw("COLZ");
  h_EvsEta_mcAll->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogx();
  c_E2d_Eta4->SetRightMargin(0.15);
  c_EvsEta_mcAll->SaveAs("plots_reco_" + decay + "/c_EvsEta_mcAll.png");
  c_EvsEta_mcAll->Close();

  auto c_EvsEta = new TCanvas("c_EvsEta", "c_EvsEta", 1200, 600);
  c_EvsEta->Divide(2, 1);
  c_EvsEta->cd(1);
  h_EvsEta_mcHits->Draw("COLZ");
  c_EvsEta->cd(2);
  h_EvsEta_mcAll->Draw("COLZ");
  c_EvsEta->SaveAs("plots_reco_" + decay + "/c_EvsEta.png");
  c_EvsEta->Close();

  if (TEfficiency::CheckConsistency(*h_EvsEta_mcHits, *h_EvsEta_mcAll)) {
    eff_EvsEta = new TEfficiency(*h_EvsEta_mcHits, *h_EvsEta_mcAll);
    // eff_EvsEta->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta->SetConfidenceLevel(0.68);
    auto c_EvsEta_eff = new TCanvas("c_EvsEta_eff", "c_EvsEta_eff", 900, 800);
    eff_EvsEta->Draw("EY");
    c_EvsEta_eff->Update();
    auto heff_EvsEta = eff_EvsEta->GetPaintedHistogram();
    // heff_EvsEta->SetMinimum(-0.0001);
    heff_EvsEta->Draw("COLZ1");
    heff_EvsEta->GetXaxis()->SetTitleOffset(1.3);
    c_EvsEta_eff->SetRightMargin(0.15);
    c_EvsEta_eff->SetBottomMargin(0.15);
    gPad->SetLogx();
    c_EvsEta_eff->Update();
    c_EvsEta_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta_eff.png");
    c_EvsEta_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EvsEta0_mcHits, *h_EvsEta0_mcAll)) {
    auto c_EvsEta0_eff =
        new TCanvas("c_EvsEta0_eff", "c_EvsEta0_eff", 900, 700);
    eff_EvsEta0 = new TEfficiency(*h_EvsEta0_mcHits, *h_EvsEta0_mcAll);
    // eff_EvsEta0->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta0->SetConfidenceLevel(0.68);
    eff_EvsEta0->Draw("APL");
    gPad->Update();
    auto graph = eff_EvsEta0->GetPaintedGraph();
    eff_EvsEta0->SetTitle(
        "; p_{T, gen}^{vis}(#tau) [GeV/c]; #tau reconstruction efficiency");

    gPad->Update();
    gPad->SetLogx();

    // Set marker properties for the graph points
    graph->SetMarkerStyle(
        20); // Change the marker style (e.g., 20 for filled circle)
    graph->SetMarkerSize(1.2); // Increase the marker size for better visibility
    graph->SetMarkerColor(kBlue); // Set the marker color (e.g., kBlue for blue)

    // Set line properties for the error bars
    graph->SetLineColor(kBlue); // Set the line color for error bars
    graph->SetLineWidth(1);     // Set the line width for error bars

    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    eff_EvsEta0->SetTitle(
        "; p_{T, gen}^{vis}(#tau) [GeV/c]; #tau reconstruction efficiency");
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetLabelOffset(0.012);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetLabelOffset(0.0125);

    graph->Draw("AP");

    TPaveText *pt2 = new TPaveText(7, 0.8, 27, 0.95);
    pt2->SetFillColor(0);
    pt2->SetMargin(0);
    pt2->SetFillStyle(0);
    pt2->SetBorderSize(0);
    pt2->SetTextAlign(11);
    // pt2->SetTextFont(50);
    pt2->AddText("Muon Collider");
    pt2->AddText("#bf{#it{Simulation}}");
    // pt2->AddText("  ");
    //  pt2->AddText("#bf{}");
    pt2->Draw();

    TPaveText *pt3 = new TPaveText(120, 0.8, 600, 0.95);
    pt3->SetFillColor(0);
    pt3->SetMargin(0);
    pt3->SetFillStyle(0);
    pt3->SetBorderSize(0);
    pt3->SetTextAlign(31);
    // pt3->SetTextFont(50);
    pt3->AddText("#bf{#tau_{h} 1-prong decays}");
pt3->AddText("#bf{0 #leq |#eta_{gen}^{vis} (#tau)| #leq 0.5}");
    // pt3->AddText("#bf{}");
    pt3->Draw();

    c_EvsEta0_eff->SetLeftMargin(0.12);
    c_EvsEta0_eff->SetRightMargin(0.03);
    c_EvsEta0_eff->SetBottomMargin(0.14);

    c_EvsEta0_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta0_eff.pdf");
    c_EvsEta0_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EvsEta1_mcHits, *h_EvsEta1_mcAll)) {
    auto c_EvsEta1_eff =
        new TCanvas("c_EvsEta1_eff", "c_EvsEta1_eff", 800, 800);
    eff_EvsEta1 = new TEfficiency(*h_EvsEta1_mcHits, *h_EvsEta1_mcAll);
    // eff_EvsEta1->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta1->SetConfidenceLevel(0.68);
    eff_EvsEta1->Draw("APL");
    gPad->Update();
    auto graph = eff_EvsEta1->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EvsEta1->SetTitle("Tau ID efficiency (0.5 < |#eta_{gen}^{vis}| < 1); "
                          "p_{T, gen}^{vis} [GeV/c]; Efficiency");
    gPad->Update();
    gPad->SetLogx();
    c_EvsEta1_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta1_eff.png");
    c_EvsEta1_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EvsEta2_mcHits, *h_EvsEta2_mcAll)) {
    auto c_EvsEta2_eff =
        new TCanvas("c_EvsEta2_eff", "c_EvsEta2_eff", 800, 800);
    eff_EvsEta2 = new TEfficiency(*h_EvsEta2_mcHits, *h_EvsEta2_mcAll);
    gStyle->SetOptStat(111);
    // eff_EvsEta2->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta2->SetConfidenceLevel(0.68);
    eff_EvsEta2->Draw("APL");
    gPad->Update();
    auto graph = eff_EvsEta2->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    gPad->SetLogx();
    eff_EvsEta2->SetTitle("Tau ID efficiency (1 < |#eta_{gen}^{vis}| < 1.5); "
                          "p_{T, gen}^{vis} [GeV/c]; Efficiency");
    gPad->Update();
    c_EvsEta2_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta2_eff.png");
    c_EvsEta2_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EvsEta3_mcHits, *h_EvsEta3_mcAll)) {
    auto c_EvsEta3_eff =
        new TCanvas("c_EvsEta3_eff", "c_EvsEta3_eff", 800, 800);
    eff_EvsEta3 = new TEfficiency(*h_EvsEta3_mcHits, *h_EvsEta3_mcAll);
    gStyle->SetOptStat(111);
    // eff_EvsEta3->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta3->SetConfidenceLevel(0.68);
    eff_EvsEta3->Draw("APL");
    gPad->Update();
    auto graph = eff_EvsEta3->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    gPad->SetLogx();
    eff_EvsEta3->SetTitle("Tau ID efficiency (1.5 < |#eta_{gen}^{vis}| < 2); "
                          "p_{T, gen}^{vis} [GeV/c]; Efficiency");
    gPad->Update();
    c_EvsEta3_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta3_eff.png");
    c_EvsEta3_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EvsEta4_mcHits, *h_EvsEta4_mcAll)) {
    auto c_EvsEta4_eff =
        new TCanvas("c_EvsEta4_eff", "c_EvsEta4_eff", 800, 800);
    eff_EvsEta4 = new TEfficiency(*h_EvsEta4_mcHits, *h_EvsEta4_mcAll);
    gStyle->SetOptStat(111);
    // eff_EvsEta4->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EvsEta4->SetConfidenceLevel(0.68);
    eff_EvsEta4->Draw("APL");
    gPad->Update();
    auto graph = eff_EvsEta4->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    gPad->SetLogx();
    eff_EvsEta4->SetTitle("Tau ID efficiency (2 < |#eta_{gen}^{vis}|); "
                          "p_{T, gen}^{vis} [GeV/c]; Efficiency");
    gPad->Update();
    c_EvsEta4_eff->SaveAs("plots_reco_" + decay + "/c_EvsEta4_eff.png");
    c_EvsEta4_eff->Close();
  }

  for (int i = 0; i != h_EtavsE0_mcHits->GetNbinsX() + 2; ++i) {
    if (h_EtavsE0_mcHits->GetBinContent(i) >
        h_EtavsE0_mcAll->GetBinContent(i)) {

      std::cout << "h_EtavsE0_mcHits->GetBinContent(i)"
                << h_EtavsE0_mcHits->GetBinContent(i) << "   bin" << i << '\n';
      std::cout << "h_EtavsE0_mcAll->GetBinContent(i)"
                << h_EtavsE0_mcAll->GetBinContent(i) << "  bin" << i << '\n';
      std::cout << "(0) Hits > All!" << '\n';
    }
    if (h_EtavsE0_mcHits->GetBinLowEdge(i) != h_EtavsE0_mcAll->GetBinLowEdge(i))
      std::cout << "(0) Wrong binning!" << '\n';
  }
  for (int i = 0; i != h_EtavsE1_mcHits->GetNbinsX() + 2; ++i) {
    if (h_EtavsE1_mcHits->GetBinContent(i) > h_EtavsE1_mcAll->GetBinContent(i))
      std::cout << "(1) Hits > All!" << '\n';
    if (h_EtavsE1_mcHits->GetBinLowEdge(i) != h_EtavsE1_mcAll->GetBinLowEdge(i))
      std::cout << "(1) Wrong binning!" << '\n';
  }
  for (int i = 0; i != h_EtavsE2_mcHits->GetNbinsX() + 2; ++i) {
    if (h_EtavsE2_mcHits->GetBinContent(i) > h_EtavsE2_mcAll->GetBinContent(i))
      std::cout << "(2) Hits > All!" << '\n';
    if (h_EtavsE2_mcHits->GetBinLowEdge(i) != h_EtavsE2_mcAll->GetBinLowEdge(i))
      std::cout << "(2) Wrong binning!" << '\n';
  }
  for (int i = 0; i != h_EtavsE3_mcHits->GetNbinsX() + 2; ++i) {
    if (h_EtavsE3_mcHits->GetBinContent(i) > h_EtavsE3_mcAll->GetBinContent(i))
      std::cout << "(3) Hits > All!" << '\n';
    if (h_EtavsE3_mcHits->GetBinLowEdge(i) != h_EtavsE3_mcAll->GetBinLowEdge(i))
      std::cout << "(3) Wrong binning!" << '\n';
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE0_mcHits, *h_EtavsE0_mcAll)) {
    auto c_EtavsE0_eff =
        new TCanvas("c_EtavsE0_eff", "c_EtavsE0_eff", 800, 800);
    eff_EtavsE0 = new TEfficiency(*h_EtavsE0_mcHits, *h_EtavsE0_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE0->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE0->SetConfidenceLevel(0.68);
    eff_EtavsE0->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE0->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE0->SetTitle("Tau ID efficiency (0 < p_{T, gen}^{vis} < 50 GeV); "
                          "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE0_eff);
    gPad->Update();
    c_EtavsE0_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE0_eff.png");
    c_EtavsE0_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE1_mcHits, *h_EtavsE1_mcAll)) {
    auto c_EtavsE1_eff =
        new TCanvas("c_EtavsE1_eff", "c_EtavsE1_eff", 900, 700);
    eff_EtavsE1 = new TEfficiency(*h_EtavsE1_mcHits, *h_EtavsE1_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE1->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE1->SetConfidenceLevel(0.68);
    eff_EtavsE1->Draw("AP");
    gPad->Update();

    auto graph = eff_EtavsE1->GetPaintedGraph();

    // Set marker properties for the graph points
    graph->SetMarkerStyle(
        20); // Change the marker style (e.g., 20 for filled circle)
    graph->SetMarkerSize(1.2); // Increase the marker size for better visibility
    graph->SetMarkerColor(kBlue); // Set the marker color (e.g., kBlue for blue)

    // Set line properties for the error bars
    graph->SetLineColor(kBlue); // Set the line color for error bars
    graph->SetLineWidth(1);     // Set the line width for error bars

    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    eff_EtavsE1->SetTitle(
        ";|#eta_{gen}^{vis}(#tau)|; #tau reconstruction efficiency");
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetLabelOffset(0.012);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetLabelOffset(0.0125);

    graph->Draw("AP");

    graph->SetTitle(
        ";|#eta_{gen}^{vis}(#tau)|; #tau reconstruction efficiency");

    gPad->Modified();
    gPad->Update();
    plot_boundaries(decay, c_EtavsE1_eff);
    gPad->Update();

    c_EtavsE1_eff->SetLeftMargin(0.12);
    c_EtavsE1_eff->SetRightMargin(0.03);
    c_EtavsE1_eff->SetBottomMargin(0.12);

    TPaveText *pt1 = new TPaveText(1.5, 0.65, 2.7, 0.97);
    pt1->SetFillColor(0);
    pt1->SetMargin(0);
    pt1->SetFillStyle(0);
    pt1->SetBorderSize(0);
    pt1->SetTextAlign(11);
    // pt1->SetTextFont(50);
    pt1->AddText("Muon Collider");
    pt1->AddText("#bf{#it{Simulation}}");
    // pt1->AddText("  ");
    pt1->AddText("#bf{#tau_{h} 1-prong decays}");
    pt1->AddText("#bf{50 GeV < p_{T, gen}^{vis} (#tau) < 100 GeV}");
    // pt1->AddText("#bf{}");
    pt1->Draw();
    gPad->Modified();
    gPad->Update();

    c_EtavsE1_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE1_eff.pdf");
    c_EtavsE1_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE2_mcHits, *h_EtavsE2_mcAll)) {
    auto c_EtavsE2_eff =
        new TCanvas("c_EtavsE2_eff", "c_EtavsE2_eff", 800, 800);
    eff_EtavsE2 = new TEfficiency(*h_EtavsE2_mcHits, *h_EtavsE2_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE2->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE2->SetConfidenceLevel(0.68);
    eff_EtavsE2->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE2->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE2->SetTitle(
        "Tau ID efficiency (100 < p_{T, gen}^{vis} < 150 GeV); "
        "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE2_eff);
    gPad->Update();
    c_EtavsE2_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE2_eff.png");
    c_EtavsE2_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE3_mcHits, *h_EtavsE3_mcAll)) {
    auto c_EtavsE3_eff =
        new TCanvas("c_EtavsE3_eff", "c_EtavsE3_eff", 800, 800);
    eff_EtavsE3 = new TEfficiency(*h_EtavsE3_mcHits, *h_EtavsE3_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE3->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE3->SetConfidenceLevel(0.68);
    eff_EtavsE3->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE3->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE3->SetTitle(
        "Tau ID efficiency (150 < p_{T, gen}^{vis} < 200 GeV); "
        "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE3_eff);
    gPad->Update();
    c_EtavsE3_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE3_eff.png");
    c_EtavsE3_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE4_mcHits, *h_EtavsE4_mcAll)) {
    auto c_EtavsE4_eff =
        new TCanvas("c_EtavsE4_eff", "c_EtavsE4_eff", 800, 800);
    eff_EtavsE4 = new TEfficiency(*h_EtavsE4_mcHits, *h_EtavsE4_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE4->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE4->SetConfidenceLevel(0.68);
    eff_EtavsE4->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE4->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE4->SetTitle(
        "Tau ID efficiency (200 < p_{T, gen}^{vis} < 300 GeV); "
        "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE4_eff);
    gPad->Update();
    c_EtavsE4_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE4_eff.png");
    c_EtavsE4_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE5_mcHits, *h_EtavsE5_mcAll)) {
    auto c_EtavsE5_eff =
        new TCanvas("c_EtavsE5_eff", "c_EtavsE5_eff", 800, 800);
    eff_EtavsE5 = new TEfficiency(*h_EtavsE5_mcHits, *h_EtavsE5_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE5->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE5->SetConfidenceLevel(0.68);
    eff_EtavsE5->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE5->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE5->SetTitle(
        "Tau ID efficiency (300 < p_{T, gen}^{vis} < 500 GeV); "
        "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE5_eff);
    gPad->Update();
    c_EtavsE5_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE5_eff.png");
    c_EtavsE5_eff->Close();
  }

  if (TEfficiency::CheckConsistency(*h_EtavsE6_mcHits, *h_EtavsE6_mcAll)) {
    auto c_EtavsE6_eff =
        new TCanvas("c_EtavsE6_eff", "c_EtavsE6_eff", 800, 800);
    eff_EtavsE6 = new TEfficiency(*h_EtavsE6_mcHits, *h_EtavsE6_mcAll);
    gStyle->SetOptStat(111);
    // eff_EtavsE6->SetStatisticOption(TEfficiency::kBBayesian);
    eff_EtavsE6->SetConfidenceLevel(0.68);
    eff_EtavsE6->Draw("APL");
    gPad->Update();
    auto graph = eff_EtavsE6->GetPaintedGraph();
    graph->SetMinimum(0.);
    graph->SetMaximum(1.);
    graph->GetXaxis()->SetTitleOffset(1.2);
    eff_EtavsE6->SetTitle("Tau ID efficiency (p_{T, gen}^{vis} > 500 GeV); "
                          "|#eta_{gen}^{vis}|; Efficiency");
    gPad->Update();
    plot_boundaries(decay, c_EtavsE6_eff);
    gPad->Update();
    c_EtavsE6_eff->SaveAs("plots_reco_" + decay + "/c_EtavsE6_eff.png");
    c_EtavsE6_eff->Close();
  }

  auto c_EvsEvis = new TCanvas("c_EvsEvis", "c_EvsEvis", 800, 800);
  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  h_EvsEvis->Draw("COLZ");
  h_EvsEvis->GetXaxis()->SetTitleOffset(1.3);
  h_EvsEvis->GetYaxis()->SetRangeUser(0., 2.1);
  c_EvsEvis->SetRightMargin(0.15);
  // gPad->SetLogz(1);
  c_EvsEvis->SetLeftMargin(0.1);
  c_EvsEvis->SaveAs("plots_reco_" + decay + "/c_EvsEvis.png");
  c_EvsEvis->Close();

  EProfile_Eta0 = h_E2d_Eta0->ProfileX("EProfile_Eta0");
  EProfile_Eta1 = h_E2d_Eta1->ProfileX("EProfile_Eta1");
  EProfile_Eta2 = h_E2d_Eta2->ProfileX("EProfile_Eta2");
  EProfile_Eta3 = h_E2d_Eta3->ProfileX("EProfile_Eta3");
  EProfile_Eta4 = h_E2d_Eta4->ProfileX("EProfile_Eta4");

  for (int i = 1; i != nlogbins_E + 1; ++i) {
    if (EProfile_Eta0->GetBinEntries(i) > 0) {
      h_ECorr_Eta0->AddBinContent(i, EProfile_Eta0->GetBinContent(i) /
                                         h_ECorr_Eta0->GetBinCenter(i));
      h_ECorr_Eta0->SetBinError(i, EProfile_Eta0->GetBinError(i) /
                                       h_ECorr_Eta0->GetBinCenter(i));
    }
    if (EProfile_Eta1->GetBinEntries(i) > 0) {
      h_ECorr_Eta1->AddBinContent(i, EProfile_Eta1->GetBinContent(i) /
                                         h_ECorr_Eta1->GetBinCenter(i));
      h_ECorr_Eta1->SetBinError(i, EProfile_Eta1->GetBinError(i) /
                                       h_ECorr_Eta1->GetBinCenter(i));
    }
    if (EProfile_Eta2->GetBinEntries(i) > 0) {
      h_ECorr_Eta2->AddBinContent(i, EProfile_Eta2->GetBinContent(i) /
                                         h_ECorr_Eta2->GetBinCenter(i));
      h_ECorr_Eta2->SetBinError(i, EProfile_Eta2->GetBinError(i) /
                                       h_ECorr_Eta2->GetBinCenter(i));
    }
    if (EProfile_Eta3->GetBinEntries(i) > 0) {
      h_ECorr_Eta3->AddBinContent(i, EProfile_Eta3->GetBinContent(i) /
                                         h_ECorr_Eta3->GetBinCenter(i));
      h_ECorr_Eta3->SetBinError(i, EProfile_Eta3->GetBinError(i) /
                                       h_ECorr_Eta3->GetBinCenter(i));
    }
  }

  auto c_ECorr_Eta0 = new TCanvas("c_ECorr_Eta0", "c_ECorr_Eta0", 1100, 800);
  gStyle->SetOptStat(0);
  h_ECorr_Eta0->Draw("E1");
  gPad->SetLogx();

  c_ECorr_Eta0->SetLeftMargin(0.18);
  c_ECorr_Eta0->SetBottomMargin(0.15);
  c_ECorr_Eta0->SetRightMargin(0.05);

  // gPad->SetLogz(1);
  // Set marker properties for the h_ECorr_Eta0 points
  h_ECorr_Eta0->SetMarkerStyle(
      20); // Change the marker style (e.g., 20 for filled circle)
  h_ECorr_Eta0->SetMarkerSize(
      1.2); // Increase the marker size for better visibility
  h_ECorr_Eta0->SetMarkerColor(
      kBlue); // Set the marker color (e.g., kBlue for blue)

  // Set line properties for the error bars
  h_ECorr_Eta0->SetLineColor(kBlue); // Set the line color for error bars
  h_ECorr_Eta0->SetLineWidth(1);     // Set the line width for error bars

  gStyle->SetStripDecimals(0);

  h_ECorr_Eta0->SetTitle("");
  h_ECorr_Eta0->SetMinimum(0.5);
  h_ECorr_Eta0->SetMaximum(3.5);
  h_ECorr_Eta0->GetXaxis()->SetTitleOffset(1.3);
  h_ECorr_Eta0->GetXaxis()->SetTitleSize(0.045);
  h_ECorr_Eta0->GetXaxis()->SetLabelSize(0.045);
  h_ECorr_Eta0->GetXaxis()->SetLabelOffset(0.012);
  h_ECorr_Eta0->GetYaxis()->SetTitleSize(0.045);
  h_ECorr_Eta0->GetYaxis()->SetTitleOffset(1.8);
  h_ECorr_Eta0->GetYaxis()->SetLabelSize(0.045);
  h_ECorr_Eta0->GetYaxis()->SetLabelOffset(0.0125);

  TPaveText *pt4 = new TPaveText(80, 1.9, 600, 3.2);
  pt4->SetFillColor(0);
  pt4->SetMargin(0);
  pt4->SetFillStyle(0);
  pt4->SetBorderSize(0);
  pt4->SetTextAlign(11);
  // pt4->SetTextFont(50);
  pt4->AddText("Muon Collider");
  pt4->AddText("#bf{#it{Simulation}}");
  pt4->AddText("  ");
  pt4->AddText("#bf{h^{#pm} + neutrals decays}");
  pt4->AddText("#bf{0 #leq |#eta_{gen}^{vis}(#tau)| #leq 0.5}");
  // pt4->AddText("#bf{}");
  pt4->Draw();
  gPad->Modified();
  gPad->Update();

  TLine *line = new TLine(h_ECorr_Eta0->GetXaxis()->GetXmin(), 1,
                          h_ECorr_Eta0->GetXaxis()->GetXmax(), 1);
  line->SetLineColor(kRed); // Set the line color to red
  line->SetLineStyle(2);    // Set the line style to dotted (2)
  line->Draw();

  c_ECorr_Eta0->SaveAs("plots_reco_" + decay + "/c_ECorr_Eta0.pdf");
  c_ECorr_Eta0->Close();

  auto c_ECorr_Eta1 = new TCanvas("c_ECorr_Eta1", "c_ECorr_Eta1", 800, 800);
  gStyle->SetOptStat(0);
  h_ECorr_Eta1->GetXaxis()->SetTitleOffset(1.2);
  h_ECorr_Eta1->SetMinimum(0.95);
  h_ECorr_Eta1->Draw("E1");
  // gPad->SetLogz(1);
  gPad->SetLogx();
  c_ECorr_Eta1->SetLeftMargin(0.15);
  c_ECorr_Eta1->SaveAs("plots_reco_" + decay + "/c_ECorr_Eta1.png");
  c_ECorr_Eta1->Close();

  auto c_ECorr_Eta2 = new TCanvas("c_ECorr_Eta2", "c_ECorr_Eta2", 800, 800);
  gStyle->SetOptStat(0);
  h_ECorr_Eta2->GetXaxis()->SetTitleOffset(1.2);
  h_ECorr_Eta2->SetMinimum(0.95);
  h_ECorr_Eta2->Draw("E1");
  // gPad->SetLogz(1);
  gPad->SetLogx();
  c_ECorr_Eta2->SetLeftMargin(0.15);
  c_ECorr_Eta2->SaveAs("plots_reco_" + decay + "/c_ECorr_Eta2.png");
  c_ECorr_Eta2->Close();

  auto c_ECorr_Eta3 = new TCanvas("c_ECorr_Eta3", "c_ECorr_Eta3", 800, 800);
  gStyle->SetOptStat(0);
  h_ECorr_Eta3->GetXaxis()->SetTitleOffset(1.2);
  h_ECorr_Eta3->SetMinimum(0.95);
  h_ECorr_Eta3->Draw("E1");
  // gPad->SetLogz(1);
  gPad->SetLogx();
  c_ECorr_Eta3->SetLeftMargin(0.15);
  c_ECorr_Eta3->SaveAs("plots_reco_" + decay + "/c_ECorr_Eta3.png");
  c_ECorr_Eta3->Close();

  auto c_EProfile_Eta0 =
      new TCanvas("c_EProfile_Eta0", "c_EProfile_Eta0", 800, 800);
  gStyle->SetOptStat(0);
  EProfile_Eta0->GetYaxis()->SetTitle("p_{T, gen}^{vis} [GeV/c]");
  EProfile_Eta0->GetXaxis()->SetTitleOffset(1.2);
  EProfile_Eta0->Draw();
  // gPad->SetLogz(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_EProfile_Eta0->SaveAs("plots_reco_" + decay + "/c_EProfile_Eta0.png");
  c_EProfile_Eta0->Close();

  auto c_EProfile_Eta1 =
      new TCanvas("c_EProfile_Eta1", "c_EProfile_Eta1", 800, 800);
  gStyle->SetOptStat(0);
  EProfile_Eta1->GetYaxis()->SetTitle("p_{T, gen}^{vis} [GeV/c]");
  EProfile_Eta1->GetXaxis()->SetTitleOffset(1.2);
  EProfile_Eta1->Draw();
  // gPad->SetLogz(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_EProfile_Eta1->SaveAs("plots_reco_" + decay + "/c_EProfile_Eta1.png");
  c_EProfile_Eta1->Close();

  auto c_EProfile_Eta2 =
      new TCanvas("c_EProfile_Eta2", "c_EProfile_Eta2", 800, 800);
  gStyle->SetOptStat(0);
  EProfile_Eta2->GetYaxis()->SetTitle("p_{T, gen}^{vis} [GeV/c]");
  EProfile_Eta2->GetXaxis()->SetTitleOffset(1.2);
  EProfile_Eta2->Draw();
  // gPad->SetLogz(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_EProfile_Eta2->SaveAs("plots_reco_" + decay + "/c_EProfile_Eta2.png");
  c_EProfile_Eta2->Close();

  auto c_EProfile_Eta3 =
      new TCanvas("c_EProfile_Eta3", "c_EProfile_Eta3", 800, 800);
  gStyle->SetOptStat(0);
  EProfile_Eta3->GetYaxis()->SetTitle("p_{T, gen}^{vis} [GeV/c]");
  EProfile_Eta3->GetXaxis()->SetTitleOffset(1.2);
  EProfile_Eta3->Draw();
  // gPad->SetLogz(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_EProfile_Eta3->SaveAs("plots_reco_" + decay + "/c_EProfile_Eta3.png");
  c_EProfile_Eta3->Close();

  auto c_EProfile_Eta4 =
      new TCanvas("c_EProfile_Eta4", "c_EProfile_Eta4", 800, 800);
  gStyle->SetOptStat(0);
  EProfile_Eta4->GetYaxis()->SetTitle("p_{T, gen}^{vis} [GeV/c]");
  EProfile_Eta4->GetXaxis()->SetTitleOffset(1.2);
  EProfile_Eta4->Draw();
  // gPad->SetLogz(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1, 1);
  c_EProfile_Eta4->SaveAs("plots_reco_" + decay + "/c_EProfile_Eta4.png");
  c_EProfile_Eta4->Close();

  TList *histo_list = new TList();
  histo_list->Add(h_ECorr_Eta0);
  histo_list->Add(h_ECorr_Eta1);
  histo_list->Add(h_ECorr_Eta2);
  histo_list->Add(h_ECorr_Eta3);
  TFile *rootfile = new TFile("./histo_corrections.root", "UPDATE");
  histo_list->Write(decay + "_corrections", TObject::kSingleKey);
  rootfile->Close();

  // Create a new .root file
  TFile *outputFile = new TFile("../efficiencies.root", "UPDATE");

  // Create a directory (folder) named "piguns_notrackfiltering" inside the root
  // file
  TDirectory *dir = outputFile->GetDirectory("tauguns_notrackfiltering_N");
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

  delete[] tau_charge;
  delete[] recE;
  delete[] mcEvis;
  delete[] mcPt;
  delete[] mcPtvis;
  delete[] recPt;
  delete[] mcEta;
  delete[] recEta;
  delete[] mcPhi;
  delete[] recPhi;
  delete[] mcDecayMode;
  delete[] recNPfos;
  delete[] pfosDeltaR;
  delete[] pfosPdg;
  delete[] pfosPt;

  h_mcDecayMode->Delete();
  h_mymcDecayMode->Delete();
  h_recDecayMode->Delete();
  h_myrecDecayMode->Delete();
  h_E2d_Eta0->Delete();
  h_E2d_Eta1->Delete();
  h_E2d_Eta2->Delete();
  h_E2d_Eta3->Delete();
  h_E2d_Eta4->Delete();

  h_EvsEta0_mcAll->Delete();
  h_EvsEta1_mcAll->Delete();
  h_EvsEta2_mcAll->Delete();
  h_EvsEta3_mcAll->Delete();
  h_EvsEta4_mcAll->Delete();
  h_EvsEta0_mcHits->Delete();
  h_EvsEta1_mcHits->Delete();
  h_EvsEta2_mcHits->Delete();
  h_EvsEta3_mcHits->Delete();
  h_EvsEta4_mcHits->Delete();
  h_EtavsE0_mcAll->Delete();
  h_EtavsE1_mcAll->Delete();
  h_EtavsE2_mcAll->Delete();
  h_EtavsE3_mcAll->Delete();
  h_EtavsE4_mcAll->Delete();
  h_EtavsE5_mcAll->Delete();
  h_EtavsE6_mcAll->Delete();
  h_EtavsE0_mcHits->Delete();
  h_EtavsE1_mcHits->Delete();
  h_EtavsE2_mcHits->Delete();
  h_EtavsE3_mcHits->Delete();
  h_EtavsE4_mcHits->Delete();
  h_EtavsE5_mcHits->Delete();
  h_EtavsE6_mcHits->Delete();
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
  eff_EtavsE6->Delete();

  h_EvsEta_mcAll->Delete();
  h_EvsEta_mcHits->Delete();
  eff_EvsEta->Delete();

  h_EvsEvis->Delete();

  h_ECorr_Eta0->Delete();
  h_ECorr_Eta1->Delete();
  h_ECorr_Eta2->Delete();
  h_ECorr_Eta3->Delete();
  h_ECorr_Eta4->Delete();
}

void TauRecoAnalysis_N(const TString filename = "output_EvalTauFinder.root") {
  // Histo_filler("all", filename);
  // Histo_filler("e", filename);
  // Histo_filler("mu", filename);
  Histo_filler("1piN", filename);
  // Histo_filler("3pi", filename);
  /*Histo_filler("1pi0N", filename);
  Histo_filler("1pi1N", filename);
  Histo_filler("1pi2N", filename);
  Histo_filler("1pi3N", filename);*/
  /*Histo_filler("3pi", filename);
  Histo_filler("other", filename);*/
}
