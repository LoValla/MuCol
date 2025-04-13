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

TH1I *h_mcDecayMode;

TH1F *h_EvsEta0_mcAll, *h_EvsEta1_mcAll, *h_EvsEta2_mcAll, *h_EvsEta3_mcAll,
    *h_EvsEta4_mcAll, *h_EvsEta0_mcHits, *h_EvsEta1_mcHits, *h_EvsEta2_mcHits,
    *h_EvsEta3_mcHits, *h_EvsEta4_mcHits;

TH1F *h_EtavsE0_mcAll, *h_EtavsE1_mcAll, *h_EtavsE2_mcAll, *h_EtavsE3_mcAll,
    *h_EtavsE4_mcAll, *h_EtavsE5_mcAll, *h_EtavsE6_mcAll, *h_EtavsE0_mcHits,
    *h_EtavsE1_mcHits, *h_EtavsE2_mcHits, *h_EtavsE3_mcHits, *h_EtavsE4_mcHits,
    *h_EtavsE5_mcHits, *h_EtavsE6_mcHits;

TH2F *h_E2d_Eta0, *h_E2d_Eta1, *h_E2d_Eta2, *h_E2d_Eta3, *h_E2d_Eta4,
    *h_EvsEta_mcAll, *h_EvsEta_mcHits, *h_EvsEvis;

TEfficiency *eff_EvsEta, *eff_EvsEta0, *eff_EvsEta1, *eff_EvsEta2, *eff_EvsEta3,
    *eff_EvsEta4, *eff_EtavsE0, *eff_EtavsE1, *eff_EtavsE2, *eff_EtavsE3,
    *eff_EtavsE4, *eff_EtavsE5, *eff_EtavsE6;

void TauAnalysis(
    const TString &filename = "Outputs/output_EvalTauFinder.root") {

  float etamax = 2.1;
  float etamin = 0.;
  int nbins_eta = 105;
  double Emax = 1000.;
  double Emin = 5.;
  int nbins_E = 50;
  int nlogbins_E = 14;
  float binsemiwidth_E = 0.5 * (Emax - Emin) / nbins_E;
  /*double logbin_edges_E[nlogbins_E + 1];

  for (int i = 0; i <= nlogbins_E; i++) {
    logbin_edges_E[i] =
        pow(10, TMath::Log10(Emin) + (TMath::Log10(Emax) - TMath::Log10(Emin)) /
                                         double(nlogbins_E) * double(i));
    // std::cout << logbin_edges_E[i] << '\n';
  }*/

  double logbin_edges_E[16] = {5.,   10.,  20.,  40.,  60.,  80.,  100., 125.,
                               150., 175., 200., 250., 300., 500., 750., 1000.};

  // histogram initialization
  {
    h_mcDecayMode =
        new TH1I("h_mcDecayMode", "MC #tau decay mode; Decay Mode; Entries", 9,
                 -1.5, 7.5);
    h_E2d_Eta0 =
        new TH2F("h_E2d_0",
                 "Tau reco energy vs Visible gen energy (0 < |#eta| < 0.5); "
                 "E_{Reco} [GeV]; E_{gen}^{vis} [GeV]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta1 =
        new TH2F("h_E2d_Eta1",
                 "Tau reco energy vs Visible gen energy (0.5 < |#eta| < 1); "
                 "E_{Reco} [GeV]; E_{gen}^{vis} [GeV]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta2 =
        new TH2F("h_E2d_Eta2",
                 "Tau reco energy vs Visible gen energy (1. < |#eta| < 1.5); "
                 "E_{Reco} [GeV]; E_{gen}^{vis} [GeV]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta3 =
        new TH2F("h_E2d_Eta3",
                 "Tau reco energy vs Visible gen energy (1.5 < |#eta| < 2); "
                 "E_{Reco} [GeV]; E_{gen}^{vis} [GeV]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);
    h_E2d_Eta4 =
        new TH2F("h_E2d_Eta4",
                 "Tau reco energy vs Visible gen energy (2 < |#eta| < 2.5); "
                 "E_{Reco} [GeV]; E_{gen}^{vis} [GeV]",
                 nlogbins_E, logbin_edges_E, nlogbins_E, logbin_edges_E);

    h_EvsEta0_mcAll = new TH1F(
        "h_EvsEta0_mcAll",
        "Tau E events (All, 0 < |#eta| < 0.5); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta1_mcAll = new TH1F(
        "h_EvsEta1_mcAll",
        "Tau events (All, 0.5 < |#eta| < 1); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta2_mcAll = new TH1F(
        "h_EvsEta2_mcAll",
        "Tau events (All, 1 < |#eta| < 1.5); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta3_mcAll = new TH1F(
        "h_EvsEta3_mcAll",
        "Tau events (All, 1.5 < |#eta| < 2); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta4_mcAll =
        new TH1F("h_EvsEta4_mcAll",
                 "Tau events (All, |#eta| > 2); E_{gen}^{vis} [GeV]; Events",
                 nlogbins_E, logbin_edges_E);

    h_EvsEta0_mcHits = new TH1F(
        "h_EvsEta0_mcHits",
        "Tau events (Hits, 0 < |#eta| < 0.5); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta1_mcHits = new TH1F(
        "h_EvsEta1_mcHits",
        "Tau events (Hits, 0.5 < |#eta| < 1); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta2_mcHits = new TH1F(
        "h_EvsEta2_mcHits",
        "Tau events (Hits, 1 < |#eta| < 1.5); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta3_mcHits = new TH1F(
        "h_EvsEta3_mcHits",
        "Tau events (Hits, 1.5 < |#eta| < 2); E_{gen}^{vis} [GeV]; Events",
        nlogbins_E, logbin_edges_E);
    h_EvsEta4_mcHits = new TH1F(
        "h_EvsEta4_mcHits",
        "Tau events (Hits, 2 < |#eta| < 2.5); E_{gen}^{vis} [GeV]; Events",
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
                 "E_{gen}^{vis} [GeV]; |#eta_{gen}|",
                 nlogbins_E, logbin_edges_E, nbins_eta, etamin, etamax);
    h_EvsEta_mcHits =
        new TH2F("h_EvsEta_mcHits",
                 "Tau visible E vs #eta (MC, reco particles); "
                 "E_{gen}^{vis} [GeV]; |#eta_{gen}|",
                 nlogbins_E, logbin_edges_E, nbins_eta, etamin, etamax);

    h_EvsEvis =
        new TH2F("h_EvsEvis",
                 "#tau E generated vs E visible (MC); "
                 "E_{gen}^{#tau} [GeV]; E_{gen}^{vis} [GeV];",
                 nbins_E + 1, Emin - binsemiwidth_E, Emax + binsemiwidth_E,
                 nbins_E + 1, Emin - binsemiwidth_E, Emax + binsemiwidth_E);
  }

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
  float *mcEvis = new float[10000000];
  float *mcE = new float[10000000];
  float *recE = new float[10000000];
  float *mcEta = new float[10000000];
  float *recEta = new float[10000000];
  int *mcDecayMode = new int[10000000];
  int *tau_charge = new int[10000000];
  int *recNPfos = new int[10000000];
  int *pfosPdg = new int[50000000];
  float *pfosPt = new float[50000000];
  float *pfosDeltaR = new float[50000000];

  tree->SetBranchAddress("nTausMC", &ntau_mc);
  tree->SetBranchAddress("nTausRec", &ntau_rec);
  tree->SetBranchAddress("EvID", &evID);
  tree->SetBranchAddress("RunID", &runID);
  tree->SetBranchAddress("mcDecayMode", mcDecayMode);
  tree->SetBranchAddress("mcE", mcE);
  tree->SetBranchAddress("mcE_vis", mcEvis);
  tree->SetBranchAddress("recE", recE);
  tree->SetBranchAddress("mcEta_vis", mcEta);
  tree->SetBranchAddress("recEta", recEta);
  tree->SetBranchAddress("recNPfos", recNPfos);
  tree->SetBranchAddress("charge", tau_charge);
  tree->SetBranchAddress("pfosPdg", pfosPdg);
  tree->SetBranchAddress("pfosPt", pfosPt);
  tree->SetBranchAddress("pfosDeltaR", pfosDeltaR);

  const long int nEntries = tree->GetEntries();

  TString plotsfolder = "plots";

  for (int ientry = 0; ientry < nEntries; ++ientry) {
    if (ientry % 10000 == 0)
      std::cout << ientry << " / " << nEntries << '\n';

    tree->GetEntry(ientry);

    // if (ntau_mc != 1 || ntau_rec != 1)
    // continue;

    int n_el = 0;
    int n_mu = 0;
    int n_pipl = 0;
    int n_pimi = 0;
    int n_gamma = 0;

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

    if (abs(tau_charge[0]) != 1)
      continue;

    if (false)
      continue;

    /*
    if (ntau_rec == 1 && ientry < 15) {
      std::cout << "recNPfos[" << ientry << "] = " << recNPfos[0] << '\n';
      for (int i = 0; i != recNPfos[0]; ++i) {
        std::cout << "pfosPt[i] = " << pfosPt[i] << '\n';
        std::cout << "pfosPdg[i] = " << pfosPdg[i] << '\n';
        std::cout << "pfosDeltaR[i] = " << pfosDeltaR[i] << '\n';
      }
    }
    */
    /*
    if (mcDecayMode[0] == -1) {
      std::cout << "Run: " << runID << " - "
                << "Ev: " << evID << '\n';
    }*/

    for (unsigned int i = 0; i < ntau_rec; ++i) {
      if (ntau_mc == ntau_rec) {
        h_mcDecayMode->Fill(mcDecayMode[i]);
      }
    }

    for (unsigned int i = 0; i < ntau_mc; ++i) {
      // fix eta, draw histos for E
      if (mcEta[i] >= 0. && mcEta[i] <= 0.5) {
        h_EvsEta0_mcAll->Fill(mcEvis[i]);
        if (ntau_rec == ntau_mc) {
          h_E2d_Eta0->Fill(recE[i], mcEvis[i]);
          h_EvsEta0_mcHits->Fill(mcEvis[i]);
        }
      } else if (mcEta[i] >= 0.5 && mcEta[i] <= 1.) {
        h_EvsEta1_mcAll->Fill(mcEvis[i]);
        if (ntau_rec == ntau_mc) {
          h_E2d_Eta1->Fill(recE[i], mcEvis[i]);
          h_EvsEta1_mcHits->Fill(mcEvis[i]);
        }
      } else if (mcEta[i] >= 1. && mcEta[i] <= 1.5) {
        h_EvsEta2_mcAll->Fill(mcEvis[i]);
        if (ntau_rec == ntau_mc) {
          h_E2d_Eta2->Fill(recE[i], mcEvis[i]);
          h_EvsEta2_mcHits->Fill(mcEvis[i]);
        }
      } else if (mcEta[i] >= 1.5 && mcEta[i] <= 2.) {
        h_EvsEta3_mcAll->Fill(mcEvis[i]);
        if (ntau_rec == ntau_mc) {
          h_E2d_Eta3->Fill(recE[i], mcEvis[i]);
          h_EvsEta3_mcHits->Fill(mcEvis[i]);
        }
      } else if (mcEta[i] >= 2. && mcEta[i] <= 2.5) {
        h_EvsEta4_mcAll->Fill(mcEvis[i]);
        if (ntau_rec == ntau_mc) {
          h_E2d_Eta4->Fill(recE[i], mcEvis[i]);
          h_EvsEta4_mcHits->Fill(mcEvis[i]);
        }
      }

      // fix E (compute efficiencies)
      if (mcEvis[i] >= 0. && mcEvis[i] < 50.) {
        h_EtavsE0_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE0_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 50. && mcEvis[i] < 100.) {
        h_EtavsE1_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE1_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 100. && mcEvis[i] < 150.) {
        h_EtavsE2_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE2_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 150. && mcEvis[i] < 200.) {
        h_EtavsE3_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE3_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 200. && mcEvis[i] < 300.) {
        h_EtavsE4_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE4_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 300. && mcEvis[i] < 500.) {
        h_EtavsE5_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE5_mcHits->Fill(mcEta[i]);
        }
      } else if (mcEvis[i] >= 500.) {
        h_EtavsE6_mcAll->Fill(mcEta[i]);
        if (ntau_rec == ntau_mc) {
          h_EtavsE6_mcHits->Fill(mcEta[i]);
        }
      }

      h_EvsEta_mcAll->Fill(mcEvis[i], mcEta[i]);
      if (ntau_rec == ntau_mc) {
        h_EvsEta_mcHits->Fill(mcEvis[i], mcEta[i]);
      }

      h_EvsEvis->Fill(mcE[i], mcEvis[i]);
    }
  }

  {
    const char *dm_labels[9] = {"#pi^{#pm}",
                                "#pi^{#pm} + #pi^{0}",
                                "#pi^{#pm} + 2#pi^{0}",
                                "#pi^{#pm} + 3#pi^{0}",
                                "3#pi^{#pm}",
                                "3#pi^{#pm} + #pi^{0}",
                                "e^{#pm}",
                                "#mu^{#pm}",
                                "Others"};

    auto c_mcDecayMode =
        new TCanvas("c_mcDecayMode", "c_mcDecayMode", 900, 800);
    h_mcDecayMode->Draw();
    // h_mcDecayMode->GetXaxis()->SetBit(TAxis::kLabelsHori);
    for (int i = 1; i <= 9; i++)
      h_mcDecayMode->GetXaxis()->SetBinLabel(i, dm_labels[i - 1]);
    h_mcDecayMode->GetXaxis()->SetLabelSize(0.045);
    c_mcDecayMode->SaveAs(plotsfolder + "/c_mcDecayMode.png");
    c_mcDecayMode->Close();

    auto c_E2d_Eta0 = new TCanvas("c_E2d_Eta0", "c_E2d_Eta0", 900, 800);
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    h_E2d_Eta0->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogy();
    c_E2d_Eta0->SaveAs(plotsfolder + "/c_E2d_Eta0.png");
    c_E2d_Eta0->Close();

    auto c_E2d_Eta1 = new TCanvas("c_E2d_Eta1", "c_E2d_Eta1", 900, 800);
    h_E2d_Eta1->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogy();
    c_E2d_Eta1->SaveAs(plotsfolder + "/c_E2d_Eta1.png");
    c_E2d_Eta1->Close();

    auto c_E2d_Eta2 = new TCanvas("c_E2d_Eta2", "c_E2d_Eta2", 900, 800);
    h_E2d_Eta2->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogy();
    c_E2d_Eta2->SaveAs(plotsfolder + "/c_E2d_Eta2.png");
    c_E2d_Eta2->Close();

    auto c_E2d_Eta3 = new TCanvas("c_E2d_Eta3", "c_E2d_Eta3", 900, 800);
    h_E2d_Eta3->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogy();
    c_E2d_Eta3->SaveAs(plotsfolder + "/c_E2d_Eta3.png");
    c_E2d_Eta3->Close();

    auto c_E2d_Eta4 = new TCanvas("c_E2d_Eta4", "c_E2d_Eta4", 900, 800);
    h_E2d_Eta4->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogy();
    c_E2d_Eta4->SaveAs(plotsfolder + "/c_E2d_Eta4.png");
    c_E2d_Eta4->Close();

    auto c_EvsEta_mcHits =
        new TCanvas("c_EvsEta_mcHits", "c_EvsEta_mcHits", 900, 800);
    h_EvsEta_mcHits->Draw("COLZ");
    gPad->SetLogx();
    c_EvsEta_mcHits->SaveAs(plotsfolder + "/c_EvsEta_mcHits.png");
    c_EvsEta_mcHits->Close();

    auto c_EvsEta_mcAll =
        new TCanvas("c_EvsEta_mcAll", "c_EvsEta_mcAll", 900, 800);
    h_EvsEta_mcAll->Draw("COLZ");
    gPad->SetLogx();
    c_EvsEta_mcAll->SaveAs(plotsfolder + "/c_EvsEta_mcAll.png");
    c_EvsEta_mcAll->Close();

    /*auto c_EvsEta = new TCanvas("c_EvsEta", "c_EvsEta", 1200, 600);
    c_EvsEta->Divide(1, 3);
    c_EvsEta->cd(1);
    h_EvsEta_mcHits->Draw("COLZ");
    c_EvsEta->cd(2);
    h_EvsEta_mcAll->Draw("COLZ");
    c_EvsEta->SaveAs(plotsfolder+"/c_EvsEta.png");
    c_EvsEta->Close();
    */

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
      c_EvsEta_eff->SaveAs(plotsfolder + "/c_EvsEta_eff.png");
      c_EvsEta_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta0_mcHits, *h_EvsEta0_mcAll)) {
      auto c_EvsEta0_eff =
          new TCanvas("c_EvsEta0_eff", "c_EvsEta0_eff", 900, 800);
      eff_EvsEta0 = new TEfficiency(*h_EvsEta0_mcHits, *h_EvsEta0_mcAll);
      // eff_EvsEta0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta0->SetConfidenceLevel(0.68);
      eff_EvsEta0->Draw("AP");
      gPad->SetLogx();
      eff_EvsEta0->SetTitle("Tau ID efficiency (0 < |#eta_{gen}| < 0.5); "
                            "E_{gen}^{vis} [GeV]; Efficiency");
      gPad->Update();
      c_EvsEta0_eff->SaveAs(plotsfolder + "/c_EvsEta0_eff.png");
      c_EvsEta0_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta1_mcHits, *h_EvsEta1_mcAll)) {
      auto c_EvsEta1_eff =
          new TCanvas("c_EvsEta1_eff", "c_EvsEta1_eff", 900, 800);
      eff_EvsEta1 = new TEfficiency(*h_EvsEta1_mcHits, *h_EvsEta1_mcAll);
      // eff_EvsEta1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta1->SetConfidenceLevel(0.68);
      eff_EvsEta1->Draw("AP");
      eff_EvsEta1->SetTitle("Tau ID efficiency (0.5 < |#eta_{gen}| < 1); "
                            "E_{gen}^{vis} [GeV]; Efficiency");
      gPad->Update();
      gPad->SetLogx();
      c_EvsEta1_eff->SaveAs(plotsfolder + "/c_EvsEta1_eff.png");
      c_EvsEta1_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta2_mcHits, *h_EvsEta2_mcAll)) {
      auto c_EvsEta2_eff =
          new TCanvas("c_EvsEta2_eff", "c_EvsEta2_eff", 900, 800);
      eff_EvsEta2 = new TEfficiency(*h_EvsEta2_mcHits, *h_EvsEta2_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta2->SetConfidenceLevel(0.68);
      eff_EvsEta2->Draw("AP");
      gPad->SetLogx();
      eff_EvsEta2->SetTitle("Tau ID efficiency (1 < |#eta_{gen}| < 1.5); "
                            "E_{gen}^{vis} [GeV]; Efficiency");
      gPad->Update();
      c_EvsEta2_eff->SaveAs(plotsfolder + "/c_EvsEta2_eff.png");
      c_EvsEta2_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta3_mcHits, *h_EvsEta3_mcAll)) {
      auto c_EvsEta3_eff =
          new TCanvas("c_EvsEta3_eff", "c_EvsEta3_eff", 900, 800);
      eff_EvsEta3 = new TEfficiency(*h_EvsEta3_mcHits, *h_EvsEta3_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta3->SetConfidenceLevel(0.68);
      eff_EvsEta3->Draw("AP");
      gPad->SetLogx();
      eff_EvsEta3->SetTitle("Tau ID efficiency (1.5 < |#eta_{gen}| < 2); "
                            "E_{gen}^{vis} [GeV]; Efficiency");
      gPad->Update();
      c_EvsEta3_eff->SaveAs(plotsfolder + "/c_EvsEta3_eff.png");
      c_EvsEta3_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EvsEta4_mcHits, *h_EvsEta4_mcAll)) {
      auto c_EvsEta4_eff =
          new TCanvas("c_EvsEta4_eff", "c_EvsEta4_eff", 900, 800);
      eff_EvsEta4 = new TEfficiency(*h_EvsEta4_mcHits, *h_EvsEta4_mcAll);
      gStyle->SetOptStat(111);
      // eff_EvsEta4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EvsEta4->SetConfidenceLevel(0.68);
      eff_EvsEta4->Draw("AP");
      gPad->SetLogx();
      eff_EvsEta4->SetTitle("Tau ID efficiency (2 < |#eta_{gen}| < 2.5); "
                            "E_{gen}^{vis} [GeV]; Efficiency");
      gPad->Update();
      c_EvsEta4_eff->SaveAs(plotsfolder + "/c_EvsEta4_eff.png");
      c_EvsEta4_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE0_mcHits, *h_EtavsE0_mcAll)) {
      auto c_EtavsE0_eff =
          new TCanvas("c_EtavsE0_eff", "c_EtavsE0_eff", 900, 800);
      eff_EtavsE0 = new TEfficiency(*h_EtavsE0_mcHits, *h_EtavsE0_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE0->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE0->SetConfidenceLevel(0.68);
      eff_EtavsE0->Draw("AP");
      eff_EtavsE0->SetTitle("Tau ID efficiency (0 < E_{gen}^{vis} < 50 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE0_eff->SaveAs(plotsfolder + "/c_EtavsE0_eff.png");
      c_EtavsE0_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE1_mcHits, *h_EtavsE1_mcAll)) {
      auto c_EtavsE1_eff =
          new TCanvas("c_EtavsE1_eff", "c_EtavsE1_eff", 900, 800);
      eff_EtavsE1 = new TEfficiency(*h_EtavsE1_mcHits, *h_EtavsE1_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE1->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE1->SetConfidenceLevel(0.68);
      eff_EtavsE1->Draw("AP");
      eff_EtavsE1->SetTitle("Tau ID efficiency (50 < E_{gen}^{vis} < 100 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE1_eff->SaveAs(plotsfolder + "/c_EtavsE1_eff.png");
      c_EtavsE1_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE2_mcHits, *h_EtavsE2_mcAll)) {
      auto c_EtavsE2_eff =
          new TCanvas("c_EtavsE2_eff", "c_EtavsE2_eff", 900, 800);
      eff_EtavsE2 = new TEfficiency(*h_EtavsE2_mcHits, *h_EtavsE2_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE2->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE2->SetConfidenceLevel(0.68);
      eff_EtavsE2->Draw("AP");
      eff_EtavsE2->SetTitle(
          "Tau ID efficiency (100 < E_{gen}^{vis} < 150 GeV); "
          "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE2_eff->SaveAs(plotsfolder + "/c_EtavsE2_eff.png");
      c_EtavsE2_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE3_mcHits, *h_EtavsE3_mcAll)) {
      auto c_EtavsE3_eff =
          new TCanvas("c_EtavsE3_eff", "c_EtavsE3_eff", 900, 800);
      eff_EtavsE3 = new TEfficiency(*h_EtavsE3_mcHits, *h_EtavsE3_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE3->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE3->SetConfidenceLevel(0.68);
      eff_EtavsE3->Draw("AP");
      eff_EtavsE3->SetTitle(
          "Tau ID efficiency (150 < E_{gen}^{vis} < 200 GeV); "
          "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE3_eff->SaveAs(plotsfolder + "/c_EtavsE3_eff.png");
      c_EtavsE3_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE4_mcHits, *h_EtavsE4_mcAll)) {
      auto c_EtavsE4_eff =
          new TCanvas("c_EtavsE4_eff", "c_EtavsE4_eff", 900, 800);
      eff_EtavsE4 = new TEfficiency(*h_EtavsE4_mcHits, *h_EtavsE4_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE4->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE4->SetConfidenceLevel(0.68);
      eff_EtavsE4->Draw("AP");
      eff_EtavsE4->SetTitle(
          "Tau ID efficiency (200 < E_{gen}^{vis} < 300 GeV); "
          "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE4_eff->SaveAs(plotsfolder + "/c_EtavsE4_eff.png");
      c_EtavsE4_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE5_mcHits, *h_EtavsE5_mcAll)) {
      auto c_EtavsE5_eff =
          new TCanvas("c_EtavsE5_eff", "c_EtavsE5_eff", 900, 800);
      eff_EtavsE5 = new TEfficiency(*h_EtavsE5_mcHits, *h_EtavsE5_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE5->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE5->SetConfidenceLevel(0.68);
      eff_EtavsE5->Draw("AP");
      eff_EtavsE5->SetTitle(
          "Tau ID efficiency (300 < E_{gen}^{vis} < 500 GeV); "
          "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE5_eff->SaveAs(plotsfolder + "/c_EtavsE5_eff.png");
      c_EtavsE5_eff->Close();
    }

    if (TEfficiency::CheckConsistency(*h_EtavsE6_mcHits, *h_EtavsE6_mcAll)) {
      auto c_EtavsE6_eff =
          new TCanvas("c_EtavsE6_eff", "c_EtavsE6_eff", 900, 800);
      eff_EtavsE6 = new TEfficiency(*h_EtavsE6_mcHits, *h_EtavsE6_mcAll);
      gStyle->SetOptStat(111);
      // eff_EtavsE6->SetStatisticOption(TEfficiency::kBBayesian);
      eff_EtavsE6->SetConfidenceLevel(0.68);
      eff_EtavsE6->Draw("AP");
      eff_EtavsE6->SetTitle("Tau ID efficiency (E_{gen}^{vis} > 500 GeV); "
                            "|#eta_{gen}|; Efficiency");
      gPad->Update();
      c_EtavsE6_eff->SaveAs(plotsfolder + "/c_EtavsE6_eff.png");
      c_EtavsE6_eff->Close();
    }

    auto c_EvsEvis = new TCanvas("c_EvsEvis", "c_EvsEvis", 900, 800);
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    h_EvsEvis->Draw("COLZ");
    // gPad->SetLogz(1);
    c_EvsEvis->SaveAs(plotsfolder + "/c_EvsEvis.png");
    c_EvsEvis->Close();
  }

  delete[] tau_charge;
  delete[] mcE;
  delete[] mcEvis;
  delete[] recE;
  delete[] mcEta;
  delete[] recEta;
  delete[] mcDecayMode;
  delete[] recNPfos;
  delete[] pfosDeltaR;
  delete[] pfosPdg;
  delete[] pfosPt;

  h_mcDecayMode->Delete();
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

  tree->Delete();
  input_file->Close();
  input_file->Delete();
}