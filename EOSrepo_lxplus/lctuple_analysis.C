#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iomanip>

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TGraph.h"

// --- Histogram declaration

TH1F *h_mc_pdg;
TH1F *h_mc_pt;
TH1F *h_reco_pdg;
TH1F *h_reco_ene;

TH1F *h_deltaeta_mc;
TH1F *h_deltaeta_reco;
TH1F *h_deltaeta_reco_pions;
TH1F *h_deltaeta_reco_nopions;

TH1F *h_deltaphi_mc;
TH1F *h_deltaphi_reco;
TH1F *h_deltaphi_reco_pions;
TH1F *h_deltaphi_reco_nopions;

TH1F *h_deltaR_mc;
TH1F *h_deltaR_reco;
TH1F *h_deltaR_reco_pions;
TH1F *h_deltaR_reco_nopions;

TH1F *h_calohits;
TH1F *h_simcalohits;
TH1F *h_calohitsratio_pions;
TH1F *h_calohitsratio_nopions;

TH1F *h_pt;
TH1F *h_pt_pions;
TH1F *h_pt_nopions;

TH1D *h_pt_unc;
TH1D *h_pt_unc_pions;
TH1D *h_pt_unc_nopions;

TH1F *h_mom_ts;
TH1F *h_mom_unc_ts;
TH1F *h_omega;
TH1F *h_omega_unc;

TH1F *h_deltaR_trst_cl;
TH1F *h_deltaR_trst_cl_ch;
TH1F *h_deltaR_trst_cl_unch;

TGraph *g_pt_pions;
TGraph *g_pt_nopions;

TGraph *g_ptvsene_ch;
TGraph *g_ptvsene_unch;

// ===========================================================================

void lctuple_analysis(const TString filename = "lctuple_gun.root", const long int maxEvents = -1)
{
  // --- Histogram booking
  h_mc_pdg = new TH1F("h_mc_pdg", "Stable MC particles; PDG code", 4500, -2250., 2250.);
  h_mc_pt = new TH1F("h_mc_pt", "p_{T} generated; p_{T} [GeV]", 50, 90., 101.0);
  h_reco_pdg = new TH1F("h_reco_pdg", "Reconstructed particles; PDG code", 4500, -2250., 2250.);
  h_reco_ene = new TH1F("h_reco_ene", "Reco #pi energy; E_{#pi} [GeV]", 40, 98.5, 101.5);

  h_deltaeta_mc = new TH1F("h_deltaeta_mc", "#Delta#eta - MC vs tst (at IP); #Delta#eta", 40, -0.0006, 0.0006);
  h_deltaeta_reco = new TH1F("h_deltaeta_reco", "#Delta#eta - reco vs tst (1st ecal hit); #Delta#eta", 50, -0.01, 0.01);
  h_deltaeta_reco_pions = new TH1F("h_deltaeta_reco_pions", "#Delta#eta - reco vs tst (1st ecal hit), only pions; #Delta#eta", 50, -0.001, 0.001);
  h_deltaeta_reco_nopions = new TH1F("h_deltaeta_reco_nopions", "#Delta#eta - reco vs tst (1st ecal hit), no pions; #Delta#eta", 50, -0.09, 0.09);

  h_deltaphi_mc = new TH1F("h_deltaphi_mc", "#Delta#phi - MC vs tst (at IP); #Delta#phi", 40, -0.01, 0.01);
  h_deltaphi_reco = new TH1F("h_deltaphi_reco", "#Delta#phi - reco vs tst (1st ecal hit); #Delta#phi", 50, -0.1, 0.1);
  h_deltaphi_reco_pions = new TH1F("h_deltaphi_reco_pions", "#Delta#phi - reco vs tst (1st ecal hit), only pions; #Delta#phi", 50, -0.1, 0.1);
  h_deltaphi_reco_nopions = new TH1F("h_deltaphi_reco_nopions", "#Delta#phi - reco vs tst (1st ecal hit), no pions; #Delta#phi", 50, -0.1, 0.1);

  h_deltaR_mc = new TH1F("h_deltaR_mc", "#DeltaR - MC vs tst (at IP); #DeltaR", 40, 0.0, 0.002);
  h_deltaR_reco = new TH1F("h_deltaR_reco", "#DeltaR - reco vs tst (1st ecal hit); #DeltaR", 50, 0.0, 0.1);
  h_deltaR_reco_pions = new TH1F("h_deltaR_reco_pions", "#DeltaR - reco vs tst (1st ecal hit), only pions; #DeltaR", 50, 0.0, 0.1);
  h_deltaR_reco_nopions = new TH1F("h_deltaR_reco_nopions", "#DeltaR - reco vs tst (1st ecal hit), no pions; #DeltaR", 50, 0.0, 0.1);

  h_calohits = new TH1F("h_calohits", "HCal Barrel hits - detected; N hits", 50, 0., 300);
  h_simcalohits = new TH1F("h_simcalohits", "HCal Barrel hits - simulated; N hits", 50, 0., 7000);
  h_calohitsratio_pions = new TH1F("h_calohitsratio_pions", "ratio hits/simhits, only pions; hits/simhits", 50, 0.0, 0.3);
  h_calohitsratio_nopions = new TH1F("h_calohitsratio_nopions", "ratio hits/simhits, no pions; hits/simhits", 50, 0.0, 0.3);

  h_pt = new TH1F("h_pt", "p_{T} of the track; p_{T} [GeV]", 50, 92., 101.0);
  h_pt_pions = new TH1F("h_pt_pions", "pT of the track (pions); p_{T} [GeV]", 50, 92., 101.);
  h_pt_nopions = new TH1F("h_pt_nopions", "pT of the track (no pions) p_{T} [GeV]", 50, 92., 101.);

  h_pt_unc = new TH1D("h_pt_unc", "p_{T} uncertainty (tracking); #sigma_{p_{T}}/^{}p_{T}", 50, 0.0015, 0.004);
  h_pt_unc_pions = new TH1D("h_pt_unc_pions", "p_{T} relative uncertainty (pions); #sigma_{p_{T}}/^{}p_{T}", 50, 0.0015, 0.004);
  h_pt_unc_nopions = new TH1D("h_pt_unc_nopions", "p_{T} uncertainty (no pions); #sigma_{p_{T}}/^{}p_{T}", 50, 0.0015, 0.004);

  h_mom_ts = new TH1F("h_mom_ts", "momentum of the track; |p| [GeV]", 100, 0., 1000.);
  h_mom_unc_ts = new TH1F("h_mom_unc_ts", "uncertainty on momentum of the track; #sigma_{p}/^{}p", 100, 0., 0.01);
  h_omega = new TH1F("h_ome", "Omega; m^{-1}", 50, 0.000, 0.012);
  h_omega_unc = new TH1F("h_ome_unc", "Omega unc; m^{-1}", 50, 0.000, 0.0001);

  h_deltaR_trst_cl = new TH1F("h_deltaR_trst_cl", "#DeltaR - cluster vs tst (1st ecal hit); #DeltaR", 50, 0.0, 0.03);
  h_deltaR_trst_cl_ch = new TH1F("h_deltaR_trst_cl_ch", "#DeltaR - cluster vs tst (1st ecal hit), charged prts; #DeltaR", 50, 0.0, 0.03);
  h_deltaR_trst_cl_unch = new TH1F("h_deltaR_trst_cl_unch", "#DeltaR - cluster vs tst (1st ecal hit), neutral prts; #DeltaR", 50, 0.0, 0.03);

  g_pt_pions = new TGraph();
  g_pt_nopions = new TGraph();

  g_ptvsene_ch = new TGraph();
  g_ptvsene_unch = new TGraph();

  //  --- Open the ntuple file and get the tree

  TFile *input_file = new TFile(filename.Data(), "read");

  TTree *myLCTuple = (TTree *)input_file->Get("MyLCTuple");

  // --- Loop over the ntuple entries

  // --- MC particles
  int n_mcp;
  int *mcp_pdg = new int[1500000];
  int *mcp_genCode = new int[1500000];
  int *mcp_siode = new int[1500000];
  float *mcp_vx = new float[1500000];
  float *mcp_vy = new float[1500000];
  float *mcp_vz = new float[1500000];
  float *mcp_px = new float[1500000];
  float *mcp_py = new float[1500000];
  float *mcp_pz = new float[1500000];
  float *mcp_ene = new float[1500000];
  float *mcp_q = new float[1500000];
  float *mcp_t = new float[1500000];

  myLCTuple->SetBranchAddress("nmcp", &n_mcp);
  myLCTuple->SetBranchAddress("mcpdg", mcp_pdg);
  myLCTuple->SetBranchAddress("mcgst", mcp_genCode);
  myLCTuple->SetBranchAddress("mcsst", mcp_siode);
  myLCTuple->SetBranchAddress("mcvtx", mcp_vx);
  myLCTuple->SetBranchAddress("mcvty", mcp_vy);
  myLCTuple->SetBranchAddress("mcvtz", mcp_vz);
  myLCTuple->SetBranchAddress("mcmox", mcp_px);
  myLCTuple->SetBranchAddress("mcmoy", mcp_py);
  myLCTuple->SetBranchAddress("mcmoz", mcp_pz);
  myLCTuple->SetBranchAddress("mcene", mcp_ene);
  myLCTuple->SetBranchAddress("mccha", mcp_q);
  myLCTuple->SetBranchAddress("mctim", mcp_t);

  // additional MC variables
  float *pt_mc = new float[100];
  float *pmag_mc = new float[100];
  float *ene_mc = new float[100];
  float *eta_mc = new float[100];
  float *phi_mc = new float[100];
  float *pdg_mc = new float[100];

  // --- RECO particles
  int n_reco;
  int *reco_type = new int[10000];
  float *reco_px = new float[10000];
  float *reco_py = new float[10000];
  float *reco_pz = new float[10000];
  float *reco_ene = new float[10000];
  float *reco_q = new float[10000];

  myLCTuple->SetBranchAddress("nrec", &n_reco);
  myLCTuple->SetBranchAddress("rctyp", reco_type);
  myLCTuple->SetBranchAddress("rcmox", reco_px);
  myLCTuple->SetBranchAddress("rcmoy", reco_py);
  myLCTuple->SetBranchAddress("rcmoz", reco_pz);
  myLCTuple->SetBranchAddress("rcene", reco_ene);
  myLCTuple->SetBranchAddress("rccha", reco_q);

  // additional reco variables
  float *pmag_rc = new float[1000];
  float *eta_rc = new float[1000];
  float *phi_rc = new float[1000];

  // trackstates
  int n_trst;
  float *trst_phi = new float[40000];
  float *trst_tnl = new float[40000];

  float *trst_ome = new float[40000];
  float *trst_cov = new float[150000];

  float *trst_rpx = new float[40000];
  float *trst_rpy = new float[40000];
  float *trst_rpz = new float[40000];

  myLCTuple->SetBranchAddress("ntrst", &n_trst);
  myLCTuple->SetBranchAddress("tsphi", trst_phi);
  myLCTuple->SetBranchAddress("tstnl", trst_tnl);
  myLCTuple->SetBranchAddress("tsome", trst_ome);
  myLCTuple->SetBranchAddress("tscov", trst_cov);

  myLCTuple->SetBranchAddress("tsrpx", trst_rpx);
  myLCTuple->SetBranchAddress("tsrpy", trst_rpy);
  myLCTuple->SetBranchAddress("tsrpz", trst_rpz);

  // additional trst variables
  float *theta_trst = new float[10000];
  float *eta_trst = new float[10000];
  float *pt_trst = new float[10000];
  float *ptunc_trst = new float[10000];
  float *mom_trst = new float[10000];
  float *mom_unc_trst = new float[10000];
  float delta_r, delta_eta, delta_phi;

  // SimCalHit, Calhit
  int n_cah, n_sch;
  float *cah_ene = new float[100000];
  myLCTuple->SetBranchAddress("ncah", &n_cah);
  myLCTuple->SetBranchAddress("nsch", &n_sch);
  myLCTuple->SetBranchAddress("caene", cah_ene);
  float calhit_ratio;

  // Pandora Cluster
  int n_cl;
  float *cl_x = new float[10000];
  float *cl_y = new float[10000];
  float *cl_z = new float[10000];
  float *cl_ene = new float[10000];
  myLCTuple->SetBranchAddress("nclu", &n_cl);
  myLCTuple->SetBranchAddress("clpox", cl_x);
  myLCTuple->SetBranchAddress("clpoy", cl_y);
  myLCTuple->SetBranchAddress("clpoz", cl_z);
  myLCTuple->SetBranchAddress("clene", cl_ene);

  float cl_distance, cl_eta, cl_phi, calorimeter_ene;

  bool matched_prt = false;

  int n_mcp_stable = 0;
  int num_matched_pions = 0;
  int num_matched_charged = 0;

  // open i/o files
  std::ofstream mcp_file;
  std::ofstream reco_file;
  std::ofstream trst_file;
  std::ofstream eff_file;

  mcp_file.open("Outputs/MCparticles.csv");
  reco_file.open("Outputs/recoparticles.csv");
  trst_file.open("Outputs/trackstates.csv");
  eff_file.open("Outputs/efficiencies.txt");

  mcp_file << "Event,Particle type,Energy,Eta,Phi" << '\n';
  reco_file << "Event,Particle type,Energy,Eta,Phi" << '\n';
  trst_file << "Event,Eta,Phi" << '\n';

  if (!mcp_file.is_open() || !reco_file.is_open() || !trst_file.is_open() || !eff_file.is_open())
  {
    std::cout << "I/O .csv file(s) can't be opened correctly" << '\n';
    return;
  }
  const long int nEntries = (maxEvents < 0 ? myLCTuple->GetEntries() : maxEvents);
  for (int ientry = 0; ientry < nEntries; ++ientry)
  {
    if (ientry % 50 == 0)
      std::cout << ientry << " / " << nEntries << '\n';

    myLCTuple->GetEntry(ientry);

    n_mcp_stable = 0;

    // --- loop over the Monte Carlo particles
    for (int i = 0; i < n_mcp; ++i)
    {
      // --- keep only the stable particles and save the information
      if (mcp_genCode[i] != 1)
        continue;

      pt_mc[n_mcp_stable] = TMath::Sqrt(pow(mcp_px[i], 2) + pow(mcp_py[i], 2));
      pmag_mc[n_mcp_stable] = TMath::Sqrt(pow(mcp_px[i], 2) + pow(mcp_py[i], 2) + pow(mcp_pz[i], 2));
      ene_mc[n_mcp_stable] = mcp_ene[i];
      eta_mc[n_mcp_stable] = 0.5 * TMath::Log((pmag_mc[i] + mcp_pz[i]) / (pmag_mc[i] - mcp_pz[i]));
      phi_mc[n_mcp_stable] = TMath::ATan2(mcp_py[i], mcp_px[i]);
      pdg_mc[n_mcp_stable] = mcp_pdg[i];

      ++n_mcp_stable;
    } // i loop

    for (int i = 0; i < n_mcp_stable; ++i){
      h_mc_pdg->Fill(pdg_mc[i]);
      h_mc_pt->Fill(pt_mc[i]);

      // save data in csv file
      mcp_file << ientry << "," << pdg_mc[i] << "," << ene_mc[i] << "," << eta_mc[i] << "," << phi_mc[i] << '\n';
    }

    // --- loop over the reconstructed particles
    for (int j = 0; j < n_reco; ++j)
    {
      h_reco_pdg->Fill(reco_type[j]);

      pmag_rc[j] = TMath::Sqrt(pow(reco_px[j], 2) + pow(reco_py[j], 2) + pow(reco_pz[j], 2));
      eta_rc[j] = 0.5 * TMath::Log((pmag_rc[j] + reco_pz[j]) / (pmag_rc[j] - reco_pz[j]));
      phi_rc[j] = TMath::ATan2(reco_py[j], reco_px[j]);

      reco_file << ientry << "," << reco_type[j] << "," << reco_ene[j] << "," << eta_rc[j] << "," << phi_rc[j] << '\n';

      if ((reco_ene[j] - ene_mc[0]) < 10.)
      {
        if (TMath::Sqrt(pow((eta_rc[j] - eta_mc[0]), 2) + pow((phi_rc[j] - phi_mc[0]), 2)) < 0.01)
        {
          if (reco_q[j] != 0)
          {
            ++num_matched_charged;
            if (reco_type[j] == 211)
            {
              ++num_matched_pions;
            }
          }
        }
      }

      // from here on, only pions
      if (reco_type[j] == 211)
      {
        h_reco_ene->Fill(reco_ene[j]);
      }

    } // j loop

    // k = 0 IP, k = 1 first tracker hit, k = 2 last tracker hit, k = 3 entrance ecal
    for (int k = 0; k < n_trst; ++k)
    {

      theta_trst[k] = TMath::ATan(1. / trst_tnl[k]);
      if (theta_trst[k] < 0.)
      {
        theta_trst[k] += TMath::Pi();
      }
      eta_trst[k] = -TMath::Log(TMath::Tan(theta_trst[k] / 2.));
      pt_trst[k] = 3.57 * 1.6022e-19 / (fabs(trst_ome[k]) * 1.e3 * 5.344286e-19);
      ptunc_trst[k] = fabs(sqrt(trst_cov[k * 15 + 5]) / trst_ome[k]);
      mom_trst[k] = fabs(pt_trst[k] / sin(theta_trst[k]));
      mom_unc_trst[k] = sqrt(pow(cos(theta_trst[k]) * sqrt(trst_cov[k * 15 + 14]) / (1 + pow(trst_tnl[k], -2)) / sin(theta_trst[k]), 2) + pow(ptunc_trst[k], 2));

      if (k % 4 == 0)
      {
        // filling track states file
        trst_file << ientry << "," << eta_trst[k] << "," << trst_phi[k] << "," << sqrt(trst_cov[k * 15 + 5]) / trst_ome[k] << '\n';
        h_pt->Fill(pt_trst[k]);
        h_pt_unc->Fill(ptunc_trst[k]);

        h_mom_ts->Fill(mom_trst[k]);
        h_mom_unc_ts->Fill(mom_unc_trst[k]);
        h_omega->Fill(trst_ome[k] * 1.e3);
        h_omega_unc->Fill(sqrt(trst_cov[k * 15 + 5]) * 1.e3);

        if (reco_type[0] == 211)
        {
          h_pt_pions->Fill(pt_trst[k]);
          h_pt_unc_pions->Fill(ptunc_trst[k]);

          g_pt_pions->AddPoint(pt_trst[k], ptunc_trst[k]);
        }
        else
        {
          h_pt_nopions->Fill(pt_trst[k]);
          h_pt_unc_nopions->Fill(ptunc_trst[k]);

          g_pt_nopions->AddPoint(pt_trst[k], ptunc_trst[k]);

          /*if (pt_trst[k] > 150.)
          {
            std::cout << "pt: " << pt_trst[k] << ", event: " << ientry << '\n';
          }
          */
        }

        // loop over mc particles to compare the direction
        for (int i = 0; i < n_mcp_stable; ++i)
        {
          delta_eta = eta_trst[k] - eta_mc[i];
          delta_phi = trst_phi[k] - phi_mc[i];
          delta_r = TMath::Sqrt(pow((trst_phi[k] - phi_mc[i]), 2) + pow((eta_trst[k] - eta_mc[i]), 2));

          h_deltaeta_mc->Fill(delta_eta);
          h_deltaphi_mc->Fill(delta_phi);
          h_deltaR_mc->Fill(delta_r);
        }
      }

      if (k % 4 == 3) // compare recoparticle and track at ecal
      {
        // loop over 1st reco particle to compare the direction
        for (int j = 0; j < 1; ++j)
        {
          delta_eta = eta_trst[k] - eta_rc[j];
          delta_phi = trst_phi[k] - phi_rc[j];
          delta_r = TMath::Sqrt(pow((trst_phi[k] - phi_rc[j]), 2) + pow((eta_trst[k] - eta_rc[j]), 2));

          h_deltaphi_reco->Fill(delta_phi);
          h_deltaeta_reco->Fill(delta_eta);
          h_deltaR_reco->Fill(delta_r);

          if (reco_type[j] == 211)
          {
            h_deltaeta_reco_pions->Fill(delta_eta);
            h_deltaphi_reco_pions->Fill(delta_phi);
            h_deltaR_reco_pions->Fill(delta_r);
          }
          else
          {
            h_deltaeta_reco_nopions->Fill(delta_eta);
            h_deltaphi_reco_nopions->Fill(delta_phi);
            h_deltaR_reco_nopions->Fill(delta_r);
          }
        }
      }

    } // k loop

    h_calohits->Fill(n_cah);
    h_simcalohits->Fill(n_sch);

    calhit_ratio = (float)n_cah / (float)n_sch;

    reco_type[0] == 211 ? h_calohitsratio_pions->Fill(calhit_ratio) : h_calohitsratio_nopions->Fill(calhit_ratio);

    calorimeter_ene = std::accumulate(cah_ene, cah_ene + n_cah, 0.);

    if (reco_q[0] != 0)
    {
      g_ptvsene_ch->AddPoint(pt_trst[0], calorimeter_ene);
    }

    else
    {
      g_ptvsene_unch->AddPoint(pt_trst[0], calorimeter_ene);
    }

    {
      if (n_reco != 1)
        continue;

      /*
      cah_x = std::inner_product(cl_x, cl_x + n_cl, cl_ene, 0.) / calorimeter_ene;
      cah_y = std::inner_product(cl_y, cl_y + n_cl, cl_ene, 0.) / calorimeter_ene;
      cah_z = std::inner_product(cl_z, cl_z + n_cl, cl_ene, 0.) / calorimeter_ene;
      */

      cl_distance = TMath::Sqrt(pow(cl_x[0], 2) + pow(cl_y[0], 2) + pow(cl_z[0], 2));
      cl_eta = 0.5 * TMath::Log((cl_distance + cl_z[0]) / (cl_distance - cl_z[0]));
      cl_phi = TMath::ATan2(cl_y[0], cl_x[0]);

      delta_r = TMath::Sqrt(pow((trst_phi[3] - cl_phi), 2) + pow((eta_trst[3] - cl_eta), 2));
      h_deltaR_trst_cl->Fill(delta_r);

      if (reco_q[0] != 0)
      {
        h_deltaR_trst_cl_ch->Fill(delta_r);
      }

      else
      {
        h_deltaR_trst_cl_unch->Fill(delta_r);
      }
    }

  } // ientry loop

  eff_file << "efficiency (pions) = " << num_matched_pions << "/" << nEntries << " = " << 
    std::setprecision(3) << float(num_matched_pions) / float(nEntries) * 100 << "%" << '\n';

  eff_file << "efficiency (charged particles) = " << num_matched_charged << "/" << nEntries << " = " << 
    std::setprecision(3) << float(num_matched_charged) / float(nEntries) * 100 << "%" << '\n';

  // close i/o files
  mcp_file.close();
  reco_file.close();
  trst_file.close();
  eff_file.close();

  // draw the histograms
  {
    TCanvas *c0 = new TCanvas("cPDG", "cPDG", 800, 800);
    c0->Divide(2, 1);
    c0->cd(1);
    h_mc_pdg->Draw();
    c0->cd(2);
    h_reco_pdg->Draw();
    c0->SaveAs("./plots/cPDG.png");
    c0->Close();

    TCanvas *c1 = new TCanvas("c_mc_pt", "c_mc_pt", 800, 800);
    h_mc_pt->Draw();
    c1->SaveAs("./plots/c_mc_pt.png");
    c1->Close();

    TCanvas *c2 = new TCanvas("c_reco", "c_reco", 800, 800);
    gStyle->SetOptStat(111111);
    h_reco_ene->Draw();
    c2->SaveAs("./plots/c_reco.png");
    c2->Close();

    TCanvas *c3 = new TCanvas("c_deltaR_mc", "c_deltaR_mc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_mc->Draw();
    c3->SaveAs("./plots/c_deltaR_mc.png");
    c3->Close();

    TCanvas *c4 = new TCanvas("c_deltaR_reco", "c_deltaR_reco", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_reco->Draw();
    c4->SaveAs("./plots/c_deltaR_reco.png");
    c4->Close();

    TCanvas *c5 = new TCanvas("c_deltaR_reco_pions", "c_deltaR_reco_pions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_reco_pions->Draw();
    c5->SaveAs("./plots/c_deltaR_reco_pions.png");
    c5->Close();

    TCanvas *c6 = new TCanvas("c_deltaR_reco_nopions", "c_deltaR_reco_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_reco_nopions->Draw();
    c6->SaveAs("./plots/c_deltaR_reco_nopions.png");
    c6->Close();

    TCanvas *c7 = new TCanvas("c_deltaeta_mc", "c_deltaeta_mc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaeta_mc->Draw();
    c7->SaveAs("./plots/c_deltaeta_mc.png");
    c7->Close();

    TCanvas *c8 = new TCanvas("c_deltaeta_reco", "c_deltaeta_reco", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaeta_reco->Draw();
    c8->SaveAs("./plots/c_deltaeta_reco.png");
    c8->Close();

    TCanvas *c9 = new TCanvas("c_deltaeta_reco_pions", "c_deltaeta_reco_pions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaeta_reco_pions->Draw();
    c9->SaveAs("./plots/c_deltaeta_reco_pions.png");
    c9->Close();

    TCanvas *c10 = new TCanvas("c_deltaeta_reco_nopions", "c_deltaeta_reco_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaeta_reco_nopions->Draw();
    c10->SaveAs("./plots/c_deltaeta_reco_nopions.png");
    c10->Close();

    TCanvas *c11 = new TCanvas("c_deltaphi_mc", "c_deltaphi_mc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaphi_mc->Draw();
    c11->SaveAs("./plots/c_deltaphi_mc.png");
    c11->Close();

    TCanvas *c12 = new TCanvas("c_deltaphi_reco", "c_deltaphi_reco", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaphi_reco->Draw();
    c12->SaveAs("./plots/c_deltaphi_reco.png");
    c12->Close();

    TCanvas *c13 = new TCanvas("c_deltaphi_reco_pions", "c_deltaphi_reco_pions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaphi_reco_pions->Draw();
    c13->SaveAs("./plots/c_deltaphi_reco_pions.png");
    c13->Close();

    TCanvas *c14 = new TCanvas("c_deltaphi_reco_nopions", "c_deltaphi_reco_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaphi_reco_nopions->Draw();
    c14->SaveAs("./plots/c_deltaphi_reco_nopions.png");
    c14->Close();

    TCanvas *c15 = new TCanvas("c_calohits", "c_calohits", 800, 800);
    gStyle->SetOptStat(111111);
    h_calohits->Draw();
    c15->SaveAs("./plots/c_calohits.png");
    c15->Close();

    TCanvas *c16 = new TCanvas("c_simcalohits", "c_simcalohits", 800, 800);
    gStyle->SetOptStat(111111);
    h_simcalohits->Draw();
    c16->SaveAs("./plots/c_simcalohits.png");
    c16->Close();

    TCanvas *c17 = new TCanvas("c_calohitsratio_pions", "c_simcalohitsratio_pions", 800, 800);
    gStyle->SetOptStat(111111);
    h_calohitsratio_pions->Draw();
    c17->SaveAs("./plots/c_calohitsratio_pions.png");
    c17->Close();

    TCanvas *c18 = new TCanvas("c_calohitsratio_nopions", "c_simcalohitsratio_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    h_calohitsratio_nopions->Draw();
    c18->SaveAs("./plots/c_calohitsratio_nopions.png");
    c18->Close();

    TCanvas *c19 = new TCanvas("c_pt", "c_pt", 800, 800);
    gStyle->SetOptStat(111111);
    h_pt->Draw();
    c19->SaveAs("./plots/c_pt.png");
    c19->Close();

    TCanvas *c20 = new TCanvas("c_pt_pions", "c_pt_pions", 800, 800);
    gStyle->SetOptStat(111111);
    h_pt_pions->Draw();
    c20->SaveAs("./plots/c_pt_pions.png");
    c20->Close();

    TCanvas *c21 = new TCanvas("c_pt_nopions", "c_pt_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    h_pt_nopions->Draw();
    c21->SaveAs("./plots/c_pt_nopions.png");
    c21->Close();

    TCanvas *c22 = new TCanvas("c_pt_unc", "c_pt_unc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_pt_unc->Draw();
    c22->SaveAs("./plots/c_pt_unc.png");
    c22->Close();

    TCanvas *c23 = new TCanvas("c_pt_unc_pions", "c_pt_unc_pions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_pt_unc_pions->Draw();
    c23->SaveAs("./plots/c_pt_unc_pions.png");
    c23->Close();

    TCanvas *c24 = new TCanvas("c_pt_unc_nopions", "c_pt_unc_nopions", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_pt_unc_nopions->Draw();
    c24->SaveAs("./plots/c_pt_unc_nopions.png");
    c24->Close();

    TCanvas *c25 = new TCanvas("c_ome", "c_ome", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_omega->Draw();
    c25->SaveAs("./plots/c_omega.png");
    c25->Close();

    TCanvas *c26 = new TCanvas("c_ome_unc", "c_ome_unc", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_omega_unc->Draw();
    c26->SaveAs("./plots/c_omega_unc.png");
    c26->Close();

    TCanvas *c27 = new TCanvas("c_mom_ts", "c_mom_ts", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_mom_ts->Draw();
    c27->SaveAs("./plots/c_mom_ts.png");
    c27->Close();

    TCanvas *c28 = new TCanvas("c_mom_unc_ts", "c_mom_unc_ts", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_mom_unc_ts->Draw();
    c28->SaveAs("./plots/c_mom_unc_ts.png");
    c28->Close();

    TCanvas *c29 = new TCanvas("c_mom_vs_unc_ts", "c_mom_unc_ts", 1100, 500);
    TGaxis::SetMaxDigits(4);
    c29->Divide(2, 1);
    c29->cd(1);
    gStyle->SetOptStat(111111);
    g_pt_pions->SetTitle("p_{T} tracks (pions);p_{T} [GeV];#sigma_{p_{T}}/^{}p_{T}");
    g_pt_pions->GetYaxis()->SetTitleOffset(0.7);
    g_pt_pions->SetMarkerStyle(26);
    g_pt_pions->SetMarkerSize(0.3);
    g_pt_pions->SetMarkerColor(kBlue);
    gPad->SetLogy();
    // g_pt_pions->GetXaxis()->SetLimits();
    g_pt_pions->SetMaximum(0.5);
    g_pt_pions->Draw("AP");
    c29->cd(2);
    gStyle->SetOptStat(111111);
    g_pt_nopions->SetTitle("p_{T} tracks (no pions);p_{T} [GeV];#sigma_{p_{T}}/^{}p_{T}");
    g_pt_nopions->GetYaxis()->SetTitleOffset(0.7);
    g_pt_nopions->SetMarkerStyle(26);
    g_pt_nopions->SetMarkerSize(0.3);
    g_pt_nopions->SetMarkerColor(kRed);
    gPad->SetLogy();
    // g_pt_nopions->GetXaxis()->SetLimits();
    g_pt_nopions->SetMaximum(0.5);
    g_pt_nopions->Draw("AP");
    c29->SaveAs("./plots/c_scatter.pdf");
    c29->Close();

    TCanvas *c30 = new TCanvas("c_deltaR_trst_cl", "c_deltaR_trst_cl", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_trst_cl->Draw();
    c30->SaveAs("./plots/c_deltaR_trst_cl.png");
    c30->Close();

    TCanvas *c31 = new TCanvas("c_deltaR_trst_cl_ch", "c_deltaR_trst_cl_ch", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_trst_cl_ch->Draw();
    c31->SaveAs("./plots/c_deltaR_trst_cl_ch.png");
    c31->Close();

    TCanvas *c32 = new TCanvas("c_deltaR_trst_cl_unch", "c_deltaR_trst_cl_unch", 800, 800);
    gStyle->SetOptStat(111111);
    TGaxis::SetMaxDigits(3);
    h_deltaR_trst_cl_unch->Draw();
    c32->SaveAs("./plots/c_deltaR_trst_cl_unch.png");
    c32->Close();

    TCanvas *c33 = new TCanvas("c_ptvsene", "c_ptvsene", 1100, 500);
    TGaxis::SetMaxDigits(4);
    c33->Divide(2, 1);
    c33->cd(1);
    gStyle->SetOptStat(111111);
    g_ptvsene_ch->SetTitle("p_{T} tracks vs clusters energies (charged prts);p_{T} [GeV]; Energy [GeV]");
    g_ptvsene_ch->GetYaxis()->SetTitleOffset(1.2);
    g_ptvsene_ch->SetMarkerStyle(26);
    g_ptvsene_ch->SetMarkerSize(0.3);
    g_ptvsene_ch->SetMarkerColor(kBlue);
    // g_ptvsene_ch->GetXaxis()->SetLimits(94., 120.);
    // g_ptvsene_unch->SetMaximum(130.);
    // g_ptvsene_unch->SetMinimum(40.);
    g_ptvsene_ch->Draw("AP");
    c33->cd(2);
    gStyle->SetOptStat(111111);
    g_ptvsene_unch->SetTitle("p_{T} tracks vs clusters energies (neutrons);p_{T} [GeV]; Energy [GeV]");
    g_ptvsene_unch->GetYaxis()->SetTitleOffset(1.2);
    g_ptvsene_unch->SetMarkerStyle(26);
    g_ptvsene_unch->SetMarkerSize(0.3);
    g_ptvsene_unch->SetMarkerColor(kRed);
    // g_ptvsene_unch->GetXaxis()->SetLimits(94., 120.);
    // g_ptvsene_unch->SetMaximum(130.);
    // g_ptvsene_unch->SetMinimum(40.);
    g_ptvsene_unch->Draw("AP");
    c33->SaveAs("./plots/c_ptvsene.pdf");
    c33->Close();
  }

  // Clean up the heap
  {
    delete[] mcp_pdg;
    delete[] mcp_genCode;
    delete[] mcp_siode;
    delete[] mcp_vx;
    delete[] mcp_vy;
    delete[] mcp_vz;
    delete[] mcp_px;
    delete[] mcp_py;
    delete[] mcp_pz;
    delete[] mcp_q;
    delete[] mcp_t;

    delete[] reco_type;
    delete[] reco_px;
    delete[] reco_py;
    delete[] reco_pz;
    delete[] reco_ene;
    delete[] reco_q;

    delete[] trst_phi;
    delete[] trst_tnl;
    delete[] trst_ome;

    delete[] trst_rpx;
    delete[] trst_rpy;
    delete[] trst_rpz;

    delete[] pt_mc;
    delete[] pmag_mc;
    delete[] eta_mc;
    delete[] phi_mc;
    delete[] ene_mc;
    delete[] pdg_mc;
    delete[] eta_rc;
    delete[] phi_rc;
    delete[] theta_trst;
    delete[] eta_trst;
    delete[] pt_trst;
    delete[] ptunc_trst;
    delete[] mom_trst;

    delete[] cl_x;
    delete[] cl_y;
    delete[] cl_z;
    delete[] cl_ene;

    delete[] cah_ene;

    h_mc_pdg->Delete();
    h_mc_pt->Delete();
    h_reco_ene->Delete();
    h_reco_pdg->Delete();
    h_deltaeta_mc->Delete();
    h_deltaeta_reco->Delete();
    h_deltaeta_reco_nopions->Delete();
    h_deltaeta_reco_pions->Delete();
    h_deltaphi_mc->Delete();
    h_deltaphi_reco->Delete();
    h_deltaphi_reco_nopions->Delete();
    h_deltaphi_reco_pions->Delete();
    h_deltaR_mc->Delete();
    h_deltaR_reco->Delete();
    h_deltaR_reco_nopions->Delete();
    h_deltaR_reco_pions->Delete();
    h_pt->Delete();
    h_pt_nopions->Delete();
    h_pt_pions->Delete();
    h_pt_unc->Delete();
    h_pt_unc_nopions->Delete();
    h_pt_unc_pions->Delete();
    h_mom_ts->Delete();
    h_mom_unc_ts->Delete();
    h_omega->Delete();
    h_omega_unc->Delete();
    h_calohits->Delete();
    h_simcalohits->Delete();
    h_calohitsratio_pions->Delete();
    h_calohitsratio_nopions->Delete();
    h_deltaR_trst_cl->Delete();
    h_deltaR_trst_cl_ch->Delete();
    h_deltaR_trst_cl_unch->Delete();

    g_pt_pions->Delete();
    g_pt_nopions->Delete();

    g_ptvsene_ch->Delete();
    g_ptvsene_unch->Delete();
  }

  // --- Close the ntuple file
  input_file->Close();
  input_file->Delete();
}
