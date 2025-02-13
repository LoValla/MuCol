#include <fstream>
#include <iostream>

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

// ===========================================================================

void update_tree(const TString filename = "lctuple_gun.root")
{
    //  --- Open the ntuple file and get the tree

    TFile *input_file = new TFile(filename, "update");

    TTree *myLCTuple = (TTree *)input_file->Get("MyLCTuple");

    // --- MC particles
    Int_t n_mcp;
    int *mcp_pdg = new int[1500000];
    int *mcp_genCode = new int[1500000];
    float *mcp_px = new float[1500000];
    float *mcp_py = new float[1500000];
    float *mcp_pz = new float[1500000];
    float *mcp_ene = new float[1500000];

    myLCTuple->SetBranchAddress("nmcp", &n_mcp);
    myLCTuple->SetBranchAddress("mcpdg", mcp_pdg);
    myLCTuple->SetBranchAddress("mcgst", mcp_genCode);
    myLCTuple->SetBranchAddress("mcmox", mcp_px);
    myLCTuple->SetBranchAddress("mcmoy", mcp_py);
    myLCTuple->SetBranchAddress("mcmoz", mcp_pz);
    myLCTuple->SetBranchAddress("mcene", mcp_ene);

    // --- RECO particles
    Int_t n_reco;
    int *reco_type = new int[100000];
    float *reco_px = new float[100000];
    float *reco_py = new float[100000];
    float *reco_pz = new float[100000];
    float *reco_ene = new float[100000];
    float *reco_q = new float[100000];

    myLCTuple->SetBranchAddress("nrec", &n_reco);
    myLCTuple->SetBranchAddress("rctyp", reco_type);
    myLCTuple->SetBranchAddress("rcmox", reco_px);
    myLCTuple->SetBranchAddress("rcmoy", reco_py);
    myLCTuple->SetBranchAddress("rcmoz", reco_pz);
    myLCTuple->SetBranchAddress("rcene", reco_ene);
    myLCTuple->SetBranchAddress("rccha", reco_q);

    // trackstates
    Int_t n_ts;
    float *ts_phi = new float[100000];
    float *ts_tnl = new float[100000];
    float *ts_ome = new float[100000];
    float *ts_cov = new float[1000000];

    myLCTuple->SetBranchAddress("ntrst", &n_ts);
    myLCTuple->SetBranchAddress("tsphi", ts_phi);
    myLCTuple->SetBranchAddress("tstnl", ts_tnl);
    myLCTuple->SetBranchAddress("tsome", ts_ome);
    myLCTuple->SetBranchAddress("tscov", ts_cov);

    // declare variables and tbranches;
    Float_t *pmom_mc = new Float_t[1000];
    Float_t *eta_mc = new Float_t[1000];
    Float_t *phi_mc = new Float_t[1000];
    Float_t *ene_mc = new Float_t[1000];
    Int_t nmcp_usr;
    myLCTuple->Branch("nmcp_usr", &nmcp_usr, "nmcp_usr/I");
    TBranch *b_mceta = myLCTuple->Branch("mceta_s", &eta_mc, "mceta_s[nmcp_usr]/F");
    // TBranch *b_mctheta = myLCTuple->Branch("mctheta", &theta_mc, "mctheta/F");
    TBranch *b_mcphi = myLCTuple->Branch("mcphi_s", &phi_mc, "mcphi_s[nmcp_usr]/F");
    TBranch *b_mcmom = myLCTuple->Branch("mcmom_s", &pmom_mc, "mcmom_s[nmcp_usr]/F");
    TBranch *b_mcene = myLCTuple->Branch("mcene_s", &ene_mc, "mcene_s[nmcp_usr]/F");

    Float_t *pmom_rc = new Float_t[1000];
    Float_t *eta_rc = new Float_t[1000];
    Float_t *phi_rc = new Float_t[1000];
    Int_t nrec_usr;
    myLCTuple->Branch("nrec_usr", &nrec_usr, "nrec_usr/I");
    TBranch *b_rceta = myLCTuple->Branch("rceta", &eta_rc, "rceta[nrec_usr]/F");
    // TBranch *b_rctheta = myLCTuple->Branch("rctheta", &theta_rc, "rctheta/F");
    TBranch *b_rcphi = myLCTuple->Branch("rcphi", &phi_rc, "rcphi[nrec_usr]/F");
    TBranch *b_rcmom = myLCTuple->Branch("rcmom", &pmom_rc, "rcmom[nrec_usr]/F");

    Float_t *theta_ts = new Float_t[1000];
    Float_t *eta_ts = new Float_t[1000];
    Float_t *pt_ts = new Float_t[1000];
    Float_t *ptunc_ts = new Float_t[1000];
    Int_t ntrst_usr;

    myLCTuple->Branch("ntrst_usr", &ntrst_usr, "ntrst_usr/I");

    TBranch *b_tspt = myLCTuple->Branch("tspt", &pt_ts, "tspt[ntrst_usr]/F");
    TBranch *b_tsptunc = myLCTuple->Branch("tsptunc", &ptunc_ts, "tsptunc[ntrst_usr]/F");
    TBranch *b_tstheta = myLCTuple->Branch("tstheta", &theta_ts, "tstheta[ntrst_usr]/F");
    TBranch *b_tseta = myLCTuple->Branch("tseta", &eta_ts, "tseta[ntrst_usr]/F");

    const long int nEntries = myLCTuple->GetEntries();

    for (int ientry = 0; ientry < nEntries; ++ientry)
    {
        /*if (ientry % 50 == 0)
        {
            std::cout << ientry << " / " << nEntries << '\n';
        }*/

        myLCTuple->GetEntry(ientry);

        nrec_usr = n_reco;
        ntrst_usr = n_ts;
        // --- loop over the Monte Carlo particles
        for (int i = 0; i < n_mcp; ++i)
        {
            // --- keep only the stable particles
            if (mcp_genCode[i] != 1)
                continue;

            pmom_mc[i] = TMath::Sqrt(pow(mcp_px[i], 2) + pow(mcp_py[i], 2) + pow(mcp_pz[i], 2));
            eta_mc[i] = 0.5 * TMath::Log((pmom_mc[i] + mcp_pz[i]) / (pmom_mc[i] - mcp_pz[i]));
            phi_mc[i] = TMath::ATan2(mcp_py[i], mcp_px[i]);
            ene_mc[i] = mcp_ene[i];
            ++nmcp_usr;

        } // i loop

        // --- loop over the reconstructed particles
        for (int j = 0; j < n_reco; ++j)
        {
            pmom_rc[j] = TMath::Sqrt(pow(reco_px[j], 2) + pow(reco_py[j], 2) + pow(reco_pz[j], 2));
            eta_rc[j] = 0.5 * TMath::Log((pmom_rc[j] + reco_pz[j]) / (pmom_rc[j] - reco_pz[j]));
            phi_rc[j] = TMath::ATan2(reco_py[j], reco_px[j]);

        } // j loop

        // --- loop over the track states (choose the reference point k = n, n = 0,1,2,3) and covariances
        // k = 0 IP, k = 1 first tracker hit, k = 2 last tracker hit, k = 3 entrance ecal
        for (int k = 0; k < n_ts; ++k)
        {

            theta_ts[k] = TMath::ATan(1. / ts_tnl[k]);

            if (theta_ts[k] < 0.)
            {
                theta_ts[k] += TMath::Pi();
            }

            eta_ts[k] = -TMath::Log(TMath::Tan(theta_ts[k] / 2.));
            pt_ts[k] = 3.57 * 1.6022e-19 / (ts_ome[k] * 1.e3 * 5.344286e-19);
            ptunc_ts[k] = sqrt(ts_cov[k * 15 + 5]) / ts_ome[k];

        } // k loop
    } // ientry loop

    // print and write the branches
    // myLCTuple->Print();
    myLCTuple->SetEntries();
    myLCTuple->Write();

    std::cout << "The TTree was updated with the high-level variables" << '\n';

    delete[] mcp_pdg;
    delete[] mcp_genCode;
    delete[] mcp_px;
    delete[] mcp_py;
    delete[] mcp_pz;
    delete[] mcp_ene;

    delete[] reco_type;
    delete[] reco_px;
    delete[] reco_py;
    delete[] reco_pz;
    delete[] reco_ene;

    delete[] ts_phi;
    delete[] ts_tnl;
    delete[] ts_ome;

    delete[] theta_ts;
    delete[] eta_ts;
    delete[] pt_ts;
    delete[] ptunc_ts;

    delete[] eta_rc;
    delete[] phi_rc;
    delete[] pmom_rc;

    delete[] eta_mc;
    delete[] phi_mc;
    delete[] pmom_mc;
    delete[] ene_mc;

    // --- Close the ntuple file

    input_file->Close();

    delete input_file;
}
