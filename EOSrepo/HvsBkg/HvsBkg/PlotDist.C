#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"

using namespace RooFit;

void display_underoverflow(TH1 *h) {
  h->GetXaxis()->SetRange(0, h->GetNbinsX() + 1);
}
void display_underflow(TH1 *h) { h->GetXaxis()->SetRange(0, h->GetNbinsX()); }
void display_overflow(TH1 *h) {
  h->GetXaxis()->SetRange(1, h->GetNbinsX() + 1);
}

void fit_plotter(const TString &filename) {

  gStyle->SetTextFont(42);

  TRandom3 *ran = new TRandom3(0);

  TFile *histo_file = new TFile(filename, "READ");
  if (!histo_file->IsOpen())
    throw std::invalid_argument("histograms filename not valid");

  TDirectory *dir = (TDirectory *)histo_file->Get("dataset/InputVariables_Id");

  TH1F *h_sig_MInv = static_cast<TH1F *>(dir->Get("recMInv__Signal_Id"));
  TH1F *h_bkg_MInv = static_cast<TH1F *>(dir->Get("recMInv__Background_Id"));
  TH1F *h_sig_PtPair = static_cast<TH1F *>(dir->Get("recPtPair__Signal_Id"));
  TH1F *h_bkg_PtPair =
      static_cast<TH1F *>(dir->Get("recPtPair__Background_Id"));
  TH1F *h_sig_EtaPair = static_cast<TH1F *>(dir->Get("recEtaPair__Signal_Id"));
  TH1F *h_bkg_EtaPair =
      static_cast<TH1F *>(dir->Get("recDeltaPhi__Background_Id"));
  TH1F *h_sig_DeltaPhi =
      static_cast<TH1F *>(dir->Get("recDeltaPhi__Signal_Id"));
  TH1F *h_bkg_DeltaPhi =
      static_cast<TH1F *>(dir->Get("recDeltaPhi__Background_Id"));
  TH1F *h_sig_SepAngle =
      static_cast<TH1F *>(dir->Get("recSepAngle__Signal_Id"));
  TH1F *h_bkg_SepAngle =
      static_cast<TH1F *>(dir->Get("recSepAngle__Background_Id"));
  TH1F *h_sig_PtSingle_0 =
      static_cast<TH1F *>(dir->Get("recPtSingle_0__Signal_Id"));
  TH1F *h_bkg_PtSingle_0 =
      static_cast<TH1F *>(dir->Get("recPtSingle_0__Background_Id"));
  TH1F *h_sig_PtSingle_1 =
      static_cast<TH1F *>(dir->Get("recPtSingle_1__Signal_Id"));
  TH1F *h_bkg_PtSingle_1 =
      static_cast<TH1F *>(dir->Get("recPtSingle_1__Background_Id"));
  TH1F *h_sig_EtaSingle_0 =
      static_cast<TH1F *>(dir->Get("recEtaSingle_0__Signal_Id"));
  TH1F *h_bkg_EtaSingle_0 =
      static_cast<TH1F *>(dir->Get("recEtaSingle_0__Background_Id"));
  TH1F *h_sig_EtaSingle_1 =
      static_cast<TH1F *>(dir->Get("recEtaSingle_1__Signal_Id"));
  TH1F *h_bkg_EtaSingle_1 =
      static_cast<TH1F *>(dir->Get("recEtaSingle_1__Background_Id"));

  std::vector<TH1F *> h_sig = {
      h_sig_MInv,       h_sig_PtPair,      h_sig_EtaPair,
      h_sig_DeltaPhi,   h_sig_SepAngle,    h_sig_PtSingle_0,
      h_sig_PtSingle_1, h_sig_EtaSingle_0, h_sig_EtaSingle_1};
  std::vector<TH1F *> h_bkg = {
      h_bkg_MInv,       h_bkg_PtPair,      h_bkg_EtaPair,
      h_bkg_DeltaPhi,   h_bkg_SepAngle,    h_bkg_PtSingle_0,
      h_bkg_PtSingle_1, h_bkg_EtaSingle_0, h_bkg_EtaSingle_1};
    
    //std::cout << "h_sig[i] = " << h_sig[i] << '\n';
    
    for(int unsigned i = 0; i != h_sig.size(); ++i){
        auto c_recMInv = new TCanvas("c_recMInv", "c_recMInv", 1000, 800);
        h_sig[i]->Scale(1 / h_sig[i]->Integral());
        h_bkg[i]->Scale(1 / h_bkg[i]->Integral());
        h_sig[i]->SetLineColor(kRed);
        h_sig[i]->SetFillColor(kRed);
        h_sig[i]->SetFillStyle(3002);
        h_bkg[i]->SetLineColor(kBlue);
        h_bkg[i]->SetFillColor(kBlue);
        h_bkg[i]->SetFillStyle(3003);
        display_underoverflow(h_sig[i]);
        display_underoverflow(h_bkg[i]);
        h_sig[i]->GetYaxis()->SetRangeUser(0., h_bkg[i]->GetMaximum() * 1.15);
        h_sig[i]->Draw("HISTO");
        h_bkg[i]->Draw("HISTO, SAME");
        h_sig[i]->GetXaxis()->SetLabelFont(42);
        h_sig[i]->GetXaxis()->SetLabelSize(0.035);
        h_sig[i]->GetXaxis()->SetTitleFont(42);
        h_sig[i]->GetXaxis()->SetTitleOffset(1.15);
        h_sig[i]->GetXaxis()->SetTitleSize(0.04);
        h_sig[i]->GetYaxis()->SetLabelFont(42);
        h_sig[i]->GetYaxis()->SetLabelSize(0.035);
        h_sig[i]->GetYaxis()->SetTitleFont(42);
        h_sig[i]->GetYaxis()->SetTitleOffset(1.38);
        h_sig[i]->GetYaxis()->SetDecimals();
        c_recMInv->SetLeftMargin(0.12);
        c_recMInv->SetBottomMargin(0.15);
        h_sig[i]->GetYaxis()->SetTitle("Density");
        auto leg0 = new TLegend(0.7, 0.78, 0.88, 0.88);
        gStyle->SetOptStat(0);
        leg0->AddEntry(h_sig[i], "H #rightarrow #tau#tau", "f");
        leg0->AddEntry(h_bkg[i], "Z #rightarrow #tau#tau", "f");
        leg0->SetTextSize(0.037);
        leg0->SetTextFont(42);
        leg0->SetLineColor(kWhite);
        leg0->Draw();
        // gPad->Modified();
        /*TPaveText *pt1 = new TPaveText(120, 0.14, 160, 0.20);
        pt1->SetFillColor(0);
        pt1->SetMargin(0);
        pt1->SetFillStyle(0);
        pt1->SetBorderSize(0);
        pt1->SetTextAlign(11);
        // pt1->SetTextFont(50);
        pt1->AddText("Muon Collider");
        pt1->AddText("#bf{#it{Simulation}}");
        pt1->AddText("#bf{#sqrt{s} = 3 TeV}");
        pt1->AddText("  ");
        pt1->AddText("#bf{p_{T}(#tau) > 10 GeV/c}");
        pt1->AddText("#bf{|#eta| < 2.5}");
        pt1->Draw();
        gStyle->SetOptTitle(0);
        // gROOT->SetStyle("Plain");
        // gROOT->ForceStyle();
        gPad->Update();*/
        c_recMInv->SaveAs(Form("./plotsDist/c.%d.png", i));
        c_recMInv->Close();
        c_recMInv->Delete();
        //pt1->Delete();
        leg0->Delete();
    }

  return;
}

void PlotDist() {
  const TString filename = "/eos/user/l/lvalla/MuColl/HvsZ/TMVACTot.root";
  fit_plotter(filename);
}
