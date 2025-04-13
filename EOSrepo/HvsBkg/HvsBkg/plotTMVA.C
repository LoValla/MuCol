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

void fit_plotter(const TString &filename) {

  gStyle->SetTextFont(42);

  TRandom3 *ran = new TRandom3();

  TFile *histo_file = new TFile(filename, "READ");
  if (!histo_file->IsOpen())
    throw std::invalid_argument("histograms filename not valid");

  TDirectory *dir = (TDirectory *)histo_file->Get("dataset/Method_BDT/BDT");

  TH1F *h_H_test = static_cast<TH1F *>(dir->Get("MVA_BDT_S"));
  TH1F *h_Z_test = static_cast<TH1F *>(dir->Get("MVA_BDT_B"));
  TH1F *h_H_train = static_cast<TH1F *>(dir->Get("MVA_BDT_Train_S"));
  TH1F *h_Z_train = static_cast<TH1F *>(dir->Get("MVA_BDT_Train_B"));

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  gStyle->SetOptStat(0000);
  gStyle->SetStripDecimals(0);

  // Style settings
  h_H_test->SetLineColor(kRed);
  h_Z_test->SetLineColor(kBlue);
  h_H_test->SetLineColor(kRed);
  h_Z_test->SetLineColor(kBlue);

  h_H_test->SetTitle(";BDT output; Events/bin");
  h_H_test->SetMaximum(7);

  // Draw histograms with "axis" option first
  h_H_test->Draw("HIST");
  h_Z_test->Draw("HIST same");
  h_H_train->Draw("E1 same");
  h_Z_train->Draw("E1 same");

  //h_H_test->Scale(1 / h_H_test->Integral());
  //h_H_test->Scale(1 / h_H_test->Integral());

  // Set styles for the hists
  h_H_test->SetLineColor(kRed);
  h_H_test->SetLineWidth(2);
  h_H_test->SetMarkerColor(kRed);
  // h_H_test->SetMarkerStyle(20);
  // h_H_test->SetMarkerSize(0.7);
  h_Z_test->SetLineColor(kBlue);
  h_Z_test->SetLineWidth(2);
  h_Z_test->SetMarkerColor(kBlue);
  // h_Z_test->SetMarkerStyle(20);
  // h_Z_test->SetMarkerSize(0.7);
  h_H_train->SetLineColor(kRed);
  h_H_train->SetMarkerColor(kRed);
  h_H_train->SetMarkerStyle(20);
  h_H_train->SetMarkerSize(1);
  h_Z_train->SetLineColor(kBlue);
  h_Z_train->SetMarkerColor(kBlue);
  h_Z_train->SetMarkerStyle(20);
  h_Z_train->SetMarkerSize(1);

  // Customize the axes
  h_H_test->GetXaxis()->SetTitleOffset(1.3);
  h_H_test->GetXaxis()->SetTitleSize(0.045);
  h_H_test->GetXaxis()->SetLabelSize(0.045);
  h_H_test->GetXaxis()->SetLabelOffset(0.012);
  h_H_test->GetYaxis()->SetTitleSize(0.045);
  h_H_test->GetYaxis()->SetTitleOffset(1.2);
  h_H_test->GetYaxis()->SetLabelSize(0.045);
  h_H_test->GetYaxis()->SetLabelOffset(0.0125);

  c1->SetLeftMargin(0.13); // Increase left margin for better space
  c1->SetBottomMargin(0.15);
  c1->SetRightMargin(0.04);

  // Create and draw the legend
  TLegend *leg = new TLegend(0.15, 0.7, 0.4, 0.88); // Top-right corner
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);   // Remove fill color
  leg->SetTextSize(0.04); // Increase text size
  leg->AddEntry(h_H_test, "Signal (test)", "l");
  leg->AddEntry(h_Z_test, "Background (test)", "l");
  leg->AddEntry(h_H_train, "Signal (train)", "lep");
  leg->AddEntry(h_Z_train, "Background (train)", "lep");
  leg->Draw();

  c1->SaveAs("plotTMVA.pdf");
  c1->Close();

  return;
}

void plotTMVA() {
  const TString filename = "/eos/user/l/lvalla/MuColl/HvsZ/TMVACTot.root";
  fit_plotter(filename);
}
