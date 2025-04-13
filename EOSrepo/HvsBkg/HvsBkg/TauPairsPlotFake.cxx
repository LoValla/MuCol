#include <fstream>
#include <iostream>

#include "Math/Vector4D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

void display_underoverflow(TH1 *h) {
  h->GetXaxis()->SetRange(0, h->GetNbinsX() + 1);
}
void display_underflow(TH1 *h) { h->GetXaxis()->SetRange(0, h->GetNbinsX()); }
void display_overflow(TH1 *h) {
  h->GetXaxis()->SetRange(1, h->GetNbinsX() + 1);
}

void HvsZ_plotter(const TString &filename, const TString &decay0 = "all",
                  const TString &decay1 = "all", const TString &prefix = "") {
  // Open and read energy correction file
  TFile *histo_file = new TFile(filename, "READ");

  if (!histo_file->IsOpen())
    throw std::invalid_argument("histograms filename not valid");

  TList *Zhisto_list =
      (TList *)histo_file->Get("Z_" + decay0 + "+" + decay1 + "_histograms");
  TList *Hhisto_list =
      (TList *)histo_file->Get("H_" + decay0 + "+" + decay1 + "_histograms");

  TH1F *Z_h_MInv = static_cast<TH1F *>(Zhisto_list->FindObject("h_recMInv"));
  TH1F *H_h_MInv = static_cast<TH1F *>(Hhisto_list->FindObject("h_recMInv"));
  TH1F *Z_h_DeltaPhi =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recDeltaPhi"));
  TH1F *H_h_DeltaPhi =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recDeltaPhi"));
  TH1F *Z_h_DeltaR =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recDeltaR"));
  TH1F *H_h_DeltaR =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recDeltaR"));
TH1F *Z_h_DeltaEta =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recDeltaEta"));
  TH1F *H_h_DeltaEta =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recDeltaEta"));
  TH1F *Z_h_PtPair =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recPtPair"));
  TH1F *H_h_PtPair =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recPtPair"));
  TH1F *Z_h_EtaPair =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recEtaPair"));
  TH1F *H_h_EtaPair =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recEtaPair"));
  TH1F *Z_h_SepAngle =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recSepAngle"));
  TH1F *H_h_SepAngle =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recSepAngle"));
  TH1F *Z_h_PtSingle =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recPtSingle"));
  TH1F *H_h_PtSingle =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recPtSingle"));
  TH1F *Z_h_ChAsymm =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recChAsymm"));
  TH1F *H_h_ChAsymm =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recChAsymm"));
  TH1F *Z_h_D0overSigma =
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recD0overSigma"));
  TH1F *H_h_D0overSigma =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recD0overSigma"));

  auto c_recMInv = new TCanvas("c_recMInv", "c_recMInv", 1000, 800);
  // H_h_MInv->Scale(1 / H_h_MInv->Integral());
  // Z_h_MInv->Scale(1 / Z_h_MInv->Integral());
  H_h_MInv->SetLineColor(kRed);
  H_h_MInv->SetFillColor(kRed);
  H_h_MInv->SetFillStyle(3002);
  Z_h_MInv->SetLineColor(kBlue);
  Z_h_MInv->SetFillColor(kBlue);
  Z_h_MInv->SetFillStyle(3003);
  display_overflow(H_h_MInv);
  display_overflow(Z_h_MInv);
  H_h_MInv->GetYaxis()->SetRangeUser(0., H_h_MInv->GetMaximum() * 1.15);
  H_h_MInv->Draw("HISTO");
  Z_h_MInv->Draw("HISTO, SAME");
  H_h_MInv->GetXaxis()->SetLabelFont(42);
  H_h_MInv->GetXaxis()->SetLabelSize(0.035);
  H_h_MInv->GetXaxis()->SetTitleFont(42);
  H_h_MInv->GetXaxis()->SetTitleOffset(1.15);
  H_h_MInv->GetXaxis()->SetTitleSize(0.04);
  H_h_MInv->GetYaxis()->SetLabelFont(42);
  H_h_MInv->GetYaxis()->SetLabelSize(0.035);
  H_h_MInv->GetYaxis()->SetTitleFont(42);
  H_h_MInv->GetYaxis()->SetTitleOffset(1.38);
  H_h_MInv->GetYaxis()->SetDecimals();
  c_recMInv->SetLeftMargin(0.12);
  c_recMInv->SetBottomMargin(0.15);
  H_h_MInv->GetYaxis()->SetTitle("Density");
  auto leg0 = new TLegend(0.7, 0.5, 0.88, 0.61);
  gStyle->SetOptStat(0);
  leg0->AddEntry(H_h_MInv, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg0->AddEntry(Z_h_MInv, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg0->SetTextSize(0.037);
  leg0->SetTextFont(42);
  leg0->SetLineColor(kWhite);
  leg0->Draw();
  // gPad->Modified();
  /*TPaveText *pt1 = new TPaveText(160, 0.10, 195, 0.15);
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
  pt1->Draw();*/
  gStyle->SetOptTitle(0);
  // gROOT->SetStyle("Plain");
  // gROOT->ForceStyle();
  gPad->Update();
  c_recMInv->SaveAs("./plotsMInv/" + prefix + "_c_recMInv_" + decay0 + decay1 +
                    ".png");
  c_recMInv->Close();

  auto c_recDeltaPhi = new TCanvas("c_recDeltaPhi", "c_recDeltaPhi", 800, 800);
  // H_h_DeltaPhi->Scale(1 / H_h_DeltaPhi->Integral());
  // Z_h_DeltaPhi->Scale(1 / Z_h_DeltaPhi->Integral());
  H_h_DeltaPhi->SetLineColor(kRed);
  H_h_DeltaPhi->SetFillColor(kRed);
  H_h_DeltaPhi->SetFillStyle(3002);
  Z_h_DeltaPhi->SetLineColor(kBlue);
  Z_h_DeltaPhi->SetFillColor(kBlue);
  Z_h_DeltaPhi->SetFillStyle(3003);
  H_h_DeltaPhi->SetMinimum(0.);
  display_underoverflow(H_h_DeltaPhi);
  display_underoverflow(Z_h_DeltaPhi);
  H_h_DeltaPhi->GetYaxis()->SetRangeUser(0., Z_h_DeltaPhi->GetMaximum() * 1.2);
  H_h_DeltaPhi->Draw("HISTO");
  Z_h_DeltaPhi->Draw("HISTO, SAME");
  auto leg8 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg8->AddEntry(H_h_DeltaPhi, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg8->AddEntry(Z_h_DeltaPhi, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg8->Draw();
  gPad->Modified();
  c_recDeltaPhi->Update();
  c_recDeltaPhi->SaveAs("./plotsDeltaPhi/" + prefix + "_c_recDeltaPhi_" +
                        decay0 + decay1 + ".png");
  c_recDeltaPhi->Close();

  auto c_recDeltaR = new TCanvas("c_recDeltaR", "c_recDeltaR", 800, 800);
  // H_h_DeltaR->Scale(1 / H_h_DeltaR->Integral());
  // Z_h_DeltaR->Scale(1 / Z_h_DeltaR->Integral());
  H_h_DeltaR->SetLineColor(kRed);
  H_h_DeltaR->SetFillColor(kRed);
  H_h_DeltaR->SetFillStyle(3002);
  Z_h_DeltaR->SetLineColor(kBlue);
  Z_h_DeltaR->SetFillColor(kBlue);
  Z_h_DeltaR->SetFillStyle(3003);
  display_overflow(H_h_DeltaR);
  display_overflow(Z_h_DeltaR);
  H_h_DeltaR->GetYaxis()->SetRangeUser(0., Z_h_DeltaR->GetMaximum() * 1.2);
  H_h_DeltaR->Draw("HISTO");
  Z_h_DeltaR->Draw("HISTO, SAME");
  auto leg10 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg10->AddEntry(H_h_DeltaR, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg10->AddEntry(Z_h_DeltaR, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg10->Draw();
  gPad->Modified();
  c_recDeltaR->Update();
  c_recDeltaR->SaveAs("./plotsDeltaR/" + prefix + "_c_recDeltaR_" + decay0 +
                      decay1 + ".png");
  c_recDeltaR->Close();

  auto c_recDeltaEta = new TCanvas("c_recDeltaEta", "c_recDeltaEta", 800, 800);
  // H_h_DeltaEta->Scale(1 / H_h_DeltaEta->Integral());
  // Z_h_DeltaEta->Scale(1 / Z_h_DeltaEta->Integral());
  H_h_DeltaEta->SetLineColor(kRed);
  H_h_DeltaEta->SetFillColor(kRed);
  H_h_DeltaEta->SetFillStyle(3002);
  Z_h_DeltaEta->SetLineColor(kBlue);
  Z_h_DeltaEta->SetFillColor(kBlue);
  Z_h_DeltaEta->SetFillStyle(3003);
  display_overflow(H_h_DeltaEta);
  display_overflow(Z_h_DeltaEta);
  H_h_DeltaEta->GetYaxis()->SetRangeUser(0., Z_h_DeltaEta->GetMaximum() * 1.2);
  H_h_DeltaEta->Draw("HISTO");
  Z_h_DeltaEta->Draw("HISTO, SAME");
  auto leg10b = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg10b->AddEntry(H_h_DeltaEta, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg10b->AddEntry(Z_h_DeltaEta, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg10b->Draw();
  gPad->Modified();
  c_recDeltaEta->Update();
  c_recDeltaEta->SaveAs("./plotsDeltaEta/c_recDeltaEta_" + decay0 + decay1 + ".png");
  c_recDeltaEta->Close();

  auto c_recEtaPair = new TCanvas("c_recEtaPair", "c_recEtaPair", 800, 800);
  // H_h_EtaPair->Scale(1 / H_h_EtaPair->Integral());
  // Z_h_EtaPair->Scale(1 / Z_h_EtaPair->Integral());
  H_h_EtaPair->SetLineColor(kRed);
  H_h_EtaPair->SetFillColor(kRed);
  H_h_EtaPair->SetFillStyle(3002);
  Z_h_EtaPair->SetLineColor(kBlue);
  Z_h_EtaPair->SetFillColor(kBlue);
  Z_h_EtaPair->SetFillStyle(3003);
  H_h_EtaPair->SetMinimum(0.);
  display_underoverflow(H_h_EtaPair);
  display_underoverflow(Z_h_EtaPair);
  H_h_EtaPair->GetYaxis()->SetRangeUser(0., H_h_EtaPair->GetMaximum() * 1.2);
  H_h_EtaPair->Draw("HISTO");
  Z_h_EtaPair->Draw("HISTO, SAME");
  auto leg4 = new TLegend(0.68, 0.78, 0.88, 0.88);
  // gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg4->AddEntry(H_h_EtaPair, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg4->AddEntry(Z_h_EtaPair, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg4->Draw();
  gPad->Modified();
  c_recEtaPair->SaveAs("./plotsEtaPair/" + prefix + "_c_recEtaPair_" + decay0 +
                       decay1 + ".png");
  c_recEtaPair->Close();

  auto c_recPtPair = new TCanvas("c_recPtPair", "c_recPtPair", 800, 800);
  // H_h_PtPair->Scale(1 / H_h_PtPair->Integral());
  // Z_h_PtPair->Scale(1 / Z_h_PtPair->Integral());
  H_h_PtPair->SetLineColor(kRed);
  H_h_PtPair->SetFillColor(kRed);
  H_h_PtPair->SetFillStyle(3002);
  Z_h_PtPair->SetLineColor(kBlue);
  Z_h_PtPair->SetFillColor(kBlue);
  Z_h_PtPair->SetFillStyle(3003);
  H_h_PtPair->GetXaxis()->SetRangeUser(0., 250.);
  display_overflow(H_h_PtPair);
  display_overflow(Z_h_PtPair);
  H_h_PtPair->GetYaxis()->SetRangeUser(0., Z_h_PtPair->GetMaximum() * 1.2);
  H_h_PtPair->Draw("HISTO");
  Z_h_PtPair->Draw("HISTO, SAME");
  auto leg2 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg2->AddEntry(H_h_PtPair, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg2->AddEntry(Z_h_PtPair, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg2->Draw();
  gPad->Modified();
  c_recPtPair->SaveAs("./plotsPtPair/" + prefix + "_c_recPtPair_" + decay0 +
                      decay1 + ".png");
  c_recPtPair->Close();

  auto c_recSepAngle = new TCanvas("c_recSepAngle", "c_recSepAngle", 800, 800);
  // H_h_SepAngle->Scale(1 / H_h_SepAngle->Integral());
  // Z_h_SepAngle->Scale(1 / Z_h_SepAngle->Integral());
  H_h_SepAngle->SetLineColor(kRed);
  H_h_SepAngle->SetFillColor(kRed);
  H_h_SepAngle->SetFillStyle(3002);
  Z_h_SepAngle->SetLineColor(kBlue);
  Z_h_SepAngle->SetFillColor(kBlue);
  Z_h_SepAngle->SetFillStyle(3003);
  H_h_SepAngle->GetYaxis()->SetRangeUser(0., Z_h_SepAngle->GetMaximum() * 1.2);
  H_h_SepAngle->Draw("HISTO");
  Z_h_SepAngle->Draw("HISTO, SAME");
  auto leg12 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg12->AddEntry(H_h_SepAngle, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg12->AddEntry(Z_h_SepAngle, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg12->Draw();
  gPad->Modified();
  c_recSepAngle->SaveAs("./plotsSepAngle/" + prefix + "_c_recSepAngle_" +
                        decay0 + decay1 + ".png");
  c_recSepAngle->Close();

  auto c_recPtSingle = new TCanvas("c_recPtSingle", "c_recPtSingle", 800, 800);
  // H_h_PtSingle->Scale(1 / H_h_PtSingle->Integral());
  // Z_h_PtSingle->Scale(1 / Z_h_PtSingle->Integral());
  H_h_PtSingle->SetLineColor(kRed);
  H_h_PtSingle->SetFillColor(kRed);
  H_h_PtSingle->SetFillStyle(3002);
  Z_h_PtSingle->SetLineColor(kBlue);
  Z_h_PtSingle->SetFillColor(kBlue);
  Z_h_PtSingle->SetFillStyle(3003);
  H_h_PtSingle->GetXaxis()->SetRangeUser(0., 250.);
  H_h_PtSingle->GetYaxis()->SetRangeUser(0., Z_h_PtSingle->GetMaximum() * 1.2);
  H_h_PtSingle->Draw("HISTO");
  Z_h_PtSingle->Draw("HISTO, SAME");
  display_overflow(H_h_PtSingle);
  display_overflow(Z_h_PtSingle);
  auto leg13 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg13->AddEntry(H_h_PtSingle, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg13->AddEntry(Z_h_PtSingle, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg13->Draw();
  gPad->Modified();
  c_recPtSingle->SaveAs("./plotsPtSingle/" + prefix + "_c_recPtSingle_" +
                        decay0 + decay1 + ".png");
  c_recPtSingle->Close();

  auto c_recChAsymm = new TCanvas("c_recChAsymm", "c_recChAsymm", 800, 800);
  // H_h_ChAsymm->Scale(1 / H_h_ChAsymm->Integral());
  // Z_h_ChAsymm->Scale(1 / Z_h_ChAsymm->Integral());
  H_h_ChAsymm->SetLineColor(kRed);
  H_h_ChAsymm->SetFillColor(kRed);
  H_h_ChAsymm->SetFillStyle(3002);
  Z_h_ChAsymm->SetLineColor(kBlue);
  Z_h_ChAsymm->SetFillColor(kBlue);
  Z_h_ChAsymm->SetFillStyle(3003);
  display_overflow(H_h_PtPair);
  display_overflow(Z_h_PtPair);
  H_h_ChAsymm->GetYaxis()->SetRangeUser(0., Z_h_ChAsymm->GetMaximum() * 1.2);
  H_h_ChAsymm->Draw("HISTO");
  Z_h_ChAsymm->Draw("HISTO, SAME");
  auto leg14 = new TLegend(0.2, 0.78, 0.38, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg14->AddEntry(H_h_ChAsymm, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg14->AddEntry(Z_h_ChAsymm, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg14->Draw();
  /*H_h_ChAsymm->SetMinimum(0.001);
  Z_h_ChAsymm->SetMinimum(0.001);
  gPad->SetLogy();*/
  gPad->Modified();
  c_recChAsymm->SaveAs("./plotsChAsymm/" + prefix + "_c_recChAsymm_" + decay0 +
                       decay1 + ".png");
  c_recChAsymm->Close();

  auto c_recD0overSigma =
      new TCanvas("c_recD0overSigma", "c_recD0overSigma", 800, 800);
  // H_h_D0overSigma->Scale(1 / H_h_D0overSigma->Integral());
  // Z_h_D0overSigma->Scale(1 / Z_h_D0overSigma->Integral());
  H_h_D0overSigma->SetLineColor(kRed);
  H_h_D0overSigma->SetFillColor(kRed);
  H_h_D0overSigma->SetFillStyle(3002);
  Z_h_D0overSigma->SetLineColor(kBlue);
  Z_h_D0overSigma->SetFillColor(kBlue);
  Z_h_D0overSigma->SetFillStyle(3003);
  display_underoverflow(H_h_D0overSigma);
  display_underoverflow(Z_h_D0overSigma);
  H_h_D0overSigma->GetYaxis()->SetRangeUser(0., Z_h_D0overSigma->GetMaximum() *
                                                    1.2);
  H_h_D0overSigma->Draw("HISTO");
  Z_h_D0overSigma->Draw("HISTO, SAME");
  gPad->SetLogx();
  auto leg15 = new TLegend(0.68, 0.78, 0.88, 0.88);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  leg15->AddEntry(H_h_D0overSigma, "#mu#mu #rightarrow H#nu_{#tau}#bar{#nu}_{#tau}", "f");
  leg15->AddEntry(Z_h_D0overSigma, "#mu#mu#rightarrow#tau#tau#mu#mu", "f");
  leg15->Draw();
  gPad->Modified();
  c_recD0overSigma->SaveAs("./plotsD0overSigma/" + prefix + "_c_recD0overSigma_" + decay0 + decay1 + ".png");
  c_recD0overSigma->Close();
}

void tau_plotter(const TString &filename, const TString &decay0 = "all",
                 const TString &decay1 = "all", const TString &prefix = "") {
  // Open and read energy correction file
  TFile *histo_file = new TFile(filename, "READ");

  if (!histo_file->IsOpen())
    throw std::invalid_argument("histograms filename not valid");

  TList *Zhisto_list =
      (TList *)histo_file->Get("Z_" + decay0 + "+" + decay1 + "_histograms");
  TList *Hhisto_list =
      (TList *)histo_file->Get("H_" + decay0 + "+" + decay1 + "_histograms");

  TH1F *h_D0 = static_cast<TH1F *>(Hhisto_list->FindObject("h_recD0Single"));
  h_D0->Add(static_cast<TH1F *>(Zhisto_list->FindObject("h_recD0Single")));
  TH1F *h_SigmaD0 =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recSigmaD0Single"));
  h_SigmaD0->Add(
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recSigmaD0Single")));
  TH1F *h_D0overSigma =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recD0overSigma"));
  h_D0overSigma->Add(
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recD0overSigma")));
  TH1F *h_ChAsymm =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recChAsymm"));
  h_ChAsymm->Add(static_cast<TH1F *>(Zhisto_list->FindObject("h_recChAsymm")));
  TH1F *h_PtSingle =
      static_cast<TH1F *>(Hhisto_list->FindObject("h_recPtSingle"));
  h_PtSingle->Add(
      static_cast<TH1F *>(Zhisto_list->FindObject("h_recPtSingle")));

  auto c_D0 = new TCanvas("c_D0", "c_D0", 800, 800);
  gStyle->SetOptStat(111111);
  h_D0->GetXaxis()->SetTitleOffset(1.2);
  h_D0->Draw("HISTO");
  display_underoverflow(h_D0);
  gPad->SetLogx();
  c_D0->SaveAs("./plotsSingleD0/" + prefix + "_c_D0_" + decay0 + decay1 +
               ".png");
  c_D0->Close();

  auto c_SigmaD0 = new TCanvas("c_SigmaD0", "c_SigmaD0", 800, 800);
  h_SigmaD0->GetXaxis()->SetTitleOffset(1.2);
  h_SigmaD0->Draw("HISTO");
  gPad->SetLogx();
  display_underoverflow(h_SigmaD0);
  c_SigmaD0->SaveAs("./plotsSingleSigmaD0/" + prefix + "_c_SigmaD0_" + decay0 +
                    decay1 + ".png");
  c_SigmaD0->Close();

  auto c_D0overSigma = new TCanvas("c_D0overSigma", "c_D0overSigma", 800, 800);
  h_D0overSigma->GetXaxis()->SetTitleOffset(1.2);
  h_D0overSigma->Draw("HISTO");
  gPad->SetLogx();
  display_underoverflow(h_D0overSigma);
  c_D0overSigma->SaveAs("./plotsD0overSigma/" + prefix + "_c_D0overSigma_" +
                        decay0 + decay1 + ".png");
  c_D0overSigma->Close();

  auto c_ChAsymm = new TCanvas("c_ChAsymm", "c_ChAsymm", 800, 800);
  h_ChAsymm->GetXaxis()->SetTitleOffset(1.2);
  h_ChAsymm->Draw("HISTO");
  display_underoverflow(h_ChAsymm);
  c_ChAsymm->SaveAs("./plotsChAsymm/" + prefix + "_c_ChAsymm_" + decay0 +
                    decay1 + ".png");
  c_ChAsymm->Close();

  auto c_PtSingle = new TCanvas("c_PtSingle", "c_PtSingle", 800, 800);
  h_PtSingle->GetXaxis()->SetTitleOffset(1.2);
  h_PtSingle->Draw("HISTO");
  display_underoverflow(h_PtSingle);
  c_PtSingle->SaveAs("./plotsPtSingle/" + prefix + "_c_PtSingle_" + decay0 +
                     decay1 + ".png");
  c_PtSingle->Close();
}

void TauPairsPlotFake(
    const TString filename =
        "/eos/user/l/lvalla/MuColl/HvsZ/HvsMuMu_Histograms.root") {

  const TString prefix = "MuMu";

  HvsZ_plotter(filename, "all", "all", prefix);
  HvsZ_plotter(filename, "had", "had", prefix);
  HvsZ_plotter(filename, "1pi", "1pi", prefix);
  HvsZ_plotter(filename, "e", "e", prefix);
  HvsZ_plotter(filename, "mu", "mu", prefix);
  HvsZ_plotter(filename, "e", "1pi", prefix);
  HvsZ_plotter(filename, "had", "e", prefix);
  HvsZ_plotter(filename, "hade", "hade", prefix);

  /*tau_plotter(filename, "all", "all", prefix);
  tau_plotter(filename, "1pi", "1pi", prefix);
  tau_plotter(filename, "had", "had", prefix);
  tau_plotter(filename, "e", "e", prefix);
  tau_plotter(filename, "mu", "mu", prefix);
  tau_plotter(filename, "e", "1pi", prefix);
  tau_plotter(filename, "had", "e", prefix);*/
}