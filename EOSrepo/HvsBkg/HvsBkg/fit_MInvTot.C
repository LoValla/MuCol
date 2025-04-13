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

void fit_plotter(const TString &decay0, const TString &decay1) {

  const TString filenameZ =
      "/eos/user/l/lvalla/MuColl/HvsZ/HvsZ_Histograms.root";
  const TString filenameMuMu =
      "/eos/user/l/lvalla/MuColl/HvsZ/HvsMuMu_Histograms.root";

  gStyle->SetTextFont(42);

  TRandom3 *ran = new TRandom3(0);

  TFile *histoZ_file = new TFile(filenameZ, "READ");

  TFile *histoMuMu_file = new TFile(filenameMuMu, "READ");

  if (!histoZ_file->IsOpen() || !histoMuMu_file->IsOpen())
    throw std::invalid_argument("histograms filename(s) not valid");

  TList *Hhisto_list =
      (TList *)histoZ_file->Get("H_" + decay0 + "+" + decay1 + "_histograms");
  TList *Zhisto_list =
      (TList *)histoZ_file->Get("Z_" + decay0 + "+" + decay1 + "_histograms");
  TList *MuMuhisto_list = (TList *)histoMuMu_file->Get("Z_" + decay0 + "+" +
                                                       decay1 + "_histograms");
  TList *MC_Hhisto_list = (TList *)histoZ_file->Get("H_MC_histograms");
  TList *MC_Zhisto_list = (TList *)histoZ_file->Get("Z_MC_histograms");
  TList *MC_MuMuhisto_list = (TList *)histoMuMu_file->Get("Z_MC_histograms");

  TH1F *h_H_MInv = static_cast<TH1F *>(Hhisto_list->FindObject("h_recMInv"));
  TH1F *h_Z_MInv = static_cast<TH1F *>(Zhisto_list->FindObject("h_recMInv"));
  TH1F *h_MuMu_MInv =
      static_cast<TH1F *>(MuMuhisto_list->FindObject("h_recMInv"));

  TH1F *h_H_MInv_test =
      new TH1F("h_H_MInv_test", "h_H_MInv_test", 1, 0., 210.);
  TH1F *h_Z_MInv_test =
      new TH1F("h_Z_MInv_test", "h_Z_MInv_test", 1, 0., 210.);
  TH1F *h_MuMu_MInv_test =
      new TH1F("h_MuMu_MInv_test", "h_MuMu_MInv_test", 1, 0., 210.);

  Int_t H_MCentries =
      static_cast<TH1F *>(MC_Hhisto_list->FindObject("h_mcMInv"))->GetEntries();
  Int_t Z_MCentries =
      static_cast<TH1F *>(MC_Zhisto_list->FindObject("h_mcMInv"))->GetEntries();
  Int_t MuMu_MCentries =
      static_cast<TH1F *>(MC_MuMuhisto_list->FindObject("h_mcMInv"))
          ->GetEntries();

  Int_t H_Recoentries = h_H_MInv->GetEntries();
  Int_t Z_Recoentries = h_Z_MInv->GetEntries();
  Int_t MuMu_Recoentries = h_MuMu_MInv->GetEntries();

  Double_t intlumi = 1.E3; // in fb-1
  Double_t Hxsec = 30.9;   // fb
  Double_t Zxsec = 73.8;   // fb
  Double_t MuMuxsec = 193;
  Double_t Heff = Double_t(H_Recoentries) / Double_t(H_MCentries);
  Double_t Zeff = Double_t(Z_Recoentries) / Double_t(Z_MCentries);
  Double_t MuMueff = Double_t(MuMu_Recoentries) / Double_t(MuMu_MCentries);

  Double_t nsig = intlumi * Hxsec * Heff;
  Double_t nbkgZ = intlumi * Zxsec * Zeff;
  Double_t nbkgMuMu = intlumi * MuMuxsec * MuMueff;
  Double_t nbkg = nbkgZ + nbkgMuMu;

  /*std::cout << "nsig" << nsig << '\n';
  std::cout << "nbkg" << nbkg << '\n';*/

  Int_t nevents = nsig + nbkg;

  /*TH1F *h_bkg_MInv = new TH1F("h_bkg_MInv", "h_bkg_MInv", 30., 0., 210.);
  TH1F *h_tot_MInv = new TH1F("h_tot_MInv", "h_tot_MInv", 30., 0., 210.);
  h_H_MInv->Sumw2();
  h_Z_MInv->Sumw2();
  h_MuMu_MInv->Sumw2();
  h_H_MInv->Scale(nsig);
  h_Z_MInv->Scale(nbkgZ);
  h_MuMu_MInv->Scale(nbkgMuMu);
  h_bkg_MInv->Add(h_Z_MInv, h_MuMu_MInv);
  h_bkg_MInv->Sumw2();
  h_tot_MInv->Add(h_H_MInv, h_bkg_MInv);*/

  TH1F *h_bkg_MInv = new TH1F("h_bkg_MInv", "h_bkg_MInv", 1, 0., 210.);
  TH1F *h_tot_MInv = new TH1F("h_tot_MInv", "h_tot_MInv", 1, 0., 210.);
  h_H_MInv->Sumw2();
  h_Z_MInv->Sumw2();
  h_MuMu_MInv->Sumw2();
  h_H_MInv->Scale(nsig);
  h_Z_MInv->Scale(nbkgZ);
  h_MuMu_MInv->Scale(nbkgMuMu);
  h_bkg_MInv->Add(h_Z_MInv_test, h_MuMu_MInv_test);
  h_bkg_MInv->Sumw2();
  h_tot_MInv->Add(h_H_MInv_test, h_bkg_MInv_test);

  // for (int i=0; i<nfits; i++){

  RooRealVar invmass("invmass", "invmass", 30., 0., 210., "GeV/c^{2}");

  RooDataHist *histo_data =
      new RooDataHist("data", "data", RooArgList(invmass), h_tot_MInv);

  RooRealVar Nsig("Nsig", "Nsig", nsig, 0, nsig + 1000);
  RooRealVar Nbkg("Nbkg", "Nbkg", nbkg, 0, nbkg + 1000);

  RooDataHist *Bkg_histo_data = new RooDataHist(
      "Bkg_histo_data", "Bkg_histo_data", RooArgList(invmass), h_bkg_MInv);

  RooDataHist *H_histo_data = new RooDataHist("H_histo_data", "H_histo_data",
                                              RooArgList(invmass), h_H_MInv);
  RooHistPdf *Bkghistpdf =
      new RooHistPdf("Bkghistpdf", "Bkghistpdf", invmass, *Bkg_histo_data);
  RooHistPdf *Hhistpdf =
      new RooHistPdf("Hhistpdf", "Hhistpdf", invmass, *H_histo_data);
  RooAddPdf model("model", "model", RooArgList(*Hhistpdf, *Bkghistpdf),
                  RooArgList(Nsig, Nbkg));

  RooMCStudy *mcstudy_ex = new RooMCStudy(
      model, RooArgSet(invmass), RooFit::Binned(false), RooFit::Silence(),
      RooFit::Extended(true),
      RooFit::FitOptions(RooFit::Save(true), RooFit::PrintEvalErrors(0)));

  mcstudy_ex->generateAndFit(10, 0, kTRUE, 0);

  RooAbsData *single_toy_data = mcstudy_ex->genData(1);
  const RooFitResult *single_toy_res = mcstudy_ex->fitResult(1);
  RooArgSet single_toy_pars = single_toy_res->floatParsFinal();

  static_cast<RooRealVar &>(single_toy_pars["Nbkg"]).Print();
  static_cast<RooRealVar &>(single_toy_pars["Nsig"]).Print();

  RooAddPdf model_toy(
      "model_toy", "model_toy", RooArgList(*Bkghistpdf, *Hhistpdf),
      RooArgList(static_cast<RooRealVar &>(single_toy_pars["Nbkg"]),
                 static_cast<RooRealVar &>(single_toy_pars["Nsig"])));
  model_toy.fitTo(*single_toy_data, RooFit::Save());
  TCanvas *can4 = new TCanvas("can4", "can4", 900, 600);
  RooPlot *frame4 = invmass.frame(Bins(30));
  single_toy_data->plotOn(frame4, Name("data"));
  model_toy.plotOn(frame4, Name("Fit"));
  // model.plotOn(frame4, Components(RooArgSet(p_sig, p_bkg)), DrawOption("L"));
  model_toy.plotOn(frame4, Components(RooArgSet(*Hhistpdf)), DrawOption("L"),
                   LineColor(kRed), Name("Higgs"));
  model_toy.plotOn(frame4, Components(RooArgSet(*Bkghistpdf)), DrawOption("L"),
                   LineColor(kGreen), Name("Background"));
  frame4->Draw();

  // TPaveText *pt2 = new TPaveText(8, 1200, 100, 1600); //FOR 1PI1PI
  TPaveText *pt2 = new TPaveText(-0.3, 390, -0.08, 630); // FOR ALLALL
  pt2->SetFillColor(0);
  pt2->SetMargin(0);
  pt2->SetFillStyle(0);
  pt2->SetBorderSize(0);
  pt2->SetTextAlign(11);
  // pt2->SetTextFont(50);
  pt2->AddText("#bf{Muon Collider}");
  pt2->AddText("#it{Simulation}");
  pt2->AddText("#sqrt{s} = 3 TeV, L = 1 ab^{-1}");
  pt2->AddText("  ");
  pt2->AddText("p_{T}(#tau) > 10 GeV/c");
  pt2->AddText("|#eta(#tau)| < 2.5");
  pt2->Draw();

  auto legend2 = new TLegend(0.665, 0.65, 0.915, 0.85);
  legend2->AddEntry("data", "Pseudo-data", "lep");
  legend2->AddEntry("fit", "Fit", "l");
  legend2->AddEntry("higgs", "H#rightarrow #tau_{h}#tau_{h}", "l");
  legend2->AddEntry("background", "Background", "l");
  legend2->SetBorderSize(0);
  legend2->SetFillColorAlpha(kBlue, 0.);
  legend2->SetTextSize(0.045);
  legend2->Draw();

  frame4->GetYaxis()->SetTitle("Events");
  frame4->GetYaxis()->SetTitleOffset(1.2);
  frame4->GetXaxis()->SetTitle("m^{vis}_{#tau#tau} [GeV/c^{2}]");
  frame4->GetXaxis()->SetTitleOffset(1.1);
  frame4->GetXaxis()->SetLabelOffset(0.013);
  frame4->GetXaxis()->SetLabelSize(0.04);
  frame4->GetYaxis()->SetLabelSize(0.04);
  frame4->GetXaxis()->SetTitleSize(0.045);
  frame4->GetYaxis()->SetTitleSize(0.045);
  frame4->SetTitle("");
  can4->SaveAs("roofit_plots/SingleMCtoy_HvsTot_" + decay0 + decay1 + ".png");
  can4->Close();

  return;

  RooMCStudy *mcstudy = new RooMCStudy(
      model, RooArgSet(invmass), RooFit::Binned(true), RooFit::Silence(),
      RooFit::Extended(true),
      RooFit::FitOptions(RooFit::Save(true), RooFit::PrintEvalErrors(0)));

  mcstudy->generateAndFit(100000);

  TH1 *h_sigevents =
      mcstudy->fitParDataSet().createHistogram("h_sigevents", Nsig);

  std::cout << "MEAN NUMBER OF SIGNAL EVENTS: " << h_sigevents->GetMean()
            << '\n';
  std::cout << "STD DEV OF SIGNAL EVENTS: " << h_sigevents->GetStdDev() << '\n';

  TH1 *h_bkgevents =
      mcstudy->fitParDataSet().createHistogram("h_bkgevents", Nbkg);

  std::cout << "MEAN NUMBER OF BKG EVENTS: " << h_bkgevents->GetMean() << '\n';
  std::cout << "STD DEV OF BKG EVENTS: " << h_bkgevents->GetStdDev() << '\n';

  RooPlot *frame1 = mcstudy->plotParam(Nsig, RooFit::Bins(64));
  RooPlot *frame2 = mcstudy->plotError(Nsig, RooFit::Bins(64));
  RooPlot *frame3 =
      mcstudy->plotPull(Nsig, RooFit::Bins(40), RooFit::FitGauss(false));
  RooRealVar pullMean("pullMean", "Mean of pull", 0, -10, 10);
  RooRealVar pullSigma("pullSigma", "Width of pull", 1, 0.1, 5);
  pullMean.setPlotLabel(
      "pull #mu"); // optional (to get nicer plot labels if you want)
  pullSigma.setPlotLabel("pull #sigma"); // optional
  RooGaussian pullGauss("pullGauss", "Gaussian of pull", *frame3->getPlotVar(),
                        pullMean, pullSigma);
  pullGauss.fitTo(const_cast<RooDataSet &>(mcstudy->fitParDataSet()),
                  RooFit::Minos(0), RooFit::PrintLevel(-1));
  pullGauss.plotOn(frame3);
  pullGauss.paramOn(frame3, RooFit::Layout(0.65, 0.9, 0.9));

  TCanvas *can1 = new TCanvas("can1", "can1", 1200, 900);
  frame1->Draw();
  can1->SaveAs("roofit_plots/Nsig_HvsTot_" + decay0 + decay1 + ".png");
  TCanvas *can2 = new TCanvas("can2", "can2", 1200, 900);
  frame2->Draw();
  can2->SaveAs("roofit_plots/NsigError_HvsTot_" + decay0 + decay1 + ".png");
  TCanvas *can3 = new TCanvas("can3", "can3", 1200, 900);
  frame3->Draw();
  can3->SaveAs("roofit_plots/NsigPull_HvsTot_" + decay0 + decay1 + ".png");

  TCanvas *can = new TCanvas("can", "can", 1200, 900);
  RooPlot *xframe = invmass.frame();

  histo_data->plotOn(xframe, Name("data"));
  model.plotOn(xframe, Name("fit"));
  // model.plotOn(xframe, Components(RooArgSet(p_sig,p_bkg)),
  DrawOption("L"),
      // LineColor(kBlue)), Name("fit") ;
      model.plotOn(xframe, Components(RooArgSet(*Hhistpdf)), DrawOption("L"),
                   LineColor(kRed), Name("higgs"));
  model.plotOn(xframe, Components(RooArgSet(*Bkghistpdf)), DrawOption("L"),
               LineColor(kGreen), Name("Background"));
  histo_data->plotOn(xframe);

  xframe->GetYaxis()->SetTitle("Events / (8 GeV)");
  xframe->GetYaxis()->SetTitleOffset(1.3);
  xframe->GetXaxis()->SetTitle("#tau#tau invariant mass [GeV]");

  xframe->Draw();

  // TPaveText *pt1 = new TPaveText(8, 1200, 100, 1600); //FOR 1PI1PI
  TPaveText *pt1 = new TPaveText(110, 300, 150, 570); // FOR ALLALL
  pt1->SetFillColor(0);
  pt1->SetMargin(0);
  pt1->SetFillStyle(0);
  pt1->SetBorderSize(0);
  pt1->SetTextAlign(11);
  // pt1->SetTextFont(50);
  pt1->AddText("#bf{Muon Collider}");
  pt1->AddText("#it{Simulation}");
  pt1->AddText("#sqrt{s} = 3 TeV, L = 1 ab^{-1}");
  pt1->AddText("  ");
  pt1->AddText("p_{T}(#tau) > 10 GeV/c");
  pt1->AddText("|#eta| < 2.5");
  pt1->Draw();

  auto legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  legend->AddEntry("data", "data", "lep");
  legend->AddEntry("fit", "fit total", "l");
  legend->AddEntry("higgs", "H#rightarrow #tau#tau", "l");
  legend->AddEntry("Background", "Background", "l");
  legend->Draw();

  can->SaveAs("roofit_plots/HvsTot_" + decay0 + decay1 + "_fit.png");
  can->Close();

  std::cout << "EXPECTED NUMBER OF EVENTS" << '\n';
  std::cout << "nbkg = " << nbkg << '\n';
  std::cout << "nsig = " << nsig << '\n';

  // res->Print();

  return;
}

void fit_MInvTot() {
  // fit_plotter(filename, "all", "all");
  fit_plotter("had", "had");
  // fit_plotter(filename, "had", "e");
  // fit_plotter(filename, "1pi", "1pi");
}
