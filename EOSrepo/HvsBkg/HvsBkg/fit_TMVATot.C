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

  TH1F *h_H_BDT = static_cast<TH1F *>(dir->Get("MVA_BDT_S"));
  TH1F *h_Z_BDT = static_cast<TH1F *>(dir->Get("MVA_BDT_B"));

  TFile *file = TFile::Open("HvsAll_Histograms.root");

  TList *HMCList = (TList *)file->Get("H_MC_histograms");
  TList *ZMCList = (TList *)file->Get("Z_MC_histograms");
  TList *MuMuMCList = (TList *)file->Get("MuMu_MC_histograms");
  TList *HList = (TList *)file->Get("H_had+had_histograms");
  TList *ZList = (TList *)file->Get("Z_had+had_histograms");
  TList *MuMuList = (TList *)file->Get("MuMu_had+had_histograms");

  Int_t H_MCentries =
      (static_cast<TH1F *>(HMCList->FindObject("h_mcPtPair")))->GetEntries();
  Int_t Z_MCentries =
      (static_cast<TH1F *>(ZMCList->FindObject("h_mcPtPair")))->GetEntries();
  Int_t ttmm_MCentries =
      (static_cast<TH1F *>(MuMuMCList->FindObject("h_mcPtPair")))->GetEntries();

  Int_t H_Recoentries =
      (static_cast<TH1F *>(HList->FindObject("h_recPtPair")))->GetEntries();
  Int_t Z_Recoentries =
      (static_cast<TH1F *>(ZList->FindObject("h_recPtPair")))->GetEntries();
  Int_t ttmm_Recoentries =
      (static_cast<TH1F *>(MuMuList->FindObject("h_recPtPair")))->GetEntries();

  file->Close();
  file->Delete();

  Double_t intlumi = 1.E3; // in fb-1
  Double_t Hxsec = 30.9;   // fb
  Double_t Zxsec = 73.8;   // fb
  Double_t ttmmxsec = 193; // fb
  Double_t Heff = Double_t(H_Recoentries) / Double_t(H_MCentries);
  Double_t Zeff = Double_t(Z_Recoentries) / Double_t(Z_MCentries);
  Double_t ttmmeff = Double_t(ttmm_Recoentries) / Double_t(ttmm_MCentries);

  std::cout << "Heff = " << Heff << '\n';
  std::cout << "Zeff = " << Zeff << '\n';
  std::cout << "MuMueff = " << ttmmeff << '\n';
  std::cout << "Hreco = " << Heff * 300000. << '\n';
  std::cout << "Zreco = " << Zeff * 300000. << '\n';
  std::cout << "MuMureco = " << ttmmeff * 300000. << '\n';

  Double_t nsig = intlumi * Hxsec * Heff; // multiply by 2 (train + test)
  Double_t nbkg = intlumi * (Zxsec * Zeff + ttmmxsec * ttmmeff);

  std::cout << "EXPECTED NUMBER OF EVENTS" << '\n';
  std::cout << "nbkgZ = " << intlumi * Zxsec * Zeff << '\n';
  std::cout << "nbkgttmm = " << intlumi * ttmmxsec * ttmmeff << '\n';
  std::cout << "nbkg = " << nbkg << '\n';
  std::cout << "nsig = " << nsig << '\n';

  Int_t nevents = nsig + nbkg;

  TH1F *h_tot_BDT = new TH1F("h_tot_BDT", "h_tot_BDT", h_H_BDT->GetNbinsX(),
                             h_H_BDT->GetBinLowEdge(1),
                             h_H_BDT->GetBinLowEdge(h_H_BDT->GetNbinsX() + 1));
  h_H_BDT->Scale(nsig);
  h_Z_BDT->Scale(nbkg);
  h_H_BDT->Sumw2();
  h_Z_BDT->Sumw2();
  h_tot_BDT->Add(h_H_BDT, h_Z_BDT);

  RooRealVar invmass("invmass", "BDT response", h_H_BDT->GetNbinsX(),
                     h_H_BDT->GetBinLowEdge(1),
                     h_H_BDT->GetBinLowEdge(h_H_BDT->GetNbinsX() + 1));

  RooDataHist *histo_data =
      new RooDataHist("data", "data", RooArgList(invmass), h_tot_BDT);

  RooRealVar Nsig("Nsig", "Nsig", nsig, 0, nevents);
  RooRealVar Nbkg("Nbkg", "Nbkg", nbkg, 0, nevents);

  RooDataHist *H_histo_data = new RooDataHist("H_histo_data", "H_histo_data",
                                              RooArgList(invmass), h_H_BDT);
  RooDataHist *Z_histo_data = new RooDataHist("Z_histo_data", "Z_histo_data",
                                              RooArgList(invmass), h_Z_BDT);

  RooHistPdf *Hhistpdf =
      new RooHistPdf("Hhistpdf", "Hhistpdf", invmass, *H_histo_data);
  RooHistPdf *Zhistpdf =
      new RooHistPdf("Zhistpdf", "Zhistpdf", invmass, *Z_histo_data);

  RooAddPdf model("model", "model", RooArgList(*Hhistpdf, *Zhistpdf),
                  RooArgList(Nsig, Nbkg));

  RooMCStudy *mcstudy_ex = new RooMCStudy(
      model, RooArgSet(invmass), RooFit::Binned(false), RooFit::Silence(),
      RooFit::Extended(true),
      RooFit::FitOptions(RooFit::Save(true), RooFit::PrintEvalErrors(0)));

  mcstudy_ex->generateAndFit(10, 0, kTRUE, 0);

  RooAbsData *single_toy_data = mcstudy_ex->genData(1);
  const RooFitResult *single_toy_res = mcstudy_ex->fitResult(1);
  RooArgSet single_toy_pars = single_toy_res->floatParsFinal();

  // plot example of mctoy fit
  RooAddPdf model_toy(
      "model_toy", "model_toy", RooArgList(*Zhistpdf, *Hhistpdf),
      RooArgList(static_cast<RooRealVar &>(single_toy_pars["Nbkg"]),
                 static_cast<RooRealVar &>(single_toy_pars["Nsig"])));
  model_toy.fitTo(*single_toy_data, RooFit::Save());
  TCanvas *can4 = new TCanvas("can4", "can4", 900, 600);
  RooPlot *frame4 = invmass.frame(Bins(h_tot_BDT->GetNbinsX()));
  single_toy_data->plotOn(frame4, Name("data"));
  model_toy.plotOn(frame4, Name("fit"));
  // model.plotOn(frame4, Components(RooArgSet(p_sig, p_bkg)), DrawOption("L"));
  model_toy.plotOn(frame4, Components(RooArgSet(*Hhistpdf)), DrawOption("L"),
                   LineColor(kRed), Name("higgs"));
  model_toy.plotOn(frame4, Components(RooArgSet(*Zhistpdf)), DrawOption("L"),
                   LineColor(kGreen), Name("background"));
  frame4->Draw();

  // TPaveText *pt1 = new TPaveText(8, 1200, 100, 1600); //FOR 1PI1PI
  TPaveText *pt1 = new TPaveText(-0.27, 265, -0.08, 400); // FOR ALLALL
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
  pt1->AddText("p_{T}(#tau) > 20 GeV/c");
  // pt1->AddText("|#eta(#tau)| < 2.5");
  pt1->Draw();

  auto legend = new TLegend(0.65, 0.65, 0.915, 0.85);
  legend->AddEntry("data", "Pseudo-data", "lep");
  legend->AddEntry("fit", "Fit", "l");
  legend->AddEntry("higgs", "H#rightarrow #tau_{h}#tau_{h}", "l");
  legend->AddEntry("background", "Background", "l");
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(kBlue, 0.);
  legend->SetTextSize(0.045);
  legend->Draw();

  frame4->GetYaxis()->SetTitle("Events/bin");
  frame4->GetYaxis()->SetTitleOffset(1.15);
  frame4->GetXaxis()->SetTitle("BDT response");
  frame4->GetXaxis()->SetTitleOffset(1.1);
  frame4->GetXaxis()->SetLabelOffset(0.013);
  frame4->GetXaxis()->SetLabelSize(0.04);
  frame4->GetYaxis()->SetLabelSize(0.04);
  frame4->GetXaxis()->SetTitleSize(0.045);
  frame4->GetYaxis()->SetTitleSize(0.045);
  frame4->SetTitle("");

  can4->SaveAs("roofit_plots/SingleMCtoy_TMVACTot.pdf");
  can4->Close();

  RooMCStudy *mcstudy = new RooMCStudy(
      model, RooArgSet(invmass), RooFit::Binned(false), RooFit::Silence(),
      RooFit::Extended(true),
      RooFit::FitOptions(RooFit::Save(true), RooFit::PrintEvalErrors(1)));

  mcstudy->generateAndFit(10000);

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
  can1->SaveAs("roofit_plots/Nsig_TMVACTot.png");
  TCanvas *can2 = new TCanvas("can2", "can2", 1200, 900);
  frame2->Draw();
  can2->SaveAs("roofit_plots/NsigError_TMVACTot.png");
  TCanvas *can3 = new TCanvas("can3", "can3", 1200, 900);
  frame3->Draw();
  can3->SaveAs("roofit_plots/NsigPull_TMVACTot.png");

  /*
TCanvas *can = new TCanvas("can", "can", 1200, 900);
RooPlot *xframe = invmass.frame();

histo_data->plotOn(xframe, Name("data"));
model.plotOn(xframe, Name("fit"));
// model.plotOn(xframe, Components(RooArgSet(p_sig,p_bkg)),
DrawOption("L"),
    // LineColor(kBlue)), Name("fit") ;
    model.plotOn(xframe, Components(RooArgSet(*Hhistpdf)), DrawOption("L"),
                 LineColor(kRed), Name("higgs"));
model.plotOn(xframe, Components(RooArgSet(*Zhistpdf)), DrawOption("L"),
             LineColor(kGreen), Name("background"));
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
legend->AddEntry("background", "Z#rightarrow #tau#tau", "l");
legend->Draw();

can->SaveAs("roofit_plots/HvsZ_TMVACTot_fit.png");
can->Close();*/

  // res->Print();

  return;
}

void fit_TMVATot() {
  const TString filename = "/eos/user/l/lvalla/MuColl/HvsZ/TMVACTot.root";
  fit_plotter(filename);
}
