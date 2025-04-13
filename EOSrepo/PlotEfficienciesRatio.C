#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

TH1F *ConvertTEfficiencyToTH1(const TEfficiency *eff) {
  if (!eff) {
    std::cerr << "Error: TEfficiency pointer is null!" << std::endl;
    return nullptr;
  }

  // Retrieve total and passed histograms
  const TH1 *totalHist = eff->GetTotalHistogram();
  const TH1 *passedHist = eff->GetPassedHistogram();

  if (!totalHist || !passedHist) {
    std::cerr << "Error: Failed to retrieve total or passed histograms!"
              << std::endl;
    return nullptr;
  }

  // Create a TH1F with the same binning as the total histogram
  TH1F *hist =
      new TH1F(Form("hist_%s", eff->GetName()), "Converted Histogram",
               totalHist->GetNbinsX(), totalHist->GetXaxis()->GetXmin(),
               totalHist->GetXaxis()->GetXmax());

  // Fill the TH1F with efficiency values
  for (int bin = 1; bin <= totalHist->GetNbinsX(); ++bin) {
    float efficiency = eff->GetEfficiency(bin);
    float errorLow = eff->GetEfficiencyErrorLow(bin);
    float errorUp = eff->GetEfficiencyErrorUp(bin);

    hist->SetBinContent(bin, efficiency);
    hist->SetBinError(bin, (errorLow + errorUp) /
                               2); // Set average error for simplicity
  }

  return hist;
}

void EfficiencyPlotter() {
  // Open the ROOT file
  TFile *inputFile = new TFile("efficiencies.root", "READ");

  if (!inputFile->IsOpen()) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return;
  }

  // Arrays of histogram names
  const char *efficiencyNames[] = {"eff_EtavsE0", "eff_EtavsE1", "eff_EtavsE2",
                                   "eff_EtavsE3", "eff_EtavsE4"};
  const int numEfficiencies =
      sizeof(efficiencyNames) / sizeof(efficiencyNames[0]);

  // Create a directory for storing output plots
  system("mkdir -p plots");

  for (int i = 0; i < numEfficiencies; ++i) {
    // Retrieve histograms from both directories

    TEfficiency *eff1 = (TEfficiency *)inputFile->Get(
        Form("tauguns/%s", efficiencyNames[i]));
    TEfficiency *eff2 = (TEfficiency *)inputFile->Get(
        Form("piguns/%s", efficiencyNames[i]));

    TH1F *hist1 = ConvertTEfficiencyToTH1(eff1);
    TH1F *hist2 = ConvertTEfficiencyToTH1(eff2);

    if (!hist1) {
      std::cerr
          << "Error: Could not retrieve histogram tauguns/"
          << efficiencyNames[i] << std::endl;
      continue;
    }

    if (!hist2) {
      std::cerr
          << "Error: Could not retrieve histogram piguns/"
          << efficiencyNames[i] << std::endl;
      continue;
    }

    // std::cout << hist1->GetNbinsX() << " vs " << hist2->GetNbinsX() << '\n';

    // Check histogram properties
    if (hist1->GetNbinsX() != hist2->GetNbinsX() ||
        hist1->GetXaxis()->GetXmin() != hist2->GetXaxis()->GetXmin() ||
        hist1->GetXaxis()->GetXmax() != hist2->GetXaxis()->GetXmax()) {
      std::cerr << "Error: Histograms " << efficiencyNames[i]
                << " have mismatched binning or ranges!" << std::endl;
      continue;
    }

    // Create a histogram to store the ratio
    TH1F *ratioHist =
        (TH1F *)hist1->Clone(Form("ratio_%s", efficiencyNames[i]));
    if (!ratioHist) {
      std::cerr << "Error: Could not clone histogram " << efficiencyNames[i]
                << std::endl;
      continue;
    }

    // Divide hist1 by hist2
    ratioHist->Divide(hist2);

    // Create a canvas
    TCanvas *c =
        new TCanvas(Form("c_ratio_%s", efficiencyNames[i]),
                    Form("Ratio for %s", efficiencyNames[i]), 1000, 800);
    gStyle->SetStripDecimals(0);
    gStyle->SetOptStat(0);

    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.06);

    // Customize the ratio histogram
    ratioHist->SetTitle(
        Form("; %s; #epsilon_{#tau}/#epsilon_{#pi}", "|#eta_{gen}^{(vis)}|"));
    ratioHist->SetLineColor(kBlue);
    ratioHist->SetMarkerColor(kBlue);
    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(1.2);
    ratioHist->SetMinimum(0.); // Adjust the range for better visibility
    ratioHist->SetMaximum(1.);

    ratioHist->GetXaxis()->SetTitleOffset(1.3);
    ratioHist->GetXaxis()->SetTitleSize(0.045);
    ratioHist->GetXaxis()->SetLabelSize(0.045);
    ratioHist->GetXaxis()->SetLabelOffset(0.012);
    ratioHist->GetYaxis()->SetTitleSize(0.045);
    ratioHist->GetYaxis()->SetLabelSize(0.045);
    ratioHist->GetYaxis()->SetLabelOffset(0.0125);

    // Draw the ratio histogram
    ratioHist->Draw("E1");

    if (efficiencyNames[i] == std::string("eff_EtavsE1")) {

      TPaveText *pt2 = new TPaveText(1.75, 0.82, 2.45, 0.96);
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

      TPaveText *pt3 = new TPaveText(0.15, 0.4, 1.2, 0.55);
      pt3->SetFillColor(0);
      pt3->SetMargin(0);
      pt3->SetFillStyle(0);
      pt3->SetBorderSize(0);
      pt3->SetTextAlign(11);
      // pt3->SetTextSize(0.035);
      //  pt3->SetTextFont(50);
      pt3->AddText("#bf{#tau_{h} 1-prong decays}");
      pt3->AddText("#bf{50 GeV #leq p_{T, gen}^{(vis)} #leq 100 GeV}");
      // pt3->AddText("#bf{}");
      pt3->Draw();
    }

    // Save the canvas as a PDF file
    c->SaveAs(Form("plots/ratio_%s.pdf", efficiencyNames[i]));

    // Cleanup
    delete c;
  }

  // Close the ROOT file
  inputFile->Close();
  delete inputFile;

  std::cout << "Ratio plots have been saved to the 'plots' directory."
            << std::endl;
  return;
}

void EfficiencyPlotter0N() {
  // Open the ROOT file
  TFile *inputFile = new TFile("efficiencies.root", "READ");

  if (!inputFile->IsOpen()) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return;
  }

  // Arrays of histogram names
  const char *efficiencyNames[] = {"eff_EtavsE0", "eff_EtavsE1", "eff_EtavsE2",
                                   "eff_EtavsE3", "eff_EtavsE4"};
  const int numEfficiencies =
      sizeof(efficiencyNames) / sizeof(efficiencyNames[0]);

  // Create a directory for storing output plots
  system("mkdir -p plots");

  for (int i = 0; i < numEfficiencies; ++i) {
    // Retrieve histograms from both directories

    TEfficiency *eff1 = (TEfficiency *)inputFile->Get(
        Form("tauguns_0N/%s", efficiencyNames[i]));
    TEfficiency *eff2 = (TEfficiency *)inputFile->Get(
        Form("piguns/%s", efficiencyNames[i]));

    TH1F *hist1 = ConvertTEfficiencyToTH1(eff1);
    TH1F *hist2 = ConvertTEfficiencyToTH1(eff2);

    if (!hist1) {
      std::cerr
          << "Error: Could not retrieve histogram tauguns/"
          << efficiencyNames[i] << std::endl;
      continue;
    }

    if (!hist2) {
      std::cerr
          << "Error: Could not retrieve histogram piguns/"
          << efficiencyNames[i] << std::endl;
      continue;
    }

    // std::cout << hist1->GetNbinsX() << " vs " << hist2->GetNbinsX() << '\n';

    // Check histogram properties
    if (hist1->GetNbinsX() != hist2->GetNbinsX() ||
        hist1->GetXaxis()->GetXmin() != hist2->GetXaxis()->GetXmin() ||
        hist1->GetXaxis()->GetXmax() != hist2->GetXaxis()->GetXmax()) {
      std::cerr << "Error: Histograms " << efficiencyNames[i]
                << " have mismatched binning or ranges!" << std::endl;
      continue;
    }

    // Create a histogram to store the ratio
    TH1F *ratioHist =
        (TH1F *)hist1->Clone(Form("ratio_%s", efficiencyNames[i]));
    if (!ratioHist) {
      std::cerr << "Error: Could not clone histogram " << efficiencyNames[i]
                << std::endl;
      continue;
    }

    // Divide hist1 by hist2
    ratioHist->Divide(hist2);

    // Create a canvas
    TCanvas *c =
        new TCanvas(Form("c_ratio_%s", efficiencyNames[i]),
                    Form("Ratio for %s", efficiencyNames[i]), 1000, 800);
    gStyle->SetStripDecimals(0);

    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.06);

    // Customize the ratio histogram
    ratioHist->SetTitle(
        Form("; %s; Ratio (tauguns/piguns)", efficiencyNames[i]));
    ratioHist->SetLineColor(kBlack);
    ratioHist->SetMarkerColor(kBlack);
    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(1.2);
    ratioHist->SetMinimum(0.5); // Adjust the range for better visibility
    ratioHist->SetMaximum(1.5);

    ratioHist->GetXaxis()->SetTitleOffset(1.2);
    ratioHist->GetXaxis()->SetTitleSize(0.045);
    ratioHist->GetXaxis()->SetLabelSize(0.045);
    ratioHist->GetXaxis()->SetLabelOffset(0.012);
    ratioHist->GetYaxis()->SetTitleSize(0.045);
    ratioHist->GetYaxis()->SetLabelSize(0.045);
    ratioHist->GetYaxis()->SetLabelOffset(0.0125);

    // Draw the ratio histogram
    ratioHist->Draw("E");

    // Add a legend
    TLegend *legend = new TLegend(0.66, 0.69, 0.96, 0.88);
    legend->AddEntry(hist1, "tauguns", "lep");
    legend->AddEntry(hist2, "piguns", "lep");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(42);
    legend->Draw();

    // Save the canvas as a PDF file
    c->SaveAs(Form("plots/ratio_%s_0N.pdf", efficiencyNames[i]));

    // Cleanup
    delete c;
  }

  // Close the ROOT file
  inputFile->Close();
  delete inputFile;

  std::cout << "Ratio plots have been saved to the 'plots' directory."
            << std::endl;
  return;
}

void EfficiencyvsTrackPlotter() {
  // Open the ROOT file
  TFile *inputFile = new TFile("efficiencies.root", "READ");

  if (!inputFile->IsOpen()) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return;
  }

  // Arrays of histogram names
  const char *efficiencyNames[] = {"eff_EtavsE0", "eff_EtavsE1", "eff_EtavsE2",
                                   "eff_EtavsE3", "eff_EtavsE4"};
  const int numEfficiencies =
      sizeof(efficiencyNames) / sizeof(efficiencyNames[0]);

  // Create a directory for storing output plots
  system("mkdir -p plots");

  for (int i = 0; i < numEfficiencies; ++i) {
    // Retrieve histograms from both directories

    TEfficiency *eff1 = (TEfficiency *)inputFile->Get(
        Form("piguns_notrackfiltering/%s", efficiencyNames[i]));
    TEfficiency *eff2 = (TEfficiency *)inputFile->Get(
        Form("piguns_tracks/%s", efficiencyNames[i]));

    TH1F *hist1 = ConvertTEfficiencyToTH1(eff1);
    TH1F *hist2 = ConvertTEfficiencyToTH1(eff2);

    if (!hist1) {
      std::cerr << "Error: Could not retrieve histogram piguns_tracks/"
                << efficiencyNames[i] << std::endl;
      continue;
    }

    if (!hist2) {
      std::cerr
          << "Error: Could not retrieve histogram piguns_notrackfiltering/"
          << efficiencyNames[i] << std::endl;
      continue;
    }

    // std::cout << hist1->GetNbinsX() << " vs " << hist2->GetNbinsX() << '\n';

    // Check histogram properties
    if (hist1->GetNbinsX() != hist2->GetNbinsX() ||
        hist1->GetXaxis()->GetXmin() != hist2->GetXaxis()->GetXmin() ||
        hist1->GetXaxis()->GetXmax() != hist2->GetXaxis()->GetXmax()) {
      std::cerr << "Error: Histograms " << efficiencyNames[i]
                << " have mismatched binning or ranges!" << std::endl;
      continue;
    }

    // Create a histogram to store the ratio
    TH1F *ratioHist =
        (TH1F *)hist1->Clone(Form("ratio_%s", efficiencyNames[i]));
    if (!ratioHist) {
      std::cerr << "Error: Could not clone histogram " << efficiencyNames[i]
                << std::endl;
      continue;
    }

    // Divide hist1 by hist2
    ratioHist->Divide(hist2);

    // Create a canvas
    TCanvas *c =
        new TCanvas(Form("c_ratio_%s", efficiencyNames[i]),
                    Form("Ratio for %s", efficiencyNames[i]), 1000, 800);
    gStyle->SetStripDecimals(0);
    gStyle->SetOptStat(0);

    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.06);

    // Customize the ratio histogram
    ratioHist->SetTitle(Form("; %s; #epsilon_{#pi}/#epsilon_{trk}",
                             "|#eta^{gen}(#pi)|"));
    ratioHist->SetLineColor(kBlue);
    ratioHist->SetMarkerColor(kBlue);
    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(1.2);
    ratioHist->SetMinimum(0.); // Adjust the range for better visibility
    ratioHist->SetMaximum(1.);

    ratioHist->GetXaxis()->SetTitleOffset(1.3);
    ratioHist->GetXaxis()->SetTitleSize(0.045);
    ratioHist->GetXaxis()->SetLabelSize(0.045);
    ratioHist->GetXaxis()->SetLabelOffset(0.012);
    ratioHist->GetYaxis()->SetTitleSize(0.045);
    ratioHist->GetYaxis()->SetLabelSize(0.045);
    ratioHist->GetYaxis()->SetLabelOffset(0.0125);

    // Draw the ratio histogram
    ratioHist->Draw("E1");

    if (efficiencyNames[i] == std::string("eff_EtavsE1")) {

      TPaveText *pt2 = new TPaveText(1.63, 0.76, 2.4, 0.94);
      pt2->SetFillColor(0);
      pt2->SetMargin(0);
      pt2->SetFillStyle(0);
      pt2->SetBorderSize(0);
      pt2->SetTextAlign(31);
      // pt2->SetTextFont(50);
      pt2->AddText("Muon Collider");
      pt2->AddText("#bf{#it{Simulation}} ");
      // pt2->AddText("  ");
      //  pt2->AddText("#bf{}");
      pt2->Draw();

      TPaveText *pt3 = new TPaveText(1.22, 0.22, 2.4, 0.3);
      pt3->SetFillColor(0);
      pt3->SetMargin(0);
      pt3->SetFillStyle(0);
      pt3->SetBorderSize(0);
      pt3->SetTextAlign(31);
      // pt3->SetTextSize(0.035);
      //  pt3->SetTextFont(50);
      // pt3->AddText("#bf{#tau_{h} 1-prong decays}");
      pt3->AddText("#bf{50 GeV #leq p_{T}^{gen}(#pi) #leq 100 GeV}");
      // pt3->AddText("#bf{}");
      pt3->Draw();
    }

    // Save the canvas as a PDF file
    c->SaveAs(Form("plots/ratio_track_%s.pdf", efficiencyNames[i]));

    // Cleanup
    delete c;
  }

  // Close the ROOT file
  inputFile->Close();
  delete inputFile;

  std::cout << "Ratio plots have been saved to the 'plots' directory."
            << std::endl;
  return;
}

int PlotEfficienciesRatio() {
  EfficiencyPlotter();
  //EfficiencyPlotter0N();

  EfficiencyvsTrackPlotter();
  return 0;
}