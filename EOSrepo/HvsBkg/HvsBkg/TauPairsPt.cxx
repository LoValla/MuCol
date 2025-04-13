#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TGraphErrors.h>

void display_overflow(TH1 *h) {
  h->GetXaxis()->SetRange(1, h->GetNbinsX() + 1);
}

void TauPairsPt(const char *filename = "HvsAll_Histograms.root") {

  gSystem->Exec("mkdir -p plots");

  // Open the file
  TFile *file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    printf("Error: cannot open file %s\n", filename);
    return;
  }

  // List of histogram names
  const char *hist_names[] = {
      "h_recPtSingle"// Add all H histograms here
  };

  TList *HList = (TList *)file->Get("H_had+had_histograms");
  TList *ZList = (TList *)file->Get("Z_had+had_histograms");
  TList *MuMuList = (TList *)file->Get("MuMu_had+had_histograms");

  // Number of histograms
  const int num_hist = sizeof(hist_names) / sizeof(hist_names[0]);

  // Loop over histograms
  for (int i = 0; i < num_hist; i++) {
    // Create a canvas for each histogram
    TCanvas *c1 = new TCanvas(Form("c1_%d", i), hist_names[i], 800, 600);

    // Get the histograms from the file for H, Z, and MuMu
    TH1D *h_hist = static_cast<TH1D *>(HList->FindObject(hist_names[i]));
    TH1D *z_hist = static_cast<TH1D *>(ZList->FindObject(hist_names[i]));
    TH1D *mumu_hist = static_cast<TH1D *>(MuMuList->FindObject(hist_names[i]));

    if (!h_hist || !z_hist || !mumu_hist) {
      printf("Error: cannot find histograms for %s\n", hist_names[i]);
      continue;
    }

    display_overflow(h_hist);
    display_overflow(z_hist);
    display_overflow(mumu_hist);

    gStyle->SetOptStat(0000);
    gStyle->SetStripDecimals(0);

    // Style settings
    h_hist->SetLineColor(kRed);
    z_hist->SetLineColor(kBlue);
    mumu_hist->SetLineColor(kGreen);

    // Determine maximum value for proper scaling
    double max_value = std::max(
        {h_hist->GetMaximum(), z_hist->GetMaximum(), mumu_hist->GetMaximum()});
    h_hist->SetMaximum(
        1.2 * max_value); // Set maximum Y to slightly above the highest peak

    h_hist->GetYaxis()->SetTitle("Probability density");
    //h_hist->GetXaxis()->SetTitle("");

    // Draw histograms with "axis" option first
    h_hist->Draw("axis");
    z_hist->Draw("axis same");
    mumu_hist->Draw("axis same");

    // Create TGraphErrors for each histogram to add horizontal error bars
    TGraphErrors *h_graph = new TGraphErrors(h_hist->GetNbinsX()+1);
    TGraphErrors *z_graph = new TGraphErrors(z_hist->GetNbinsX()+1);
    TGraphErrors *mumu_graph = new TGraphErrors(mumu_hist->GetNbinsX()+1);

    for (int j = 0; j < h_hist->GetNbinsX()+1; ++j) {
      double x = h_hist->GetBinCenter(j + 1);
      double y = h_hist->GetBinContent(j + 1);
      double xerr = h_hist->GetBinWidth(j + 1) / 2.0;
      double yerr = 0.0; // No vertical error

      h_graph->SetPoint(j, x, y);
      h_graph->SetPointError(j, xerr, yerr);

      z_graph->SetPoint(j, x, z_hist->GetBinContent(j + 1));
      z_graph->SetPointError(j, xerr, yerr);

      mumu_graph->SetPoint(j, x, mumu_hist->GetBinContent(j + 1));
      mumu_graph->SetPointError(j, xerr, yerr);
    }

    // Set styles for the graphs
    h_graph->SetLineColor(kRed);
    h_graph->SetMarkerColor(kRed);
    h_graph->SetMarkerStyle(20);
    h_graph->SetMarkerSize(0.7);
    z_graph->SetLineColor(kBlue);
    z_graph->SetMarkerColor(kBlue);
    z_graph->SetMarkerStyle(20);
    z_graph->SetMarkerSize(0.7);
    mumu_graph->SetLineColor(kGreen);
    mumu_graph->SetMarkerColor(kGreen);
    mumu_graph->SetMarkerStyle(20);
    mumu_graph->SetMarkerSize(0.7);

    // Draw the graphs with horizontal error bars
    h_graph->Draw("P same");
    z_graph->Draw("P same");
    mumu_graph->Draw("P same");

    // Customize the axes
    h_hist->GetXaxis()->SetTitleOffset(1.3);
    h_hist->GetXaxis()->SetTitleSize(0.045);
    h_hist->GetXaxis()->SetLabelSize(0.045);
    h_hist->GetXaxis()->SetLabelOffset(0.012);
    h_hist->GetYaxis()->SetTitleSize(0.045);
    h_hist->GetYaxis()->SetTitleOffset(1.35);
    h_hist->GetYaxis()->SetLabelSize(0.045);
    h_hist->GetYaxis()->SetLabelOffset(0.0125);

    c1->SetLeftMargin(0.13); // Increase left margin for better space
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.04);

    // Create and draw the legend
    TLegend *leg = new TLegend(0.71, 0.7, 0.95, 0.88); // Top-right corner
    leg->SetBorderSize(1);
    leg->SetFillStyle(0); // Remove fill color
    leg->SetTextSize(0.04); // Increase text size
    leg->AddEntry(h_graph, "H #rightarrow #tau_{h}#tau_{h} (VBF)", "lep");
    leg->AddEntry(z_graph, "Z #rightarrow #tau_{h}#tau_{h} (VBF)", "lep");
    leg->AddEntry(mumu_graph, "#tau#tau #rightarrow #tau_{h}#tau_{h}#mu#mu", "lep");
    leg->Draw();

    // Create and draw TPaveText with increased size and adjusted position
    TPaveText *pt1 = new TPaveText(0.33, 0.73, 0.75, 0.86, "NDC");
    pt1->SetFillColor(0);
    pt1->SetMargin(0.1); // Increase margin for better spacing
    pt1->SetFillStyle(0);
    pt1->SetBorderSize(0);
    pt1->SetTextAlign(11);
    pt1->SetTextSize(0.043); // Adjust text size
    pt1->AddText("Tau pairs distributions");
    pt1->AddText("#bf{#it{Simulation}} ");
    pt1->AddText("#bf{#sqrt{s} = 3 TeV} ");
    pt1->Draw();

    TPaveText *pt2 = new TPaveText(0.71, 0.6, 0.9, 0.7, "NDC"); // Positioned under the legend
    pt2->SetFillColor(0);
    pt2->SetMargin(0.1);
    pt2->SetFillStyle(0);
    pt2->SetBorderSize(0);
    pt2->SetTextAlign(11);
    pt2->SetTextSize(0.042); // Adjust text size
    pt2->AddText("#bf{p_{T}(#tau) > 20 GeV/c}");
    //pt2->AddText("#bf{|#eta| < 2.5}");
    pt2->Draw();

    // Save each canvas as a PDF file in the "plots" directory
    TString output_filename = Form("plots/%s.pdf", hist_names[i]);
    c1->SaveAs(output_filename);

    // Clean up
    delete c1;
    delete h_graph;
    delete z_graph;
    delete mumu_graph;
  }

  // Close the file
  file->Close();
  delete file;
}
