#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void PlotAndCompareEfficiencies() {
  // Open the ROOT file
  TFile *inputFile = new TFile("efficiencies.root", "READ");

  if (!inputFile->IsOpen()) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return;
  }

  // Arrays of TEfficiency objects' names for "EvsEta" and "EtavsE" categories
  const char *efficiencyEvsEta[] = {"eff_EvsEta0", "eff_EvsEta1", "eff_EvsEta2",
                                    "eff_EvsEta3", "eff_EvsEta4"};
  const char *efficiencyEtavsE[] = {"eff_EtavsE0", "eff_EtavsE1", "eff_EtavsE2",
                                    "eff_EtavsE3", "eff_EtavsE4", "eff_EtavsE5",
                                    "eff_EtavsE6"};

  // Number of histograms in each category
  const int numEvsEta = sizeof(efficiencyEvsEta) / sizeof(efficiencyEvsEta[0]);
  const int numEtavsE = sizeof(efficiencyEtavsE) / sizeof(efficiencyEtavsE[0]);

  // Create a directory for storing output plots
  system("mkdir -p plots");

  // Function to process and plot histograms
  auto processHistograms = [&](const char *dirName,
                               const char *efficiencyNames[],
                               int numEfficiencies, const char *xAxisTitle,
                               const char *label1, const char *label2) {
    for (int i = 0; i < numEfficiencies; ++i) {
      // Retrieve TEfficiency objects from both directories
      TEfficiency *eff1 = (TEfficiency *)inputFile->Get(
          Form("%s/%s", dirName, efficiencyNames[i]));
      TEfficiency *eff2 = (TEfficiency *)inputFile->Get(
          Form("piguns/%s", efficiencyNames[i]));

      if (!eff1 || !eff2) {
        std::cerr << "Error: Could not retrieve TEfficiency objects "
                  << efficiencyNames[i] << std::endl;
        continue;
      }

      // Create a canvas
      TCanvas *c =
          new TCanvas(Form("c_%s", efficiencyNames[i]),
                      Form("Canvas for %s", efficiencyNames[i]), 1000, 800);
      gStyle->SetStripDecimals(0);

      c->SetLeftMargin(0.12);
      c->SetBottomMargin(0.13);
      c->SetRightMargin(0.06);

      // Removed the grid
      c->SetGrid(0);

      // Draw the first efficiency
      eff1->Draw("AP");
      // Draw the second efficiency on the same canvas
      eff2->Draw("P SAME");

      // Set X-axis title and Y-axis title
      eff1->SetTitle(
          Form("; %s;#pi  reconstruction efficiency", xAxisTitle)); // No title for the plot

      gPad->Update();

      TGraph *graph1 = eff1->GetPaintedGraph();
      graph1->SetMinimum(0.);
      graph1->SetMaximum(1.);
      TGraph *graph2 = eff2->GetPaintedGraph();
      graph2->SetMinimum(0.);
      graph2->SetMaximum(1.);

      graph1->GetXaxis()->SetTitleOffset(1.2);
      graph1->GetXaxis()->SetTitleSize(0.045);
      graph1->GetXaxis()->SetLabelSize(0.045);
      graph1->GetXaxis()->SetLabelOffset(0.012);
      graph1->GetYaxis()->SetTitleSize(0.045);
      graph1->GetYaxis()->SetLabelSize(0.045);
      graph1->GetYaxis()->SetLabelOffset(0.0125);

      eff1->SetLineColor(kRed);
      eff1->SetMarkerColor(kRed);
      eff1->SetMarkerStyle(20);
      eff1->SetMarkerSize(1.2);

      eff2->SetLineColor(kBlue);
      eff2->SetMarkerColor(kBlue);
      eff2->SetMarkerStyle(21);
      eff2->SetMarkerSize(1.2);

      // Set Y-axis range from 0.0 to 1.0
      gPad->SetTicks(); // Add ticks to the axes
      gPad->Update();   // Update the pad to apply the changes
      gPad->SetLogy(0); // Ensure Y-axis is not in log scale

      if (TString(efficiencyNames[i]) == TString("eff_EtavsE1")) {

        eff1->SetTitle(
          Form(";|#eta^{gen}(#pi)|;#pi  reconstruction efficiency"));

        TLegend *legend = new TLegend(0.66, 0.69, 0.96, 0.88);
        legend->AddEntry(eff2, label2, "lep");
        legend->AddEntry(eff1, label1, "lep");
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0); // Solid fill style
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);

        // legend->SetMargin(0.1);  // Reduce the space around legend
        legend->SetEntrySeparation(0.03); // Reduce the space between entries

        legend->Draw();

        TPaveText *pt2 = new TPaveText(0.65, 0.83, 1.3, 0.94);
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

        TPaveText *pt3 = new TPaveText(1.55, 0.63, 2.63, 0.75);
        pt3->SetFillColor(0);
        pt3->SetMargin(0);
        pt3->SetFillStyle(0);
        pt3->SetBorderSize(0);
        pt3->SetTextAlign(31);
        // pt3->SetTextSize(0.035);
        //  pt3->SetTextFont(50);
        //pt3->AddText("#bf{#tau_{h} 1-prong decays}");
        pt3->AddText("#bf{50 GeV #leq p_{T}^{gen} (#pi) #leq 100 GeV}");
        // pt3->AddText("#bf{}");
        pt3->Draw();
      }

      if (TString(efficiencyNames[i]) == TString("eff_EvsEta0")) {

        eff1->SetTitle(
          Form(";p_{T}^{gen}(#pi) [GeV/c];#pi  reconstruction efficiency"));

        gPad->SetLogx();

        graph1->SetMinimum(0.);
        graph1->SetMaximum(1.4);
        graph2->SetMinimum(0.);
        graph2->SetMaximum(1.4);
        gPad->Update();

        TLegend *legend = new TLegend(0.62, 0.67, 0.95, 0.87);
        legend->AddEntry(eff2, label2, "lep");
        legend->AddEntry(eff1, label1, "lep");
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0); // Solid fill style
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);

        // legend->SetMargin(0.1);  // Reduce the space around legend
        legend->SetEntrySeparation(0.03); // Reduce the space between entries

        legend->Draw();

        graph1->GetXaxis()->SetLabelOffset(0.006);

        TPaveText *pt4 = new TPaveText(6, 1.03, 27, 1.33);
        pt4->SetFillColor(0);
        pt4->SetMargin(0);
        pt4->SetFillStyle(0);
        pt4->SetBorderSize(0);
        pt4->SetTextAlign(11);
        pt4->SetTextSize(0.043);
        pt4->AddText("Muon Collider");
        pt4->AddText("#bf{#it{Simulation}}");
        //pt4->AddText("  ");
        //pt4->AddText("#bf{#pi_{h} 1-prong decays}");
        pt4->AddText("#bf{0 #leq |#eta^{gen}(#pi)| #leq 0.5}");
        // pt4->AddText("#bf{}");
        pt4->Draw();
        gPad->Modified();
        gPad->Update();
      }

      gPad->SetGridx(0);
      gPad->SetGridy(0);

      // Save the canvas as a PDF file
      c->SaveAs(Form("plots/%s_comparison.pdf", efficiencyNames[i]));

      // Cleanup
      delete c;
    }
  };

  // Process "EvsEta" histograms with specific X-axis title and legend labels
  processHistograms("piguns", efficiencyEvsEta, numEvsEta,
                    "p_{T, mc}", "250 keV HCAL", "2 MeV HCAL");

  // Process "EtavsE" histograms with specific X-axis title and legend labels
  processHistograms("piguns", efficiencyEtavsE, numEtavsE,
                    "|#eta_{mc}|", "250 keV HCAL", "2 MeV HCAL");

  // Close the ROOT file
  inputFile->Close();
  delete inputFile;

  std::cout << "Plots have been saved to the 'plots' directory." << std::endl;
}

int main() {
  PlotAndCompareEfficiencies();
  return 0;
}
