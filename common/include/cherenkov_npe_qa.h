#ifndef CHERENKOV_NPE_QA_H
#define CHERENKOV_NPE_QA_H

#include <algorithm>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>

namespace ana_helper {
//_____________________________________________________________________________
inline void
append_cherenkov_npe_page(
    TFile* input,
    TCanvas* canvas,
    const TString& pdf,
    const std::vector<TString>& names,
    const TString& particle,
    int columns,
    int rows
)
{
    canvas->cd();
    canvas->Clear();
    canvas->Divide(columns, rows);
    for (std::size_t index = 0; index < names.size(); ++index) {
        canvas->cd(index+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(.13);
        gPad->SetRightMargin(.04);
        gPad->SetTopMargin(.10);
        gPad->SetBottomMargin(.13);
        const TString suffix = particle.IsNull() ? "" : "_"+particle;
        const std::vector<TString> keys{names[index], names[index]+"_AwT", names[index]+suffix, names[index]+"_AwT"+suffix};
        const std::vector<TString> labels{"No PID cut", "No PID + AwT", particle, particle+" + AwT"};
        const int colors[]{kGray+2, kGreen+2, kBlue, kOrange};
        double maximum = 0.;
        for (const auto& key : keys) {
            auto* histogram = dynamic_cast<TH1*>(input->Get(key));
            if (histogram) maximum = std::max(maximum, histogram->GetMaximum());
        }
        auto* legend = new TLegend(.38, .66, .95, .90);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(.030);
        bool drawn = false;
        for (std::size_t selection = 0; selection < keys.size(); ++selection) {
            if (selection >= 2 && suffix.IsNull()) continue;
            auto* histogram = dynamic_cast<TH1*>(input->Get(keys[selection]));
            if (!histogram) {
                legend->AddEntry(static_cast<TObject*>(nullptr), labels[selection]+": missing", "");
                std::cerr << "[Warning] NPE QA histogram missing: " << keys[selection] << std::endl;
                continue;
            }
            // Draw a copy so the stored spectrum and its binning stay intact.
            auto* displayed = histogram->DrawCopy(drawn ? "hist same" : "hist");
            displayed->SetTitle(names[index]+" ("+particle+")");
            displayed->GetXaxis()->SetTitle("NPE");
            displayed->GetYaxis()->SetTitle("Counts");
            displayed->SetLineColor(colors[selection]);
            displayed->SetLineWidth(selection == 0 ? 1 : 2);
            displayed->SetStats(false);
            displayed->SetMinimum(.5);
            displayed->SetMaximum(std::max(2., maximum*2.));
            legend->AddEntry(displayed, Form("%s (N=%.0f)", labels[selection].Data(), histogram->GetEntries()), "l");
            drawn = true;
        }
        if (drawn) {
            legend->Draw();
        } else {
            delete legend;
            auto* note = new TLatex();
            note->SetNDC();
            note->SetTextSize(.045);
            note->DrawLatex(.1, .6, names[index]);
            note->DrawLatex(.1, .45, "NPE histogram is not in the input ROOT");
        }
    }
    canvas->Print(pdf);
}
}
#endif
