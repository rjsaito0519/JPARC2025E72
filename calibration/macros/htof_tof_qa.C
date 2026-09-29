// QA of the scattered-particle TOF (TPC helix x HTOF) from DstTPCHelixHTOF outputs.
// Tracks: non-beam, non-accidental, HTOF-matched, valid TOF. Read-only; writes one PDF.
//   root -l -b -q 'calibration/macros/htof_tof_qa.C+("<out.pdf>", "<DATA_DIR>", "02601,02602,...")'
// Pages: (1) dt_pi, m2, m2 vs p*q, dt_pi vs seg; (2) m2 by dE/dx tag, dt_p, dt_pi per seg,
//        dt_pi vs L_sec; (3) m2 vs p for q>0 / q<0, nsigma_pi, nsigma_p.
// Note: ACLiC ("+") writes htof_tof_qa_C.so/.d/.pcm next to this file (not git-ignored);
// e.g. call gSystem->SetBuildDir("/tmp/<you>/aclic", kTRUE) first to keep them out of the repo.
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <cstdio>
#include <vector>
#include <TObjArray.h>
#include <TObjString.h>

void htof_tof_qa(const char* out = "htof_tof_qa.pdf",
                 const char* data_dir = "/group/had/sks/Users/sryuta/JPARC2025E72",
                 const char* run_list = "02601,02602,02603,02604,02606,02607") {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetPadRightMargin(0.13);
    std::vector<TString> runs;
    { TString rl(run_list); TObjArray* tok = rl.Tokenize(","); for (auto* o : *tok) runs.push_back(((TObjString*)o)->GetString()); delete tok; }
    const char* D = data_dir;
    const double mpi2 = 0.13957*0.13957, mk2 = 0.493677*0.493677, mp2 = 0.938272*0.938272;

    TH1D* h_dtpi   = new TH1D("h_dtpi", "dt_{#pi} (pion-tagged);t_{sec} - t_{calc,#pi} [ns];tracks", 400, -10, 10);
    TH1D* h_m2     = new TH1D("h_m2", "m^{2} (all);m^{2} [GeV^{2}];tracks", 500, -0.5, 2.0);
    TH2D* h_m2p    = new TH2D("h_m2p", "m^{2} vs p #times charge;p #times q [GeV/c];m^{2} [GeV^{2}]", 200, -1.5, 1.5, 250, -0.5, 2.0);
    TH2D* h_dtseg  = new TH2D("h_dtseg", "dt_{#pi} vs HTOF seg (pion-tagged);HTOF seg;dt_{#pi} [ns]", 34, -0.5, 33.5, 200, -3, 3);
    TH1D* h_m2_pi  = new TH1D("h_m2_pi", "m^{2} by dE/dx tag;m^{2} [GeV^{2}];tracks", 500, -0.5, 2.0);
    TH1D* h_m2_k   = new TH1D("h_m2_k", "", 500, -0.5, 2.0);
    TH1D* h_m2_pr  = new TH1D("h_m2_pr", "", 500, -0.5, 2.0);
    TH1D* h_dtp    = new TH1D("h_dtp", "dt_{p} (proton-tagged);t_{sec} - t_{calc,p} [ns];tracks", 400, -10, 10);
    TH2D* h_dtL    = new TH2D("h_dtL", "dt_{#pi} vs L_{sec} (pion-tagged);L_{sec} [mm];dt_{#pi} [ns]", 100, 0, 1000, 200, -3, 3);
    TH2D* h_m2p_pos = new TH2D("h_m2p_pos", "m^{2} vs p (q > 0);p [GeV/c];m^{2} [GeV^{2}]", 150, 0, 1.5, 250, -0.5, 2.0);
    TH2D* h_m2p_neg = new TH2D("h_m2p_neg", "m^{2} vs p (q < 0);p [GeV/c];m^{2} [GeV^{2}]", 150, 0, 1.5, 250, -0.5, 2.0);
    TH1D* h_nspi   = new TH1D("h_nspi", "n#sigma_{#pi} (pion-tagged);n#sigma_{#pi};tracks", 200, -10, 10);
    TH1D* h_nsp    = new TH1D("h_nsp", "n#sigma_{p} (proton-tagged);n#sigma_{p};tracks", 200, -10, 10);
    long n_tr = 0;
    for (const auto& r : runs) {
        TFile f(Form("%s/run%s_TPCHelixHTOF.root", D, r.Data()));
        if (f.IsZombie()) { printf("Error: cannot open run%s_TPCHelixHTOF.root\n", r.Data()); return; }
        TTreeReader rd("htof", &f);
        TTreeReaderValue<std::vector<Int_t>> pid(rd, "pid"), seg_match(rd, "seg_match"), isb(rd, "is_beam"), isa(rd, "is_accidental"), chg(rd, "charge");
        TTreeReaderValue<std::vector<Double_t>> seg(rd, "extrap_seg"), dtpi(rd, "dt_pi"), dtp(rd, "dt_p"), m2(rd, "m2"),
            pv(rd, "p_vtx"), L(rd, "L_sec"), nspi(rd, "nsigma_pi"), nsp(rd, "nsigma_p");
        while (rd.Next()) {
            for (size_t i = 0; i < pid->size(); i++) {
                if ((*seg_match)[i] != 1 || (*isb)[i] || (*isa)[i] || TMath::IsNaN((*dtpi)[i])) continue;
                n_tr++;
                const int p = (*pid)[i], q = (*chg)[i];
                const bool is_pi = (p & 1) && !(p & 4), is_pr = (p & 4) && !(p & 1), is_k = (p & 2);
                const double m = (*m2)[i], mom = (*pv)[i];
                h_m2->Fill(m); h_m2p->Fill(mom*q, m);
                (q > 0 ? h_m2p_pos : h_m2p_neg)->Fill(mom, m);
                if (is_pi) { h_dtpi->Fill((*dtpi)[i]); h_dtseg->Fill((*seg)[i], (*dtpi)[i]); h_m2_pi->Fill(m); h_dtL->Fill((*L)[i], (*dtpi)[i]); h_nspi->Fill((*nspi)[i]); }
                if (is_pr) { h_dtp->Fill((*dtp)[i]); h_m2_pr->Fill(m); h_nsp->Fill((*nsp)[i]); }
                if (is_k)  h_m2_k->Fill(m);
            }
        }
    }
    auto mline = [](double x, double ymax, int col) { TLine* l = new TLine(x, 0, x, ymax); l->SetLineColor(col); l->SetLineStyle(2); l->Draw(); };
    auto gfit = [](TH1D* h, double w) { double c = h->GetBinCenter(h->GetMaximumBin()); TF1* g = new TF1(Form("g_%s", h->GetName()), "gaus", c - w, c + w);
        h->Fit(g, "0QR"); g->SetLineColor(kRed); g->Draw("same"); TLatex t; t.SetNDC(); t.SetTextSize(0.045);
        t.DrawLatex(0.15, 0.84, Form("core mean %+.3f ns, #sigma %.3f ns", g->GetParameter(1), std::fabs(g->GetParameter(2))));
        t.DrawLatex(0.15, 0.78, Form("entries %.0f", h->GetEntries())); };

    TCanvas c("c", "", 1600, 1200);
    c.Print(Form("%s[", out));
    // page 1
    c.Clear(); c.Divide(2, 2);
    c.cd(1); h_dtpi->Draw("HIST"); gfit(h_dtpi, 0.6);
    c.cd(2); gPad->SetLogy(); h_m2->Draw("HIST"); { double y = h_m2->GetMaximum(); mline(mpi2, y, kBlue); mline(mk2, y, kGreen+2); mline(mp2, y, kRed);
        TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.55, 0.84, "#color[4]{#pi}  #color[418]{K}  #color[2]{p}  (PDG m^{2})"); }
    c.cd(3); gPad->SetLogz(); h_m2p->Draw("colz");
    for (double m : {mpi2, mk2, mp2}) { TLine* l = new TLine(-1.5, m, 1.5, m); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
    c.cd(4); gPad->SetLogz(); h_dtseg->Draw("colz"); { TLine* l = new TLine(-0.5, 0, 33.5, 0); l->SetLineStyle(2); l->Draw(); }
    c.Print(out);
    // page 2
    c.Clear(); c.Divide(2, 2);
    c.cd(1); gPad->SetLogy();
    h_m2_pi->SetLineColor(kBlue); h_m2_k->SetLineColor(kGreen+2); h_m2_pr->SetLineColor(kRed);
    h_m2_pi->SetMaximum(2*std::max({h_m2_pi->GetMaximum(), h_m2_k->GetMaximum(), h_m2_pr->GetMaximum()}));
    h_m2_pi->Draw("HIST"); h_m2_k->Draw("HIST same"); h_m2_pr->Draw("HIST same");
    { TLegend* lg = new TLegend(0.55, 0.65, 0.85, 0.85); lg->AddEntry(h_m2_pi, "dE/dx #pi (no p bit)", "l"); lg->AddEntry(h_m2_k, "dE/dx K bit", "l"); lg->AddEntry(h_m2_pr, "dE/dx p (no #pi bit)", "l"); lg->Draw(); }
    c.cd(2); h_dtp->Draw("HIST"); gfit(h_dtp, 0.8);
    c.cd(3);
    {
        TGraphErrors* gm = new TGraphErrors(); TGraphErrors* gs = new TGraphErrors();
        for (int s = 0; s < 34; s++) {
            TH1D* pj = h_dtseg->ProjectionY(Form("pj%d", s), s + 1, s + 1);
            if (pj->GetEntries() < 50) { delete pj; continue; }
            double cc = pj->GetBinCenter(pj->GetMaximumBin()); TF1 g("gs", "gaus", cc - 0.6, cc + 0.6); pj->Fit(&g, "0QR");
            int k = gm->GetN(); gm->SetPoint(k, s, g.GetParameter(1)); gm->SetPointError(k, 0, g.GetParError(1));
            gs->SetPoint(k, s, std::fabs(g.GetParameter(2))); delete pj;
        }
        gm->SetTitle("pion dt_{#pi} per seg (black: core mean, red: #sigma);HTOF seg;[ns]");
        gm->SetMarkerStyle(20); gm->SetMinimum(-0.6); gm->SetMaximum(0.8); gm->Draw("AP");
        gs->SetMarkerStyle(24); gs->SetMarkerColor(kRed); gs->Draw("P same");
        TLine* l = new TLine(-0.5, 0, 33.5, 0); l->SetLineStyle(2); l->Draw();
    }
    c.cd(4); gPad->SetLogz(); h_dtL->Draw("colz"); { TLine* l = new TLine(0, 0, 1000, 0); l->SetLineStyle(2); l->Draw(); }
    c.Print(out);
    // page 3
    c.Clear(); c.Divide(2, 2);
    c.cd(1); gPad->SetLogz(); h_m2p_pos->Draw("colz");
    for (double m : {mpi2, mk2, mp2}) { TLine* l = new TLine(0, m, 1.5, m); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
    c.cd(2); gPad->SetLogz(); h_m2p_neg->Draw("colz");
    for (double m : {mpi2, mk2, mp2}) { TLine* l = new TLine(0, m, 1.5, m); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
    c.cd(3); h_nspi->Draw("HIST"); gfit(h_nspi, 1.5);
    c.cd(4); h_nsp->Draw("HIST"); gfit(h_nsp, 1.5);
    c.Print(out);
    c.Print(Form("%s]", out));
    printf("tracks used: %ld -> %s\n", n_tr, out);
}
