// HTOF beam-window check with field-off BcOut tracks: beam profile at the HTOF plane-0 z and
// whether plane-0 HTOF (seg 0-5) fired. Used to check the beam-window rectangle assumed by
// TPCAnalyzer::ExtrapolateToHTOF (analyzer src/TPCAnalyzer.cc: TPC frame |x| < 71 mm,
// |y - 4| < 56 mm on plane 0). Read-only; writes one PDF.
//   root -l -b -q 'calibration/macros/htof_beam_window.C+("<out.pdf>", "<DATA_DIR>", "02569,02570")'
// Inputs per run: <DATA_DIR>/run<run>_BcOut.root (tree bcout: ntrack, chisqr, x0, y0, u0, v0) and
// <DATA_DIR>/run<run>_Hodo.root (tree hodo: htof_raw_seg, htof_tdc_u/d, htof_hit_seg), same entries.
// Track: the smallest-chisqr BcOut track of the event (any ntrack >= 1), straight line to z_htof
// (BcOut local z = global z in DCGEO: L = Z for BLC2), x = x0 + u0 z, y = y0 + v0 z.
// "Fired" (raw): any U or D TDC of the segment inside [tdc_lo, tdc_hi]; "fired" (hit): the segment
// is in htof_hit_seg (analyzer hit). Transmission = fraction of tracks with no plane-0 segment fired.
// Pages: (1) chisqr, raw TDC of seg 0-5 with the gate, fired-segment counts; (2) track maps: all,
// no plane-0 fired, plane-0 fired; (3) per-segment fired maps seg 0-5; (4) transmission map and its
// x / y slices fitted with a smoothed box (edges x1, x2 and resolution sigma), code rectangle overlaid;
// (5) beam profile (all tracks) x / y slices (log) fitted with (narrow + wide Gaussian) x smoothed box:
// the edges where the beam tails are cut (what constrains the window size if the beam tails are cut there).
// Assumes the TPC frame x/y used by the code equal the global x/y (HypTPC X = Y = 0 in DCGEO);
// check this before comparing with the code values.
// Note: ACLiC ("+") writes htof_beam_window_C.so/.d/.pcm next to this file (not git-ignored);
// e.g. call gSystem->SetBuildDir("/tmp/<you>/aclic", kTRUE) first to keep them out of the repo.
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <cstdio>
#include <vector>

namespace {
// smoothed box: plateau p0 between edges p1 < p2, Gaussian edge width p3, constant p4
Double_t smooth_box(Double_t* x, Double_t* p) {
    const Double_t s = std::max(std::fabs(p[3]), 1e-3) * TMath::Sqrt2();
    return 0.5 * p[0] * (TMath::Erf((x[0] - p[1]) / s) - TMath::Erf((x[0] - p[2]) / s)) + p[4];
}

// beam profile cut by an aperture: (narrow + wide Gaussian) x smoothed box (unit plateau) + constant
// p0 amp1, p1 mean, p2 sigma1, p3 amp2, p4 sigma2, p5 edge_lo, p6 edge_hi, p7 edge sigma, p8 const
Double_t cut_profile(Double_t* x, Double_t* p) {
    Double_t box_par[5] = {1.0, p[5], p[6], p[7], 0.0};
    const Double_t g = p[0] * TMath::Gaus(x[0], p[1], p[2]) + p[3] * TMath::Gaus(x[0], p[1], p[4]);
    return g * smooth_box(x, box_par) + p[8];
}

// fit a transmission slice with the smoothed box; returns the TF1 (drawn by the caller)
TF1* fit_box(TH1D* h, Double_t lo_guess, Double_t hi_guess, Double_t range_lo, Double_t range_hi) {
    TF1* f = new TF1(Form("f_%s", h->GetName()), smooth_box, range_lo, range_hi, 5);
    f->SetParNames("plateau", "edge_lo", "edge_hi", "sigma", "const");
    f->SetParameters(0.9, lo_guess, hi_guess, 3.0, 0.02);
    f->SetParLimits(0, 0.0, 1.2);
    f->SetParLimits(3, 0.1, 30.0);
    f->SetParLimits(4, 0.0, 0.5);
    f->SetLineColor(kRed);
    h->Fit(f, "0QR");
    return f;
}
} // namespace

void htof_beam_window(const char* out = "htof_beam_window.pdf",
                      const char* data_dir = "/group/had/sks/Users/sryuta/JPARC2025E72",
                      const char* run_list = "02569,02570",
                      Double_t z_htof = -499.7,       // global z of HTOF plane 0 [mm]
                      Double_t tdc_lo = 700000.,      // raw TDC gate (check page 1)
                      Double_t tdc_hi = 730000.,
                      Double_t chisqr_max = 1e9) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetPadRightMargin(0.13);
    std::vector<TString> runs;
    { TString rl(run_list); TObjArray* tok = rl.Tokenize(","); for (auto* o : *tok) runs.push_back(((TObjString*)o)->GetString()); delete tok; }

    // code rectangle (TPC frame, see header)
    const Double_t code_x_half = 71.0, code_y_center = 4.0, code_y_half = 56.0;
    const Int_t n_seg_p0 = 6;
    const Int_t nb = 150; const Double_t xy_lim = 150.;

    TH1D* h_chi   = new TH1D("h_chi", "BcOut best-track #chi^{2};#chi^{2};events", 200, 0, 20);
    TH1D* h_ntr   = new TH1D("h_ntr", "BcOut ntrack;ntrack;events", 10, -0.5, 9.5);
    TH1D* h_tdc   = new TH1D("h_tdc", "HTOF raw TDC (U and D), seg 0-5;TDC;hits", 400, 600000, 800000);
    TH1D* h_nfire = new TH1D("h_nfire", "fired plane-0 segments per track (raw);n seg;tracks", 7, -0.5, 6.5);
    TH2D* h_all   = new TH2D("h_all", Form("all tracks at z = %.1f;x [mm];y [mm]", z_htof), nb, -xy_lim, xy_lim, nb, -xy_lim, xy_lim);
    TH2D* h_none  = new TH2D("h_none", "no plane-0 segment fired (raw);x [mm];y [mm]", nb, -xy_lim, xy_lim, nb, -xy_lim, xy_lim);
    TH2D* h_fire  = new TH2D("h_fire", "any plane-0 segment fired (raw);x [mm];y [mm]", nb, -xy_lim, xy_lim, nb, -xy_lim, xy_lim);
    TH2D* h_none_hit = new TH2D("h_none_hit", "no plane-0 segment in htof_hit_seg;x [mm];y [mm]", nb, -xy_lim, xy_lim, nb, -xy_lim, xy_lim);
    TH2D* h_seg[n_seg_p0];
    for (Int_t s = 0; s < n_seg_p0; s++)
        h_seg[s] = new TH2D(Form("h_seg%d", s), Form("seg%d fired (raw);x [mm];y [mm]", s), nb / 2, -xy_lim, xy_lim, nb / 2, -xy_lim, xy_lim);

    Long64_t n_ev = 0, n_trk = 0, n_mismatch = 0;
    for (const auto& r : runs) {
        TFile fb(Form("%s/run%s_BcOut.root", data_dir, r.Data()));
        TFile fh(Form("%s/run%s_Hodo.root", data_dir, r.Data()));
        if (fb.IsZombie() || fh.IsZombie()) { printf("Error: cannot open run %s inputs\n", r.Data()); return; }
        TTreeReader rb("bcout", &fb), rh("hodo", &fh);
        TTreeReaderValue<UInt_t> ev_b(rb, "event_number"), ev_h(rh, "event_number");
        TTreeReaderValue<Int_t> ntrack(rb, "ntrack");
        TTreeReaderValue<std::vector<Double_t>> chi(rb, "chisqr"), x0(rb, "x0"), y0(rb, "y0"), u0(rb, "u0"), v0(rb, "v0");
        TTreeReaderValue<std::vector<Double_t>> raw_seg(rh, "htof_raw_seg"), hit_seg(rh, "htof_hit_seg");
        TTreeReaderValue<std::vector<std::vector<Double_t>>> tdc_u(rh, "htof_tdc_u"), tdc_d(rh, "htof_tdc_d");
        while (rb.Next()) {
            if (!rh.Next()) { printf("Error: hodo tree shorter than bcout in run %s\n", r.Data()); break; }
            n_ev++;
            if (*ev_b != *ev_h) { n_mismatch++; continue; }
            h_ntr->Fill(*ntrack);
            // best track
            Int_t ib = -1;
            for (std::size_t i = 0; i < chi->size(); i++)
                if (ib < 0 || (*chi)[i] < (*chi)[ib]) ib = static_cast<Int_t>(i);
            if (ib < 0) continue;
            h_chi->Fill((*chi)[ib]);
            if ((*chi)[ib] > chisqr_max) continue;
            n_trk++;
            const Double_t x = (*x0)[ib] + (*u0)[ib] * z_htof, y = (*y0)[ib] + (*v0)[ib] * z_htof;
            // plane-0 fired segments
            Bool_t fired[n_seg_p0] = {};
            for (std::size_t j = 0; j < raw_seg->size(); j++) {
                const Int_t s = TMath::Nint((*raw_seg)[j]);
                if (s < 0 || s >= n_seg_p0) continue;
                for (const auto* v : {&(*tdc_u)[j], &(*tdc_d)[j]})
                    for (Double_t t : *v) {
                        h_tdc->Fill(t);
                        if (t >= tdc_lo && t <= tdc_hi) fired[s] = kTRUE;
                    }
            }
            Bool_t fired_hit = kFALSE;
            for (Double_t s : *hit_seg) if (TMath::Nint(s) >= 0 && TMath::Nint(s) < n_seg_p0) fired_hit = kTRUE;
            Int_t nf = 0;
            for (Int_t s = 0; s < n_seg_p0; s++) if (fired[s]) { nf++; h_seg[s]->Fill(x, y); }
            h_nfire->Fill(nf);
            h_all->Fill(x, y);
            (nf ? h_fire : h_none)->Fill(x, y);
            if (!fired_hit) h_none_hit->Fill(x, y);
        }
    }
    printf("events %lld, event_number mismatches %lld, tracks used %lld\n", n_ev, n_mismatch, n_trk);

    // transmission maps and slices
    TH2D* h_tr = (TH2D*)h_none->Clone("h_tr");
    h_tr->Divide(h_none, h_all, 1, 1, "B");
    h_tr->SetTitle("transmission (no plane-0 fired, raw) / all;x [mm];y [mm]");
    TH2D* h_tr_hit = (TH2D*)h_none_hit->Clone("h_tr_hit");
    h_tr_hit->Divide(h_none_hit, h_all, 1, 1, "B");
    // x slice: |y - y_c| < 30 mm; y slice: |x| < 40 mm (inside the code window)
    auto slice = [&](TH2D* num, const char* name, Bool_t along_x, Double_t lo, Double_t hi) {
        TH1D* n = along_x ? num->ProjectionX(Form("%s_n", name), num->GetYaxis()->FindBin(lo), num->GetYaxis()->FindBin(hi))
                          : num->ProjectionY(Form("%s_n", name), num->GetXaxis()->FindBin(lo), num->GetXaxis()->FindBin(hi));
        TH1D* d = along_x ? h_all->ProjectionX(Form("%s_d", name), h_all->GetYaxis()->FindBin(lo), h_all->GetYaxis()->FindBin(hi))
                          : h_all->ProjectionY(Form("%s_d", name), h_all->GetXaxis()->FindBin(lo), h_all->GetXaxis()->FindBin(hi));
        TH1D* t = (TH1D*)n->Clone(name);
        t->Divide(n, d, 1, 1, "B");
        for (Int_t b = 1; b <= t->GetNbinsX(); b++) if (d->GetBinContent(b) < 5) { t->SetBinContent(b, 0); t->SetBinError(b, 0); }
        t->SetTitle(Form("transmission, %s;%s [mm];fraction", along_x ? Form("%.0f < y < %.0f", lo, hi) : Form("%.0f < x < %.0f", lo, hi), along_x ? "x" : "y"));
        return t;
    };
    const Double_t yc = code_y_center;
    TH1D* t_x     = slice(h_none, "t_x", kTRUE, yc - 30., yc + 30.);
    TH1D* t_y     = slice(h_none, "t_y", kFALSE, -40., 40.);
    TH1D* t_x_hit = slice(h_none_hit, "t_x_hit", kTRUE, yc - 30., yc + 30.);
    TH1D* t_y_hit = slice(h_none_hit, "t_y_hit", kFALSE, -40., 40.);
    TF1* f_x     = fit_box(t_x, -code_x_half, code_x_half, -xy_lim, xy_lim);
    TF1* f_y     = fit_box(t_y, yc - code_y_half, yc + code_y_half, -xy_lim, xy_lim);
    TF1* f_x_hit = fit_box(t_x_hit, -code_x_half, code_x_half, -xy_lim, xy_lim);
    TF1* f_y_hit = fit_box(t_y_hit, yc - code_y_half, yc + code_y_half, -xy_lim, xy_lim);
    auto report = [](const char* tag, TF1* fx, TF1* fy) {
        printf("%-4s x: edges %7.1f %7.1f (+-%.1f %.1f)  center %6.1f  width %6.1f  sigma %.1f | "
               "y: edges %7.1f %7.1f (+-%.1f %.1f)  center %6.1f  height %6.1f  sigma %.1f\n", tag,
               fx->GetParameter(1), fx->GetParameter(2), fx->GetParError(1), fx->GetParError(2),
               0.5 * (fx->GetParameter(1) + fx->GetParameter(2)), fx->GetParameter(2) - fx->GetParameter(1), std::fabs(fx->GetParameter(3)),
               fy->GetParameter(1), fy->GetParameter(2), fy->GetParError(1), fy->GetParError(2),
               0.5 * (fy->GetParameter(1) + fy->GetParameter(2)), fy->GetParameter(2) - fy->GetParameter(1), std::fabs(fy->GetParameter(3)));
    };
    // beam profile slices (counts): Gaussian core/tails cut by a smoothed box (edges from the tail cut)
    TH1D* p_x = h_all->ProjectionX("p_x", h_all->GetYaxis()->FindBin(yc - 30.), h_all->GetYaxis()->FindBin(yc + 30.));
    TH1D* p_y = h_all->ProjectionY("p_y", h_all->GetXaxis()->FindBin(-40.), h_all->GetXaxis()->FindBin(40.));
    p_x->SetTitle(Form("beam profile, %.0f < y < %.0f;x [mm];tracks", yc - 30., yc + 30.));
    p_y->SetTitle("beam profile, -40 < x < 40;y [mm];tracks");
    auto fit_profile = [&](TH1D* h, Double_t lo, Double_t hi) {
        TF1* f = new TF1(Form("f_%s", h->GetName()), cut_profile, -xy_lim, xy_lim, 9);
        f->SetParNames("amp1", "mean", "sigma1", "amp2", "sigma2", "edge_lo", "edge_hi", "edge_sigma", "const");
        const Double_t top = h->GetBinContent(h->GetMaximumBin());
        f->SetParameters(top, h->GetMean(), 0.7 * h->GetStdDev(), 0.05 * top, 2.5 * h->GetStdDev(), lo, hi, 2.0, 0.1);
        f->SetParLimits(2, 1.0, 100.0); f->SetParLimits(3, 0.0, top); f->SetParLimits(4, 5.0, 300.0);
        f->SetParLimits(5, lo - 40., lo + 40.); f->SetParLimits(6, hi - 40., hi + 40.);
        f->SetParLimits(7, 0.1, 20.0); f->SetParLimits(8, 0.0, 0.01 * top);
        f->SetLineColor(kRed); f->SetNpx(1000);
        h->Fit(f, "0QRL");
        return f;
    };
    TF1* g_x = fit_profile(p_x, -code_x_half, code_x_half);
    TF1* g_y = fit_profile(p_y, yc - code_y_half, yc + code_y_half);
    printf("code rectangle: x %.1f %.1f, y %.1f %.1f (TPC frame)\n", -code_x_half, code_x_half, yc - code_y_half, yc + code_y_half);
    report("raw", f_x, f_y);
    report("hit", f_x_hit, f_y_hit);
    printf("prof x: edges %7.1f %7.1f (+-%.1f %.1f)  center %6.1f  width %6.1f  sigma %.1f | "
           "y: edges %7.1f %7.1f (+-%.1f %.1f)  center %6.1f  height %6.1f  sigma %.1f\n",
           g_x->GetParameter(5), g_x->GetParameter(6), g_x->GetParError(5), g_x->GetParError(6),
           0.5 * (g_x->GetParameter(5) + g_x->GetParameter(6)), g_x->GetParameter(6) - g_x->GetParameter(5), g_x->GetParameter(7),
           g_y->GetParameter(5), g_y->GetParameter(6), g_y->GetParError(5), g_y->GetParError(6),
           0.5 * (g_y->GetParameter(5) + g_y->GetParameter(6)), g_y->GetParameter(6) - g_y->GetParameter(5), g_y->GetParameter(7));

    auto draw_code_box = [&]() {
        TBox* b = new TBox(-code_x_half, yc - code_y_half, code_x_half, yc + code_y_half);
        b->SetFillStyle(0); b->SetLineColor(kRed); b->SetLineWidth(2); b->SetLineStyle(2); b->Draw();
    };
    auto draw_fit_box = [&](TF1* fx, TF1* fy) {
        TBox* b = new TBox(fx->GetParameter(1), fy->GetParameter(1), fx->GetParameter(2), fy->GetParameter(2));
        b->SetFillStyle(0); b->SetLineColor(kMagenta); b->SetLineWidth(2); b->Draw();
    };
    auto draw_edges = [](TF1* f, Double_t lo, Double_t hi) {
        for (Double_t e : {lo, hi}) { TLine* l = new TLine(e, 0, e, 1.1); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw(); }
        TLatex t; t.SetNDC(); t.SetTextSize(0.04);
        t.DrawLatex(0.15, 0.93 - 0.0, Form("fit edges %.1f / %.1f (#pm%.1f/%.1f), #sigma %.1f", f->GetParameter(1), f->GetParameter(2), f->GetParError(1), f->GetParError(2), std::fabs(f->GetParameter(3))));
    };

    TCanvas c("c", "", 1600, 1200);
    c.Print(Form("%s[", out));
    // page 1
    c.Clear(); c.Divide(2, 2);
    c.cd(1); h_chi->Draw("HIST");
    c.cd(2); gPad->SetLogy(); h_tdc->Draw("HIST");
    for (Double_t t : {tdc_lo, tdc_hi}) { TLine* l = new TLine(t, 0, t, h_tdc->GetMaximum()); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw(); }
    c.cd(3); gPad->SetLogy(); h_ntr->Draw("HIST");
    c.cd(4); gPad->SetLogy(); h_nfire->Draw("HIST");
    c.Print(out);
    // page 2
    c.Clear(); c.Divide(2, 2);
    c.cd(1); gPad->SetLogz(); h_all->Draw("colz"); draw_code_box();
    c.cd(2); gPad->SetLogz(); h_none->Draw("colz"); draw_code_box();
    c.cd(3); gPad->SetLogz(); h_fire->Draw("colz"); draw_code_box();
    c.cd(4); gPad->SetLogz(); h_none_hit->Draw("colz"); draw_code_box();
    c.Print(out);
    // page 3
    c.Clear(); c.Divide(3, 2);
    for (Int_t s = 0; s < n_seg_p0; s++) { c.cd(s + 1); h_seg[s]->Draw("colz"); draw_code_box(); }
    c.Print(out);
    // page 4
    c.Clear(); c.Divide(2, 2);
    c.cd(1); h_tr->SetMinimum(0); h_tr->SetMaximum(1); h_tr->Draw("colz"); draw_code_box(); draw_fit_box(f_x, f_y);
    { TLatex t; t.SetNDC(); t.SetTextSize(0.035); t.DrawLatex(0.12, 0.93, "red dashed: code rectangle, magenta: fit (raw)"); }
    c.cd(2); h_tr_hit->SetMinimum(0); h_tr_hit->SetMaximum(1); h_tr_hit->SetTitle("transmission (no plane-0 in htof_hit_seg) / all;x [mm];y [mm]");
    h_tr_hit->Draw("colz"); draw_code_box(); draw_fit_box(f_x_hit, f_y_hit);
    c.cd(3); t_x->SetMinimum(0); t_x->SetMaximum(1.1); t_x->Draw("E"); f_x->Draw("same"); draw_edges(f_x, -code_x_half, code_x_half);
    t_x_hit->SetLineColor(kGray + 1); t_x_hit->SetMarkerColor(kGray + 1); t_x_hit->Draw("E same");
    c.cd(4); t_y->SetMinimum(0); t_y->SetMaximum(1.1); t_y->Draw("E"); f_y->Draw("same"); draw_edges(f_y, yc - code_y_half, yc + code_y_half);
    t_y_hit->SetLineColor(kGray + 1); t_y_hit->SetMarkerColor(kGray + 1); t_y_hit->Draw("E same");
    c.Print(out);
    // page 5
    c.Clear(); c.Divide(2, 2);
    c.cd(1); gPad->SetLogz(); h_all->Draw("colz"); draw_code_box();
    { TBox* b = new TBox(g_x->GetParameter(5), g_y->GetParameter(5), g_x->GetParameter(6), g_y->GetParameter(6));
      b->SetFillStyle(0); b->SetLineColor(kMagenta); b->SetLineWidth(2); b->Draw(); }
    { TLatex t; t.SetNDC(); t.SetTextSize(0.035); t.DrawLatex(0.12, 0.93, "red dashed: code rectangle, magenta: profile-edge fit"); }
    c.cd(3); gPad->SetLogy(); p_x->SetMinimum(0.5); p_x->Draw("E"); g_x->Draw("same");
    for (Double_t e : {-code_x_half, code_x_half}) { TLine* l = new TLine(e, 0, e, p_x->GetMaximum()); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw(); }
    { TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.15, 0.93, Form("cut edges %.1f / %.1f (#pm%.1f/%.1f), #sigma %.1f", g_x->GetParameter(5), g_x->GetParameter(6), g_x->GetParError(5), g_x->GetParError(6), g_x->GetParameter(7))); }
    c.cd(4); gPad->SetLogy(); p_y->SetMinimum(0.5); p_y->Draw("E"); g_y->Draw("same");
    for (Double_t e : {yc - code_y_half, yc + code_y_half}) { TLine* l = new TLine(e, 0, e, p_y->GetMaximum()); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw(); }
    { TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.15, 0.93, Form("cut edges %.1f / %.1f (#pm%.1f/%.1f), #sigma %.1f", g_y->GetParameter(5), g_y->GetParameter(6), g_y->GetParError(5), g_y->GetParError(6), g_y->GetParameter(7))); }
    c.Print(out);
    c.Print(Form("%s]", out));
    printf("-> %s\n", out);
}
