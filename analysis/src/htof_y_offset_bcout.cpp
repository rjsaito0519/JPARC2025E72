// BcOut + Hodo: HTOF upstream y-offset / beam-window review -> PDF
// Usage: htof_y_offset_bcout <run_number> [-o <out.pdf>] [--pion | --kaon]
// Default BcOut/Hodo: run%05d_{BcOut,Hodo}.root under DATA_DIR or DECODE_DIR
// Default PDF: OUTPUT_DIR/img/runXXXXX/runXXXXX_HTOF_y_offset_bcout.pdf

#include "ana_helper.h"
#include "paths.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

namespace
{

constexpr Double_t kHtofZ = -499.7; // mm, FF frame
constexpr Int_t kUpstreamSegLo = 0;
constexpr Int_t kUpstreamSegHi = 5;

// Beam-window guide (draw only; matches TPCAnalyzer upstream hole)
constexpr Double_t kWinXHalf = 71.0;
constexpr Double_t kWinYLo = -50.0;
constexpr Double_t kWinYHi = 62.0;
constexpr Double_t kWinHeightNom = 112.0; // |design| height; soft reference only
constexpr Double_t kWinHeightSoftHalf = 20.0; // allow H in [112-20, 112+20]

constexpr Int_t kXYBins = 200;
constexpr Double_t kXYLo = -200.0;
constexpr Double_t kXYHi = 200.0;
constexpr Int_t kYBins = 200;
constexpr Double_t kYLo = -150.0;
constexpr Double_t kYHi = 150.0;

// Erf edge-fit windows (mm); centered near design window edges
constexpr Double_t kEdgeFitLoMin = -120.0;
constexpr Double_t kEdgeFitLoMax = 20.0;
constexpr Double_t kEdgeFitHiMin = 0.0;
constexpr Double_t kEdgeFitHiMax = 120.0;
constexpr Double_t kMinAllForRatio = 20.0; // skip bins with too few denom entries

constexpr Int_t kBeamFlagNoFilter = -1;
constexpr Int_t kBeamFlagPion = 1;
constexpr Int_t kBeamFlagKaon = 2;

struct RunStats
{
  Long64_t n_hodo = 0;
  Long64_t n_synced = 0;
  Long64_t n_filled_hit = 0;
  Long64_t n_filled_pass = 0;
  Long64_t n_filled_in_x = 0;
  Long64_t n_skipped_beam_flag = 0;
  Long64_t n_no_track = 0;
};

struct EdgeFit
{
  bool ok = false;
  Double_t y0 = 0.0;      // edge center (50% point)
  Double_t y0_err = 0.0;
  Double_t width = 0.0;   // transition width
  Double_t width_err = 0.0;
  Double_t amp = 0.0;
  Double_t offset = 0.0;
};

struct EdgeResult
{
  EdgeFit lo;
  EdgeFit hi;
  bool ok = false;          // usable center available
  bool soft_ok = false;     // joint plateau with soft H~112 succeeded
  Double_t y_center = 0.0;  // preferred (soft if available)
  Double_t y_center_err = 0.0;
  Double_t height = 0.0;    // preferred (soft if available)
  Double_t height_err = 0.0;
  Double_t y_center_free = 0.0; // independent (ylo+yhi)/2
  Double_t height_free = 0.0;
};

struct HistSet
{
  TH2D* xy_hit = nullptr;
  TH2D* xy_pass = nullptr;
  TH1D* y_all_in_x = nullptr;   // |x|<kWinXHalf, all
  TH1D* y_pass_in_x = nullptr;  // |x|<kWinXHalf, no upstream TDC hit
  TH1D* y_hit_in_x = nullptr;   // |x|<kWinXHalf, upstream TDC hit
  TH1D* y_hit = nullptr;
  TH1D* y_pass = nullptr;
  TH1D* r_pass = nullptr;       // pass/all vs y
};

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " <run_number> [-o <out.pdf>] [--pion | --kaon]\n"
            << "  Default BcOut/Hodo: run%05d_{BcOut,Hodo}.root under DATA_DIR or DECODE_DIR\n"
            << "  Upstream HTOF hit: raw TDC (seg " << kUpstreamSegLo << "-" << kUpstreamSegHi
            << ") inside UserParam HTOF_TDC window\n"
            << "  Default PDF: " << OUTPUT_DIR << "/img/runXXXXX/"
            << "runXXXXX_HTOF_y_offset_bcout.pdf\n";
}

void
setup_style()
{
  gROOT->SetBatch();
  gStyle->SetOptStat(110);
  gStyle->SetLabelSize(0.042, "XYZ");
  gStyle->SetTitleSize(0.042, "XYZ");
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.08);
}

Int_t
parse_run_number(const char* arg)
{
  if (!arg || !*arg) return -1;
  Int_t run = -1;
  if (std::sscanf(arg, "%d", &run) != 1 || run < 0) return -1;
  return run;
}

TString
beam_flag_filter_label(Int_t beam_flag_filter)
{
  if (beam_flag_filter == kBeamFlagPion) return "pion (beam_flag=1)";
  if (beam_flag_filter == kBeamFlagKaon) return "kaon (beam_flag=2)";
  return "none";
}

TString
default_pdf_path(Int_t run_num)
{
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
  return Form("%s/run%05d_HTOF_y_offset_bcout.pdf", img_dir.Data(), run_num);
}

TString
resolve_input_path(Int_t run_num, const char* tag)
{
  const TString candidates[] = {
    Form("%s/run%05d_%s.root", DATA_DIR.Data(), run_num, tag),
    Form("%s/run%05d/run%05d_%s.root", DECODE_DIR.Data(), run_num, run_num, tag),
  };
  for (const auto& path : candidates) {
    if (gSystem->AccessPathName(path) == 0) return path;
  }
  return candidates[0];
}

std::pair<Double_t, Double_t>
load_htof_tdc_window(Int_t run_num)
{
  auto range = ana_helper::get_user_param("HTOF_TDC", run_num);
  if ((range.first != 0.0 || range.second != 0.0) && range.first < range.second) {
    return range;
  }

  const TString candidates[] = {
    Form("%s/param/USER/e72/UserParam_run%05d", ANALYZER_DIR.Data(), run_num),
    Form("%s/param/USER/e72/UserParam_run%05d", WORK_DIR.Data(), run_num),
  };
  for (const auto& param_file : candidates) {
    struct stat st;
    if (stat(param_file.Data(), &st) != 0) continue;
    std::ifstream ifs(param_file.Data());
    if (!ifs.is_open()) continue;
    std::string line;
    while (std::getline(ifs, line)) {
      TString tline = line;
      if (!tline.BeginsWith("HTOF_TDC")) continue;
      TObjArray* tokens = tline.Tokenize(" \t");
      if (tokens->GetEntries() >= 3) {
        const Double_t tmin = ((TObjString*)tokens->At(1))->GetString().Atof();
        const Double_t tmax = ((TObjString*)tokens->At(2))->GetString().Atof();
        delete tokens;
        if (tmin < tmax) return {tmin, tmax};
      } else {
        delete tokens;
      }
    }
  }
  return {0.0, 0.0};
}

std::unordered_map<UInt_t, Long64_t>
build_event_index(TTree* tree)
{
  std::unordered_map<UInt_t, Long64_t> index;
  if (!tree) return index;

  UInt_t event_number = 0;
  tree->SetBranchAddress("event_number", &event_number);
  const Long64_t nent = tree->GetEntries();
  index.reserve(static_cast<std::size_t>(nent));
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tree->GetEntry(ie);
    index[event_number] = ie;
  }
  return index;
}

Double_t
blc_x_at_z(Double_t x0, Double_t u0, Double_t z)
{
  return x0 + u0 * z;
}

Double_t
blc_y_at_z(Double_t y0, Double_t v0, Double_t z)
{
  return y0 + v0 * z;
}

bool
tdc_in_window(const std::vector<Double_t>& tdcs, Double_t tmin, Double_t tmax)
{
  for (Double_t t : tdcs) {
    if (t >= tmin && t <= tmax) return true;
  }
  return false;
}

bool
has_upstream_htof_hit(const std::vector<Double_t>* raw_seg,
                      const std::vector<std::vector<Double_t>>* tdc_u,
                      const std::vector<std::vector<Double_t>>* tdc_d,
                      const std::vector<std::vector<Double_t>>* tdc_s,
                      Double_t tmin, Double_t tmax)
{
  if (!raw_seg) return false;
  for (std::size_t i = 0; i < raw_seg->size(); ++i) {
    const Int_t seg = static_cast<Int_t>((*raw_seg)[i]);
    if (seg < kUpstreamSegLo || seg > kUpstreamSegHi) continue;

    if (tdc_u && i < tdc_u->size() && tdc_in_window((*tdc_u)[i], tmin, tmax)) return true;
    if (tdc_d && i < tdc_d->size() && tdc_in_window((*tdc_d)[i], tmin, tmax)) return true;
    if (tdc_s && i < tdc_s->size() && tdc_in_window((*tdc_s)[i], tmin, tmax)) return true;
  }
  return false;
}

HistSet
make_hists()
{
  HistSet h;
  h.xy_hit = new TH2D(
    "xy_hit",
    Form("BcOut @ z=%.1f mm (upstream HTOF TDC hit);x [mm];y [mm]", kHtofZ),
    kXYBins, kXYLo, kXYHi, kXYBins, kXYLo, kXYHi);
  h.xy_pass = new TH2D(
    "xy_pass",
    Form("BcOut @ z=%.1f mm (no upstream HTOF TDC hit);x [mm];y [mm]", kHtofZ),
    kXYBins, kXYLo, kXYHi, kXYBins, kXYLo, kXYHi);
  h.y_all_in_x = new TH1D(
    "y_all_in_x",
    Form("all, |x|<%.0f @ z=%.1f;y [mm];entries", kWinXHalf, kHtofZ),
    kYBins, kYLo, kYHi);
  h.y_pass_in_x = new TH1D(
    "y_pass_in_x",
    Form("pass, |x|<%.0f @ z=%.1f;y [mm];entries", kWinXHalf, kHtofZ),
    kYBins, kYLo, kYHi);
  h.y_hit_in_x = new TH1D(
    "y_hit_in_x",
    Form("hit, |x|<%.0f @ z=%.1f;y [mm];entries", kWinXHalf, kHtofZ),
    kYBins, kYLo, kYHi);
  h.y_hit = new TH1D(
    "y_hit",
    Form("y @ z=%.1f (upstream HTOF TDC hit);y [mm];entries", kHtofZ),
    kYBins, kYLo, kYHi);
  h.y_pass = new TH1D(
    "y_pass",
    Form("y @ z=%.1f (no upstream HTOF TDC hit);y [mm];entries", kHtofZ),
    kYBins, kYLo, kYHi);
  h.r_pass = nullptr;
  return h;
}

void
prepare_hist_pad(TH2D* h)
{
  if (!h) return;
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.20);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.08);
  h->GetXaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleOffset(1.45);
  h->GetZaxis()->SetTitleOffset(1.25);
}

void
draw_placeholder(const char* msg)
{
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.5, msg);
}

void
draw_hist2_or_placeholder(TH2D* h, const char* opt = "COLZ")
{
  if (!h || h->GetEntries() <= 0) {
    draw_placeholder(h ? "empty histogram" : "histogram missing");
    return;
  }
  prepare_hist_pad(h);
  h->Draw(opt);
}

void
draw_hist1_or_placeholder(TH1D* h, const char* opt = "")
{
  if (!h || h->GetEntries() <= 0) {
    draw_placeholder(h ? "empty histogram" : "histogram missing");
    return;
  }
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.08);
  h->GetXaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleOffset(1.45);
  h->Draw(opt);
}

void
draw_window_box()
{
  auto* box = new TBox(-kWinXHalf, kWinYLo, kWinXHalf, kWinYHi);
  box->SetFillStyle(0);
  box->SetLineColor(kRed);
  box->SetLineWidth(2);
  box->SetLineStyle(2);
  box->Draw("same");
}

TH1D*
make_pass_ratio(TH1D* y_pass, TH1D* y_all)
{
  if (!y_pass || !y_all) return nullptr;
  auto* ratio = dynamic_cast<TH1D*>(y_pass->Clone("r_pass"));
  ratio->SetTitle(Form("pass/all (|x|<%.0f) @ z=%.1f;y [mm];N_{pass}/N_{all}",
                       kWinXHalf, kHtofZ));
  ratio->Divide(y_pass, y_all, 1.0, 1.0, "B");
  for (Int_t ib = 1; ib <= ratio->GetNbinsX(); ++ib) {
    if (y_all->GetBinContent(ib) < kMinAllForRatio) {
      ratio->SetBinContent(ib, 0.0);
      ratio->SetBinError(ib, 0.0);
    }
  }
  ratio->SetMinimum(0.0);
  ratio->SetMaximum(1.05);
  return ratio;
}

EdgeFit
fit_edge_erf(TH1D* ratio, bool lower_edge, Double_t xmin, Double_t xmax)
{
  EdgeFit fit;
  if (!ratio || ratio->GetEntries() <= 0) return fit;

  // lower: rise  amp*0.5*(1+Erf((x-y0)/w))+off
  // upper: fall  amp*0.5*(1+Erf((y0-x)/w))+off
  const char* formula = lower_edge
    ? "[0]*0.5*(1.+TMath::Erf((x-[1])/[2]))+[3]"
    : "[0]*0.5*(1.+TMath::Erf(([1]-x)/[2]))+[3]";
  auto* f = new TF1(Form("erf_%s_%s", ratio->GetName(), lower_edge ? "lo" : "hi"),
                    formula, xmin, xmax);
  f->SetParNames("amp", "y0", "width", "offset");
  f->SetParameters(0.8, lower_edge ? kWinYLo : kWinYHi, 5.0, 0.05);
  f->SetParLimits(0, 0.05, 1.5);
  f->SetParLimits(1, xmin, xmax);
  f->SetParLimits(2, 0.5, 40.0);
  f->SetParLimits(3, -0.1, 0.5);
  f->SetLineColor(lower_edge ? kBlue + 1 : kMagenta + 1);
  f->SetLineWidth(2);

  if (ratio->Fit(f, "RSQ0", "", xmin, xmax) != 0) return fit;

  fit.ok = true;
  fit.amp = f->GetParameter(0);
  fit.y0 = f->GetParameter(1);
  fit.width = f->GetParameter(2);
  fit.offset = f->GetParameter(3);
  fit.y0_err = f->GetParError(1);
  fit.width_err = f->GetParError(2);
  return fit;
}

EdgeResult
fit_window_edges(TH1D* ratio)
{
  EdgeResult res;
  res.lo = fit_edge_erf(ratio, true, kEdgeFitLoMin, kEdgeFitLoMax);
  res.hi = fit_edge_erf(ratio, false, kEdgeFitHiMin, kEdgeFitHiMax);

  if (res.lo.ok && res.hi.ok && res.hi.y0 > res.lo.y0) {
    res.y_center_free = 0.5 * (res.lo.y0 + res.hi.y0);
    res.height_free = res.hi.y0 - res.lo.y0;
    res.y_center = res.y_center_free;
    res.height = res.height_free;
    res.y_center_err = 0.5 * std::sqrt(res.lo.y0_err * res.lo.y0_err +
                                       res.hi.y0_err * res.hi.y0_err);
    res.ok = true;
  }

  if (!ratio || ratio->GetEntries() <= 0) return res;

  // Soft joint plateau: yc free, H softly near 112 mm (not fixed).
  // pass ~ amp * rise(y_lo) * fall(y_hi) + offset
  // y_lo = yc - H/2, y_hi = yc + H/2
  const Double_t xmin = kEdgeFitLoMin;
  const Double_t xmax = kEdgeFitHiMax;
  auto* f = new TF1(Form("erf_%s_soft", ratio->GetName()),
                    "[0]*0.5*(1.+TMath::Erf((x-([1]-[2]/2.))/[3]))"
                    "*0.5*(1.+TMath::Erf((([1]+[2]/2.)-x)/[4]))+[5]",
                    xmin, xmax);
  f->SetParNames("amp", "yc", "H", "w_lo", "w_hi", "offset");

  Double_t yc0 = res.ok ? res.y_center_free : 0.5 * (kWinYLo + kWinYHi);
  Double_t H0 = res.ok ? res.height_free : kWinHeightNom;
  if (H0 < kWinHeightNom - kWinHeightSoftHalf) H0 = kWinHeightNom;
  if (H0 > kWinHeightNom + kWinHeightSoftHalf) H0 = kWinHeightNom;
  const Double_t wlo0 = res.lo.ok ? res.lo.width : 5.0;
  const Double_t whi0 = res.hi.ok ? res.hi.width : 5.0;
  const Double_t amp0 = res.lo.ok ? res.lo.amp : 0.8;
  const Double_t off0 = res.lo.ok ? res.lo.offset : 0.05;

  f->SetParameters(amp0, yc0, H0, wlo0, whi0, off0);
  f->SetParLimits(0, 0.05, 1.5);
  f->SetParLimits(1, -80.0, 80.0);
  f->SetParLimits(2, kWinHeightNom - kWinHeightSoftHalf,
                  kWinHeightNom + kWinHeightSoftHalf);
  f->SetParLimits(3, 0.5, 40.0);
  f->SetParLimits(4, 0.5, 40.0);
  f->SetParLimits(5, -0.1, 0.5);
  f->SetLineColor(kOrange + 7);
  f->SetLineWidth(2);

  if (ratio->Fit(f, "RSQ0", "", xmin, xmax) == 0) {
    res.soft_ok = true;
    res.ok = true;
    res.y_center = f->GetParameter(1);
    res.y_center_err = f->GetParError(1);
    res.height = f->GetParameter(2);
    res.height_err = f->GetParError(2);
    // Update displayed edge markers from soft model
    res.lo.y0 = res.y_center - 0.5 * res.height;
    res.hi.y0 = res.y_center + 0.5 * res.height;
    res.lo.width = f->GetParameter(3);
    res.hi.width = f->GetParameter(4);
    res.lo.ok = true;
    res.hi.ok = true;
  }

  return res;
}

void
print_page(TCanvas* c, const TString& pdf_path)
{
  c->Print(pdf_path);
}

void
draw_title_page(TCanvas* c, const TString& pdf_path, Int_t run_num,
                const TString& path_bcout, const TString& path_hodo,
                Double_t tdc_min, Double_t tdc_max, const RunStats& stats,
                Int_t beam_flag_filter)
{
  c->Clear();
  c->cd();

  TTimeStamp now;
  now.Set();

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.92, "HTOF y offset from BcOut (+ Hodo raw TDC)");
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.5, 0.84, Form("run_number = %d", run_num));
  if (beam_flag_filter >= 0) {
    tex->SetTextColor(kBlue + 1);
    tex->DrawLatex(0.5, 0.78,
                   Form("beam_flag filter: %s", beam_flag_filter_label(beam_flag_filter).Data()));
    tex->SetTextColor(kBlack);
  }
  tex->SetTextAlign(12);
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.06, 0.70, Form("BcOut: %s", path_bcout.Data()));
  tex->DrawLatex(0.06, 0.64, Form("Hodo : %s", path_hodo.Data()));
  tex->SetTextAlign(22);
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.5, 0.54, Form("z_HTOF = %.1f mm (FF), upstream seg = %d-%d",
                                 kHtofZ, kUpstreamSegLo, kUpstreamSegHi));
  tex->DrawLatex(0.5, 0.48, Form("HTOF_TDC window: [%.0f, %.0f]", tdc_min, tdc_max));
  tex->DrawLatex(0.5, 0.42, Form("Hodo entries: %lld, synced+track: %lld",
                                 stats.n_hodo, stats.n_synced));
  tex->DrawLatex(0.5, 0.36, Form("fills hit / pass: %lld / %lld (|x|<%.0f for ratio: %lld)",
                                 stats.n_filled_hit, stats.n_filled_pass,
                                 kWinXHalf, stats.n_filled_in_x));
  if (beam_flag_filter >= 0) {
    tex->DrawLatex(0.5, 0.30,
                   Form("Skipped by beam_flag: %lld", stats.n_skipped_beam_flag));
  }
  tex->DrawLatex(0.5, 0.22,
                 Form("y-offset: pass/all edges; soft H~%.0f#pm%.0f mm for center",
                      kWinHeightNom, kWinHeightSoftHalf));
  tex->DrawLatex(0.5, 0.12, Form("created: %s", now.AsString("s")));
  print_page(c, pdf_path);
}

void
draw_xy_page(TCanvas* c, const TString& pdf_path, TH2D* h, const char* note)
{
  c->Clear();
  c->cd();
  draw_hist2_or_placeholder(h);
  if (h && h->GetEntries() > 0) draw_window_box();

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.028);
  tex->SetTextColor(kRed + 1);
  tex->DrawLatex(0.15, 0.92, note);
  print_page(c, pdf_path);
}

void
draw_edge_ratio_page(TCanvas* c, const TString& pdf_path, TH1D* ratio,
                     const EdgeResult& edges)
{
  c->Clear();
  c->cd();
  if (!ratio || ratio->GetEntries() <= 0) {
    draw_placeholder("empty pass/all ratio");
    print_page(c, pdf_path);
    return;
  }

  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.16);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.5);
  ratio->SetLineColor(kBlack);
  ratio->Draw("E1");

  auto* fsoft = ratio->GetFunction(Form("erf_%s_soft", ratio->GetName()));
  if (fsoft) fsoft->Draw("same");

  if (edges.lo.ok) {
    auto* flo = ratio->GetFunction(Form("erf_%s_lo", ratio->GetName()));
    if (flo && !edges.soft_ok) flo->Draw("same");
    auto* llo = new TLine(edges.lo.y0, 0.0, edges.lo.y0, 1.05);
    llo->SetLineColor(kBlue + 1);
    llo->SetLineStyle(2);
    llo->Draw("same");
  }
  if (edges.hi.ok) {
    auto* fhi = ratio->GetFunction(Form("erf_%s_hi", ratio->GetName()));
    if (fhi && !edges.soft_ok) fhi->Draw("same");
    auto* lhi = new TLine(edges.hi.y0, 0.0, edges.hi.y0, 1.05);
    lhi->SetLineColor(kMagenta + 1);
    lhi->SetLineStyle(2);
    lhi->Draw("same");
  }
  if (edges.ok) {
    auto* lc = new TLine(edges.y_center, 0.0, edges.y_center, 1.05);
    lc->SetLineColor(kRed);
    lc->SetLineWidth(2);
    lc->Draw("same");
  }

  // design edges
  auto* dlo = new TLine(kWinYLo, 0.0, kWinYLo, 1.05);
  dlo->SetLineColor(kGray + 2);
  dlo->SetLineStyle(3);
  dlo->Draw("same");
  auto* dhi = new TLine(kWinYHi, 0.0, kWinYHi, 1.05);
  dhi->SetLineColor(kGray + 2);
  dhi->SetLineStyle(3);
  dhi->Draw("same");

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.026);
  tex->DrawLatex(0.15, 0.955,
                 Form("pass/all; soft joint H#approx%.0f#pm%.0f mm (orange); gray=design",
                      kWinHeightNom, kWinHeightSoftHalf));
  if (edges.lo.ok) {
    tex->DrawLatex(0.15, 0.915,
                   Form("y_lo = %.2f mm", edges.lo.y0));
  }
  if (edges.hi.ok) {
    tex->DrawLatex(0.45, 0.915,
                   Form("y_hi = %.2f mm", edges.hi.y0));
  }
  if (edges.ok) {
    tex->SetTextColor(kRed + 1);
    tex->DrawLatex(0.15, 0.875,
                   Form("y_center = %.2f #pm %.2f mm  (H=%.1f mm%s)",
                        edges.y_center, edges.y_center_err, edges.height,
                        edges.soft_ok ? ", soft" : ", free"));
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.15, 0.835,
                   Form("vs y=0 => offset ~ %.2f mm; free (ylo+yhi)/2=%.2f (H_free=%.1f)",
                        edges.y_center,
                        edges.y_center_free, edges.height_free));
  } else {
    tex->DrawLatex(0.15, 0.875, "edge fit incomplete");
  }
  print_page(c, pdf_path);
}

void
draw_counts_page(TCanvas* c, const TString& pdf_path, TH1D* y_all, TH1D* y_pass,
                 TH1D* y_hit)
{
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  if (y_all && y_all->GetEntries() > 0) {
    y_all->SetLineColor(kBlack);
    y_all->SetLineWidth(2);
    y_all->Draw();
    if (y_pass && y_pass->GetEntries() > 0) {
      y_pass->SetLineColor(kBlue + 1);
      y_pass->SetLineWidth(2);
      y_pass->Draw("same");
    }
    if (y_hit && y_hit->GetEntries() > 0) {
      y_hit->SetLineColor(kRed + 1);
      y_hit->SetLineWidth(2);
      y_hit->Draw("same");
    }
    auto* leg = new TLegend(0.55, 0.70, 0.90, 0.88);
    if (y_all) leg->AddEntry(y_all, "all", "l");
    if (y_pass && y_pass->GetEntries() > 0) leg->AddEntry(y_pass, "pass", "l");
    if (y_hit && y_hit->GetEntries() > 0) leg->AddEntry(y_hit, "hit", "l");
    leg->Draw();
  } else {
    draw_placeholder("empty |x| counts");
  }
  c->cd(2);
  draw_hist1_or_placeholder(y_pass);
  print_page(c, pdf_path);
}

void
draw_compare_page(TCanvas* c, const TString& pdf_path, TH1D* y_hit, TH1D* y_pass,
                  const EdgeResult& edges)
{
  c->Clear();
  c->cd();

  if ((!y_hit || y_hit->GetEntries() <= 0) && (!y_pass || y_pass->GetEntries() <= 0)) {
    draw_placeholder("empty histograms");
    print_page(c, pdf_path);
    return;
  }

  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.12);

  TH1D* first = nullptr;
  if (y_hit && y_hit->GetEntries() > 0) {
    y_hit->SetLineColor(kRed + 1);
    y_hit->SetLineWidth(2);
    y_hit->SetTitle(Form("y @ z=%.1f: hit vs pass;y [mm];entries", kHtofZ));
    y_hit->Draw();
    first = y_hit;
  }
  if (y_pass && y_pass->GetEntries() > 0) {
    y_pass->SetLineColor(kBlue + 1);
    y_pass->SetLineWidth(2);
    if (first) y_pass->Draw("same");
    else {
      y_pass->SetTitle(Form("y @ z=%.1f: hit vs pass;y [mm];entries", kHtofZ));
      y_pass->Draw();
    }
  }

  auto* leg = new TLegend(0.62, 0.72, 0.92, 0.88);
  if (y_hit && y_hit->GetEntries() > 0) leg->AddEntry(y_hit, "TDC hit (seg0-5)", "l");
  if (y_pass && y_pass->GetEntries() > 0) leg->AddEntry(y_pass, "no TDC hit", "l");
  leg->Draw();

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.14, 0.94,
                 "Hit: hole in xy. Soft H~112 joint edges -> window y center.");
  if (edges.ok) {
    tex->DrawLatex(0.14, 0.90,
                   Form("y_center = %.2f #pm %.2f mm (offset vs y=0)",
                        edges.y_center, edges.y_center_err));
  }
  print_page(c, pdf_path);
}

Int_t
analyze(Int_t run_num, const TString& output_path_in, Int_t beam_flag_filter)
{
  setup_style();

  const auto tdc_win = load_htof_tdc_window(run_num);
  if (!(tdc_win.first < tdc_win.second)) {
    std::cerr << "Error: could not load HTOF_TDC window for run " << run_num
              << " (UserParam_runXXXXX)." << std::endl;
    return 1;
  }
  const Double_t tdc_min = tdc_win.first;
  const Double_t tdc_max = tdc_win.second;

  const TString path_bcout = resolve_input_path(run_num, "BcOut");
  const TString path_hodo = resolve_input_path(run_num, "Hodo");

  auto* f_bcout = TFile::Open(path_bcout);
  auto* f_hodo = TFile::Open(path_hodo);
  if (!f_bcout || f_bcout->IsZombie()) {
    std::cerr << "Error: could not open " << path_bcout << std::endl;
    return 1;
  }
  if (!f_hodo || f_hodo->IsZombie()) {
    std::cerr << "Error: could not open " << path_hodo << std::endl;
    f_bcout->Close();
    return 1;
  }

  auto* tree_bcout = dynamic_cast<TTree*>(f_bcout->Get("bcout"));
  auto* tree_hodo = dynamic_cast<TTree*>(f_hodo->Get("hodo"));
  if (!tree_bcout || !tree_hodo) {
    std::cerr << "Error: required trees not found (bcout, hodo)." << std::endl;
    f_bcout->Close();
    f_hodo->Close();
    return 1;
  }

  std::cout << "BcOut: " << path_bcout << std::endl;
  std::cout << "Hodo : " << path_hodo << std::endl;
  std::cout << "HTOF_TDC: [" << tdc_min << ", " << tdc_max << "]" << std::endl;
  std::cout << "z_HTOF: " << kHtofZ << " mm" << std::endl;
  if (beam_flag_filter >= 0) {
    std::cout << "beam_flag filter: " << beam_flag_filter_label(beam_flag_filter)
              << std::endl;
  }

  const auto bcout_idx = build_event_index(tree_bcout);
  HistSet h = make_hists();
  RunStats stats;
  stats.n_hodo = tree_hodo->GetEntries();

  UInt_t event_number = 0;
  Int_t beam_flag = 0;
  std::vector<Double_t>* htof_raw_seg = nullptr;
  std::vector<std::vector<Double_t>>* htof_tdc_u = nullptr;
  std::vector<std::vector<Double_t>>* htof_tdc_d = nullptr;
  std::vector<std::vector<Double_t>>* htof_tdc_s = nullptr;
  tree_hodo->SetBranchAddress("event_number", &event_number);
  tree_hodo->SetBranchAddress("beam_flag", &beam_flag);
  tree_hodo->SetBranchAddress("htof_raw_seg", &htof_raw_seg);
  tree_hodo->SetBranchAddress("htof_tdc_u", &htof_tdc_u);
  tree_hodo->SetBranchAddress("htof_tdc_d", &htof_tdc_d);
  tree_hodo->SetBranchAddress("htof_tdc_s", &htof_tdc_s);

  Int_t ntrack = 0;
  std::vector<Double_t>* x0 = nullptr;
  std::vector<Double_t>* y0 = nullptr;
  std::vector<Double_t>* u0 = nullptr;
  std::vector<Double_t>* v0 = nullptr;
  tree_bcout->SetBranchAddress("ntrack", &ntrack);
  tree_bcout->SetBranchAddress("x0", &x0);
  tree_bcout->SetBranchAddress("y0", &y0);
  tree_bcout->SetBranchAddress("u0", &u0);
  tree_bcout->SetBranchAddress("v0", &v0);

  for (Long64_t ie = 0; ie < stats.n_hodo; ++ie) {
    tree_hodo->GetEntry(ie);
    if (beam_flag_filter >= 0 && beam_flag != beam_flag_filter) {
      ++stats.n_skipped_beam_flag;
      continue;
    }

    const auto it = bcout_idx.find(event_number);
    if (it == bcout_idx.end()) continue;
    tree_bcout->GetEntry(it->second);

    if (ntrack < 1 || !x0 || !y0 || !u0 || !v0 ||
        x0->empty() || y0->empty() || u0->empty() || v0->empty()) {
      ++stats.n_no_track;
      continue;
    }
    ++stats.n_synced;

    const Double_t x = blc_x_at_z((*x0)[0], (*u0)[0], kHtofZ);
    const Double_t y = blc_y_at_z((*y0)[0], (*v0)[0], kHtofZ);
    const bool hit = has_upstream_htof_hit(htof_raw_seg, htof_tdc_u, htof_tdc_d, htof_tdc_s,
                                           tdc_min, tdc_max);

    if (hit) {
      h.xy_hit->Fill(x, y);
      h.y_hit->Fill(y);
      ++stats.n_filled_hit;
    } else {
      h.xy_pass->Fill(x, y);
      h.y_pass->Fill(y);
      ++stats.n_filled_pass;
    }

    if (std::fabs(x) < kWinXHalf) {
      h.y_all_in_x->Fill(y);
      ++stats.n_filled_in_x;
      if (hit) h.y_hit_in_x->Fill(y);
      else h.y_pass_in_x->Fill(y);
    }
  }

  TString output_path = output_path_in;
  if (output_path.IsNull()) {
    output_path = default_pdf_path(run_num);
    std::cout << "Default output: " << output_path << std::endl;
  }

  h.r_pass = make_pass_ratio(h.y_pass_in_x, h.y_all_in_x);
  const EdgeResult edges = fit_window_edges(h.r_pass);
  if (edges.lo.ok) {
    std::cout << "edge lo: y0=" << edges.lo.y0 << " +/- " << edges.lo.y0_err
              << " mm, w=" << edges.lo.width << std::endl;
  }
  if (edges.hi.ok) {
    std::cout << "edge hi: y0=" << edges.hi.y0 << " +/- " << edges.hi.y0_err
              << " mm, w=" << edges.hi.width << std::endl;
  }
  if (edges.ok) {
    std::cout << "window y_center=" << edges.y_center << " +/- " << edges.y_center_err
              << " mm, height=" << edges.height
              << (edges.soft_ok ? " mm (soft H~112)" : " mm (free)")
              << ", free_center=" << edges.y_center_free
              << " free_H=" << edges.height_free
              << ", offset_vs_y0=" << edges.y_center << " mm" << std::endl;
  }

  auto* canvas = new TCanvas("htof_y_offset_bcout", output_path, 1200, 900);
  canvas->Print(output_path + "[");

  draw_title_page(canvas, output_path, run_num, path_bcout, path_hodo,
                  tdc_min, tdc_max, stats, beam_flag_filter);
  draw_xy_page(canvas, output_path, h.xy_hit,
               "red dashed: beam-window guide; expect hole if TDC hits on seg0-5");
  draw_xy_page(canvas, output_path, h.xy_pass,
               "no upstream TDC hit: beam through window");
  draw_counts_page(canvas, output_path, h.y_all_in_x, h.y_pass_in_x, h.y_hit_in_x);
  draw_edge_ratio_page(canvas, output_path, h.r_pass, edges);
  draw_compare_page(canvas, output_path, h.y_hit, h.y_pass, edges);

  canvas->Print(output_path + "]");
  delete canvas;

  f_bcout->Close();
  f_hodo->Close();

  std::cout << "Wrote " << output_path
            << " (run=" << run_num
            << ", hit=" << stats.n_filled_hit
            << ", pass=" << stats.n_filled_pass
            << ", in_x=" << stats.n_filled_in_x
            << ")" << std::endl;
  return 0;
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  Int_t run_num = -1;
  TString output_path;
  bool has_output = false;
  Int_t beam_flag_filter = kBeamFlagNoFilter;

  for (Int_t i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    if (std::strcmp(arg, "-o") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      output_path = argv[++i];
      has_output = true;
      continue;
    }
    if (std::strcmp(arg, "--pion") == 0) {
      if (beam_flag_filter != kBeamFlagNoFilter) {
        std::cerr << "Specify only one of --pion or --kaon." << std::endl;
        usage(argv[0]);
        return 1;
      }
      beam_flag_filter = kBeamFlagPion;
      continue;
    }
    if (std::strcmp(arg, "--kaon") == 0) {
      if (beam_flag_filter != kBeamFlagNoFilter) {
        std::cerr << "Specify only one of --pion or --kaon." << std::endl;
        usage(argv[0]);
        return 1;
      }
      beam_flag_filter = kBeamFlagKaon;
      continue;
    }
    if (arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
    if (run_num >= 0) {
      std::cerr << "Multiple run numbers specified." << std::endl;
      usage(argv[0]);
      return 1;
    }
    run_num = parse_run_number(arg);
    if (run_num < 0) {
      std::cerr << "Invalid run number: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
  }

  if (run_num < 0) {
    usage(argv[0]);
    return 1;
  }

  if (has_output) gSystem->ExpandPathName(output_path);
  return analyze(run_num, has_output ? output_path : TString(), beam_flag_filter);
}
