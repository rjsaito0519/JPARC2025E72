// D5 wire residual review (UserD5Tracking histograms) -> multi-page PDF
// Usage: d5_wire_resi_review [-o <out.pdf>] [--dcgeo <DCGEO>] <decode.root>
// Default PDF: OUTPUT_DIR/img/runXXXXX/runXXXXX_D5_wire_resi_review.pdf
// Default DCGEO: ANALYZER_DIR/param/DCGEO/DCGeomParam_e72_example
//
// Extra (D5Fit only): project plane residuals to (rx,ry) and least-squares
// rigid (dx,dy) with cos*dx+sin*dy ≈ -mu_s (same sign as update_residual Ofs).
// TA per layer from DCGEO (fallback: U=-45, V=+45).

#include "ana_helper.h"
#include "config.h"
#include "paths.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

Config& conf = Config::getInstance();

const char* kDetectors[] = {"BLC1a", "BLC1b", "BLC2a", "BLC2b"};
const Int_t kNumDetectors = 4;

const char* kPlaneLabels[] = {
  "U1", "UP1", "V1", "VP1", "U2", "UP2", "V2", "VP2"
};

constexpr Double_t kFitNSigma = 1.3; // same as ana_helper::residual_fit
constexpr Int_t kProjBins = 200;
constexpr Double_t kProjLo = -5.0;   // per a/b chamber
constexpr Double_t kProjHi = 5.0;
constexpr Double_t kGroupProjLo = -10.0; // BLC1/BLC2 a+b combined (wider)
constexpr Double_t kGroupProjHi = 10.0;

struct GaussFit {
  bool ok = false;
  Double_t mean = 0.0;
  Double_t mean_err = 0.0;
  Double_t sigma = 0.0;
};

struct DetectorShift {
  bool ok = false;
  Double_t dx = 0.0;
  Double_t dy = 0.0;
  Double_t mu_s[8] = {};
  bool mu_ok[8] = {};
};

// BLC1 = a+b (16 planes), BLC2 = a+b (16 planes)
struct GroupShift {
  bool ok = false;
  Double_t dx = 0.0;
  Double_t dy = 0.0;
  const char* group = nullptr;   // "BLC1" or "BLC2"
  const char* flag = nullptr;    // "--bcin" or "--bcout"
  const char* det_a = nullptr;
  const char* det_b = nullptr;
  Double_t mu_s[16] = {};
  bool mu_ok[16] = {};
};

// layer name (e.g. BLC1a-U1) -> TA [deg]
std::map<std::string, Double_t> g_ta_deg;
TString g_dcgeo_path;

TString
default_dcgeo_path()
{
  return Form("%s/param/DCGEO/DCGeomParam_e72_example", ANALYZER_DIR.Data());
}

bool
load_dcgeo_ta(const TString& path)
{
  g_ta_deg.clear();
  g_dcgeo_path = path;
  std::ifstream ifs(path.Data());
  if (!ifs) {
    std::cerr << "Cannot open DCGEO: " << path << std::endl;
    return false;
  }
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#' || line[0] == '/') continue;
    // strip trailing comments
    const auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    std::istringstream iss(line);
    std::string id, name;
    Double_t x = 0, y = 0, z = 0, ta = 0;
    if (!(iss >> id >> name >> x >> y >> z >> ta)) continue;
    if (name.find("BLC") != 0) continue;
    if (name.find('-') == std::string::npos) continue;
    g_ta_deg[name] = ta;
  }
  std::cout << "DCGEO TA: " << path << " (" << g_ta_deg.size()
            << " BLC layers)" << std::endl;
  return !g_ta_deg.empty();
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0
            << " [-o <out.pdf>] [--dcgeo <DCGEO>] <decode.root>\n"
            << "  Reads D5WireResi_* histograms from D5 Tracking decode ROOT.\n"
            << "  Also estimates rigid (dx,dy) from D5Fit plane residuals.\n"
            << "  TA from DCGEO (default: " << default_dcgeo_path() << ").\n"
            << "  Default PDF: " << OUTPUT_DIR << "/img/runXXXXX/"
            << "runXXXXX_D5_wire_resi_review.pdf\n";
}

TObject*
get_object(TFile* f, const char* name)
{
  if (!f || !name) return nullptr;
  return f->Get(name);
}

Long64_t
hist_entries(TFile* f, const char* name)
{
  auto* h = dynamic_cast<TH1*>(get_object(f, name));
  if (!h) return 0;
  return static_cast<Long64_t>(h->GetEntries());
}

void
setup_style()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
}

Int_t
parse_run_from_path(const TString& path)
{
  const std::string s(path.Data());
  const std::size_t pos = s.find("run");
  if (pos == std::string::npos || pos + 8 > s.size()) return -1;
  Int_t run = -1;
  if (std::sscanf(s.c_str() + pos + 3, "%5d", &run) == 1) return run;
  return -1;
}

TString
default_pdf_path(Int_t run_num)
{
  const Int_t run = run_num >= 0 ? run_num : 0;
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run);
  return Form("%s/run%05d_D5_wire_resi_review.pdf", img_dir.Data(), run);
}

void
print_page(TCanvas* c, const TString& pdf_path)
{
  c->Print(pdf_path);
}

void
draw_placeholder(const char* message)
{
  auto* tex = new TLatex(0.5, 0.5, message);
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->Draw();
}

void
draw_plane_label(const char* label)
{
  auto* tex = new TLatex(0.18, 0.88, label);
  tex->SetNDC();
  tex->SetTextAlign(11);
  tex->SetTextSize(0.045);
  tex->Draw();
}

void
draw_gauss_mean_vline(TH1* h, Double_t mean)
{
  if (!h || !std::isfinite(mean)) return;
  const Double_t ymax = h->GetMaximum();
  if (!(ymax > 0.0)) return;
  auto* line = new TLine(mean, 0.0, mean, ymax);
  line->SetLineStyle(2); // dashed
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("SAME");
}

// Prefer DCGEO TA for detector-plane (e.g. BLC1a-U1); else nominal ±45.
Double_t
plane_ta_deg(const char* detector, Int_t plane)
{
  if (plane < 0 || plane >= 8 || !detector) return 0.0;
  const char* lab = kPlaneLabels[plane];
  const std::string key = std::string(detector) + "-" + lab;
  const auto it = g_ta_deg.find(key);
  if (it != g_ta_deg.end()) return it->second;
  if (lab[0] == 'U') return -45.0;
  if (lab[0] == 'V') return 45.0;
  return 0.0;
}

GaussFit
fit_gauss_1p3(TH1* h)
{
  GaussFit fit;
  if (!h || h->GetEntries() < 20) return fit;

  const Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
  const Double_t width = h->GetStdDev();
  if (!(width > 0.0) || !std::isfinite(width)) return fit;

  const TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";
  const Double_t axis_min = h->GetXaxis()->GetXmin();
  const Double_t axis_max = h->GetXaxis()->GetXmax();

  Double_t lower = std::max(peak_pos - kFitNSigma * width, axis_min);
  Double_t upper = std::min(peak_pos + kFitNSigma * width, axis_max);
  if (!(upper > lower)) return fit;

  auto* f_pre = new TF1(Form("pre_gaus_%s", h->GetName()), "gausn", lower, upper);
  f_pre->SetParameter(1, peak_pos);
  f_pre->SetParameter(2, width * 0.9);
  h->Fit(f_pre, fit_option.Data(), "", lower, upper);
  const Double_t pre_mean = f_pre->GetParameter(1);
  const Double_t pre_sigma = std::fabs(f_pre->GetParameter(2));
  const Double_t pre_amp = f_pre->GetParameter(0);
  delete f_pre;
  if (!(pre_sigma > 0.0) || !std::isfinite(pre_mean) || !std::isfinite(pre_sigma))
    return fit;

  lower = std::max(pre_mean - kFitNSigma * pre_sigma, axis_min);
  upper = std::min(pre_mean + kFitNSigma * pre_sigma, axis_max);
  if (!(upper > lower)) return fit;

  auto* gaus = new TF1(Form("gaus_%s", h->GetName()), "gausn", lower, upper);
  gaus->SetLineColor(kOrange + 1);
  gaus->SetLineWidth(2);
  gaus->SetNpx(500);
  gaus->SetParameter(0, pre_amp);
  gaus->SetParameter(1, pre_mean);
  gaus->SetParameter(2, pre_sigma * 0.9);
  h->Fit(gaus, fit_option.Data(), "", lower, upper);

  const Double_t mean = gaus->GetParameter(1);
  const Double_t sigma = std::fabs(gaus->GetParameter(2));
  if (!(sigma > 0.0) || !std::isfinite(mean) || !std::isfinite(sigma))
    return fit;

  fit.ok = true;
  fit.mean = mean;
  fit.sigma = sigma;
  fit.mean_err = gaus->GetParError(1);
  gaus->SetRange(lower, upper);
  return fit;
}

Double_t
plane_mu_s(TH1* h, GaussFit* fit_out = nullptr)
{
  const GaussFit fit = fit_gauss_1p3(h);
  if (fit_out) *fit_out = fit;
  if (fit.ok) return fit.mean;
  if (h && h->GetEntries() > 0) return h->GetMean();
  return 0.0;
}

bool
solve_dx_dy(const Double_t* mu, const bool* ok, const Double_t* cta,
            const Double_t* sta, Int_t nplane, Double_t& dx, Double_t& dy)
{
  // Match update_residual.py: Ofs correction uses -mu_s
  //   cos(TA)*dx + sin(TA)*dy ≈ -mu_s
  Double_t a11 = 0.0, a12 = 0.0, a22 = 0.0, b1 = 0.0, b2 = 0.0;
  Int_t n = 0;
  for (Int_t i = 0; i < nplane; ++i) {
    if (!ok[i]) continue;
    const Double_t rhs = -mu[i];
    a11 += cta[i] * cta[i];
    a12 += cta[i] * sta[i];
    a22 += sta[i] * sta[i];
    b1 += cta[i] * rhs;
    b2 += sta[i] * rhs;
    ++n;
  }
  if (n < 2) return false;
  const Double_t det = a11 * a22 - a12 * a12;
  if (!(std::fabs(det) > 1.0e-12)) return false;
  dx = (a22 * b1 - a12 * b2) / det;
  dy = (-a12 * b1 + a11 * b2) / det;
  return std::isfinite(dx) && std::isfinite(dy);
}

void
fill_rx_ry_from_planes(TFile* f, const char* detector, TH1D* h_rx, TH1D* h_ry)
{
  const Int_t nplane = conf.num_of_ch.at("blc");
  for (Int_t i = 0; i < nplane && i < 8; ++i) {
    auto* h = dynamic_cast<TH1*>(
      get_object(f, Form("D5WireResi_%s_D5Fit_plane%d", detector, i)));
    if (!h || h->GetEntries() <= 0) continue;
    const Double_t ta = plane_ta_deg(detector, i) * TMath::DegToRad();
    const Double_t c = TMath::Cos(ta);
    const Double_t s = TMath::Sin(ta);
    const Int_t nb = h->GetNbinsX();
    for (Int_t ib = 1; ib <= nb; ++ib) {
      const Double_t w = h->GetBinContent(ib);
      if (!(w > 0.0)) continue;
      const Double_t r = h->GetBinCenter(ib);
      h_rx->Fill(r * c, w);
      h_ry->Fill(r * s, w);
    }
  }
}

DetectorShift
estimate_detector_shift(TFile* f, const char* detector)
{
  DetectorShift out;
  const Int_t nplane = conf.num_of_ch.at("blc");
  Double_t cta[8] = {};
  Double_t sta[8] = {};
  for (Int_t i = 0; i < nplane && i < 8; ++i) {
    auto* h = dynamic_cast<TH1*>(
      get_object(f, Form("D5WireResi_%s_D5Fit_plane%d", detector, i)));
    const Double_t ta = plane_ta_deg(detector, i) * TMath::DegToRad();
    cta[i] = TMath::Cos(ta);
    sta[i] = TMath::Sin(ta);
    if (!h || h->GetEntries() < 20) {
      out.mu_ok[i] = false;
      out.mu_s[i] = 0.0;
      continue;
    }
    GaussFit gf;
    out.mu_s[i] = plane_mu_s(h, &gf);
    out.mu_ok[i] = true;
  }
  out.ok = solve_dx_dy(out.mu_s, out.mu_ok, cta, sta, 8, out.dx, out.dy);
  return out;
}

GroupShift
estimate_group_shift(TFile* f, const char* group, const char* flag,
                     const char* det_a, const char* det_b)
{
  GroupShift out;
  out.group = group;
  out.flag = flag;
  out.det_a = det_a;
  out.det_b = det_b;
  Double_t cta[16] = {};
  Double_t sta[16] = {};
  const char* dets[2] = {det_a, det_b};
  for (Int_t id = 0; id < 2; ++id) {
    for (Int_t i = 0; i < 8; ++i) {
      const Int_t k = id * 8 + i;
      auto* h = dynamic_cast<TH1*>(
        get_object(f, Form("D5WireResi_%s_D5Fit_plane%d", dets[id], i)));
      const Double_t ta = plane_ta_deg(dets[id], i) * TMath::DegToRad();
      cta[k] = TMath::Cos(ta);
      sta[k] = TMath::Sin(ta);
      if (!h || h->GetEntries() < 20) {
        out.mu_ok[k] = false;
        out.mu_s[k] = 0.0;
        continue;
      }
      GaussFit gf;
      out.mu_s[k] = plane_mu_s(h, &gf);
      out.mu_ok[k] = true;
    }
  }
  out.ok = solve_dx_dy(out.mu_s, out.mu_ok, cta, sta, 16, out.dx, out.dy);
  return out;
}

void
draw_title_page(TCanvas* c, const TString& pdf_path, TFile* f, Int_t run_num)
{
  c->Clear();
  c->cd();
  TTimeStamp now;
  now.Set();

  Long64_t total = 0;
  for (Int_t id = 0; id < kNumDetectors; ++id) {
    total += hist_entries(f, Form("D5WireResi_%s_LocalFit_vs_plane", kDetectors[id]));
  }

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.84, "D5 wire residual review (rank0 D5 pairs)");
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.74, f->GetName());
  tex->DrawLatex(0.5, 0.68, Form("run_number = %d", run_num));
  tex->DrawLatex(0.5, 0.62,
                 Form("D5WireResi LocalFit_vs_plane entries (sum) = %lld",
                      static_cast<long long>(total)));
  tex->DrawLatex(0.5, 0.54, "LocalFit = chamber DoFit | D5Fit = joint-fit prediction");
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.5, 0.46, "Policy: BLC1a / BLC1b independent Ofs; BLC2 fixed; TA fixed");
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.40, "Use --residual-offset per BLC1a/BLC1b (not --bcin). Do not apply BLC2.");
  tex->DrawLatex(0.5, 0.34, Form("TA from: %s (read-only for LS; not tuned here)", g_dcgeo_path.Data()));
  tex->DrawLatex(0.5, 0.28, "also: D5Fit r vs s_{perp}/s_{para} (monitor; #DeltaTA not applied)");
  tex->DrawLatex(0.5, 0.22, "end pages: D5Track #Delta / p / #chi^{2} / residuals (QA)");
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.14, Form("created: %s", now.AsString("s")));
  print_page(c, pdf_path);
}

void
draw_vs_plane_page(TCanvas* c, const TString& pdf_path, TFile* f,
                   const char* detector)
{
  c->Clear();
  c->Divide(2, 1);
  const char* names[] = {"LocalFit", "D5Fit"};
  for (Int_t i = 0; i < 2; ++i) {
    c->cd(i + 1);
    const TString hname =
      Form("D5WireResi_%s_%s_vs_plane", detector, names[i]);
    auto* h = get_object(f, hname);
    if (!h) {
      draw_placeholder(Form("NOT FOUND:\n%s", hname.Data()));
      continue;
    }
    gPad->SetLogz();
    h->Draw("colz");
    auto* tex = new TLatex(0.18, 0.88, Form("%s %s", detector, names[i]));
    tex->SetNDC();
    tex->SetTextAlign(11);
    tex->SetTextSize(0.045);
    tex->Draw();
  }
  print_page(c, pdf_path);
}

void
draw_plane_overlay_page(TCanvas* c, const TString& pdf_path, TFile* f,
                        const char* detector)
{
  const Int_t nplane = conf.num_of_ch.at("blc");
  c->Clear();
  c->Divide(4, 2);
  for (Int_t i = 0; i < nplane; ++i) {
    c->cd(i + 1);
    const TString h_local =
      Form("D5WireResi_%s_LocalFit_plane%d", detector, i);
    const TString h_d5 =
      Form("D5WireResi_%s_D5Fit_plane%d", detector, i);
    auto* hl = dynamic_cast<TH1*>(get_object(f, h_local));
    auto* hd = dynamic_cast<TH1*>(get_object(f, h_d5));
    if (!hl && !hd) {
      draw_placeholder(Form("NOT FOUND:\nplane %d", i));
      continue;
    }
    if (hl) {
      hl->SetLineColor(kBlack);
      hl->Draw("hist");
    }
    if (hd) {
      hd->SetLineColor(kBlue + 1);
      if (hl) hd->Draw("hist same");
      else hd->Draw("hist");
      // D5Fit Gauss (same as LS input) + red dashed line at mean
      const GaussFit gf = fit_gauss_1p3(hd);
      if (gf.ok) {
        auto* gaus = hd->GetFunction(Form("gaus_%s", hd->GetName()));
        if (gaus) gaus->Draw("SAME");
        draw_gauss_mean_vline(hd, gf.mean);
      }
    }
    if (hl || hd) {
      auto* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
      if (hl) leg->AddEntry(hl, "LocalFit", "l");
      if (hd) leg->AddEntry(hd, "D5Fit", "l");
      leg->Draw();
    }
    if (i < 8) draw_plane_label(kPlaneLabels[i]);
  }
  print_page(c, pdf_path);
}

void
draw_resi_vs_s_page(TCanvas* c, const TString& pdf_path, TFile* f,
                    const char* detector, const char* s_tag, const char* s_label)
{
  // s_tag: "sperp" or "spara"; same r as D5Fit wire residual
  const Int_t nplane = conf.num_of_ch.at("blc");
  c->Clear();
  c->Divide(4, 2);
  for (Int_t i = 0; i < nplane && i < 8; ++i) {
    c->cd(i + 1);
    const TString hname =
      Form("D5WireResi_%s_D5Fit_plane%d_vs_%s", detector, i, s_tag);
    auto* h = dynamic_cast<TH2*>(get_object(f, hname));
    if (!h || h->GetEntries() <= 0) {
      draw_placeholder(Form("NOT FOUND / empty:\n%s\n(re-decode needed)",
                            hname.Data()));
      if (i < 8) draw_plane_label(kPlaneLabels[i]);
      continue;
    }
    gPad->SetLogz();
    h->Draw("colz");
    auto* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextAlign(11);
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.14, 0.88,
                   Form("%s D5Fit r vs %s", detector, s_label));
    tex->SetTextSize(0.028);
    tex->DrawLatex(0.14, 0.82, "(same r as D5Fit wire residual)");
    draw_plane_label(kPlaneLabels[i]);
  }
  print_page(c, pdf_path);
}

void
draw_xy_proj_and_summary(TCanvas* c, const TString& pdf_path, TFile* f,
                         const char* detector, const DetectorShift& sh)
{
  auto* h_rx = new TH1D(
    Form("D5WireResi_%s_D5Fit_rx", detector),
    Form("%s D5Fit r_{x}=r_{s} cos TA;r_{x} [mm];entries", detector),
    kProjBins, kProjLo, kProjHi);
  auto* h_ry = new TH1D(
    Form("D5WireResi_%s_D5Fit_ry", detector),
    Form("%s D5Fit r_{y}=r_{s} sin TA;r_{y} [mm];entries", detector),
    kProjBins, kProjLo, kProjHi);
  h_rx->Sumw2(kFALSE);
  h_ry->Sumw2(kFALSE);
  fill_rx_ry_from_planes(f, detector, h_rx, h_ry);

  // Page: rx, ry
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  if (h_rx->GetEntries() > 0) {
    h_rx->SetLineColor(kBlue + 1);
    h_rx->Draw("hist");
  } else {
    draw_placeholder("empty r_{x}");
  }
  c->cd(2);
  if (h_ry->GetEntries() > 0) {
    h_ry->SetLineColor(kBlue + 1);
    h_ry->Draw("hist");
  } else {
    draw_placeholder("empty r_{y}");
  }
  print_page(c, pdf_path);

  // Page: summary + recommended CLI
  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.10, 0.90, Form("%s D5Fit #rightarrow rigid (dx,dy) LS", detector));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.10, 0.84, "cos(TA) dx + sin(TA) dy #approx -#mu_{s}  (TA from DCGEO; same sign as update_residual)");

  Double_t y = 0.76;
  for (Int_t i = 0; i < 8; ++i) {
    if (!sh.mu_ok[i]) {
      tex->DrawLatex(0.10, y, Form("  %s: n/a", kPlaneLabels[i]));
    } else {
      tex->DrawLatex(0.10, y, Form("  %s: #mu_{s} = %+.4f mm", kPlaneLabels[i], sh.mu_s[i]));
    }
    y -= 0.045;
  }

  tex->SetTextSize(0.035);
  if (sh.ok) {
    tex->DrawLatex(0.10, 0.32, Form("LS: dx = %+.4f mm,  dy = %+.4f mm", sh.dx, sh.dy));
    tex->SetTextSize(0.028);
    const bool is_blc1 = (std::strncmp(detector, "BLC1", 4) == 0);
    const bool is_blc2 = (std::strncmp(detector, "BLC2", 4) == 0);
    if (is_blc1) {
      tex->DrawLatex(0.10, 0.24, "APPLY (a/b independent; TA fixed):");
      tex->DrawLatex(0.10, 0.18, "apply_dcgeo_offset.py <DCGEO> --skip-xy \\");
      tex->DrawLatex(0.10, 0.12,
                     Form("  --residual-offset %s:%.4f,%.4f", detector, sh.dx, sh.dy));
      tex->DrawLatex(0.10, 0.06, "(dry-run first; then --write). Do not use --bcin.");
    } else if (is_blc2) {
      tex->DrawLatex(0.10, 0.24, "BLC2 FIXED — do not apply this offset.");
      tex->DrawLatex(0.10, 0.18, "(values below are monitor-only)");
      tex->SetTextSize(0.024);
      tex->DrawLatex(0.10, 0.12,
                     Form("  [monitor] --residual-offset %s:%.4f,%.4f",
                          detector, sh.dx, sh.dy));
    } else {
      tex->DrawLatex(0.10, 0.24, "apply_dcgeo_offset.py <DCGEO> --skip-xy \\");
      tex->DrawLatex(0.10, 0.18,
                     Form("  --residual-offset %s:%.4f,%.4f", detector, sh.dx, sh.dy));
      tex->DrawLatex(0.10, 0.12, "(dry-run first; then --write)");
    }
  } else {
    tex->DrawLatex(0.10, 0.32, "LS failed (need >=2 valid planes / non-singular).");
  }
  print_page(c, pdf_path);
}

void
print_shift_log(const char* detector, const DetectorShift& sh)
{
  std::cout << "=== " << detector << " D5Fit -> (dx,dy) ===" << std::endl;
  for (Int_t i = 0; i < 8; ++i) {
    if (!sh.mu_ok[i])
      std::cout << "  " << kPlaneLabels[i] << ": n/a" << std::endl;
    else
      std::cout << "  " << kPlaneLabels[i] << ": mu_s=" << sh.mu_s[i] << " mm"
                << std::endl;
  }
  if (sh.ok) {
    std::cout << "  LS: dx=" << sh.dx << " mm, dy=" << sh.dy << " mm" << std::endl;
    if (std::strncmp(detector, "BLC1", 4) == 0) {
      std::cout << "  => APPLY: apply_dcgeo_offset.py <DCGEO> --skip-xy"
                << " --residual-offset " << detector << ":" << sh.dx << "," << sh.dy
                << std::endl;
    } else if (std::strncmp(detector, "BLC2", 4) == 0) {
      std::cout << "  => BLC2 FIXED (monitor only): --residual-offset "
                << detector << ":" << sh.dx << "," << sh.dy << std::endl;
    } else {
      std::cout << "  => apply_dcgeo_offset.py <DCGEO> --skip-xy"
                << " --residual-offset " << detector << ":" << sh.dx << "," << sh.dy
                << std::endl;
    }
  } else {
    std::cout << "  LS: failed" << std::endl;
  }
}

void
print_group_log(const GroupShift& g)
{
  std::cout << "=== " << g.group << " (a+b) D5Fit -> (dx,dy) ===" << std::endl;
  for (Int_t id = 0; id < 2; ++id) {
    const char* det = (id == 0) ? g.det_a : g.det_b;
    for (Int_t i = 0; i < 8; ++i) {
      const Int_t k = id * 8 + i;
      if (!g.mu_ok[k])
        std::cout << "  " << det << "-" << kPlaneLabels[i] << ": n/a" << std::endl;
      else
        std::cout << "  " << det << "-" << kPlaneLabels[i]
                  << ": mu_s=" << g.mu_s[k] << " mm" << std::endl;
    }
  }
  if (g.ok) {
    std::cout << "  LS(a+b): dx=" << g.dx << " mm, dy=" << g.dy << " mm" << std::endl;
    if (std::strcmp(g.group, "BLC1") == 0) {
      std::cout << "  => do NOT use --bcin; apply BLC1a and BLC1b"
                << " --residual-offset separately" << std::endl;
    } else {
      std::cout << "  => BLC2 FIXED; do NOT use --bcout" << std::endl;
    }
  } else {
    std::cout << "  LS(a+b): failed" << std::endl;
  }
}

void
draw_group_xy_and_summary(TCanvas* c, const TString& pdf_path, TFile* f,
                          const GroupShift& g)
{
  auto* h_rx = new TH1D(
    Form("D5WireResi_%s_D5Fit_rx", g.group),
    Form("%s (a+b) D5Fit r_{x};r_{x} [mm];entries", g.group),
    kProjBins, kGroupProjLo, kGroupProjHi);
  auto* h_ry = new TH1D(
    Form("D5WireResi_%s_D5Fit_ry", g.group),
    Form("%s (a+b) D5Fit r_{y};r_{y} [mm];entries", g.group),
    kProjBins, kGroupProjLo, kGroupProjHi);
  h_rx->Sumw2(kFALSE);
  h_ry->Sumw2(kFALSE);
  fill_rx_ry_from_planes(f, g.det_a, h_rx, h_ry);
  fill_rx_ry_from_planes(f, g.det_b, h_rx, h_ry);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  if (h_rx->GetEntries() > 0) {
    h_rx->SetLineColor(kBlue + 1);
    h_rx->Draw("hist");
  } else {
    draw_placeholder("empty r_{x}");
  }
  c->cd(2);
  if (h_ry->GetEntries() > 0) {
    h_ry->SetLineColor(kBlue + 1);
    h_ry->Draw("hist");
  } else {
    draw_placeholder("empty r_{y}");
  }
  print_page(c, pdf_path);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.08, 0.93,
                 Form("%s (a+b) D5Fit #rightarrow rigid (dx,dy) LS", g.group));
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.08, 0.88, "16 planes (a+b). cos dx + sin dy #approx -#mu_{s} (TA from DCGEO)");

  // Two columns: a | b  so text does not overflow the page
  for (Int_t id = 0; id < 2; ++id) {
    const char* det = (id == 0) ? g.det_a : g.det_b;
    const Double_t x0 = (id == 0) ? 0.08 : 0.52;
    Double_t y = 0.82;
    tex->SetTextSize(0.028);
    tex->DrawLatex(x0, y, Form("%s:", det));
    y -= 0.038;
    tex->SetTextSize(0.024);
    for (Int_t i = 0; i < 8; ++i) {
      const Int_t k = id * 8 + i;
      if (!g.mu_ok[k])
        tex->DrawLatex(x0, y, Form("  %s: n/a", kPlaneLabels[i]));
      else
        tex->DrawLatex(x0, y, Form("  %s: #mu_{s}=%+.4f mm", kPlaneLabels[i], g.mu_s[k]));
      y -= 0.038;
    }
  }

  tex->SetTextSize(0.032);
  if (g.ok) {
    tex->DrawLatex(0.08, 0.28, Form("LS(a+b) [info only]: dx = %+.4f mm,  dy = %+.4f mm",
                                    g.dx, g.dy));
    tex->SetTextSize(0.026);
    if (std::strcmp(g.group, "BLC1") == 0) {
      tex->DrawLatex(0.08, 0.20, "Do NOT use --bcin (a/b must be independent).");
      tex->DrawLatex(0.08, 0.14, "Apply --residual-offset on BLC1a and BLC1b pages / final CLI page.");
      tex->DrawLatex(0.08, 0.08, "TA fixed (no #DeltaTA from this review).");
    } else {
      tex->DrawLatex(0.08, 0.20, "BLC2 FIXED — do NOT use --bcout / residual-offset here.");
      tex->DrawLatex(0.08, 0.14, "This page is monitor-only (mean / slope checks).");
      tex->DrawLatex(0.08, 0.08, "TA fixed.");
    }
  } else {
    tex->DrawLatex(0.08, 0.28, "LS(a+b) failed.");
  }
  print_page(c, pdf_path);
}

void
draw_recommended_cli_page(TCanvas* c, const TString& pdf_path,
                          const DetectorShift& blc1a, const DetectorShift& blc1b,
                          const DetectorShift& blc2a, const DetectorShift& blc2b)
{
  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.08, 0.92, "Recommended apply (current policy)");
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.08, 0.86, "BLC1a / BLC1b independent Ofs | BLC2 fixed | TA fixed");

  tex->SetTextSize(0.032);
  tex->DrawLatex(0.08, 0.78, "Run (dry-run first; then --write):");
  tex->SetTextSize(0.026);
  if (blc1a.ok && blc1b.ok) {
    tex->DrawLatex(0.08, 0.70, "apply_dcgeo_offset.py <DCGEO> --skip-xy \\");
    tex->DrawLatex(0.08, 0.64,
                   Form("  --residual-offset BLC1a:%.4f,%.4f \\", blc1a.dx, blc1a.dy));
    tex->DrawLatex(0.08, 0.58,
                   Form("  --residual-offset BLC1b:%.4f,%.4f", blc1b.dx, blc1b.dy));
  } else {
    tex->DrawLatex(0.08, 0.70, "BLC1a and/or BLC1b LS failed — see per-detector pages.");
  }

  tex->SetTextSize(0.028);
  tex->DrawLatex(0.08, 0.48, "Do NOT use:");
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.08, 0.42, "  --bcin / --bcout  (forces a+b common shift)");
  tex->DrawLatex(0.08, 0.36, "  BLC2 --residual-offset  (BLC2 fixed)");
  tex->DrawLatex(0.08, 0.30, "  TA / RA edits from this review  (rotation fixed for now)");

  tex->SetTextSize(0.024);
  tex->DrawLatex(0.08, 0.20, Form("BLC2 monitor LS:  a dx=%+.4f dy=%+.4f | b dx=%+.4f dy=%+.4f",
                                  blc2a.ok ? blc2a.dx : 0.0, blc2a.ok ? blc2a.dy : 0.0,
                                  blc2b.ok ? blc2b.dx : 0.0, blc2b.ok ? blc2b.dy : 0.0));
  tex->DrawLatex(0.08, 0.14, "(BLC2 numbers are for watching drift only)");
  print_page(c, pdf_path);
}

void
draw_one_hist_pad(TFile* f, const char* name)
{
  auto* obj = get_object(f, name);
  if (!obj) {
    draw_placeholder(Form("NOT FOUND:\n%s", name));
    return;
  }
  if (auto* h2 = dynamic_cast<TH2*>(obj)) {
    if (h2->GetEntries() <= 0) {
      draw_placeholder(Form("empty:\n%s", name));
      return;
    }
    gPad->SetLogz();
    h2->Draw("colz");
    return;
  }
  if (auto* h1 = dynamic_cast<TH1*>(obj)) {
    if (h1->GetEntries() <= 0) {
      draw_placeholder(Form("empty:\n%s", name));
      return;
    }
    h1->SetLineColor(kBlue + 1);
    h1->Draw("hist");
    return;
  }
  draw_placeholder(Form("unsupported:\n%s", name));
}

void
draw_d5_track_qa_pages(TCanvas* c, const TString& pdf_path, TFile* f)
{
  // Basic D5Track_* QA (beam flag off -> no suffix). Same ROOT as wire residuals.
  const char* page1[] = {
    "D5Track_Delta",
    "D5Track_Momentum",
    "D5Track_D5Chi2",
    "D5Track_D5Chi2NDF",
  };
  c->Clear();
  c->Divide(2, 2);
  for (Int_t i = 0; i < 4; ++i) {
    c->cd(i + 1);
    draw_one_hist_pad(f, page1[i]);
  }
  print_page(c, pdf_path);

  const char* page2[] = {
    "D5Track_ResidualX",
    "D5Track_ResidualY",
    "D5Track_ResidualXY",
    "D5Track_FitX0",
  };
  c->Clear();
  c->Divide(2, 2);
  for (Int_t i = 0; i < 4; ++i) {
    c->cd(i + 1);
    draw_one_hist_pad(f, page2[i]);
  }
  print_page(c, pdf_path);

  const char* page3[] = {
    "D5Track_FitY0",
    "D5Track_MtxoutX",
    "D5Track_MtxoutY",
    "D5Track_FitY0VsBLC1y",
  };
  c->Clear();
  c->Divide(2, 2);
  for (Int_t i = 0; i < 4; ++i) {
    c->cd(i + 1);
    draw_one_hist_pad(f, page3[i]);
  }
  print_page(c, pdf_path);

  // FitX0 vs BLC1x on its own half-page with note
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_one_hist_pad(f, "D5Track_FitX0VsBLC1x");
  c->cd(2);
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.12, 0.70, "D5Track QA (end of review)");
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.12, 0.58, "Delta / Momentum / #chi^{2} / residuals / fit coords");
  tex->DrawLatex(0.12, 0.48, "Filled by UserD5Tracking (rank0 pairs).");
  tex->DrawLatex(0.12, 0.38, "Independent of wire-residual alignment pages.");
  print_page(c, pdf_path);
}

Int_t
analyze(const TString& input_path, const TString& output_path_in,
        const TString& dcgeo_path_in)
{
  setup_style();

  TString dcgeo_path = dcgeo_path_in.IsNull() ? default_dcgeo_path() : dcgeo_path_in;
  gSystem->ExpandPathName(dcgeo_path);
  if (!load_dcgeo_ta(dcgeo_path)) {
    std::cerr << "Warning: DCGEO TA load failed; falling back to nominal +/-45 deg.\n";
  }

  auto* f = TFile::Open(input_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open: " << input_path << std::endl;
    return 1;
  }

  Int_t run_num = -1;
  auto* d5 = dynamic_cast<TTree*>(f->Get("d5"));
  if (d5 && d5->GetEntries() > 0) {
    UInt_t run_number = 0;
    if (d5->GetBranch("run_number"))
      d5->SetBranchAddress("run_number", &run_number);
    d5->GetEntry(0);
    run_num = static_cast<Int_t>(run_number);
  }
  if (run_num < 0) run_num = parse_run_from_path(input_path);

  TString output_path = output_path_in;
  if (output_path.IsNull()) {
    output_path = default_pdf_path(run_num);
    std::cout << "Default output: " << output_path << std::endl;
  }

  if (hist_entries(f, "D5WireResi_BLC1a_LocalFit_vs_plane") <= 0) {
    std::cerr << "Warning: D5WireResi_* histograms not found or empty.\n"
              << "  Re-decode with updated UserD5Tracking.\n";
  }

  DetectorShift shifts[kNumDetectors];
  for (Int_t id = 0; id < kNumDetectors; ++id) {
    shifts[id] = estimate_detector_shift(f, kDetectors[id]);
    print_shift_log(kDetectors[id], shifts[id]);
  }

  const GroupShift g_blc1 =
    estimate_group_shift(f, "BLC1", "--bcin", "BLC1a", "BLC1b");
  const GroupShift g_blc2 =
    estimate_group_shift(f, "BLC2", "--bcout", "BLC2a", "BLC2b");
  print_group_log(g_blc1);
  print_group_log(g_blc2);

  auto* canvas = new TCanvas("d5_wire_resi_review", output_path, 1600, 1200);
  canvas->Print(output_path + "[");
  draw_title_page(canvas, output_path, f, run_num);
  for (Int_t id = 0; id < kNumDetectors; ++id) {
    draw_vs_plane_page(canvas, output_path, f, kDetectors[id]);
    draw_plane_overlay_page(canvas, output_path, f, kDetectors[id]);
    draw_resi_vs_s_page(canvas, output_path, f, kDetectors[id], "sperp",
                        "s_{perp}");
    draw_resi_vs_s_page(canvas, output_path, f, kDetectors[id], "spara",
                        "s_{para}");
    draw_xy_proj_and_summary(canvas, output_path, f, kDetectors[id], shifts[id]);
  }
  draw_group_xy_and_summary(canvas, output_path, f, g_blc1);
  draw_group_xy_and_summary(canvas, output_path, f, g_blc2);
  draw_recommended_cli_page(canvas, output_path,
                            shifts[0], shifts[1], shifts[2], shifts[3]);
  draw_d5_track_qa_pages(canvas, output_path, f);
  canvas->Print(output_path + "]");
  delete canvas;
  f->Close();

  std::cout << "Wrote " << output_path << " (run=" << run_num << ")" << std::endl;
  return 0;
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  TString input_path;
  TString output_path;
  TString dcgeo_path;
  bool has_output = false;

  for (Int_t i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    if (strcmp(arg, "-o") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      output_path = argv[++i];
      has_output = true;
      continue;
    }
    if (strcmp(arg, "--dcgeo") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      dcgeo_path = argv[++i];
      continue;
    }
    if (arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
    if (!input_path.IsNull()) {
      std::cerr << "Multiple input files specified." << std::endl;
      usage(argv[0]);
      return 1;
    }
    input_path = arg;
  }

  if (input_path.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  gSystem->ExpandPathName(input_path);
  return analyze(input_path, has_output ? output_path : TString(), dcgeo_path);
}
