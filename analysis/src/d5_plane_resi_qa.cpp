// Per-plane residual QA for the D5 alignment (BLC1a/b, BLC2a/b), K window versus pi
// Usage: d5_plane_resi_qa <run> [--bcin <BcIn.root>] [--bcout <BcOut.root>] [--d5 <D5.root>]
//                         [-o <out.pdf>]
// Defaults: BcIn  DECODE_DIR/runXXXXX/runXXXXX_BcIn_resi_K.root
//           BcOut DECODE_DIR/runXXXXX/runXXXXX_BcOut_resi_K.root
//           D5    DECODE_DIR/runXXXXX/runXXXXX_D5_K.root
//           pdf   OUTPUT_DIR/img/runXXXXX/d5_plane_resi_qa_runXXXXX.pdf
//
// Chamber-fit residuals (BcIn: BLC1a+b, BcOut: BLC2a+b) are drawn per plane for the K
// window (*_K) and pi (*_Pi), pi scaled to the K entries. <|r|> (within +-kRange) gives
// the half distance of the two peaks when the residual is split by the left/right side. The D5 joint-fit residuals
// (D5WireResi_*_D5Fit_plane*) are filled for every beam flag in UserD5Tracking and are
// drawn as they are (all beam) for comparison. The first pages show residual versus plane
// (2D, linear colour scale, profile mean overlaid) for K, pi and the D5 fit.

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{

const std::vector<TString> kDets = { "BLC1a", "BLC1b", "BLC2a", "BLC2b" };
const std::vector<TString> kPlanes = { "U1", "UP1", "V1", "VP1", "U2", "UP2", "V2", "VP2" };
constexpr Double_t kRange = 1.2; // [mm] drawn residual range

struct Core
{
  Double_t mean = NAN;
  Double_t sigma = NAN;
};

// Gauss core in +-1.5 sigma, iterated from the peak
Core
fit_core(TH1* h, Double_t sigma0)
{
  Core c;
  if (!h || h->GetEntries() < 100) return c;
  Double_t m = h->GetBinCenter(h->GetMaximumBin()), s = sigma0;
  auto* f = new TF1(Form("g_%s", h->GetName()), "gaus", m - 1.5 * s, m + 1.5 * s);
  f->SetLineColor(kBlack);
  f->SetLineWidth(1);
  for (Int_t it = 0; it < 5; ++it) {
    f->SetRange(m - 1.5 * s, m + 1.5 * s);
    f->SetParameters(h->GetMaximum(), m, s);
    h->Fit(f, "QR0");
    m = f->GetParameter(1);
    s = std::fabs(f->GetParameter(2));
  }
  c.mean = m;
  c.sigma = s;
  return c;
}

// <|r|> within +-kRange: half distance of the two peaks of a left/right split residual
Double_t
abs_mean(const TH1* h)
{
  Double_t sw = 0.0, sx = 0.0;
  for (Int_t i = 1; i <= h->GetNbinsX(); ++i) {
    const Double_t x = h->GetBinCenter(i);
    if (std::fabs(x) > kRange) continue;
    sw += h->GetBinContent(i);
    sx += h->GetBinContent(i) * std::fabs(x);
  }
  return sw > 0.0 ? sx / sw : NAN;
}

TH1*
get_clone(TFile* f, const TString& name)
{
  if (!f) return nullptr;
  auto* h = dynamic_cast<TH1*>(f->Get(name));
  if (!h) return nullptr;
  auto* c = dynamic_cast<TH1*>(h->Clone(name + "_qa"));
  c->SetDirectory(nullptr);
  return c;
}

// residual vs plane from the per-plane 1D histograms (x: plane 0-7, y: residual)
TH2D*
stack_planes(TFile* f, const TString& det, const TString& suffix, Double_t range)
{
  auto* h0 = dynamic_cast<TH1*>(f->Get(Form("%s_Track_Residual_plane0%s", det.Data(), suffix.Data())));
  if (!h0) return nullptr;
  const Int_t lo = h0->FindBin(-range), hi = h0->FindBin(range);
  auto* h2 = new TH2D(Form("h2_%s%s", det.Data(), suffix.Data()), "", 8, -0.5, 7.5, hi - lo + 1,
                      h0->GetBinLowEdge(lo), h0->GetBinLowEdge(hi + 1));
  h2->SetDirectory(nullptr);
  for (Int_t p = 0; p < 8; ++p) {
    auto* h = dynamic_cast<TH1*>(f->Get(Form("%s_Track_Residual_plane%d%s", det.Data(), p, suffix.Data())));
    if (!h) continue;
    for (Int_t b = lo; b <= hi; ++b)
      h2->SetBinContent(p + 1, b - lo + 1, h->GetBinContent(b));
  }
  for (Int_t p = 0; p < 8; ++p) h2->GetXaxis()->SetBinLabel(p + 1, kPlanes[p]);
  return h2;
}

void
draw_vs_plane(TH2* h2, const TString& title)
{
  h2->SetTitle(title + ";plane;residual [mm]");
  h2->GetXaxis()->SetLabelSize(0.06);
  h2->Draw("colz");
  auto* pr = h2->ProfileX(Form("%s_pfx", h2->GetName()));
  pr->SetLineColor(kBlack);
  pr->SetMarkerColor(kBlack);
  pr->SetMarkerStyle(20);
  pr->SetMarkerSize(0.6);
  pr->Draw("same");
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " <run> [--bcin <file>] [--bcout <file>] [--d5 <file>]"
            << " [-o <out.pdf>]\n";
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  if (argc < 2) { usage(argv[0]); return 1; }
  const Int_t run = std::atoi(argv[1]);
  TString path_bcin = Form("%s/run%05d/run%05d_BcIn_resi_K.root", DECODE_DIR.Data(), run, run);
  TString path_bcout = Form("%s/run%05d/run%05d_BcOut_resi_K.root", DECODE_DIR.Data(), run, run);
  TString path_d5 = Form("%s/run%05d/run%05d_D5_K.root", DECODE_DIR.Data(), run, run);
  TString path_pdf = Form("%s/d5_plane_resi_qa_run%05d.pdf",
                          ana_helper::get_img_dir(OUTPUT_DIR, run).Data(), run);
  for (Int_t i = 2; i < argc; ++i) {
    const TString a = argv[i];
    if (a == "--bcin" && i + 1 < argc) path_bcin = argv[++i];
    else if (a == "--bcout" && i + 1 < argc) path_bcout = argv[++i];
    else if (a == "--d5" && i + 1 < argc) path_d5 = argv[++i];
    else if (a == "-o" && i + 1 < argc) path_pdf = argv[++i];
    else { usage(argv[0]); return 1; }
  }

  auto* f_bcin = TFile::Open(path_bcin);
  auto* f_bcout = TFile::Open(path_bcout);
  auto* f_d5 = TFile::Open(path_d5);
  if (!f_bcin || !f_bcout || !f_d5) {
    std::cerr << "cannot open input files" << std::endl;
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetPaperSize(29.0, 20.5);

  // core[det][plane]: K, pi, D5Fit(all)
  Core c_k[4][8], c_d5[4][8];
  Double_t abs_k[4][8], abs_pi[4][8];
  for (auto& a : abs_k) for (auto& v : a) v = NAN;
  for (auto& a : abs_pi) for (auto& v : a) v = NAN;

  auto* c = new TCanvas("c", "c", 1600, 1130);
  c->Print(path_pdf + "[");

  // residual vs plane (2D): K window, pi (chamber fit) and D5 fit (all beam)
  for (const TString& kind : { TString("_K"), TString("_Pi"), TString("D5Fit") }) {
    c->Clear();
    c->Divide(2, 2, 0.002, 0.002);
    for (std::size_t d = 0; d < kDets.size(); ++d) {
      c->cd(d + 1);
      gPad->SetRightMargin(0.12);
      TH2* h2 = nullptr;
      TString title;
      if (kind == "D5Fit") {
        auto* src = dynamic_cast<TH2*>(f_d5->Get(Form("D5WireResi_%s_D5Fit_vs_plane", kDets[d].Data())));
        if (!src) continue;
        h2 = dynamic_cast<TH2*>(src->Clone(Form("h2_%s_d5", kDets[d].Data())));
        h2->SetDirectory(nullptr);
        h2->GetYaxis()->SetRangeUser(-1.5, 1.5);
        for (Int_t p = 0; p < 8; ++p) h2->GetXaxis()->SetBinLabel(p + 1, kPlanes[p]);
        title = Form("%s D5 fit (all beam)", kDets[d].Data());
      } else {
        h2 = stack_planes((d < 2) ? f_bcin : f_bcout, kDets[d], kind, 0.8);
        title = Form("%s chamber fit, %s", kDets[d].Data(), kind == "_K" ? "K window" : "#pi");
      }
      if (h2) draw_vs_plane(h2, title);
    }
    c->cd();
    c->Print(path_pdf, Form("Title:residual vs plane %s", kind.Data()));
  }

  for (std::size_t d = 0; d < kDets.size(); ++d) {
    TFile* f_track = (d < 2) ? f_bcin : f_bcout;
    c->Clear();
    c->Divide(4, 2, 0.002, 0.002);
    for (Int_t p = 0; p < 8; ++p) {
      c->cd(p + 1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.03);
      auto* hk = get_clone(f_track, Form("%s_Track_Residual_plane%d_K", kDets[d].Data(), p));
      auto* hp = get_clone(f_track, Form("%s_Track_Residual_plane%d_Pi", kDets[d].Data(), p));
      if (!hk || !hp) continue;
      c_k[d][p] = fit_core(hk, 0.1);
      abs_k[d][p] = abs_mean(hk);
      abs_pi[d][p] = abs_mean(hp);
      hk->Rebin(2);
      hp->Rebin(2);
      if (hp->Integral() > 0) hp->Scale(hk->Integral() / hp->Integral());
      hk->SetTitle(Form("%s-%s (chamber fit);residual [mm];count", kDets[d].Data(), kPlanes[p].Data()));
      hk->GetXaxis()->SetRangeUser(-kRange, kRange);
      hk->SetMaximum(1.15 * std::max(hk->GetMaximum(), hp->GetMaximum()));
      hk->SetLineColor(kBlack);
      hk->SetFillColorAlpha(kGray, 0.5);
      hp->SetLineColor(kRed + 1);
      hp->SetLineWidth(2);
      hk->Draw("hist");
      hp->Draw("hist same");
      TLatex t;
      t.SetNDC();
      t.SetTextSize(0.05);
      t.DrawLatex(0.15, 0.84, Form("K: core %+.3f / #sigma %.3f", c_k[d][p].mean, c_k[d][p].sigma));
      t.SetTextColor(kRed + 1);
      t.DrawLatex(0.15, 0.78, Form("#pi: <|r|> %.3f (K %.3f)", abs_pi[d][p], abs_k[d][p]));
      if (p == 0) {
        auto* leg = new TLegend(0.62, 0.70, 0.97, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hk, "K window", "f");
        leg->AddEntry(hp, "#pi (scaled)", "l");
        leg->Draw();
      }
    }
    c->cd();
    c->Print(path_pdf, Form("Title:%s chamber fit K vs pi", kDets[d].Data()));
  }

  for (std::size_t d = 0; d < kDets.size(); ++d) {
    c->Clear();
    c->Divide(4, 2, 0.002, 0.002);
    for (Int_t p = 0; p < 8; ++p) {
      c->cd(p + 1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.03);
      auto* hd = get_clone(f_d5, Form("D5WireResi_%s_D5Fit_plane%d", kDets[d].Data(), p));
      if (!hd) continue;
      c_d5[d][p] = fit_core(hd, 0.3);
      hd->SetTitle(Form("%s-%s (D5 fit, all beam);residual [mm];count", kDets[d].Data(), kPlanes[p].Data()));
      hd->GetXaxis()->SetRangeUser(-2.0, 2.0);
      hd->SetLineColor(kBlue + 1);
      hd->SetFillColorAlpha(kAzure - 9, 0.5);
      hd->Draw("hist");
      TLatex t;
      t.SetNDC();
      t.SetTextSize(0.05);
      t.DrawLatex(0.15, 0.84, Form("core %+.3f / #sigma %.3f", c_d5[d][p].mean, c_d5[d][p].sigma));
    }
    c->cd();
    c->Print(path_pdf, Form("Title:%s D5 fit", kDets[d].Data()));
  }

  // summary of the core means
  c->Clear();
  c->cd();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(82);
  t.SetTextSize(0.016);
  Double_t y = 0.95;
  t.DrawLatex(0.03, y, Form("run %05d  per-plane residual core mean / sigma [mm]", run));
  y -= 0.03;
  t.DrawLatex(0.03, y, Form("chamber fit: %s, %s", gSystem->BaseName(path_bcin), gSystem->BaseName(path_bcout)));
  y -= 0.025;
  t.DrawLatex(0.03, y, Form("D5 fit (all beam): %s", gSystem->BaseName(path_d5)));
  y -= 0.035;
  t.DrawLatex(0.03, y, "plane          K mean  K sig   K <|r|> pi <|r|> D5fit mean D5fit sig");
  y -= 0.022;
  for (std::size_t d = 0; d < kDets.size(); ++d) {
    for (Int_t p = 0; p < 8; ++p) {
      t.DrawLatex(0.03, y, Form("%-5s-%-4s    %+.3f  %.3f   %.3f   %.3f    %+.3f     %.3f",
                                kDets[d].Data(), kPlanes[p].Data(), c_k[d][p].mean, c_k[d][p].sigma,
                                abs_k[d][p], abs_pi[d][p], c_d5[d][p].mean, c_d5[d][p].sigma));
      y -= 0.0255;
    }
  }
  c->Print(path_pdf, "Title:summary");
  c->Print(path_pdf + "]");

  std::cout << "wrote " << path_pdf << std::endl;
  return 0;
}
