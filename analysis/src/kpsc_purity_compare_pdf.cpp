// -*- C++ -*-
// Kp scattering purity review: multipage PDF from DstKpScattering output.
//
// Usage:
//   kpsc_purity_compare_pdf <kpsc.root> [-o <out.pdf>] [-r <run>]
//
// Input:
//   TTree "kpsc" from bin/DstKpScattering output.
//
// Output (default):
//   OUTPUT_DIR/img/runXXXXX/kpsc_purity_compare_runXXXXX.pdf

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr Double_t kMK = 0.494;
constexpr Double_t kMmassWindowDiag = 0.15;
constexpr Double_t kMmassWindowTight = 0.10;
constexpr Double_t kAngleMissKGood = 0.3;
constexpr Double_t kSidebandLo = 0.65;
constexpr Double_t kSidebandHi = 0.85;
constexpr Double_t kNsigmaPionWindow = 3.0;
constexpr Double_t kNsigmaKaonReject = -2.0;
constexpr Double_t kNsigmaProtonReject = 3.0;

constexpr const char* kCutBase = "mmass>0";
constexpr const char* kCutB2 = "mmass>0 && effective_ntTpc==2";
constexpr const char* kCutB3 = "mmass>0 && effective_ntTpc<=3";
constexpr const char* kCutB4 = "mmass>0 && effective_ntTpc<=4";
constexpr const char* kCutC2 =
  "mmass>0 && effective_ntTpc==2 && abs(mmass-0.494)<0.15";
constexpr const char* kCutC3 =
  "mmass>0 && effective_ntTpc<=3 && abs(mmass-0.494)<0.15";
constexpr const char* kCutC4 =
  "mmass>0 && effective_ntTpc<=4 && abs(mmass-0.494)<0.15";
constexpr const char* kCutGood =
  "mmass>0 && effective_ntTpc==2 && angle_miss_k<0.3 && abs(mmass-0.494)<0.1";
constexpr const char* kCutWorking =
  "mmass>0 && effective_ntTpc==2 && abs(mmass-0.494)<0.15"
  " && nsigma_k_sel>-2 && nsigma_k_sel<2"
  " && !(abs(nsigma_pi_sel)<3 && nsigma_k_sel<-2)"
  " && nsigma_p_sel>-1.5";
constexpr const char* kCutSidebandNt2 =
  "mmass>0 && effective_ntTpc==2 && mmass>0.65 && mmass<0.85";
constexpr const char* kCutSidebandNt4 =
  "mmass>0 && effective_ntTpc<=4 && mmass>0.65 && mmass<0.85";
constexpr const char* kCutPiVeto =
  "mmass>0 && effective_ntTpc==2 && abs(mmass-0.494)<0.15"
  " && !(abs(nsigma_pi_sel)<3 && nsigma_k_sel<-2)"
  " && !(abs(nsigma_pi_p_sel)<3 && nsigma_p_sel>3)";

struct CutSet
{
  const char* name;
  const char* label;
  const char* cut;
  Color_t color;
};

const CutSet kSummaryCuts[] = {
  {"A", "A: kinematics OK", kCutBase, kBlack},
  {"B2", "B2: Nt==2", kCutB2, kBlue + 1},
  {"B3", "B3: Nt<=3", kCutB3, kAzure + 2},
  {"B4", "B4: Nt<=4", kCutB4, kCyan + 2},
  {"C2", "C2: Nt2+MMKW", kCutC2, kRed + 1},
  {"C3", "C3: Nt<=3+MMKW", kCutC3, kOrange + 7},
  {"C4", "C4: Nt<=4+MMKW", kCutC4, kMagenta + 1},
  {"G", "G: good (AM+MM0.1)", kCutGood, kOrange + 2},
  {"D", "D: working", kCutWorking, kGreen + 2},
};

struct Var1D
{
  const char* var;
  const char* title;
  Int_t nb;
  Double_t xlo;
  Double_t xhi;
};

struct Var2D
{
  const char* var;
  const char* xtitle;
  Int_t nbx;
  Double_t xlo;
  Double_t xhi;
};

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " <kpsc.root> [-o <out.pdf>] [-r <run>]\n"
            << "  Reads TTree kpsc.\n"
            << "  Default PDF: OUTPUT_DIR/img/runXXXXX/kpsc_purity_compare_runXXXXX.pdf\n";
}

Int_t
ParseRunFromTree(TTree* t)
{
  if (!t || t->GetEntries() <= 0)
    return -1;
  if (!t->GetBranch("run_number"))
    return -1;
  t->GetEntry(0);
  return static_cast<Int_t>(t->GetMaximum("run_number"));
}

Int_t
ParseRunFromPath(const TString& path)
{
  TString s = path;
  s.ToLower();
  const Ssiz_t n = s.Length();
  for (Ssiz_t i = 0; i + 3 < n; ++i) {
    if (s(i, 3) != "run")
      continue;
    const Ssiz_t j = i + 3;
    if (j >= n)
      continue;
    const char c = s[j];
    if (c < '0' || c > '9')
      continue;
    Ssiz_t k = j;
    while (k < n) {
      const char d = s[k];
      if (d < '0' || d > '9')
        break;
      ++k;
    }
    if (k - j > 6)
      continue;
    return TString(s(j, k - j)).Atoi();
  }
  return 0;
}

void
StyleCanvas(TCanvas& c, Bool_t wideRight = false)
{
  c.SetLeftMargin(0.12);
  c.SetRightMargin(wideRight ? 0.15 : 0.05);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.10);
}

TH1D*
Book1D(const char* name, const char* title, Int_t nb, Double_t xlo, Double_t xhi)
{
  auto* h = new TH1D(name, title, nb, xlo, xhi);
  h->SetDirectory(nullptr);
  h->Sumw2();
  return h;
}

TH1D*
Fill1D(TTree* t, const char* var, const char* cut,
       const char* name, const char* title, Int_t nb, Double_t xlo, Double_t xhi)
{
  if (auto* old = gDirectory->Get(name))
    delete old;
  auto* h = new TH1D(name, title, nb, xlo, xhi);
  t->Project(name, var, cut);
  auto* hc = static_cast<TH1D*>(h->Clone(Form("%s_c", name)));
  hc->SetDirectory(nullptr);
  hc->SetTitle(title);
  delete h;
  return hc;
}

void
PreparePalette()
{
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
}

void
DrawTree2D(TTree* t, TCanvas& c, const char* yvar, const char* xvar,
             const char* cut, const char* name, const char* title,
             Int_t nbx, Double_t xlo, Double_t xhi,
             Int_t nby, Double_t ylo, Double_t yhi, const char* opt)
{
  if (auto* old = gDirectory->Get(name))
    delete old;
  PreparePalette();
  const TString drawSpec = Form("%s:%s>>%s(%d,%g,%g,%d,%g,%g)",
                              yvar, xvar, name, nbx, xlo, xhi, nby, ylo, yhi);
  t->Draw(drawSpec, cut, opt);
  auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(name));
  if (h) {
    h->SetTitle(title);
    if (h->GetEntries() > 0. && h->GetMaximum() > 0.) {
      h->SetMinimum(0.);
      h->SetMaximum(h->GetMaximum() * 1.05);
    }
  }
  gPad->Modified();
  gPad->Update();
}

Double_t
HistMax(const std::vector<TH1D*>& hs)
{
  Double_t ymax = 0.;
  for (const auto* h : hs) {
    if (!h)
      continue;
    ymax = std::max(ymax, h->GetMaximum());
  }
  return ymax;
}

Double_t
CountEvents(TTree* t, const char* cut)
{
  t->Draw(">>h_count_tmp", cut, "goff");
  return static_cast<Double_t>(t->GetSelectedRows());
}

class PdfWriter
{
public:
  explicit PdfWriter(TString path) : path_(std::move(path)) {}

  void Print(TCanvas& c)
  {
    ++page_;
    if (page_ == 1)
      c.Print(path_ + "(");
    else
      c.Print(path_);
  }

  void Close(TCanvas& c)
  {
    if (page_ <= 0)
      return;
    if (page_ == 1)
      c.Print(path_);
    else
      c.Print(path_ + ")");
  }

  Int_t Pages() const { return page_; }

private:
  TString path_;
  Int_t page_ = 0;
};

void
DrawMmassBands(const TH1* ref)
{
  if (!ref)
    return;
  const Double_t xmin = ref->GetXaxis()->GetXmin();
  const Double_t xmax = ref->GetXaxis()->GetXmax();
  const Double_t peakLo = kMK - kMmassWindowDiag;
  const Double_t peakHi = kMK + kMmassWindowDiag;
  auto* lPeakLo = new TLine(xmin, peakLo, xmax, peakLo);
  auto* lPeakHi = new TLine(xmin, peakHi, xmax, peakHi);
  auto* lSbLo = new TLine(xmin, kSidebandLo, xmax, kSidebandLo);
  auto* lSbHi = new TLine(xmin, kSidebandHi, xmax, kSidebandHi);
  for (auto* l : {lPeakLo, lPeakHi, lSbLo, lSbHi}) {
    l->SetLineColor(kRed + 1);
    l->SetLineStyle(l == lSbLo || l == lSbHi ? 7 : 1);
    l->SetLineWidth(1);
    l->Draw("same");
  }
}

void
DrawSummaryTable(TTree* t, TCanvas& c, PdfWriter& writer)
{
  c.Clear();
  c.cd();
  StyleCanvas(c);

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.032);
  tx->DrawLatex(0.08, 0.94, "Kp purity cut summary");

  Double_t y = 0.88;
  for (const auto& cs : kSummaryCuts) {
    const Double_t n = CountEvents(t, cs.cut);
    tx->SetTextColor(cs.color);
    tx->SetTextSize(0.028);
    tx->DrawLatex(0.08, y, Form("%s : %.0f events", cs.label, n));
    y -= 0.055;
  }

  tx->SetTextColor(kGray + 2);
  tx->DrawLatex(0.08, y, Form("sideband Nt2 : %.0f", CountEvents(t, kCutSidebandNt2)));
  y -= 0.045;
  tx->DrawLatex(0.08, y, Form("sideband Nt<=4 : %.0f", CountEvents(t, kCutSidebandNt4)));
  y -= 0.045;
  tx->SetTextColor(kGreen + 2);
  tx->DrawLatex(0.08, y, Form("Nt2+MMKW + pi veto : %.0f", CountEvents(t, kCutPiVeto)));

  tx->SetTextColor(kBlack);
  tx->SetTextSize(0.022);
  tx->DrawLatex(0.08, 0.16, Form("MMKW: |mmass - %.3f| < %.2f GeV (diag) / %.2f (good)",
                                 kMK, kMmassWindowDiag, kMmassWindowTight));
  tx->DrawLatex(0.08, 0.12, Form("good cut: Nt==2 && angle_miss_k<%.1f && |mmass-m_K|<%.2f",
                                 kAngleMissKGood, kMmassWindowTight));
  tx->DrawLatex(0.08, 0.08, Form("sideband: %.2f < mmass < %.2f GeV", kSidebandLo, kSidebandHi));
  tx->DrawLatex(0.08, 0.04, Form("#pi-like K: |n#sigma_{#pi}|<%.0f && n#sigma_{K}<%.0f",
                                 kNsigmaPionWindow, kNsigmaKaonReject));

  writer.Print(c);
}

void
DrawOverlay1D(TTree* t, TCanvas& c, PdfWriter& writer,
                const char* var, const char* xtitle, Int_t nb,
                Double_t xlo, Double_t xhi, const char* pageTitle,
                const std::vector<CutSet>& cuts)
{
  c.Clear();
  c.cd();
  StyleCanvas(c);

  std::vector<TH1D*> hs;
  hs.reserve(cuts.size());
  for (const auto& cs : cuts) {
    const TString hname = Form("h_%s_%s", var, cs.name);
    hs.push_back(Fill1D(t, var, cs.cut, hname, xtitle, nb, xlo, xhi));
    hs.back()->SetLineColor(cs.color);
    hs.back()->SetLineWidth(2);
  }

  const Double_t ymax = HistMax(hs);
  hs.front()->SetMaximum(ymax > 0. ? ymax * 1.30 : 1.);
  hs.front()->SetTitle(pageTitle);
  hs.front()->GetXaxis()->SetTitle(xtitle);
  hs.front()->GetYaxis()->SetTitle("Counts");
  hs.front()->Draw("hist");
  for (std::size_t i = 1; i < hs.size(); ++i)
    hs[i]->Draw("hist same");

  auto* leg = new TLegend(0.48, 0.55, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  for (std::size_t i = 0; i < hs.size(); ++i)
    leg->AddEntry(hs[i], Form("%s (%0.f)", cuts[i].label, hs[i]->GetEntries()), "l");
  leg->Draw();

  writer.Print(c);
}

void
DrawPair1D(TTree* t, TCanvas& c, PdfWriter& writer,
           const char* var, const char* xtitle,
           const char* cutLeft, const char* labelLeft, Color_t colorLeft,
           const char* cutRight, const char* labelRight, Color_t colorRight,
           Int_t nb, Double_t xlo, Double_t xhi, const char* pageTitle)
{
  static Int_t pairIdx = 0;
  ++pairIdx;

  c.Clear();
  c.cd();
  StyleCanvas(c);

  const TString h1name = Form("h_pair_%d_%s_L", pairIdx, var);
  const TString h2name = Form("h_pair_%d_%s_R", pairIdx, var);
  auto* h1 = Fill1D(t, var, cutLeft, h1name, xtitle, nb, xlo, xhi);
  auto* h2 = Fill1D(t, var, cutRight, h2name, xtitle, nb, xlo, xhi);
  h1->SetLineColor(colorLeft);
  h2->SetLineColor(colorRight);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);

  const Double_t ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
  h1->SetMaximum(ymax > 0. ? ymax * 1.25 : 1.);
  h1->SetTitle(pageTitle);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->GetYaxis()->SetTitle("Counts");
  h1->Draw("hist");
  h2->Draw("hist same");

  auto* leg = new TLegend(0.52, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h1, Form("%s (%0.f)", labelLeft, h1->GetEntries()), "l");
  leg->AddEntry(h2, Form("%s (%0.f)", labelRight, h2->GetEntries()), "l");
  leg->Draw();

  writer.Print(c);
}

void
DrawMmass2D(TTree* t, TCanvas& c, PdfWriter& writer,
            const Var2D& v, const char* baseCut, const char* pageTitle)
{
  static Int_t mm2dIdx = 0;
  ++mm2dIdx;

  c.Clear();
  c.cd();
  StyleCanvas(c, true);

  const TString hname = Form("h_mmass_vs_%s_%d", v.var, mm2dIdx);
  const TString title = Form("%s; %s;missing mass [GeV]", pageTitle, v.xtitle);
  DrawTree2D(t, c, "mmass", v.var, baseCut, hname, title,
             v.nbx, v.xlo, v.xhi, 120, 0., 1.2, "COL");
  if (auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(hname)))
    DrawMmassBands(h);

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.028);
  if (auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(hname)))
    tx->DrawLatex(0.14, 0.92, Form("entries = %.0f", h->GetEntries()));
  tx->DrawLatex(0.14, 0.87, baseCut);

  writer.Print(c);
}

void
DrawNtMmassOverlay(TTree* t, TCanvas& c, PdfWriter& writer)
{
  c.Clear();
  c.cd();
  StyleCanvas(c);

  struct NtCut { const char* name; const char* cut; Color_t color; };
  const NtCut ntCuts[] = {
    {"nt2", "mmass>0 && effective_ntTpc==2", kBlue + 1},
    {"nt3", "mmass>0 && effective_ntTpc==3", kOrange + 7},
    {"nt4", "mmass>0 && effective_ntTpc==4", kMagenta + 1},
    {"ntle4", "mmass>0 && effective_ntTpc<=4", kGreen + 2},
  };

  std::vector<TH1D*> hs;
  for (const auto& nc : ntCuts) {
    auto* h = Fill1D(t, "mmass", nc.cut, Form("h_mmass_%s", nc.name),
                     "missing mass [GeV]", 120, 0., 1.2);
    h->SetLineColor(nc.color);
    h->SetLineWidth(2);
    hs.push_back(h);
  }

  const Double_t ymax = HistMax(hs);
  hs.front()->SetMaximum(ymax > 0. ? ymax * 1.25 : 1.);
  hs.front()->SetTitle("Missing mass by effective_ntTpc");
  hs.front()->GetXaxis()->SetTitle("missing mass [GeV]");
  hs.front()->Draw("hist");
  for (std::size_t i = 1; i < hs.size(); ++i)
    hs[i]->Draw("hist same");

  auto* leg = new TLegend(0.52, 0.62, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hs[0], Form("Nt==2 (%0.f)", hs[0]->GetEntries()), "l");
  leg->AddEntry(hs[1], Form("Nt==3 (%0.f)", hs[1]->GetEntries()), "l");
  leg->AddEntry(hs[2], Form("Nt==4 (%0.f)", hs[2]->GetEntries()), "l");
  leg->AddEntry(hs[3], Form("Nt<=4 (%0.f)", hs[3]->GetEntries()), "l");
  leg->Draw();

  writer.Print(c);
}

void
DrawNsigma2D(TTree* t, TCanvas& c, PdfWriter& writer,
             const char* cut, const char* title)
{
  static Int_t ns2dIdx = 0;
  ++ns2dIdx;

  c.Clear();
  c.cd();
  StyleCanvas(c, true);

  const TString hname = Form("h_nsigma_k_vs_pi_%s_%d", title, ns2dIdx);
  const TString pageTitle = Form("%s;n#sigma_{#pi}^{sel};n#sigma_{K}^{sel}", title);
  DrawTree2D(t, c, "nsigma_k_sel", "nsigma_pi_sel", cut, hname, pageTitle,
             100, -10., 10., 100, -10., 10., "COL");

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.030);
  if (auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(hname)))
    tx->DrawLatex(0.12, 0.92, Form("entries = %.0f", h->GetEntries()));
  tx->DrawLatex(0.12, 0.87, cut);

  writer.Print(c);
}

void
DrawPidDedxVsQp(TTree* t, TCanvas& c, PdfWriter& writer,
                const char* cut, const char* pageTitle, const char* htag,
                Bool_t allTracks)
{
  static Int_t pidIdx = 0;
  ++pidIdx;

  c.Clear();
  c.cd();
  StyleCanvas(c, true);
  PreparePalette();

  // Input TFile is READ; book TH2D ourselves (TTree::Draw creates TH2F by default,
  // so dynamic_cast<TH2D*> after >>h(...) always failed).
  // Axis: analyzer TPC_BINS_SIGNED_P / TPC_BINS_DE range; 150 bins for PDF.
  gROOT->cd();
  const TString hname = Form("h_pid_dedx_%s_%d", htag, pidIdx);
  if (auto* old = gROOT->FindObject(hname))
    delete old;

  auto* h = new TH2D(hname,
                     Form("%s;q#timesp [GeV/c];dE/dx (a.u.)", pageTitle),
                     150, -1.5, 1.5, 150, 0., 800.);
  const char* cutUse = (cut && cut[0] != '\0') ? cut : "";
  if (allTracks) {
    // All TPC track vector elements; no event selection cut.
    t->Project(hname, "dEdx:charge*mom0", cutUse);
  } else {
    t->Project(hname, "dEdx[i_kcand]:charge[i_kcand]*mom0[i_kcand]", cutUse);
    t->Draw(Form("dEdx[i_p]:charge[i_p]*mom0[i_p]>>+%s", hname.Data()),
            cutUse, "goff");
  }

  c.cd();
  if (h->GetEntries() <= 0. || h->Integral() <= 0.) {
    auto* tx = new TLatex(0.15, 0.5, "PID hist fill failed");
    tx->Draw();
    writer.Print(c);
    return;
  }

  if (h->GetMaximum() > 0.) {
    h->SetMinimum(0.);
    h->SetMaximum(h->GetMaximum() * 1.05);
  }
  h->Draw("COL");
  gPad->SetLogz(1);

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.028);
  if (allTracks) {
    tx->DrawLatex(0.14, 0.92,
                  Form("entries = %.0f  integral = %.0f  (all tracks)",
                       h->GetEntries(), h->Integral()));
    tx->DrawLatex(0.14, 0.87, "(no cut)");
  } else {
    tx->DrawLatex(0.14, 0.92,
                  Form("entries = %.0f  integral = %.0f  (K cand + proton)",
                       h->GetEntries(), h->Integral()));
    tx->DrawLatex(0.14, 0.87, (cutUse[0] != '\0') ? cutUse : "(no cut)");
    tx->SetTextSize(0.024);
    tx->DrawLatex(0.14, 0.82, "left: K cand (q<0)   right: proton (q>0)");
  }

  writer.Print(c);
  gPad->SetLogz(0);
}

// qp window where good-cut selected K concentrate (~mean±~2 rms of good sample).
constexpr Double_t kKaonQpLo = -0.85;
constexpr Double_t kKaonQpHi = -0.30;

void
DrawKaonRegionDedx(TTree* t, TCanvas& c, PdfWriter& writer)
{
  // Selected good K peak near qp≈-0.55; window covers that band on the K side.
  const TString cutQpAll =
    Form("charge*mom0>%g && charge*mom0<%g", kKaonQpLo, kKaonQpHi);
  const TString cutQpSel =
    Form("i_kcand>=0 && charge[i_kcand]*mom0[i_kcand]>%g"
         " && charge[i_kcand]*mom0[i_kcand]<%g",
         kKaonQpLo, kKaonQpHi);
  const TString cutQpGood =
    Form("%s && %s", cutQpSel.Data(), kCutGood);

  auto* hAll = Fill1D(t, "dEdx", cutQpAll, "h_dedx_kwin_all",
                      "dE/dx (a.u.)", 100, 0., 200.);
  auto* hSel = Fill1D(t, "dEdx[i_kcand]", cutQpSel, "h_dedx_kwin_sel",
                      "dE/dx (a.u.)", 100, 0., 200.);
  auto* hGood = Fill1D(t, "dEdx[i_kcand]", cutQpGood, "h_dedx_kwin_good",
                       "dE/dx (a.u.)", 100, 0., 200.);

  const TString winLabel =
    Form("%.2f < q#timesp < %.2f [GeV/c]", kKaonQpLo, kKaonQpHi);

  auto drawOne = [&](TH1D* h, Color_t col, const char* title, const char* note) {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    h->SetLineColor(col);
    h->SetLineWidth(2);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("dE/dx (a.u.)");
    h->GetYaxis()->SetTitle("Counts");
    if (h->GetMaximum() > 0.)
      h->SetMaximum(h->GetMaximum() * 1.25);
    h->Draw("hist");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.028);
    tx->DrawLatex(0.14, 0.92, Form("entries = %.0f", h->GetEntries()));
    tx->DrawLatex(0.14, 0.87, winLabel);
    tx->SetTextSize(0.024);
    tx->DrawLatex(0.14, 0.82, note);
    writer.Print(c);
  };

  drawOne(hAll, kBlack,
          "dE/dx in K qp window (all tracks)",
          "all tracks in window, no event cut");
  drawOne(hSel, kBlue + 1,
          "dE/dx in K qp window (selected K cand)",
          "dEdx[i_kcand], i_kcand>=0");
  drawOne(hGood, kOrange + 2,
          "dE/dx in K qp window (good-cut K)",
          kCutGood);

  // Overlay (unit-normalized) — stats differ a lot; compare shapes.
  c.Clear();
  c.cd();
  StyleCanvas(c);
  auto* nAll = static_cast<TH1D*>(hAll->Clone("h_dedx_kwin_all_n"));
  auto* nSel = static_cast<TH1D*>(hSel->Clone("h_dedx_kwin_sel_n"));
  auto* nGood = static_cast<TH1D*>(hGood->Clone("h_dedx_kwin_good_n"));
  nAll->SetDirectory(nullptr);
  nSel->SetDirectory(nullptr);
  nGood->SetDirectory(nullptr);
  const Double_t iAll = nAll->Integral();
  const Double_t iSel = nSel->Integral();
  const Double_t iGood = nGood->Integral();
  if (iAll > 0.)
    nAll->Scale(1. / iAll);
  if (iSel > 0.)
    nSel->Scale(1. / iSel);
  if (iGood > 0.)
    nGood->Scale(1. / iGood);

  nAll->SetLineColor(kBlack);
  nAll->SetLineWidth(2);
  nSel->SetLineColor(kBlue + 1);
  nSel->SetLineWidth(2);
  nGood->SetLineColor(kOrange + 2);
  nGood->SetLineWidth(2);

  const Double_t ymax = HistMax({nAll, nSel, nGood});
  nAll->SetMaximum(ymax > 0. ? ymax * 1.35 : 1.);
  nAll->SetTitle("dE/dx in K qp window (unit-normalized overlay)");
  nAll->GetXaxis()->SetTitle("dE/dx (a.u.)");
  nAll->GetYaxis()->SetTitle("Normalized counts");
  nAll->Draw("hist");
  nSel->Draw("hist same");
  nGood->Draw("hist same");

  auto* leg = new TLegend(0.48, 0.62, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(nAll, Form("all tracks (N=%.0f)", hAll->GetEntries()), "l");
  leg->AddEntry(nSel, Form("selected K (N=%.0f)", hSel->GetEntries()), "l");
  leg->AddEntry(nGood, Form("good-cut K (N=%.0f)", hGood->GetEntries()), "l");
  leg->Draw();

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.026);
  tx->DrawLatex(0.14, 0.92, winLabel);
  writer.Print(c);
}

void
DrawPidNsigmaOverlay(TTree* t, TCanvas& c, PdfWriter& writer,
                     const char* cutNone, const char* cutGood)
{
  c.Clear();
  c.cd();
  StyleCanvas(c);

  auto* hp0 = Fill1D(t, "nsigma_p_sel", cutNone, "h_nsp_none", "n#sigma", 100, -10., 10.);
  auto* hk0 = Fill1D(t, "nsigma_k_sel", cutNone, "h_nsk_none", "n#sigma", 100, -10., 10.);
  auto* hp1 = Fill1D(t, "nsigma_p_sel", cutGood, "h_nsp_good", "n#sigma", 100, -10., 10.);
  auto* hk1 = Fill1D(t, "nsigma_k_sel", cutGood, "h_nsk_good", "n#sigma", 100, -10., 10.);

  hp0->SetLineColor(kRed);       hp0->SetLineStyle(2); hp0->SetLineWidth(2);
  hk0->SetLineColor(kBlue);      hk0->SetLineStyle(2); hk0->SetLineWidth(2);
  hp1->SetLineColor(kRed + 1);   hp1->SetLineWidth(2);
  hk1->SetLineColor(kBlue + 1);  hk1->SetLineWidth(2);

  const Double_t ymax = HistMax({hp0, hk0, hp1, hk1});
  hp0->SetMaximum(ymax > 0. ? ymax * 1.30 : 1.);
  hp0->SetTitle("n#sigma PID: no cut vs good cut");
  hp0->GetXaxis()->SetTitle("n#sigma");
  hp0->GetYaxis()->SetTitle("Counts");
  hp0->Draw("hist");
  hk0->Draw("hist same");
  hp1->Draw("hist same");
  hk1->Draw("hist same");

  auto* leg = new TLegend(0.48, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.026);
  leg->AddEntry(hp0, Form("p n#sigma (mmass>0) %.0f", hp0->GetEntries()), "l");
  leg->AddEntry(hk0, Form("K n#sigma (mmass>0) %.0f", hk0->GetEntries()), "l");
  leg->AddEntry(hp1, Form("p n#sigma (good) %.0f", hp1->GetEntries()), "l");
  leg->AddEntry(hk1, Form("K n#sigma (good) %.0f", hk1->GetEntries()), "l");
  leg->Draw();

  writer.Print(c);
}

} // namespace

int
main(int argc, char** argv)
{
  Int_t run = -1;
  TString inPath;
  TString outPdf;

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-r" && i + 1 < argc) {
      run = TString(argv[++i]).Atoi();
    } else if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.Length() > 0 && a[0] == '-') {
      std::cerr << "Unknown option: " << a << std::endl;
      usage(argv[0]);
      return 1;
    } else if (inPath.IsNull()) {
      inPath = a;
    } else {
      std::cerr << "Extra argument: " << a << std::endl;
      usage(argv[0]);
      return 1;
    }
  }

  if (inPath.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  Int_t runOut = run;
  if (runOut < 0)
    runOut = 0;

  gROOT->SetBatch(kTRUE);
  PreparePalette();
  gStyle->SetOptStat(1110);

  TFile f(inPath, "READ");
  if (!f.IsOpen() || f.IsZombie()) {
    std::cerr << "Error: cannot open " << inPath << std::endl;
    return 1;
  }

  auto* t = dynamic_cast<TTree*>(f.Get("kpsc"));
  if (!t) {
    std::cerr << "Error: tree kpsc not found in " << inPath << std::endl;
    return 1;
  }

  if (run < 0) {
    const Int_t runTree = ParseRunFromTree(t);
    if (runTree > 0)
      runOut = runTree;
    else
      runOut = ParseRunFromPath(inPath);
  }
  if (runOut < 0)
    runOut = 0;

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
  if (outPdf.IsNull())
    outPdf = Form("%s/kpsc_purity_compare_run%05d.pdf", imgDir.Data(), runOut);

  TCanvas c("c_kpsc_purity", "kpsc purity compare", 900, 700);
  PdfWriter writer(outPdf);

  const std::vector<CutSet> mmassCuts = {
    kSummaryCuts[0], kSummaryCuts[1], kSummaryCuts[4], kSummaryCuts[7], kSummaryCuts[8]};

  DrawSummaryTable(t, c, writer);
  DrawNtMmassOverlay(t, c, writer);
  DrawOverlay1D(t, c, writer, "mmass", "missing mass [GeV]", 120, 0., 1.2,
                "Missing mass vs cut stage", mmassCuts);
  DrawOverlay1D(t, c, writer, "effective_ntTpc", "effective_ntTpc", 11, -0.5, 10.5,
                "effective_ntTpc vs cut stage", mmassCuts);

  DrawPidDedxVsQp(t, c, writer, "",
                  "PID dE/dx vs qp (all tracks, no cut)", "all", true);
  DrawPidDedxVsQp(t, c, writer, kCutGood,
                  "PID dE/dx vs qp (good: Nt2+AM+MM0.1)", "good", false);
  DrawKaonRegionDedx(t, c, writer);
  DrawPidNsigmaOverlay(t, c, writer, kCutBase, kCutGood);
  DrawNsigma2D(t, c, writer, kCutGood, "good");

  const Var2D mmass2dVars[] = {
    {"effective_ntTpc", "effective_ntTpc", 11, -0.5, 10.5},
    {"dp_miss", "dp_{miss} [GeV/c]", 100, 0., 1.},
    {"angle_miss_k", "#angle(p_{miss},p_{K}) [rad]", 100, 0., TMath::Pi()},
    {"nsigma_k_sel", "n#sigma_{K}^{sel}", 100, -10., 10.},
    {"nsigma_pi_sel", "n#sigma_{#pi}^{sel}", 100, -10., 10.},
    {"nsigma_p_sel", "n#sigma_{p}^{sel}", 100, -10., 10.},
    {"close_dist", "close_{dist} [mm]", 100, 0., 50.},
    {"dphi", "d#phi [rad]", 64, 0., TMath::Pi()},
    {"n_proton", "n_{proton}", 11, -0.5, 10.5},
    {"n_qminus", "n_{q-}", 11, -0.5, 10.5},
  };
  for (const auto& v : mmass2dVars)
    DrawMmass2D(t, c, writer, v, kCutBase, Form("mmass vs %s (kinematics OK)", v.var));

  const Var1D peakVars[] = {
    {"nsigma_k_sel", "n#sigma_{K}^{sel}", 100, -10., 10.},
    {"nsigma_pi_sel", "n#sigma_{#pi}^{sel}", 100, -10., 10.},
    {"nsigma_p_sel", "n#sigma_{p}^{sel}", 100, -10., 10.},
    {"angle_miss_k", "#angle(p_{miss},p_{K}) [rad]", 100, 0., TMath::Pi()},
    {"dp_miss", "dp_{miss} [GeV/c]", 100, 0., 1.},
    {"close_dist", "close_{dist} [mm]", 100, 0., 50.},
    {"dphi", "d#phi [rad]", 64, 0., TMath::Pi()},
  };
  for (const auto& v : peakVars) {
    DrawPair1D(t, c, writer, v.var, v.title,
               kCutC2, "C2: Nt2+MMKW peak", kRed + 1,
               kCutSidebandNt2, "Nt2 sideband", kGray + 2,
               v.nb, v.xlo, v.xhi,
               Form("Peak vs sideband (Nt2): %s", v.var));
    DrawPair1D(t, c, writer, v.var, v.title,
               kCutC4, "C4: Nt<=4+MMKW peak", kGreen + 2,
               kCutSidebandNt4, "Nt<=4 sideband", kGray + 2,
               v.nb, v.xlo, v.xhi,
               Form("Peak vs sideband (Nt<=4): %s", v.var));
  }

  DrawNsigma2D(t, c, writer, kCutC2, "Nt2+MMKW");
  DrawNsigma2D(t, c, writer, kCutC4, "Nt<=4+MMKW");

  DrawPair1D(t, c, writer, "dp_miss", "dp_{miss} [GeV/c]",
             kCutC2, "C2: Nt2+MMKW", kRed + 1,
             kCutWorking, "D: working", kGreen + 2,
             100, 0., 1., "dp_{miss}: Nt2+MMKW vs working cut");
  DrawPair1D(t, c, writer, "angle_miss_k", "#angle(p_{miss},p_{K}) [rad]",
             kCutC2, "C2: Nt2+MMKW", kRed + 1,
             kCutGood, "G: good", kOrange + 2,
             100, 0., TMath::Pi(), "angle_miss_k: Nt2+MMKW vs good cut");
  DrawPair1D(t, c, writer, "cos_theta_k_cm", "cos#theta_{K}^{CM}",
             kCutGood, "G: good", kOrange + 2,
             kCutWorking, "D: working", kGreen + 2,
             100, -1., 1., "cos#theta_{K}^{CM}: good vs working");

  writer.Close(c);
  std::cout << "Wrote " << outPdf << " (" << writer.Pages() << " pages)" << std::endl;
  return 0;
}
