// -*- C++ -*-
// Kp scattering purity review: multipage PDF from DstKpScattering output.
//
// Usage:
//   kpsc_purity_compare_pdf <kpsc.root> [-o <out.pdf>] [-r <run>] [--forslide]
//
// Input:
//   TTree "kpsc" from bin/DstKpScattering output.
//
// Output (default):
//   OUTPUT_DIR/img/runXXXXX/kpsc_purity_compare_runXXXXX.pdf
//   --forslide: tip-cut figures only (mmass 2D / tip mmass / PID), slide styling.

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
#include <TPad.h>
#include <TProfile.h>
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
constexpr Double_t kAngleMissKTip = 0.1;
constexpr Double_t kCloseDistTip = 5.0;
constexpr Double_t kPidYMaxSlide = 500.;
constexpr Double_t kSlideMmassHi = 0.8;
constexpr Double_t kSlideCloseDistHi = 30.;
constexpr Double_t kSlideDedxHi = 500.;
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
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.3 && abs(mmass-0.494)<0.1";
// Tip cut (no mmass window): Nt==2 + tight angle + close_dist.
constexpr const char* kCutTip =
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.1 && close_dist<5";
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
  {"C2", "C2: Nt2+MMass_Peak_Wide", kCutC2, kRed + 1},
  {"C3", "C3: Nt<=3+MMass_Peak_Wide", kCutC3, kOrange + 7},
  {"C4", "C4: Nt<=4+MMass_Peak_Wide", kCutC4, kMagenta + 1},
  {"G", "G: good (AM+MM0.1)", kCutGood, kOrange + 2},
  {"T", "T: tip (Nt2+AM0.1+cd5)", kCutTip, kViolet + 1},
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
  std::cerr << "Usage: " << argv0
            << " <kpsc.root> [-o <out.pdf>] [-r <run>] [--forslide]\n"
            << "  Reads TTree kpsc.\n"
            << "  Default PDF: OUTPUT_DIR/img/runXXXXX/kpsc_purity_compare_runXXXXX.pdf\n"
            << "  --forslide: tip-cut slide pages only (no legends/titles).\n";
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
  auto* h = dynamic_cast<TH1*>(gPad->GetPrimitive(name));
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
  tx->DrawLatex(0.08, y, Form("Nt2+MMass_Peak_Wide + pi veto : %.0f", CountEvents(t, kCutPiVeto)));

  tx->SetTextColor(kBlack);
  tx->SetTextSize(0.022);
  tx->DrawLatex(0.08, 0.16, Form("MMass_Peak_Wide: |mmass - %.3f| < %.2f GeV (diag) / %.2f (good)",
                                 kMK, kMmassWindowDiag, kMmassWindowTight));
  tx->DrawLatex(0.08, 0.12, Form("good cut: Nt==2 && diff_angle<%.1f && |mmass-m_K|<%.2f",
                                 kAngleMissKGood, kMmassWindowTight));
  tx->DrawLatex(0.08, 0.08, Form("tip cut: Nt==2 && diff_angle<%.2f && close_dist<%.0f (no mmass window)",
                                 kAngleMissKTip, kCloseDistTip));
  tx->DrawLatex(0.08, 0.04, Form("sideband: %.2f < mmass < %.2f GeV; #pi-like K: |n#sigma_{#pi}|<%.0f && n#sigma_{K}<%.0f",
                                 kSidebandLo, kSidebandHi, kNsigmaPionWindow, kNsigmaKaonReject));

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
            const Var2D& v, const char* baseCut, const char* pageTitle,
            Int_t nby = 120, Double_t ylo = 0., Double_t yhi = 1.2,
            Bool_t drawBands = true, Bool_t markGoodAngle = false)
{
  static Int_t mm2dIdx = 0;
  ++mm2dIdx;

  c.Clear();
  c.cd();
  StyleCanvas(c, true);

  const TString hname = Form("h_mmass_vs_%s_%d", v.var, mm2dIdx);
  const TString title = Form("%s; %s;missing mass [GeV]", pageTitle, v.xtitle);
  DrawTree2D(t, c, "mmass", v.var, baseCut, hname, title,
             v.nbx, v.xlo, v.xhi, nby, ylo, yhi, "COL");
  auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(hname));
  if (h && drawBands)
    DrawMmassBands(h);
  if (h && markGoodAngle) {
    auto* lAng = new TLine(kAngleMissKGood, ylo, kAngleMissKGood, yhi);
    lAng->SetLineColor(kOrange + 2);
    lAng->SetLineStyle(2);
    lAng->SetLineWidth(2);
    lAng->Draw("same");
  }

  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.028);
  if (h)
    tx->DrawLatex(0.14, 0.92, Form("entries = %.0f", h->GetEntries()));
  tx->DrawLatex(0.14, 0.87, baseCut);

  writer.Print(c);
}

// mmass–diff_angle V: zoom / slices / profile / |p| correlation (Nt2 focus).
void
DrawDiffAngleVStudy(TTree* t, TCanvas& c, PdfWriter& writer)
{
  constexpr Double_t kZoomLo = 0.30;
  constexpr Double_t kZoomHi = 0.70;
  constexpr Int_t kNbAngle = 64;
  constexpr Double_t kAngHi = TMath::Pi();
  const Var2D vAng{"diff_angle", "#angle(p_{miss},p_{K}) [rad]",
                   kNbAngle, 0., kAngHi};

  // --- summary numbers ---
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    const Double_t nNt2 = CountEvents(t, kCutB2);
    const Double_t nTip = CountEvents(
      t, "mmass>0 && effective_ntTpc==2 && diff_angle<0.3");
    const Double_t nTipPeak = CountEvents(t, kCutGood);
    const Double_t nWing = CountEvents(
      t, "mmass>0 && effective_ntTpc==2 && diff_angle>=0.3");
    const Double_t nWingPeak = CountEvents(
      t,
      "mmass>0 && effective_ntTpc==2 && diff_angle>=0.3"
      " && abs(mmass-0.494)<0.1");

    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.032);
    tx->DrawLatex(0.08, 0.92, "diff_angle V study (kinematics consistency)");
    tx->SetTextSize(0.026);
    tx->DrawLatex(0.08, 0.84,
                  "mmass = miss from beam+target-p; diff_angle = #angle(p_{miss},p_{K})");
    tx->DrawLatex(0.08, 0.78,
                  "V near m_{K}: correlated proton resolution (not global mom flip)");
    tx->SetTextSize(0.028);
    Double_t y = 0.68;
    tx->DrawLatex(0.08, y, Form("Nt2 (B2)                 : %.0f", nNt2));
    y -= 0.05;
    tx->DrawLatex(0.08, y, Form("Nt2 && diff_angle<0.3    : %.0f  (%.1f%% of Nt2)",
                                nTip, nNt2 > 0. ? 100. * nTip / nNt2 : 0.));
    y -= 0.05;
    tx->DrawLatex(0.08, y, Form("good (tip + |mmass-m_K|<0.1): %.0f  (%.1f%% of Nt2)",
                                nTipPeak, nNt2 > 0. ? 100. * nTipPeak / nNt2 : 0.));
    y -= 0.05;
    tx->DrawLatex(0.08, y, Form("Nt2 && diff_angle>=0.3   : %.0f", nWing));
    y -= 0.05;
    tx->DrawLatex(0.08, y,
                  Form("wing && |mmass-m_K|<0.1  : %.0f  (wrong-K-like band?)",
                       nWingPeak));
    y -= 0.08;
    tx->SetTextSize(0.024);
    tx->SetTextColor(kGray + 2);
    tx->DrawLatex(0.08, y, "Next pages: zoom 2D (Base/Nt2/sideband), profile, angle slices, |p| corr.");
    writer.Print(c);
  }

  // --- zoomed 2D ---
  DrawMmass2D(t, c, writer, vAng, kCutBase,
              "ZOOM: mmass vs diff_angle (Base)", 80, kZoomLo, kZoomHi, true, true);
  DrawMmass2D(t, c, writer, vAng, kCutB2,
              "ZOOM: mmass vs diff_angle (Nt2)", 80, kZoomLo, kZoomHi, true, true);
  DrawMmass2D(t, c, writer, vAng, kCutSidebandNt2,
              "ZOOM: mmass vs diff_angle (Nt2 sideband)", 80, kZoomLo, kZoomHi,
              true, true);

  // --- profile <mmass> vs diff_angle ---
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    static Int_t profIdx = 0;
    ++profIdx;
    const TString pname = Form("h_mmass_prof_dang_%d", profIdx);
    if (auto* old = gDirectory->Get(pname))
      delete old;
    auto* prof = new TProfile(pname,
                              "Nt2: <mmass> vs diff_angle; #angle(p_{miss},p_{K}) [rad];"
                              "<missing mass> [GeV]",
                              kNbAngle, 0., kAngHi, kZoomLo, kZoomHi, "s");
    t->Project(pname, "mmass:diff_angle", kCutB2, "prof");
    auto* pc = static_cast<TProfile*>(prof->Clone(Form("%s_c", pname.Data())));
    pc->SetDirectory(nullptr);
    delete prof;
    pc->SetLineColor(kBlue + 1);
    pc->SetMarkerColor(kBlue + 1);
    pc->SetMarkerStyle(20);
    pc->SetMarkerSize(0.7);
    pc->GetYaxis()->SetRangeUser(kZoomLo, kZoomHi);
    pc->Draw("E1");
    auto* lMk = new TLine(0., kMK, kAngHi, kMK);
    lMk->SetLineColor(kRed + 1);
    lMk->SetLineWidth(2);
    lMk->Draw("same");
    auto* lAng = new TLine(kAngleMissKGood, kZoomLo, kAngleMissKGood, kZoomHi);
    lAng->SetLineColor(kOrange + 2);
    lAng->SetLineStyle(2);
    lAng->SetLineWidth(2);
    lAng->Draw("same");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.028);
    tx->DrawLatex(0.14, 0.92, Form("entries = %.0f  (%s)", pc->GetEntries(), kCutB2));
    tx->DrawLatex(0.14, 0.87, "red: m_{K}; orange dashed: good cut angle");
    writer.Print(c);
  }

  // --- mmass 1D in angle slices (Nt2) ---
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    struct Slice {
      const char* name;
      const char* cut;
      const char* label;
      Color_t color;
    };
    const Slice slices[] = {
      {"tip",
       "mmass>0 && effective_ntTpc==2 && diff_angle<0.1",
       "diff_angle<0.1", kBlue + 1},
      {"mid",
       "mmass>0 && effective_ntTpc==2 && diff_angle>=0.1 && diff_angle<0.3",
       "0.1#leq diff_angle<0.3", kOrange + 7},
      {"wing",
       "mmass>0 && effective_ntTpc==2 && diff_angle>=0.3",
       "diff_angle#geq0.3", kGray + 2},
    };
    std::vector<TH1D*> hs;
    for (const auto& s : slices) {
      auto* h = Fill1D(t, "mmass", s.cut, Form("h_mmass_dang_%s", s.name),
                       "missing mass [GeV]", 100, kZoomLo, kZoomHi);
      h->SetLineColor(s.color);
      h->SetLineWidth(2);
      hs.push_back(h);
    }
    const Double_t ymax = HistMax(hs);
    hs.front()->SetMaximum(ymax > 0. ? ymax * 1.25 : 1.);
    hs.front()->SetTitle("Nt2: mmass in diff_angle slices (zoom)");
    hs.front()->GetXaxis()->SetTitle("missing mass [GeV]");
    hs.front()->Draw("hist");
    for (std::size_t i = 1; i < hs.size(); ++i)
      hs[i]->Draw("hist same");
    auto* leg = new TLegend(0.50, 0.62, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (std::size_t i = 0; i < hs.size(); ++i)
      leg->AddEntry(hs[i], Form("%s (%.0f)", slices[i].label, hs[i]->GetEntries()), "l");
    leg->Draw();
    writer.Print(c);
  }

  // --- |p| correlation (Nt2, zoomed mmass) ---
  {
    const Var2D vSigned{"diff_p_miss_k_signed",
                        "|p_{miss}|-|p_{K}| [GeV/c]", 100, -0.5, 0.5};
    DrawMmass2D(t, c, writer, vSigned, kCutB2,
                "ZOOM: mmass vs Diff_P_MissK_Signed (Nt2)",
                80, kZoomLo, kZoomHi, true, false);
  }
  {
    c.Clear();
    c.cd();
    StyleCanvas(c, true);
    static Int_t corrIdx = 0;
    ++corrIdx;
    const TString hname = Form("h_dang_vs_dpsigned_%d", corrIdx);
    DrawTree2D(t, c, "diff_angle", "diff_p_miss_k_signed", kCutB2, hname,
               "Nt2: diff_angle vs Diff_P_MissK_Signed;"
               "|p_{miss}|-|p_{K}| [GeV/c];#angle(p_{miss},p_{K}) [rad]",
               100, -0.5, 0.5, kNbAngle, 0., kAngHi, "COL");
    auto* lAng = new TLine(-0.5, kAngleMissKGood, 0.5, kAngleMissKGood);
    lAng->SetLineColor(kOrange + 2);
    lAng->SetLineStyle(2);
    lAng->SetLineWidth(2);
    lAng->Draw("same");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.028);
    if (auto* h = dynamic_cast<TH2D*>(gPad->GetPrimitive(hname)))
      tx->DrawLatex(0.14, 0.92, Form("entries = %.0f", h->GetEntries()));
    tx->DrawLatex(0.14, 0.87, kCutB2);
    writer.Print(c);
  }
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

void
StyleSlideAxes(TH1* h, const char* xtitle, const char* ytitle)
{
  if (!h)
    return;
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->CenterTitle(false);
  h->GetYaxis()->CenterTitle(false);
  h->GetXaxis()->SetTitleSize(0.050);
  h->GetYaxis()->SetTitleSize(0.050);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.25);
}

// Plot ~70% (left, fixed aspect), cut note on the right.
struct SlideLayout
{
  TPad* plot = nullptr;
  TPad* note = nullptr;
};

SlideLayout
BeginSlidePage(TCanvas& c, Bool_t zPalette)
{
  c.Clear();
  c.cd();
  // Left ~70%: plot. Extra top margin so any title band stays above axes
  // and can be cropped out of a screenshot of the figure body.
  auto* plot = new TPad("slide_plot", "", 0.01, 0.04, 0.70, 0.92);
  plot->SetFillStyle(0);
  plot->SetLeftMargin(0.14);
  plot->SetRightMargin(zPalette ? 0.15 : 0.05);
  plot->SetBottomMargin(0.14);
  plot->SetTopMargin(0.12);
  plot->Draw();

  auto* note = new TPad("slide_note", "", 0.71, 0.04, 0.99, 0.92);
  note->SetFillStyle(0);
  note->SetLeftMargin(0.05);
  note->SetRightMargin(0.05);
  note->Draw();
  return {plot, note};
}

void
DrawSlideCutNote(TPad* note, const std::vector<TString>& lines)
{
  if (!note)
    return;
  note->cd();
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextAlign(13);
  tx->SetTextFont(42);
  tx->SetTextSize(0.070);
  Double_t y = 0.92;
  for (const auto& line : lines) {
    if (line.IsNull()) {
      y -= 0.04;
      continue;
    }
    tx->DrawLatex(0.04, y, line);
    y -= 0.085;
  }
}

const std::vector<TString> kSlideNoteBase = {
  "cut:",
  "M_{miss} > 0",
};

const std::vector<TString> kSlideNoteTip = {
  "cut:",
  "M_{miss} > 0",
  "N_{TPC}^{eff} = 2",
  "diff_angle < 0.1",
  "d_{close} < 5 mm",
};

void
DrawSlideMmass2D(TTree* t, TCanvas& c, PdfWriter& writer, const Var2D& v)
{
  static Int_t slide2dIdx = 0;
  ++slide2dIdx;

  auto pads = BeginSlidePage(c, true);
  pads.plot->cd();

  const TString hname = Form("h_slide_mmass_vs_%s_%d", v.var, slide2dIdx);
  const TString title = Form(";%s;M_{miss} [GeV/c^{2}]", v.xtitle);
  DrawTree2D(t, c, "mmass", v.var, kCutBase, hname, title,
             v.nbx, v.xlo, v.xhi, 100, 0., kSlideMmassHi, "COLZ");
  if (auto* h = dynamic_cast<TH1*>(gPad->GetPrimitive(hname))) {
    h->SetTitle("");
    StyleSlideAxes(h, v.xtitle, "M_{miss} [GeV/c^{2}]");
  }

  DrawSlideCutNote(pads.note, kSlideNoteBase);
  c.cd();
  writer.Print(c);
}

void
DrawSlide1D(TTree* t, TCanvas& c, PdfWriter& writer,
            const char* var, const char* xtitle, const char* cut,
            Int_t nb, Double_t xlo, Double_t xhi, Color_t color,
            const std::vector<TString>& noteLines)
{
  static Int_t slide1dIdx = 0;
  ++slide1dIdx;

  auto pads = BeginSlidePage(c, false);
  pads.plot->cd();

  const TString hname = Form("h_slide_%s_%d", var, slide1dIdx);
  auto* h = Fill1D(t, var, cut, hname, "", nb, xlo, xhi);
  h->SetLineColor(color);
  h->SetLineWidth(2);
  StyleSlideAxes(h, xtitle, "Counts");
  if (h->GetMaximum() > 0.)
    h->SetMaximum(h->GetMaximum() * 1.15);
  h->Draw("hist");

  DrawSlideCutNote(pads.note, noteLines);
  c.cd();
  writer.Print(c);
}

void
DrawSlideMmassOverlay(TTree* t, TCanvas& c, PdfWriter& writer)
{
  static Int_t ovIdx = 0;
  ++ovIdx;

  auto pads = BeginSlidePage(c, false);
  pads.plot->cd();

  auto* hAll = Fill1D(t, "mmass", kCutBase,
                      Form("h_slide_mmass_all_%d", ovIdx), "", 100, 0., kSlideMmassHi);
  auto* hTip = Fill1D(t, "mmass", kCutTip,
                      Form("h_slide_mmass_tip_%d", ovIdx), "", 100, 0., kSlideMmassHi);
  hAll->SetLineColor(kGray + 2);
  hAll->SetLineWidth(2);
  hTip->SetLineColor(kRed + 1);
  hTip->SetLineWidth(2);
  StyleSlideAxes(hAll, "M_{miss} [GeV/c^{2}]", "Counts");
  const Double_t ymax = std::max(hAll->GetMaximum(), hTip->GetMaximum());
  hAll->SetMaximum(ymax > 0. ? ymax * 1.20 : 1.);
  hAll->Draw("hist");
  hTip->Draw("hist same");

  const std::vector<TString> note = {
    "gray: M_{miss} > 0",
    "",
    "red (tip):",
    "M_{miss} > 0",
    "N_{TPC}^{eff} = 2",
    "diff_angle < 0.1",
    "d_{close} < 5 mm",
  };
  DrawSlideCutNote(pads.note, note);
  c.cd();
  writer.Print(c);
}

void
DrawSlideDedxVsQp(TTree* t, TCanvas& c, PdfWriter& writer,
                  const char* cut, const char* htag, Bool_t allTracks,
                  const std::vector<TString>& noteLines)
{
  static Int_t slidePidIdx = 0;
  ++slidePidIdx;

  auto pads = BeginSlidePage(c, true);
  pads.plot->cd();
  PreparePalette();

  gROOT->cd();
  const TString hname = Form("h_slide_dedx_%s_%d", htag, slidePidIdx);
  if (auto* old = gROOT->FindObject(hname))
    delete old;

  auto* h = new TH2D(hname, ";q#timesp [GeV/c];dE/dx [a.u.]",
                     150, -1.5, 1.5, 150, 0., kSlideDedxHi);
  const char* cutUse = (cut && cut[0] != '\0') ? cut : "";
  if (allTracks) {
    t->Project(hname, "dEdx:charge*mom0", cutUse);
  } else {
    t->Project(hname, "dEdx[i_kcand]:charge[i_kcand]*mom0[i_kcand]", cutUse);
    t->Draw(Form("dEdx[i_p]:charge[i_p]*mom0[i_p]>>+%s", hname.Data()),
            cutUse, "goff");
  }

  pads.plot->cd();
  StyleSlideAxes(h, "q#timesp [GeV/c]", "dE/dx [a.u.]");
  if (h->GetEntries() > 0. && h->Integral() > 0. && h->GetMaximum() > 0.) {
    h->SetMinimum(0.);
    h->SetMaximum(h->GetMaximum() * 1.05);
  }
  h->Draw("COLZ");
  gPad->SetLogz(1);

  DrawSlideCutNote(pads.note, noteLines);
  c.cd();
  writer.Print(c);
  pads.plot->cd();
  gPad->SetLogz(0);
}

void
DrawForSlide(TTree* t, TCanvas& c, PdfWriter& writer)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c.SetCanvasSize(1000, 700);

  const Var2D slide2d[] = {
    {"effective_ntTpc", "N_{TPC}^{eff}", 11, -0.5, 10.5},
    {"diff_angle", "diff_angle", 100, 0., TMath::Pi()},
    {"close_dist", "d_{close} [mm]", 100, 0., kSlideCloseDistHi},
  };
  for (const auto& v : slide2d)
    DrawSlideMmass2D(t, c, writer, v);

  DrawSlide1D(t, c, writer, "mmass", "M_{miss} [GeV/c^{2}]", kCutBase,
              100, 0., kSlideMmassHi, kGray + 2, kSlideNoteBase);
  DrawSlide1D(t, c, writer, "mmass", "M_{miss} [GeV/c^{2}]", kCutTip,
              100, 0., kSlideMmassHi, kRed + 1, kSlideNoteTip);
  DrawSlideMmassOverlay(t, c, writer);

  const std::vector<TString> noteAllTracks = {
    "cut: (none)",
    "all TPC tracks",
  };
  const std::vector<TString> noteTipTracks = {
    "cut (tip):",
    "M_{miss} > 0",
    "N_{TPC}^{eff} = 2",
    "diff_angle < 0.1",
    "d_{close} < 5 mm",
    "",
    "tracks: K cand + p",
  };
  DrawSlideDedxVsQp(t, c, writer, "", "all", true, noteAllTracks);
  DrawSlideDedxVsQp(t, c, writer, kCutTip, "tip", false, noteTipTracks);
}

void
DrawTipPidPages(TTree* t, TCanvas& c, PdfWriter& writer)
{
  const Int_t oldStat = gStyle->GetOptStat();
  gStyle->SetOptStat(0);
  const Var1D pidVars[] = {
    {"nsigma_p_sel", "n#sigma_{p}^{sel}", 100, -10., 10.},
    {"nsigma_pi_p_sel", "n#sigma_{#pi}^{p,sel}", 100, -10., 10.},
    {"nsigma_k_sel", "n#sigma_{K}^{sel}", 100, -10., 10.},
    {"nsigma_pi_sel", "n#sigma_{#pi}^{sel}", 100, -10., 10.},
  };
  for (const auto& v : pidVars) {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    const TString hname = Form("h_tip_%s", v.var);
    auto* h = Fill1D(t, v.var, kCutTip, hname, v.title, v.nb, v.xlo, v.xhi);
    h->SetLineColor(kViolet + 1);
    h->SetLineWidth(2);
    h->SetTitle(Form("tip cut: %s", v.title));
    h->GetXaxis()->SetTitle(v.title);
    h->GetYaxis()->SetTitle("Counts");
    h->SetMaximum(kPidYMaxSlide);
    h->Draw("hist");
    writer.Print(c);
  }
  gStyle->SetOptStat(oldStat);
}

} // namespace

int
main(int argc, char** argv)
{
  Int_t run = -1;
  TString inPath;
  TString outPdf;
  Bool_t forSlide = false;

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-r" && i + 1 < argc) {
      run = TString(argv[++i]).Atoi();
    } else if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "--forslide") {
      forSlide = true;
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
  gStyle->SetOptStat(forSlide ? 0 : 1110);

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
  if (outPdf.IsNull()) {
    if (forSlide)
      outPdf = Form("%s/kpsc_purity_slide_run%05d.pdf", imgDir.Data(), runOut);
    else
      outPdf = Form("%s/kpsc_purity_compare_run%05d.pdf", imgDir.Data(), runOut);
  }

  TCanvas c("c_kpsc_purity", "kpsc purity compare", 900, 700);
  PdfWriter writer(outPdf);

  if (forSlide) {
    DrawForSlide(t, c, writer);
    writer.Close(c);
    std::cout << "Wrote " << outPdf << " (" << writer.Pages() << " pages, forslide)"
              << std::endl;
    return 0;
  }

  // A, B2, C2, G, T, D
  const std::vector<CutSet> mmassCuts = {
    kSummaryCuts[0], kSummaryCuts[1], kSummaryCuts[4],
    kSummaryCuts[7], kSummaryCuts[8], kSummaryCuts[9]};

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
  DrawTipPidPages(t, c, writer);

  const Var2D mmass2dVars[] = {
    {"effective_ntTpc", "effective_ntTpc", 11, -0.5, 10.5},
    {"diff_p_miss_k", "dp_{miss} [GeV/c]", 100, 0., 1.},
    {"diff_angle", "#angle(p_{miss},p_{K}) [rad]", 100, 0., TMath::Pi()},
    {"nsigma_k_sel", "n#sigma_{K}^{sel}", 100, -10., 10.},
    {"nsigma_pi_sel", "n#sigma_{#pi}^{sel}", 100, -10., 10.},
    {"nsigma_p_sel", "n#sigma_{p}^{sel}", 100, -10., 10.},
    {"close_dist", "close_{dist} [mm]", 100, 0., 50.},
    {"delta_phi_scat", "d#phi [rad]", 64, 0., TMath::Pi()},
    {"n_proton", "n_{proton}", 11, -0.5, 10.5},
    {"n_qminus", "n_{q-}", 11, -0.5, 10.5},
  };
  for (const auto& v : mmass2dVars)
    DrawMmass2D(t, c, writer, v, kCutBase, Form("mmass vs %s (kinematics OK)", v.var));

  DrawDiffAngleVStudy(t, c, writer);

  const Var1D peakVars[] = {
    {"nsigma_k_sel", "n#sigma_{K}^{sel}", 100, -10., 10.},
    {"nsigma_pi_sel", "n#sigma_{#pi}^{sel}", 100, -10., 10.},
    {"nsigma_p_sel", "n#sigma_{p}^{sel}", 100, -10., 10.},
    {"diff_angle", "#angle(p_{miss},p_{K}) [rad]", 100, 0., TMath::Pi()},
    {"diff_p_miss_k", "dp_{miss} [GeV/c]", 100, 0., 1.},
    {"close_dist", "close_{dist} [mm]", 100, 0., 50.},
    {"delta_phi_scat", "d#phi [rad]", 64, 0., TMath::Pi()},
  };
  for (const auto& v : peakVars) {
    DrawPair1D(t, c, writer, v.var, v.title,
               kCutC2, "C2: Nt2+MMass_Peak_Wide peak", kRed + 1,
               kCutSidebandNt2, "Nt2 sideband", kGray + 2,
               v.nb, v.xlo, v.xhi,
               Form("Peak vs sideband (Nt2): %s", v.var));
    DrawPair1D(t, c, writer, v.var, v.title,
               kCutC4, "C4: Nt<=4+MMass_Peak_Wide peak", kGreen + 2,
               kCutSidebandNt4, "Nt<=4 sideband", kGray + 2,
               v.nb, v.xlo, v.xhi,
               Form("Peak vs sideband (Nt<=4): %s", v.var));
  }

  DrawNsigma2D(t, c, writer, kCutC2, "Nt2+MMass_Peak_Wide");
  DrawNsigma2D(t, c, writer, kCutC4, "Nt<=4+MMass_Peak_Wide");

  DrawPair1D(t, c, writer, "diff_p_miss_k", "dp_{miss} [GeV/c]",
             kCutC2, "C2: Nt2+MMass_Peak_Wide", kRed + 1,
             kCutWorking, "D: working", kGreen + 2,
             100, 0., 1., "dp_{miss}: Nt2+MMass_Peak_Wide vs working cut");
  DrawPair1D(t, c, writer, "diff_angle", "#angle(p_{miss},p_{K}) [rad]",
             kCutC2, "C2: Nt2+MMass_Peak_Wide", kRed + 1,
             kCutGood, "G: good", kOrange + 2,
             100, 0., TMath::Pi(), "diff_angle: Nt2+MMass_Peak_Wide vs good cut");
  DrawPair1D(t, c, writer, "cos_theta_cm", "cos#theta_{K}^{CM}",
             kCutGood, "G: good", kOrange + 2,
             kCutWorking, "D: working", kGreen + 2,
             100, -1., 1., "cos#theta_{K}^{CM}: good vs working");
  DrawPair1D(t, c, writer, "mmass", "missing mass [GeV]",
             kCutGood, "G: good", kOrange + 2,
             kCutTip, "T: tip", kViolet + 1,
             120, 0., 1.2, "mmass: good vs tip cut");

  writer.Close(c);
  std::cout << "Wrote " << outPdf << " (" << writer.Pages() << " pages)" << std::endl;
  return 0;
}
