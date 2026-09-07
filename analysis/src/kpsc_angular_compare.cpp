// -*- C++ -*-
// Kp angular distribution: acceptance-corrected shape vs existing dσ data.
//
// Usage:
//   kpsc_angular_compare <kpsc.root> [-o <out.pdf>] [-a <acceptance.root>] [-r <run>]
//
// Input:
//   TTree "kpsc" from DstKpScattering.
// Default acceptance:
//   ANALYZER_DIR/trigger_study_E72/results/root/acceptance_kp_flat_9999.root
// Output (default):
//   OUTPUT_DIR/img/runXXXXX/kpsc_angular_compare_runXXXXX.pdf

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
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

#include <Math/SpecFuncMathMore.h>

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr Double_t kAccMomStart = 600.0;   // MeV/c
constexpr Double_t kAccMomStep = 0.5;      // MeV/c
constexpr Int_t kAccMomN = 800;
constexpr Double_t kAccMin = 0.20;         // ignore / hide A <= 20%
constexpr Int_t kNCosBins = 50;
constexpr Int_t kNLegendre = 6;
// Sparse display: Legendre CSV is only A_ell(p), not dense cosθ points.
constexpr Double_t kRefCosStep = 0.20;
constexpr Double_t kMK = 0.494;
constexpr Double_t kMmassWindow = 0.10;
constexpr Double_t kMassKaon = 493.677;    // MeV
constexpr Double_t kMassProton = 938.27208816; // MeV

constexpr const char* kCutSelect =
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.3"
  " && abs(mmass-0.494)<0.1 && close_dist<=5";

struct LegendreRow
{
  Double_t mom = 0.;
  Double_t a[kNLegendre] = {};
  Double_t a_err[kNLegendre] = {};
};

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
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0
            << " <kpsc.root> [-o <out.pdf>] [-a <acceptance.root>] [-r <run>]\n"
            << "  Cut: " << kCutSelect << "\n"
            << "  Beam mom from d5_momentum mean -> acceptance index + nearest dσ.\n"
            << "  Bins with acceptance <= " << kAccMin << " are ignored.\n";
}

void
StyleCanvas(TCanvas& c)
{
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.10);
}

void
PrintProgress(Long64_t i, Long64_t n, const char* label)
{
  if (n <= 0)
    return;
  const Long64_t step = TMath::Max(n / 200LL, 1LL);
  if (i != n && (i % step) != 0)
    return;
  constexpr Int_t W = 40;
  const Double_t frac = Double_t(i) / Double_t(n);
  const Int_t filled = TMath::Min(W, Int_t(W * frac + 0.5));
  char bar[W + 1];
  for (Int_t k = 0; k < W; ++k)
    bar[k] = (k < filled) ? '#' : '-';
  bar[W] = '\0';
  std::fprintf(stderr, "\r%s [%s] %5.1f%% (%lld/%lld)",
               label, bar, 100. * frac,
               static_cast<long long>(i), static_cast<long long>(n));
  if (i >= n)
    std::fprintf(stderr, "\n");
  std::fflush(stderr);
}

Int_t
ParseRunFromTree(TTree* t)
{
  if (!t || t->GetEntries() <= 0 || !t->GetBranch("run_number"))
    return -1;
  return static_cast<Int_t>(t->GetMaximum("run_number"));
}

Int_t
ParseRunFromPath(const TString& path)
{
  TString s = path;
  s.ToLower();
  if (s.Contains("735"))
    return 735;
  const Ssiz_t n = s.Length();
  for (Ssiz_t i = 0; i + 3 < n; ++i) {
    if (s(i, 3) != "run")
      continue;
    Ssiz_t k = i + 3;
    if (k >= n || s[k] < '0' || s[k] > '9')
      continue;
    Ssiz_t e = k;
    while (e < n && s[e] >= '0' && s[e] <= '9')
      ++e;
    return TString(s(k, e - k)).Atoi();
  }
  return -1;
}

TH1D*
Book1D(const char* name, const char* title, Int_t nb, Double_t xlo, Double_t xhi)
{
  if (auto* old = gROOT->FindObject(name))
    delete old;
  auto* h = new TH1D(name, title, nb, xlo, xhi);
  h->SetDirectory(gROOT);
  h->Sumw2(kFALSE);
  return h;
}

Int_t
AcceptanceIndex(Double_t p_mev)
{
  if (p_mev < kAccMomStart
      || p_mev >= kAccMomStart + kAccMomN * kAccMomStep)
    return -1;
  return static_cast<Int_t>((p_mev - kAccMomStart) / kAccMomStep);
}

std::vector<LegendreRow>
LoadLegendreCsv(const TString& path)
{
  std::vector<LegendreRow> rows;
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "Cannot open Legendre CSV: " << path << "\n";
    return rows;
  }
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream iss(line);
    LegendreRow r;
    Double_t mom_err = 0.;
    if (!(iss >> r.mom >> mom_err))
      continue;
    bool ok = true;
    for (Int_t ell = 0; ell < kNLegendre; ++ell) {
      if (!(iss >> r.a[ell] >> r.a_err[ell])) {
        ok = false;
        break;
      }
    }
    if (ok)
      rows.push_back(r);
  }
  return rows;
}

const LegendreRow*
NearestLegendre(const std::vector<LegendreRow>& rows, Double_t p_mev)
{
  const LegendreRow* best = nullptr;
  Double_t best_dp = std::numeric_limits<Double_t>::max();
  for (const auto& r : rows) {
    const Double_t dp = TMath::Abs(r.mom - p_mev);
    if (dp < best_dp) {
      best_dp = dp;
      best = &r;
    }
  }
  return best;
}

void
EvalLegendre(const LegendreRow& row, Double_t cos_th, Double_t& y, Double_t& ey)
{
  y = 0.;
  Double_t ey2 = 0.;
  for (Int_t ell = 0; ell < kNLegendre; ++ell) {
    const Double_t pl = ROOT::Math::legendre(ell, cos_th);
    y += row.a[ell] * pl;
    ey2 += TMath::Power(row.a_err[ell] * pl, 2);
  }
  ey = TMath::Sqrt(ey2);
}

// d5_momentum is GeV/c; √s returned in MeV (same as cusp_study fit_functions).
Double_t
SqrtSFromD5GeV(Double_t p_gev)
{
  const Double_t p = p_gev * 1000.;
  const Double_t ek = TMath::Sqrt(p * p + kMassKaon * kMassKaon);
  return TMath::Sqrt(TMath::Power(ek + kMassProton, 2.) - p * p);
}

void
NormalizeUnitDensity(TH1D* h)
{
  if (!h)
    return;
  const Double_t integ = h->Integral("width");
  if (integ > 0.)
    h->Scale(1. / integ);
}

TH1D*
MakeCorrected(const TH1D* h_raw, const TH1D* h_acc, const char* name)
{
  auto* h = static_cast<TH1D*>(h_raw->Clone(name));
  h->SetDirectory(nullptr);
  h->SetTitle("acceptance-corrected cos#theta_{K}^{CM};cos#theta_{K}^{CM};normalized");
  h->Reset("ICES");
  for (Int_t i = 1; i <= h_raw->GetNbinsX(); ++i) {
    const Double_t a = h_acc->GetBinContent(i);
    const Double_t n = h_raw->GetBinContent(i);
    if (a <= kAccMin)
      continue;
    h->SetBinContent(i, n / a);
  }
  NormalizeUnitDensity(h);
  return h;
}

TH1D*
MaskLowAcceptance(const TH1D* h_acc, const char* name)
{
  auto* h = static_cast<TH1D*>(h_acc->Clone(name));
  h->SetDirectory(nullptr);
  for (Int_t i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) <= kAccMin)
      h->SetBinContent(i, 0.);
  }
  return h;
}

// Sparse markers from Legendre (A_ell only). Norm uses all valid bins (fine grid).
TGraphErrors*
MakeRefGraph(const char* name, const LegendreRow& row, const TH1D* h_acc,
             Color_t col, Style_t mstyle)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  const Double_t dx = h_acc->GetXaxis()->GetBinWidth(1);

  Double_t integ = 0.;
  for (Int_t i = 1; i <= h_acc->GetNbinsX(); ++i) {
    if (h_acc->GetBinContent(i) <= kAccMin)
      continue;
    Double_t y = 0., ey = 0.;
    EvalLegendre(row, h_acc->GetXaxis()->GetBinCenter(i), y, ey);
    if (y < 0.)
      y = 0.;
    integ += y * dx;
  }
  const Double_t scale = (integ > 0.) ? (1. / integ) : 1.;

  Double_t last_x = -999.;
  Int_t ip = 0;
  for (Int_t i = 1; i <= h_acc->GetNbinsX(); ++i) {
    if (h_acc->GetBinContent(i) <= kAccMin)
      continue;
    const Double_t x = h_acc->GetXaxis()->GetBinCenter(i);
    if (ip > 0 && (x - last_x) < kRefCosStep - 1e-9)
      continue;
    Double_t y = 0., ey = 0.;
    EvalLegendre(row, x, y, ey);
    if (y < 0.)
      y = 0.;
    g->SetPoint(ip, x, y * scale);
    g->SetPointError(ip, 0.5 * kRefCosStep, ey * scale);
    last_x = x;
    ++ip;
  }

  g->SetMarkerStyle(mstyle);
  g->SetMarkerColor(col);
  g->SetLineColor(col);
  g->SetMarkerSize(1.1);
  return g;
}

struct FillResult
{
  TH1D* h_p = nullptr;
  TH1D* h_raw = nullptr;
  TH1D* h_mmass_nt2 = nullptr;
  TH1D* h_mmass_good = nullptr;
  TH1D* h_mmass_sel = nullptr;
  TH2D* h_sqrts_costh = nullptr;
  Long64_t n_sel = 0;
};

TH2D*
Book2D(const char* name, const char* title,
       Int_t nx, Double_t xlo, Double_t xhi,
       Int_t ny, Double_t ylo, Double_t yhi)
{
  if (auto* old = gROOT->FindObject(name))
    delete old;
  auto* h = new TH2D(name, title, nx, xlo, xhi, ny, ylo, yhi);
  h->SetDirectory(gROOT);
  h->Sumw2(kFALSE);
  return h;
}

FillResult
FillHistograms(TTree* t, Int_t nb_cos, Double_t xlo, Double_t xhi)
{
  FillResult out;
  gROOT->cd();
  out.h_p = Book1D("h_d5mom_sel",
                   "d5_momentum (selected);p_{D5} [GeV/c];Counts",
                   200, 0.5, 1.0);
  out.h_raw = Book1D("h_costh_raw",
                     "cos#theta_{K}^{CM} (selected);cos#theta_{K}^{CM};Counts",
                     nb_cos, xlo, xhi);
  out.h_mmass_nt2 = Book1D("h_mmass_nt2",
                           "missing mass; M_{miss} [GeV];Counts",
                           120, 0., 1.2);
  out.h_mmass_good = Book1D("h_mmass_good",
                            "missing mass; M_{miss} [GeV];Counts",
                            120, 0., 1.2);
  out.h_mmass_sel = Book1D("h_mmass_sel",
                           "missing mass; M_{miss} [GeV];Counts",
                           120, 0., 1.2);
  out.h_sqrts_costh = Book2D(
    "h_sqrts_costh_sel",
    "selected: #sqrt{s} vs cos#theta_{K}^{CM};#sqrt{s} [MeV];cos#theta_{K}^{CM}",
    80, 1600., 1760., nb_cos, xlo, xhi);

  Double_t mmass = 0., diff_angle = 0., close_dist = 0.;
  Double_t d5_momentum = 0., cos_theta_cm = 0.;
  Int_t effective_ntTpc = 0;

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("mmass", 1);
  t->SetBranchStatus("effective_ntTpc", 1);
  t->SetBranchStatus("diff_angle", 1);
  t->SetBranchStatus("close_dist", 1);
  t->SetBranchStatus("d5_momentum", 1);
  t->SetBranchStatus("cos_theta_cm", 1);

  t->SetBranchAddress("mmass", &mmass);
  t->SetBranchAddress("effective_ntTpc", &effective_ntTpc);
  t->SetBranchAddress("diff_angle", &diff_angle);
  t->SetBranchAddress("close_dist", &close_dist);
  t->SetBranchAddress("d5_momentum", &d5_momentum);
  t->SetBranchAddress("cos_theta_cm", &cos_theta_cm);

  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    PrintProgress(i + 1, n, "fill");

    if (!(mmass > 0.))
      continue;

    const bool nt2 = (effective_ntTpc == 2);
    if (nt2)
      out.h_mmass_nt2->Fill(mmass);

    const bool good = nt2 && (diff_angle < 0.3)
                      && (TMath::Abs(mmass - kMK) < kMmassWindow);
    if (good)
      out.h_mmass_good->Fill(mmass);

    const bool sel = good && (close_dist <= 5.);
    if (!sel)
      continue;
    out.h_mmass_sel->Fill(mmass);
    out.h_p->Fill(d5_momentum);
    out.h_raw->Fill(cos_theta_cm);
    out.h_sqrts_costh->Fill(SqrtSFromD5GeV(d5_momentum), cos_theta_cm);
    ++out.n_sel;
  }

  t->SetBranchStatus("*", 1);
  t->ResetBranchAddresses();
  return out;
}

} // namespace

int
main(int argc, char** argv)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }

  TString inPath;
  TString outPdf;
  TString accPath =
    Form("%s/trigger_study_E72/results/root/acceptance_kp_flat_9999.root",
         ANALYZER_DIR.Data());
  Int_t runArg = -1;

  for (Int_t i = 1; i < argc; ++i) {
    const TString a = argv[i];
    if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "-a" && i + 1 < argc) {
      accPath = argv[++i];
    } else if (a == "-r" && i + 1 < argc) {
      runArg = TString(argv[++i]).Atoi();
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.BeginsWith("-")) {
      std::cerr << "Unknown option: " << a << "\n";
      usage(argv[0]);
      return 1;
    } else if (inPath.IsNull()) {
      inPath = a;
    } else {
      std::cerr << "Unexpected argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }

  if (inPath.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  auto* fin = TFile::Open(inPath, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open " << inPath << "\n";
    return 1;
  }
  auto* t = dynamic_cast<TTree*>(fin->Get("kpsc"));
  if (!t) {
    std::cerr << "No TTree kpsc in " << inPath << "\n";
    return 1;
  }

  Int_t runOut = runArg;
  if (runOut < 0)
    runOut = ParseRunFromTree(t);
  if (runOut < 0)
    runOut = ParseRunFromPath(inPath);
  if (runOut < 0)
    runOut = 0;

  if (outPdf.IsNull()) {
    const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
    outPdf = Form("%s/kpsc_angular_compare_run%05d.pdf", imgDir.Data(), runOut);
  }

  // Need acceptance binning before fill (cos hist).
  // Temporary: open acceptance with a placeholder index after fill for p mean.
  // Fill first with default [-1,1]/50, then re-project if needed — keep 50 bins.
  auto fill = FillHistograms(t, kNCosBins, -1., 1.);
  if (fill.n_sel <= 0) {
    std::cerr << "No events after selection cut.\n";
    return 1;
  }

  auto* h_p = fill.h_p;
  auto* h_raw = fill.h_raw;
  const Double_t p_gev = h_p->GetMean();
  const Double_t p_mev = p_gev * 1000.;
  const Int_t accIdx = AcceptanceIndex(p_mev);
  if (accIdx < 0) {
    std::cerr << "Beam mean p=" << p_mev << " MeV/c outside acceptance grid.\n";
    return 1;
  }

  auto* facc = TFile::Open(accPath, "READ");
  if (!facc || facc->IsZombie()) {
    std::cerr << "Cannot open acceptance file " << accPath << "\n";
    return 1;
  }
  auto* h_acc = dynamic_cast<TH1D*>(facc->Get(Form("acceptance%d", accIdx)));
  if (!h_acc) {
    std::cerr << "Missing acceptance" << accIdx << " in " << accPath << "\n";
    return 1;
  }
  h_acc = static_cast<TH1D*>(h_acc->Clone("h_acceptance_use"));
  h_acc->SetDirectory(nullptr);

  // Rebin raw cos to match acceptance edges if needed (same 50 on [-1,1]).
  if (h_acc->GetNbinsX() != h_raw->GetNbinsX()
      || h_acc->GetXaxis()->GetXmin() != h_raw->GetXaxis()->GetXmin()
      || h_acc->GetXaxis()->GetXmax() != h_raw->GetXaxis()->GetXmax()) {
    std::cerr << "Warning: acceptance binning differs from fill hist; "
                 "using acceptance edges via clone contents best-effort.\n";
  }

  const TString sparkCsv =
    Form("%s/cusp_study_E72/data/Kp_spark_chamber_legendre.csv",
         ANALYZER_DIR.Data());
  const TString bubbleCsv =
    Form("%s/cusp_study_E72/data/Kp_bubble_chamber1970_legendre.csv",
         ANALYZER_DIR.Data());
  const auto sparkRows = LoadLegendreCsv(sparkCsv);
  const auto bubbleRows = LoadLegendreCsv(bubbleCsv);
  const LegendreRow* spark = NearestLegendre(sparkRows, p_mev);
  const LegendreRow* bubble = NearestLegendre(bubbleRows, p_mev);
  if (!spark || !bubble) {
    std::cerr << "Failed to load nearest Legendre rows.\n";
    return 1;
  }

  auto* h_corr = MakeCorrected(h_raw, h_acc, "h_costh_corr");
  auto* h_acc_draw = MaskLowAcceptance(h_acc, "h_acceptance_draw");
  auto* g_spark = MakeRefGraph("g_spark", *spark, h_acc, kRed + 1, 20);
  auto* g_bubble = MakeRefGraph("g_bubble", *bubble, h_acc, kBlue + 1, 21);

  TCanvas c("c_kpsc_ang", "kpsc angular compare", 900, 700);
  PdfWriter writer(outPdf);

  // Page 1: cut summary
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.032);
    tx->DrawLatex(0.10, 0.90, "Kp angular compare: selection");
    tx->SetTextSize(0.024);
    tx->DrawLatex(0.10, 0.82, Form("input: %s", inPath.Data()));
    tx->DrawLatex(0.10, 0.76, Form("cut: %s", kCutSelect));
    tx->DrawLatex(0.10, 0.70, Form("selected entries = %lld",
                                   static_cast<long long>(fill.n_sel)));
    tx->DrawLatex(0.10, 0.64,
                  Form("d5_momentum mean = %.4f GeV/c (%.1f MeV/c), rms = %.4f",
                       p_gev, p_mev, h_p->GetRMS()));
    tx->DrawLatex(0.10, 0.58, Form("acceptance index = %d  (%s)",
                                   accIdx, h_acc->GetTitle()));
    tx->DrawLatex(0.10, 0.52,
                  Form("spark ref mom = %.0f MeV/c (%d pts)  bubble = %.0f (%d pts)",
                       spark->mom, g_spark->GetN(), bubble->mom, g_bubble->GetN()));
    tx->SetTextSize(0.022);
    tx->DrawLatex(0.10, 0.40,
                  Form("Bins with A #leq %.0f%% ignored. Ref markers every #Delta cos#theta=%.2f.",
                       100. * kAccMin, kRefCosStep));
    tx->DrawLatex(0.10, 0.34,
                  "Acceptance = G4 trigger/geometry only (not reconstruction/PID).");
    tx->DrawLatex(0.10, 0.28,
                  "Refs from Legendre A_{#ell}#pmerr (not digitized dense cos#theta data).");
    writer.Print(c);
  }

  // Page 2: missing mass comparison
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    auto* h_nt2 = fill.h_mmass_nt2;
    auto* h_good = fill.h_mmass_good;
    auto* h_sel = fill.h_mmass_sel;
    h_nt2->SetLineColor(kGray + 2);
    h_nt2->SetLineWidth(2);
    h_good->SetLineColor(kBlue + 1);
    h_good->SetLineWidth(2);
    h_sel->SetLineColor(kRed + 1);
    h_sel->SetLineWidth(2);
    h_nt2->SetTitle("missing mass vs cut; M_{miss} [GeV];Counts");
    h_nt2->Draw("hist");
    h_good->Draw("hist same");
    h_sel->Draw("hist same");

    auto* lLo = new TLine(kMK - kMmassWindow, 0., kMK - kMmassWindow,
                          h_nt2->GetMaximum());
    auto* lHi = new TLine(kMK + kMmassWindow, 0., kMK + kMmassWindow,
                          h_nt2->GetMaximum());
    for (auto* l : {lLo, lHi}) {
      l->SetLineColor(kOrange + 7);
      l->SetLineStyle(2);
      l->Draw("same");
    }

    auto* leg = new TLegend(0.50, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_nt2, Form("Nt==2 (%.0f)", h_nt2->GetEntries()), "l");
    leg->AddEntry(h_good,
                  Form("good (no close_dist) (%.0f)", h_good->GetEntries()), "l");
    leg->AddEntry(h_sel, Form("selected (%.0f)", h_sel->GetEntries()), "l");
    leg->Draw();
    writer.Print(c);
  }

  // Page 3: √s vs cosθ (selected)
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    c.SetRightMargin(0.14);
    auto* h2 = fill.h_sqrts_costh;
    h2->Draw("colz");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.024);
    tx->DrawLatex(0.12, 0.92,
                  Form("selected N=%.0f  (#sqrt{s} from d5_momentum, LH2 rest frame)",
                       h2->GetEntries()));
    writer.Print(c);
    c.SetRightMargin(0.05);
  }

  // Page 4: beam mom
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    h_p->SetLineColor(kBlack);
    h_p->SetLineWidth(2);
    h_p->Draw("hist");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.028);
    tx->DrawLatex(0.14, 0.92,
                  Form("mean=%.4f GeV/c  -> acceptance%d", p_gev, accIdx));
    writer.Print(c);
  }

  // Page 5: raw cos theta
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    h_raw->SetLineColor(kBlack);
    h_raw->SetLineWidth(2);
    h_raw->Draw("hist");
    writer.Print(c);
  }

  // Page 6: acceptance
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    h_acc_draw->SetTitle(Form("acceptance%d (A>%.0f%%);cos#theta_{CM};A",
                              accIdx, 100. * kAccMin));
    h_acc_draw->SetLineColor(kBlack);
    h_acc_draw->SetLineWidth(2);
    h_acc_draw->GetYaxis()->SetRangeUser(0., 1.05);
    h_acc_draw->Draw("hist");
    auto* cutLine = new TLine(-1., kAccMin, 1., kAccMin);
    cutLine->SetLineColor(kRed + 1);
    cutLine->SetLineStyle(2);
    cutLine->Draw("same");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.024);
    tx->DrawLatex(0.14, 0.92, Form("%s", h_acc->GetTitle()));
    writer.Print(c);
  }

  // Page 7: corrected overlay
  {
    c.Clear();
    c.cd();
    StyleCanvas(c);
    h_corr->SetLineColor(kBlack);
    h_corr->SetLineWidth(2);
    Double_t ymax = h_corr->GetMaximum();
    for (auto* g : {g_spark, g_bubble}) {
      for (Int_t i = 0; i < g->GetN(); ++i) {
        Double_t x = 0., y = 0.;
        g->GetPoint(i, x, y);
        ymax = TMath::Max(ymax, y + g->GetErrorY(i));
      }
    }
    h_corr->SetMaximum(ymax * 1.35);
    h_corr->Draw("hist");
    g_spark->Draw("P SAME");
    g_bubble->Draw("P SAME");

    auto* leg = new TLegend(0.14, 0.68, 0.62, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_corr, "data (acc. corrected, unit norm.)", "l");
    leg->AddEntry(g_spark,
                  Form("spark %.0f MeV/c (#Delta cos=%.2f)", spark->mom, kRefCosStep),
                  "ep");
    leg->AddEntry(g_bubble,
                  Form("bubble %.0f MeV/c (#Delta cos=%.2f)", bubble->mom, kRefCosStep),
                  "ep");
    leg->Draw();

    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.022);
    tx->DrawLatex(0.14, 0.92,
                  Form("p_{beam}=%.1f MeV/c  A_{min}=%.2f  N_{sel}=%lld",
                       p_mev, kAccMin, static_cast<long long>(fill.n_sel)));
    writer.Print(c);
  }

  writer.Close(c);
  std::cout << "Wrote " << outPdf << " (" << writer.Pages() << " pages)\n"
            << "  p_beam mean = " << p_mev << " MeV/c, acceptance" << accIdx
            << ", spark=" << spark->mom << " (" << g_spark->GetN() << " pts)"
            << ", bubble=" << bubble->mom << " (" << g_bubble->GetN() << " pts)"
            << ", A_min=" << kAccMin << "\n";

  delete g_spark;
  delete g_bubble;
  delete h_corr;
  delete h_acc_draw;
  delete h_acc;
  delete fin;
  delete facc;
  return 0;
}
