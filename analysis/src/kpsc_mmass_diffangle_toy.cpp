// -*- C++ -*-
// Toy MC: K+p elastic kinematics + component-wise resolution smear.
// Check whether mmass vs diff_angle V-shape is expected from beam/proton/K smearing.
//
// Usage:
//   kpsc_mmass_diffangle_toy [-i <kpsc.root>] [-o <out.pdf>] [-r <run>]
//                            [-n <nevt>] [--pbeam <GeV/c>]
//                            [--sig-p-rel <σ>] [--sig-p-ang <rad>]
//                            [--sig-k-rel <σ>] [--sig-k-ang <rad>]
//                            [--sig-beam-rel <σ>] [--sig-beam-ang <rad>]
//                            [--mom-bias-p <σ>] [--mom-bias-beam <σ>] [--mom-bias-k <σ>]
//                            [--seed <n>]
//
// Output (default, tentative):
//   OUTPUT_DIR/img/run00000/kpsc_mmass_diffangle_toy_run00000.pdf
//
// Definitions match DstKpScattering::ComputeKinematics:
//   mmass      = (beam + target - proton).M()
//   diff_angle = angle(p_miss, p_K)   // consistency angle, not scattering angle
//   PDF x-axis for diff_angle: 0 .. 0.5 rad

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{

constexpr Double_t kMK = 0.493677;   // GeV (PDG K+)
constexpr Double_t kMP = 0.938272;   // GeV
constexpr Double_t kMKline = 0.494;  // display line (same as purity PDF)

constexpr Double_t kZoomMLo = 0.30;
constexpr Double_t kZoomMHi = 0.70;
constexpr Double_t kAngHi = 0.5;       // display range for diff_angle
constexpr Double_t kAngTip = 0.1;      // tip-like slice (matches offline tip cut)
constexpr Int_t kNbAng = 50;
constexpr Int_t kNbM = 80;

constexpr Double_t kDefPBeam = 1.0;       // GeV/c
constexpr Long64_t kDefNevt = 80000;
// Larger default smears (tendency study, not tuned to data).
constexpr Double_t kDefSigPRel = 0.15;
constexpr Double_t kDefSigPAng = 0.08;    // rad
constexpr Double_t kDefSigKRel = 0.15;
constexpr Double_t kDefSigKAng = 0.08;
constexpr Double_t kDefSigBeamRel = 0.05;
constexpr Double_t kDefSigBeamAng = 0.02;
// Multiplicative |p| bias width for "wrong momentum" (inconsistent) scenarios.
constexpr Double_t kDefMomBiasP = 0.30;
constexpr Double_t kDefMomBiasBeam = 0.20;
constexpr Double_t kDefMomBiasK = 0.30;

void
usage(const char* argv0)
{
  std::cerr
    << "Usage: " << argv0 << " [options]\n"
    << "  -i <kpsc.root>     optional data overlay (tree kpsc)\n"
    << "  -o <out.pdf>       output PDF\n"
    << "  -r <run>           run label for default path (default 0 -> run00000)\n"
    << "  -n <nevt>          toy events (default " << kDefNevt << ")\n"
    << "  --pbeam <GeV/c>    beam |p| (default " << kDefPBeam << ")\n"
    << "  --sig-p-rel/ang    proton |p| relative / angle smear\n"
    << "  --sig-k-rel/ang    kaon smear\n"
    << "  --sig-beam-rel/ang beam smear\n"
    << "  --mom-bias-p/beam/k  |p| wrong-scale Gaussian width (inconsistent)\n"
    << "  --seed <n>         RNG seed (default 1)\n"
    << "Default PDF: OUTPUT_DIR/img/run00000/kpsc_mmass_diffangle_toy_run00000.pdf\n"
    << "diff_angle axis: 0 .. " << kAngHi << " rad\n";
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
StyleCanvas(TCanvas& c)
{
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.14);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.10);
}

struct TrueEvt
{
  TVector3 p_beam;
  TVector3 p_p;
  TVector3 p_k;
};

struct KinObs
{
  Double_t mmass = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t diff_angle = std::numeric_limits<Double_t>::quiet_NaN();
};

KinObs
ComputeKinematics(const TVector3& p_beam, const TVector3& mom_p, const TVector3& mom_k)
{
  KinObs o;
  const Double_t e_beam = TMath::Sqrt(p_beam.Mag2() + kMK * kMK);
  const Double_t e_p = TMath::Sqrt(mom_p.Mag2() + kMP * kMP);
  const TLorentzVector lv_beam(p_beam.X(), p_beam.Y(), p_beam.Z(), e_beam);
  const TLorentzVector lv_tgt(0., 0., 0., kMP);
  const TLorentzVector lv_p(mom_p.X(), mom_p.Y(), mom_p.Z(), e_p);
  const TLorentzVector lv_miss = lv_beam + lv_tgt - lv_p;
  o.mmass = lv_miss.M();
  const TVector3 p_miss = lv_miss.Vect();
  o.diff_angle
    = (p_miss.Mag() > 1.e-12 && mom_k.Mag() > 1.e-12) ? p_miss.Angle(mom_k)
                                                      : std::numeric_limits<Double_t>::quiet_NaN();
  return o;
}

Bool_t
GenerateElastic(Double_t p_beam_mag, TRandom3& rng, TrueEvt& out)
{
  if (!(p_beam_mag > 0.))
    return false;

  const TLorentzVector lv_beam(0., 0., p_beam_mag,
                               TMath::Sqrt(p_beam_mag * p_beam_mag + kMK * kMK));
  const TLorentzVector lv_tgt(0., 0., 0., kMP);
  const TLorentzVector lv_cms = lv_beam + lv_tgt;
  const TVector3 beta = lv_cms.BoostVector();

  TLorentzVector lv_beam_cm = lv_beam;
  lv_beam_cm.Boost(-beta);
  const Double_t p_star = lv_beam_cm.P();
  if (!(p_star > 0.))
    return false;

  const Double_t cos_th = rng.Uniform(-1., 1.);
  const Double_t sin_th = TMath::Sqrt(TMath::Max(0., 1. - cos_th * cos_th));
  const Double_t phi = rng.Uniform(0., TMath::TwoPi());
  const Double_t e_k_cm = TMath::Sqrt(p_star * p_star + kMK * kMK);

  TLorentzVector lv_k_cm(p_star * sin_th * TMath::Cos(phi),
                         p_star * sin_th * TMath::Sin(phi),
                         p_star * cos_th,
                         e_k_cm);
  TLorentzVector lv_k = lv_k_cm;
  lv_k.Boost(beta);
  const TLorentzVector lv_p = lv_cms - lv_k;

  out.p_beam = lv_beam.Vect();
  out.p_p = lv_p.Vect();
  out.p_k = lv_k.Vect();
  return out.p_p.Mag() > 0. && out.p_k.Mag() > 0.;
}

void
SmearMomentum(TVector3& p, Double_t sig_rel, Double_t sig_ang, TRandom3& rng)
{
  if (p.Mag() < 1.e-12)
    return;
  if (sig_rel > 0.) {
    const Double_t f = 1. + rng.Gaus(0., sig_rel);
    if (f > 0.05)
      p *= f;
  }
  if (sig_ang > 0.) {
    const Double_t mag = p.Mag();
    const TVector3 n = p.Unit();
    TVector3 a = n.Orthogonal();
    if (a.Mag2() < 1.e-24)
      a = TVector3(1., 0., 0.).Cross(n);
    a = a.Unit();
    const TVector3 b = n.Cross(a).Unit();
    const Double_t dx = rng.Gaus(0., sig_ang);
    const Double_t dy = rng.Gaus(0., sig_ang);
    const TVector3 nd = (n + dx * a + dy * b).Unit();
    p = nd * mag;
  }
}

TVector3
RandomDirection(TRandom3& rng)
{
  const Double_t cos_th = rng.Uniform(-1., 1.);
  const Double_t sin_th = TMath::Sqrt(TMath::Max(0., 1. - cos_th * cos_th));
  const Double_t phi = rng.Uniform(0., TMath::TwoPi());
  return TVector3(sin_th * TMath::Cos(phi), sin_th * TMath::Sin(phi), cos_th);
}

enum class Scenario
{
  POnly,
  BeamOnly,
  KOnly,
  All,
  PMomWrong,
  BeamMomWrong,
  KMomWrong,
  SwapPK,
  WrongK
};

const char*
ScenarioName(Scenario s)
{
  switch (s) {
  case Scenario::POnly: return "p_only";
  case Scenario::BeamOnly: return "beam_only";
  case Scenario::KOnly: return "k_only";
  case Scenario::All: return "all";
  case Scenario::PMomWrong: return "p_mom_wrong";
  case Scenario::BeamMomWrong: return "beam_mom_wrong";
  case Scenario::KMomWrong: return "k_mom_wrong";
  case Scenario::SwapPK: return "swap_pk";
  case Scenario::WrongK: return "wrong_k";
  }
  return "unknown";
}

const char*
ScenarioTitle(Scenario s)
{
  switch (s) {
  case Scenario::POnly: return "toy: proton smear only (large; expect V)";
  case Scenario::BeamOnly: return "toy: beam smear only (large; expect V-like)";
  case Scenario::KOnly: return "toy: K smear only (large; expect horizontal)";
  case Scenario::All: return "toy: beam+proton+K smear (large)";
  case Scenario::PMomWrong: return "inconsistent: wrong |p_p| (dir fixed)";
  case Scenario::BeamMomWrong: return "inconsistent: wrong |p_beam| (dir fixed)";
  case Scenario::KMomWrong: return "inconsistent: wrong |p_K| (dir fixed; angle~unchanged)";
  case Scenario::SwapPK: return "inconsistent: swap p <-> K momentum vectors";
  case Scenario::WrongK: return "inconsistent: wrong K direction (true |p_K|)";
  }
  return "toy";
}

struct SmearCfg
{
  Double_t p_rel = kDefSigPRel;
  Double_t p_ang = kDefSigPAng;
  Double_t k_rel = kDefSigKRel;
  Double_t k_ang = kDefSigKAng;
  Double_t beam_rel = kDefSigBeamRel;
  Double_t beam_ang = kDefSigBeamAng;
  Double_t mom_bias_p = kDefMomBiasP;
  Double_t mom_bias_beam = kDefMomBiasBeam;
  Double_t mom_bias_k = kDefMomBiasK;
};

Double_t
RandomMomScale(TRandom3& rng, Double_t bias_width)
{
  if (!(bias_width > 0.))
    return 1.;
  // Keep scale positive and not near zero.
  return TMath::Max(0.25, 1. + rng.Gaus(0., bias_width));
}

void
ScaleMomentumMag(TVector3& p, Double_t factor)
{
  if (factor > 0.05)
    p *= factor;
}

KinObs
ApplyScenario(const TrueEvt& truth, Scenario s, const SmearCfg& cfg, TRandom3& rng)
{
  TVector3 p_beam = truth.p_beam;
  TVector3 p_p = truth.p_p;
  TVector3 p_k = truth.p_k;

  switch (s) {
  case Scenario::POnly:
    SmearMomentum(p_p, cfg.p_rel, cfg.p_ang, rng);
    break;
  case Scenario::BeamOnly:
    SmearMomentum(p_beam, cfg.beam_rel, cfg.beam_ang, rng);
    break;
  case Scenario::KOnly:
    SmearMomentum(p_k, cfg.k_rel, cfg.k_ang, rng);
    break;
  case Scenario::All:
    SmearMomentum(p_beam, cfg.beam_rel, cfg.beam_ang, rng);
    SmearMomentum(p_p, cfg.p_rel, cfg.p_ang, rng);
    SmearMomentum(p_k, cfg.k_rel, cfg.k_ang, rng);
    break;
  case Scenario::PMomWrong:
    ScaleMomentumMag(p_p, RandomMomScale(rng, cfg.mom_bias_p));
    break;
  case Scenario::BeamMomWrong:
    ScaleMomentumMag(p_beam, RandomMomScale(rng, cfg.mom_bias_beam));
    break;
  case Scenario::KMomWrong:
    ScaleMomentumMag(p_k, RandomMomScale(rng, cfg.mom_bias_k));
    break;
  case Scenario::SwapPK:
    p_p = truth.p_k;
    p_k = truth.p_p;
    break;
  case Scenario::WrongK:
    p_k = RandomDirection(rng) * truth.p_k.Mag();
    break;
  }
  return ComputeKinematics(p_beam, p_p, p_k);
}

TH2D*
Book2D(const char* name, const char* title)
{
  auto* h = new TH2D(name, title, kNbAng, 0., kAngHi, kNbM, kZoomMLo, kZoomMHi);
  h->SetDirectory(nullptr);
  h->Sumw2(false);
  return h;
}

void
DrawMmassMkLine(TH2* /*h*/)
{
  auto* l = new TLine(0., kMKline, kAngHi, kMKline);
  l->SetLineColor(kRed + 1);
  l->SetLineWidth(2);
  l->Draw("same");
}

struct PageLayout
{
  TPad* plot = nullptr;
  TPad* note = nullptr;
};

PageLayout
BeginNotePage(TCanvas& c, Bool_t zPalette)
{
  c.Clear();
  c.cd();
  auto* plot = new TPad("toy_plot", "", 0.01, 0.04, 0.68, 0.96);
  plot->SetFillStyle(0);
  plot->SetLeftMargin(0.14);
  plot->SetRightMargin(zPalette ? 0.14 : 0.05);
  plot->SetBottomMargin(0.14);
  plot->SetTopMargin(0.10);
  plot->Draw();

  auto* note = new TPad("toy_note", "", 0.69, 0.04, 0.99, 0.96);
  note->SetFillStyle(0);
  note->SetLeftMargin(0.04);
  note->SetRightMargin(0.02);
  note->Draw();
  return {plot, note};
}

void
DrawNoteLines(TPad* note, const std::vector<TString>& lines)
{
  if (!note)
    return;
  note->cd();
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextAlign(13);
  tx->SetTextFont(42);
  Double_t y = 0.96;
  for (const auto& line : lines) {
    if (line.IsNull()) {
      y -= 0.025;
      continue;
    }
    // Section headers slightly larger.
    const Bool_t header = line.EndsWith(":") || line.BeginsWith("#bf");
    tx->SetTextSize(header ? 0.055 : 0.045);
    tx->DrawLatex(0.02, y, line);
    y -= header ? 0.055 : 0.048;
    if (y < 0.04)
      break;
  }
}

std::vector<TString>
CommonKinNotes()
{
  return {
    "#bf{definitions}",
    "M_{miss}=(p_{b}+p_{t}-p_{p}).M()",
    "p_{miss}=p_{b}+p_{t}-p_{p}",
    "diff_angle=#angle(p_{miss},p_{K})",
    "(not scatter angle)",
    "",
    Form("axis: 0-%.2f rad", kAngHi),
    Form("red line: m_{K}=%.3f", kMKline),
  };
}

std::vector<TString>
ScenarioNoteLines(Scenario s, const SmearCfg& cfg)
{
  std::vector<TString> lines;
  lines.push_back(TString("#bf{") + ScenarioName(s) + "}");
  lines.push_back("");
  {
    const std::vector<TString> common = CommonKinNotes();
    lines.insert(lines.end(), common.begin(), common.end());
  }
  lines.push_back("");
  lines.push_back("#bf{this page}");

  switch (s) {
  case Scenario::POnly:
    lines.push_back("smear: proton only");
    lines.push_back(Form("|p|: 1+Gaus(0,%.2f)", cfg.p_rel));
    lines.push_back(Form("ang: Gaus(0,%.2f) rad", cfg.p_ang));
    lines.push_back("beam, K: truth");
    lines.push_back("");
    lines.push_back("expect: V near m_{K}");
    lines.push_back("(p moves p_{miss})");
    break;
  case Scenario::BeamOnly:
    lines.push_back("smear: beam only");
    lines.push_back(Form("|p|: 1+Gaus(0,%.2f)", cfg.beam_rel));
    lines.push_back(Form("ang: Gaus(0,%.2f) rad", cfg.beam_ang));
    lines.push_back("p, K: truth");
    lines.push_back("");
    lines.push_back("expect: V-like");
    break;
  case Scenario::KOnly:
    lines.push_back("smear: K only");
    lines.push_back(Form("|p|: 1+Gaus(0,%.2f)", cfg.k_rel));
    lines.push_back(Form("ang: Gaus(0,%.2f) rad", cfg.k_ang));
    lines.push_back("beam, p: truth");
    lines.push_back("");
    lines.push_back("M_{miss} stays ~m_{K}");
    lines.push_back("expect: horizontal");
    break;
  case Scenario::All:
    lines.push_back("smear: beam+p+K");
    lines.push_back(Form("p: rel %.2f ang %.2f", cfg.p_rel, cfg.p_ang));
    lines.push_back(Form("beam: %.2f / %.2f", cfg.beam_rel, cfg.beam_ang));
    lines.push_back(Form("K: %.2f / %.2f", cfg.k_rel, cfg.k_ang));
    lines.push_back("");
    lines.push_back("closest to data-like");
    break;
  case Scenario::PMomWrong:
    lines.push_back("inconsistent |p_p|");
    lines.push_back("direction fixed");
    lines.push_back(Form("|p_p|*=1+Gaus(0,%.2f)", cfg.mom_bias_p));
    lines.push_back("no ang smear");
    lines.push_back("");
    lines.push_back("expect: V from |p|");
    break;
  case Scenario::BeamMomWrong:
    lines.push_back("inconsistent |p_beam|");
    lines.push_back("direction fixed");
    lines.push_back(Form("|p_b|*=1+Gaus(0,%.2f)", cfg.mom_bias_beam));
    lines.push_back("no ang smear");
    lines.push_back("");
    lines.push_back("expect: V-like");
    break;
  case Scenario::KMomWrong:
    lines.push_back("inconsistent |p_K|");
    lines.push_back("direction fixed");
    lines.push_back(Form("|p_K|*=1+Gaus(0,%.2f)", cfg.mom_bias_k));
    lines.push_back("");
    lines.push_back("M_{miss} untouched");
    lines.push_back("diff_angle ~unchanged");
    lines.push_back("(angle uses dir)");
    break;
  case Scenario::SwapPK:
    lines.push_back("inconsistent assign");
    lines.push_back("p <-> K vectors");
    lines.push_back("p_p(use)=p_K(true)");
    lines.push_back("p_K(use)=p_p(true)");
    lines.push_back("");
    lines.push_back("expect: off tip");
    break;
  case Scenario::WrongK:
    lines.push_back("inconsistent K dir");
    lines.push_back("p, beam: truth");
    lines.push_back("|p_K| kept");
    lines.push_back("dir: isotropic");
    lines.push_back("");
    lines.push_back("M_{miss}~m_{K}");
    lines.push_back("expect: horizontal");
    break;
  }
  return lines;
}

std::vector<TString>
DataNoteLines()
{
  std::vector<TString> lines;
  lines.push_back("#bf{DATA Nt2}");
  lines.push_back("");
  {
    const std::vector<TString> common = CommonKinNotes();
    lines.insert(lines.end(), common.begin(), common.end());
  }
  lines.push_back("");
  lines.push_back("#bf{selection}");
  lines.push_back("M_{miss} > 0");
  lines.push_back("N_{TPC}^{eff} = 2");
  lines.push_back("");
  lines.push_back("compare to toy V");
  return lines;
}

std::vector<TString>
TruthNoteLines()
{
  std::vector<TString> lines;
  lines.push_back("#bf{truth (no smear)}");
  lines.push_back("");
  {
    const std::vector<TString> common = CommonKinNotes();
    lines.insert(lines.end(), common.begin(), common.end());
  }
  lines.push_back("");
  lines.push_back("#bf{this page}");
  lines.push_back("ideal elastic 2->2");
  lines.push_back("no resolution");
  lines.push_back("");
  lines.push_back("expect: tip only");
  lines.push_back("diff_angle~0");
  lines.push_back("M_{miss}~m_{K}");
  return lines;
}

std::vector<TString>
SliceNoteLines(const char* toyLabel, Bool_t hasData)
{
  std::vector<TString> lines = {
    "#bf{angle slices}",
    Form("toy: %s", toyLabel),
    "",
    "#bf{definitions}",
    "M_{miss}=(p_{b}+p_{t}-p_{p}).M()",
    "diff_angle=#angle(p_{miss},p_{K})",
    "",
    "#bf{slices}",
    Form("blue: #Delta#theta < %.2f", kAngTip),
    Form("orange: #Delta#theta #geq %.2f", kAngTip),
    "(normalized 1D)",
    "",
  };
  if (hasData) {
    lines.push_back("markers: data Nt2");
    lines.push_back("same angle slices");
  } else {
    lines.push_back("no data overlay");
  }
  return lines;
}

void
DrawTitlePage(TCanvas& c, PdfWriter& writer, Double_t pbeam, Long64_t nevt,
              const SmearCfg& cfg, Bool_t hasData)
{
  c.Clear();
  c.cd();
  StyleCanvas(c);
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.035);
  tx->DrawLatex(0.08, 0.90, "Kp elastic toy: mmass vs diff_angle V study");
  tx->SetTextSize(0.026);
  tx->DrawLatex(0.08, 0.82,
                "M_{miss}=(p_{b}+p_{t}-p_{p}).M();  diff_angle=#angle(p_{miss},p_{K})");
  tx->DrawLatex(0.08, 0.76,
                "diff_angle = consistency angle (not lab / CM scattering angle)");
  tx->SetTextSize(0.028);
  Double_t y = 0.66;
  tx->DrawLatex(0.08, y, Form("p_{beam} = %.3f GeV/c,  N_{toy} = %lld%s",
                              pbeam, nevt, hasData ? ",  + data kpsc overlay" : ""));
  y -= 0.05;
  tx->DrawLatex(0.08, y,
                Form("#sigma_p: rel=%.3f  ang=%.3f rad", cfg.p_rel, cfg.p_ang));
  y -= 0.045;
  tx->DrawLatex(0.08, y,
                Form("#sigma_{beam}: rel=%.3f  ang=%.3f rad", cfg.beam_rel, cfg.beam_ang));
  y -= 0.045;
  tx->DrawLatex(0.08, y,
                Form("#sigma_K: rel=%.3f  ang=%.3f rad", cfg.k_rel, cfg.k_ang));
  y -= 0.045;
  tx->DrawLatex(0.08, y,
                Form("|p| bias (inconsistent): p=%.2f  beam=%.2f  K=%.2f",
                     cfg.mom_bias_p, cfg.mom_bias_beam, cfg.mom_bias_k));
  y -= 0.055;
  tx->SetTextSize(0.024);
  tx->DrawLatex(0.08, y,
                Form("Each plot page: left figure, right note (what/how much).  x: 0-%.2f rad",
                     kAngHi));
  y -= 0.045;
  tx->DrawLatex(0.08, y, "Expect V: p/beam smear or wrong |p|;  horizontal: K dir / wrong_k");
  y -= 0.045;
  tx->DrawLatex(0.08, y, "cos#theta_{CM} flat (shape tendency only; not d#sigma)");
  writer.Print(c);
}

void
Draw2DPage(TCanvas& c, PdfWriter& writer, TH2D* h, const char* title,
           const std::vector<TString>& noteLines)
{
  auto pads = BeginNotePage(c, true);
  pads.plot->cd();
  h->SetTitle(title);
  h->GetXaxis()->SetTitle("#angle(p_{miss}, p_{K}) [rad]");
  h->GetYaxis()->SetTitle("M_{miss} [GeV]");
  h->GetXaxis()->CenterTitle(false);
  h->GetYaxis()->CenterTitle(false);
  h->Draw("COLZ");
  DrawMmassMkLine(h);
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.035);
  tx->DrawLatex(0.16, 0.93, Form("N=%.0f", h->GetEntries()));
  DrawNoteLines(pads.note, noteLines);
  c.cd();
  writer.Print(c);
}

void
DrawSlicePage(TCanvas& c, PdfWriter& writer, TH2D* hToy, TH2D* hData, const char* label)
{
  auto pads = BeginNotePage(c, false);
  pads.plot->cd();

  // Project mmass for tip-like and wing-like angle slices (within 0..kAngHi).
  static Int_t sliceIdx = 0;
  ++sliceIdx;
  auto* hTip = hToy->ProjectionY(Form("h_slice_tip_%d", sliceIdx),
                                 1, hToy->GetXaxis()->FindBin(kAngTip - 1.e-6));
  auto* hWing = hToy->ProjectionY(Form("h_slice_wing_%d", sliceIdx),
                                  hToy->GetXaxis()->FindBin(kAngTip),
                                  hToy->GetNbinsX());
  hTip->SetDirectory(nullptr);
  hWing->SetDirectory(nullptr);
  hTip->SetLineColor(kBlue + 1);
  hWing->SetLineColor(kOrange + 7);
  hTip->SetLineWidth(2);
  hWing->SetLineWidth(2);
  if (hTip->Integral() > 0)
    hTip->Scale(1. / hTip->Integral());
  if (hWing->Integral() > 0)
    hWing->Scale(1. / hWing->Integral());

  hTip->SetTitle(Form("%s; M_{miss} [GeV]; normalized", label));
  hTip->GetYaxis()->SetRangeUser(0., 1.15 * TMath::Max(hTip->GetMaximum(), hWing->GetMaximum()));
  hTip->Draw("HIST");
  hWing->Draw("HIST SAME");

  const Bool_t hasData = (hData && hData->GetEntries() > 0);
  if (hasData) {
    auto* dTip = hData->ProjectionY(Form("h_data_tip_%d", sliceIdx),
                                    1, hData->GetXaxis()->FindBin(kAngTip - 1.e-6));
    auto* dWing = hData->ProjectionY(Form("h_data_wing_%d", sliceIdx),
                                     hData->GetXaxis()->FindBin(kAngTip),
                                     hData->GetNbinsX());
    dTip->SetDirectory(nullptr);
    dWing->SetDirectory(nullptr);
    dTip->SetMarkerStyle(20);
    dWing->SetMarkerStyle(24);
    dTip->SetMarkerColor(kBlue + 1);
    dWing->SetMarkerColor(kOrange + 7);
    dTip->SetMarkerSize(0.7);
    dWing->SetMarkerSize(0.7);
    if (dTip->Integral() > 0)
      dTip->Scale(1. / dTip->Integral());
    if (dWing->Integral() > 0)
      dWing->Scale(1. / dWing->Integral());
    dTip->Draw("E1 SAME");
    dWing->Draw("E1 SAME");
  }

  auto* leg = new TLegend(0.45, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hTip, Form("#Delta#theta<%.2f", kAngTip), "l");
  leg->AddEntry(hWing, Form("#Delta#theta#geq%.2f", kAngTip), "l");
  leg->Draw();

  auto* lMk = new TLine(kMKline, 0., kMKline, hTip->GetMaximum());
  lMk->SetLineColor(kRed + 1);
  lMk->SetLineStyle(2);
  lMk->Draw("same");

  DrawNoteLines(pads.note, SliceNoteLines(label, hasData));
  c.cd();
  writer.Print(c);
}

Bool_t
FillDataHist(TTree* t, TH2D* h)
{
  if (!t || !h)
    return false;
  Double_t mmass = 0., diff_angle = 0.;
  Int_t effective_ntTpc = 0;
  t->SetBranchStatus("*", 0);
  if (!t->GetBranch("mmass") || !t->GetBranch("diff_angle")
      || !t->GetBranch("effective_ntTpc")) {
    std::cerr << "Warning: kpsc missing required branches; skip data overlay\n";
    return false;
  }
  t->SetBranchStatus("mmass", 1);
  t->SetBranchStatus("diff_angle", 1);
  t->SetBranchStatus("effective_ntTpc", 1);
  t->SetBranchAddress("mmass", &mmass);
  t->SetBranchAddress("diff_angle", &diff_angle);
  t->SetBranchAddress("effective_ntTpc", &effective_ntTpc);
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (!(mmass > 0.) || effective_ntTpc != 2)
      continue;
    if (!std::isfinite(mmass) || !std::isfinite(diff_angle))
      continue;
    h->Fill(diff_angle, mmass);
  }
  t->SetBranchStatus("*", 1);
  return h->GetEntries() > 0;
}

} // namespace

Int_t
main(Int_t argc, Char_t** argv)
{
  TString inPath;
  TString outPdf;
  Int_t run = 0; // tentative default -> run00000
  Long64_t nevt = kDefNevt;
  Double_t pbeam = kDefPBeam;
  UInt_t seed = 1;
  SmearCfg cfg;

  for (Int_t i = 1; i < argc; ++i) {
    const TString a = argv[i];
    auto need = [&](const char* opt) -> TString {
      if (i + 1 >= argc) {
        std::cerr << "Error: " << opt << " needs a value\n";
        usage(argv[0]);
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a == "-i") {
      inPath = need("-i");
    } else if (a == "-o") {
      outPdf = need("-o");
    } else if (a == "-r") {
      run = need("-r").Atoi();
    } else if (a == "-n") {
      nevt = need("-n").Atoll();
    } else if (a == "--pbeam") {
      pbeam = need("--pbeam").Atof();
    } else if (a == "--sig-p-rel") {
      cfg.p_rel = need("--sig-p-rel").Atof();
    } else if (a == "--sig-p-ang") {
      cfg.p_ang = need("--sig-p-ang").Atof();
    } else if (a == "--sig-k-rel") {
      cfg.k_rel = need("--sig-k-rel").Atof();
    } else if (a == "--sig-k-ang") {
      cfg.k_ang = need("--sig-k-ang").Atof();
    } else if (a == "--sig-beam-rel") {
      cfg.beam_rel = need("--sig-beam-rel").Atof();
    } else if (a == "--sig-beam-ang") {
      cfg.beam_ang = need("--sig-beam-ang").Atof();
    } else if (a == "--mom-bias-p") {
      cfg.mom_bias_p = need("--mom-bias-p").Atof();
    } else if (a == "--mom-bias-beam") {
      cfg.mom_bias_beam = need("--mom-bias-beam").Atof();
    } else if (a == "--mom-bias-k") {
      cfg.mom_bias_k = need("--mom-bias-k").Atof();
    } else if (a == "--seed") {
      seed = static_cast<UInt_t>(need("--seed").Atoi());
    } else if (a.BeginsWith("-")) {
      std::cerr << "Unknown option: " << a << "\n";
      usage(argv[0]);
      return 1;
    } else if (inPath.IsNull()) {
      inPath = a; // positional kpsc.root
    } else {
      std::cerr << "Unexpected argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }

  if (nevt <= 0 || !(pbeam > 0.)) {
    std::cerr << "Error: need positive -n and --pbeam\n";
    return 1;
  }
  if (run < 0)
    run = 0;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TRandom3 rng(seed);

  std::vector<TrueEvt> truths;
  truths.reserve(static_cast<std::size_t>(nevt));
  Long64_t nGen = 0;
  while (static_cast<Long64_t>(truths.size()) < nevt && nGen < nevt * 20) {
    ++nGen;
    TrueEvt ev;
    if (!GenerateElastic(pbeam, rng, ev))
      continue;
    // Sanity: true kinematics should sit at tip.
    const KinObs kin0 = ComputeKinematics(ev.p_beam, ev.p_p, ev.p_k);
    if (!std::isfinite(kin0.mmass) || !std::isfinite(kin0.diff_angle))
      continue;
    truths.push_back(ev);
  }
  if (truths.empty()) {
    std::cerr << "Error: failed to generate elastic events\n";
    return 1;
  }

  TH2D* hData = Book2D("h_data_nt2", "data Nt2; #angle(p_{miss},p_{K}) [rad]; missing mass [GeV]");
  Bool_t hasData = false;
  TFile* fData = nullptr;
  if (!inPath.IsNull()) {
    fData = TFile::Open(inPath, "READ");
    if (!fData || fData->IsZombie()) {
      std::cerr << "Warning: cannot open " << inPath << "; continue without data\n";
    } else {
      auto* t = dynamic_cast<TTree*>(fData->Get("kpsc"));
      if (!t) {
        std::cerr << "Warning: tree kpsc not found; continue without data\n";
      } else {
        hasData = FillDataHist(t, hData);
        if (hasData)
          std::cout << "Data Nt2 entries filled: " << hData->GetEntries() << std::endl;
      }
    }
  }

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, run);
  if (outPdf.IsNull())
    outPdf = Form("%s/kpsc_mmass_diffangle_toy_run%05d.pdf", imgDir.Data(), run);

  TCanvas c("c_kpsc_toy", "kpsc mmass diff_angle toy", 900, 700);
  PdfWriter writer(outPdf);

  DrawTitlePage(c, writer, pbeam, static_cast<Long64_t>(truths.size()), cfg, hasData);

  if (hasData)
    Draw2DPage(c, writer, hData, "DATA Nt2", DataNoteLines());

  const Scenario scenarios[] = {
    Scenario::POnly,
    Scenario::BeamOnly,
    Scenario::KOnly,
    Scenario::All,
    Scenario::PMomWrong,
    Scenario::BeamMomWrong,
    Scenario::KMomWrong,
    Scenario::SwapPK,
    Scenario::WrongK,
  };

  TH2D* hPOnly = nullptr;
  TH2D* hAll = nullptr;

  for (const Scenario s : scenarios) {
    // Independent RNG stream per scenario from same truths (re-seed by scenario id).
    TRandom3 rngS(seed + 17u * (1u + static_cast<UInt_t>(s)));
    auto* h = Book2D(Form("h_toy_%s", ScenarioName(s)),
                     Form("%s; #angle(p_{miss},p_{K}) [rad]; M_{miss} [GeV]",
                          ScenarioTitle(s)));
    for (const TrueEvt& ev : truths) {
      const KinObs o = ApplyScenario(ev, s, cfg, rngS);
      if (!std::isfinite(o.mmass) || !std::isfinite(o.diff_angle))
        continue;
      h->Fill(o.diff_angle, o.mmass);
    }
    Draw2DPage(c, writer, h, ScenarioTitle(s), ScenarioNoteLines(s, cfg));
    if (s == Scenario::POnly)
      hPOnly = h;
    if (s == Scenario::All)
      hAll = h;
  }

  // Angle-slice comparison: prefer all-smear; fall back to p_only.
  TH2D* hSlice = hAll ? hAll : hPOnly;
  if (hSlice)
    DrawSlicePage(c, writer, hSlice, hasData ? hData : nullptr,
                  hAll ? "all-smear" : "p_only");

  // Truth check page: unsmeared should peak at tip.
  {
    auto* hTrue = Book2D("h_toy_truth",
                         "toy truth (no smear); #angle(p_{miss},p_{K}) [rad]; M_{miss} [GeV]");
    for (const TrueEvt& ev : truths) {
      const KinObs o = ComputeKinematics(ev.p_beam, ev.p_p, ev.p_k);
      if (!std::isfinite(o.mmass) || !std::isfinite(o.diff_angle))
        continue;
      hTrue->Fill(o.diff_angle, o.mmass);
    }
    Draw2DPage(c, writer, hTrue, "toy truth (no smear)", TruthNoteLines());
  }

  writer.Close(c);
  std::cout << "Wrote " << outPdf << " (" << writer.Pages() << " pages)" << std::endl;

  if (fData)
    fData->Close();
  return 0;
}
