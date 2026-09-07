// -*- C++ -*-
// Offline TPC dE/dx PID calibration from DstTpcDedxPidSamples tree.
//
// Usage:
//   tpc_dedx_pid_calib <root>... [-o out.pdf] [--out-user snip.txt] [-r <run>]
//     [--close-dist-max 10] [--diff-angle-max 0.3]
//     [--mass-lo X --mass-hi Y]  # Lambda tag_mass window (default 1.115±0.1)
//     [--kp-mass-lo X --kp-mass-hi Y]  # Kp missing-mass window (default m_K±0.1)
//
// Input: TTree "dedx_pid_sample" (one row = one tagged track).
// Fit/PDF: pool kp+lambda by species (pi, K, p) only — no reaction split.
// Fit: proton median dEdx/Bethe(conv=1) → TpcDedxConversionFactor;
//      residual RMS vs |p| → CalcTPCdEdxSigma 5-param (E72 init).
// Defaults:
//   PDF:  OUTPUT_DIR/img/runXXXXX/tpc_dedx_pid_calib_runXXXXX.pdf
//   User: OUTPUT_DIR/tpc_dedx_pid/tpc_dedx_pid_user_runXXXXX.txt

#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

// Copied from analysis/src/tpc_pid_dedx_boundaries.cpp (E72 / Kinematics).
constexpr Double_t kConvE72 = 7388.11;
constexpr Double_t kSigmaPiE72[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
constexpr Double_t kSigmaKE72[5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
constexpr Double_t kSigmaPE72[5] = {12.9717, -8.43799, 3.10608, 166.494, -6.56123};

// Current HypTPCdEdxPID nσ windows (Kinematics E72 defaults).
constexpr Double_t kWinPiE72[2] = {-3., 3.};
constexpr Double_t kWinKE72[2] = {-2., 2.};
constexpr Double_t kWinPE72[2] = {-1.5, 6.};
// Candidate tightening discussed for π (not applied unless UserParam changes).
constexpr Double_t kWinPiTight[2] = {-2., 2.};

constexpr Double_t kDefaultCloseDistMax = 10.;   // mm
constexpr Double_t kDefaultDiffAngleMax = 0.3;   // rad

constexpr Int_t kNSpecies = 3; // 0=pi, 1=K, 2=p

constexpr Int_t kNpBins = 48; // ~25 MeV bins on [kPLo, kPHi]
constexpr Double_t kPLo = 0.05;
constexpr Double_t kPHi = 1.25;
constexpr Int_t kMinBinEntries = 40;
// Skip ultra-low |p| for sigma RMS fit (exp term otherwise overfits the tip).
constexpr Double_t kSigmaFitPMin = 0.15;

constexpr Int_t kDedxBins = 200;
constexpr Double_t kDedxMax = 400.;
constexpr Int_t kNsigmaBins = 120;
constexpr Double_t kNsigmaMax = 12.;

// Lambda invariant-mass window (GeV/c2): 1.115 ± 0.1
constexpr Double_t kLambdaMassNom = 1.115;
constexpr Double_t kLambdaMassHalfWin = 0.1;
constexpr Double_t kDefaultLambdaMassLo = kLambdaMassNom - kLambdaMassHalfWin;
constexpr Double_t kDefaultLambdaMassHi = kLambdaMassNom + kLambdaMassHalfWin;

// Kp missing mass (tag_mass = mmass) window: m_K ± 0.1 GeV/c2 (filled at runtime).
constexpr Double_t kKpMassHalfWin = 0.1;

const char* kSpeciesName[kNSpecies] = {"pi", "K", "p"};

Double_t PdgMassGeV(int pdg_code)
{
  TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pdg_code);
  return p ? p->Mass() : -1.;
}

Double_t PionMass() { return PdgMassGeV(kPiMinus); }
Double_t KaonMass() { return PdgMassGeV(kKMinus); }
Double_t ProtonMass() { return PdgMassGeV(kProton); }
Double_t ElectronMass() { return PdgMassGeV(kElectron); }

// Analyzer HypTPCdEdxElectron (Kinematics.cc). Defaults from E72.
constexpr Double_t kElectronNsigmaPiMax = -3.5; // require nsigma_pi < this
constexpr Double_t kElectronNsigmaAbsMax = 3.5; // |nsigma_e| < this
constexpr Double_t kElectronPoqAbsMax = 0.1;    // |poq| < this [GeV/c]
constexpr Double_t kSigmaDedxEpE72 = 7.8511;
constexpr Double_t kSigmaDedxEmE72 = 8.45029;

Double_t SpeciesMassGeV(Int_t species)
{
  if (species == 0)
    return PionMass();
  if (species == 1)
    return KaonMass();
  if (species == 2)
    return ProtonMass();
  return -1.;
}

Double_t BetaFromEandP(Double_t energy, Double_t momentum)
{
  return momentum / energy;
}

Double_t DensityEffectCorrection(Double_t betagamma, Double_t* par)
{
  const Double_t constant = 2. * TMath::Log(10.);
  Double_t delta = 0.;
  const Double_t X = TMath::Log10(betagamma);
  if (X <= par[2])
    delta = par[5] * TMath::Power(10., 2. * (X - par[2]));
  else if (par[2] < X && X < par[3])
    delta = constant * X - par[4] + par[0] * TMath::Power((par[3] - X), par[1]);
  else if (X >= par[3])
    delta = constant * X - par[4];
  return delta;
}

// P10 only (materialid==0), same as Kinematics::HypTPCdEdx for gas.
Double_t HypTPCdEdxP10(Double_t mass /* MeV/c2 */, Double_t beta)
{
  const Double_t rho = TMath::Power(10., -3) * (0.9 * 1.662 + 0.1 * 0.6672);
  const Double_t ZoverA = 17.2 / 37.6;
  const Double_t I = 0.9 * 188.0 + 0.1 * 41.7;
  Double_t density_effect_par[6] = {0.};
  density_effect_par[0] = 0.9 * 0.19714 + 0.1 * 0.09253;
  density_effect_par[1] = 0.9 * 2.9618 + 0.1 * 3.6257;
  density_effect_par[2] = 0.9 * 1.7635 + 0.1 * 1.6263;
  density_effect_par[3] = 0.9 * 4.4855 + 0.1 * 3.9716;
  density_effect_par[4] = 0.9 * 11.9480 + 0.1 * 9.5243;
  density_effect_par[5] = 0.;

  const Double_t Z = 1.;
  const Double_t me = 0.5109989461;
  const Double_t K = 0.307075;
  const Double_t constant = rho * K * ZoverA;
  const Double_t I2 = I * I;
  const Double_t beta2 = beta * beta;
  const Double_t gamma2 = 1. / (1. - beta2);
  const Double_t MeVToeV = TMath::Power(10., 6);
  const Double_t Wmax =
    2. * me * beta2 * gamma2 / (TMath::Sq(me / mass + 1.) + 2. * (me / mass) * (TMath::Sqrt(gamma2) - 1.));
  const Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2 * gamma2), density_effect_par);
  const Double_t dedx =
    constant * Z * Z / beta2 *
    (0.5 * TMath::Log(2. * me * beta2 * gamma2 * Wmax * MeVToeV * MeVToeV / I2) - beta2 - 0.5 * delta);
  return dedx;
}

Double_t HypTPCBethe(Double_t poq, Double_t mass_mev, Double_t conv)
{
  const Double_t momentum = 1000. * TMath::Abs(poq);
  const Double_t energy = TMath::Hypot(mass_mev, momentum);
  const Double_t beta = BetaFromEandP(energy, momentum);
  return conv * HypTPCdEdxP10(mass_mev, beta);
}

Double_t CalcTPCdEdxSigma(const Double_t sigma_par[5], Double_t poq)
{
  const Double_t abspoq = TMath::Abs(poq);
  return sigma_par[0] + sigma_par[1] * abspoq + sigma_par[2] * poq * poq +
         sigma_par[3] * TMath::Exp(sigma_par[4] * abspoq);
}

// Same logic as Kinematics::HypTPCdEdxElectron (hardcoded nσ / |p| cuts).
Bool_t IsHypTPCdEdxElectron(Double_t dedx, Double_t poq, Double_t conv, Double_t offset_pi,
                            const Double_t sigma_pi[5])
{
  if (!(std::isfinite(dedx) && std::isfinite(poq)))
    return kFALSE;
  if (!(TMath::Abs(poq) < kElectronPoqAbsMax))
    return kFALSE;
  const Double_t m_pi = 1000. * PionMass();
  const Double_t m_e = 1000. * ElectronMass();
  const Double_t mean_pi = HypTPCBethe(poq, m_pi, conv) + offset_pi;
  const Double_t mean_e = HypTPCBethe(poq, m_e, conv);
  const Double_t sig_pi = CalcTPCdEdxSigma(sigma_pi, poq);
  const Double_t sig_e = (poq > 0.) ? kSigmaDedxEpE72 : kSigmaDedxEmE72;
  if (!(sig_pi > 0.) || !(sig_e > 0.) || !std::isfinite(mean_pi) || !std::isfinite(mean_e))
    return kFALSE;
  const Double_t nsigma_pi = (dedx - mean_pi) / sig_pi;
  const Double_t nsigma_e = (dedx - mean_e) / sig_e;
  return (nsigma_pi < kElectronNsigmaPiMax) && (TMath::Abs(nsigma_e) < kElectronNsigmaAbsMax);
}

void CopySigma(Double_t* dst, const Double_t src[5])
{
  for (Int_t i = 0; i < 5; ++i)
    dst[i] = src[i];
}

Double_t MedianSorted(std::vector<Double_t> v)
{
  if (v.empty())
    return 0.;
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  if (n % 2 == 1)
    return v[n / 2];
  return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct Sample
{
  Int_t species = -1;
  Int_t charge = 0;
  Double_t mom0 = 0.;
  Double_t dEdx = 0.;

  Double_t AbsP() const { return TMath::Abs(mom0); }
  // |poq| = |mom0|; charge*mom0 for Bethe (Abs inside HypTPCBethe).
  Double_t Poq() const { return static_cast<Double_t>(charge) * mom0; }
};

struct CutCfg
{
  Double_t close_dist_max = kDefaultCloseDistMax;
  Double_t diff_angle_max = kDefaultDiffAngleMax;
  Double_t lambda_mass_lo = kDefaultLambdaMassLo;
  Double_t lambda_mass_hi = kDefaultLambdaMassHi;
  Double_t kp_mass_lo = 0.; // set from KaonMass()±half unless CLI overrides
  Double_t kp_mass_hi = 0.;
};

bool PassCuts(Int_t source, Int_t effective_ntTpc, Double_t close_dist, Double_t diff_angle, Double_t tag_mass,
              const CutCfg& cfg)
{
  if (!(std::isfinite(close_dist) && close_dist < cfg.close_dist_max))
    return false;
  if (source == 0) {
    // Kp: quality + missing-mass (tag_mass = mmass) around m_K.
    if (effective_ntTpc != 2)
      return false;
    if (!(std::isfinite(diff_angle) && diff_angle < cfg.diff_angle_max))
      return false;
    if (!(std::isfinite(tag_mass) && tag_mass >= cfg.kp_mass_lo && tag_mass <= cfg.kp_mass_hi))
      return false;
  } else if (source == 1) {
    // Lambda: close_dist + invariant-mass window.
    if (!(std::isfinite(tag_mass) && tag_mass >= cfg.lambda_mass_lo && tag_mass <= cfg.lambda_mass_hi))
      return false;
  } else {
    return false;
  }
  return true;
}

const Double_t* SigmaE72(Int_t species)
{
  if (species == 0)
    return kSigmaPiE72;
  if (species == 1)
    return kSigmaKE72;
  return kSigmaPE72;
}

TGraph MakeBetheCurve(Double_t mass_mev, Double_t conv, Double_t offset, Color_t col)
{
  constexpr Int_t n = 200;
  std::vector<Double_t> x, y;
  x.reserve(n);
  y.reserve(n);
  for (Int_t i = 0; i < n; ++i) {
    const Double_t t = static_cast<Double_t>(i) / static_cast<Double_t>(n - 1);
    const Double_t p = kPLo + t * (kPHi - kPLo);
    if (p < 1e-4)
      continue;
    const Double_t mean = HypTPCBethe(p, mass_mev, conv) + offset;
    if (!std::isfinite(mean) || mean <= 0.)
      continue;
    x.push_back(p);
    y.push_back(mean);
  }
  TGraph g(static_cast<Int_t>(x.size()), x.data(), y.data());
  g.SetLineColor(col);
  g.SetLineWidth(2);
  return g;
}

// mean(p) + nsig * σ(p) for drawing PID window edges (dotted).
TGraph MakeBetheNsigmaEdge(Double_t mass_mev, Double_t conv, Double_t offset, const Double_t sigma_par[5],
                           Double_t nsig, Color_t col)
{
  constexpr Int_t n = 200;
  std::vector<Double_t> x, y;
  x.reserve(n);
  y.reserve(n);
  for (Int_t i = 0; i < n; ++i) {
    const Double_t t = static_cast<Double_t>(i) / static_cast<Double_t>(n - 1);
    const Double_t p = kPLo + t * (kPHi - kPLo);
    if (p < 1e-4)
      continue;
    const Double_t mean = HypTPCBethe(p, mass_mev, conv) + offset;
    const Double_t sig = CalcTPCdEdxSigma(sigma_par, p);
    if (!std::isfinite(mean) || !(sig > 0.))
      continue;
    const Double_t yy = mean + nsig * sig;
    if (!(yy > 0.))
      continue;
    x.push_back(p);
    y.push_back(yy);
  }
  TGraph g(static_cast<Int_t>(x.size()), x.data(), y.data());
  g.SetLineColor(col);
  g.SetLineWidth(1);
  g.SetLineStyle(3); // dotted
  return g;
}

const Double_t* WinE72(Int_t species)
{
  if (species == 0)
    return kWinPiE72;
  if (species == 1)
    return kWinKE72;
  return kWinPE72;
}

void StyleCanvas(TCanvas& c)
{
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.14);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.08);
}

struct PdfWriter
{
  TString path;
  Int_t n_printed = 0;
  Bool_t closed = kFALSE;

  explicit PdfWriter(const TString& p) : path(p)
  {
    gSystem->mkdir(gSystem->DirName(path), kTRUE);
    std::remove(path.Data());
  }

  void Print(TCanvas& c)
  {
    if (closed)
      return;
    if (n_printed == 0)
      c.Print(path + "(");
    else
      c.Print(path);
    ++n_printed;
  }

  void Close(TCanvas& c)
  {
    if (closed)
      return;
    if (n_printed == 0) {
      c.Clear();
      auto* tx = new TLatex(0.2, 0.5, "no pages");
      tx->SetNDC();
      tx->Draw();
      c.Print(path);
    } else {
      c.Print(path + ")");
    }
    closed = kTRUE;
  }
};

Int_t ParseRunFromPath(const TString& path)
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

Int_t ParseRunFromChain(TChain& chain)
{
  if (chain.GetEntries() <= 0)
    return -1;
  if (!chain.GetBranch("run_number"))
    return -1;
  return static_cast<Int_t>(chain.GetMaximum("run_number"));
}

void usage(const char* a0)
{
  std::cerr << "Usage: " << a0 << " <root>... [-o out.pdf] [--out-user snip.txt] [-r <run>]\n"
            << "  [--close-dist-max 10] [--diff-angle-max 0.3]\n"
            << "  [--mass-lo X --mass-hi Y]  # Lambda tag_mass window (default 1.115±0.1)\n"
            << "  [--kp-mass-lo X --kp-mass-hi Y]  # Kp mmass window (default m_K±0.1)\n"
            << "  Tree: dedx_pid_sample from DstTpcDedxPidSamples\n"
            << "  Fit/PDF: pool kp+lambda by species (pi/K/p) only\n"
            << "  Default PDF:  OUTPUT_DIR/img/runXXXXX/tpc_dedx_pid_calib_runXXXXX.pdf\n"
            << "  Default User: OUTPUT_DIR/tpc_dedx_pid/tpc_dedx_pid_user_runXXXXX.txt\n"
            << "  Default cuts: close_dist < 10 mm;\n"
            << "                kp: effective_ntTpc==2 && diff_angle < 0.3 && mmass in m_K±0.1;\n"
            << "                lambda: tag_mass in [1.015, 1.215]; no delta_phi\n";
}

Double_t MedianScaleToGraph(TGraphErrors& gr, const Double_t shape[5])
{
  std::vector<Double_t> scales;
  scales.reserve(static_cast<std::size_t>(gr.GetN()));
  for (Int_t i = 0; i < gr.GetN(); ++i) {
    Double_t x = 0., y = 0.;
    gr.GetPoint(i, x, y);
    const Double_t s0 = CalcTPCdEdxSigma(shape, x);
    if (!(s0 > 1e-6) || !(y > 0.) || !std::isfinite(y))
      continue;
    scales.push_back(y / s0);
  }
  return scales.empty() ? 1. : MedianSorted(scales);
}

void ApplySigmaSoftLimits(TF1& f)
{
  f.SetParLimits(0, 0., 200.);
  f.SetParLimits(3, 0., 2000.);
  f.SetParLimits(4, -40., 0.);
}

Bool_t SigmaParamsAtHardLimit(const Double_t p[5], TString& why)
{
  constexpr Double_t eps = 1e-3;
  if (!(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) && std::isfinite(p[3]) &&
        std::isfinite(p[4]))) {
    why = "nonfinite";
    return true;
  }
  if (p[0] <= eps) {
    why = "[0]~0";
    return true;
  }
  if (p[0] >= 200. - eps) {
    why = "[0]~200";
    return true;
  }
  if (p[3] >= 2000. - 1.) {
    why = "[3]~2000";
    return true;
  }
  if (p[4] <= -40. + eps) {
    why = "[4]~-40";
    return true;
  }
  if (p[4] >= -eps) {
    why = "[4]~0";
    return true;
  }
  return false;
}

Double_t SigmaFitRelMedianAbsDev(TGraphErrors& gr, const Double_t par[5])
{
  std::vector<Double_t> rel;
  rel.reserve(static_cast<std::size_t>(gr.GetN()));
  for (Int_t i = 0; i < gr.GetN(); ++i) {
    Double_t x = 0., y = 0.;
    gr.GetPoint(i, x, y);
    if (!(y > 1e-6) || !std::isfinite(y))
      continue;
    const Double_t f = CalcTPCdEdxSigma(par, x);
    if (!(f > 0.) || !std::isfinite(f))
      continue;
    rel.push_back(TMath::Abs(f / y - 1.));
  }
  return rel.empty() ? 1e9 : MedianSorted(rel);
}

void CopySigmaFromTF1(Double_t out[5], const TF1& f)
{
  for (Int_t i = 0; i < 5; ++i)
    out[i] = f.GetParameter(i);
}

// Staged fit: (1) median scale of E72 shape → (2) free 5-param → (3) retry with
// [1],[2] fixed then released. Reject hard-limit sticking or poor agreement vs
// RMS points; else E72×scale fallback.
Bool_t FitSigma5(TGraphErrors& gr, const Double_t init[5], Double_t out[5], Double_t& scale_fallback,
                 Bool_t& used_fallback, TString& detail)
{
  used_fallback = kFALSE;
  scale_fallback = 1.;
  CopySigma(out, init);
  detail = "fallback";
  if (gr.GetN() < 3) {
    used_fallback = kTRUE;
    detail = "fallback (Nbins<3)";
    return false;
  }

  // Stage1: global scale of E72 shape (median of y/sigma_E72).
  const Double_t sc = MedianScaleToGraph(gr, init);
  scale_fallback = sc;
  Double_t scaled[5];
  for (Int_t i = 0; i < 5; ++i)
    scaled[i] = sc * init[i];
  constexpr Double_t kMaxRelMed = 0.25; // reject 5-param if median |fit/y-1| exceeds this

  auto try_accept = [&](TF1& f, Int_t fit_status, const char* stage_tag) -> Bool_t {
    Double_t cand[5];
    CopySigmaFromTF1(cand, f);
    TString why;
    if (SigmaParamsAtHardLimit(cand, why)) {
      detail = Form("%s rejected (%s, fitStatus=%d)", stage_tag, why.Data(), fit_status);
      return false;
    }
    const Double_t rel = SigmaFitRelMedianAbsDev(gr, cand);
    if (!(rel <= kMaxRelMed)) {
      detail = Form("%s rejected (relMed=%.2f>%.2f, fitStatus=%d)", stage_tag, rel, kMaxRelMed, fit_status);
      return false;
    }
    CopySigma(out, cand);
    detail = Form("OK %s (scale0=%.3g, Nbins=%d, fitStatus=%d, relMed=%.3f)", stage_tag, sc, gr.GetN(),
                  fit_status, rel);
    return true;
  };

  TF1 f("f_sigma", "[0]+[1]*TMath::Abs(x)+[2]*x*x+[3]*TMath::Exp([4]*TMath::Abs(x))", kPLo, kPHi);

  // Stage2: free all 5 from scaled E72 start.
  f.SetParameters(scaled[0], scaled[1], scaled[2], scaled[3], scaled[4]);
  ApplySigmaSoftLimits(f);
  Int_t status = gr.Fit(&f, "QRMN");
  if (try_accept(f, status, "stage2"))
    return true;

  // Stage3: fit [0],[3],[4] with [1],[2] fixed, then release and refit.
  f.SetParameters(scaled[0], scaled[1], scaled[2], scaled[3], scaled[4]);
  ApplySigmaSoftLimits(f);
  f.FixParameter(1, scaled[1]);
  f.FixParameter(2, scaled[2]);
  status = gr.Fit(&f, "QRMN");
  f.ReleaseParameter(1);
  f.ReleaseParameter(2);
  ApplySigmaSoftLimits(f);
  status = gr.Fit(&f, "QRMN");
  if (try_accept(f, status, "stage3"))
    return true;

  // Fallback: keep E72 shape × median scale.
  for (Int_t i = 0; i < 5; ++i)
    out[i] = scaled[i];
  used_fallback = kTRUE;
  if (!detail.BeginsWith("stage") && !detail.Contains("rejected"))
    detail = Form("fallback scale=%.3g", sc);
  else
    detail = Form("fallback scale=%.3g after %s", sc, detail.Data());
  return false;
}

void FormatSigmaLine(std::ostream& os, const char* key, const Double_t s[5])
{
  os << key << "  " << s[0] << "  " << s[1] << "  " << s[2] << "  " << s[3] << "  " << s[4] << "\n";
}

} // namespace

int main(int argc, char** argv)
{
  std::vector<TString> inFiles;
  TString outPdf;
  TString outUser;
  Int_t run = -1;
  CutCfg cuts;
  Bool_t mass_lo_set = kFALSE;
  Bool_t mass_hi_set = kFALSE;
  Bool_t kp_mass_lo_set = kFALSE;
  Bool_t kp_mass_hi_set = kFALSE;

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "--out-user" && i + 1 < argc) {
      outUser = argv[++i];
    } else if (a == "-r" && i + 1 < argc) {
      run = TString(argv[++i]).Atoi();
    } else if (a == "--close-dist-max" && i + 1 < argc) {
      cuts.close_dist_max = TString(argv[++i]).Atof();
    } else if (a == "--diff-angle-max" && i + 1 < argc) {
      cuts.diff_angle_max = TString(argv[++i]).Atof();
    } else if (a == "--mass-lo" && i + 1 < argc) {
      cuts.lambda_mass_lo = TString(argv[++i]).Atof();
      mass_lo_set = kTRUE;
    } else if (a == "--mass-hi" && i + 1 < argc) {
      cuts.lambda_mass_hi = TString(argv[++i]).Atof();
      mass_hi_set = kTRUE;
    } else if (a == "--kp-mass-lo" && i + 1 < argc) {
      cuts.kp_mass_lo = TString(argv[++i]).Atof();
      kp_mass_lo_set = kTRUE;
    } else if (a == "--kp-mass-hi" && i + 1 < argc) {
      cuts.kp_mass_hi = TString(argv[++i]).Atof();
      kp_mass_hi_set = kTRUE;
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.Length() > 0 && a[0] == '-') {
      std::cerr << "Unknown option: " << a << "\n";
      usage(argv[0]);
      return 1;
    } else {
      inFiles.push_back(a);
    }
  }

  if (inFiles.empty()) {
    usage(argv[0]);
    return 1;
  }
  if (mass_lo_set != mass_hi_set) {
    std::cerr << "Error: set both --mass-lo and --mass-hi, or neither.\n";
    return 1;
  }
  if (kp_mass_lo_set != kp_mass_hi_set) {
    std::cerr << "Error: set both --kp-mass-lo and --kp-mass-hi, or neither.\n";
    return 1;
  }
  if (!kp_mass_lo_set) {
    const Double_t mk = KaonMass();
    cuts.kp_mass_lo = mk - kKpMassHalfWin;
    cuts.kp_mass_hi = mk + kKpMassHalfWin;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);

  TChain chain("dedx_pid_sample");
  for (const auto& f : inFiles) {
    const Int_t n = chain.Add(f);
    if (n <= 0)
      std::cerr << "Warning: failed to add " << f << "\n";
  }
  if (chain.GetNtrees() <= 0) {
    std::cerr << "Error: no input trees loaded.\n";
    return 1;
  }

  Int_t runOut = run;
  if (runOut < 0)
    runOut = 0;
  if (run < 0) {
    const Int_t runTree = ParseRunFromChain(chain);
    if (runTree > 0)
      runOut = runTree;
    else
      runOut = ParseRunFromPath(inFiles.front());
  }
  if (runOut < 0)
    runOut = 0;

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
  if (outPdf.IsNull())
    outPdf = Form("%s/tpc_dedx_pid_calib_run%05d.pdf", imgDir.Data(), runOut);

  if (outUser.IsNull()) {
    const TString userDir = Form("%s/tpc_dedx_pid", OUTPUT_DIR.Data());
    gSystem->mkdir(userDir.Data(), kTRUE);
    outUser = Form("%s/tpc_dedx_pid_user_run%05d.txt", userDir.Data(), runOut);
  }

  UInt_t run_number = 0;
  UInt_t event_number = 0;
  Int_t source = -1;
  Int_t species = -1;
  Int_t charge = 0;
  Double_t mom0 = 0.;
  Double_t dEdx = 0.;
  Int_t pid = 0;
  Int_t itag = -1;
  Int_t effective_ntTpc = -1;
  Double_t close_dist = 0.;
  Double_t delta_phi_scat = 0.;
  Double_t diff_angle = 0.;
  Double_t tag_mass = 0.;

  chain.SetBranchAddress("run_number", &run_number);
  chain.SetBranchAddress("event_number", &event_number);
  chain.SetBranchAddress("source", &source);
  chain.SetBranchAddress("species", &species);
  chain.SetBranchAddress("charge", &charge);
  chain.SetBranchAddress("mom0", &mom0);
  chain.SetBranchAddress("dEdx", &dEdx);
  chain.SetBranchAddress("pid", &pid);
  chain.SetBranchAddress("itag", &itag);
  chain.SetBranchAddress("effective_ntTpc", &effective_ntTpc);
  chain.SetBranchAddress("close_dist", &close_dist);
  chain.SetBranchAddress("delta_phi_scat", &delta_phi_scat);
  chain.SetBranchAddress("diff_angle", &diff_angle);
  chain.SetBranchAddress("tag_mass", &tag_mass);

  // Pre / post cut distributions.
  TH1D h_close_pre("h_close_pre", "close_dist (all);close_dist [mm];Entries", 200, 0., 50.);
  TH1D h_close_post("h_close_post", "close_dist (pass);close_dist [mm];Entries", 200, 0., 50.);
  TH1D h_dangle_kp_pre("h_dangle_kp_pre", "diff_angle kp (all);diff_angle [rad];Entries", 200, 0., 1.5);
  TH1D h_dangle_kp_post("h_dangle_kp_post", "diff_angle kp (pass);diff_angle [rad];Entries", 200, 0., 1.5);
  TH1D h_mass_kp_pre("h_mass_kp_pre", "tag_mass kp (all);tag_mass [GeV/c^{2}];Entries", 200, 0., 1.2);
  TH1D h_mass_kp_post("h_mass_kp_post", "tag_mass kp (pass);tag_mass [GeV/c^{2}];Entries", 200, 0., 1.2);
  TH1D h_mass_lam_pre("h_mass_lam_pre", "tag_mass lambda (all);tag_mass [GeV/c^{2}];Entries", 200, 1.0, 1.3);
  TH1D h_mass_lam_post("h_mass_lam_post", "tag_mass lambda (pass);tag_mass [GeV/c^{2}];Entries", 200, 1.0, 1.3);
  for (TH1D* h : {&h_close_pre, &h_close_post, &h_dangle_kp_pre, &h_dangle_kp_post, &h_mass_kp_pre, &h_mass_kp_post,
                   &h_mass_lam_pre, &h_mass_lam_post}) {
    h->SetDirectory(nullptr);
    h->Sumw2(kFALSE);
  }

  std::vector<Sample> passed;
  passed.reserve(100000);

  const Long64_t nent = chain.GetEntries();
  std::cout << "Reading " << nent << " entries from " << chain.GetNtrees() << " file(s)...\n";
  for (Long64_t ie = 0; ie < nent; ++ie) {
    displayProgressBar(static_cast<Int_t>(ie + 1), static_cast<Int_t>(nent));
    chain.GetEntry(ie);

    if (std::isfinite(close_dist))
      h_close_pre.Fill(close_dist);
    if (source == 0) {
      if (std::isfinite(diff_angle))
        h_dangle_kp_pre.Fill(diff_angle);
      if (std::isfinite(tag_mass))
        h_mass_kp_pre.Fill(tag_mass);
    } else if (source == 1) {
      if (std::isfinite(tag_mass))
        h_mass_lam_pre.Fill(tag_mass);
    }

    if (!PassCuts(source, effective_ntTpc, close_dist, diff_angle, tag_mass, cuts))
      continue;
    if (species < 0 || species >= kNSpecies)
      continue;
    if (!(std::isfinite(mom0) && std::isfinite(dEdx) && dEdx > 0. && TMath::Abs(mom0) > 1e-4))
      continue;

    if (std::isfinite(close_dist))
      h_close_post.Fill(close_dist);
    if (source == 0) {
      if (std::isfinite(diff_angle))
        h_dangle_kp_post.Fill(diff_angle);
      if (std::isfinite(tag_mass))
        h_mass_kp_post.Fill(tag_mass);
    } else if (source == 1) {
      if (std::isfinite(tag_mass))
        h_mass_lam_post.Fill(tag_mass);
    }

    Sample s;
    s.species = species;
    s.charge = charge;
    s.mom0 = mom0;
    s.dEdx = dEdx;
    passed.push_back(s);
  }

  std::cout << "Passed cuts: " << passed.size() << " / " << nent << "\n";

  // --- Conversion: proton median dEdx/Bethe(conv=1); apply to mean + sigma fit ---
  const Double_t mp_mev = 1000. * ProtonMass();
  std::vector<Double_t> proton_ratios;
  proton_ratios.reserve(passed.size());
  for (const auto& s : passed) {
    if (s.species != 2)
      continue;
    const Double_t bethe1 = HypTPCBethe(s.Poq(), mp_mev, 1.0);
    if (!(bethe1 > 0.) || !std::isfinite(bethe1))
      continue;
    proton_ratios.push_back(s.dEdx / bethe1);
  }
  Double_t conv_fit = kConvE72;
  if (proton_ratios.empty()) {
    std::cout << "TpcDedxConversionFactor fallback E72=" << conv_fit
              << " (no proton samples for median)\n";
  } else {
    conv_fit = MedianSorted(proton_ratios);
    std::cout << "TpcDedxConversionFactor=" << conv_fit
              << "  (proton median dEdx/Bethe_1, N=" << proton_ratios.size()
              << "; E72 was " << kConvE72 << ")\n";
  }

  // --- Per-species additive offset (median residual vs conv*Bethe); enables π/K/p together ---
  Double_t offset_fit[kNSpecies] = {0., 0., 0.};
  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    const Double_t mass_mev = 1000. * SpeciesMassGeV(sp);
    std::vector<Double_t> resid;
    resid.reserve(passed.size() / 3 + 8);
    for (const auto& s : passed) {
      if (s.species != sp)
        continue;
      const Double_t mean = HypTPCBethe(s.Poq(), mass_mev, conv_fit);
      if (!std::isfinite(mean))
        continue;
      resid.push_back(s.dEdx - mean);
    }
    if (!resid.empty()) {
      offset_fit[sp] = MedianSorted(resid);
      std::cout << "TpcDedxOffset" << kSpeciesName[sp] << "=" << offset_fit[sp]
                << "  (median residual, N=" << resid.size() << ")\n";
    }
  }

  // --- Per-species residual RMS vs |p| → sigma (means = conv*Bethe + offset) ---
  Double_t sigma_new[kNSpecies][5];
  Bool_t sigma_fallback[kNSpecies] = {kFALSE, kFALSE, kFALSE};
  Double_t sigma_scale[kNSpecies] = {1., 1., 1.};
  TGraphErrors gr_rms[kNSpecies];

  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    CopySigma(sigma_new[sp], SigmaE72(sp));
    const Double_t mass_mev = 1000. * SpeciesMassGeV(sp);
    std::vector<std::vector<Double_t>> bins(static_cast<std::size_t>(kNpBins));
    for (const auto& s : passed) {
      if (s.species != sp)
        continue;
      const Double_t p = s.AbsP();
      if (p < kSigmaFitPMin || p >= kPHi)
        continue;
      const Int_t ib = static_cast<Int_t>((p - kPLo) / (kPHi - kPLo) * kNpBins);
      if (ib < 0 || ib >= kNpBins)
        continue;
      const Double_t mean = HypTPCBethe(s.Poq(), mass_mev, conv_fit) + offset_fit[sp];
      if (!std::isfinite(mean))
        continue;
      bins[static_cast<std::size_t>(ib)].push_back(s.dEdx - mean);
    }

    std::vector<Double_t> xs, ys, exs, eys;
    for (Int_t ib = 0; ib < kNpBins; ++ib) {
      const auto& v = bins[static_cast<std::size_t>(ib)];
      if (static_cast<Int_t>(v.size()) < kMinBinEntries)
        continue;
      Double_t sum = 0., sum2 = 0.;
      for (Double_t r : v) {
        sum += r;
        sum2 += r * r;
      }
      const Double_t n = static_cast<Double_t>(v.size());
      const Double_t mean = sum / n;
      const Double_t rms = TMath::Sqrt(TMath::Max(0., sum2 / n - mean * mean));
      if (!(rms > 0.) || !std::isfinite(rms))
        continue;
      const Double_t pcen = kPLo + (static_cast<Double_t>(ib) + 0.5) * (kPHi - kPLo) / kNpBins;
      const Double_t dp = 0.5 * (kPHi - kPLo) / kNpBins;
      xs.push_back(pcen);
      ys.push_back(rms);
      exs.push_back(dp);
      eys.push_back(rms / TMath::Sqrt(2. * n));
    }
    gr_rms[sp] = TGraphErrors(static_cast<Int_t>(xs.size()), xs.data(), ys.data(), exs.data(), eys.data());
    gr_rms[sp].SetName(Form("gr_rms_%s", kSpeciesName[sp]));
    gr_rms[sp].SetTitle(Form("%s residual RMS vs |p|;|p| [GeV/c];RMS(dE/dx-Bethe)", kSpeciesName[sp]));
    gr_rms[sp].SetMarkerStyle(20);
    gr_rms[sp].SetMarkerSize(0.8);

    TString fit_detail;
    FitSigma5(gr_rms[sp], SigmaE72(sp), sigma_new[sp], sigma_scale[sp], sigma_fallback[sp], fit_detail);
    std::cout << "Sigma " << kSpeciesName[sp] << ": " << fit_detail << "\n";
    if (sigma_fallback[sp]) {
      std::cout << "  → using E72 shape * " << sigma_scale[sp] << " (sigma[0]: " << SigmaE72(sp)[0] << " → "
                << sigma_new[sp][0] << ")\n";
    }
  }

  // --- Species histos (kp+lambda pooled) ---
  TH2D* h2_dedx_all =
    new TH2D("h2_dedx_all", "dE/dx vs |p| (all tagged #pi/K/p);|p| [GeV/c];dE/dx (a.u.)", 80, kPLo, kPHi,
             kDedxBins, 0., kDedxMax);
  h2_dedx_all->SetDirectory(nullptr);
  h2_dedx_all->Sumw2(kFALSE);
  TH2D* h2_dedx_all_noe =
    new TH2D("h2_dedx_all_noe",
             "dE/dx vs |p| (HypTPCdEdxElectron veto);|p| [GeV/c];dE/dx (a.u.)", 80, kPLo, kPHi, kDedxBins, 0.,
             kDedxMax);
  h2_dedx_all_noe->SetDirectory(nullptr);
  h2_dedx_all_noe->Sumw2(kFALSE);
  Long64_t n_e_id_analyzer = 0;

  TH2D* h2_dedx[kNSpecies] = {};
  TH2D* h2_resid[kNSpecies] = {};
  TH2D* h2_ns_old[kNSpecies] = {};
  TH2D* h2_ns_new[kNSpecies] = {};
  TH1D* h1_ns_old[kNSpecies] = {};
  TH1D* h1_ns_new[kNSpecies] = {};
  Long64_t n_sp[kNSpecies] = {};

  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    const char* sn = kSpeciesName[sp];
    h2_dedx[sp] =
      new TH2D(Form("h2_dedx_%s", sn), Form("dE/dx vs |p| (%s);|p| [GeV/c];dE/dx (a.u.)", sn), 80, kPLo, kPHi,
               kDedxBins, 0., kDedxMax);
    h2_resid[sp] =
      new TH2D(Form("h2_resid_%s", sn), Form("residual vs |p| (%s);|p| [GeV/c];dE/dx - Bethe_new", sn), 80, kPLo,
               kPHi, 200, -100., 100.);
    h2_ns_old[sp] =
      new TH2D(Form("h2_ns_old_%s", sn), Form("n#sigma old (%s);|p| [GeV/c];n#sigma", sn), 80, kPLo, kPHi,
               kNsigmaBins, -kNsigmaMax, kNsigmaMax);
    h2_ns_new[sp] =
      new TH2D(Form("h2_ns_new_%s", sn), Form("n#sigma new (%s);|p| [GeV/c];n#sigma", sn), 80, kPLo, kPHi,
               kNsigmaBins, -kNsigmaMax, kNsigmaMax);
    h1_ns_old[sp] = new TH1D(Form("h1_ns_old_%s", sn), Form("n#sigma old 1D (%s);n#sigma;Entries", sn),
                             kNsigmaBins, -kNsigmaMax, kNsigmaMax);
    h1_ns_new[sp] = new TH1D(Form("h1_ns_new_%s", sn), Form("n#sigma new 1D (%s);n#sigma;Entries", sn),
                             kNsigmaBins, -kNsigmaMax, kNsigmaMax);
    for (TH1* h : {static_cast<TH1*>(h2_dedx[sp]), static_cast<TH1*>(h2_resid[sp]),
                   static_cast<TH1*>(h2_ns_old[sp]), static_cast<TH1*>(h2_ns_new[sp]),
                   static_cast<TH1*>(h1_ns_old[sp]), static_cast<TH1*>(h1_ns_new[sp])}) {
      h->SetDirectory(nullptr);
      h->Sumw2(kFALSE);
    }
  }

  for (const auto& s : passed) {
    const Int_t sp = s.species;
    const Double_t p = s.AbsP();
    const Double_t poq = s.Poq();
    const Double_t mass_mev = 1000. * SpeciesMassGeV(sp);
    const Double_t mean_new = HypTPCBethe(poq, mass_mev, conv_fit) + offset_fit[sp];
    const Double_t mean_old = HypTPCBethe(poq, mass_mev, kConvE72);
    const Double_t sig_old = CalcTPCdEdxSigma(SigmaE72(sp), poq);
    const Double_t sig_new = CalcTPCdEdxSigma(sigma_new[sp], poq);
    if (!(sig_old > 0.) || !(sig_new > 0.))
      continue;
    const Double_t ns_old = (s.dEdx - mean_old) / sig_old;
    const Double_t ns_new = (s.dEdx - mean_new) / sig_new;

    ++n_sp[sp];
    h2_dedx_all->Fill(p, s.dEdx);
    // Same track-level condition as Kinematics::HypTPCdEdxElectron (new conv/offset + new σ_π).
    if (IsHypTPCdEdxElectron(s.dEdx, poq, conv_fit, offset_fit[0], sigma_new[0])) {
      ++n_e_id_analyzer;
    } else {
      h2_dedx_all_noe->Fill(p, s.dEdx);
    }
    h2_dedx[sp]->Fill(p, s.dEdx);
    if (std::isfinite(mean_new))
      h2_resid[sp]->Fill(p, s.dEdx - mean_new);
    if (std::isfinite(ns_old)) {
      h2_ns_old[sp]->Fill(p, ns_old);
      h1_ns_old[sp]->Fill(ns_old);
    }
    if (std::isfinite(ns_new)) {
      h2_ns_new[sp]->Fill(p, ns_new);
      h1_ns_new[sp]->Fill(ns_new);
    }
  }

  std::cout << "HypTPCdEdxElectron veto refill: removed " << n_e_id_analyzer << " / "
            << h2_dedx_all->GetEntries() << "  (kept " << h2_dedx_all_noe->GetEntries() << ")\n";

  // --- UserParam snip ---
  std::ostringstream snip;
  snip << "# tpc_dedx_pid_calib recommended (Npass=" << passed.size()
       << ", Nproton_ratio=" << proton_ratios.size() << ")\n";
  snip << "# cuts: close_dist<" << cuts.close_dist_max;
  snip << "; kp: effective_ntTpc==2 && diff_angle<" << cuts.diff_angle_max
       << " && mmass in [" << cuts.kp_mass_lo << "," << cuts.kp_mass_hi << "]";
  snip << "; lambda: tag_mass in [" << cuts.lambda_mass_lo << "," << cuts.lambda_mass_hi << "]";
  snip << "; no delta_phi cut; PDF pooled by species\n";
  snip << "# conversion = proton median; offsets = per-species residual median; then sigma\n";
  snip << "# windows E72: pi [" << kWinPiE72[0] << "," << kWinPiE72[1] << "]  K [" << kWinKE72[0]
       << "," << kWinKE72[1] << "]  p [" << kWinPE72[0] << "," << kWinPE72[1]
       << "]  (candidate: pi -> [" << kWinPiTight[0] << "," << kWinPiTight[1] << "])\n";
  snip << "TpcDedxConversionFactor  " << conv_fit << "\n";
  snip << "TpcDedxOffsetPion  " << offset_fit[0] << "\n";
  snip << "TpcDedxOffsetKaon  " << offset_fit[1] << "\n";
  snip << "TpcDedxOffsetProton  " << offset_fit[2] << "\n";
  FormatSigmaLine(snip, "TpcDedxSigmaPion", sigma_new[0]);
  FormatSigmaLine(snip, "TpcDedxSigmaKaon", sigma_new[1]);
  FormatSigmaLine(snip, "TpcDedxSigmaProton", sigma_new[2]);
  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    if (sigma_fallback[sp])
      snip << "# note: Sigma" << kSpeciesName[sp] << " from E72 shape * " << sigma_scale[sp]
           << " (5-param fit failed)\n";
  }

  std::cout << "\n----- UserParam snip -----\n" << snip.str() << "--------------------------\n";
  {
    std::ofstream ofs(outUser.Data());
    if (!ofs) {
      std::cerr << "Warning: cannot write UserParam snip to " << outUser << "\n";
    } else {
      ofs << snip.str();
      std::cout << "Wrote UserParam snip to " << outUser << "\n";
    }
  }

  // --- PDF ---
  PdfWriter pdf(outPdf);
  TCanvas c("c_dedx_calib", "", 1200, 900);
  auto draw_prepost = [&](TH1D& pre, TH1D& post, const char* title) {
    c.Clear();
    StyleCanvas(c);
    c.cd();
    pre.SetLineColor(kBlack);
    post.SetLineColor(kRed + 1);
    pre.SetTitle(title);
    if (pre.GetMaximum() < post.GetMaximum())
      pre.SetMaximum(post.GetMaximum() * 1.15);
    pre.Draw("hist");
    post.Draw("hist same");
    auto* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(&pre, "pre-cut", "l");
    leg->AddEntry(&post, "post-cut", "l");
    leg->SetBit(kCanDelete);
    leg->Draw();
    c.Update();
    pdf.Print(c);
  };

  draw_prepost(h_close_pre, h_close_post, "close_dist pre/post;close_dist [mm];Entries");
  draw_prepost(h_dangle_kp_pre, h_dangle_kp_post, "diff_angle (kp) pre/post;diff_angle [rad];Entries");
  draw_prepost(h_mass_kp_pre, h_mass_kp_post, "tag_mass (kp) pre/post;tag_mass [GeV/c^{2}];Entries");
  draw_prepost(h_mass_lam_pre, h_mass_lam_post, "tag_mass (lambda) pre/post;tag_mass [GeV/c^{2}];Entries");

  // Sigma RMS fit pages
  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    if (gr_rms[sp].GetN() <= 0)
      continue;
    c.Clear();
    StyleCanvas(c);
    c.cd();
    gr_rms[sp].Draw("AP");
    TF1 fdraw(Form("fdraw_%d", sp), "[0]+[1]*TMath::Abs(x)+[2]*x*x+[3]*TMath::Exp([4]*TMath::Abs(x))", kPLo, kPHi);
    fdraw.SetParameters(sigma_new[sp]);
    fdraw.SetLineColor(kRed + 1);
    fdraw.Draw("same");
    TF1 fold(Form("fold_%d", sp), "[0]+[1]*TMath::Abs(x)+[2]*x*x+[3]*TMath::Exp([4]*TMath::Abs(x))", kPLo, kPHi);
    fold.SetParameters(SigmaE72(sp));
    fold.SetLineColor(kGray + 2);
    fold.SetLineStyle(2);
    fold.Draw("same");
    auto* tx = new TLatex();
    tx->SetNDC();
    tx->SetTextSize(0.03);
    tx->SetBit(kCanDelete);
    tx->DrawLatex(0.15, 0.92,
                  Form("%s sigma: solid=new  dashed=E72%s", kSpeciesName[sp],
                       sigma_fallback[sp] ? Form("  (fallback scale=%.3g)", sigma_scale[sp]) : ""));
    c.Update();
    pdf.Print(c);
  }

  Color_t bethe_col[kNSpecies] = {kRed, static_cast<Color_t>(kGreen + 2), kBlue};
  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    if (n_sp[sp] <= 0)
      continue;
    const Double_t mass_mev = 1000. * SpeciesMassGeV(sp);
    TGraph g_new = MakeBetheCurve(mass_mev, conv_fit, offset_fit[sp], bethe_col[sp]);
    const Double_t* win = WinE72(sp);
    // Dotted: current E72 window × E72 σ. Dashed: same window × fitted σ.
    TGraph g_lo = MakeBetheNsigmaEdge(mass_mev, conv_fit, offset_fit[sp], SigmaE72(sp), win[0], bethe_col[sp]);
    TGraph g_hi = MakeBetheNsigmaEdge(mass_mev, conv_fit, offset_fit[sp], SigmaE72(sp), win[1], bethe_col[sp]);
    TGraph g_lo_new =
      MakeBetheNsigmaEdge(mass_mev, conv_fit, offset_fit[sp], sigma_new[sp], win[0], bethe_col[sp]);
    TGraph g_hi_new =
      MakeBetheNsigmaEdge(mass_mev, conv_fit, offset_fit[sp], sigma_new[sp], win[1], bethe_col[sp]);
    g_lo_new.SetLineStyle(2);
    g_hi_new.SetLineStyle(2);

    auto draw2 = [&](TH2D* h, const char* note, Bool_t overlay_bethe, Bool_t logz) {
      c.Clear();
      StyleCanvas(c);
      c.cd();
      if (logz)
        gPad->SetLogz();
      else
        gPad->SetLogz(0);
      h->Draw("COLZ");
      if (overlay_bethe) {
        g_new.Draw("L same");
        g_lo.Draw("L same");
        g_hi.Draw("L same");
        g_lo_new.Draw("L same");
        g_hi_new.Draw("L same");
        auto* leg = new TLegend(0.48, 0.68, 0.88, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.022);
        leg->AddEntry(&g_new, Form("mean (conv+offset=%.3g)", offset_fit[sp]), "l");
        leg->AddEntry(&g_lo, Form("E72 #sigma window [%.2g,%.2g] dotted", win[0], win[1]), "l");
        leg->AddEntry(&g_lo_new, "fit #sigma same window dashed", "l");
        leg->SetBit(kCanDelete);
        leg->Draw();
      }
      auto* tx2 = new TLatex();
      tx2->SetNDC();
      tx2->SetTextSize(0.028);
      tx2->SetBit(kCanDelete);
      tx2->DrawLatex(0.12, 0.93, Form("%s  N=%lld", note, static_cast<long long>(n_sp[sp])));
      c.Update();
      pdf.Print(c);
      gPad->SetLogz(0);
    };

    draw2(h2_dedx[sp], Form("%s : dE/dx vs |p|", kSpeciesName[sp]), kTRUE, kTRUE);
    draw2(h2_resid[sp], Form("%s : residual vs |p|", kSpeciesName[sp]), kFALSE, kFALSE);
    draw2(h2_ns_old[sp], Form("%s : nsigma OLD", kSpeciesName[sp]), kFALSE, kFALSE);
    draw2(h2_ns_new[sp], Form("%s : nsigma NEW", kSpeciesName[sp]), kFALSE, kFALSE);

    c.Clear();
    StyleCanvas(c);
    c.SetRightMargin(0.08);
    c.cd();
    h1_ns_old[sp]->SetLineColor(kGray + 2);
    h1_ns_new[sp]->SetLineColor(bethe_col[sp]);
    h1_ns_old[sp]->SetTitle(Form("n#sigma 1D (%s);n#sigma;Entries", kSpeciesName[sp]));
    const Double_t ymax = 1.15 * TMath::Max(h1_ns_old[sp]->GetMaximum(), h1_ns_new[sp]->GetMaximum());
    h1_ns_old[sp]->SetMaximum(TMath::Max(ymax, 1.));
    h1_ns_old[sp]->Draw("hist");
    h1_ns_new[sp]->Draw("hist same");
    auto* leg = new TLegend(0.62, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h1_ns_old[sp],
                  Form("old mean=%.2f rms=%.2f", h1_ns_old[sp]->GetMean(), h1_ns_old[sp]->GetRMS()), "l");
    leg->AddEntry(h1_ns_new[sp],
                  Form("new mean=%.2f rms=%.2f", h1_ns_new[sp]->GetMean(), h1_ns_new[sp]->GetRMS()), "l");
    leg->SetBit(kCanDelete);
    leg->Draw();
    TLine z0(0., 0., 0., h1_ns_old[sp]->GetMaximum());
    z0.SetLineStyle(2);
    z0.Draw();
    c.Update();
    pdf.Print(c);
  }

  // Overview: Bethe (fitted conv) + E72-σ window (dotted) + fit-σ window (dashed)
  auto draw_overview = [&](TH2D* h, const char* title_note, Long64_t nshow) {
    c.Clear();
    StyleCanvas(c);
    c.cd();
    gPad->SetLogz();
    h->Draw("COLZ");
    const Double_t m_pi = 1000. * PionMass();
    const Double_t m_k = 1000. * KaonMass();
    const Double_t m_p = 1000. * ProtonMass();
    const Double_t m_e = 1000. * ElectronMass();
    TGraph g_pi = MakeBetheCurve(m_pi, conv_fit, offset_fit[0], kRed);
    TGraph g_k = MakeBetheCurve(m_k, conv_fit, offset_fit[1], static_cast<Color_t>(kGreen + 2));
    TGraph g_p = MakeBetheCurve(m_p, conv_fit, offset_fit[2], kBlue);
    TGraph g_e = MakeBetheCurve(m_e, conv_fit, 0., kMagenta + 1);
    g_e.SetLineStyle(3);
    g_e.SetLineWidth(2);
    TGraph g_pi_lo = MakeBetheNsigmaEdge(m_pi, conv_fit, offset_fit[0], kSigmaPiE72, kWinPiE72[0], kRed);
    TGraph g_pi_hi = MakeBetheNsigmaEdge(m_pi, conv_fit, offset_fit[0], kSigmaPiE72, kWinPiE72[1], kRed);
    TGraph g_k_lo =
      MakeBetheNsigmaEdge(m_k, conv_fit, offset_fit[1], kSigmaKE72, kWinKE72[0], static_cast<Color_t>(kGreen + 2));
    TGraph g_k_hi =
      MakeBetheNsigmaEdge(m_k, conv_fit, offset_fit[1], kSigmaKE72, kWinKE72[1], static_cast<Color_t>(kGreen + 2));
    TGraph g_p_lo = MakeBetheNsigmaEdge(m_p, conv_fit, offset_fit[2], kSigmaPE72, kWinPE72[0], kBlue);
    TGraph g_p_hi = MakeBetheNsigmaEdge(m_p, conv_fit, offset_fit[2], kSigmaPE72, kWinPE72[1], kBlue);
    TGraph g_pi_lo_n = MakeBetheNsigmaEdge(m_pi, conv_fit, offset_fit[0], sigma_new[0], kWinPiE72[0], kRed);
    TGraph g_pi_hi_n = MakeBetheNsigmaEdge(m_pi, conv_fit, offset_fit[0], sigma_new[0], kWinPiE72[1], kRed);
    TGraph g_k_lo_n =
      MakeBetheNsigmaEdge(m_k, conv_fit, offset_fit[1], sigma_new[1], kWinKE72[0], static_cast<Color_t>(kGreen + 2));
    TGraph g_k_hi_n =
      MakeBetheNsigmaEdge(m_k, conv_fit, offset_fit[1], sigma_new[1], kWinKE72[1], static_cast<Color_t>(kGreen + 2));
    TGraph g_p_lo_n = MakeBetheNsigmaEdge(m_p, conv_fit, offset_fit[2], sigma_new[2], kWinPE72[0], kBlue);
    TGraph g_p_hi_n = MakeBetheNsigmaEdge(m_p, conv_fit, offset_fit[2], sigma_new[2], kWinPE72[1], kBlue);
    for (TGraph* g : {&g_pi_lo_n, &g_pi_hi_n, &g_k_lo_n, &g_k_hi_n, &g_p_lo_n, &g_p_hi_n})
      g->SetLineStyle(2);
    g_pi.Draw("L same");
    g_k.Draw("L same");
    g_p.Draw("L same");
    g_e.Draw("L same");
    g_pi_lo.Draw("L same");
    g_pi_hi.Draw("L same");
    g_k_lo.Draw("L same");
    g_k_hi.Draw("L same");
    g_p_lo.Draw("L same");
    g_p_hi.Draw("L same");
    g_pi_lo_n.Draw("L same");
    g_pi_hi_n.Draw("L same");
    g_k_lo_n.Draw("L same");
    g_k_hi_n.Draw("L same");
    g_p_lo_n.Draw("L same");
    g_p_hi_n.Draw("L same");
    auto* leg = new TLegend(0.45, 0.55, 0.88, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.020);
    leg->AddEntry(&g_pi, Form("mean (conv=%.4g + offset)", conv_fit), "l");
    leg->AddEntry(&g_e, "e Bethe (ref)", "l");
    leg->AddEntry(&g_pi_lo, "E72 #sigma #times window (dotted)", "l");
    leg->AddEntry(&g_pi_lo_n, "fit #sigma #times window (dashed)", "l");
    leg->SetBit(kCanDelete);
    leg->Draw();
    auto* txo = new TLatex();
    txo->SetNDC();
    txo->SetTextSize(0.024);
    txo->SetBit(kCanDelete);
    txo->DrawLatex(0.12, 0.93, Form("%s  N=%lld", title_note, static_cast<long long>(nshow)));
    c.Update();
    pdf.Print(c);
    gPad->SetLogz(0);
  };

  draw_overview(h2_dedx_all, "all tagged #pi/K/p (log Z)", static_cast<Long64_t>(passed.size()));
  draw_overview(h2_dedx_all_noe,
                Form("HypTPCdEdxElectron veto refill: n#sigma_#pi<%.2g & |n#sigma_e|<%.2g & |p/q|<%.2g (removed %lld)",
                     kElectronNsigmaPiMax, kElectronNsigmaAbsMax, kElectronPoqAbsMax,
                     static_cast<long long>(n_e_id_analyzer)),
                static_cast<Long64_t>(h2_dedx_all_noe->GetEntries()));

  // Cuts / electron note page
  {
    c.Clear();
    c.cd();
    auto* txc = new TLatex();
    txc->SetNDC();
    txc->SetTextFont(82);
    txc->SetTextSize(0.028);
    Double_t yc = 0.92;
    auto line = [&](const char* s) {
      txc->DrawLatex(0.08, yc, s);
      yc -= 0.045;
    };
    line("tpc_dedx_pid_calib — selection & electron note");
    txc->SetTextSize(0.022);
    yc -= 0.02;
    line(Form("Sample cuts (applied here): close_dist < %.3g mm", cuts.close_dist_max));
    line(Form("  kp: effective_ntTpc==2 && diff_angle < %.3g && mmass in [%.4g, %.4g] (m_K #pm 0.1)",
              cuts.diff_angle_max, cuts.kp_mass_lo, cuts.kp_mass_hi));
    line(Form("  lambda: tag_mass in [%.4g, %.4g] GeV/c^{2}  (nom 1.115 #pm 0.1)", cuts.lambda_mass_lo,
              cuts.lambda_mass_hi));
    line("  no delta_phi cut; PDF/fit pool kp+lambda by species only");
    yc -= 0.02;
    line("Species in this tool: 0=#pi, 1=K, 2=p  (reaction-tagged; no e sample)");
    line("Electron is excluded from conversion/sigma fits (not a tagged species).");
    yc -= 0.02;
    line("Analyzer HypTPCdEdxElectron (Kinematics; track-level refill on overview page):");
    line(Form("  n#sigma_#pi < %.2g  &&  |n#sigma_e| < %.2g  &&  |p/q| < %.2g GeV/c",
              kElectronNsigmaPiMax, kElectronNsigmaAbsMax, kElectronPoqAbsMax));
    line(Form("  #sigma_e constant: e+ = %.4g,  e- = %.4g  (UserParam TpcDedxSigmaElectron+/-)",
              kSigmaDedxEpE72, kSigmaDedxEmE72));
    line(Form("  Overview noe hist: refill excluding that flag (removed %lld; |p|<0.1 only).",
              static_cast<long long>(n_e_id_analyzer)));
    line("  High-|p| MIP band is mostly #pi (e/#pi Bethe overlap) — not removed by this ID.");
    line("  This calib does not retune electron sigma / e ID cuts.");
    c.Update();
    pdf.Print(c);
  }

  // Summary page
  c.Clear();
  c.cd();
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextFont(82);
  tx->SetTextSize(0.025);
  Double_t y = 0.92;
  tx->DrawLatex(0.08, y, "tpc_dedx_pid_calib — recommended UserParam");
  y -= 0.05;
  tx->SetTextSize(0.022);
  tx->DrawLatex(0.08, y,
                Form("Npass=%zu  conv=%.4g (proton median; E72 was %.4g)  Nproton=%zu", passed.size(),
                     conv_fit, kConvE72, proton_ratios.size()));
  y -= 0.04;
  tx->DrawLatex(0.08, y,
                Form("cuts: close_dist<%.3g; kp Nt==2 && dA<%.3g; lam M=[%.3g,%.3g]", cuts.close_dist_max,
                     cuts.diff_angle_max, cuts.lambda_mass_lo, cuts.lambda_mass_hi));
  y -= 0.06;
  std::string line;
  std::istringstream iss(snip.str());
  while (std::getline(iss, line)) {
    if (y < 0.08)
      break;
    tx->DrawLatex(0.08, y, line.c_str());
    y -= 0.035;
  }
  c.Update();
  pdf.Print(c);
  pdf.Close(c);

  for (Int_t sp = 0; sp < kNSpecies; ++sp) {
    delete h2_dedx[sp];
    delete h2_resid[sp];
    delete h2_ns_old[sp];
    delete h2_ns_new[sp];
    delete h1_ns_old[sp];
    delete h1_ns_new[sp];
  }
  delete h2_dedx_all;
  delete h2_dedx_all_noe;

  std::cout << "Wrote " << pdf.n_printed << " page(s) to " << outPdf << "\n";
  return 0;
}
