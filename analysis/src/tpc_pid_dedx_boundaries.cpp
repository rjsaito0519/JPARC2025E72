// -*- C++ -*-
// Plot dE/dx vs signed p and 50 MeV/c dE/dx slices from DstTPCHelixTracking.
//
// PID curves are a manual copy of JPARC2025E72/src/Kinematics.cc (E72).
// Do not edit the analyzer to tune PDF overlays: change the copy / CLI here,
// then (later) port confirmed values into Kinematics.cc.
//
// Usage:
//   tpc_pid_dedx_boundaries <input.root> [-o <out.pdf>] [-r <run>]
//     [--conv <factor>]
//     [--pi-window lo,hi] [--p-window lo,hi] [--k-window lo,hi]
//     [--sep <power>]
//     [--sigma-pi a,b,c,d,e] [--sigma-k ...] [--sigma-p ...]
//
// Any --conv / window / sep / sigma flag turns on dual overlay
// (solid = E72 copy, dashed = trial). Default is E72 only.
//
// Input: TTree "tpc" (ntTpc, charge, mom0, dEdx); fallback TH2 PID_dEdx_vs_SignedMom.
// Output: multipage ROOT PDF — 2D (lin/log dE/dx, |qp|<=1) then 2x2 slices
// (|qp|<=1 GeV/c, 50 MeV/c, log Y, dE/dx to 250).

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPDGCode.h>
#include <TPad.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace
{

// Same binning as HistTools.cc BuildTPCHelixTracking (fill range).
constexpr Double_t kBinsSignedP[3] = {600., -3.0, 3.0};
constexpr Double_t kBinsDe[3] = {400., 0., 800.};
// Slice / display scan: |qp| <= 1 GeV/c is enough.
constexpr Double_t kScanQp[2] = {-1.0, 1.0};
constexpr Double_t kSliceDp = 0.05; // 50 MeV/c
constexpr Double_t kSliceDeMax = 250.;

// Copied from Kinematics.cc (E72). Keep in sync manually.
constexpr Double_t kConvE72 = 7388.11;
constexpr Double_t kSigmaPiE72[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
constexpr Double_t kSigmaKE72[5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
constexpr Double_t kSigmaPE72[5] = {12.9717, -8.43799, 3.10608, 166.494, -6.56123};
constexpr Double_t kWinPiE72[2] = {-3., 3.};
constexpr Double_t kWinPE72[2] = {-1.5, 6.};
constexpr Double_t kWinKE72[2] = {-2., 2.};
constexpr Double_t kSepE72 = 3.;

constexpr Double_t kQpEps = 1e-4;

void SanitizePdfTitlesSignedQ(TH1* h)
{
  if (!h)
    return;
  h->SetTitle("PID dE/dx (signed qp);qp [GeV/c];dE/dx (a.u.)");
  h->GetXaxis()->SetTitle("qp [GeV/c]");
  h->GetYaxis()->SetTitle("dE/dx (a.u.)");
}

Double_t PdgMassGeV(int pdg_code)
{
  TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pdg_code);
  return p ? p->Mass() : -1.;
}

Double_t PionMass() { return PdgMassGeV(kPiMinus); }
Double_t KaonMass() { return PdgMassGeV(kKMinus); }
Double_t ProtonMass() { return PdgMassGeV(kProton); }

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

struct PidParams {
  Double_t conv = kConvE72;
  Double_t sigma_pi[5] = {kSigmaPiE72[0], kSigmaPiE72[1], kSigmaPiE72[2], kSigmaPiE72[3], kSigmaPiE72[4]};
  Double_t sigma_k[5] = {kSigmaKE72[0], kSigmaKE72[1], kSigmaKE72[2], kSigmaKE72[3], kSigmaKE72[4]};
  Double_t sigma_p[5] = {kSigmaPE72[0], kSigmaPE72[1], kSigmaPE72[2], kSigmaPE72[3], kSigmaPE72[4]};
  Double_t win_pi[2] = {kWinPiE72[0], kWinPiE72[1]};
  Double_t win_p[2] = {kWinPE72[0], kWinPE72[1]};
  Double_t win_k[2] = {kWinKE72[0], kWinKE72[1]};
  Double_t sep_thr = kSepE72;

  Double_t MeanPi(Double_t poq) const { return HypTPCBethe(poq, 1000.0 * PionMass(), conv); }
  Double_t MeanK(Double_t poq) const { return HypTPCBethe(poq, 1000.0 * KaonMass(), conv); }
  Double_t MeanP(Double_t poq) const { return HypTPCBethe(poq, 1000.0 * ProtonMass(), conv); }
  Double_t SigmaPi(Double_t poq) const { return CalcTPCdEdxSigma(sigma_pi, poq); }
  Double_t SigmaK(Double_t poq) const { return CalcTPCdEdxSigma(sigma_k, poq); }
  Double_t SigmaP(Double_t poq) const { return CalcTPCdEdxSigma(sigma_p, poq); }

  Double_t SepCut(Double_t poq) const
  {
    const Double_t dedx_pi = MeanPi(poq);
    const Double_t dedx_p = MeanP(poq);
    const Double_t spi = SigmaPi(poq);
    const Double_t spp = SigmaP(poq);
    const Double_t avg_sigma = 0.5 * (spi + spp);
    const Double_t dedx_diff = TMath::Abs(dedx_pi - dedx_p);
    const Double_t separation_power = dedx_diff / avg_sigma;
    return dedx_pi < dedx_p ? dedx_pi + 0.5 * separation_power * spi
                            : dedx_p + 0.5 * separation_power * spp;
  }
};

bool ParseNDoubles(const char* s, Double_t* out, Int_t n)
{
  TString str(s);
  TObjArray* parts = str.Tokenize(",");
  if (!parts || parts->GetEntries() != n) {
    delete parts;
    return false;
  }
  for (Int_t i = 0; i < n; ++i) {
    auto* os = dynamic_cast<TObjString*>(parts->At(i));
    if (!os) {
      delete parts;
      return false;
    }
    TString t = os->GetString();
    t = t.Strip(TString::kBoth);
    out[i] = t.Atof();
  }
  delete parts;
  return true;
}

Int_t NSlice()
{
  return static_cast<Int_t>(TMath::Nint((kScanQp[1] - kScanQp[0]) / kSliceDp));
}

Int_t SliceIndex(Double_t qp)
{
  if (qp < kScanQp[0] || qp >= kScanQp[1])
    return -1;
  const Int_t i = static_cast<Int_t>((qp - kScanQp[0]) / kSliceDp);
  const Int_t n = NSlice();
  if (i < 0 || i >= n)
    return -1;
  return i;
}

Double_t SliceLo(Int_t i) { return kScanQp[0] + kSliceDp * static_cast<Double_t>(i); }
Double_t SliceHi(Int_t i) { return SliceLo(i) + kSliceDp; }
Double_t SliceCenter(Int_t i) { return 0.5 * (SliceLo(i) + SliceHi(i)); }

TH2* MakeHistFromBins(const char* name, const char* title, const Double_t bx[3], const Double_t by[3])
{
  return new TH2D(name, title, static_cast<Int_t>(bx[0]), bx[1], bx[2], static_cast<Int_t>(by[0]), by[1], by[2]);
}

TH1D* MakeSliceHist(Int_t i)
{
  auto* h = new TH1D(Form("h_dedx_slice_%d", i), "", static_cast<Int_t>(kBinsDe[0]), kBinsDe[1], kBinsDe[2]);
  h->SetDirectory(nullptr);
  h->Sumw2(kFALSE);
  return h;
}

void FillFromTree(TFile& f, TH2* h2, std::vector<TH1D*>& slices)
{
  TTree* t = nullptr;
  f.GetObject("tpc", t);
  if (!t) {
    std::cerr << "No TTree 'tpc' in " << f.GetName() << std::endl;
    return;
  }
  Int_t ntTpc = 0;
  std::vector<Int_t>* charge = nullptr;
  std::vector<Double_t>* mom0 = nullptr;
  std::vector<Double_t>* dEdx = nullptr;
  t->SetBranchAddress("ntTpc", &ntTpc);
  t->SetBranchAddress("charge", &charge);
  t->SetBranchAddress("mom0", &mom0);
  t->SetBranchAddress("dEdx", &dEdx);
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (!charge || !mom0 || !dEdx)
      continue;
    const Int_t nt = ntTpc;
    for (Int_t j = 0; j < nt; ++j) {
      if (j >= static_cast<Int_t>(charge->size()) || j >= static_cast<Int_t>(mom0->size()) ||
          j >= static_cast<Int_t>(dEdx->size()))
        break;
      const Double_t signed_p = static_cast<Double_t>((*charge)[j]) * (*mom0)[j];
      const Double_t de = (*dEdx)[j];
      h2->Fill(signed_p, de);
      const Int_t is = SliceIndex(signed_p);
      if (is >= 0)
        slices[static_cast<std::size_t>(is)]->Fill(de);
    }
  }
}

void FillSlicesFromTH2(TH2* h2, std::vector<TH1D*>& slices)
{
  if (!h2)
    return;
  const TAxis* xa = h2->GetXaxis();
  const TAxis* ya = h2->GetYaxis();
  for (Int_t ix = 1; ix <= h2->GetNbinsX(); ++ix) {
    const Double_t qp = xa->GetBinCenter(ix);
    const Int_t is = SliceIndex(qp);
    if (is < 0)
      continue;
    TH1D* hs = slices[static_cast<std::size_t>(is)];
    for (Int_t iy = 1; iy <= h2->GetNbinsY(); ++iy) {
      const Double_t w = h2->GetBinContent(ix, iy);
      if (w == 0.)
        continue;
      hs->Fill(ya->GetBinCenter(iy), w);
    }
  }
}

bool LoadData(const TString& path, TH2*& h2, std::vector<TH1D*>& slices)
{
  TFile f(path.Data(), "READ");
  if (!f.IsOpen() || f.IsZombie()) {
    std::cerr << "Cannot open " << path << std::endl;
    return false;
  }

  h2 = MakeHistFromBins("h_pid", "PID dE/dx (signed qp);qp [GeV/c];dE/dx (a.u.)", kBinsSignedP, kBinsDe);
  h2->SetDirectory(nullptr);
  slices.resize(static_cast<std::size_t>(NSlice()), nullptr);
  for (Int_t i = 0; i < NSlice(); ++i)
    slices[static_cast<std::size_t>(i)] = MakeSliceHist(i);

  TTree* t = nullptr;
  f.GetObject("tpc", t);
  if (t) {
    FillFromTree(f, h2, slices);
  } else {
    TObject* o = f.Get("PID_dEdx_vs_SignedMom");
    TH2* src = dynamic_cast<TH2*>(o);
    if (!src) {
      std::cerr << "No TTree 'tpc' and no TH2 PID_dEdx_vs_SignedMom in " << path << std::endl;
      delete h2;
      h2 = nullptr;
      for (auto* h : slices)
        delete h;
      slices.clear();
      return false;
    }
    TH2* hc = static_cast<TH2*>(src->Clone("h_pid"));
    hc->SetDirectory(nullptr);
    delete h2;
    h2 = hc;
    FillSlicesFromTH2(src, slices);
  }

  SanitizePdfTitlesSignedQ(h2);
  if (h2->GetEntries() <= 0.) {
    std::cerr << "No entries in " << path << std::endl;
    delete h2;
    h2 = nullptr;
    for (auto* h : slices)
      delete h;
    slices.clear();
    return false;
  }
  return true;
}

void usage(const char* a0)
{
  std::cerr << "Usage: " << a0 << " <input.root> [-o <out.pdf>] [-r <run>]\n"
            << "  [--conv <factor>]\n"
            << "  [--pi-window lo,hi] [--p-window lo,hi] [--k-window lo,hi]\n"
            << "  [--sep <power>]\n"
            << "  [--sigma-pi a,b,c,d,e] [--sigma-k ...] [--sigma-p ...]\n"
            << "  Default curves: E72 copy from Kinematics.cc (conv=" << kConvE72 << ").\n"
            << "  Any of the flags above enables dual overlay (solid=E72, dashed=trial).\n"
            << "  Reads TTree tpc (charge,mom0,dEdx) or TH2 PID_dEdx_vs_SignedMom.\n"
            << "  PDF: 2D (lin/log Y, |qp|<=1) then 50 MeV/c slices in 2x2\n"
            << "  (|qp|<=1 GeV/c, log counts, dE/dx xmax=250).\n"
            << "  OUTPUT_DIR/img/runXXXXX/ (run from -r or first run##### in path; else 0)\n";
}

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

void StylePidCanvas(TCanvas& c)
{
  c.SetLeftMargin(0.10);
  c.SetRightMargin(0.15);
  c.SetBottomMargin(0.10);
  c.SetTopMargin(0.08);
}

void ConfigureLogZForColz(TH2* h)
{
  if (!h)
    return;
  const Double_t zmax = h->GetMaximum();
  if (zmax <= 0.)
    return;
  Double_t zminPos = zmax;
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      const Double_t v = h->GetBinContent(ix, iy);
      if (v > 0. && v < zminPos)
        zminPos = v;
    }
  }
  if (zminPos >= zmax)
    zminPos = zmax * 1e-6;
  h->SetMinimum(TMath::Max(1e-18, zminPos * 0.35));
  h->SetMaximum(zmax * 1.08);
}

Double_t LogYAxisMin(TH2* h)
{
  if (!h)
    return 0.5;
  const TAxis* ya = h->GetYaxis();
  Double_t loBest = ya->GetXmax();
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      if (h->GetBinContent(ix, iy) <= 0.)
        continue;
      const Double_t e = ya->GetBinLowEdge(iy);
      if (e < loBest)
        loBest = e;
    }
  }
  if (loBest >= ya->GetXmax() - 1e-9)
    return TMath::Max(0.5, ya->GetXmin() + 1e-3);
  return TMath::Max(0.08, loBest * 0.35);
}

struct OverlayGraphs {
  TGraph pi, k, p, cut;
  TGraph pi_lo, pi_hi, k_lo, k_hi, p_lo, p_hi;
};

void StyleMean(TGraph& g, Color_t col, Style_t sty, Width_t w)
{
  g.SetLineColor(col);
  g.SetLineStyle(sty);
  g.SetLineWidth(w);
}

void StyleBand(TGraph& g, Color_t col, Style_t sty)
{
  g.SetLineColor(col);
  g.SetLineStyle(sty);
  g.SetLineWidth(1);
}

bool AppendPoint(std::vector<Double_t>& x, std::vector<Double_t>& y, Double_t poq, Double_t val)
{
  if (!std::isfinite(val) || TMath::Abs(poq) < kQpEps)
    return false;
  x.push_back(poq);
  y.push_back(val);
  return true;
}

OverlayGraphs BuildOverlay(const PidParams& par, Bool_t trial)
{
  OverlayGraphs og;
  const Int_t ngrid = 601;
  const Double_t pmin = kScanQp[0];
  const Double_t pmax = kScanQp[1];
  std::vector<Double_t> x_pi, y_pi, x_k, y_k, x_p, y_p, x_cut, y_cut;
  std::vector<Double_t> x_pilo, y_pilo, x_pihi, y_pihi;
  std::vector<Double_t> x_klo, y_klo, x_khi, y_khi;
  std::vector<Double_t> x_plo, y_plo, x_phi, y_phi;
  x_pi.reserve(ngrid);
  y_pi.reserve(ngrid);

  for (int i = 0; i < ngrid; ++i) {
    const Double_t t = static_cast<Double_t>(i) / static_cast<Double_t>(ngrid - 1);
    const Double_t poq = pmin + t * (pmax - pmin);
    const Double_t mpi = par.MeanPi(poq);
    const Double_t mk = par.MeanK(poq);
    const Double_t mp = par.MeanP(poq);
    const Double_t spi = par.SigmaPi(poq);
    const Double_t spk = par.SigmaK(poq);
    const Double_t spp = par.SigmaP(poq);
    AppendPoint(x_pi, y_pi, poq, mpi);
    AppendPoint(x_k, y_k, poq, mk);
    AppendPoint(x_p, y_p, poq, mp);
    AppendPoint(x_cut, y_cut, poq, par.SepCut(poq));
    AppendPoint(x_pilo, y_pilo, poq, mpi + par.win_pi[0] * spi);
    AppendPoint(x_pihi, y_pihi, poq, mpi + par.win_pi[1] * spi);
    AppendPoint(x_klo, y_klo, poq, mk + par.win_k[0] * spk);
    AppendPoint(x_khi, y_khi, poq, mk + par.win_k[1] * spk);
    AppendPoint(x_plo, y_plo, poq, mp + par.win_p[0] * spp);
    AppendPoint(x_phi, y_phi, poq, mp + par.win_p[1] * spp);
  }

  auto make = [](const std::vector<Double_t>& x, const std::vector<Double_t>& y) {
    return TGraph(static_cast<Int_t>(x.size()), x.data(), y.data());
  };
  og.pi = make(x_pi, y_pi);
  og.k = make(x_k, y_k);
  og.p = make(x_p, y_p);
  og.cut = make(x_cut, y_cut);
  og.pi_lo = make(x_pilo, y_pilo);
  og.pi_hi = make(x_pihi, y_pihi);
  og.k_lo = make(x_klo, y_klo);
  og.k_hi = make(x_khi, y_khi);
  og.p_lo = make(x_plo, y_plo);
  og.p_hi = make(x_phi, y_phi);

  const Style_t mean_sty = trial ? 2 : 1;
  const Style_t band_sty = trial ? 2 : 3;
  const Width_t mean_w = trial ? 2 : 2;
  StyleMean(og.pi, kRed, mean_sty, mean_w);
  StyleMean(og.k, kGreen + 2, mean_sty, mean_w);
  StyleMean(og.p, kBlue, mean_sty, mean_w);
  StyleMean(og.cut, kMagenta, trial ? 2 : 2, mean_w);
  StyleBand(og.pi_lo, kRed, band_sty);
  StyleBand(og.pi_hi, kRed, band_sty);
  StyleBand(og.k_lo, kGreen + 2, band_sty);
  StyleBand(og.k_hi, kGreen + 2, band_sty);
  StyleBand(og.p_lo, kBlue, band_sty);
  StyleBand(og.p_hi, kBlue, band_sty);
  return og;
}

void DrawOverlay(OverlayGraphs& og)
{
  og.pi.Draw("L same");
  og.k.Draw("L same");
  og.p.Draw("L same");
  og.cut.Draw("L same");
  og.pi_lo.Draw("L same");
  og.pi_hi.Draw("L same");
  og.k_lo.Draw("L same");
  og.k_hi.Draw("L same");
  og.p_lo.Draw("L same");
  og.p_hi.Draw("L same");
}

TLegend* BuildPidLegend(OverlayGraphs& cur, Bool_t has_trial)
{
  const Double_t y1 = has_trial ? 0.62 : 0.72;
  auto* leg = new TLegend(0.58, y1, 0.93, 0.93);
  leg->SetFillColorAlpha(kWhite, 0.92);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.017);
  leg->AddEntry(&cur.pi, "pi mean (E72)", "l");
  leg->AddEntry(&cur.k, "K mean (E72)", "l");
  leg->AddEntry(&cur.p, "p mean (E72)", "l");
  leg->AddEntry(&cur.cut, "pi/p separation cut", "l");
  leg->AddEntry(&cur.pi_lo, "pi nsigma window", "l");
  leg->AddEntry(&cur.k_lo, "K nsigma window", "l");
  leg->AddEntry(&cur.p_lo, "p nsigma window", "l");
  if (has_trial)
    leg->AddEntry(static_cast<TObject*>(nullptr), "dashed = trial params", "");
  return leg;
}

void DrawPidPageColz(TCanvas& cv, TH2* hist, Bool_t logy, OverlayGraphs& cur, OverlayGraphs* trial,
                     Bool_t has_trial, Double_t sep_e72, Double_t sep_trial)
{
  cv.Clear();
  StylePidCanvas(cv);
  cv.cd();
  ConfigureLogZForColz(hist);
  hist->Draw("COLZ");
  cv.SetLogx(kFALSE);
  cv.SetLogy(logy);
  cv.SetLogz(kTRUE);
  DrawOverlay(cur);
  if (has_trial && trial)
    DrawOverlay(*trial);
  TLegend* leg = BuildPidLegend(cur, has_trial);
  leg->SetBit(kCanDelete);
  leg->Draw();
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.018);
  tx->SetBit(kCanDelete);
  if (has_trial)
    tx->DrawLatex(0.10, 0.96, Form("sep threshold E72=%.2g  trial=%.2g (low sep: nsigma windows)", sep_e72, sep_trial));
  else
    tx->DrawLatex(0.10, 0.96, Form("E72 sep threshold = %.2g (low sep: nsigma windows)", sep_e72));
  cv.Update();
}

TLine* VLine(Double_t x, Double_t ymin, Double_t ymax, Color_t col, Style_t sty, Width_t w)
{
  if (!std::isfinite(x) || x <= 0. || x >= kSliceDeMax)
    return nullptr;
  auto* ln = new TLine(x, ymin, x, ymax);
  ln->SetLineColor(col);
  ln->SetLineStyle(sty);
  ln->SetLineWidth(w);
  ln->SetBit(kCanDelete);
  return ln;
}

void DrawSliceMarkers(Double_t poq, Double_t ymin, Double_t ymax, const PidParams& par, Bool_t trial)
{
  const Style_t mean_sty = trial ? 2 : 1;
  const Style_t win_sty = trial ? 2 : 3;
  const Width_t mean_w = 2;
  const Double_t mpi = par.MeanPi(poq);
  const Double_t mk = par.MeanK(poq);
  const Double_t mp = par.MeanP(poq);
  const Double_t spi = par.SigmaPi(poq);
  const Double_t spk = par.SigmaK(poq);
  const Double_t spp = par.SigmaP(poq);
  if (auto* ln = VLine(mpi, ymin, ymax, kRed, mean_sty, mean_w))
    ln->Draw();
  if (auto* ln = VLine(mk, ymin, ymax, kGreen + 2, mean_sty, mean_w))
    ln->Draw();
  if (auto* ln = VLine(mp, ymin, ymax, kBlue, mean_sty, mean_w))
    ln->Draw();
  if (auto* ln = VLine(mpi + par.win_pi[0] * spi, ymin, ymax, kRed, win_sty, 1))
    ln->Draw();
  if (auto* ln = VLine(mpi + par.win_pi[1] * spi, ymin, ymax, kRed, win_sty, 1))
    ln->Draw();
  if (auto* ln = VLine(mk + par.win_k[0] * spk, ymin, ymax, kGreen + 2, win_sty, 1))
    ln->Draw();
  if (auto* ln = VLine(mk + par.win_k[1] * spk, ymin, ymax, kGreen + 2, win_sty, 1))
    ln->Draw();
  if (auto* ln = VLine(mp + par.win_p[0] * spp, ymin, ymax, kBlue, win_sty, 1))
    ln->Draw();
  if (auto* ln = VLine(mp + par.win_p[1] * spp, ymin, ymax, kBlue, win_sty, 1))
    ln->Draw();
}

void DrawOneSlice(TPad* pad, TH1D* h, Int_t islice, const PidParams& e72, const PidParams& trial, Bool_t has_trial)
{
  pad->cd();
  pad->SetLeftMargin(0.14);
  pad->SetRightMargin(0.04);
  pad->SetBottomMargin(0.14);
  pad->SetTopMargin(0.12);
  pad->SetLogy(kTRUE);

  const Double_t lo = SliceLo(islice);
  const Double_t hi = SliceHi(islice);
  const Double_t qc = SliceCenter(islice);
  h->SetTitle(Form("qp [%.3f, %.3f)  N=%.0f;dE/dx (a.u.);Counts", lo, hi, h->GetEntries()));
  h->GetXaxis()->SetTitle("dE/dx (a.u.)");
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetRangeUser(0., kSliceDeMax);
  h->SetLineColor(kBlack);
  h->SetStats(kFALSE);

  Double_t ymax_vis = 1.;
  const Int_t ib0 = TMath::Max(1, h->FindBin(1e-6));
  const Int_t ib1 = TMath::Min(h->GetNbinsX(), h->FindBin(kSliceDeMax - 1e-6));
  for (Int_t ib = ib0; ib <= ib1; ++ib)
    ymax_vis = TMath::Max(ymax_vis, h->GetBinContent(ib));
  const Double_t ymin = 0.5;
  const Double_t ymax = TMath::Max(ymin * 10., ymax_vis) * 1.20;
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->Draw("hist");
  DrawSliceMarkers(qc, ymin, ymax, e72, kFALSE);
  if (has_trial)
    DrawSliceMarkers(qc, ymin, ymax, trial, kTRUE);
  auto* tx = new TLatex();
  tx->SetNDC();
  tx->SetTextSize(0.035);
  tx->SetBit(kCanDelete);
  if (has_trial)
    tx->DrawLatex(0.16, 0.93, "solid=E72  dashed=trial  red=pi  green=K  blue=p");
  else
    tx->DrawLatex(0.16, 0.93, "red=pi  green=K  blue=p  (E72)");
}

void PrintCanvas(TCanvas& c, const TString& path, Int_t page, Int_t ntotal)
{
  if (ntotal <= 1) {
    c.Print(path);
    return;
  }
  if (page == 1)
    c.Print(path + "(");
  else if (page == ntotal)
    c.Print(path + ")");
  else
    c.Print(path);
}

void PrintPidParams(const char* tag, const PidParams& p)
{
  std::cerr << tag << ": conv=" << p.conv
            << " pi_win=[" << p.win_pi[0] << "," << p.win_pi[1] << "]"
            << " p_win=[" << p.win_p[0] << "," << p.win_p[1] << "]"
            << " k_win=[" << p.win_k[0] << "," << p.win_k[1] << "]"
            << " sep=" << p.sep_thr << "\n";
  std::cerr << "  sigma_pi=" << p.sigma_pi[0] << "," << p.sigma_pi[1] << "," << p.sigma_pi[2] << ","
            << p.sigma_pi[3] << "," << p.sigma_pi[4] << "\n";
  std::cerr << "  sigma_k=" << p.sigma_k[0] << "," << p.sigma_k[1] << "," << p.sigma_k[2] << ","
            << p.sigma_k[3] << "," << p.sigma_k[4] << "\n";
  std::cerr << "  sigma_p=" << p.sigma_p[0] << "," << p.sigma_p[1] << "," << p.sigma_p[2] << ","
            << p.sigma_p[3] << "," << p.sigma_p[4] << "\n";
}

} // namespace

int main(int argc, char** argv)
{
  Int_t run = -1;
  TString inPath;
  TString outPdf;
  PidParams e72;
  PidParams trial = e72;
  Bool_t has_trial = kFALSE;

  auto mark_trial = [&]() { has_trial = kTRUE; };

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-r" && i + 1 < argc) {
      run = TString(argv[++i]).Atoi();
    } else if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "--conv" && i + 1 < argc) {
      trial.conv = TString(argv[++i]).Atof();
      mark_trial();
    } else if (a == "--pi-window" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.win_pi, 2)) {
        std::cerr << "Bad --pi-window (want lo,hi)\n";
        return 1;
      }
      mark_trial();
    } else if (a == "--p-window" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.win_p, 2)) {
        std::cerr << "Bad --p-window (want lo,hi)\n";
        return 1;
      }
      mark_trial();
    } else if (a == "--k-window" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.win_k, 2)) {
        std::cerr << "Bad --k-window (want lo,hi)\n";
        return 1;
      }
      mark_trial();
    } else if (a == "--sep" && i + 1 < argc) {
      trial.sep_thr = TString(argv[++i]).Atof();
      mark_trial();
    } else if (a == "--sigma-pi" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.sigma_pi, 5)) {
        std::cerr << "Bad --sigma-pi (want 5 comma-separated values)\n";
        return 1;
      }
      mark_trial();
    } else if (a == "--sigma-k" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.sigma_k, 5)) {
        std::cerr << "Bad --sigma-k (want 5 comma-separated values)\n";
        return 1;
      }
      mark_trial();
    } else if (a == "--sigma-p" && i + 1 < argc) {
      if (!ParseNDoubles(argv[++i], trial.sigma_p, 5)) {
        std::cerr << "Bad --sigma-p (want 5 comma-separated values)\n";
        return 1;
      }
      mark_trial();
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
    runOut = ParseRunFromPath(inPath);
  if (runOut < 0)
    runOut = 0;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  PrintPidParams("E72", e72);
  if (has_trial)
    PrintPidParams("trial", trial);

  TH2* h2 = nullptr;
  std::vector<TH1D*> slices;
  if (!LoadData(inPath, h2, slices))
    return 1;

  OverlayGraphs g_e72 = BuildOverlay(e72, kFALSE);
  OverlayGraphs g_trial;
  if (has_trial)
    g_trial = BuildOverlay(trial, kTRUE);

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
  if (outPdf.IsNull())
    outPdf = Form("%s/PID_dEdx_vs_SignedMom_run%05d.pdf", imgDir.Data(), runOut);

  std::vector<Int_t> nonempty;
  nonempty.reserve(slices.size());
  for (Int_t i = 0; i < static_cast<Int_t>(slices.size()); ++i) {
    if (slices[static_cast<std::size_t>(i)]->GetEntries() > 0.)
      nonempty.push_back(i);
  }
  const Int_t nslice_pages = nonempty.empty() ? 0 : (static_cast<Int_t>(nonempty.size()) + 3) / 4;
  const Int_t ntotal = 2 + nslice_pages;

  TH2* h_p1 = static_cast<TH2*>(h2->Clone("h_pid_p1"));
  h_p1->SetDirectory(nullptr);
  SanitizePdfTitlesSignedQ(h_p1);
  h_p1->GetXaxis()->SetRangeUser(kScanQp[0], kScanQp[1]);
  if (has_trial) {
    TString t(h_p1->GetTitle());
    const Ssiz_t sep = t.First(';');
    const TString tag = " [solid=E72, dashed=trial]";
    if (sep >= 0)
      t.Insert(sep, tag);
    else
      t += tag;
    h_p1->SetTitle(t);
  }

  TH2* h_p2 = static_cast<TH2*>(h2->Clone("h_pid_p2"));
  h_p2->SetDirectory(nullptr);
  SanitizePdfTitlesSignedQ(h_p2);
  h_p2->SetTitle("PID dE/dx (signed qp log Y);qp [GeV/c];dE/dx (a.u.)");
  h_p2->GetXaxis()->SetRangeUser(kScanQp[0], kScanQp[1]);
  if (has_trial) {
    TString t(h_p2->GetTitle());
    const Ssiz_t sep = t.First(';');
    const TString tag = " [solid=E72, dashed=trial]";
    if (sep >= 0)
      t.Insert(sep, tag);
    else
      t += tag;
    h_p2->SetTitle(t);
  }
  const Double_t yLogMax = kBinsDe[2];
  const Double_t yLogMin = TMath::Max(1e-6, LogYAxisMin(h_p2));
  h_p2->GetYaxis()->SetRangeUser(yLogMin, yLogMax);

  constexpr Int_t kCanvasW = 1400;
  constexpr Int_t kCanvasH = 1000;
  if (!outPdf.IsNull())
    std::remove(outPdf.Data());

  Int_t page = 0;
  TCanvas c2d("c_pid_2d", "", kCanvasW, kCanvasH);
  DrawPidPageColz(c2d, h_p1, kFALSE, g_e72, has_trial ? &g_trial : nullptr, has_trial, e72.sep_thr, trial.sep_thr);
  PrintCanvas(c2d, outPdf, ++page, ntotal);

  DrawPidPageColz(c2d, h_p2, kTRUE, g_e72, has_trial ? &g_trial : nullptr, has_trial, e72.sep_thr, trial.sep_thr);
  PrintCanvas(c2d, outPdf, ++page, ntotal);

  TCanvas csl("c_pid_slices", "", kCanvasW, kCanvasH);
  for (Int_t ip = 0; ip < nslice_pages; ++ip) {
    csl.Clear();
    csl.Divide(2, 2, 0.002, 0.002);
    for (Int_t k = 0; k < 4; ++k) {
      const Int_t idx = ip * 4 + k;
      if (idx >= static_cast<Int_t>(nonempty.size()))
        break;
      const Int_t islice = nonempty[static_cast<std::size_t>(idx)];
      DrawOneSlice(static_cast<TPad*>(csl.cd(k + 1)), slices[static_cast<std::size_t>(islice)], islice, e72, trial,
                   has_trial);
    }
    csl.Update();
    PrintCanvas(csl, outPdf, ++page, ntotal);
  }

  delete h_p2;
  delete h_p1;
  delete h2;
  for (auto* h : slices)
    delete h;

  std::cout << "Wrote " << ntotal << " page(s) to " << outPdf << "\n";
  return 0;
}
