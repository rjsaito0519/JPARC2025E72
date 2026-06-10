// -*- C++ -*-
// Plot PID_dEdx_vs_SignedMom with HypTPC dE/dx boundaries (same formulas as analyzer
// JPARC2025E72/src/Kinematics.cc — copied here; keep in sync manually).
//
// HypTPCdEdxPID logic (bit0=pi, bit1=K, bit2=p; non-exclusive in low-separation):
// - Electron: HypTPCdEdxElectron true -> return 0 (no bits).
// - If separation_power = |dedx_pi-dedx_p|/avg_sigma < 6: pi if nsigma_pi in (-3,3);
//   p if nsigma_p > -4 (upper 6 in array unused in code).
// - Else: pi if dedx < ppi_separation_cut, p if dedx > cut
//   (ppi_separation_cut = min(dedx_pi,dedx_p)+0.5*separation_power*sigma on that side).
// - K bit if nsigma_K in (-3,3), added on top of pi/p.
// PDF: single TCanvas pad — TH2 Draw("COLZ") then TGraph Draw("L same") + legend; SaveAs PNG then
// img2pdf (else ROOT Print). Page 1: linear dE/dx + log Z; page 2: log dE/dx + log Z; both signed qp.
//
// Usage: tpc_pid_dedx_boundaries <input.root> [-o <out.pdf>] [-r <run>]
//   [--pi-baseline] [--p-baseline] [--p-align-sigma <n>]
//   Default: pi mean = Bethe_pi - 3*sigma_pi; p mean = Bethe_p - 1*sigma_p (modest).
//   --pi-baseline / --p-baseline: E42 Bethe centers for that species.
//   ROOT: TH2 "PID_dEdx_vs_SignedMom" or TTree "tpc" (ntTpc, charge, mom0, dEdx)
//   Plot-only; does not change Kinematics.cc / DST pid

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH2.h>
#include <TH2D.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1.h>
#include <TTree.h>

#include <sys/wait.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace
{

/** Acrobat / some viewers choke on TLatex '#' in PDF text from ROOT. Use plain titles. */
void SanitizePdfTitlesSignedQ(TH1* h)
{
  if (!h)
    return;
  h->SetTitle("PID dE/dx (signed qp);qp [GeV/c];dE/dx (a.u.)");
  h->GetXaxis()->SetTitle("qp [GeV/c]");
  h->GetYaxis()->SetTitle("dE/dx (a.u.)");
}

// --- Same binning as HistTools.cc BuildTPCHelixTracking ---
constexpr Double_t kBinsSignedP[3] = {600., -3.0, 3.0};
constexpr Double_t kBinsDe[3] = {400., 0., 800.};

// Mass helpers: same convention as Kinematics.cc (pdg::ProtonMass() [GeV/c^2], Bethe: 1000.0 * …)
constexpr int kPdgDeuteron = 1000010020;

Double_t PdgMassGeV(int pdg_code)
{
  if (pdg_code == kPdgDeuteron) return 1.875613;
  TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(pdg_code);
  return p ? p->Mass() : -1.;
}

Double_t PionMass() { return PdgMassGeV(kPiMinus); }
Double_t KaonMass() { return PdgMassGeV(kKMinus); }
Double_t ProtonMass() { return PdgMassGeV(kProton); }

constexpr Double_t conversion_factor = 12171.3;
constexpr Double_t sigma_dedx_pi[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
constexpr Double_t sigma_dedx_k[5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
constexpr Double_t sigma_dedx_p[5] = {12.9717, -8.43799, 3.10608, 166.494, -6.56123};

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

// P10 only (materialid==0), same as Kinematics::HypTPCdEdx for gas
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

Double_t HypTPCBethe(Double_t poq /* GeV/c */, Double_t mass_mev)
{
  Double_t x = poq;
  Double_t p[2] = {conversion_factor, mass_mev};
  const Double_t momentum = 1000. * TMath::Abs(x);
  const Double_t energy = TMath::Hypot(p[1], momentum);
  const Double_t beta = BetaFromEandP(energy, momentum);
  return p[0] * HypTPCdEdxP10(p[1], beta);
}

Double_t CalcTPCdEdxSigma(const Double_t sigma_par[5], Double_t poq)
{
  const Double_t abspoq = TMath::Abs(poq);
  return sigma_par[0] + sigma_par[1] * abspoq + sigma_par[2] * poq * poq +
         sigma_par[3] * TMath::Exp(sigma_par[4] * abspoq);
}

Double_t HypTPCdEdxPionCurve(Double_t poq) { return HypTPCBethe(poq, 1000.0 * PionMass()); }

Double_t HypTPCdEdxKaonCurve(Double_t poq) { return HypTPCBethe(poq, 1000.0 * KaonMass()); }

Double_t HypTPCdEdxProtonCurve(Double_t poq) { return HypTPCBethe(poq, 1000.0 * ProtonMass()); }

/** Plot-only pi alignment: put legacy -3sigma band on the mean curve. */
struct PiAlign {
  bool center_at_legacy_minus3sigma = true;

  bool IsAligned() const { return center_at_legacy_minus3sigma; }

  TString Summary() const
  {
    return center_at_legacy_minus3sigma ? "pi center = legacy -3sigma (Bethe - 3*sigma_pi)"
                                        : "pi E42 baseline (Bethe center)";
  }

  void LegacyPiExpectation(Double_t poq, Double_t& dedx_bethe, Double_t& sigma_pi) const
  {
    dedx_bethe = HypTPCdEdxPionCurve(poq);
    sigma_pi = CalcTPCdEdxSigma(sigma_dedx_pi, poq);
  }

  Double_t PiCenter(Double_t poq) const
  {
    Double_t dedx_bethe = 0.;
    Double_t sigma_pi = 0.;
    LegacyPiExpectation(poq, dedx_bethe, sigma_pi);
    return center_at_legacy_minus3sigma ? dedx_bethe - 3. * sigma_pi : dedx_bethe;
  }
};

void PrintPiAlignReference(const PiAlign& align)
{
  if (!align.IsAligned())
    return;
  std::cerr << "Pi alignment (plot-only): center(p) = Bethe_pi(p) - 3*sigma_pi(p); "
            << "new +/-3sigma band around that center (Kinematics.cc unchanged)\n";
  for (const double absp : {0.3, 0.5, 1.0, 1.5}) {
    for (const double sign : {+1., -1.}) {
      const double poq = sign * absp;
      Double_t dedx_bethe = 0.;
      Double_t sigma_pi = 0.;
      align.LegacyPiExpectation(poq, dedx_bethe, sigma_pi);
      const double center = dedx_bethe - 3. * sigma_pi;
      std::cerr << "  qp=" << poq << " GeV/c: Bethe=" << dedx_bethe << " sigma=" << sigma_pi
                << " legacy_-3sigma=" << center << " (new center)\n";
    }
  }
}

/** Plot-only proton alignment: modest downward shift vs pi (default 1*sigma_p). */
struct ProtonAlign {
  bool enabled = true;
  Double_t sigma_shift = 1.0;

  bool IsAligned() const { return enabled && sigma_shift > 1e-12; }

  TString Summary() const
  {
    if (!IsAligned())
      return "p E42 baseline (Bethe center)";
    return Form("p center = Bethe - %.2g*sigma_p", sigma_shift);
  }

  void LegacyProtonExpectation(Double_t poq, Double_t& dedx_bethe, Double_t& sigma_p) const
  {
    dedx_bethe = HypTPCdEdxProtonCurve(poq);
    sigma_p = CalcTPCdEdxSigma(sigma_dedx_p, poq);
  }

  Double_t ProtonCenter(Double_t poq) const
  {
    Double_t dedx_bethe = 0.;
    Double_t sigma_p = 0.;
    LegacyProtonExpectation(poq, dedx_bethe, sigma_p);
    return IsAligned() ? dedx_bethe - sigma_shift * sigma_p : dedx_bethe;
  }
};

void PrintProtonAlignReference(const ProtonAlign& align)
{
  if (!align.IsAligned())
    return;
  std::cerr << "Proton alignment (plot-only): center(p) = Bethe_p(p) - " << align.sigma_shift
            << "*sigma_p(p); bands still -4/+6 sigma around center\n";
  for (const double absp : {0.3, 0.5, 1.0, 1.5}) {
    for (const double sign : {+1., -1.}) {
      const double poq = sign * absp;
      Double_t dedx_bethe = 0.;
      Double_t sigma_p = 0.;
      align.LegacyProtonExpectation(poq, dedx_bethe, sigma_p);
      const double center = dedx_bethe - align.sigma_shift * sigma_p;
      std::cerr << "  qp=" << poq << " GeV/c: Bethe=" << dedx_bethe << " sigma=" << sigma_p
                << " new_center=" << center << "\n";
    }
  }
}

Double_t PpiSeparationCut(Double_t poq)
{
  const Double_t dedx_pi = HypTPCdEdxPionCurve(poq);
  const Double_t dedx_p = HypTPCdEdxProtonCurve(poq);
  const Double_t sigma_pi = CalcTPCdEdxSigma(sigma_dedx_pi, poq);
  const Double_t sigma_p = CalcTPCdEdxSigma(sigma_dedx_p, poq);
  const Double_t avg_sigma = 0.5 * (sigma_pi + sigma_p);
  const Double_t dedx_diff = TMath::Abs(dedx_pi - dedx_p);
  const Double_t separation_power = dedx_diff / avg_sigma;
  return dedx_pi < dedx_p ? dedx_pi + 0.5 * separation_power * sigma_pi
                         : dedx_p + 0.5 * separation_power * sigma_p;
}

TH2* MakeHistFromBins(const char* name, const char* title, const Double_t bx[3], const Double_t by[3])
{
  return new TH2D(name, title, static_cast<Int_t>(bx[0]), bx[1], bx[2], static_cast<Int_t>(by[0]), by[1], by[2]);
}

void FillFromTree(TFile& f, TH2* h)
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
      h->Fill(signed_p, (*dEdx)[j]);
    }
  }
}

TH2* LoadOrBuildHistogram(const TString& path)
{
  TFile f(path.Data(), "READ");
  if (!f.IsOpen() || f.IsZombie()) {
    std::cerr << "Cannot open " << path << std::endl;
    return nullptr;
  }
  TObject* o = f.Get("PID_dEdx_vs_SignedMom");
  if (o) {
    TH2* h2 = dynamic_cast<TH2*>(o);
    if (h2) {
      TH2* hc = static_cast<TH2*>(h2->Clone("h_pid"));
      hc->SetDirectory(nullptr);
      SanitizePdfTitlesSignedQ(hc);
      return hc;
    }
  }
  TH2* h = MakeHistFromBins("h_pid", "PID dE/dx (signed qp);qp [GeV/c];dE/dx (a.u.)", kBinsSignedP, kBinsDe);
  FillFromTree(f, h);
  if (h->GetEntries() <= 0.) {
    delete h;
    return nullptr;
  }
  SanitizePdfTitlesSignedQ(h);
  return h;
}

void usage(const char* a0)
{
  std::cerr << "Usage: " << a0 << " <input.root> [-o <out.pdf>] [-r <run>]\n"
            << "  [--pi-baseline] [--p-baseline] [--p-align-sigma <n>]\n"
            << "  Default: pi center = Bethe_pi - 3*sigma_pi; p center = Bethe_p - 1*sigma_p.\n"
            << "  --pi-baseline / --p-baseline: E42 Bethe center for that species.\n"
            << "  --p-align-sigma: proton downward shift in units of sigma_p (default 1).\n"
            << "  Reads TH2 PID_dEdx_vs_SignedMom or TTree tpc (charge,mom0,dEdx).\n"
            << "  Output: 2-page PDF (PNG->img2pdf when available). COLZ then curves Draw(\"L same\") on one pad per page.\n"
            << "  OUTPUT_DIR/img/runXXXXX/ (run from -r or first run##### in path; else 0)\n";
}

void AppendAlignSubtitle(TH1* h, const PiAlign& pi_align, const ProtonAlign& p_align)
{
  if (!h || (!pi_align.IsAligned() && !p_align.IsAligned()))
    return;
  TString tag;
  if (pi_align.IsAligned())
    tag = pi_align.Summary();
  if (p_align.IsAligned()) {
    if (!tag.IsNull())
      tag += "; ";
    tag += p_align.Summary();
  }
  TString t(h->GetTitle());
  const Ssiz_t sep = t.First(';');
  const TString bracket = Form(" [%s]", tag.Data());
  if (sep >= 0)
    t.Insert(sep, bracket);
  else
    t += bracket;
  h->SetTitle(t);
}

/** First run followed by digits in path/filename (e.g. .../run01234/...), else 0. */
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

void DrawBoundaryCurves(TGraph& gr_pi, TGraph& gr_k, TGraph& gr_p, TGraph& gr_cut, TGraph& gr_pilo,
                        TGraph& gr_pihi, TGraph& gr_klo, TGraph& gr_khi, TGraph& gr_plo, TGraph& gr_phi)
{
  TGraph* gs[] = {&gr_pi,  &gr_k,   &gr_p,   &gr_cut, &gr_pilo,
                  &gr_pihi, &gr_klo, &gr_khi, &gr_plo, &gr_phi};
  for (auto* g : gs)
    g->Draw("L same");
}

void ConfigureLogZForColz(TH2* h);

/** Upper-right legend (NDC). Heap-allocated so it survives until PDF flush (stack TLegend vanishes). */
TLegend* BuildPidLegendUpperRight(TGraph& gr_pi, TGraph& gr_k, TGraph& gr_p, TGraph& gr_cut, TGraph& gr_pilo,
                                  TGraph& gr_klo, TGraph& gr_plo, const PiAlign& pi_align, const ProtonAlign& p_align)
{
  const Bool_t extra = pi_align.IsAligned() || p_align.IsAligned();
  const Double_t y1 = extra ? 0.66 : 0.78;
  auto* leg = new TLegend(0.62, y1, 0.93, 0.93);
  leg->SetFillColorAlpha(kWhite, 0.92);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.017);
  leg->SetEntrySeparation(0.20);
  // ASCII only (Acrobat-safe with vector PDF)
  leg->AddEntry(&gr_pi, pi_align.IsAligned() ? "pi mean (legacy -3sigma)" : "pi mean (Bethe)", "l");
  leg->AddEntry(&gr_k, "K mean (Bethe)", "l");
  leg->AddEntry(&gr_p, p_align.IsAligned() ? "p mean (aligned)" : "p mean (Bethe)", "l");
  leg->AddEntry(&gr_cut, "pi/p separation cut", "l");
  leg->AddEntry(&gr_pilo, "pi +/- 3 sigma", "l");
  leg->AddEntry(&gr_klo, "K +/- 3 sigma", "l");
  leg->AddEntry(&gr_plo, "p band (-4..+6 sigma)", "l");
  if (pi_align.IsAligned())
    leg->AddEntry(static_cast<TObject*>(nullptr), pi_align.Summary(), "");
  if (p_align.IsAligned())
    leg->AddEntry(static_cast<TObject*>(nullptr), p_align.Summary(), "");
  return leg;
}

/** One pad: COLZ then curves + legend with "same" (export this canvas to PNG/PDF). */
void DrawPidPageColzSamePad(TCanvas& cv, TH2* hist, Bool_t logx, Bool_t logy, Bool_t logz, TGraph& gr_pi,
                            TGraph& gr_k, TGraph& gr_p, TGraph& gr_cut, TGraph& gr_pilo, TGraph& gr_pihi,
                            TGraph& gr_klo, TGraph& gr_khi, TGraph& gr_plo, TGraph& gr_phi, const PiAlign& pi_align,
                            const ProtonAlign& p_align)
{
  cv.Clear();
  StylePidCanvas(cv);
  cv.cd();
  if (logz)
    ConfigureLogZForColz(hist);
  hist->Draw("COLZ");
  cv.SetLogx(logx);
  cv.SetLogy(logy);
  cv.SetLogz(logz);
  DrawBoundaryCurves(gr_pi, gr_k, gr_p, gr_cut, gr_pilo, gr_pihi, gr_klo, gr_khi, gr_plo, gr_phi);
  TLegend* leg = BuildPidLegendUpperRight(gr_pi, gr_k, gr_p, gr_cut, gr_pilo, gr_klo, gr_plo, pi_align, p_align);
  leg->Draw();
  cv.Update();
}

/** Lowest Y bin edge that has content>0; fallback for empty hist. */
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

/** Minimum positive bin content for COLZ log palette (Z = counts / weight). */
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

} // namespace

static bool ShellCommandSucceeded(const char* cmd)
{
  const int st = std::system(cmd);
  if (st == -1)
    return false;
  if (!WIFEXITED(st))
    return false;
  return WEXITSTATUS(st) == 0;
}

/** Raster multipage PDF from two PNGs (matches on-screen); false => caller may use ROOT vector Print. */
static bool MergeTwoPngsToPdf(const TString& pdfPath, const TString& p1, const TString& p2)
{
  TString cmd;
  cmd.Form("img2pdf -o '%s' '%s' '%s'", pdfPath.Data(), p1.Data(), p2.Data());
  if (ShellCommandSucceeded(cmd.Data()))
    return true;
  cmd.Form("magick '%s' '%s' '%s'", p1.Data(), p2.Data(), pdfPath.Data());
  if (ShellCommandSucceeded(cmd.Data()))
    return true;
  cmd.Form("convert '%s' '%s' '%s'", p1.Data(), p2.Data(), pdfPath.Data());
  return ShellCommandSucceeded(cmd.Data());
}

int main(int argc, char** argv)
{
  Int_t run = -1; // -1 = not set by -r; use parse or 0
  TString inPath;
  TString outPdf;
  PiAlign pi_align;
  ProtonAlign proton_align;

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-r" && i + 1 < argc) {
      run = TString(argv[++i]).Atoi();
    } else if (a == "-o" && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (a == "--pi-baseline") {
      pi_align.center_at_legacy_minus3sigma = false;
    } else if (a == "--p-baseline") {
      proton_align.enabled = false;
    } else if (a == "--p-align-sigma" && i + 1 < argc) {
      proton_align.sigma_shift = TString(argv[++i]).Atof();
      proton_align.enabled = true;
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

  PrintPiAlignReference(pi_align);
  PrintProtonAlignReference(proton_align);

  TH2* h = LoadOrBuildHistogram(inPath);
  if (!h) {
    std::cerr << "Could not load histogram or fill from tree in " << inPath << std::endl;
    return 1;
  }

  const Int_t ngrid = 601;
  std::vector<Double_t> gx(ngrid), gpi(ngrid), gk(ngrid), gp(ngrid), gcut(ngrid);
  std::vector<Double_t> gpi_lo(ngrid), gpi_hi(ngrid), gk_lo(ngrid), gk_hi(ngrid), gp_lo(ngrid), gp_hi(ngrid);
  const Double_t pmin = kBinsSignedP[1];
  const Double_t pmax = kBinsSignedP[2];
  for (int i = 0; i < ngrid; ++i) {
    const Double_t t = static_cast<Double_t>(i) / static_cast<Double_t>(ngrid - 1);
    const Double_t poq = pmin + t * (pmax - pmin);
    gx[i] = poq;
    Double_t dedx_bethe = 0.;
    Double_t spi = 0.;
    pi_align.LegacyPiExpectation(poq, dedx_bethe, spi);
    gpi[i] = pi_align.PiCenter(poq);
    gk[i] = HypTPCdEdxKaonCurve(poq);
    gp[i] = proton_align.ProtonCenter(poq);
    gcut[i] = PpiSeparationCut(poq);
    const Double_t spk = CalcTPCdEdxSigma(sigma_dedx_k, poq);
    const Double_t spp = CalcTPCdEdxSigma(sigma_dedx_p, poq);
    gpi_lo[i] = gpi[i] - 3. * spi;
    gpi_hi[i] = gpi[i] + 3. * spi;
    gk_lo[i] = gk[i] - 3. * spk;
    gk_hi[i] = gk[i] + 3. * spk;
    gp_lo[i] = gp[i] - 4. * spp;
    gp_hi[i] = gp[i] + 6. * spp;
  }

  TGraph gr_pi(ngrid, gx.data(), gpi.data());
  gr_pi.SetLineColor(kRed);
  gr_pi.SetLineWidth(2);
  TGraph gr_k(ngrid, gx.data(), gk.data());
  gr_k.SetLineColor(kGreen + 2);
  gr_k.SetLineWidth(2);
  TGraph gr_p(ngrid, gx.data(), gp.data());
  gr_p.SetLineColor(kBlue);
  gr_p.SetLineWidth(2);
  TGraph gr_cut(ngrid, gx.data(), gcut.data());
  gr_cut.SetLineColor(kMagenta);
  gr_cut.SetLineStyle(2);
  gr_cut.SetLineWidth(2);
  TGraph gr_pilo(ngrid, gx.data(), gpi_lo.data());
  gr_pilo.SetLineColor(kRed);
  gr_pilo.SetLineStyle(3);
  TGraph gr_pihi(ngrid, gx.data(), gpi_hi.data());
  gr_pihi.SetLineColor(kRed);
  gr_pihi.SetLineStyle(3);
  TGraph gr_klo(ngrid, gx.data(), gk_lo.data());
  gr_klo.SetLineColor(kGreen + 2);
  gr_klo.SetLineStyle(3);
  TGraph gr_khi(ngrid, gx.data(), gk_hi.data());
  gr_khi.SetLineColor(kGreen + 2);
  gr_khi.SetLineStyle(3);
  TGraph gr_plo(ngrid, gx.data(), gp_lo.data());
  gr_plo.SetLineColor(kBlue);
  gr_plo.SetLineStyle(3);
  TGraph gr_phi(ngrid, gx.data(), gp_hi.data());
  gr_phi.SetLineColor(kBlue);
  gr_phi.SetLineStyle(3);

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
  const Bool_t any_align = pi_align.IsAligned() || proton_align.IsAligned();
  if (outPdf.IsNull()) {
    if (any_align)
      outPdf = Form("%s/PID_dEdx_vs_SignedMom_boundaries_aligned_run%05d.pdf", imgDir.Data(), runOut);
    else
      outPdf = Form("%s/PID_dEdx_vs_SignedMom_boundaries_run%05d.pdf", imgDir.Data(), runOut);
  }

  TH2* h_p1 = static_cast<TH2*>(h->Clone("h_pid_p1"));
  h_p1->SetDirectory(nullptr);
  SanitizePdfTitlesSignedQ(h_p1);
  AppendAlignSubtitle(h_p1, pi_align, proton_align);

  TH2* h_p2 = static_cast<TH2*>(h->Clone("h_pid_p2"));
  h_p2->SetDirectory(nullptr);
  SanitizePdfTitlesSignedQ(h_p2);
  h_p2->SetTitle("PID dE/dx (signed qp log Y);qp [GeV/c];dE/dx (a.u.)");
  AppendAlignSubtitle(h_p2, pi_align, proton_align);
  const Double_t yLogMax = kBinsDe[2];
  const Double_t yLogMin = TMath::Max(1e-6, LogYAxisMin(h_p2));
  h_p2->GetYaxis()->SetRangeUser(yLogMin, yLogMax);

  constexpr Int_t kCanvasW = 1400;
  constexpr Int_t kCanvasH = 1000;

  if (!outPdf.IsNull())
    std::remove(outPdf.Data());

  TCanvas c1("c_pid_p1", "", kCanvasW, kCanvasH);
  DrawPidPageColzSamePad(c1, h_p1, kFALSE, kFALSE, kTRUE, gr_pi, gr_k, gr_p, gr_cut, gr_pilo, gr_pihi, gr_klo, gr_khi,
                        gr_plo, gr_phi, pi_align, proton_align);

  TCanvas c2("c_pid_p2", "", kCanvasW, kCanvasH);
  DrawPidPageColzSamePad(c2, h_p2, kFALSE, kTRUE, kTRUE, gr_pi, gr_k, gr_p, gr_cut, gr_pilo, gr_pihi, gr_klo, gr_khi,
                        gr_plo, gr_phi, pi_align, proton_align);

  const TString tmpP1 = Form("/tmp/tpc_pid_dedx_boundaries_%d_p1.png", static_cast<int>(getpid()));
  const TString tmpP2 = Form("/tmp/tpc_pid_dedx_boundaries_%d_p2.png", static_cast<int>(getpid()));
  c1.SaveAs(tmpP1.Data());
  c2.SaveAs(tmpP2.Data());

  bool rasterPdf = MergeTwoPngsToPdf(outPdf, tmpP1, tmpP2);
  std::remove(tmpP1.Data());
  std::remove(tmpP2.Data());
  if (!rasterPdf) {
    std::cerr << "img2pdf / ImageMagick not available or merge failed; writing ROOT vector PDF instead.\n";
    c1.Print(outPdf + "(");
    c2.Print(outPdf + ")");
  }

  delete h_p2;
  delete h_p1;
  delete h;
  if (rasterPdf)
    std::cout << "Wrote 2-page raster PDF (PNG embedded): " << outPdf << "\n";
  else
    std::cout << "Wrote 2-page vector PDF (fallback): " << outPdf << "\n";
  return 0;
}
