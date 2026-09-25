// -*- C++ -*-
// Kp systematics isolation + tip offset QA (merged).
//
// Usage:
//   kpsc_systematics_isolation <baseline.root> [-o out.pdf] [-t out.txt] [-r run]
//     [--p0-tag label --p0-file file.root]...
//
// Run number: default from tree branch run_number (minimum = youngest).
//   Fallback: parse "runNNNNN" from path. Override with -r.
// Output default:
//   OUTPUT_DIR/img/runXXXXX/kpsc_systematics_isolation_runXXXXX.{pdf,txt}
//
// Tip cut:
//   mmass>0 && effective_ntTpc==2 && diff_angle<0.1 && close_dist<5

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
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace
{

constexpr Double_t kMK = 0.494;
constexpr Double_t kMP = 0.938272;
constexpr const char* kCutTip =
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.1 && close_dist<5";
constexpr const char* kCutTipPxP =
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.1 && close_dist<5 && px_p>0";
constexpr const char* kCutTipPxM =
  "mmass>0 && effective_ntTpc==2 && diff_angle<0.1 && close_dist<5 && px_p<0";

struct P0Sample
{
  TString tag;
  TString path;
};

void
usage(const char* a0)
{
  std::cerr
    << "Usage: " << a0
    << " <baseline.root> [-o out.pdf] [-t out.txt] [-r run]\n"
    << "  [--p0-tag <label> --p0-file <kpsc.root>]...\n"
    << "  Run: tree run_number minimum (youngest); else path; else -r.\n"
    << "  PDF/TXT default: OUTPUT_DIR/img/runXXXXX/kpsc_systematics_isolation_runXXXXX.*\n";
}

Int_t
ParseRunFromTree(TTree* t)
{
  if (!t || t->GetEntries() <= 0)
    return -1;
  if (!t->GetBranch("run_number"))
    return -1;
  const Double_t mn = t->GetMinimum("run_number");
  if (!std::isfinite(mn) || mn < 0.)
    return -1;
  return static_cast<Int_t>(mn);
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
    if (s[j] < '0' || s[j] > '9')
      continue;
    Ssiz_t k = j;
    while (k < n && s[k] >= '0' && s[k] <= '9')
      ++k;
    if (k - j > 6)
      continue;
    return TString(s(j, k - j)).Atoi();
  }
  return -1;
}

Double_t
MedianOfHist(TH1* h)
{
  if (!h || h->GetEntries() <= 0)
    return std::numeric_limits<Double_t>::quiet_NaN();
  Double_t q = 0.5;
  Double_t x = 0.;
  h->GetQuantiles(1, &x, &q);
  return x;
}

Double_t
PeakX(TH1* h)
{
  if (!h || h->GetEntries() <= 0)
    return std::numeric_limits<Double_t>::quiet_NaN();
  return h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
}

void
StyleHist(TH1* h, Color_t c)
{
  h->SetLineColor(c);
  h->SetLineWidth(2);
  h->SetStats(0);
}

void
PrintPage(TCanvas& c, const TString& pdf)
{
  c.Print(pdf);
  c.Clear();
  c.cd();
}

struct TipStats
{
  Long64_t n = 0;
  Double_t med_lam_p = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_lam_pk = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_cbm = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_d5 = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t peak_mmass = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_abs_dm = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_diff_angle = std::numeric_limits<Double_t>::quiet_NaN();
  Double_t med_theta_p = std::numeric_limits<Double_t>::quiet_NaN();
  Bool_t has_new = false;
};

TipStats
FillTipStats(TTree* t, const char* cut, const char* tag, Bool_t need_new)
{
  TipStats s;
  s.n = t->GetEntries(cut);
  s.has_new = need_new && t->GetBranch("p_p_calc") && t->GetBranch("calc_beam_mom")
              && t->GetBranch("px_p");

  gROOT->cd();
  auto* hmm = new TH1D(Form("hmm_%s", tag), ";M_{miss} [GeV];Counts", 120, 0., 1.2);
  auto* hdm = new TH1D(Form("hdm_%s", tag), ";|M_{miss}-m_K|;Counts", 100, 0., 0.5);
  auto* hda = new TH1D(Form("hda_%s", tag), ";diff_angle;Counts", 100, 0., 0.5);
  auto* hd5 = new TH1D(Form("hd5_%s", tag), ";d5_momentum;Counts", 100, 0.4, 1.2);
  t->Project(hmm->GetName(), "mmass", cut);
  t->Project(hdm->GetName(), "abs(mmass-0.494)", cut);
  t->Project(hda->GetName(), "diff_angle", cut);
  t->Project(hd5->GetName(), "d5_momentum", cut);
  s.peak_mmass = PeakX(hmm);
  s.med_abs_dm = MedianOfHist(hdm);
  s.med_diff_angle = MedianOfHist(hda);
  s.med_d5 = MedianOfHist(hd5);

  if (s.has_new) {
    auto* hlp = new TH1D(Form("hlp_%s", tag), ";#lambda_p;Counts", 120, 0.5, 1.5);
    auto* hlk = new TH1D(Form("hlk_%s", tag), ";#lambda_{pK};Counts", 120, 0.5, 1.5);
    auto* hcb = new TH1D(Form("hcb_%s", tag), ";calc_beam_mom;Counts", 100, 0.4, 1.2);
    auto* htp = new TH1D(Form("htp_%s", tag), ";theta_p_lab;Counts", 90, 0., TMath::Pi());
    const TString cLamP = Form("%s && p_p>0 && p_p_calc>0 && TMath::Finite(p_p_calc)", cut);
    const TString cLamK =
      Form("%s && calc_beam_mom>0 && d5_momentum>0 && TMath::Finite(calc_beam_mom)", cut);
    t->Project(hlp->GetName(), "p_p_calc/p_p", cLamP);
    t->Project(hlk->GetName(), "d5_momentum/calc_beam_mom", cLamK);
    t->Project(hcb->GetName(), "calc_beam_mom", cut);
    if (t->GetBranch("theta_p_lab"))
      t->Project(htp->GetName(), "theta_p_lab", cut);
    s.med_lam_p = MedianOfHist(hlp);
    s.med_lam_pk = MedianOfHist(hlk);
    s.med_cbm = MedianOfHist(hcb);
    s.med_theta_p = MedianOfHist(htp);
  }
  return s;
}

void
Draw2D(TTree* t, TCanvas& c, const TString& pdf, const char* varexp, const char* cut,
       const char* hname, const char* title, Int_t nbx, Double_t xlo, Double_t xhi,
       Int_t nby, Double_t ylo, Double_t yhi)
{
  gROOT->cd();
  if (auto* old = gROOT->FindObject(hname))
    delete old;
  auto* h = new TH2D(hname, title, nbx, xlo, xhi, nby, ylo, yhi);
  t->Project(hname, varexp, cut);
  h->SetStats(0);
  h->Draw("COLZ");
  c.SetLogz(h->GetEntries() > 0);
  PrintPage(c, pdf);
  c.SetLogz(0);
}

void
DrawOffsetQa(TTree* t, TCanvas& c, const TString& pdf, const TString& inPath,
             const TipStats& sall, const TipStats& spxp, const TipStats& spxm)
{
  if (!sall.has_new)
    return;

  gROOT->cd();
  auto* hLam = new TH1D("qa_hLam", ";#lambda_p = p_p_calc/p_p;Counts", 100, 0.5, 1.5);
  auto* hLamK = new TH1D("qa_hLamK", ";#lambda_{pK} = d5/calc_beam_mom;Counts", 100, 0.5, 1.5);
  auto* hLamP = new TH1D("qa_hLamP", "px>0;#lambda_p;Counts", 100, 0.5, 1.5);
  auto* hLamM = new TH1D("qa_hLamM", "px<0;#lambda_p;Counts", 100, 0.5, 1.5);
  auto* hPP = new TH2D("qa_hPP", "tip;p_p_calc [GeV/c];p_p [GeV/c]", 80, 0., 1.2, 80, 0., 1.2);
  auto* hCB = new TH2D("qa_hCB", "tip;d5_momentum;calc_beam_mom", 80, 0.4, 1.2, 80, 0.4, 1.2);

  const TString cLam = Form("%s && p_p>0 && p_p_calc>0", kCutTip);
  const TString cLamK = Form("%s && calc_beam_mom>0 && d5_momentum>0", kCutTip);
  t->Project("qa_hLam", "p_p_calc/p_p", cLam);
  t->Project("qa_hLamK", "d5_momentum/calc_beam_mom", cLamK);
  t->Project("qa_hLamP", "p_p_calc/p_p", Form("%s && px_p>0", cLam.Data()));
  t->Project("qa_hLamM", "p_p_calc/p_p", Form("%s && px_p<0", cLam.Data()));
  t->Project("qa_hPP", "p_p:p_p_calc", cLam);
  t->Project("qa_hCB", "calc_beam_mom:d5_momentum", cLamK);

  const Double_t medLp = sall.med_lam_p;
  const Double_t medLk = sall.med_lam_pk;

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.055);
  tx.DrawLatex(0.10, 0.88, "QA: tip offsets (#lambda)");
  tx.SetTextSize(0.028);
  tx.DrawLatex(0.10, 0.80, Form("input: %s", inPath.Data()));
  tx.DrawLatex(0.10, 0.74,
               Form("tip N = %lld   (px>0: %lld,  px<0: %lld)", sall.n, spxp.n, spxm.n));
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.10, 0.62,
               Form("#lambda_p median = %.3f   (%.1f %%)", medLp, 100. * (medLp - 1.)));
  tx.DrawLatex(0.10, 0.54,
               Form("   px>0: %.3f (%.1f %%)    px<0: %.3f (%.1f %%)", spxp.med_lam_p,
                    100. * (spxp.med_lam_p - 1.), spxm.med_lam_p,
                    100. * (spxm.med_lam_p - 1.)));
  tx.DrawLatex(0.10, 0.42,
               Form("#lambda_pK median = %.3f   (%.1f %%)", medLk, 100. * (medLk - 1.)));
  tx.SetTextSize(0.028);
  tx.DrawLatex(0.10, 0.28, "Red dashed on 1D = median; black dashed = 1 (no offset)");
  tx.DrawLatex(0.10, 0.22, "#lambda_pK >> #lambda_p => not a single common B-scale");
  PrintPage(c, pdf);

  c.Divide(2, 1);
  c.cd(1);
  StyleHist(hLam, kViolet + 1);
  hLam->Draw("hist");
  auto* l1 = new TLine(1., 0., 1., hLam->GetMaximum() * 1.05);
  l1->SetLineStyle(2);
  l1->Draw();
  if (std::isfinite(medLp)) {
    auto* l1m = new TLine(medLp, 0., medLp, hLam->GetMaximum() * 1.05);
    l1m->SetLineColor(kRed + 1);
    l1m->Draw();
  }
  c.cd(2);
  StyleHist(hLamK, kTeal + 2);
  hLamK->Draw("hist");
  auto* l2 = new TLine(1., 0., 1., hLamK->GetMaximum() * 1.05);
  l2->SetLineStyle(2);
  l2->Draw();
  if (std::isfinite(medLk)) {
    auto* l2m = new TLine(medLk, 0., medLk, hLamK->GetMaximum() * 1.05);
    l2m->SetLineColor(kRed + 1);
    l2m->Draw();
  }
  PrintPage(c, pdf);

  c.Divide(2, 1);
  c.cd(1);
  StyleHist(hLamP, kBlue + 1);
  hLamP->Draw("hist");
  auto* lp = new TLine(1., 0., 1., std::max(hLamP->GetMaximum(), 1.) * 1.05);
  lp->SetLineStyle(2);
  lp->Draw();
  c.cd(2);
  StyleHist(hLamM, kOrange + 7);
  hLamM->Draw("hist");
  auto* lm = new TLine(1., 0., 1., std::max(hLamM->GetMaximum(), 1.) * 1.05);
  lm->SetLineStyle(2);
  lm->Draw();
  PrintPage(c, pdf);

  c.Divide(2, 1);
  c.cd(1);
  hPP->SetStats(0);
  hPP->Draw("COLZ");
  gPad->SetLogz(hPP->GetEntries() > 0);
  auto* d1 = new TLine(0., 0., 1.2, 1.2);
  d1->SetLineColor(kRed);
  d1->SetLineStyle(2);
  d1->Draw();
  c.cd(2);
  hCB->SetStats(0);
  hCB->Draw("COLZ");
  gPad->SetLogz(hCB->GetEntries() > 0);
  auto* d2 = new TLine(0.4, 0.4, 1.2, 1.2);
  d2->SetLineColor(kRed);
  d2->SetLineStyle(2);
  d2->Draw();
  PrintPage(c, pdf);
  c.cd();
  gPad->SetLogz(0);
}

Double_t
OpeningAngleFromMags(Double_t cbm, Double_t pp, Double_t pk)
{
  if (!(cbm > 0.) || !(pp > 0.) || !(pk > 0.))
    return std::numeric_limits<Double_t>::quiet_NaN();
  Double_t cosg = (cbm * cbm - pp * pp - pk * pk) / (2. * pp * pk);
  if (cosg > 1.)
    cosg = 1.;
  if (cosg < -1.)
    cosg = -1.;
  return TMath::ACos(cosg);
}

Double_t
Pearson(const std::vector<Double_t>& x, const std::vector<Double_t>& y);

void
DrawLamPkDiagnosis(TTree* t, TCanvas& c, const TString& pdf, std::ostream& log)
{
  if (!t->GetBranch("calc_beam_mom") || !t->GetBranch("d5_momentum")
      || !t->GetBranch("p_p") || !t->GetBranch("p_k") || !t->GetBranch("px_p"))
    return;

  const Bool_t has_thp = t->GetBranch("theta_p_lab");
  const Bool_t has_thk = t->GetBranch("theta_k_lab");
  const Bool_t has_pm = t->GetBranch("p_miss");

  gROOT->cd();
  auto* hdAll = new TH1D("diag_dcb_all", "tip;calc_beam-d5 [GeV/c];Counts", 80, -0.5, 0.2);
  auto* hdP = new TH1D("diag_dcb_pxp", "tip px>0;calc_beam-d5;Counts", 80, -0.5, 0.2);
  auto* hdM = new TH1D("diag_dcb_pxm", "tip px<0;calc_beam-d5;Counts", 80, -0.5, 0.2);
  auto* hGam = new TH1D("diag_gamma", "tip;#gamma(p,K) from |p|+|k| [rad];Counts", 80, 0., TMath::Pi());
  auto* h2Thp = new TH2D("diag_dcb_thp", "tip;#theta_{p lab};calc_beam-d5", 60, 0., TMath::Pi(),
                         60, -0.5, 0.2);
  auto* h2Thk = new TH2D("diag_dcb_thk", "tip;#theta_{K lab};calc_beam-d5", 60, 0., TMath::Pi(),
                         60, -0.5, 0.2);
  auto* h2Pp = new TH2D("diag_dcb_pp", "tip;p_p [GeV/c];calc_beam-d5", 60, 0., 1.2, 60, -0.5, 0.2);
  auto* h2Pk = new TH2D("diag_dcb_pk", "tip;p_k [GeV/c];calc_beam-d5", 60, 0., 1.2, 60, -0.5, 0.2);
  auto* h2Da = new TH2D("diag_dcb_da", "tip;diff_angle;calc_beam-d5", 50, 0., 0.5, 60, -0.5, 0.2);
  auto* h2Gam = new TH2D("diag_dcb_gam", "tip;#gamma(p,K);calc_beam-d5", 60, 0., TMath::Pi(),
                         60, -0.5, 0.2);
  auto* h2PmPk = new TH2D("diag_pm_pk", "tip;p_k [GeV/c];p_miss [GeV/c]", 60, 0., 1.2, 60, 0., 1.2);
  auto* h2Cbd5 = new TH2D("diag_cb_d5", "tip;d5_momentum;calc_beam_mom", 60, 0.4, 1.2, 60, 0.4, 1.2);

  Double_t mmass = 0., da = 0., cd = 0., d5 = 0., cbm = 0., pp = 0., pk = 0., px = 0.;
  Double_t thp = 0., thk = 0., pm = 0.;
  Int_t nt = 0;
  t->SetBranchAddress("mmass", &mmass);
  t->SetBranchAddress("diff_angle", &da);
  t->SetBranchAddress("close_dist", &cd);
  t->SetBranchAddress("effective_ntTpc", &nt);
  t->SetBranchAddress("d5_momentum", &d5);
  t->SetBranchAddress("calc_beam_mom", &cbm);
  t->SetBranchAddress("p_p", &pp);
  t->SetBranchAddress("p_k", &pk);
  t->SetBranchAddress("px_p", &px);
  if (has_thp)
    t->SetBranchAddress("theta_p_lab", &thp);
  if (has_thk)
    t->SetBranchAddress("theta_k_lab", &thk);
  if (has_pm)
    t->SetBranchAddress("p_miss", &pm);

  std::vector<Double_t> v_dcb, v_thp, v_thk, v_pp, v_pk, v_da, v_gam;
  v_dcb.reserve(5000);

  const Long64_t nent = t->GetEntries();
  for (Long64_t i = 0; i < nent; ++i) {
    t->GetEntry(i);
    if (!(mmass > 0.) || nt != 2 || !(da < 0.1) || !(cd < 5.))
      continue;
    if (!(cbm > 0.) || !(d5 > 0.) || !std::isfinite(cbm) || !std::isfinite(d5))
      continue;
    const Double_t dcb = cbm - d5;
    hdAll->Fill(dcb);
    if (px > 0.)
      hdP->Fill(dcb);
    else if (px < 0.)
      hdM->Fill(dcb);
    h2Da->Fill(da, dcb);
    h2Cbd5->Fill(d5, cbm);
    if (pp > 0.)
      h2Pp->Fill(pp, dcb);
    if (pk > 0.)
      h2Pk->Fill(pk, dcb);
    if (has_thp && std::isfinite(thp))
      h2Thp->Fill(thp, dcb);
    if (has_thk && std::isfinite(thk))
      h2Thk->Fill(thk, dcb);
    if (has_pm && std::isfinite(pm) && pk > 0.)
      h2PmPk->Fill(pk, pm);
    const Double_t gam = OpeningAngleFromMags(cbm, pp, pk);
    if (std::isfinite(gam)) {
      hGam->Fill(gam);
      h2Gam->Fill(gam, dcb);
    }
    v_dcb.push_back(dcb);
    v_da.push_back(da);
    v_pp.push_back(pp);
    v_pk.push_back(pk);
    v_thp.push_back(has_thp && std::isfinite(thp) ? thp : std::numeric_limits<Double_t>::quiet_NaN());
    v_thk.push_back(has_thk && std::isfinite(thk) ? thk : std::numeric_limits<Double_t>::quiet_NaN());
    v_gam.push_back(gam);
  }
  t->ResetBranchAddresses();

  auto pairCorr = [](const std::vector<Double_t>& a, const std::vector<Double_t>& b) {
    std::vector<Double_t> xa, ya;
    xa.reserve(a.size());
    ya.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (!std::isfinite(a[i]) || !std::isfinite(b[i]))
        continue;
      xa.push_back(a[i]);
      ya.push_back(b[i]);
    }
    return Pearson(xa, ya);
  };

  const Double_t medAll = MedianOfHist(hdAll);
  const Double_t medP = MedianOfHist(hdP);
  const Double_t medM = MedianOfHist(hdM);
  const Double_t medGam = MedianOfHist(hGam);
  log << "## StepA lam_pK / calc_beam-d5\n";
  log << "  N_tip_finite=" << v_dcb.size() << "\n";
  log << "  med(cbm-d5)=" << medAll << "  px>0=" << medP << "  px<0=" << medM << "\n";
  log << "  med(gamma)=" << medGam << "\n";
  log << "  corr(cbm-d5, theta_p)=" << pairCorr(v_dcb, v_thp) << "\n";
  log << "  corr(cbm-d5, theta_k)=" << pairCorr(v_dcb, v_thk) << "\n";
  log << "  corr(cbm-d5, p_p)=" << pairCorr(v_dcb, v_pp) << "\n";
  log << "  corr(cbm-d5, p_k)=" << pairCorr(v_dcb, v_pk) << "\n";
  log << "  corr(cbm-d5, diff_angle)=" << pairCorr(v_dcb, v_da) << "\n";
  log << "  corr(cbm-d5, gamma)=" << pairCorr(v_dcb, v_gam) << "\n\n";

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.12, 0.80, "StepA: calc_beam - d5 diagnosis");
  tx.SetTextSize(0.03);
  tx.DrawLatex(0.12, 0.70, Form("med(cbm-d5)=%.4f  (px>0 %.4f, px<0 %.4f)", medAll, medP, medM));
  tx.DrawLatex(0.12, 0.62, Form("med(#gamma)=%.3f rad", medGam));
  tx.DrawLatex(0.12, 0.52, Form("corr vs #theta_p / #theta_K / p_p / p_k / #gamma:"));
  tx.DrawLatex(0.12, 0.44,
               Form("  %.3f / %.3f / %.3f / %.3f / %.3f", pairCorr(v_dcb, v_thp),
                    pairCorr(v_dcb, v_thk), pairCorr(v_dcb, v_pp), pairCorr(v_dcb, v_pk),
                    pairCorr(v_dcb, v_gam)));
  PrintPage(c, pdf);

  c.Divide(3, 1);
  c.cd(1);
  StyleHist(hdAll, kBlack);
  hdAll->Draw("hist");
  c.cd(2);
  StyleHist(hdP, kBlue + 1);
  hdP->Draw("hist");
  c.cd(3);
  StyleHist(hdM, kOrange + 7);
  hdM->Draw("hist");
  PrintPage(c, pdf);

  auto draw2 = [&](TH2* h) {
    c.cd();
    h->SetStats(0);
    h->Draw("COLZ");
    c.SetLogz(h->GetEntries() > 0);
    PrintPage(c, pdf);
    c.SetLogz(0);
  };
  draw2(h2Thp);
  draw2(h2Thk);
  draw2(h2Pp);
  draw2(h2Pk);
  draw2(h2Da);
  draw2(h2Gam);

  c.Divide(2, 1);
  c.cd(1);
  h2PmPk->SetStats(0);
  h2PmPk->Draw("COLZ");
  gPad->SetLogz(h2PmPk->GetEntries() > 0);
  auto* diag = new TLine(0., 0., 1.2, 1.2);
  diag->SetLineColor(kRed);
  diag->SetLineStyle(2);
  diag->Draw();
  c.cd(2);
  h2Cbd5->SetStats(0);
  h2Cbd5->Draw("COLZ");
  gPad->SetLogz(h2Cbd5->GetEntries() > 0);
  auto* d2 = new TLine(0.4, 0.4, 1.2, 1.2);
  d2->SetLineColor(kRed);
  d2->SetLineStyle(2);
  d2->Draw();
  PrintPage(c, pdf);
  c.cd();
  gPad->SetLogz(0);

  StyleHist(hGam, kTeal + 2);
  hGam->Draw("hist");
  PrintPage(c, pdf);
}

void
DrawLamPHemisphere(TTree* t, TCanvas& c, const TString& pdf, std::ostream& log)
{
  if (!t->GetBranch("p_p_calc") || !t->GetBranch("p_p") || !t->GetBranch("px_p"))
    return;

  const Bool_t has_thp = t->GetBranch("theta_p_lab");

  gROOT->cd();
  auto* h2Th = new TH2D("diag_lam_th", "tip;#theta_{p lab};#lambda_p", 60, 0., TMath::Pi(), 60,
                        0.5, 1.5);
  auto* h2ThP = new TH2D("diag_lam_th_p", "tip px>0;#theta_{p lab};#lambda_p", 60, 0., TMath::Pi(),
                         60, 0.5, 1.5);
  auto* h2ThM = new TH2D("diag_lam_th_m", "tip px<0;#theta_{p lab};#lambda_p", 60, 0., TMath::Pi(),
                         60, 0.5, 1.5);
  auto* h2Px = new TH2D("diag_lam_px", "tip;px_p [GeV/c];#lambda_p", 60, -1., 1., 60, 0.5, 1.5);
  auto* h2Pp = new TH2D("diag_lam_pp", "tip;p_p [GeV/c];#lambda_p", 60, 0., 1.2, 60, 0.5, 1.5);

  constexpr Int_t kNThetaBin = 6;
  const Double_t thLo = 0.;
  const Double_t thHi = TMath::Pi();
  const Double_t dth = (thHi - thLo) / kNThetaBin;
  std::vector<std::vector<Double_t>> lamBins(kNThetaBin), lamBinsP(kNThetaBin), lamBinsM(kNThetaBin);

  Double_t mmass = 0., da = 0., cd = 0., pp = 0., ppc = 0., px = 0., thp = 0.;
  Int_t nt = 0;
  t->SetBranchAddress("mmass", &mmass);
  t->SetBranchAddress("diff_angle", &da);
  t->SetBranchAddress("close_dist", &cd);
  t->SetBranchAddress("effective_ntTpc", &nt);
  t->SetBranchAddress("p_p", &pp);
  t->SetBranchAddress("p_p_calc", &ppc);
  t->SetBranchAddress("px_p", &px);
  if (has_thp)
    t->SetBranchAddress("theta_p_lab", &thp);

  const Long64_t nent = t->GetEntries();
  for (Long64_t i = 0; i < nent; ++i) {
    t->GetEntry(i);
    if (!(mmass > 0.) || nt != 2 || !(da < 0.1) || !(cd < 5.))
      continue;
    if (!(pp > 0.) || !(ppc > 0.) || !std::isfinite(ppc))
      continue;
    const Double_t lam = ppc / pp;
    h2Px->Fill(px, lam);
    h2Pp->Fill(pp, lam);
    if (has_thp && std::isfinite(thp)) {
      h2Th->Fill(thp, lam);
      if (px > 0.)
        h2ThP->Fill(thp, lam);
      else if (px < 0.)
        h2ThM->Fill(thp, lam);
      if (thp >= thLo && thp < thHi) {
        Int_t ib = static_cast<Int_t>((thp - thLo) / dth);
        if (ib < 0)
          ib = 0;
        if (ib >= kNThetaBin)
          ib = kNThetaBin - 1;
        lamBins[ib].push_back(lam);
        if (px > 0.)
          lamBinsP[ib].push_back(lam);
        else if (px < 0.)
          lamBinsM[ib].push_back(lam);
      }
    }
  }
  t->ResetBranchAddresses();

  auto medVec = [](std::vector<Double_t> v) {
    if (v.empty())
      return std::numeric_limits<Double_t>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2)
      return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
  };

  log << "## StepB lam_p hemisphere\n";
  for (Int_t ib = 0; ib < kNThetaBin; ++ib) {
    const Double_t t0 = thLo + ib * dth;
    const Double_t t1 = t0 + dth;
    log << "  theta_bin[" << t0 << "," << t1 << ") N=" << lamBins[ib].size()
        << " med_lam-1=" << (medVec(lamBins[ib]) - 1.)
        << " px>0=" << (medVec(lamBinsP[ib]) - 1.)
        << " px<0=" << (medVec(lamBinsM[ib]) - 1.) << "\n";
  }
  log << "\n";

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.12, 0.78, "StepB: #lambda_p hemisphere vs #theta_p / p");
  tx.SetTextSize(0.028);
  for (Int_t ib = 0; ib < kNThetaBin; ++ib) {
    const Double_t t0 = thLo + ib * dth;
    tx.DrawLatex(0.12, 0.68 - 0.07 * ib,
                 Form("[%4.2f,%4.2f) med(#lambda-1)=%.3f  (+ %.3f / - %.3f)  N=%zu", t0, t0 + dth,
                      medVec(lamBins[ib]) - 1., medVec(lamBinsP[ib]) - 1.,
                      medVec(lamBinsM[ib]) - 1., lamBins[ib].size()));
  }
  PrintPage(c, pdf);

  auto draw2 = [&](TH2* h) {
    c.cd();
    h->SetStats(0);
    h->Draw("COLZ");
    c.SetLogz(h->GetEntries() > 0);
    PrintPage(c, pdf);
    c.SetLogz(0);
  };
  draw2(h2Th);
  draw2(h2ThP);
  draw2(h2ThM);
  draw2(h2Px);
  draw2(h2Pp);
}

void
DrawScaledMmass(TTree* t, TCanvas& c, const TString& pdf, Double_t lam,
                const char* cut, const char* tag, std::ostream& log)
{
  Double_t mmass, d5, px, py, pz, pp, da, cd;
  Int_t nt;
  t->SetBranchAddress("mmass", &mmass);
  t->SetBranchAddress("d5_momentum", &d5);
  t->SetBranchAddress("px_p", &px);
  t->SetBranchAddress("py_p", &py);
  t->SetBranchAddress("pz_p", &pz);
  t->SetBranchAddress("p_p", &pp);
  t->SetBranchAddress("diff_angle", &da);
  t->SetBranchAddress("effective_ntTpc", &nt);
  t->SetBranchAddress("close_dist", &cd);

  gROOT->cd();
  auto* h0 = new TH1D(Form("hmm0_%s", tag), ";M_{miss} [GeV];Counts", 100, 0., 0.8);
  auto* h1 = new TH1D(Form("hmm1_%s", tag), ";M_{miss} [GeV];Counts", 100, 0., 0.8);
  auto* h2 = new TH2D(Form("hv0_%s", tag), "baseline tip;diff_angle;M_{miss}",
                      50, 0., 0.5, 80, 0.2, 0.8);
  auto* h3 = new TH2D(Form("hv1_%s", tag),
                      Form("p_p #times %.4f (beam#parallel z approx);diff_angle;M_{miss}", lam),
                      50, 0., 0.5, 80, 0.2, 0.8);

  const Bool_t req_pxp = (TString(cut).Contains("px_p>0"));
  const Bool_t req_pxm = (TString(cut).Contains("px_p<0"));
  Long64_t nuse = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (!(mmass > 0.) || nt != 2 || !(da < 0.1) || !(cd < 5.))
      continue;
    if (req_pxp && !(px > 0.))
      continue;
    if (req_pxm && !(px < 0.))
      continue;
    if (!(pp > 1.e-6) || !std::isfinite(px) || !std::isfinite(d5) || !(d5 > 0.))
      continue;
    const TVector3 mom0(px, py, pz);
    if (mom0.Mag() < 1.e-6)
      continue;
    const TVector3 mom1 = mom0.Unit() * (lam * pp);
    const Double_t e_b = TMath::Sqrt(d5 * d5 + kMK * kMK);
    const TLorentzVector lv_b(0., 0., d5, e_b);
    const TLorentzVector lv_t(0., 0., 0., kMP);
    const Double_t e0 = TMath::Sqrt(mom0.Mag2() + kMP * kMP);
    const Double_t e1 = TMath::Sqrt(mom1.Mag2() + kMP * kMP);
    const TLorentzVector lv_p0(mom0.X(), mom0.Y(), mom0.Z(), e0);
    const TLorentzVector lv_p1(mom1.X(), mom1.Y(), mom1.Z(), e1);
    h0->Fill((lv_b + lv_t - lv_p0).M());
    h1->Fill((lv_b + lv_t - lv_p1).M());
    h2->Fill(da, (lv_b + lv_t - lv_p0).M());
    h3->Fill(da, (lv_b + lv_t - lv_p1).M());
    ++nuse;
  }

  auto* hd0 = new TH1D(Form("hdm0_%s", tag), "", 100, 0., 0.5);
  auto* hd1 = new TH1D(Form("hdm1_%s", tag), "", 100, 0., 0.5);
  for (Int_t b = 1; b <= h0->GetNbinsX(); ++b) {
    if (h0->GetBinContent(b) > 0.)
      hd0->Fill(TMath::Abs(h0->GetBinCenter(b) - kMK), h0->GetBinContent(b));
    if (h1->GetBinContent(b) > 0.)
      hd1->Fill(TMath::Abs(h1->GetBinCenter(b) - kMK), h1->GetBinContent(b));
  }
  log << "Step1 virtual mmass (" << tag << "): N=" << nuse << " lam=" << lam
      << " peak0=" << PeakX(h0) << " peak1=" << PeakX(h1) << "\n";
  log << "  med|mmass-mK| baseline~=" << MedianOfHist(hd0)
      << " scaled~=" << MedianOfHist(hd1) << "\n";

  c.Divide(2, 2);
  c.cd(1);
  StyleHist(h0, kGray + 2);
  StyleHist(h1, kRed + 1);
  h0->Draw("hist");
  h1->Draw("hist same");
  auto* leg = new TLegend(0.55, 0.7, 0.88, 0.88);
  leg->AddEntry(h0, "baseline (z-beam approx)", "l");
  leg->AddEntry(h1, Form("p_p#times%.3f", lam), "l");
  leg->Draw();
  c.cd(2);
  h2->SetStats(0);
  h2->Draw("COLZ");
  c.cd(3);
  h3->SetStats(0);
  h3->Draw("COLZ");
  c.cd(4);
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.1, 0.7, "Note: beam #parallel z approx;");
  tx.DrawLatex(0.1, 0.6, "diff_angle from tree");
  PrintPage(c, pdf);
  c.Clear();
  c.Divide(1, 1);
  c.cd();
}

Double_t
Pearson(const std::vector<Double_t>& x, const std::vector<Double_t>& y)
{
  if (x.size() < 2 || x.size() != y.size())
    return std::numeric_limits<Double_t>::quiet_NaN();
  const Double_t mx = std::accumulate(x.begin(), x.end(), 0.) / x.size();
  const Double_t my = std::accumulate(y.begin(), y.end(), 0.) / y.size();
  Double_t num = 0., dx = 0., dy = 0.;
  for (std::size_t i = 0; i < x.size(); ++i) {
    num += (x[i] - mx) * (y[i] - my);
    dx += (x[i] - mx) * (x[i] - mx);
    dy += (y[i] - my) * (y[i] - my);
  }
  if (!(dx > 0.) || !(dy > 0.))
    return std::numeric_limits<Double_t>::quiet_NaN();
  return num / TMath::Sqrt(dx * dy);
}

} // namespace

int
main(int argc, char** argv)
{
  TString inPath, outPdf, outTxt;
  Int_t run = -1;
  std::vector<P0Sample> p0s;

  for (int i = 1; i < argc; ++i) {
    TString a(argv[i]);
    if (a == "-o" && i + 1 < argc)
      outPdf = argv[++i];
    else if (a == "-t" && i + 1 < argc)
      outTxt = argv[++i];
    else if (a == "-r" && i + 1 < argc)
      run = TString(argv[++i]).Atoi();
    else if (a == "--p0-tag" && i + 1 < argc) {
      P0Sample s;
      s.tag = argv[++i];
      if (i + 2 < argc && TString(argv[i + 1]) == "--p0-file") {
        i += 2;
        s.path = argv[i];
        p0s.push_back(s);
      } else {
        std::cerr << "--p0-tag needs --p0-file\n";
        return 1;
      }
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.BeginsWith("-")) {
      usage(argv[0]);
      return 1;
    } else if (inPath.IsNull()) {
      inPath = a;
    } else {
      usage(argv[0]);
      return 1;
    }
  }
  if (inPath.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile fin(inPath, "READ");
  if (!fin.IsOpen() || fin.IsZombie()) {
    std::cerr << "cannot open " << inPath << "\n";
    return 1;
  }
  auto* t = dynamic_cast<TTree*>(fin.Get("kpsc"));
  if (!t) {
    std::cerr << "no kpsc tree\n";
    return 1;
  }

  Int_t runOut = run;
  if (runOut < 0) {
    runOut = ParseRunFromTree(t);
    if (runOut < 0)
      runOut = ParseRunFromPath(inPath);
  }
  if (runOut < 0)
    runOut = 0;

  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, runOut);
  if (outPdf.IsNull())
    outPdf = Form("%s/kpsc_systematics_isolation_run%05d.pdf", imgDir.Data(), runOut);
  if (outTxt.IsNull())
    outTxt = Form("%s/kpsc_systematics_isolation_run%05d.txt", imgDir.Data(), runOut);

  std::ofstream log(outTxt.Data());
  if (!log) {
    std::cerr << "cannot write " << outTxt << "\n";
    return 1;
  }

  log << "# Kp systematics isolation (+ tip offset QA)\n";
  log << "baseline: " << inPath << "\n";
  log << "run_out: " << runOut << "\n";
  log << "tip: " << kCutTip << "\n";
  log << "entries_all: " << t->GetEntries() << "\n";

  const Bool_t has_new = t->GetBranch("p_p_calc") && t->GetBranch("px_p")
                         && t->GetBranch("calc_beam_mom");
  log << "has_new_branches: " << (has_new ? 1 : 0) << "\n\n";

  TipStats sall = FillTipStats(t, kCutTip, "all", has_new);
  TipStats spxp, spxm;
  if (has_new) {
    spxp = FillTipStats(t, kCutTipPxP, "pxp", true);
    spxm = FillTipStats(t, kCutTipPxM, "pxm", true);
  }

  auto dump = [&](const char* name, const TipStats& s) {
    log << "[" << name << "]\n  N=" << s.n << "\n";
    log << "  peak_mmass=" << s.peak_mmass << " med|mmass-mK|=" << s.med_abs_dm << "\n";
    log << "  med_diff_angle=" << s.med_diff_angle << " med_d5=" << s.med_d5 << "\n";
    if (s.has_new) {
      log << "  med_lam_p-1=" << (s.med_lam_p - 1.)
          << " med_lam_pK-1=" << (s.med_lam_pk - 1.) << "\n";
      log << "  med_calc_beam_mom=" << s.med_cbm
          << " cbm-d5=" << (s.med_cbm - s.med_d5) << "\n";
    }
    log << "\n";
  };
  log << "## Step0\n";
  dump("tip_all", sall);
  if (has_new) {
    dump("tip_px_p>0", spxp);
    dump("tip_px_p<0", spxm);
  }

  TCanvas c("c", "kpsc systematics", 1000, 700);
  c.Print(outPdf + "[");

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.045);
  tx.DrawLatex(0.12, 0.75, "Kp systematics isolation + QA");
  tx.SetTextSize(0.03);
  tx.DrawLatex(0.12, 0.65, Form("%s", inPath.Data()));
  tx.DrawLatex(0.12, 0.58, Form("run %05d  tip N=%lld", runOut, sall.n));
  PrintPage(c, outPdf);

  DrawOffsetQa(t, c, outPdf, inPath, sall, spxp, spxm);

  if (has_new) {
    DrawLamPkDiagnosis(t, c, outPdf, log);
    DrawLamPHemisphere(t, c, outPdf, log);
  }

  Draw2D(t, c, outPdf, "mmass:diff_angle", kCutTip, "h_mm_da",
         "tip; diff_angle; M_{miss} [GeV]", 50, 0., 0.5, 80, 0., 0.8);
  if (has_new) {
    Draw2D(t, c, outPdf, "mmass:diff_angle", kCutTipPxP, "h_mm_da_p",
           "tip px>0; diff_angle; M_{miss}", 50, 0., 0.5, 80, 0., 0.8);
    Draw2D(t, c, outPdf, "mmass:diff_angle", kCutTipPxM, "h_mm_da_m",
           "tip px<0; diff_angle; M_{miss}", 50, 0., 0.5, 80, 0., 0.8);
  }

  if (auto* h = dynamic_cast<TH1*>(fin.Get("KpScat_Proton_PVtx_Minus_Mom0_Nt2"))) {
    h->SetStats(0);
    h->Draw("hist");
    PrintPage(c, outPdf);
  }
  if (auto* h = dynamic_cast<TH1*>(fin.Get("KpScat_Resi_RkTpc"))) {
    h->SetStats(0);
    h->Draw("hist");
    PrintPage(c, outPdf);
  }
  if (auto* h = dynamic_cast<TH1*>(fin.Get("KpScat_N_Veto"))) {
    h->SetStats(0);
    h->Draw("hist");
    PrintPage(c, outPdf);
  }

  log << "## Step1\n";
  if (has_new && std::isfinite(sall.med_lam_p) && sall.med_lam_p > 0.) {
    DrawScaledMmass(t, c, outPdf, sall.med_lam_p, kCutTip, "all", log);
    if (std::isfinite(spxp.med_lam_p) && spxp.med_lam_p > 0.)
      DrawScaledMmass(t, c, outPdf, spxp.med_lam_p, kCutTipPxP, "pxp", log);
    if (std::isfinite(spxm.med_lam_p) && spxm.med_lam_p > 0.)
      DrawScaledMmass(t, c, outPdf, spxm.med_lam_p, kCutTipPxM, "pxm", log);
  }

  log << "## Step2/3 p0 scan\n";
  struct Row { TString tag; TipStats st; };
  std::vector<Row> rows{{"baseline", sall}};
  for (const auto& s : p0s) {
    TFile f(s.path, "READ");
    if (!f.IsOpen() || f.IsZombie())
      continue;
    auto* tp = dynamic_cast<TTree*>(f.Get("kpsc"));
    if (!tp)
      continue;
    const Bool_t neu = tp->GetBranch("p_p_calc") && tp->GetBranch("calc_beam_mom");
    TipStats st = FillTipStats(tp, kCutTip, Form("p0_%s", s.tag.Data()), neu);
    rows.push_back({s.tag, st});
    dump(Form("p0_%s", s.tag.Data()), st);
  }

  std::vector<Double_t> ys_dm, ys_dp, ys_dth;
  for (const auto& r : rows) {
    ys_dm.push_back(r.st.med_abs_dm);
    ys_dp.push_back(r.st.med_d5);
    ys_dth.push_back(r.st.med_diff_angle);
    log << "  row " << r.tag << " med|mmass-mK|=" << r.st.med_abs_dm
        << " med_d5=" << r.st.med_d5 << " N=" << r.st.n << "\n";
  }
  if (rows.size() >= 2) {
    std::vector<Double_t> ddm, ddp, ddth;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      ddm.push_back(ys_dm[i] - ys_dm[0]);
      ddp.push_back(ys_dp[i] - ys_dp[0]);
      ddth.push_back(ys_dth[i] - ys_dth[0]);
    }
    log << "corr(Delta m, Delta d5)=" << Pearson(ddm, ddp) << "\n";
    log << "corr(Delta m, Delta diff_angle)=" << Pearson(ddm, ddth) << "\n";
  }

  c.cd();
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.1, 0.75, Form("Done. run=%05d", runOut));
  tx.DrawLatex(0.1, 0.65, outPdf.Data());
  PrintPage(c, outPdf);
  c.Print(outPdf + "]");

  log << "\nWrote PDF: " << outPdf << "\n";
  std::cout << "Wrote " << outPdf << "\nWrote " << outTxt << "\n";
  return 0;
}
