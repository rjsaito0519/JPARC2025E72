// Figures for docs/htof/manual/htof_tof_manual.tex. Run it inside docs/htof/manual/fig/:
//   root -l -b -q '<myanalysis>/macro/htof_manual_make_figs.C("<DATA_DIR>", "<calib_root_dir>")'
// DATA_DIR: run0XXXX_Hodo.root / run0XXXX_TPCHelixHTOF.root (analyzer outputs);
// calib_root_dir: run02601_HTOF_{ADC,PHC}_Pi.root (HTOF calibration outputs, for the weak-side and
// PHC figures). Writes tdc_example.pdf, adc_example.pdf, weak_example.pdf, phc_example.pdf,
// qa_dtpi_m2.pdf and phc_fit_steps.pdf in the current directory.
// phc_fit_steps.pdf alone: root -l -b -q -e '.L <myanalysis>/macro/htof_manual_make_figs.C+' -e 'make_phcfit_fig("<calib_root_dir>")'
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TProfile.h>
#include <TBox.h>
#include <TH1F.h>
#include <TPad.h>
#include <TLegend.h>
#include <vector>
#include <algorithm>

namespace {
const char* kRuns[] = {"02601", "02602", "02603", "02604", "02606", "02607"};

void label(Double_t x, Double_t y, const char* s, Double_t size = 0.045)
{
  TLatex t; t.SetNDC(); t.SetTextSize(size); t.DrawLatex(x, y, s);
}

TH1D* sum_hist(const char* data_dir, const char* name)
{
  TH1D* out = nullptr;
  for (const char* r : kRuns) {
    TFile f(Form("%s/run%s_Hodo.root", data_dir, r));
    TH1D* h = (TH1D*)f.Get(name);
    if (!h) continue;
    if (!out) { out = (TH1D*)h->Clone(Form("%s_sum", name)); out->SetDirectory(nullptr); }
    else out->Add(h);
  }
  return out;
}
// Steps of the HTOF PHC fit (same procedure as ana_helper::htof_phc_fit) for one channel, drawn
// in three pads starting at pad0: ridge, fit range, band + profile + fit.
// ref = stage-1 parameters (p1 fixed to it, as in stage 2), dense_ratio as in stage 2.
struct PhcSteps { Double_t fit_lo, fit_hi, w, p0; };
PhcSteps draw_phc_steps(TCanvas& c, Int_t pad0, TH2D* h, const Double_t* ref, Double_t dense_ratio, const char* tag)
{
  const Double_t de_min = 0.2, de_max = 5.0;
  const Int_t nx = h->GetNbinsX(), slice = 2, n_sl = nx / slice;
  const Double_t ybw = h->GetYaxis()->GetBinWidth(1);
  const Double_t ylo = h->GetYaxis()->GetXmin(), yhi = h->GetYaxis()->GetXmax();
  struct Pt { Double_t x = 0, xl = 0, xh = 0, y = 0, ey = 0, s = 0, n = 0; Bool_t ok = false; };
  std::vector<Pt> pts(n_sl);
  auto slice_fit = [&](Int_t isl, Double_t yc, Double_t yh) {
    Pt& p = pts[isl];
    const Int_t b0 = isl * slice + 1, b1 = b0 + slice - 1;
    p.xl = h->GetXaxis()->GetBinLowEdge(b0); p.xh = h->GetXaxis()->GetBinUpEdge(b1); p.x = 0.5 * (p.xl + p.xh);
    TH1D* py = h->ProjectionY(Form("%s_sl%d_fig", h->GetName(), isl), b0, b1);
    const Double_t lo = yh < 0 ? ylo : std::max(ylo, yc - yh), hi = yh < 0 ? yhi : std::min(yhi, yc + yh);
    Int_t blo = py->FindBin(lo), bhi = py->FindBin(hi), pb = blo;
    auto s3 = [&](Int_t k) { return py->GetBinContent(k - 1) + py->GetBinContent(k) + py->GetBinContent(k + 1); };
    for (Int_t b = blo; b <= bhi; ++b) if (s3(b) > s3(pb)) pb = b;
    const Double_t pk = py->GetBinCenter(pb);
    const Double_t ncore = py->Integral(py->FindBin(pk - 1.0), py->FindBin(pk + 1.0));
    p.ok = false; p.n = ncore;
    if (ncore >= 30.0) {
      TF1 g("g_sl_fig", "gaus", pk - 0.8, pk + 0.8);
      g.SetParameters(py->GetBinContent(pb), pk, 0.4);
      g.SetParLimits(1, pk - 0.8, pk + 0.8); g.SetParLimits(2, 0.5 * ybw, 1.5);
      py->Fit(&g, "0QRB", "", pk - 0.8, pk + 0.8);
      p.y = g.GetParameter(1); p.ey = std::max(g.GetParError(1), 1e-3); p.s = g.GetParameter(2); p.ok = true;
    }
    delete py;
    return p.ok;
  };
  Int_t seed = -1; Double_t seed_n = -1;
  for (Int_t isl = 0; isl < n_sl; ++isl) {
    const Double_t xc = h->GetXaxis()->GetBinUpEdge(isl * slice + 1);
    if (xc < 0.8 || xc > 1.3) continue;
    const Double_t n = h->Integral(isl * slice + 1, isl * slice + slice, 1, h->GetNbinsY());
    if (n > seed_n) { seed_n = n; seed = isl; }
  }
  PhcSteps out{de_min, de_max, 0, 0};
  if (seed < 0 || !slice_fit(seed, 0, -1)) return out;
  for (Int_t dir : {+1, -1}) {
    Double_t yp = pts[seed].y; Int_t nf = 0;
    for (Int_t isl = seed + dir; isl >= 0 && isl < n_sl; isl += dir) {
      if (slice_fit(isl, yp, 1.0)) { yp = pts[isl].y; nf = 0; } else if (++nf >= 3) break;
    }
  }
  Double_t nmax = 0; for (const auto& p : pts) if (p.ok) nmax = std::max(nmax, p.n);
  auto dense = [&](Int_t i) { return pts[i].ok && pts[i].n >= dense_ratio * nmax; };
  Int_t il = seed, ir = seed;
  while (il - 1 >= 0 && dense(il - 1)) --il;
  while (ir + 1 < n_sl && dense(ir + 1)) ++ir;
  const Double_t flo = std::max(pts[il].xl, de_min), fhi = std::min(pts[ir].xh, de_max);
  std::vector<Double_t> sig; for (Int_t i = il; i <= ir; ++i) if (pts[i].ok) sig.push_back(pts[i].s);
  std::sort(sig.begin(), sig.end());
  const Double_t w = std::max(2.0 * (sig.empty() ? 0.4 : sig[sig.size() / 2]), 3.0 * ybw);
  const char* form = "-[0]/TMath::Sqrt(TMath::Abs(x-[1]))+[2]";
  TF1* f = new TF1(Form("fphc_%s", tag), form, flo, fhi);
  const Double_t lim[3][2] = {{0.001, 15.0}, {-5.0, flo}, {-10.0, 10.0}};
  for (Int_t i = 0; i < 3; ++i) { f->SetParameter(i, (3 * lim[i][0] + lim[i][1]) / 4); f->SetParLimits(i, lim[i][0], lim[i][1]); }
  f->SetParameter(1, flo - (fhi - flo) * 0.005);
  for (Int_t i = 0; i < 3; ++i) f->SetParameter(i, std::min(std::max(ref[i], lim[i][0]), lim[i][1]));
  f->FixParameter(1, ref[1]);
  TGraphErrors* gr = new TGraphErrors();      // ridge points inside the fit range
  TGraphErrors* gr_all = new TGraphErrors();  // all ridge points
  for (const auto& p : pts) {
    if (!p.ok) continue;
    gr_all->SetPoint(gr_all->GetN(), p.x, p.y); gr_all->SetPointError(gr_all->GetN() - 1, 0, p.ey);
    if (p.x < flo || p.x > fhi) continue;
    gr->SetPoint(gr->GetN(), p.x, p.y); gr->SetPointError(gr->GetN() - 1, 0, p.ey);
  }
  if (gr->GetN() >= 4) gr->Fit(f, "0QR", "", flo, fhi);
  TF1* fpre = (TF1*)f->Clone(Form("fpre_%s", tag));
  TProfile* pf = new TProfile(Form("pf_%s", tag), "", nx, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  pf->SetDirectory(nullptr);
  for (Int_t bx = h->GetXaxis()->FindBin(flo); bx <= h->GetXaxis()->FindBin(fhi); ++bx) {
    const Double_t xc = h->GetXaxis()->GetBinCenter(bx), yr = f->Eval(xc);
    for (Int_t by = 1; by <= h->GetNbinsY(); ++by) {
      const Double_t yc = h->GetYaxis()->GetBinCenter(by), cnt = h->GetBinContent(bx, by);
      if (cnt > 0 && std::fabs(yc - yr) < w) pf->Fill(xc, yc, cnt);
    }
  }
  pf->Fit(f, "0QEMR", "", flo, fhi);
  out = {flo, fhi, w, f->GetParameter(0)};

  const Double_t yc0 = pts[seed].y;
  const Double_t y0 = std::max(ylo, yc0 - 6.0), y1 = std::min(yhi, yc0 + 4.0);
  // (a) ridge
  c.cd(pad0); gPad->SetLogz();
  TH2D* h1 = (TH2D*)h->Clone(Form("%s_a", tag)); h1->SetDirectory(nullptr);
  h1->SetTitle(Form("%s  (1) ridge;#DeltaE [MIP];Y [ns]", tag));
  h1->GetXaxis()->SetRangeUser(0, 4); h1->GetYaxis()->SetRangeUser(y0, y1); h1->Draw("colz");
  gr_all->SetMarkerStyle(20); gr_all->SetMarkerSize(0.6); gr_all->SetLineWidth(1); gr_all->Draw("P same");
  TGraph* gseed = new TGraph(1); gseed->SetPoint(0, pts[seed].x, pts[seed].y);
  gseed->SetMarkerStyle(29); gseed->SetMarkerSize(3.0); gseed->SetMarkerColor(kMagenta); gseed->Draw("P same");
  // (b) fit range
  c.cd(pad0 + 1);
  TGraph* gn = new TGraph();
  Double_t nmin = dense_ratio * nmax;
  for (const auto& p : pts) if (p.n > 0) { gn->SetPoint(gn->GetN(), p.x, p.n); nmin = std::min(nmin, p.n); }
  gPad->SetLogy();
  TH1F* fr = gPad->DrawFrame(0, 0.5 * nmin, 4, 2.0 * nmax,
                             Form("%s  (2) fit range;#DeltaE [MIP];entries within #pm1 ns of the ridge", tag));
  fr->GetYaxis()->SetTitleOffset(1.3);
  TBox* bx = new TBox(flo, 0.5 * nmin, fhi, 2.0 * nmax); bx->SetFillColorAlpha(kGreen + 1, 0.25); bx->Draw();
  gn->SetMarkerStyle(20); gn->SetMarkerSize(0.8); gn->Draw("P same");
  TLine* lt = new TLine(0, dense_ratio * nmax, 4, dense_ratio * nmax); lt->SetLineColor(kRed); lt->SetLineStyle(2); lt->SetLineWidth(2); lt->Draw();
  TLine* lm = new TLine(de_min, 0.5 * nmin, de_min, 2.0 * nmax); lm->SetLineColor(kGray + 2); lm->SetLineStyle(3); lm->SetLineWidth(2); lm->Draw();
  label(0.45, 0.84, Form("fit range %.2f-%.2f MIP", flo, fhi), 0.05);
  // (c) band, profile, fit
  c.cd(pad0 + 2); gPad->SetLogz();
  TH2D* h3 = (TH2D*)h->Clone(Form("%s_c", tag)); h3->SetDirectory(nullptr);
  h3->SetTitle(Form("%s  (3) band and fit;#DeltaE [MIP];Y [ns]", tag));
  h3->GetXaxis()->SetRangeUser(0, 4); h3->GetYaxis()->SetRangeUser(y0, y1); h3->Draw("colz");
  for (Double_t sgn : {+1.0, -1.0}) {
    TF1* fb = new TF1(Form("fb_%s_%+d", tag, (Int_t)sgn), Form("%s%+f", form, sgn * w), flo, fhi);
    for (Int_t i = 0; i < 3; ++i) fb->SetParameter(i, fpre->GetParameter(i));
    fb->SetLineColor(kGray + 2); fb->SetLineStyle(2); fb->SetLineWidth(1); fb->SetNpx(300); fb->Draw("same");
  }
  pf->SetMarkerStyle(20); pf->SetMarkerSize(0.5); pf->SetMarkerColor(kRed); pf->SetLineColor(kRed); pf->Draw("same");
  f->SetLineColor(kOrange + 1); f->SetLineWidth(2); f->SetNpx(300); f->Draw("same");
  for (Double_t x : {flo, fhi}) { TLine* l = new TLine(x, y0, x, y1); l->SetLineStyle(3); l->Draw(); }
  label(0.40, 0.84, Form("w = %.2f ns", w), 0.05);
  return out;
}
} // namespace

void make_phcfit_fig(const char* calib_dir = "/home/had/sryuta/JPARC2025E72/results/root/run02601")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetTitleFontSize(0.06);
  TFile f(Form("%s/run02601_HTOF_PHC_Pi.root", calib_dir));
  TTree* t = (TTree*)f.Get("tree");
  Int_t ch; std::vector<Double_t> *p0a = nullptr, *p1a = nullptr, *p2a = nullptr, *p0 = nullptr, *flo = nullptr, *fhi = nullptr;
  t->SetBranchAddress("ch", &ch); t->SetBranchAddress("p0_all_val", &p0a); t->SetBranchAddress("p1_all_val", &p1a);
  t->SetBranchAddress("p2_all_val", &p2a); t->SetBranchAddress("p0_val", &p0);
  t->SetBranchAddress("fit_lo", &flo); t->SetBranchAddress("fit_hi", &fhi);
  TCanvas c("c_phcfit", "", 1800, 1100);
  c.Divide(3, 2);
  Int_t row = 0;
  for (Int_t seg : {20, 7}) {
    t->GetEntry(seg);
    const Double_t ref[3] = {(*p0a)[0], (*p1a)[0], (*p2a)[0]};
    TH2D* h = (TH2D*)f.Get(Form("HTOF_seg%dU_PHC_Y_sel", seg));
    if (!h) continue;
    const PhcSteps r = draw_phc_steps(c, 3 * row + 1, h, ref, 0.04, Form("seg%d U", seg));
    printf("seg%dU: fit range %.3f-%.3f (stored %.3f-%.3f), w %.3f, p0 %.4f (stored %.4f)\n",
           seg, r.fit_lo, r.fit_hi, (*flo)[0], (*fhi)[0], r.w, r.p0, (*p0)[0]);
    ++row;
  }
  c.Print("phc_fit_steps.pdf");
}

void make_figs(const char* data_dir = "/group/had/sks/Users/sryuta/JPARC2025E72",
               const char* calib_dir = "/home/had/sryuta/JPARC2025E72/results/root/run02601")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetTitleFontSize(0.05);

  // ---- TDC: seg20 U, all events, main peak with a gauss fit
  {
    TH1D* h = sum_hist(data_dir, "HTOF_TDC_seg20U");
    TCanvas c("c_tdc", "", 900, 600);
    if (h) {
      h->SetTitle("HTOF seg20 U, TDC (all events, 6 runs);TDC [count];entries");
      const Double_t pk = h->GetBinCenter(h->GetMaximumBin());
      h->GetXaxis()->SetRangeUser(pk - 20000, pk + 20000);
      TF1 g("g_tdc", "gaus", pk - 3000, pk + 3000);
      h->Fit(&g, "0QR");
      h->Draw("HIST"); g.SetLineColor(kRed); g.DrawCopy("same");
      label(0.58, 0.82, Form("offset = %.0f count", g.GetParameter(1)), 0.04);
      label(0.58, 0.76, "(mean of the gauss fit)", 0.04);
      label(0.58, 0.70, "gain = -9.765625#times10^{-4} ns/count", 0.035);
    }
    c.Print("tdc_example.pdf");
  }

  // ---- ADC: seg20 U, raw (all events) and pion-tagged (from the DST)
  {
    TH1D* raw = sum_hist(data_dir, "HTOF_ADC_seg20U");
    TH1D* sel = new TH1D("sel20U", "HTOF seg20 U, ADC of pion-tagged matched tracks;ADC;entries", 1024, 0, 4096);
    sel->SetDirectory(nullptr);
    for (const char* r : kRuns) {
      // pion bit set, proton bit not set, matched; vertex tracks + no-vertex tracks with dseg = 0
      TFile f(Form("%s/run%s_TPCHelixHTOF.root", data_dir, r));
      TTreeReader rd("htof", &f);
      TTreeReaderValue<std::vector<Int_t>> ok(rd, "seg_match"), pid(rd, "pid"), vs(rd, "vertex_source");
      TTreeReaderValue<std::vector<Double_t>> seg(rd, "extrap_seg"), cls(rd, "match_cl_seg"), adc(rd, "htof_adc_u");
      while (rd.Next()) for (size_t i = 0; i < ok->size(); i++) {
        if ((*ok)[i] != 1 || !((*pid)[i] & 1) || ((*pid)[i] & 4) || TMath::Nint((*seg)[i]) != 20) continue;
        if ((*vs)[i] == 0 && TMath::Abs((*seg)[i] - (*cls)[i]) > 0.5) continue;
        if (!TMath::IsNaN((*adc)[i])) sel->Fill((*adc)[i]);
      }
    }
    TCanvas c("c_adc", "", 1500, 500);
    c.Divide(3, 1);
    c.cd(1); gPad->SetLogy();
    Double_t ped = 0;
    if (raw) {
      raw->SetTitle("seg20 U, ADC (all events): pedestal;ADC;entries");
      raw->GetXaxis()->SetRangeUser(0, 1200);
      ped = raw->GetBinCenter(raw->GetMaximumBin());
      TF1 gp("gp", "gaus", ped - 5, ped + 5); raw->Fit(&gp, "0QR");
      TF1 gp2("gp2", "gaus", gp.GetParameter(1) - 2 * gp.GetParameter(2), gp.GetParameter(1) + 2 * gp.GetParameter(2));
      raw->Fit(&gp2, "0QR"); ped = gp2.GetParameter(1);
      raw->Draw("HIST");
      TLine* l = new TLine(ped, 1, ped, raw->GetMaximum()); l->SetLineColor(kMagenta); l->SetLineStyle(2); l->Draw();
      label(0.40, 0.82, Form("pedestal P = %.1f", ped));
    }
    c.cd(2); gPad->SetLogy();
    sel->GetXaxis()->SetRangeUser(0, 1200); sel->Draw("HIST");
    TLine* l2 = new TLine(ped, 1, ped, sel->GetMaximum()); l2->SetLineColor(kMagenta); l2->SetLineStyle(2); l2->Draw();
    label(0.35, 0.82, "pedestal-height entries remain");
    c.cd(3);
    TH1D* z = (TH1D*)sel->Rebin(2, "sel20U_z");
    z->SetTitle("seg20 U, pion-tagged (close-up) with the MIP fit;ADC;entries");
    z->GetXaxis()->SetRangeUser(ped + 20, ped + 400);
    const Double_t mpk = z->GetBinCenter(z->GetMaximumBin());
    const Double_t w = 0.35 * (mpk - ped);
    TF1 fl("fl", "landaun(0)+pol0(3)", mpk - 0.6 * w, mpk + 1.2 * w);
    fl.SetParameters(z->Integral(z->FindBin(mpk - w), z->FindBin(mpk + w), "width"), mpk, 0.3 * w, 0);
    fl.SetParLimits(3, 0, z->GetMaximum());
    z->Fit(&fl, "0QR");
    z->Draw("HIST"); fl.SetLineColor(kRed); fl.DrawCopy("same");
    label(0.45, 0.82, Form("MIP M = %.1f", fl.GetMaximumX(mpk - 0.6 * w, mpk + 1.2 * w)));
    c.Print("adc_example.pdf");
  }

  // ---- weak side (beam-window short segment): seg1 D, tracks with strong side near its MIP
  {
    TFile f(Form("%s/run02601_HTOF_ADC_Pi.root", calib_dir));
    TH1D* hc = (TH1D*)f.Get("HTOF_ADC_seg1D_Pi_cond");
    TH1D* hs = (TH1D*)f.Get("HTOF_ADC_seg1D_Pi");
    TCanvas c("c_weak", "", 1400, 500);
    c.Divide(2, 1);
    c.cd(1);
    if (hs) {
      TH1D* a = (TH1D*)hs->Rebin(4, "w_all"); a->SetTitle("seg1 D (weak side), pion-tagged, all;ADC;entries");
      a->GetXaxis()->SetRangeUser(100, 900); a->Draw("HIST");
      label(0.40, 0.82, "spike just above the pedestal");
    }
    c.cd(2);
    if (hc) {
      TH1D* b = (TH1D*)hc->Rebin(4, "w_cond"); b->SetTitle("seg1 D, strong side (U) within #pm30% of its MIP;ADC;entries");
      b->GetXaxis()->SetRangeUser(100, 900);
      const Double_t ped = 171.1, u = 355.0 - 132.6;
      TF1 f2("f2", "gaus(0)+landau(3)+[6]", ped + 5, ped + 3 * u);
      b->GetXaxis()->SetRangeUser(ped + 5, ped + 80); const Double_t sp = b->GetBinCenter(b->GetMaximumBin()), spa = b->GetMaximum();
      b->GetXaxis()->SetRangeUser(100, 900);
      f2.SetParameters(spa, sp, 12, 0.3 * spa, ped + 0.6 * u, 0.15 * u, 0);
      f2.SetParLimits(1, ped, ped + 80); f2.SetParLimits(2, 2, 40); f2.SetParLimits(4, sp + 30, ped + 2 * u); f2.SetParLimits(5, 5, u); f2.SetParLimits(6, 0, 0.05 * spa);
      b->Fit(&f2, "0QRL");
      b->Draw("HIST"); f2.SetLineColor(kOrange + 1); f2.DrawCopy("same");
      TF1 fs("fs", "gaus", ped + 5, ped + 3 * u); fs.SetParameters(f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2)); fs.SetLineColor(kBlue); fs.SetLineStyle(2); fs.SetNpx(500); fs.DrawCopy("same");
      TF1 fh("fh", "landau", ped + 5, ped + 3 * u); fh.SetParameters(f2.GetParameter(3), f2.GetParameter(4), f2.GetParameter(5)); fh.SetLineColor(kGreen + 2); fh.SetLineStyle(2); fh.SetNpx(500); fh.DrawCopy("same");
      label(0.38, 0.82, Form("spike (blue) + hump (green): MIP = %.0f", fh.GetMaximumX(ped + 5, ped + 3 * u)));
    }
    c.Print("weak_example.pdf");
  }

  // ---- PHC: seg20 U, Y_U vs dE_U before and (t_c mean - time0 - expected) after
  {
    TFile f(Form("%s/run02601_HTOF_PHC_Pi.root", calib_dir));
    TH2D* y = (TH2D*)f.Get("HTOF_seg20U_PHC_Y_sel");
    TH2D* a = (TH2D*)f.Get("HTOF_seg20_dTPi_vs_DeltaEU_calib");
    TCanvas c("c_phc", "", 1400, 550);
    c.Divide(2, 1);
    c.cd(1); gPad->SetLogz();
    if (y) { y->SetTitle("seg20 U: Y_{U} = (time0 - t_{U}) + (time0 - t_{c,D}) vs #DeltaE_{U};#DeltaE_{U} [MIP];Y_{U} [ns]"); y->GetXaxis()->SetRangeUser(0, 4); y->Draw("colz"); }
    c.cd(2); gPad->SetLogz();
    if (a) { a->SetTitle("seg20: dt_{#pi} after the correction vs #DeltaE_{U};#DeltaE_{U} [MIP];dt_{#pi} [ns]"); a->GetXaxis()->SetRangeUser(0, 4); a->Draw("colz"); TLine* l = new TLine(0, 0, 4, 0); l->SetLineStyle(2); l->Draw(); }
    c.Print("phc_example.pdf");
  }

  // ---- QA: dt_pi and dt_pi per segment (pions with p_vtx >= 0.20 GeV/c, the absolute-offset sample), m2 vs p (q > 0)
  {
    TH1D* hd = new TH1D("qa_dtpi", "dt_{#pi}, pion-tagged, p_{vtx} #geq 0.20 GeV/c (non-beam, TOF valid);dt_{#pi} [ns];tracks", 400, -5, 5);
    TH2D* hm = new TH2D("qa_m2p", "m^{2} vs p, q > 0 (non-beam, TOF valid);p_{vtx} [GeV/c];m^{2} [GeV^{2}]", 150, 0, 1.5, 250, -0.5, 2.0);
    TH2D* hs = new TH2D("qa_seg", "", 34, -0.5, 33.5, 200, -3, 3);
    for (const char* r : kRuns) {
      TFile f(Form("%s/run%s_TPCHelixHTOF.root", data_dir, r));
      TTreeReader rd("htof", &f);
      TTreeReaderValue<std::vector<Int_t>> tv(rd, "tof_valid"), isb(rd, "is_beam"), isa(rd, "is_accidental"), pid(rd, "pid"), q(rd, "charge");
      TTreeReaderValue<std::vector<Double_t>> dt(rd, "dt_pi"), m2(rd, "m2"), pv(rd, "p_vtx"), seg(rd, "extrap_seg");
      while (rd.Next()) for (size_t i = 0; i < tv->size(); i++) {
        if ((*tv)[i] != 1 || (*isb)[i] || (*isa)[i]) continue;
        const Bool_t is_pi = ((*pid)[i] & 1) && !((*pid)[i] & 4);
        if (is_pi && (*pv)[i] >= 0.20) { hd->Fill((*dt)[i]); hs->Fill((*seg)[i], (*dt)[i]); }
        if ((*q)[i] > 0) hm->Fill((*pv)[i], (*m2)[i]);
      }
    }
    TCanvas c("c_qa", "", 1500, 500);
    c.Divide(3, 1);
    c.cd(1);
    const Double_t pk = hd->GetBinCenter(hd->GetMaximumBin());
    TF1 g("g_qa", "gaus", pk - 0.6, pk + 0.6); hd->Fit(&g, "0QR");
    hd->Draw("HIST"); g.SetLineColor(kRed); g.DrawCopy("same");
    label(0.15, 0.82, Form("core %+.3f ns, #sigma %.3f ns", g.GetParameter(1), g.GetParameter(2)));
    c.cd(2); gPad->SetLogz(); hm->Draw("colz");
    for (Double_t m : {0.0195, 0.2437, 0.8804}) { TLine* l = new TLine(0, m, 1.5, m); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
    c.cd(3);
    TGraphErrors* gm = new TGraphErrors();
    TGraphErrors* gs = new TGraphErrors();
    for (Int_t s = 0; s < 34; s++) {
      TH1D* p = hs->ProjectionY(Form("qa_p%d", s), s + 1, s + 1);
      if (p->GetEntries() < 50) { delete p; continue; }
      const Double_t cc = p->GetBinCenter(p->GetMaximumBin());
      TF1 gg("gg", "gaus", cc - 0.6, cc + 0.6); p->Fit(&gg, "0QR");
      const Int_t k = gm->GetN();
      gm->SetPoint(k, s, gg.GetParameter(1)); gm->SetPointError(k, 0, gg.GetParError(1));
      gs->SetPoint(k, s, std::fabs(gg.GetParameter(2)));
      delete p;
    }
    gm->SetTitle("dt_{#pi} per segment, p_{vtx} #geq 0.20 GeV/c (black: core mean, red: #sigma);HTOF seg;[ns]");
    gm->SetMarkerStyle(20); gm->SetMinimum(-0.6); gm->SetMaximum(0.8); gm->Draw("AP");
    gs->SetMarkerStyle(24); gs->SetMarkerColor(kRed); gs->Draw("P same");
    TLine* l0 = new TLine(-0.5, 0, 33.5, 0); l0->SetLineStyle(2); l0->Draw();
    c.Print("qa_dtpi_m2.pdf");
  }

  make_phcfit_fig(calib_dir);
}
