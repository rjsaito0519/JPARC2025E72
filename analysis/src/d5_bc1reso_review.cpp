// BLC1a/b offsets against BLC2 from the BLC2-referenced D5 fit, with beam-flag selection
// Usage: d5_bc1reso_review <run> [--pion | --kaon | --all] [--d5 <D5_bc1reso.root>]
//                          [--dcgeo <DCGEO>] [-o <out.pdf>] [--txt <out.txt>]
// Defaults: D5   DECODE_DIR/runXXXXX/runXXXXX_D5_K_bc1reso.root  (BLC1 Res inflated)
//           geo  ANALYZER_DIR/param/DCGEO/e72/DCGeomParam_runXXXXX_K  (BLC1a/b z)
//           pdf  OUTPUT_DIR/img/runXXXXX/d5_bc1reso_review_runXXXXX.pdf
//           txt  kTxtDir/runXXXXX/d5_bc1reso_review_runXXXXX.txt  (numbers for run_d5.py; not in img)
//           beam flag: --kaon (beam_flag == 2)
//
// With the BLC1 Res inflated the D5 fit line at BLC1 (d5_fit_*0) is set by BLC2 alone.
// For each rank0 pair the offset r(z) = fit line - BLC1 track is taken at the mean z of
// BLC1a and of BLC1b, and the core mean of r is the --residual-offset of that chamber
// (same sign convention as the LS result of d5_wire_resi_review). The per-plane
// histograms of d5_wire_resi_review mix all beam flags; this tool selects them.
// The common x shift of BLC1 is degenerate with delta: x is shown for information only.

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

const TString kTxtDir = "/group/had/sks/Users/sryuta/tmp/d5_align";
constexpr Int_t kBeamFlagPion = 1;
constexpr Int_t kBeamFlagKaon = 2;
constexpr Double_t kCoreWindow = 1.5; // +-sigma window of the iterated Gauss core
constexpr Int_t kCoreIter = 6;

struct GaussResult
{
  Double_t mean = NAN;
  Double_t mean_err = NAN;
  Double_t sigma = NAN;
};

// Mean Z of the BLC1a / BLC1b layers from DCGEO (column Z).
Bool_t
load_blc1_z(const TString& path, Double_t& z_a, Double_t& z_b)
{
  std::ifstream ifs(path.Data());
  if (!ifs.is_open()) return kFALSE;
  Double_t sum_a = 0.0, sum_b = 0.0;
  Int_t n_a = 0, n_b = 0;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string id, name;
    Double_t x, y, z;
    if (!(iss >> id >> name >> x >> y >> z)) continue;
    if (name.rfind("BLC1a-", 0) == 0) { sum_a += z; ++n_a; }
    if (name.rfind("BLC1b-", 0) == 0) { sum_b += z; ++n_b; }
  }
  if (n_a == 0 || n_b == 0) return kFALSE;
  z_a = sum_a / n_a;
  z_b = sum_b / n_b;
  return kTRUE;
}

// Gauss fit in +-kCoreWindow sigma, iterated from the histogram peak.
// r(z) is wide (multiple scattering between BLC1 and BLC2) with non-Gaussian tails.
GaussResult
fit_core(TH1D* h, Double_t sigma0)
{
  GaussResult g;
  if (!h || h->GetEntries() < 100) return g;
  Double_t m = h->GetBinCenter(h->GetMaximumBin()), s = sigma0;
  auto* f = new TF1(Form("g_%s", h->GetName()), "gaus", m - kCoreWindow * s, m + kCoreWindow * s);
  f->SetLineColor(kRed);
  for (Int_t it = 0; it < kCoreIter; ++it) {
    f->SetRange(m - kCoreWindow * s, m + kCoreWindow * s);
    f->SetParameters(h->GetMaximum(), m, s);
    h->Fit(f, "QR0");
    m = f->GetParameter(1);
    s = std::fabs(f->GetParameter(2));
  }
  g.mean = m;
  g.mean_err = f->GetParError(1);
  g.sigma = s;
  return g;
}

void
draw_with_fit(TH1* h)
{
  h->Draw();
  if (auto* f = h->GetFunction(Form("g_%s", h->GetName()))) {
    f->ResetBit(TF1::kNotDraw);
    f->Draw("same");
  }
}

TString
flag_label(Int_t flag)
{
  if (flag == kBeamFlagPion) return "pi (beam_flag=1)";
  if (flag == kBeamFlagKaon) return "K window (beam_flag=2)";
  return "all beam_flag";
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0
            << " <run> [--pion | --kaon | --all] [--d5 <D5_bc1reso.root>]"
            << " [--dcgeo <file>] [-o <out.pdf>] [--txt <out.txt>]\n";
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  if (argc < 2) { usage(argv[0]); return 1; }
  const Int_t run = std::atoi(argv[1]);
  Int_t flag = kBeamFlagKaon;
  TString path_d5 = Form("%s/run%05d/run%05d_D5_K_bc1reso.root", DECODE_DIR.Data(), run, run);
  TString path_geo = Form("%s/param/DCGEO/e72/DCGeomParam_run%05d_K", ANALYZER_DIR.Data(), run);
  TString path_pdf = Form("%s/d5_bc1reso_review_run%05d.pdf",
                          ana_helper::get_img_dir(OUTPUT_DIR, run).Data(), run);
  TString path_txt;
  for (Int_t i = 2; i < argc; ++i) {
    const TString a = argv[i];
    if (a == "--pion") flag = kBeamFlagPion;
    else if (a == "--kaon") flag = kBeamFlagKaon;
    else if (a == "--all") flag = -1;
    else if (a == "--d5" && i + 1 < argc) path_d5 = argv[++i];
    else if (a == "--dcgeo" && i + 1 < argc) path_geo = argv[++i];
    else if (a == "-o" && i + 1 < argc) path_pdf = argv[++i];
    else if (a == "--txt" && i + 1 < argc) path_txt = argv[++i];
    else { usage(argv[0]); return 1; }
  }
  if (path_txt.IsNull()) {
    const TString dir = Form("%s/run%05d", kTxtDir.Data(), run);
    gSystem->mkdir(dir, kTRUE);
    path_txt = Form("%s/%s", dir.Data(), gSystem->BaseName(path_pdf));
    path_txt.ReplaceAll(".pdf", ".txt");
  }

  Double_t z_a = 0.0, z_b = 0.0;
  if (!load_blc1_z(path_geo, z_a, z_b)) {
    std::cerr << "cannot read BLC1 z from: " << path_geo << std::endl;
    return 1;
  }

  auto* f = TFile::Open(path_d5);
  if (!f || f->IsZombie()) { std::cerr << "cannot open " << path_d5 << std::endl; return 1; }
  auto* t = dynamic_cast<TTree*>(f->Get("d5"));
  if (!t) { std::cerr << "tree d5 not found" << std::endl; return 1; }

  Int_t beam_flag = 0;
  std::vector<Double_t> *x1 = nullptr, *u1 = nullptr, *y1 = nullptr, *v1 = nullptr;
  std::vector<Double_t> *fx = nullptr, *fu = nullptr, *fy = nullptr, *fv = nullptr;
  std::vector<Double_t>* d5_delta = nullptr;
  std::vector<Int_t> *rank = nullptr, *trk1 = nullptr;
  t->SetBranchAddress("beam_flag", &beam_flag);
  t->SetBranchAddress("blc1_x", &x1); t->SetBranchAddress("blc1_u", &u1);
  t->SetBranchAddress("blc1_y", &y1); t->SetBranchAddress("blc1_v", &v1);
  t->SetBranchAddress("d5_fit_x0", &fx); t->SetBranchAddress("d5_fit_u0", &fu);
  t->SetBranchAddress("d5_fit_y0", &fy); t->SetBranchAddress("d5_fit_v0", &fv);
  t->SetBranchAddress("d5_delta", &d5_delta);
  t->SetBranchAddress("d5_pair_rank", &rank);
  t->SetBranchAddress("d5_blc1_trk", &trk1);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  auto* h_ya = new TH1D("h_ya", "BLC1a: r_{y} = y_{fit} - y_{BLC1}; r_{y} [mm]; count", 600, -30, 30);
  auto* h_yb = new TH1D("h_yb", "BLC1b: r_{y} = y_{fit} - y_{BLC1}; r_{y} [mm]; count", 600, -30, 30);
  auto* h_xa = new TH1D("h_xa", "BLC1a: r_{x} = x_{fit} - x_{BLC1} (info); r_{x} [mm]; count", 400, -4, 4);
  auto* h_xb = new TH1D("h_xb", "BLC1b: r_{x} = x_{fit} - x_{BLC1} (info); r_{x} [mm]; count", 400, -4, 4);
  auto* p_ya_delta = new TProfile("p_ya_delta", "BLC1a r_{y} vs #delta (|r_{y}|<10); d5_delta [%]; r_{y} [mm]", 30, -5, 10);
  auto* p_yb_delta = new TProfile("p_yb_delta", "BLC1b r_{y} vs #delta (|r_{y}|<10); d5_delta [%]; r_{y} [mm]", 30, -5, 10);

  Long64_t n_used = 0;
  for (Long64_t e = 0, ne = t->GetEntries(); e < ne; ++e) {
    t->GetEntry(e);
    if (flag >= 0 && beam_flag != flag) continue;
    for (std::size_t k = 0; k < rank->size(); ++k) {
      if ((*rank)[k] != 0) continue;
      const Int_t a = (*trk1)[k];
      const Double_t delta = (*d5_delta)[k];
      const Double_t rya = ((*fy)[k] + (*fv)[k] * z_a) - ((*y1)[a] + (*v1)[a] * z_a);
      const Double_t ryb = ((*fy)[k] + (*fv)[k] * z_b) - ((*y1)[a] + (*v1)[a] * z_b);
      const Double_t rxa = ((*fx)[k] + (*fu)[k] * z_a) - ((*x1)[a] + (*u1)[a] * z_a);
      const Double_t rxb = ((*fx)[k] + (*fu)[k] * z_b) - ((*x1)[a] + (*u1)[a] * z_b);
      h_ya->Fill(rya);
      h_yb->Fill(ryb);
      h_xa->Fill(rxa);
      h_xb->Fill(rxb);
      if (std::fabs(rya) < 10.0) p_ya_delta->Fill(delta, rya);
      if (std::fabs(ryb) < 10.0) p_yb_delta->Fill(delta, ryb);
      ++n_used;
    }
  }

  const GaussResult g_ya = fit_core(h_ya, h_ya->GetRMS());
  const GaussResult g_yb = fit_core(h_yb, h_yb->GetRMS());
  const GaussResult g_xa = fit_core(h_xa, 0.1);
  const GaussResult g_xb = fit_core(h_xb, 0.1);

  std::ofstream ofs(path_txt.Data());
  ofs << "run " << run << "\n"
      << "beam_flag " << flag_label(flag) << "\n"
      << "d5_root " << path_d5 << "\n"
      << "dcgeo " << path_geo << "\n"
      << "n_pairs " << n_used << "\n"
      << Form("blc1a_z %.3f\nblc1b_z %.3f\n", z_a, z_b)
      << Form("dy_a_mm %.5f\ndy_a_err_mm %.5f\ndy_a_sigma_mm %.4f\n", g_ya.mean, g_ya.mean_err, g_ya.sigma)
      << Form("dy_b_mm %.5f\ndy_b_err_mm %.5f\ndy_b_sigma_mm %.4f\n", g_yb.mean, g_yb.mean_err, g_yb.sigma)
      << Form("dx_a_mm_info %.5f\ndx_b_mm_info %.5f\n", g_xa.mean, g_xb.mean)
      << Form("apply_y --residual-offset BLC1a:0,%.4f --residual-offset BLC1b:0,%.4f\n", g_ya.mean, g_yb.mean);
  ofs.close();

  auto* c = new TCanvas("c", "c", 1200, 900);
  c->Print(path_pdf + "[");
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  Double_t yy = 0.92;
  auto line = [&](const TString& s) { tex->DrawLatex(0.06, yy, s); yy -= 0.045; };
  line(Form("D5 BLC2-referenced review: run %05d, %s, %lld pairs (rank0)", run, flag_label(flag).Data(), n_used));
  line(Form("D5 ROOT: %s", path_d5.Data()));
  line(Form("DCGEO: %s  (BLC1a z=%.1f, BLC1b z=%.1f mm)", path_geo.Data(), z_a, z_b));
  line("r(z) = D5 fit line (BLC1 Res inflated: set by BLC2) - BLC1 track, core Gauss mean");
  line(Form("  BLC1a dy = %+.4f #pm %.4f mm  (#sigma %.2f mm)", g_ya.mean, g_ya.mean_err, g_ya.sigma));
  line(Form("  BLC1b dy = %+.4f #pm %.4f mm  (#sigma %.2f mm)", g_yb.mean, g_yb.mean_err, g_yb.sigma));
  line(Form("  #rightarrow --residual-offset BLC1a:0,%.4f --residual-offset BLC1b:0,%.4f", g_ya.mean, g_yb.mean));
  line(Form("  x (info, common x degenerate with #delta): BLC1a %+.4f, BLC1b %+.4f mm", g_xa.mean, g_xb.mean));
  c->Print(path_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_ya);
  c->cd(2); draw_with_fit(h_yb);
  c->cd(3); p_ya_delta->Draw();
  c->cd(4); p_yb_delta->Draw();
  c->Print(path_pdf);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); draw_with_fit(h_xa);
  c->cd(2); draw_with_fit(h_xb);
  c->Print(path_pdf);
  c->Print(path_pdf + "]");

  std::cout << "wrote " << path_pdf << "\n      " << path_txt << std::endl;
  std::cout << Form("dy_a=%+.4f+-%.4f dy_b=%+.4f+-%.4f mm (%s)", g_ya.mean, g_ya.mean_err,
                    g_yb.mean, g_yb.mean_err, flag_label(flag).Data()) << std::endl;
  return 0;
}
