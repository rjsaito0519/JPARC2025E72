// D5 fit-free consistency review (BLC1 track transported by the D5 matrix vs BLC2 track)
// Usage: d5_fitfree_review <run> [--pion | --kaon | --all] [--d5 <D5.root>]
//                          [--mtx <D5TransferMatrix>] [--dcgeo <DCGEO>] [-o <out.pdf>]
//                          [--txt <out.txt>]
// Defaults: D5   DECODE_DIR/runXXXXX/runXXXXX_D5_K.root
//           mtx  ANALYZER_DIR/param/TM/e72/D5TransferMatrix_runXXXXX
//           geo  ANALYZER_DIR/param/DCGEO/e72/DCGeomParam_runXXXXX_K
//           pdf  OUTPUT_DIR/img/runXXXXX/d5_fitfree_review_runXXXXX.pdf
//           txt  kTxtDir/runXXXXX/d5_fitfree_review_runXXXXX.txt  (numbers for run_d5align.py; not in img)
//           beam flag: --kaon (beam_flag == 2)
//
// Horizontal: for each rank0 pair, solve delta so that the matrix reproduces
// the BLC2 x at D5ZOut (delta_x) and the BLC2 slope u (delta_u) from the BLC1
// track alone. <delta_x - delta_u> = 0 for consistent geometry and matrix.
// An opposite x shift of BLC1a (+s) and BLC1b (-s) changes it by k*s with
//   k = (-R12/R16 + R22/R26) * (-2000 / (z_b - z_a))   [%/mm]
// so s = -<delta_x - delta_u> / k is the BLC1a/b correction
// (apply_dcgeo_offset.py --residual-offset BLC1a:+s,0 --residual-offset BLC1b:-s,0).
// The common x shift of BLC1 is degenerate with delta and is not determined here.
//
// Vertical: dy = y_BLC2(D5ZOut) - y_mtx(BLC1), dv likewise [mrad]; information only.

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
constexpr Double_t kZExtrapBLC1 = 0.0;
constexpr Double_t kDeltaFitLo = -2.0; // [%] range for linear fits vs delta
constexpr Double_t kDeltaFitHi = 6.0;

struct D5Matrix
{
  Double_t r[6][6] = {};
  Double_t t[5][6][6] = {};
  Double_t z_in = 0.0;
  Double_t z_out = 0.0;
  Bool_t ready = false;
};

struct GaussResult
{
  Double_t mean = 0.0;
  Double_t mean_err = 0.0;
  Double_t sigma = 0.0;
};

Bool_t
load_matrix(const TString& path, D5Matrix& m)
{
  std::ifstream ifs(path.Data());
  if (!ifs.is_open()) return kFALSE;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string key;
    iss >> key;
    if (key == "D5ZIn:") {
      iss >> m.z_in;
    } else if (key == "D5ZOut:") {
      iss >> m.z_out;
    } else if (key == "D5Matrix:") {
      Int_t order = 0;
      iss >> order;
      if (order == 1) {
        for (Int_t i = 0; i < 6; ++i)
          for (Int_t j = 0; j < 6; ++j) ifs >> m.r[i][j];
      } else if (order >= 21 && order <= 25) {
        for (Int_t i = 0; i < 6; ++i)
          for (Int_t j = 0; j < 6; ++j) ifs >> m.t[order - 21][i][j];
      }
    }
  }
  m.ready = kTRUE;
  return kTRUE;
}

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

// in: x,y [mm] at z_in, u,v [dx/dz], delta [%]; out: x,y [mm], u,v [mrad] at z_out
void
transport(const D5Matrix& m, Double_t x, Double_t u, Double_t y, Double_t v,
          Double_t delta, Double_t* out)
{
  const Double_t p[6] = { (x + u * m.z_in) / 10.0, u * 1000.0,
                          (y + v * m.z_in) / 10.0, v * 1000.0, 0.0, delta };
  for (Int_t i = 0; i < 4; ++i) {
    Double_t s = 0.0;
    for (Int_t j = 0; j < 6; ++j) s += m.r[i][j] * p[j];
    for (Int_t j = 0; j < 6; ++j)
      for (Int_t k = 0; k < 6; ++k) s += m.t[i][j][k] * p[j] * p[k];
    out[i] = s;
  }
  out[0] *= 10.0;
  out[2] *= 10.0;
}

// delta [%] such that transported out[idx] equals target (Newton)
Double_t
solve_delta(const D5Matrix& m, Double_t x, Double_t u, Double_t y, Double_t v,
            Int_t idx, Double_t target)
{
  Double_t d = 0.0;
  for (Int_t it = 0; it < 20; ++it) {
    Double_t o0[4], o1[4];
    transport(m, x, u, y, v, d, o0);
    transport(m, x, u, y, v, d + 1.0e-3, o1);
    const Double_t der = (o1[idx] - o0[idx]) / 1.0e-3;
    if (der == 0.0) return NAN;
    const Double_t step = (target - o0[idx]) / der;
    d += step;
    if (std::fabs(step) < 1.0e-7) return d;
  }
  return d;
}

GaussResult
fit_core(TH1D* h)
{
  GaussResult g;
  if (!h || h->GetEntries() < 20) return g;
  const Double_t m = h->GetMean(), r = h->GetRMS();
  auto* f = new TF1(Form("g_%s", h->GetName()), "gaus", m - 1.5 * r, m + 1.5 * r);
  f->SetLineColor(kRed);
  h->Fit(f, "QR0");
  g.mean = f->GetParameter(1);
  g.mean_err = f->GetParError(1);
  g.sigma = f->GetParameter(2);
  return g;
}

Double_t
fit_line(TProfile* p, Double_t lo, Double_t hi, Double_t* slope = nullptr)
{
  auto* f = new TF1(Form("l_%s", p->GetName()), "pol1", lo, hi);
  f->SetLineColor(kRed);
  p->Fit(f, "QR0");
  if (slope) *slope = f->GetParameter(1);
  return f->GetParameter(0);
}

void
draw_with_fit(TH1* h)
{
  h->Draw();
  auto* f = h->GetFunction(Form("g_%s", h->GetName()));
  if (!f) f = h->GetFunction(Form("l_%s", h->GetName()));
  if (f) { f->ResetBit(TF1::kNotDraw); f->Draw("same"); }
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
            << " <run> [--pion | --kaon | --all] [--d5 <D5.root>] [--mtx <file>]"
            << " [--dcgeo <file>] [-o <out.pdf>] [--txt <out.txt>]\n";
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  if (argc < 2) { usage(argv[0]); return 1; }
  const Int_t run = std::atoi(argv[1]);
  Int_t flag = kBeamFlagKaon;
  TString path_d5 = Form("%s/run%05d/run%05d_D5_K.root", DECODE_DIR.Data(), run, run);
  TString path_mtx = Form("%s/param/TM/e72/D5TransferMatrix_run%05d", ANALYZER_DIR.Data(), run);
  TString path_geo = Form("%s/param/DCGEO/e72/DCGeomParam_run%05d_K", ANALYZER_DIR.Data(), run);
  TString path_pdf = Form("%s/d5_fitfree_review_run%05d.pdf",
                          ana_helper::get_img_dir(OUTPUT_DIR, run).Data(), run);
  TString path_txt;
  for (Int_t i = 2; i < argc; ++i) {
    const TString a = argv[i];
    if (a == "--pion") flag = kBeamFlagPion;
    else if (a == "--kaon") flag = kBeamFlagKaon;
    else if (a == "--all") flag = -1;
    else if (a == "--d5" && i + 1 < argc) path_d5 = argv[++i];
    else if (a == "--mtx" && i + 1 < argc) path_mtx = argv[++i];
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

  D5Matrix mtx;
  if (!load_matrix(path_mtx, mtx)) {
    std::cerr << "cannot open matrix: " << path_mtx << std::endl;
    return 1;
  }
  Double_t z_a = 0.0, z_b = 0.0;
  if (!load_blc1_z(path_geo, z_a, z_b)) {
    std::cerr << "cannot read BLC1 z from: " << path_geo << std::endl;
    return 1;
  }
  const Double_t ddu_ddelta_slope =
    -mtx.r[0][1] / mtx.r[0][5] + mtx.r[1][1] / mtx.r[1][5]; // [%/mrad]
  const Double_t k_s = ddu_ddelta_slope * (-2000.0 / (z_b - z_a)); // [%/mm]

  auto* f = TFile::Open(path_d5);
  if (!f || f->IsZombie()) { std::cerr << "cannot open " << path_d5 << std::endl; return 1; }
  auto* t = dynamic_cast<TTree*>(f->Get("d5"));
  if (!t) { std::cerr << "tree d5 not found" << std::endl; return 1; }

  Int_t beam_flag = 0;
  std::vector<Double_t> *x1 = nullptr, *u1 = nullptr, *y1 = nullptr, *v1 = nullptr;
  std::vector<Double_t> *x2 = nullptr, *u2 = nullptr, *y2 = nullptr, *v2 = nullptr;
  std::vector<Double_t>* d5_delta = nullptr;
  std::vector<Int_t> *rank = nullptr, *trk1 = nullptr, *trk2 = nullptr;
  t->SetBranchAddress("beam_flag", &beam_flag);
  t->SetBranchAddress("blc1_x", &x1); t->SetBranchAddress("blc1_u", &u1);
  t->SetBranchAddress("blc1_y", &y1); t->SetBranchAddress("blc1_v", &v1);
  t->SetBranchAddress("blc2_x", &x2); t->SetBranchAddress("blc2_u", &u2);
  t->SetBranchAddress("blc2_y", &y2); t->SetBranchAddress("blc2_v", &v2);
  t->SetBranchAddress("d5_delta", &d5_delta);
  t->SetBranchAddress("d5_pair_rank", &rank);
  t->SetBranchAddress("d5_blc1_trk", &trk1);
  t->SetBranchAddress("d5_blc2_trk", &trk2);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  auto* h_dd = new TH1D("h_dd", "#delta_{x} - #delta_{u}; #delta_{x} - #delta_{u} [%]; count", 400, -3, 3);
  auto* h_dy = new TH1D("h_dy", "#Deltay = y_{BLC2} - y_{mtx}(BLC1) @ D5ZOut; #Deltay [mm]; count", 400, -20, 20);
  auto* h_dv = new TH1D("h_dv", "#Deltav = v_{BLC2} - v_{mtx}(BLC1); #Deltav [mrad]; count", 400, -10, 10);
  auto* p_dd_delta = new TProfile("p_dd_delta", "#delta_{x}-#delta_{u} vs #delta; d5_delta [%]; #delta_{x}-#delta_{u} [%]", 40, -5, 10);
  auto* p_dd_u1 = new TProfile("p_dd_u1", "#delta_{x}-#delta_{u} vs BLC1 u; u_{BLC1} [mrad]; #delta_{x}-#delta_{u} [%]", 40, -20, 20);
  auto* p_dd_x1 = new TProfile("p_dd_x1", "#delta_{x}-#delta_{u} vs BLC1 x; x_{BLC1} [mm]; #delta_{x}-#delta_{u} [%]", 40, -60, 60);
  auto* p_dy_delta = new TProfile("p_dy_delta", "#Deltay vs #delta; d5_delta [%]; #Deltay [mm]", 40, -5, 10);
  auto* p_dv_delta = new TProfile("p_dv_delta", "#Deltav vs #delta; d5_delta [%]; #Deltav [mrad]", 40, -5, 10);
  auto* p_dy_y1 = new TProfile("p_dy_y1", "#Deltay vs BLC1 y; y_{BLC1} [mm]; #Deltay [mm]", 40, -30, 30);

  Long64_t n_used = 0;
  for (Long64_t e = 0, ne = t->GetEntries(); e < ne; ++e) {
    t->GetEntry(e);
    if (flag >= 0 && beam_flag != flag) continue;
    for (std::size_t k = 0; k < rank->size(); ++k) {
      if ((*rank)[k] != 0) continue;
      const Int_t a = (*trk1)[k], b = (*trk2)[k];
      const Double_t x2z = (*x2)[b] + (*u2)[b] * mtx.z_out;
      const Double_t y2z = (*y2)[b] + (*v2)[b] * mtx.z_out;
      const Double_t u2m = std::atan((*u2)[b]) * 1000.0;
      const Double_t v2m = std::atan((*v2)[b]) * 1000.0;
      const Double_t dx = solve_delta(mtx, (*x1)[a], (*u1)[a], (*y1)[a], (*v1)[a], 0, x2z);
      const Double_t du = solve_delta(mtx, (*x1)[a], (*u1)[a], (*y1)[a], (*v1)[a], 1, u2m);
      if (!std::isfinite(dx) || !std::isfinite(du)) continue;
      Double_t out[4];
      transport(mtx, (*x1)[a], (*u1)[a], (*y1)[a], (*v1)[a], (*d5_delta)[k], out);
      const Double_t dy = y2z - out[2];
      const Double_t dv = v2m - out[3];
      const Double_t delta = (*d5_delta)[k];
      h_dd->Fill(dx - du);
      h_dy->Fill(dy);
      h_dv->Fill(dv);
      p_dd_delta->Fill(delta, dx - du);
      p_dd_u1->Fill((*u1)[a] * 1000.0, dx - du);
      p_dd_x1->Fill((*x1)[a] + (*u1)[a] * kZExtrapBLC1, dx - du);
      if (std::fabs(dy) < 20.0) {
        p_dy_delta->Fill(delta, dy);
        p_dy_y1->Fill((*y1)[a], dy);
      }
      if (std::fabs(dv) < 10.0) p_dv_delta->Fill(delta, dv);
      ++n_used;
    }
  }

  const GaussResult g_dd = fit_core(h_dd);
  const GaussResult g_dy = fit_core(h_dy);
  const GaussResult g_dv = fit_core(h_dv);
  Double_t dy_slope = 0.0, dv_slope = 0.0;
  const Double_t dy_at0 = fit_line(p_dy_delta, kDeltaFitLo, kDeltaFitHi, &dy_slope);
  const Double_t dv_at0 = fit_line(p_dv_delta, kDeltaFitLo, kDeltaFitHi, &dv_slope);
  const Double_t s_rec = (k_s != 0.0) ? -g_dd.mean / k_s : 0.0;

  // Text summary (key value)
  std::ofstream ofs(path_txt.Data());
  ofs << "run " << run << "\n"
      << "beam_flag " << flag_label(flag) << "\n"
      << "d5_root " << path_d5 << "\n"
      << "matrix " << path_mtx << "\n"
      << "dcgeo " << path_geo << "\n"
      << "n_pairs " << n_used << "\n"
      << Form("blc1a_z %.3f\nblc1b_z %.3f\n", z_a, z_b)
      << Form("k_s_percent_per_mm %.5f\n", k_s)
      << Form("ddxu_mean_percent %.5f\nddxu_mean_err %.5f\nddxu_sigma %.5f\n", g_dd.mean, g_dd.mean_err, g_dd.sigma)
      << Form("blc1_ab_opposite_shift_s_mm %.5f\n", s_rec)
      << Form("dy_mean_mm %.4f\ndy_sigma_mm %.4f\ndy_at_delta0_mm %.4f\ndy_slope_mm_per_percent %.4f\n", g_dy.mean, g_dy.sigma, dy_at0, dy_slope)
      << Form("dv_mean_mrad %.4f\ndv_sigma_mrad %.4f\ndv_at_delta0_mrad %.4f\ndv_slope_mrad_per_percent %.4f\n", g_dv.mean, g_dv.sigma, dv_at0, dv_slope)
      << Form("apply_x --residual-offset BLC1a:%.4f,0 --residual-offset BLC1b:%.4f,0\n", s_rec, -s_rec);
  for (Int_t i = 1; i <= p_dd_x1->GetNbinsX(); ++i)
    if (p_dd_x1->GetBinEntries(i) > 0)
      ofs << Form("ddxu_vs_x1 %.2f %.5f %.5f %.0f\n", p_dd_x1->GetBinCenter(i), p_dd_x1->GetBinContent(i),
                  p_dd_x1->GetBinError(i), p_dd_x1->GetBinEntries(i));
  for (Int_t i = 1; i <= p_dd_delta->GetNbinsX(); ++i)
    if (p_dd_delta->GetBinEntries(i) > 0)
      ofs << Form("ddxu_vs_delta %.2f %.5f %.5f %.0f\n", p_dd_delta->GetBinCenter(i), p_dd_delta->GetBinContent(i),
                  p_dd_delta->GetBinError(i), p_dd_delta->GetBinEntries(i));
  ofs.close();

  // PDF
  auto* c = new TCanvas("c", "c", 1200, 900);
  c->Print(path_pdf + "[");
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  Double_t yy = 0.92;
  auto line = [&](const TString& s) { tex->DrawLatex(0.06, yy, s); yy -= 0.045; };
  line(Form("D5 fit-free review: run %05d, %s, %lld pairs (rank0)", run, flag_label(flag).Data(), n_used));
  line(Form("D5 ROOT: %s", path_d5.Data()));
  line(Form("matrix: %s", path_mtx.Data()));
  line(Form("DCGEO: %s  (BLC1a z=%.1f, BLC1b z=%.1f mm)", path_geo.Data(), z_a, z_b));
  line("Horizontal: #delta_{x} (BLC2 x) - #delta_{u} (BLC2 u), BLC1 track transported alone");
  line(Form("  <#delta_{x}-#delta_{u}> = %+.4f #pm %.4f %%   (#sigma = %.3f %%)", g_dd.mean, g_dd.mean_err, g_dd.sigma));
  line(Form("  k = d<#delta_{x}-#delta_{u}>/ds = %.4f %%/mm  (BLC1a +s, BLC1b -s)", k_s));
  line(Form("  s = %+.4f mm  #rightarrow  --residual-offset BLC1a:%.4f,0 --residual-offset BLC1b:%.4f,0", s_rec, s_rec, -s_rec));
  line("  (BLC1 common x shift is degenerate with #delta; not determined here)");
  line("Vertical (info): #Deltay, #Deltav at D5ZOut, BLC1 transported with fitted #delta");
  line(Form("  <#Deltay> = %+.3f mm (#sigma %.2f), at #delta=0: %+.3f mm, slope %+.3f mm/%%", g_dy.mean, g_dy.sigma, dy_at0, dy_slope));
  line(Form("  <#Deltav> = %+.4f mrad (#sigma %.2f), at #delta=0: %+.4f mrad, slope %+.4f mrad/%%", g_dv.mean, g_dv.sigma, dv_at0, dv_slope));
  c->Print(path_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_dd);
  c->cd(2); p_dd_delta->Draw();
  c->cd(3); p_dd_u1->Draw();
  c->cd(4); p_dd_x1->Draw();
  c->Print(path_pdf);

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_dy);
  c->cd(2); draw_with_fit(h_dv);
  c->cd(3); draw_with_fit(p_dy_delta);
  c->cd(4); draw_with_fit(p_dv_delta);
  c->Print(path_pdf);

  c->Clear();
  p_dy_y1->Draw();
  c->Print(path_pdf);
  c->Print(path_pdf + "]");

  std::cout << "wrote " << path_pdf << "\n      " << path_txt << std::endl;
  std::cout << Form("ddxu=%+.4f%% k=%.4f%%/mm s=%+.4f mm | dy=%+.3f (d0 %+.3f) dv=%+.4f mrad",
                    g_dd.mean, k_s, s_rec, g_dy.mean, dy_at0, g_dv.mean) << std::endl;
  return 0;
}
