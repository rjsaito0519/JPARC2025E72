// BLC1a/b position against BLC2 through the D5 matrix (measurement for run_d5.py pos)
// Usage: d5_pos_review <run> [--pion | --kaon | --all] [--d5 <D5.root>] [--d5ref <D5_bc1reso.root>]
//                      [--mtx <D5TransferMatrix>] [--dcgeo <DCGEO>] [-o <out.pdf>] [--txt <out.txt>]
// Defaults (S = Pi for --pion, K for --kaon / --all):
//           D5     DECODE_DIR/runXXXXX/runXXXXX_D5_S.root            (normal D5 fit)
//           D5ref  SCRATCH_DIR/runXXXXX/runXXXXX_D5_S_bc1reso.root   (BLC1 Res inflated)
//           mtx    ANALYZER_DIR/param/TM/e72/D5TransferMatrix_runXXXXX
//           geo    ANALYZER_DIR/param/DCGEO/e72/DCGeomParam_runXXXXX_S
//           pdf    OUTPUT_DIR/img/runXXXXX/d5_pos_review_runXXXXX.pdf
//           txt    kTxtDir/runXXXXX/d5_pos_review_runXXXXX.txt  (numbers for run_d5.py; not in img)
//           beam flag: --pion (beam_flag == 1)
//
// dy (BLC2-referenced fit, D5ref): with the BLC1 Res inflated the D5 fit line at BLC1 (d5_fit_*0)
// is set by BLC2 alone. r(z) = fit line - BLC1 track at the mean z of BLC1a and of BLC1b; the
// core mean of r is the --residual-offset of that chamber. x is shown for information only.
//
// s (fit-free, D5): for each rank0 pair, solve delta so that the matrix reproduces the BLC2 x at
// D5ZOut (delta_x) and the BLC2 slope u (delta_u) from the BLC1 track alone;
// <delta_x - delta_u> = 0 for consistent geometry and matrix. An opposite x shift of BLC1a (+s)
// and BLC1b (-s) changes it by k*s with k = (-R12/R16 + R22/R26) * (-2000 / (z_b - z_a)) [%/mm],
// so s = -<delta_x - delta_u> / k. The common x shift of BLC1 is degenerate with delta.
// Also dy, dv at D5ZOut (BLC2 - BLC1 transported with the fitted delta), information only.

#include "ana_helper.h"
#include "d5_helper.h"
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
constexpr Double_t kDeltaFitLo = -2.0; // [%] range for linear fits vs delta
constexpr Double_t kDeltaFitHi = 6.0;

struct D5Matrix
{
  Double_t r[6][6] = {};
  Double_t t[5][6][6] = {};
  Double_t z_in = 0.0;
  Double_t z_out = 0.0;
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

// branches of the d5 tree used here
struct D5Tree
{
  TFile* file = nullptr;
  TTree* tree = nullptr;
  Int_t beam_flag = 0;
  std::vector<Double_t> *x1 = nullptr, *u1 = nullptr, *y1 = nullptr, *v1 = nullptr;
  std::vector<Double_t> *x2 = nullptr, *u2 = nullptr, *y2 = nullptr, *v2 = nullptr;
  std::vector<Double_t> *fx = nullptr, *fu = nullptr, *fy = nullptr, *fv = nullptr;
  std::vector<Double_t>* delta = nullptr;
  std::vector<Int_t> *rank = nullptr, *trk1 = nullptr, *trk2 = nullptr;

  Bool_t open(const TString& path)
  {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return kFALSE;
    tree = dynamic_cast<TTree*>(file->Get("d5"));
    if (!tree) return kFALSE;
    tree->SetBranchAddress("beam_flag", &beam_flag);
    tree->SetBranchAddress("blc1_x", &x1); tree->SetBranchAddress("blc1_u", &u1);
    tree->SetBranchAddress("blc1_y", &y1); tree->SetBranchAddress("blc1_v", &v1);
    tree->SetBranchAddress("blc2_x", &x2); tree->SetBranchAddress("blc2_u", &u2);
    tree->SetBranchAddress("blc2_y", &y2); tree->SetBranchAddress("blc2_v", &v2);
    tree->SetBranchAddress("d5_fit_x0", &fx); tree->SetBranchAddress("d5_fit_u0", &fu);
    tree->SetBranchAddress("d5_fit_y0", &fy); tree->SetBranchAddress("d5_fit_v0", &fv);
    tree->SetBranchAddress("d5_delta", &delta);
    tree->SetBranchAddress("d5_pair_rank", &rank);
    tree->SetBranchAddress("d5_blc1_trk", &trk1);
    tree->SetBranchAddress("d5_blc2_trk", &trk2);
    return kTRUE;
  }
};

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0
            << " <run> [--pion | --kaon | --all] [--d5 <file>] [--d5ref <file>] [--mtx <file>]"
            << " [--dcgeo <file>] [-o <out.pdf>] [--txt <out.txt>]\n";
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  if (argc < 2) { usage(argv[0]); return 1; }
  const Int_t run = std::atoi(argv[1]);
  d5_helper::BeamSel beam;
  TString path_d5, path_ref, path_geo, path_txt;
  TString path_mtx = Form("%s/param/TM/e72/D5TransferMatrix_run%05d", ANALYZER_DIR.Data(), run);
  TString path_pdf = Form("%s/d5_pos_review_run%05d.pdf",
                          ana_helper::get_img_dir(OUTPUT_DIR, run).Data(), run);
  for (Int_t i = 2; i < argc; ++i) {
    const TString a = argv[i];
    if (d5_helper::parse_beam_option(a, beam)) continue;
    if (a == "--d5" && i + 1 < argc) path_d5 = argv[++i];
    else if (a == "--d5ref" && i + 1 < argc) path_ref = argv[++i];
    else if (a == "--mtx" && i + 1 < argc) path_mtx = argv[++i];
    else if (a == "--dcgeo" && i + 1 < argc) path_geo = argv[++i];
    else if (a == "-o" && i + 1 < argc) path_pdf = argv[++i];
    else if (a == "--txt" && i + 1 < argc) path_txt = argv[++i];
    else { usage(argv[0]); return 1; }
  }
  const TString& sfx = beam.param_suffix;
  if (path_d5.IsNull())
    path_d5 = Form("%s/run%05d/run%05d_D5_%s.root", DECODE_DIR.Data(), run, run, sfx.Data());
  if (path_ref.IsNull())
    path_ref = Form("%s/run%05d/run%05d_D5_%s_bc1reso.root", SCRATCH_DIR.Data(), run, run, sfx.Data());
  if (path_geo.IsNull())
    path_geo = Form("%s/param/DCGEO/e72/DCGeomParam_run%05d_%s", ANALYZER_DIR.Data(), run, sfx.Data());
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
  const auto geo = d5_helper::load_dcgeo_blc(path_geo);
  Double_t z_a = 0.0, z_b = 0.0;
  if (!d5_helper::chamber_mean_z(geo, "BLC1a", z_a) || !d5_helper::chamber_mean_z(geo, "BLC1b", z_b)) {
    std::cerr << "cannot read BLC1 z from: " << path_geo << std::endl;
    return 1;
  }
  const Double_t ddu_ddelta_slope =
    -mtx.r[0][1] / mtx.r[0][5] + mtx.r[1][1] / mtx.r[1][5]; // [%/mrad]
  const Double_t k_s = ddu_ddelta_slope * (-2000.0 / (z_b - z_a)); // [%/mm]

  D5Tree ref, nor;
  if (!ref.open(path_ref)) { std::cerr << "cannot read d5 tree of " << path_ref << std::endl; return 1; }
  if (!nor.open(path_d5)) { std::cerr << "cannot read d5 tree of " << path_d5 << std::endl; return 1; }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);

  // --- dy from the BLC2-referenced fit ---
  auto* h_ya = new TH1D("h_ya", "BLC1a: r_{y} = y_{fit} - y_{BLC1} (BLC2-referenced); r_{y} [mm]; count", 600, -30, 30);
  auto* h_yb = new TH1D("h_yb", "BLC1b: r_{y} = y_{fit} - y_{BLC1} (BLC2-referenced); r_{y} [mm]; count", 600, -30, 30);
  auto* h_xa = new TH1D("h_xa", "BLC1a: r_{x} (info); r_{x} [mm]; count", 400, -4, 4);
  auto* h_xb = new TH1D("h_xb", "BLC1b: r_{x} (info); r_{x} [mm]; count", 400, -4, 4);
  auto* p_ya_delta = new TProfile("p_ya_delta", "BLC1a r_{y} vs #delta (|r_{y}|<10); d5_delta [%]; r_{y} [mm]", 30, -5, 10);
  auto* p_yb_delta = new TProfile("p_yb_delta", "BLC1b r_{y} vs #delta (|r_{y}|<10); d5_delta [%]; r_{y} [mm]", 30, -5, 10);
  Long64_t n_ref = 0;
  for (Long64_t e = 0, ne = ref.tree->GetEntries(); e < ne; ++e) {
    ref.tree->GetEntry(e);
    if (beam.flag >= 0 && ref.beam_flag != beam.flag) continue;
    for (std::size_t k = 0; k < ref.rank->size(); ++k) {
      if ((*ref.rank)[k] != 0) continue;
      const Int_t a = (*ref.trk1)[k];
      const Double_t rya = ((*ref.fy)[k] + (*ref.fv)[k] * z_a) - ((*ref.y1)[a] + (*ref.v1)[a] * z_a);
      const Double_t ryb = ((*ref.fy)[k] + (*ref.fv)[k] * z_b) - ((*ref.y1)[a] + (*ref.v1)[a] * z_b);
      h_ya->Fill(rya);
      h_yb->Fill(ryb);
      h_xa->Fill(((*ref.fx)[k] + (*ref.fu)[k] * z_a) - ((*ref.x1)[a] + (*ref.u1)[a] * z_a));
      h_xb->Fill(((*ref.fx)[k] + (*ref.fu)[k] * z_b) - ((*ref.x1)[a] + (*ref.u1)[a] * z_b));
      if (std::fabs(rya) < 10.0) p_ya_delta->Fill((*ref.delta)[k], rya);
      if (std::fabs(ryb) < 10.0) p_yb_delta->Fill((*ref.delta)[k], ryb);
      ++n_ref;
    }
  }
  auto peak = [](TH1* h) { return h->GetBinCenter(h->GetMaximumBin()); };
  const auto g_ya = d5_helper::fit_core(h_ya, peak(h_ya), h_ya->GetRMS(), 6);
  const auto g_yb = d5_helper::fit_core(h_yb, peak(h_yb), h_yb->GetRMS(), 6);
  const auto g_xa = d5_helper::fit_core(h_xa, peak(h_xa), 0.1, 6);
  const auto g_xb = d5_helper::fit_core(h_xb, peak(h_xb), 0.1, 6);

  // --- s (and dy, dv at D5ZOut) from the BLC1 track transported alone ---
  auto* h_dd = new TH1D("h_dd", "#delta_{x} - #delta_{u}; #delta_{x} - #delta_{u} [%]; count", 400, -3, 3);
  auto* h_dy = new TH1D("h_dy", "#Deltay = y_{BLC2} - y_{mtx}(BLC1) @ D5ZOut; #Deltay [mm]; count", 400, -20, 20);
  auto* h_dv = new TH1D("h_dv", "#Deltav = v_{BLC2} - v_{mtx}(BLC1); #Deltav [mrad]; count", 400, -10, 10);
  auto* p_dd_delta = new TProfile("p_dd_delta", "#delta_{x}-#delta_{u} vs #delta; d5_delta [%]; #delta_{x}-#delta_{u} [%]", 40, -5, 10);
  auto* p_dd_u1 = new TProfile("p_dd_u1", "#delta_{x}-#delta_{u} vs BLC1 u; u_{BLC1} [mrad]; #delta_{x}-#delta_{u} [%]", 40, -20, 20);
  auto* p_dd_x1 = new TProfile("p_dd_x1", "#delta_{x}-#delta_{u} vs BLC1 x; x_{BLC1} [mm]; #delta_{x}-#delta_{u} [%]", 40, -60, 60);
  auto* p_dy_delta = new TProfile("p_dy_delta", "#Deltay vs #delta; d5_delta [%]; #Deltay [mm]", 40, -5, 10);
  auto* p_dv_delta = new TProfile("p_dv_delta", "#Deltav vs #delta; d5_delta [%]; #Deltav [mrad]", 40, -5, 10);
  Long64_t n_nor = 0;
  for (Long64_t e = 0, ne = nor.tree->GetEntries(); e < ne; ++e) {
    nor.tree->GetEntry(e);
    if (beam.flag >= 0 && nor.beam_flag != beam.flag) continue;
    for (std::size_t k = 0; k < nor.rank->size(); ++k) {
      if ((*nor.rank)[k] != 0) continue;
      const Int_t a = (*nor.trk1)[k], b = (*nor.trk2)[k];
      const Double_t x1 = (*nor.x1)[a], u1 = (*nor.u1)[a], y1 = (*nor.y1)[a], v1 = (*nor.v1)[a];
      const Double_t x2z = (*nor.x2)[b] + (*nor.u2)[b] * mtx.z_out;
      const Double_t y2z = (*nor.y2)[b] + (*nor.v2)[b] * mtx.z_out;
      const Double_t u2m = std::atan((*nor.u2)[b]) * 1000.0;
      const Double_t v2m = std::atan((*nor.v2)[b]) * 1000.0;
      const Double_t dx = solve_delta(mtx, x1, u1, y1, v1, 0, x2z);
      const Double_t du = solve_delta(mtx, x1, u1, y1, v1, 1, u2m);
      if (!std::isfinite(dx) || !std::isfinite(du)) continue;
      const Double_t delta = (*nor.delta)[k];
      Double_t out[4];
      transport(mtx, x1, u1, y1, v1, delta, out);
      const Double_t dy = y2z - out[2];
      const Double_t dv = v2m - out[3];
      h_dd->Fill(dx - du);
      h_dy->Fill(dy);
      h_dv->Fill(dv);
      p_dd_delta->Fill(delta, dx - du);
      p_dd_u1->Fill(u1 * 1000.0, dx - du);
      p_dd_x1->Fill(x1, dx - du);
      if (std::fabs(dy) < 20.0) p_dy_delta->Fill(delta, dy);
      if (std::fabs(dv) < 10.0) p_dv_delta->Fill(delta, dv);
      ++n_nor;
    }
  }
  const auto g_dd = d5_helper::fit_core(h_dd, h_dd->GetMean(), h_dd->GetRMS(), 1, 20);
  const auto g_dy = d5_helper::fit_core(h_dy, h_dy->GetMean(), h_dy->GetRMS(), 1, 20);
  const auto g_dv = d5_helper::fit_core(h_dv, h_dv->GetMean(), h_dv->GetRMS(), 1, 20);
  Double_t dy_slope = 0.0, dv_slope = 0.0;
  const Double_t dy_at0 = fit_line(p_dy_delta, kDeltaFitLo, kDeltaFitHi, &dy_slope);
  const Double_t dv_at0 = fit_line(p_dv_delta, kDeltaFitLo, kDeltaFitHi, &dv_slope);
  const Double_t s_rec = (k_s != 0.0) ? -g_dd.mean / k_s : 0.0;

  // --- text (key value) ---
  std::ofstream ofs(path_txt.Data());
  ofs << "run " << run << "\n"
      << "beam_flag " << beam.label << "\n"
      << "d5_root " << path_d5 << "\n"
      << "d5ref_root " << path_ref << "\n"
      << "matrix " << path_mtx << "\n"
      << "dcgeo " << path_geo << "\n"
      << "n_pairs_ref " << n_ref << "\n"
      << "n_pairs " << n_nor << "\n"
      << Form("blc1a_z %.3f\nblc1b_z %.3f\n", z_a, z_b)
      << Form("dy_a_mm %.5f\ndy_a_err_mm %.5f\ndy_a_sigma_mm %.4f\n", g_ya.mean, g_ya.mean_err, g_ya.sigma)
      << Form("dy_b_mm %.5f\ndy_b_err_mm %.5f\ndy_b_sigma_mm %.4f\n", g_yb.mean, g_yb.mean_err, g_yb.sigma)
      << Form("dx_a_mm_info %.5f\ndx_b_mm_info %.5f\n", g_xa.mean, g_xb.mean)
      << Form("k_s_percent_per_mm %.5f\n", k_s)
      << Form("ddxu_mean_percent %.5f\nddxu_mean_err %.5f\nddxu_sigma %.5f\n", g_dd.mean, g_dd.mean_err, g_dd.sigma)
      << Form("blc1_ab_opposite_shift_s_mm %.5f\n", s_rec)
      << Form("dy_mean_mm %.4f\ndy_sigma_mm %.4f\ndy_at_delta0_mm %.4f\ndy_slope_mm_per_percent %.4f\n", g_dy.mean, g_dy.sigma, dy_at0, dy_slope)
      << Form("dv_mean_mrad %.4f\ndv_sigma_mrad %.4f\ndv_at_delta0_mrad %.4f\ndv_slope_mrad_per_percent %.4f\n", g_dv.mean, g_dv.sigma, dv_at0, dv_slope)
      << Form("apply --residual-offset BLC1a:%.4f,%.4f --residual-offset BLC1b:%.4f,%.4f\n",
              s_rec, g_ya.mean, -s_rec, g_yb.mean);
  for (Int_t i = 1; i <= p_dd_x1->GetNbinsX(); ++i)
    if (p_dd_x1->GetBinEntries(i) > 0)
      ofs << Form("ddxu_vs_x1 %.2f %.5f %.5f %.0f\n", p_dd_x1->GetBinCenter(i), p_dd_x1->GetBinContent(i),
                  p_dd_x1->GetBinError(i), p_dd_x1->GetBinEntries(i));
  for (Int_t i = 1; i <= p_dd_delta->GetNbinsX(); ++i)
    if (p_dd_delta->GetBinEntries(i) > 0)
      ofs << Form("ddxu_vs_delta %.2f %.5f %.5f %.0f\n", p_dd_delta->GetBinCenter(i), p_dd_delta->GetBinContent(i),
                  p_dd_delta->GetBinError(i), p_dd_delta->GetBinEntries(i));
  ofs.close();

  // --- PDF ---
  auto* c = new TCanvas("c", "c", 1200, 900);
  c->Print(path_pdf + "[");
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.026);
  Double_t yy = 0.93;
  auto line = [&](const TString& s) { tex->DrawLatex(0.05, yy, s); yy -= 0.042; };
  line(Form("BLC1 position against BLC2: run %05d, %s", run, beam.label.Data()));
  line(Form("D5 (normal): %s  (%lld rank0 pairs)", gSystem->BaseName(path_d5), n_nor));
  line(Form("D5 (BLC2-referenced): %s  (%lld rank0 pairs)", gSystem->BaseName(path_ref), n_ref));
  line(Form("matrix: %s;  DCGEO: %s", gSystem->BaseName(path_mtx), gSystem->BaseName(path_geo)));
  line(Form("BLC1a z = %.1f mm, BLC1b z = %.1f mm", z_a, z_b));
  line("dy: D5 fit line (BLC1 Res inflated: set by BLC2) - BLC1 track, core Gauss mean");
  line(Form("  BLC1a dy = %+.4f #pm %.4f mm  (#sigma %.2f mm)", g_ya.mean, g_ya.mean_err, g_ya.sigma));
  line(Form("  BLC1b dy = %+.4f #pm %.4f mm  (#sigma %.2f mm)", g_yb.mean, g_yb.mean_err, g_yb.sigma));
  line("s: <#delta_{x}-#delta_{u}> of the BLC1 track transported alone (fit-free)");
  line(Form("  <#delta_{x}-#delta_{u}> = %+.4f #pm %.4f %%,  k = %.4f %%/mm  #rightarrow  s = %+.4f mm (BLC1a +s, BLC1b -s)",
            g_dd.mean, g_dd.mean_err, k_s, s_rec));
  line(Form("#rightarrow --residual-offset BLC1a:%.4f,%.4f --residual-offset BLC1b:%.4f,%.4f",
            s_rec, g_ya.mean, -s_rec, g_yb.mean));
  line("(BLC1 common x is degenerate with #delta; not determined)");
  line(Form("info: #Deltay = %+.3f mm (at #delta=0 %+.3f), #Deltav = %+.4f mrad (at #delta=0 %+.4f)",
            g_dy.mean, dy_at0, g_dv.mean, dv_at0));
  c->Print(path_pdf, "Title:summary");

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_ya);
  c->cd(2); draw_with_fit(h_yb);
  c->cd(3); p_ya_delta->Draw();
  c->cd(4); p_yb_delta->Draw();
  c->Print(path_pdf, "Title:dy (BLC2-referenced)");

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_dd);
  c->cd(2); p_dd_delta->Draw();
  c->cd(3); p_dd_u1->Draw();
  c->cd(4); p_dd_x1->Draw();
  c->Print(path_pdf, "Title:s (fit-free)");

  c->Clear();
  c->Divide(2, 2);
  c->cd(1); draw_with_fit(h_dy);
  c->cd(2); draw_with_fit(h_dv);
  c->cd(3); draw_with_fit(p_dy_delta);
  c->cd(4); draw_with_fit(p_dv_delta);
  c->Print(path_pdf, "Title:dy dv at D5ZOut (info)");

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); draw_with_fit(h_xa);
  c->cd(2); draw_with_fit(h_xb);
  c->Print(path_pdf, "Title:x (info)");
  c->Print(path_pdf + "]");

  std::cout << "wrote " << path_pdf << "\n      " << path_txt << std::endl;
  std::cout << Form("dy_a=%+.4f+-%.4f dy_b=%+.4f+-%.4f mm, s=%+.4f mm (%s)", g_ya.mean, g_ya.mean_err,
                    g_yb.mean, g_yb.mean_err, s_rec, beam.label.Data()) << std::endl;
  return 0;
}
