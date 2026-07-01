// BcIn / BcOut / Hodo (BH2 + BHT) correlation review + D5 transport check -> PDF
// Usage: blc_bh2_correlation_review <run_number> [--bcin <bcin.root>] [-o <out.pdf>] [--swap-bcin-xy]
// BcIn/BcOut/Hodo default: run%05d_{BcIn,BcOut,Hodo}.root under DATA_DIR or DECODE_DIR

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

constexpr Int_t kNumBh2Seg = 15;
constexpr Double_t kBh2Width = 14.0;
constexpr Double_t kBh2Z = -745.7;
constexpr Double_t kBh2SegLo = -0.5;
constexpr Double_t kBh2SegHi = 14.5;

constexpr Int_t kNumBhtSeg = 63;
constexpr Double_t kBhtZ = 0.0;
constexpr Double_t kBhtSegLo = -0.5;
constexpr Double_t kBhtSegHi = 62.5;

constexpr Double_t kPosRangeLo = -120.0;
constexpr Double_t kPosRangeHi = 120.0;

constexpr Double_t kD5Delta = 0.0;
constexpr Int_t kD5XYBins = 200;
constexpr Double_t kD5XYLo = -200.0;
constexpr Double_t kD5XYHi = 200.0;

constexpr const char* kD5ParamRel = "param/D5TransferMatrix.param";

struct D5MatrixData
{
  bool ready = false;
  Double_t z_in  = 0.0;
  Double_t z_out = 0.0;
  Double_t matrix_1st[6][6] = {};
  Double_t matrix_2nd[5][6][6] = {};
};

struct Bh2HistSet
{
  TH2D* x0_vs_seg = nullptr;
  TH2D* y0_vs_seg = nullptr;
  TH2D* u0_vs_seg = nullptr;
  TH2D* v0_vs_seg = nullptr;
};

struct BhtHistSet
{
  TH2D* x0_vs_seg = nullptr;
  TH2D* y0_vs_seg = nullptr;
  TH2D* u0_vs_seg = nullptr;
  TH2D* v0_vs_seg = nullptr;
};

struct Bh2D5HistSet
{
  TH2D* xatbh2_d5_vs_xmm = nullptr;
  TH2D* yatbh2_d5_vs_xmm = nullptr;
};

struct D5Blc2HistSet
{
  TH2D* mtx_x_vs_blc2x = nullptr;
  TH2D* mtx_y_vs_blc2y = nullptr;
  TH2D* mtx_x_vs_blc2y = nullptr;
  TH2D* mtx_y_vs_blc2x = nullptr;
};

struct RunStats
{
  Long64_t n_hodo = 0;
  Long64_t n_bh2_hit = 0;
  Long64_t n_bht_hit = 0;
  Long64_t n_synced_bh2 = 0;
  Long64_t n_synced_bht = 0;
  Long64_t n_filled_bh2 = 0;
  Long64_t n_filled_bht = 0;
  Long64_t n_filled_d5_bh2 = 0;
  Long64_t n_filled_d5_blc2 = 0;
};

Double_t
bh2_x_pos(Int_t seg)
{
  const Double_t offset = (kNumBh2Seg - 1) / 2.0;
  return (seg - offset) * kBh2Width;
}

TString
d5_param_path()
{
  return Form("%s/%s", ANALYZER_DIR.Data(), kD5ParamRel);
}

Bool_t
LoadD5Matrix(const TString& filename, D5MatrixData& mtx)
{
  std::ifstream ifs(filename.Data());
  if (!ifs.is_open()) {
    std::cerr << "LoadD5Matrix: file open fail: " << filename << std::endl;
    return kFALSE;
  }

  mtx = D5MatrixData{};

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string key;
    iss >> key;

    if (key == "D5ZIn:") {
      iss >> mtx.z_in;
    } else if (key == "D5ZOut:") {
      iss >> mtx.z_out;
    } else if (key == "D5Matrix:") {
      Int_t order = 0;
      iss >> order;
      if (order == 1) {
        for (Int_t i = 0; i < 6; ++i) {
          for (Int_t j = 0; j < 6; ++j) {
            ifs >> mtx.matrix_1st[i][j];
          }
        }
      } else if (order >= 21 && order <= 25) {
        const Int_t idx = order - 21;
        for (Int_t i = 0; i < 6; ++i) {
          for (Int_t j = 0; j < 6; ++j) {
            ifs >> mtx.matrix_2nd[idx][i][j];
          }
        }
      }
    }
  }

  mtx.ready = true;
  return kTRUE;
}

void
BuildTransportInput(Double_t x0, Double_t u0, Double_t y0, Double_t v0,
                    Double_t delta, Double_t z_ref_in, Double_t* in)
{
  in[0] = x0 + u0 * z_ref_in;
  in[1] = u0 * 1000.0;
  in[2] = y0 + v0 * z_ref_in;
  in[3] = v0 * 1000.0;
  in[4] = delta;
}

Bool_t
TransportFull(const D5MatrixData& mtx, const Double_t* in, Double_t* out)
{
  if (!mtx.ready) return kFALSE;

  const Double_t unit = 10.0;
  Double_t parin[6] = { in[0] / unit, in[1], in[2] / unit, in[3], 0.0, in[4] };

  Double_t out1st[6] = {};
  for (Int_t i = 0; i < 6; ++i) {
    for (Int_t j = 0; j < 6; ++j) {
      out1st[i] += mtx.matrix_1st[i][j] * parin[j];
    }
  }

  Double_t out2nd[6] = {};
  for (Int_t i = 0; i < 5; ++i) {
    for (Int_t j = 0; j < 6; ++j) {
      for (Int_t k = 0; k < 6; ++k) {
        out2nd[i] += mtx.matrix_2nd[i][j][k] * parin[j] * parin[k];
      }
    }
  }

  out[0] = (out1st[0] + out2nd[0]) * unit;
  out[1] = (out1st[1] + out2nd[1]);
  out[2] = (out1st[2] + out2nd[2]) * unit;
  out[3] = (out1st[3] + out2nd[3]);
  return kTRUE;
}

Double_t
ExtrapPosition(Double_t x_at_z, Double_t u_mrad, Double_t z_from, Double_t z_to)
{
  return x_at_z + (u_mrad / 1000.0) * (z_to - z_from);
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " <run_number> [--bcin <bcin.root>] [-o <out.pdf>] [--swap-bcin-xy]\n"
            << "  Default BcIn/BcOut/Hodo: run%05d_{BcIn,BcOut,Hodo}.root under DATA_DIR or DECODE_DIR\n"
            << "  --bcin: override BcIn ROOT (UserBcInTracking bcin tree), e.g. bcin_ta_plus90.root\n"
            << "  --swap-bcin-xy: swap BcIn x0<->y0, u0<->v0 (post-fit diagnostic)\n"
            << "  Default PDF: " << OUTPUT_DIR << "/img/runXXXXX/"
            << "runXXXXX_BLC_BH2_BHT_correlation[_ta_plus90].pdf\n";
}

void
setup_style()
{
  gROOT->SetBatch();
  gStyle->SetOptStat(110);
  gStyle->SetLabelSize(0.042, "XYZ");
  gStyle->SetTitleSize(0.042, "XYZ");
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadTopMargin(0.08);
}

void
prepare_hist_pad(TH2D* h)
{
  if (!h) return;
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.20);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.08);
  h->GetXaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleOffset(1.45);
  h->GetZaxis()->SetTitleOffset(1.25);
}

Int_t
parse_run_number(const char* arg)
{
  if (!arg || !*arg) return -1;
  Int_t run = -1;
  if (std::sscanf(arg, "%d", &run) != 1 || run < 0) return -1;
  return run;
}

TString
bcin_config_label(const TString& bcin_path)
{
  if (bcin_path.Contains("Uta180") || bcin_path.Contains("uta180")) return "BLC1 Uta180";
  if (bcin_path.Contains("ta_plus90")) return "BLC1 ta_plus90";
  if (bcin_path.Contains("ta_minus90")) return "BLC1 ta_minus90";
  if (bcin_path.Contains("_Pi") || bcin_path.Contains("Pi")) return "BLC1 Pi";
  if (bcin_path.Contains("blctest")) return "BLC1 debug";
  return "BLC1 default";
}

TString
hist_title_prefix(const char* detector, const TString& bcin_cfg)
{
  if (strcmp(detector, "BcIn") == 0) {
    return Form("%s [%s]", detector, bcin_cfg.Data());
  }
  return detector;
}

void
fill_bh2_x_edges(Double_t* edges)
{
  for (Int_t s = 0; s < kNumBh2Seg; ++s) {
    edges[s] = bh2_x_pos(s) - kBh2Width / 2.0;
  }
  edges[kNumBh2Seg] = bh2_x_pos(kNumBh2Seg - 1) + kBh2Width / 2.0;
}

TH2D*
make_bh2_d5_vs_xmm_hist(const char* name, const char* title, const char* ytitle)
{
  Double_t xedges[kNumBh2Seg + 1];
  fill_bh2_x_edges(xedges);
  return new TH2D(name, Form("%s;%s;%s", title, "BH2 seg x [mm]", ytitle),
                  kNumBh2Seg, xedges,
                  kD5XYBins, kD5XYLo, kD5XYHi);
}

TString
pdf_suffix_from_bcin(const TString& bcin_path)
{
  if (bcin_path.Contains("Uta180") || bcin_path.Contains("uta180")) return "_Uta180";
  if (bcin_path.Contains("ta_plus90")) return "_ta_plus90";
  if (bcin_path.Contains("ta_minus90")) return "_ta_minus90";
  return "";
}

TString
default_pdf_path(Int_t run_num, const TString& bcin_path)
{
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
  const TString suffix = pdf_suffix_from_bcin(bcin_path);
  return Form("%s/run%05d_BLC_BH2_BHT_correlation%s.pdf",
              img_dir.Data(), run_num, suffix.Data());
}

TString
resolve_input_path(Int_t run_num, const char* tag)
{
  const TString candidates[] = {
    Form("%s/run%05d_%s.root", DATA_DIR.Data(), run_num, tag),
    Form("%s/run%05d/run%05d_%s.root", DECODE_DIR.Data(), run_num, run_num, tag),
  };
  for (const auto& path : candidates) {
    if (gSystem->AccessPathName(path) == 0) return path;
  }
  return candidates[0];
}

std::unordered_map<UInt_t, Long64_t>
build_event_index(TTree* tree)
{
  std::unordered_map<UInt_t, Long64_t> index;
  if (!tree) return index;

  UInt_t event_number = 0;
  tree->SetBranchAddress("event_number", &event_number);
  const Long64_t nent = tree->GetEntries();
  index.reserve(static_cast<std::size_t>(nent));
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tree->GetEntry(ie);
    index[event_number] = ie;
  }
  return index;
}

TH2D*
make_d5_hist(const char* name, const char* title, const char* xtitle, const char* ytitle)
{
  return new TH2D(name, Form("%s;%s;%s", title, xtitle, ytitle),
                  kD5XYBins, kD5XYLo, kD5XYHi,
                  kD5XYBins, kD5XYLo, kD5XYHi);
}

Bh2HistSet
make_bh2_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  Bh2HistSet h;
  h.x0_vs_seg = new TH2D(
    Form("%s_X0_vs_BH2seg", detector),
    Form("%s x0 vs BH2 seg;BH2 seg;x0 [mm]", pfx.Data()),
    kNumBh2Seg, kBh2SegLo, kBh2SegHi,
    120, kPosRangeLo, kPosRangeHi);
  h.y0_vs_seg = new TH2D(
    Form("%s_Y0_vs_BH2seg", detector),
    Form("%s y0 vs BH2 seg;BH2 seg;y0 [mm]", pfx.Data()),
    kNumBh2Seg, kBh2SegLo, kBh2SegHi,
    120, kPosRangeLo, kPosRangeHi);
  h.u0_vs_seg = new TH2D(
    Form("%s_U0_vs_BH2seg", detector),
    Form("%s u0 vs BH2 seg;BH2 seg;u0 [rad]", pfx.Data()),
    kNumBh2Seg, kBh2SegLo, kBh2SegHi,
    120, -0.05, 0.05);
  h.v0_vs_seg = new TH2D(
    Form("%s_V0_vs_BH2seg", detector),
    Form("%s v0 vs BH2 seg;BH2 seg;v0 [rad]", pfx.Data()),
    kNumBh2Seg, kBh2SegLo, kBh2SegHi,
    120, -0.05, 0.05);
  return h;
}

BhtHistSet
make_bht_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  BhtHistSet h;
  h.x0_vs_seg = new TH2D(
    Form("%s_X0_vs_BHTseg", detector),
    Form("%s x0 vs BHT seg;BHT seg;x0 [mm]", pfx.Data()),
    kNumBhtSeg, kBhtSegLo, kBhtSegHi,
    120, kPosRangeLo, kPosRangeHi);
  h.y0_vs_seg = new TH2D(
    Form("%s_Y0_vs_BHTseg", detector),
    Form("%s y0 vs BHT seg;BHT seg;y0 [mm]", pfx.Data()),
    kNumBhtSeg, kBhtSegLo, kBhtSegHi,
    120, kPosRangeLo, kPosRangeHi);
  h.u0_vs_seg = new TH2D(
    Form("%s_U0_vs_BHTseg", detector),
    Form("%s u0 vs BHT seg;BHT seg;u0 [rad]", pfx.Data()),
    kNumBhtSeg, kBhtSegLo, kBhtSegHi,
    120, -0.05, 0.05);
  h.v0_vs_seg = new TH2D(
    Form("%s_V0_vs_BHTseg", detector),
    Form("%s v0 vs BHT seg;BHT seg;v0 [rad]", pfx.Data()),
    kNumBhtSeg, kBhtSegLo, kBhtSegHi,
    120, -0.05, 0.05);
  return h;
}

Bh2D5HistSet
make_bh2_d5_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  Bh2D5HistSet h;
  h.xatbh2_d5_vs_xmm = make_bh2_d5_vs_xmm_hist(
    Form("%s_XatBH2D5_vs_Bh2Xmm", detector),
    Form("%s D5 transport: x0+u0*z @ z_BH2 vs BH2 seg x (15 seg)", pfx.Data()),
    Form("x0+u0*z (D5 @ z=%.1f) [mm]", kBh2Z));
  h.yatbh2_d5_vs_xmm = make_bh2_d5_vs_xmm_hist(
    Form("%s_YatBH2D5_vs_Bh2Xmm", detector),
    Form("%s D5 transport: y0+v0*z @ z_BH2 vs BH2 seg x (15 seg)", pfx.Data()),
    Form("y0+v0*z (D5 @ z=%.1f) [mm]", kBh2Z));
  return h;
}

D5Blc2HistSet
make_d5_blc2_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  D5Blc2HistSet h;
  h.mtx_x_vs_blc2x = make_d5_hist(
    Form("%s_D5_MtxX_vs_BLC2x", detector),
    Form("%s D5 x0+u0*z vs BLC2 x0+u0*z @ z_out", pfx.Data()),
    "x0+u0*z (D5 @ z_out) [mm]", "x0+u0*z (BLC2 @ z_out) [mm]");
  h.mtx_y_vs_blc2y = make_d5_hist(
    Form("%s_D5_MtxY_vs_BLC2y", detector),
    Form("%s D5 y0+v0*z vs BLC2 y0+v0*z @ z_out", pfx.Data()),
    "y0+v0*z (D5 @ z_out) [mm]", "y0+v0*z (BLC2 @ z_out) [mm]");
  h.mtx_x_vs_blc2y = make_d5_hist(
    Form("%s_D5_MtxX_vs_BLC2y", detector),
    Form("%s D5 x0+u0*z vs BLC2 y0+v0*z @ z_out (cross)", pfx.Data()),
    "x0+u0*z (D5 @ z_out) [mm]", "y0+v0*z (BLC2 @ z_out) [mm]");
  h.mtx_y_vs_blc2x = make_d5_hist(
    Form("%s_D5_MtxY_vs_BLC2x", detector),
    Form("%s D5 y0+v0*z vs BLC2 x0+u0*z @ z_out (cross)", pfx.Data()),
    "y0+v0*z (D5 @ z_out) [mm]", "x0+u0*z (BLC2 @ z_out) [mm]");
  return h;
}

void
fill_bh2_hist(Bh2HistSet& h, Double_t x0, Double_t y0, Double_t u0, Double_t v0,
              Int_t seg)
{
  if (seg < 0 || seg >= kNumBh2Seg) return;

  h.x0_vs_seg->Fill(seg, x0);
  h.y0_vs_seg->Fill(seg, y0);
  h.u0_vs_seg->Fill(seg, u0);
  h.v0_vs_seg->Fill(seg, v0);
}

void
fill_bht_hist(BhtHistSet& h, Double_t x0, Double_t y0, Double_t u0, Double_t v0,
              Int_t seg)
{
  if (seg < 0 || seg >= kNumBhtSeg) return;

  h.x0_vs_seg->Fill(seg, x0);
  h.y0_vs_seg->Fill(seg, y0);
  h.u0_vs_seg->Fill(seg, u0);
  h.v0_vs_seg->Fill(seg, v0);
}

void
fill_bh2_tracks(Bh2HistSet& h, Int_t ntrack,
                const std::vector<Double_t>& x0,
                const std::vector<Double_t>& y0,
                const std::vector<Double_t>& u0,
                const std::vector<Double_t>& v0,
                const std::vector<Double_t>& hit_seg,
                RunStats& stats)
{
  if (ntrack < 1) return;

  for (std::size_t it = 0; it < static_cast<std::size_t>(ntrack); ++it) {
    if (it >= x0.size() || it >= y0.size() || it >= u0.size() || it >= v0.size())
      break;
    for (const Double_t seg_d : hit_seg) {
      const Int_t seg = static_cast<Int_t>(seg_d + 0.5);
      fill_bh2_hist(h, x0[it], y0[it], u0[it], v0[it], seg);
      ++stats.n_filled_bh2;
    }
  }
}

void
fill_bht_tracks(BhtHistSet& h, Int_t ntrack,
                const std::vector<Double_t>& x0,
                const std::vector<Double_t>& y0,
                const std::vector<Double_t>& u0,
                const std::vector<Double_t>& v0,
                const std::vector<Double_t>& hit_seg,
                RunStats& stats)
{
  if (ntrack < 1) return;

  for (std::size_t it = 0; it < static_cast<std::size_t>(ntrack); ++it) {
    if (it >= x0.size() || it >= y0.size() || it >= u0.size() || it >= v0.size())
      break;
    for (const Double_t seg_d : hit_seg) {
      const Int_t seg = static_cast<Int_t>(seg_d + 0.5);
      fill_bht_hist(h, x0[it], y0[it], u0[it], v0[it], seg);
      ++stats.n_filled_bht;
    }
  }
}

void
fill_d5_bcin_bh2(Bh2D5HistSet& h, const D5MatrixData& mtx,
                 Double_t x0, Double_t y0, Double_t u0, Double_t v0,
                 const std::vector<Double_t>& hit_seg, RunStats& stats)
{
  Double_t in[5] = {};
  Double_t out[4] = {};
  BuildTransportInput(x0, u0, y0, v0, kD5Delta, mtx.z_in, in);
  if (!TransportFull(mtx, in, out)) return;

  const Double_t x_bh2_d5 = ExtrapPosition(out[0], out[1], mtx.z_out, kBh2Z);
  const Double_t y_bh2_d5 = ExtrapPosition(out[2], out[3], mtx.z_out, kBh2Z);

  for (const Double_t seg_d : hit_seg) {
    const Int_t seg = static_cast<Int_t>(seg_d + 0.5);
    if (seg < 0 || seg >= kNumBh2Seg) continue;
    const Double_t bh2_xmm = bh2_x_pos(seg);
    h.xatbh2_d5_vs_xmm->Fill(bh2_xmm, x_bh2_d5);
    h.yatbh2_d5_vs_xmm->Fill(bh2_xmm, y_bh2_d5);
    ++stats.n_filled_d5_bh2;
  }
}

void
fill_d5_bcin_blc2(D5Blc2HistSet& h, const D5MatrixData& mtx,
                  Double_t x0, Double_t y0, Double_t u0, Double_t v0,
                  Double_t x2, Double_t y2, Double_t u2, Double_t v2,
                  RunStats& stats)
{
  Double_t in[5] = {};
  Double_t out[4] = {};
  BuildTransportInput(x0, u0, y0, v0, kD5Delta, mtx.z_in, in);
  if (!TransportFull(mtx, in, out)) return;

  const Double_t blc2_x = x2 + u2 * mtx.z_out;
  const Double_t blc2_y = y2 + v2 * mtx.z_out;

  h.mtx_x_vs_blc2x->Fill(out[0], blc2_x);
  h.mtx_y_vs_blc2y->Fill(out[2], blc2_y);
  h.mtx_x_vs_blc2y->Fill(out[0], blc2_y);
  h.mtx_y_vs_blc2x->Fill(out[2], blc2_x);
  ++stats.n_filled_d5_blc2;
}

void
swap_xy_track_params(Int_t ntrack,
                     const std::vector<Double_t>& x0_in,
                     const std::vector<Double_t>& y0_in,
                     const std::vector<Double_t>& u0_in,
                     const std::vector<Double_t>& v0_in,
                     std::vector<Double_t>& x0_out,
                     std::vector<Double_t>& y0_out,
                     std::vector<Double_t>& u0_out,
                     std::vector<Double_t>& v0_out)
{
  const std::size_t n = static_cast<std::size_t>(ntrack);
  x0_out.resize(n);
  y0_out.resize(n);
  u0_out.resize(n);
  v0_out.resize(n);
  for (std::size_t it = 0; it < n; ++it) {
    if (it >= x0_in.size() || it >= y0_in.size() || it >= u0_in.size() || it >= v0_in.size())
      break;
    x0_out[it] = y0_in[it];
    y0_out[it] = x0_in[it];
    u0_out[it] = v0_in[it];
    v0_out[it] = u0_in[it];
  }
}

void
draw_placeholder(const char* message)
{
  auto* tex = new TLatex(0.5, 0.5, message);
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.06);
  tex->Draw();
}

void
draw_hist_or_placeholder(TH2D* h, const char* opt = "COLZ")
{
  if (!h || h->GetEntries() <= 0) {
    draw_placeholder(h ? "empty histogram" : "histogram missing");
    return;
  }
  prepare_hist_pad(h);
  h->Draw(opt);
}

void
print_page(TCanvas* c, const TString& pdf_path)
{
  c->Print(pdf_path);
}

void
draw_title_page(TCanvas* c, const TString& pdf_path, Int_t run_num,
                const TString& path_bcin, const TString& path_bcout,
                const TString& path_hodo, const RunStats& stats,
                bool swap_bcin_xy, const D5MatrixData& mtx,
                const TString& bcin_cfg)
{
  c->Clear();
  c->cd();

  TTimeStamp now;
  now.Set();

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.92, "BLC-BH2/BHT correlation review (+ D5)");
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.5, 0.84, Form("run_number = %d  |  BcIn config: %s",
                                 run_num, bcin_cfg.Data()));
  if (swap_bcin_xy) {
    tex->SetTextColor(kRed);
    tex->DrawLatex(0.5, 0.78, "BcIn: x0<->y0, u0<->v0 swapped (diagnostic)");
    tex->SetTextColor(kBlack);
  }
  tex->SetTextAlign(12);
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.06, 0.72, Form("BcIn : %s", path_bcin.Data()));
  tex->DrawLatex(0.06, 0.66, Form("BcOut: %s", path_bcout.Data()));
  tex->DrawLatex(0.06, 0.60, Form("Hodo : %s", path_hodo.Data()));
  tex->SetTextAlign(22);
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.5, 0.50, Form("Hodo entries: %lld", stats.n_hodo));
  tex->DrawLatex(0.5, 0.44, Form("Synced BH2/BHT fills: %lld / %lld",
                                 stats.n_filled_bh2, stats.n_filled_bht));
  tex->DrawLatex(0.5, 0.38, Form("D5@BH2 fills: %lld, D5 vs BLC2: %lld",
                                 stats.n_filled_d5_bh2, stats.n_filled_d5_blc2));
  tex->DrawLatex(0.5, 0.30, Form("D5 z_in=%.1f mm, z_out=%.1f mm", mtx.z_in, mtx.z_out));
  tex->DrawLatex(0.5, 0.24, Form("z_BH2=%.1f mm, z_BHT=%.1f mm", kBh2Z, kBhtZ));
  tex->DrawLatex(0.5, 0.16, Form("created: %s", now.AsString("s")));
  print_page(c, pdf_path);
}

void
draw_pair_page(TCanvas* c, const TString& pdf_path, TH2D* left, TH2D* right)
{
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_hist_or_placeholder(left);
  c->cd(2);
  draw_hist_or_placeholder(right);
  print_page(c, pdf_path);
}

void
draw_bh2_pages(TCanvas* c, const TString& pdf_path, const char* label,
               const Bh2HistSet& h)
{
  draw_pair_page(c, pdf_path, h.x0_vs_seg, h.y0_vs_seg);
  draw_pair_page(c, pdf_path, h.u0_vs_seg, h.v0_vs_seg);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.88, Form("%s BH2 summary", label));
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.5, 0.76,
                 "BH2 is segmented along X. Expect x0 band vs BH2 segment.");
  tex->DrawLatex(0.5, 0.70,
                 "If y0 or v0 correlates with BH2 segment, suspect x/y swap.");
  tex->DrawLatex(0.5, 0.62,
                 "BLC1 x0,y0,u0,v0 only (no extrapolation to z_BH2).");
  tex->DrawLatex(0.5, 0.54, Form("%s_X0_vs_BH2seg entries: %.0f",
                                 label, h.x0_vs_seg ? h.x0_vs_seg->GetEntries() : 0));
  print_page(c, pdf_path);
}

void
draw_bht_pages(TCanvas* c, const TString& pdf_path, const char* label,
               const BhtHistSet& h)
{
  draw_pair_page(c, pdf_path, h.x0_vs_seg, h.y0_vs_seg);
  draw_pair_page(c, pdf_path, h.u0_vs_seg, h.v0_vs_seg);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.88, Form("%s BHT summary", label));
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.5, 0.76, Form("BHT: 63 segments at z=%.1f mm.", kBhtZ));
  tex->DrawLatex(0.5, 0.70, "x0,y0,u0,v0 at BLC fit reference (no z extrap).");
  tex->DrawLatex(0.5, 0.60, Form("%s_X0_vs_BHTseg entries: %.0f",
                                 label, h.x0_vs_seg ? h.x0_vs_seg->GetEntries() : 0));
  print_page(c, pdf_path);
}

void
draw_d5_bh2_pages(TCanvas* c, const TString& pdf_path, const Bh2D5HistSet& h,
                  const TString& bcin_cfg)
{
  draw_pair_page(c, pdf_path, h.xatbh2_d5_vs_xmm, h.yatbh2_d5_vs_xmm);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.88, Form("BcIn D5 transport @ BH2 [%s]", bcin_cfg.Data()));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.76,
                 "D5: BLC1 -> z_in, matrix to z_out, then x0+u0*z / y0+v0*z to z_BH2.");
  tex->DrawLatex(0.5, 0.69, Form("z_BH2=%.1f mm. Expect slope ~ +1 vs BH2 seg x.", kBh2Z));
  tex->DrawLatex(0.5, 0.62, "No BLC1-only extrapolation to z_BH2 (not measurable there).");
  print_page(c, pdf_path);
}

void
draw_d5_blc2_pages(TCanvas* c, const TString& pdf_path, const D5Blc2HistSet& h,
                   const TString& bcin_cfg)
{
  draw_pair_page(c, pdf_path, h.mtx_x_vs_blc2x, h.mtx_y_vs_blc2y);
  draw_pair_page(c, pdf_path, h.mtx_x_vs_blc2y, h.mtx_y_vs_blc2x);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.88, Form("BcIn D5 transport vs BLC2 [%s]", bcin_cfg.Data()));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.76,
                 "D5 x0+u0*z / y0+v0*z @ z_out vs BLC2 x0+u0*z / y0+v0*z @ z_out.");
  tex->DrawLatex(0.5, 0.70, "Cross terms flat => matrix does not mix x/y.");
  print_page(c, pdf_path);
}

void
draw_readme_page(TCanvas* c, const TString& pdf_path)
{
  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.90, "How to read (TA +/-90 check)");
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.82, Form("BH2 (15 seg, along X at z=%.1f mm):", kBh2Z));
  tex->DrawLatex(0.5, 0.77, "  Normal: x0 shows band vs BH2 segment.");
  tex->DrawLatex(0.5, 0.72, "  Normal: y0 vs BH2 segment is mostly flat.");
  tex->DrawLatex(0.5, 0.67, "  Do NOT use BLC1 x0+u0*z at z_BH2 (not measurable).");
  tex->DrawLatex(0.5, 0.60, "BHT (63 seg): x0,u0 vs segment; y0,v0 mostly flat.");
  tex->DrawLatex(0.5, 0.52,
                 "D5@BH2: matrix z_in->z_out, then x0+u0*z / y0+v0*z to z_BH2.");
  tex->DrawLatex(0.5, 0.45, "Compare ta_plus90 vs ta_minus90 PDFs side by side.");
  tex->DrawLatex(0.5, 0.38, "Swap suspect: y0 or v0 correlates with X-segmented hodoscope.");
  tex->DrawLatex(0.5, 0.30, Form("z_BH2=%.1f mm, z_BHT=%.1f mm", kBh2Z, kBhtZ));
  tex->DrawLatex(0.5, 0.18, "Bh2Xmm = (seg - 7) * 14 mm for seg 0..14");
  print_page(c, pdf_path);
}

Int_t
analyze(Int_t run_num, const TString& path_bcin_in, const TString& output_path_in,
        bool swap_bcin_xy)
{
  setup_style();

  TString path_bcin = path_bcin_in;
  if (path_bcin.IsNull()) {
    path_bcin = resolve_input_path(run_num, "BcIn");
  }
  gSystem->ExpandPathName(path_bcin);
  if (path_bcin.IsNull() || gSystem->AccessPathName(path_bcin) != 0) {
    std::cerr << "Error: BcIn ROOT not found: " << path_bcin << std::endl;
    return 1;
  }

  D5MatrixData mtx;
  if (!LoadD5Matrix(d5_param_path(), mtx)) {
    return 1;
  }

  const TString path_bcout = resolve_input_path(run_num, "BcOut");
  const TString path_hodo = resolve_input_path(run_num, "Hodo");

  auto* f_bcin = TFile::Open(path_bcin);
  auto* f_bcout = TFile::Open(path_bcout);
  auto* f_hodo = TFile::Open(path_hodo);

  if (!f_bcin || f_bcin->IsZombie()) {
    std::cerr << "Error: could not open " << path_bcin << std::endl;
    return 1;
  }
  if (!f_bcout || f_bcout->IsZombie()) {
    std::cerr << "Error: could not open " << path_bcout << std::endl;
    f_bcin->Close();
    return 1;
  }
  if (!f_hodo || f_hodo->IsZombie()) {
    std::cerr << "Error: could not open " << path_hodo << std::endl;
    f_bcin->Close();
    f_bcout->Close();
    return 1;
  }

  auto* tree_bcin = dynamic_cast<TTree*>(f_bcin->Get("bcin"));
  auto* tree_bcout = dynamic_cast<TTree*>(f_bcout->Get("bcout"));
  auto* tree_hodo = dynamic_cast<TTree*>(f_hodo->Get("hodo"));
  if (!tree_bcin || !tree_bcout || !tree_hodo) {
    std::cerr << "Error: required trees not found (bcin, bcout, hodo)." << std::endl;
    f_bcin->Close();
    f_bcout->Close();
    f_hodo->Close();
    return 1;
  }

  std::cout << "BcIn : " << path_bcin << std::endl;
  std::cout << "BcOut: " << path_bcout << std::endl;
  std::cout << "Hodo : " << path_hodo << std::endl;
  std::cout << "D5   : " << d5_param_path()
            << " z_in=" << mtx.z_in << " z_out=" << mtx.z_out << std::endl;

  const TString bcin_cfg = bcin_config_label(path_bcin);

  const auto bcin_idx = build_event_index(tree_bcin);
  const auto bcout_idx = build_event_index(tree_bcout);

  Bh2HistSet h_bcin_bh2 = make_bh2_hist_set("BcIn", bcin_cfg);
  Bh2HistSet h_bcout_bh2 = make_bh2_hist_set("BcOut", bcin_cfg);
  BhtHistSet h_bcin_bht = make_bht_hist_set("BcIn", bcin_cfg);
  BhtHistSet h_bcout_bht = make_bht_hist_set("BcOut", bcin_cfg);
  Bh2D5HistSet h_bcin_d5_bh2 = make_bh2_d5_hist_set("BcIn", bcin_cfg);
  D5Blc2HistSet h_bcin_d5_blc2 = make_d5_blc2_hist_set("BcIn", bcin_cfg);
  RunStats stats;
  stats.n_hodo = tree_hodo->GetEntries();

  UInt_t event_number = 0;
  std::vector<Double_t>* bh2_hit_seg = nullptr;
  std::vector<Double_t>* bht_hit_seg = nullptr;
  tree_hodo->SetBranchAddress("event_number", &event_number);
  tree_hodo->SetBranchAddress("bh2_hit_seg", &bh2_hit_seg);
  tree_hodo->SetBranchAddress("bht_hit_seg", &bht_hit_seg);

  Int_t ntrack_bcin = 0;
  Int_t ntrack_bcout = 0;
  std::vector<Double_t>* x0_bcin = nullptr;
  std::vector<Double_t>* y0_bcin = nullptr;
  std::vector<Double_t>* u0_bcin = nullptr;
  std::vector<Double_t>* v0_bcin = nullptr;
  std::vector<Double_t>* x0_bcout = nullptr;
  std::vector<Double_t>* y0_bcout = nullptr;
  std::vector<Double_t>* u0_bcout = nullptr;
  std::vector<Double_t>* v0_bcout = nullptr;
  tree_bcin->SetBranchAddress("ntrack", &ntrack_bcin);
  tree_bcin->SetBranchAddress("x0", &x0_bcin);
  tree_bcin->SetBranchAddress("y0", &y0_bcin);
  tree_bcin->SetBranchAddress("u0", &u0_bcin);
  tree_bcin->SetBranchAddress("v0", &v0_bcin);
  tree_bcout->SetBranchAddress("ntrack", &ntrack_bcout);
  tree_bcout->SetBranchAddress("x0", &x0_bcout);
  tree_bcout->SetBranchAddress("y0", &y0_bcout);
  tree_bcout->SetBranchAddress("u0", &u0_bcout);
  tree_bcout->SetBranchAddress("v0", &v0_bcout);

  std::vector<Double_t> x0_bcin_sw, y0_bcin_sw, u0_bcin_sw, v0_bcin_sw;

  for (Long64_t ie = 0; ie < stats.n_hodo; ++ie) {
    tree_hodo->GetEntry(ie);

    const bool has_bh2 = bh2_hit_seg && !bh2_hit_seg->empty();
    const bool has_bht = bht_hit_seg && !bht_hit_seg->empty();
    if (!has_bh2 && !has_bht) continue;

    const auto it_in = bcin_idx.find(event_number);
    const auto it_out = bcout_idx.find(event_number);
    if (it_in == bcin_idx.end() || it_out == bcout_idx.end()) continue;

    tree_bcin->GetEntry(it_in->second);
    tree_bcout->GetEntry(it_out->second);

    const std::vector<Double_t>* x0_in = x0_bcin;
    const std::vector<Double_t>* y0_in = y0_bcin;
    const std::vector<Double_t>* u0_in = u0_bcin;
    const std::vector<Double_t>* v0_in = v0_bcin;
    if (swap_bcin_xy && x0_bcin && y0_bcin && u0_bcin && v0_bcin) {
      swap_xy_track_params(ntrack_bcin, *x0_bcin, *y0_bcin, *u0_bcin, *v0_bcin,
                           x0_bcin_sw, y0_bcin_sw, u0_bcin_sw, v0_bcin_sw);
      x0_in = &x0_bcin_sw;
      y0_in = &y0_bcin_sw;
      u0_in = &u0_bcin_sw;
      v0_in = &v0_bcin_sw;
    }

    if (ntrack_bcin >= 1 && ntrack_bcout >= 1 &&
        x0_in && y0_in && u0_in && v0_in &&
        x0_bcout && y0_bcout && u0_bcout && v0_bcout) {
      fill_d5_bcin_blc2(h_bcin_d5_blc2, mtx,
                        (*x0_in)[0], (*y0_in)[0], (*u0_in)[0], (*v0_in)[0],
                        (*x0_bcout)[0], (*y0_bcout)[0], (*u0_bcout)[0], (*v0_bcout)[0],
                        stats);
    }

    if (has_bh2) {
      ++stats.n_bh2_hit;
      ++stats.n_synced_bh2;
      if (x0_in && y0_in && u0_in && v0_in) {
        fill_bh2_tracks(h_bcin_bh2, ntrack_bcin, *x0_in, *y0_in, *u0_in, *v0_in,
                        *bh2_hit_seg, stats);
        if (ntrack_bcin >= 1) {
          fill_d5_bcin_bh2(h_bcin_d5_bh2, mtx,
                           (*x0_in)[0], (*y0_in)[0], (*u0_in)[0], (*v0_in)[0],
                           *bh2_hit_seg, stats);
        }
      }
      if (x0_bcout && y0_bcout && u0_bcout && v0_bcout) {
        fill_bh2_tracks(h_bcout_bh2, ntrack_bcout, *x0_bcout, *y0_bcout, *u0_bcout, *v0_bcout,
                        *bh2_hit_seg, stats);
      }
    }

    if (has_bht) {
      ++stats.n_bht_hit;
      ++stats.n_synced_bht;
      if (x0_in && y0_in && u0_in && v0_in) {
        fill_bht_tracks(h_bcin_bht, ntrack_bcin, *x0_in, *y0_in, *u0_in, *v0_in,
                        *bht_hit_seg, stats);
      }
      if (x0_bcout && y0_bcout && u0_bcout && v0_bcout) {
        fill_bht_tracks(h_bcout_bht, ntrack_bcout, *x0_bcout, *y0_bcout, *u0_bcout, *v0_bcout,
                        *bht_hit_seg, stats);
      }
    }
  }

  TString output_path = output_path_in;
  if (output_path.IsNull()) {
    output_path = default_pdf_path(run_num, path_bcin);
    std::cout << "Default output: " << output_path << std::endl;
  }

  auto* canvas = new TCanvas("blc_bh2_correlation_review", output_path, 1800, 900);
  canvas->Print(output_path + "[");

  draw_title_page(canvas, output_path, run_num, path_bcin, path_bcout, path_hodo, stats,
                  swap_bcin_xy, mtx, bcin_cfg);
  draw_bh2_pages(canvas, output_path, Form("BcIn [%s]", bcin_cfg.Data()), h_bcin_bh2);
  draw_bht_pages(canvas, output_path, Form("BcIn [%s]", bcin_cfg.Data()), h_bcin_bht);
  draw_d5_bh2_pages(canvas, output_path, h_bcin_d5_bh2, bcin_cfg);
  draw_d5_blc2_pages(canvas, output_path, h_bcin_d5_blc2, bcin_cfg);
  draw_bh2_pages(canvas, output_path, "BcOut", h_bcout_bh2);
  draw_bht_pages(canvas, output_path, "BcOut", h_bcout_bht);
  draw_readme_page(canvas, output_path);

  canvas->Print(output_path + "]");
  delete canvas;

  f_bcin->Close();
  f_bcout->Close();
  f_hodo->Close();

  std::cout << "Wrote " << output_path << " (run=" << run_num
            << ", fills_bh2=" << stats.n_filled_bh2
            << ", fills_d5_bh2=" << stats.n_filled_d5_bh2
            << ", fills_d5_blc2=" << stats.n_filled_d5_blc2 << ")" << std::endl;
  return 0;
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  Int_t run_num = -1;
  TString output_path;
  TString path_bcin;
  bool has_output = false;
  bool has_bcin = false;
  bool swap_bcin_xy = false;

  for (Int_t i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    if (strcmp(arg, "-o") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      output_path = argv[++i];
      has_output = true;
      continue;
    }
    if (strcmp(arg, "--bcin") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      path_bcin = argv[++i];
      has_bcin = true;
      continue;
    }
    if (strcmp(arg, "--swap-bcin-xy") == 0) {
      swap_bcin_xy = true;
      continue;
    }
    if (arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
    if (run_num >= 0) {
      std::cerr << "Multiple run numbers specified." << std::endl;
      usage(argv[0]);
      return 1;
    }
    run_num = parse_run_number(arg);
    if (run_num < 0) {
      std::cerr << "Invalid run number: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
  }

  if (run_num < 0) {
    usage(argv[0]);
    return 1;
  }

  if (has_output) gSystem->ExpandPathName(output_path);
  return analyze(run_num, path_bcin, has_output ? output_path : TString(), swap_bcin_xy);
}
