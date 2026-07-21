// BcIn / BcOut / Hodo (BH2 + BHT) correlation review + D5 transport check -> PDF
// Usage: blc_bh2_correlation_review <run_number> [--bcin <bcin.root>] [-o <out.pdf>]
//        [--swap-bcin-xy] [--pion | --kaon]
// BcIn/BcOut/Hodo default: run%05d_{BcIn,BcOut,Hodo}.root under DATA_DIR or DECODE_DIR

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cmath>
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

struct BlcYResidualHistSet
{
  TH1D* dy_at_center = nullptr;
  TH1D* dy_at_center_lowv = nullptr;
  TH2D* v0_vs_dy = nullptr;
};

struct BlcYResidualGaussFit
{
  bool ok = false;
  Double_t mean = 0.0;
  Double_t mean_err = 0.0;
  Double_t sigma = 0.0;
  Double_t sigma_err = 0.0;
  Double_t amp = 0.0;
};

struct BlcBeamProfileHistSet
{
  TH2D* blc1_xy = nullptr;
  TH2D* blc2_xy = nullptr;
};

struct BlcSlopeHistSet
{
  TH2D* v0_blc1_vs_v0_blc2 = nullptr;
  TH1D* dy_extrap = nullptr;
  TH2D* v0_vs_dy_extrap = nullptr;
  TH2D* v0_vs_dy_extrap_zoom = nullptr;
  TH2D* dy_vs_blc1_extrap_z = nullptr;      // BLC2 fixed; scan BLC1 extrap z
  TProfile* dy_mean_vs_blc1_extrap_z = nullptr;
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
  Long64_t n_filled_blc_y_residual = 0;
  Long64_t n_filled_blc_y_residual_lowv = 0;
  Long64_t n_filled_blc_slope = 0;
  Long64_t n_filled_blc1_beam_profile = 0;
  Long64_t n_filled_blc2_beam_profile = 0;
  Long64_t n_skipped_beam_flag = 0;
};

constexpr Double_t kBlcChamberCenterZ = 0.0; // track fit ref z=0 (chamber center)
constexpr Double_t kBlc1V0FlatCut = 0.01;  // |BLC1 v0| cut for flat-track residual
constexpr Double_t kBlc1BeamProfileZ = 0.0;
constexpr Double_t kBlc2BeamProfileZ = -1299.815; // BLC2 chamber center (fixed)
constexpr Double_t kBlc1ExtrapZ = 2900.0; // 3 m from BLC1 chamber center
constexpr Double_t kExtrapDyRange = 500.0;
constexpr Int_t kExtrapDyBins = 1000; // 1 mm/bin over ±500 mm
constexpr Double_t kExtrapZoomV0 = 0.02;
constexpr Double_t kExtrapZoomDy = 100.0;
constexpr Int_t kExtrapZoomBins = 200; // 1 mm/bin over ±100 mm for Δy
// BLC1->BLC2: BLC2 y fixed at chamber center; scan BLC1 extrap z around 2900
constexpr Double_t kBlc1ExtrapScanCenter = kBlc1ExtrapZ;
constexpr Double_t kBlc1ExtrapScanHalfRange = 300.0;
constexpr Double_t kBlc1ExtrapScanLo = kBlc1ExtrapScanCenter - kBlc1ExtrapScanHalfRange;
constexpr Double_t kBlc1ExtrapScanHi = kBlc1ExtrapScanCenter + kBlc1ExtrapScanHalfRange;
constexpr Int_t kBlc1ExtrapScanBins = 60; // 10 mm/bin over ±300 mm
constexpr Double_t kGaussDrawNSigma = 5.0;
constexpr Double_t kGaussDrawMinHalfRange = 15.0; // mm

constexpr Int_t kBeamFlagNoFilter = -1;
constexpr Int_t kBeamFlagPion = 1;
constexpr Int_t kBeamFlagKaon = 2;

TString
beam_flag_filter_label(Int_t beam_flag_filter)
{
  if (beam_flag_filter == kBeamFlagPion) return "pi only (beam_flag=1)";
  if (beam_flag_filter == kBeamFlagKaon) return "K only (beam_flag=2)";
  return "all beam_flag";
}

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
            << "        [--pion | --kaon]\n"
            << "  Default BcIn/BcOut/Hodo: run%05d_{BcIn,BcOut,Hodo}.root under DATA_DIR or DECODE_DIR\n"
            << "  --bcin: override BcIn ROOT (UserBcInTracking bcin tree), e.g. bcin_ta_plus90.root\n"
            << "  --swap-bcin-xy: swap BcIn x0<->y0, u0<->v0 (post-fit diagnostic)\n"
            << "  --pion: use beam_flag=1 events only; --kaon: use beam_flag=2 events only\n"
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

BlcYResidualHistSet
make_blc_y_residual_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  BlcYResidualHistSet h;
  h.dy_at_center = new TH1D(
    Form("%s_BLC2y_minus_BLC1y_at_center", detector),
    Form("%s BLC2-BLC1 #Deltay @ chamber center (z=0);"
         "#Deltay = (y0+v0*z)_{BLC2} - y0_{BLC1} [mm];entries",
         pfx.Data()),
    240, kPosRangeLo, kPosRangeHi);
  h.dy_at_center_lowv = new TH1D(
    Form("%s_BLC2y_minus_BLC1y_at_center_lowv", detector),
    Form("%s BLC2-BLC1 #Deltay @ chamber center (|BLC1 v0|#leq%.3f);"
         "#Deltay = (y0+v0*z)_{BLC2} - y0_{BLC1} [mm];entries",
         pfx.Data(), kBlc1V0FlatCut),
    240, kPosRangeLo, kPosRangeHi);
  h.v0_vs_dy = new TH2D(
    Form("%s_BLC1v0_vs_BLC2y_residual", detector),
    Form("%s BLC1 v0 vs BLC2-BLC1 #Deltay @ chamber center;"
         "BLC1 v0;#Deltay = (y0+v0*z)_{BLC2} - y0_{BLC1} [mm]",
         pfx.Data()),
    120, -0.05, 0.05,
    240, kPosRangeLo, kPosRangeHi);
  return h;
}

BlcSlopeHistSet
make_blc_slope_hist_set(const char* detector, const TString& bcin_cfg)
{
  const TString pfx = hist_title_prefix(detector, bcin_cfg);
  BlcSlopeHistSet h;
  h.v0_blc1_vs_v0_blc2 = new TH2D(
    Form("%s_BLC1v0_vs_BLC2v0", detector),
    Form("%s BLC1 v0 vs BLC2 v0;BLC1 v0;BLC2 v0", pfx.Data()),
    120, -0.05, 0.05,
    120, -0.05, 0.05);
  h.dy_extrap = new TH1D(
    Form("%s_BLC2y_minus_BLC1y_extrap", detector),
    Form("%s BLC2@z=%.3f - BLC1@z=%.0f #Deltay;"
         "#Deltay = (y0+v0*z)_{BLC2} - (y0+v0*z)_{BLC1} [mm];entries",
         pfx.Data(), kBlc2BeamProfileZ, kBlc1ExtrapZ),
    kExtrapDyBins, -kExtrapDyRange, kExtrapDyRange);
  h.v0_vs_dy_extrap = new TH2D(
    Form("%s_BLC1v0_vs_extrap_dy", detector),
    Form("%s BLC1 v0 vs extrap #Deltay;"
         "BLC1 v0;#Deltay = (y0+v0*%.3f)_{BLC2} - (y0+v0*%.0f)_{BLC1} [mm]",
         pfx.Data(), kBlc2BeamProfileZ, kBlc1ExtrapZ),
    120, -0.05, 0.05,
    kExtrapDyBins, -kExtrapDyRange, kExtrapDyRange);
  h.v0_vs_dy_extrap_zoom = new TH2D(
    Form("%s_BLC1v0_vs_extrap_dy_zoom", detector),
    Form("%s BLC1 v0 vs extrap #Deltay (close-up);"
         "BLC1 v0;#Deltay = (y0+v0*%.3f)_{BLC2} - (y0+v0*%.0f)_{BLC1} [mm]",
         pfx.Data(), kBlc2BeamProfileZ, kBlc1ExtrapZ),
    kExtrapZoomBins, -kExtrapZoomV0, kExtrapZoomV0,
    kExtrapZoomBins, -kExtrapZoomDy, kExtrapZoomDy);
  h.dy_vs_blc1_extrap_z = new TH2D(
    Form("%s_dy_vs_blc1_extrap_z", detector),
    Form("%s #Deltay vs BLC1 extrap z (BLC2 fixed @ %.3f, z=%.0f #pm %.0f);"
         "z_{BLC1 extrap} [mm];#Deltay = y_{BLC2}(z_{BLC2}) - y_{BLC1}(z) [mm]",
         pfx.Data(), kBlc2BeamProfileZ, kBlc1ExtrapScanCenter, kBlc1ExtrapScanHalfRange),
    kBlc1ExtrapScanBins, kBlc1ExtrapScanLo, kBlc1ExtrapScanHi,
    kExtrapDyBins, -kExtrapDyRange, kExtrapDyRange);
  h.dy_mean_vs_blc1_extrap_z = new TProfile(
    Form("%s_dy_mean_vs_blc1_extrap_z", detector),
    Form("%s <#Deltay> vs BLC1 extrap z (BLC2 fixed @ %.3f, z=%.0f #pm %.0f);"
         "z_{BLC1 extrap} [mm];<#Deltay> [mm]",
         pfx.Data(), kBlc2BeamProfileZ, kBlc1ExtrapScanCenter, kBlc1ExtrapScanHalfRange),
    kBlc1ExtrapScanBins, kBlc1ExtrapScanLo, kBlc1ExtrapScanHi,
    -kExtrapDyRange, kExtrapDyRange);
  return h;
}

BlcBeamProfileHistSet
make_blc_beam_profile_hist_set(const TString& bcin_cfg)
{
  const TString pfx_bcin = hist_title_prefix("BcIn", bcin_cfg);
  BlcBeamProfileHistSet h;
  h.blc1_xy = make_d5_hist(
    "BcIn_BLC1_beam_xy",
    Form("%s BLC1 beam profile @ z=%.0f mm", pfx_bcin.Data(), kBlc1BeamProfileZ),
    "x [mm]", "y [mm]");
  h.blc2_xy = make_d5_hist(
    "BcOut_BLC2_beam_xy",
    Form("BcOut BLC2 beam profile @ z=%.3f mm", kBlc2BeamProfileZ),
    "x [mm]", "y [mm]");
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

Double_t
blc_y_at_chamber_center(Double_t y0, Double_t v0)
{
  return y0 + v0 * kBlcChamberCenterZ;
}

Double_t
blc_x_at_z(Double_t x0, Double_t u0, Double_t z)
{
  return x0 + u0 * z;
}

Double_t
blc_y_at_z(Double_t y0, Double_t v0, Double_t z)
{
  return y0 + v0 * z;
}

void
fill_blc1_beam_profile(BlcBeamProfileHistSet& h,
                       Double_t x0, Double_t y0, Double_t u0, Double_t v0,
                       RunStats& stats)
{
  const Double_t x = blc_x_at_z(x0, u0, kBlc1BeamProfileZ);
  const Double_t y = blc_y_at_z(y0, v0, kBlc1BeamProfileZ);
  h.blc1_xy->Fill(x, y);
  ++stats.n_filled_blc1_beam_profile;
}

void
fill_blc2_beam_profile(BlcBeamProfileHistSet& h,
                       Double_t x0, Double_t y0, Double_t u0, Double_t v0,
                       RunStats& stats)
{
  const Double_t x = blc_x_at_z(x0, u0, kBlc2BeamProfileZ);
  const Double_t y = blc_y_at_z(y0, v0, kBlc2BeamProfileZ);
  h.blc2_xy->Fill(x, y);
  ++stats.n_filled_blc2_beam_profile;
}

void
fill_blc_y_residual(BlcYResidualHistSet& h,
                    Double_t y0_blc1, Double_t v0_blc1,
                    Double_t y0_blc2, Double_t v0_blc2,
                    RunStats& stats)
{
  const Double_t y_blc1 = y0_blc1;
  const Double_t y_blc2 = blc_y_at_chamber_center(y0_blc2, v0_blc2);
  const Double_t dy = y_blc2 - y_blc1;

  h.dy_at_center->Fill(dy);
  h.v0_vs_dy->Fill(v0_blc1, dy);
  if (std::fabs(v0_blc1) <= kBlc1V0FlatCut) {
    h.dy_at_center_lowv->Fill(dy);
    ++stats.n_filled_blc_y_residual_lowv;
  }
  ++stats.n_filled_blc_y_residual;
}

void
fill_blc_slope(BlcSlopeHistSet& h,
               Double_t y0_blc1, Double_t v0_blc1,
               Double_t y0_blc2, Double_t v0_blc2,
               RunStats& stats)
{
  // Point residual: BLC1@kBlc1ExtrapZ vs BLC2@chamber center (legacy)
  const Double_t y_blc1 = blc_y_at_z(y0_blc1, v0_blc1, kBlc1ExtrapZ);
  const Double_t y_blc2 = blc_y_at_z(y0_blc2, v0_blc2, kBlc2BeamProfileZ);
  const Double_t dy = y_blc2 - y_blc1;

  h.v0_blc1_vs_v0_blc2->Fill(v0_blc1, v0_blc2);
  h.dy_extrap->Fill(dy);
  h.v0_vs_dy_extrap->Fill(v0_blc1, dy);
  h.v0_vs_dy_extrap_zoom->Fill(v0_blc1, dy);

  // BLC1->BLC2: BLC2 y fixed at chamber center; scan BLC1 extrapolation z
  if (h.dy_vs_blc1_extrap_z) {
    for (Int_t iz = 1; iz <= h.dy_vs_blc1_extrap_z->GetNbinsX(); ++iz) {
      const Double_t z = h.dy_vs_blc1_extrap_z->GetXaxis()->GetBinCenter(iz);
      const Double_t dy_z = y_blc2 - blc_y_at_z(y0_blc1, v0_blc1, z);
      h.dy_vs_blc1_extrap_z->Fill(z, dy_z);
      if (h.dy_mean_vs_blc1_extrap_z) h.dy_mean_vs_blc1_extrap_z->Fill(z, dy_z);
    }
  }

  ++stats.n_filled_blc_slope;
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
draw_hist1_or_placeholder(TH1D* h, const char* opt = "")
{
  if (!h || h->GetEntries() <= 0) {
    draw_placeholder(h ? "empty histogram" : "histogram missing");
    return;
  }
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.08);
  h->GetXaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleOffset(1.45);
  h->Draw(opt);
}

BlcYResidualGaussFit
fit_dy_gauss(TH1D* h)
{
  BlcYResidualGaussFit fit;
  if (!h || h->GetEntries() < 20) return fit;

  const Double_t mean = h->GetMean();
  const Double_t rms  = h->GetRMS();
  if (rms <= 0.0) return fit;

  const Double_t xmin = mean - 3.0 * rms;
  const Double_t xmax = mean + 3.0 * rms;
  auto* gaus = new TF1(Form("gaus_%s", h->GetName()), "gaus", xmin, xmax);
  gaus->SetLineColor(kRed);
  gaus->SetParameters(h->GetMaximum(), mean, rms);

  if (h->Fit(gaus, "RSQ0") != 0) return fit;

  fit.ok = true;
  fit.amp = gaus->GetParameter(0);
  fit.mean = gaus->GetParameter(1);
  fit.sigma = gaus->GetParameter(2);
  fit.mean_err = gaus->GetParError(1);
  fit.sigma_err = gaus->GetParError(2);
  return fit;
}

void
draw_dy_with_gauss_fit(TH1D* h, const BlcYResidualGaussFit& fit, Int_t color = kBlue + 1,
                       const char* third_line = "BLC1 y shift #approx +%.3f mm (fix BLC2)")
{
  if (!h || h->GetEntries() <= 0) {
    draw_placeholder(h ? "empty histogram" : "histogram missing");
    return;
  }

  h->SetLineColor(color);
  h->SetLineWidth(2);

  if (fit.ok && fit.sigma > 0.0) {
    Double_t half = kGaussDrawNSigma * std::fabs(fit.sigma);
    if (half < kGaussDrawMinHalfRange) half = kGaussDrawMinHalfRange;
    const Double_t xmin = fit.mean - half;
    const Double_t xmax = fit.mean + half;
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    auto* gaus = h->GetFunction(Form("gaus_%s", h->GetName()));
    if (gaus) gaus->SetRange(xmin, xmax);
  }

  draw_hist1_or_placeholder(h);

  if (!fit.ok) return;

  auto* gaus = h->GetFunction(Form("gaus_%s", h->GetName()));
  if (!gaus) return;
  gaus->SetLineColor(kRed);
  gaus->SetLineWidth(2);
  gaus->Draw("same");

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(12);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.16, 0.82,
                 Form("Gauss fit: #mu = %.3f #pm %.3f mm", fit.mean, fit.mean_err));
  tex->DrawLatex(0.16, 0.76,
                 Form("#sigma = %.3f #pm %.3f mm", fit.sigma, fit.sigma_err));
  if (third_line && third_line[0] != '\0') {
    tex->DrawLatex(0.16, 0.70, Form(third_line, fit.mean));
  }
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
                const TString& bcin_cfg, Int_t beam_flag_filter)
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
  Double_t title_y = 0.78;
  if (beam_flag_filter >= 0) {
    tex->SetTextColor(kBlue + 1);
    tex->DrawLatex(0.5, title_y,
                   Form("beam_flag filter: %s", beam_flag_filter_label(beam_flag_filter).Data()));
    tex->SetTextColor(kBlack);
    title_y -= 0.06;
  }
  if (swap_bcin_xy) {
    tex->SetTextColor(kRed);
    tex->DrawLatex(0.5, title_y, "BcIn: x0<->y0, u0<->v0 swapped (diagnostic)");
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
  if (beam_flag_filter >= 0) {
    tex->DrawLatex(0.5, 0.44,
                   Form("Skipped by beam_flag filter: %lld", stats.n_skipped_beam_flag));
  }
  tex->DrawLatex(0.5, 0.38, Form("Synced BH2/BHT fills: %lld / %lld",
                                 stats.n_filled_bh2, stats.n_filled_bht));
  tex->DrawLatex(0.5, 0.32, Form("D5@BH2 fills: %lld, D5 vs BLC2: %lld",
                                 stats.n_filled_d5_bh2, stats.n_filled_d5_blc2));
  tex->DrawLatex(0.5, 0.24, Form("D5 z_in=%.1f mm, z_out=%.1f mm", mtx.z_in, mtx.z_out));
  tex->DrawLatex(0.5, 0.18, Form("z_BH2=%.1f mm, z_BHT=%.1f mm", kBh2Z, kBhtZ));
  tex->DrawLatex(0.5, 0.10, Form("created: %s", now.AsString("s")));
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
draw_blc_y_residual_pages(TCanvas* c, const TString& pdf_path,
                          const BlcYResidualHistSet& h, const TString& bcin_cfg)
{
  const BlcYResidualGaussFit fit_all  = fit_dy_gauss(h.dy_at_center);
  const BlcYResidualGaussFit fit_lowv = fit_dy_gauss(h.dy_at_center_lowv);

  if (fit_all.ok) {
    std::cout << "BLC y residual gauss fit [all]: mean=" << fit_all.mean
              << " +/- " << fit_all.mean_err << " mm"
              << " => shift BLC1 y by +" << fit_all.mean << " mm (fix BLC2)" << std::endl;
  }
  if (fit_lowv.ok) {
    std::cout << "BLC y residual gauss fit [|BLC1 v0|<=" << kBlc1V0FlatCut
              << "]: mean=" << fit_lowv.mean
              << " +/- " << fit_lowv.mean_err << " mm"
              << " => shift BLC1 y by +" << fit_lowv.mean << " mm (fix BLC2)" << std::endl;
  }

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_dy_with_gauss_fit(h.dy_at_center, fit_all, kBlack);
  c->cd(2);
  draw_hist_or_placeholder(h.v0_vs_dy);
  print_page(c, pdf_path);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_dy_with_gauss_fit(h.dy_at_center, fit_all, kBlack);
  c->cd(2);
  draw_dy_with_gauss_fit(h.dy_at_center_lowv, fit_lowv, kBlue + 1);
  print_page(c, pdf_path);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.90, Form("BLC1 vs BLC2 y residual @ chamber center [%s]",
                                 bcin_cfg.Data()));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.82,
                 Form("#Deltay = (y0+v0*z)_{BLC2} - y0_{BLC1} at z=%.1f mm (chamber center).",
                      kBlcChamberCenterZ));
  tex->DrawLatex(0.5, 0.76,
                 "Gauss #mu: offset of BLC1 y when BLC2 is fixed (#Deltay = y_{BLC2}-y_{BLC1}).");
  tex->DrawLatex(0.5, 0.70,
                 Form("Apply BLC1 y shift #approx +#mu (same sign as fit mean)."));
  tex->DrawLatex(0.5, 0.64,
                 Form("Flat-track subset: |BLC1 v0| #leq %.3f", kBlc1V0FlatCut));
  if (fit_all.ok) {
    tex->DrawLatex(0.5, 0.56,
                   Form("All tracks: #mu = %.3f #pm %.3f mm, #sigma = %.3f #pm %.3f mm",
                        fit_all.mean, fit_all.mean_err, fit_all.sigma, fit_all.sigma_err));
  } else {
    tex->DrawLatex(0.5, 0.56, "All tracks: gauss fit failed / too few entries");
  }
  if (fit_lowv.ok) {
    tex->DrawLatex(0.5, 0.48,
                   Form("|v0| cut: #mu = %.3f #pm %.3f mm, #sigma = %.3f #pm %.3f mm",
                        fit_lowv.mean, fit_lowv.mean_err, fit_lowv.sigma, fit_lowv.sigma_err));
  } else {
    tex->DrawLatex(0.5, 0.48, Form("|v0| cut: gauss fit failed / too few entries"));
  }
  tex->DrawLatex(0.5, 0.38, Form("entries: all=%.0f, |v0|#leq%.3f=%.0f, v0_vs_dy=%.0f",
                                 h.dy_at_center ? h.dy_at_center->GetEntries() : 0,
                                 kBlc1V0FlatCut,
                                 h.dy_at_center_lowv ? h.dy_at_center_lowv->GetEntries() : 0,
                                 h.v0_vs_dy ? h.v0_vs_dy->GetEntries() : 0));
  print_page(c, pdf_path);
}

void
draw_blc_slope_pages(TCanvas* c, const TString& pdf_path,
                     const BlcSlopeHistSet& h, const TString& bcin_cfg)
{
  const BlcYResidualGaussFit fit_extrap = fit_dy_gauss(h.dy_extrap);
  if (fit_extrap.ok) {
    std::cout << "BLC extrap y residual gauss fit [BLC1@z=" << kBlc1ExtrapZ
              << ", BLC2@z=" << kBlc2BeamProfileZ << "]: mean=" << fit_extrap.mean
              << " +/- " << fit_extrap.mean_err << " mm"
              << " => shift BLC1 y by +" << fit_extrap.mean << " mm (fix BLC2)"
              << std::endl;
  }

  Double_t dy_at_scan_center = 0.0;
  bool dy_at_scan_center_ok = false;
  if (h.dy_mean_vs_blc1_extrap_z && h.dy_mean_vs_blc1_extrap_z->GetEntries() > 0) {
    const Int_t ib = h.dy_mean_vs_blc1_extrap_z->FindBin(kBlc1ExtrapScanCenter);
    if (ib >= 1 && ib <= h.dy_mean_vs_blc1_extrap_z->GetNbinsX()
        && h.dy_mean_vs_blc1_extrap_z->GetBinEntries(ib) > 0) {
      dy_at_scan_center = h.dy_mean_vs_blc1_extrap_z->GetBinContent(ib);
      dy_at_scan_center_ok = true;
      std::cout << "BLC1->BLC2 <#Deltay> @ BLC1 extrap z=" << kBlc1ExtrapScanCenter
                << " mm (BLC2 fixed @ " << kBlc2BeamProfileZ << "): "
                << dy_at_scan_center << " mm" << std::endl;
    }
  }

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_hist_or_placeholder(h.v0_blc1_vs_v0_blc2);
  c->cd(2);
  draw_dy_with_gauss_fit(h.dy_extrap, fit_extrap, kBlack);
  print_page(c, pdf_path);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_hist_or_placeholder(h.v0_vs_dy_extrap);
  c->cd(2);
  draw_dy_with_gauss_fit(h.dy_extrap, fit_extrap, kBlack);
  print_page(c, pdf_path);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_hist_or_placeholder(h.v0_vs_dy_extrap_zoom);
  c->cd(2);
  draw_dy_with_gauss_fit(h.dy_extrap, fit_extrap, kBlack);
  print_page(c, pdf_path);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_hist_or_placeholder(h.dy_vs_blc1_extrap_z);
  c->cd(2);
  if (!h.dy_mean_vs_blc1_extrap_z || h.dy_mean_vs_blc1_extrap_z->GetEntries() <= 0) {
    draw_placeholder(h.dy_mean_vs_blc1_extrap_z ? "empty histogram" : "histogram missing");
  } else {
    h.dy_mean_vs_blc1_extrap_z->SetMarkerStyle(20);
    h.dy_mean_vs_blc1_extrap_z->SetMarkerSize(0.7);
    h.dy_mean_vs_blc1_extrap_z->SetLineColor(kBlue + 1);
    h.dy_mean_vs_blc1_extrap_z->Draw("E1");
  }
  print_page(c, pdf_path);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.90, Form("BLC1/BLC2 v0 correlation and extrap y residual [%s]",
                                 bcin_cfg.Data()));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.80,
                 Form("Point extrap: y_{BLC1}@z=%.0f mm; y_{BLC2}@z=%.3f mm.",
                      kBlc1ExtrapZ, kBlc2BeamProfileZ));
  tex->DrawLatex(0.5, 0.74,
                 "#Deltay = y_{BLC2} - y_{BLC1} (same sign as chamber-center residual).");
  tex->DrawLatex(0.5, 0.68, "Gauss #mu on #Deltay: BLC1 y offset when BLC2 is fixed.");
  tex->DrawLatex(0.5, 0.62,
                 Form("Close-up 2D: |v0| #leq %.3f, |#Deltay| #leq %.0f mm.",
                      kExtrapZoomV0, kExtrapZoomDy));
  tex->DrawLatex(0.5, 0.56,
                 Form("BLC1->BLC2 scan: BLC2 fixed @ %.3f; vary z_{BLC1} #in [%.0f, %.0f] (%.0f#pm%.0f).",
                      kBlc2BeamProfileZ, kBlc1ExtrapScanLo, kBlc1ExtrapScanHi,
                      kBlc1ExtrapScanCenter, kBlc1ExtrapScanHalfRange));
  if (dy_at_scan_center_ok) {
    tex->DrawLatex(0.5, 0.50,
                   Form("<#Deltay> @ z_{BLC1}=%.0f: %.3f mm",
                        kBlc1ExtrapScanCenter, dy_at_scan_center));
  } else {
    tex->DrawLatex(0.5, 0.50, "<#Deltay> @ z_{BLC1}=2900: n/a");
  }
  if (fit_extrap.ok) {
    tex->DrawLatex(0.5, 0.44,
                   Form("Point extrap #Deltay: #mu = %.3f #pm %.3f mm, #sigma = %.3f #pm %.3f mm",
                        fit_extrap.mean, fit_extrap.mean_err,
                        fit_extrap.sigma, fit_extrap.sigma_err));
    tex->DrawLatex(0.5, 0.38,
                   Form("BLC1 y shift #approx +%.3f mm (fix BLC2)", fit_extrap.mean));
  } else {
    tex->DrawLatex(0.5, 0.44, "Point extrap #Deltay: gauss fit failed / too few entries");
  }
  tex->DrawLatex(0.5, 0.30, Form("entries: v0_vs_v0=%.0f, dy=%.0f, v0_vs_dy=%.0f, zoom=%.0f, dy_vs_z=%.0f",
                                 h.v0_blc1_vs_v0_blc2 ? h.v0_blc1_vs_v0_blc2->GetEntries() : 0,
                                 h.dy_extrap ? h.dy_extrap->GetEntries() : 0,
                                 h.v0_vs_dy_extrap ? h.v0_vs_dy_extrap->GetEntries() : 0,
                                 h.v0_vs_dy_extrap_zoom ? h.v0_vs_dy_extrap_zoom->GetEntries() : 0,
                                 h.dy_vs_blc1_extrap_z ? h.dy_vs_blc1_extrap_z->GetEntries() : 0));
  print_page(c, pdf_path);
}

void
draw_blc_beam_profile_pages(TCanvas* c, const TString& pdf_path,
                            const BlcBeamProfileHistSet& h, const TString& bcin_cfg)
{
  draw_pair_page(c, pdf_path, h.blc1_xy, h.blc2_xy);

  TH1D* blc1_x = h.blc1_xy ? h.blc1_xy->ProjectionX("blc1_beam_px") : nullptr;
  TH1D* blc1_y = h.blc1_xy ? h.blc1_xy->ProjectionY("blc1_beam_py") : nullptr;
  TH1D* blc2_x = h.blc2_xy ? h.blc2_xy->ProjectionX("blc2_beam_px") : nullptr;
  TH1D* blc2_y = h.blc2_xy ? h.blc2_xy->ProjectionY("blc2_beam_py") : nullptr;

  if (blc1_x) {
    blc1_x->SetTitle(Form("BLC1 x @ z=%.0f mm;x [mm];entries", kBlc1BeamProfileZ));
  }
  if (blc1_y) {
    blc1_y->SetTitle(Form("BLC1 y @ z=%.0f mm;y [mm];entries", kBlc1BeamProfileZ));
  }
  if (blc2_x) {
    blc2_x->SetTitle(Form("BLC2 x @ z=%.3f mm;x [mm];entries", kBlc2BeamProfileZ));
  }
  if (blc2_y) {
    blc2_y->SetTitle(Form("BLC2 y @ z=%.3f mm;y [mm];entries", kBlc2BeamProfileZ));
  }

  const BlcYResidualGaussFit fit_blc1_x = fit_dy_gauss(blc1_x);
  const BlcYResidualGaussFit fit_blc1_y = fit_dy_gauss(blc1_y);
  const BlcYResidualGaussFit fit_blc2_x = fit_dy_gauss(blc2_x);
  const BlcYResidualGaussFit fit_blc2_y = fit_dy_gauss(blc2_y);

  if (fit_blc1_x.ok || fit_blc1_y.ok) {
    std::cout << "BLC1 beam center @ z=" << kBlc1BeamProfileZ << " mm:";
    if (fit_blc1_x.ok) {
      std::cout << " x=" << fit_blc1_x.mean << " +/- " << fit_blc1_x.mean_err << " mm";
    }
    if (fit_blc1_y.ok) {
      std::cout << " y=" << fit_blc1_y.mean << " +/- " << fit_blc1_y.mean_err << " mm";
    }
    std::cout << std::endl;
  }
  if (fit_blc2_x.ok || fit_blc2_y.ok) {
    std::cout << "BLC2 beam center @ z=" << kBlc2BeamProfileZ << " mm:";
    if (fit_blc2_x.ok) {
      std::cout << " x=" << fit_blc2_x.mean << " +/- " << fit_blc2_x.mean_err << " mm";
    }
    if (fit_blc2_y.ok) {
      std::cout << " y=" << fit_blc2_y.mean << " +/- " << fit_blc2_y.mean_err << " mm";
    }
    std::cout << std::endl;
  }

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_dy_with_gauss_fit(blc1_x, fit_blc1_x, kBlack, "beam center x #approx %.3f mm");
  c->cd(2);
  draw_dy_with_gauss_fit(blc1_y, fit_blc1_y, kBlue + 1, "beam center y #approx %.3f mm");
  print_page(c, pdf_path);

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  draw_dy_with_gauss_fit(blc2_x, fit_blc2_x, kBlack, "beam center x #approx %.3f mm");
  c->cd(2);
  draw_dy_with_gauss_fit(blc2_y, fit_blc2_y, kBlue + 1, "beam center y #approx %.3f mm");
  print_page(c, pdf_path);

  c->Clear();
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.038);
  tex->DrawLatex(0.5, 0.90, Form("BLC1 / BLC2 beam profile + center fit [%s]",
                                 bcin_cfg.Data()));
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.5, 0.82,
                 Form("BLC1 (BcIn): x0+u0*z, y0+v0*z @ z=%.0f mm.",
                      kBlc1BeamProfileZ));
  tex->DrawLatex(0.5, 0.76,
                 Form("BLC2 (BcOut): x0+u0*z, y0+v0*z @ z=%.3f mm.",
                      kBlc2BeamProfileZ));
  tex->DrawLatex(0.5, 0.70,
                 "Beam center from Gauss fit on ProjectionX / ProjectionY.");
  if (fit_blc1_x.ok && fit_blc1_y.ok) {
    tex->DrawLatex(0.5, 0.62,
                   Form("BLC1 center: (%.3f #pm %.3f, %.3f #pm %.3f) mm",
                        fit_blc1_x.mean, fit_blc1_x.mean_err,
                        fit_blc1_y.mean, fit_blc1_y.mean_err));
    tex->DrawLatex(0.5, 0.56,
                   Form("BLC1 #sigma_{x}=%.3f, #sigma_{y}=%.3f mm",
                        fit_blc1_x.sigma, fit_blc1_y.sigma));
  } else {
    tex->DrawLatex(0.5, 0.62, "BLC1 center: gauss fit failed / too few entries");
  }
  if (fit_blc2_x.ok && fit_blc2_y.ok) {
    tex->DrawLatex(0.5, 0.48,
                   Form("BLC2 center: (%.3f #pm %.3f, %.3f #pm %.3f) mm",
                        fit_blc2_x.mean, fit_blc2_x.mean_err,
                        fit_blc2_y.mean, fit_blc2_y.mean_err));
    tex->DrawLatex(0.5, 0.42,
                   Form("BLC2 #sigma_{x}=%.3f, #sigma_{y}=%.3f mm",
                        fit_blc2_x.sigma, fit_blc2_y.sigma));
  } else {
    tex->DrawLatex(0.5, 0.48, "BLC2 center: gauss fit failed / too few entries");
  }
  tex->DrawLatex(0.5, 0.32, Form("entries: BLC1=%.0f, BLC2=%.0f",
                                 h.blc1_xy ? h.blc1_xy->GetEntries() : 0,
                                 h.blc2_xy ? h.blc2_xy->GetEntries() : 0));
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
        bool swap_bcin_xy, Int_t beam_flag_filter)
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
  if (beam_flag_filter >= 0) {
    std::cout << "beam_flag filter: " << beam_flag_filter_label(beam_flag_filter)
              << std::endl;
  }

  const TString bcin_cfg = bcin_config_label(path_bcin);

  const auto bcin_idx = build_event_index(tree_bcin);
  const auto bcout_idx = build_event_index(tree_bcout);

  Bh2HistSet h_bcin_bh2 = make_bh2_hist_set("BcIn", bcin_cfg);
  Bh2HistSet h_bcout_bh2 = make_bh2_hist_set("BcOut", bcin_cfg);
  BhtHistSet h_bcin_bht = make_bht_hist_set("BcIn", bcin_cfg);
  BhtHistSet h_bcout_bht = make_bht_hist_set("BcOut", bcin_cfg);
  Bh2D5HistSet h_bcin_d5_bh2 = make_bh2_d5_hist_set("BcIn", bcin_cfg);
  D5Blc2HistSet h_bcin_d5_blc2 = make_d5_blc2_hist_set("BcIn", bcin_cfg);
  BlcYResidualHistSet h_blc_y_residual = make_blc_y_residual_hist_set("BcIn", bcin_cfg);
  BlcSlopeHistSet h_blc_slope = make_blc_slope_hist_set("BcIn", bcin_cfg);
  BlcBeamProfileHistSet h_blc_beam_profile = make_blc_beam_profile_hist_set(bcin_cfg);
  RunStats stats;
  stats.n_hodo = tree_hodo->GetEntries();

  UInt_t event_number = 0;
  Int_t beam_flag = 0;
  std::vector<Double_t>* bh2_hit_seg = nullptr;
  std::vector<Double_t>* bht_hit_seg = nullptr;
  tree_hodo->SetBranchAddress("event_number", &event_number);
  tree_hodo->SetBranchAddress("beam_flag", &beam_flag);
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

    if (beam_flag_filter >= 0 && beam_flag != beam_flag_filter) {
      ++stats.n_skipped_beam_flag;
      continue;
    }

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
      fill_blc_y_residual(h_blc_y_residual,
                          (*y0_in)[0], (*v0_in)[0],
                          (*y0_bcout)[0], (*v0_bcout)[0],
                          stats);
      fill_blc_slope(h_blc_slope,
                     (*y0_in)[0], (*v0_in)[0],
                     (*y0_bcout)[0], (*v0_bcout)[0],
                     stats);
    }

    if (ntrack_bcin >= 1 && x0_in && y0_in && u0_in && v0_in) {
      fill_blc1_beam_profile(h_blc_beam_profile,
                             (*x0_in)[0], (*y0_in)[0], (*u0_in)[0], (*v0_in)[0],
                             stats);
    }
    if (ntrack_bcout >= 1 && x0_bcout && y0_bcout && u0_bcout && v0_bcout) {
      fill_blc2_beam_profile(h_blc_beam_profile,
                             (*x0_bcout)[0], (*y0_bcout)[0],
                             (*u0_bcout)[0], (*v0_bcout)[0],
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
                  swap_bcin_xy, mtx, bcin_cfg, beam_flag_filter);
  draw_bh2_pages(canvas, output_path, Form("BcIn [%s]", bcin_cfg.Data()), h_bcin_bh2);
  draw_bht_pages(canvas, output_path, Form("BcIn [%s]", bcin_cfg.Data()), h_bcin_bht);
  draw_d5_bh2_pages(canvas, output_path, h_bcin_d5_bh2, bcin_cfg);
  draw_d5_blc2_pages(canvas, output_path, h_bcin_d5_blc2, bcin_cfg);
  draw_blc_y_residual_pages(canvas, output_path, h_blc_y_residual, bcin_cfg);
  draw_blc_slope_pages(canvas, output_path, h_blc_slope, bcin_cfg);
  draw_blc_beam_profile_pages(canvas, output_path, h_blc_beam_profile, bcin_cfg);
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
            << ", fills_d5_blc2=" << stats.n_filled_d5_blc2
            << ", fills_blc_y_residual=" << stats.n_filled_blc_y_residual
            << ", fills_blc_y_residual_lowv=" << stats.n_filled_blc_y_residual_lowv
            << ", fills_blc_slope=" << stats.n_filled_blc_slope
            << ", fills_blc1_beam=" << stats.n_filled_blc1_beam_profile
            << ", fills_blc2_beam=" << stats.n_filled_blc2_beam_profile << ")" << std::endl;
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
  Int_t beam_flag_filter = kBeamFlagNoFilter;

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
    if (strcmp(arg, "--pion") == 0) {
      if (beam_flag_filter != kBeamFlagNoFilter) {
        std::cerr << "Specify only one of --pion or --kaon." << std::endl;
        usage(argv[0]);
        return 1;
      }
      beam_flag_filter = kBeamFlagPion;
      continue;
    }
    if (strcmp(arg, "--kaon") == 0) {
      if (beam_flag_filter != kBeamFlagNoFilter) {
        std::cerr << "Specify only one of --pion or --kaon." << std::endl;
        usage(argv[0]);
        return 1;
      }
      beam_flag_filter = kBeamFlagKaon;
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
  return analyze(run_num, path_bcin, has_output ? output_path : TString(), swap_bcin_xy,
                 beam_flag_filter);
}
