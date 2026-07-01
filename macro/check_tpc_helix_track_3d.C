// -*- C++ -*-
// DstTPCHelixTracking 出力（ツリー "tpc"）の 3D 目視レビュー。ベースは check_tpc_track_3d.C。
//
// Usage:
//   root -l
//   .L check_tpc_helix_track_3d.C+
//   helix_set_path("file.root")
//   helix_event()              // 1 イベント（ランダム）
//   helix_event(123)           // entry 固定
//   helix_status()             // 現在のフィルタ設定
//   helix_coords_help()        // 座標系と ROOT 3D 軸の対応表
//   .L 読み込み時に空の TPC フレーム窓を開く（helix_open_gui）。データは helix_set_path 後に helix_event()
//   SetMomRange(0.3, 0.8); SetMomTrackIndex(1);
//   SetPidFilter(0x6, 0x2);    // (pid & 0x6) == 0x2  → K ビット必須
//   SetRequireK0(kTRUE);
//   SetCloseDistMax(20); SetNtTpcRange(2, 2);
//   helix_next(50000)          // フィルタに合う次の entry を表示
//   helix_browse_filtered()    // Enter=helix_next(), q=終了
//   helix_reset_config()       // フィルタをデフォルトに
//   SetHelixViewMargin(0.25)   // 3D 引き（余白。既定 0.22）
//   SetHelixPreferHitTheta(kTRUE)  // 螺旋 θ 区間を helix_t 優先（ビームの長い線を抑制）
//   SetHelixDisplayTargetOrigin(kTRUE)  // 縦軸 z をターゲット中心 (z'=0) で表示
//   SetHelixDrawCalpos(kTRUE)      // calpos_* を小十字で重ね描画（座標検証）
//   SetHelixDrawLayerRings(kTRUE)  // 底面 (y=-310) に TPC レイヤー円（既定 ON）
//
// 必須枝: run_number, event_number, ntTpc, nhtrack, pid,
//         helix_cx, helix_cy, helix_z0, helix_r, helix_dz, helix_t,
//         hitpos_x, hitpos_y, hitpos_z
//
// 座標（3 系 + ROOT 描画の軸入れ替え）— helix_coords_help() でも表示
//
//  [A] 螺旋フィット内部 (helix_cx, helix_cy, …)  … DST の helix_* パラメータ
//      水平面はターゲット中心基準 (internal y = z_tpc - Z_TARGET)。そのまま描かない。
//
//  [B] TPC 実座標 (x, y, z)  … hitpos_*, calpos_*, vtxTpc（アナライザの LocalToGlobal 後）
//      x: パッド平面の横方向 [mm]
//      y: ドリフト方向（行深さ） [mm]
//      z: ビーム方向 [mm]、原点 z=0 は TPC 機械中心、ターゲットは z = Z_TARGET = -143 mm
//      螺旋線は [A]→[B] 変換後に描く: x=-xi, y=zi, z=yi+Z_TARGET
//
//  [C] GEANT / DC の実験グローバル … この DST 枝には無い（本マクロの対象外）
//
//  ROOT TPolyLine3D::SetPoint(i, X, Y, Z) への渡し方（check_tpc_track_3d.C と同じ慣習）:
//      X = z_tpc  （オプションで z' = z_tpc - Z_TARGET、SetHelixDisplayTargetOrigin）
//      Y = x_tpc
//      Z = y_tpc
//  角の矢印ラベルは TPC 軸名: 赤=z, 緑=x, 青=y（ROOT 標準の X/Y/Z ラベルではない）
//
//  フレーム: 八角形は TPC 機械中心 (x,z)=(0,0) 基準。ターゲット円は z=Z_TARGET に固定。
// PID 色: K=青, p=赤, pi=緑（複数ビット時 K > p > pi）。

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

static constexpr Double_t kTPC_Z_TARGET = -143.0;

TFile* gHlxFile = nullptr;
TTreeReader* gHlxReader = nullptr;
TRandom3* gHlxRandom = nullptr;
TCanvas* gHlxCanvas = nullptr;
TView* gHlxView = nullptr;
Long64_t gHlxCurrentEntry = -1;

static Int_t gHlxHelixTRangeMode = 0;
// TView 余白（各軸方向に range の何割を足すか）。大きいほど引きの構図。
static Double_t gHlxViewMargin = 0.22;
// 螺旋描画: kTRUE なら helix_t（ヒット列）から θ 区間を優先（ビーム等の長い extrap 線を避ける）
static Bool_t gHlxHelixPreferHitTheta = kTRUE;
// kTRUE: 表示の X 軸に z' = z_tpc - Z_TARGET を使う（既定は kFALSE = z_tpc そのまま）
static Bool_t gHlxDisplayTargetOrigin = kFALSE;
// calpos_* 枝があれば螺旋上の参照点を重ね描画（座標変換の目視検証）
static Bool_t gHlxDrawCalpos = kTRUE;
// 底面に TPC パッドレイヤー半径の円を描画（ターゲット中心 z=Z_TARGET 基準）
static Bool_t gHlxDrawLayerRings = kTRUE;

// TPCPadHelper::padParameter[][kRadius] と同じ（マクロ単体で .L 可能にするため埋め込み）
static const Int_t kHelixNumTpcLayers = 32;
static const Double_t kHelixTpcLayerRadius[kHelixNumTpcLayers] = {
  14.75, 24.25, 33.75, 43.25, 52.75, 62.25, 71.75, 81.25, 90.75, 100.25, 111.5, 124.5, 137.5,
  150.5, 163.5, 176.5, 189.5, 202.5, 215.5, 228.5, 241.5, 254.5, 267.5, 280.5, 293.5, 306.5,
  319.5, 332.5, 345.5, 358.5, 371.5, 384.5};

struct HelixBranchFlags {
  Bool_t has_is_beam = kFALSE;
  Bool_t has_is_accidental = kFALSE;
  Bool_t has_charge = kFALSE;
  Bool_t has_mom0 = kFALSE;
  Bool_t has_dEdx = kFALSE;
  Bool_t has_beam_flag = kFALSE;
  Bool_t has_effective_ntTpc = kFALSE;
  Bool_t has_closeDistTpc = kFALSE;
  Bool_t has_vtxTpc = kFALSE;
  Bool_t has_calpos = kFALSE;
  Bool_t has_k0_mass = kFALSE;
  Bool_t has_lambda_mass = kFALSE;
};

HelixBranchFlags gHlxBranches;

struct HelixReviewConfig {
  Double_t mom_min = -1.0;
  Double_t mom_max = -1.0;
  Int_t track_for_mom = -1;
  UInt_t pid_mask = 0;
  UInt_t pid_require = 0;
  Bool_t skip_accidental_tracks = kFALSE;
  Bool_t require_beam_track = kFALSE;
  Int_t min_ntTpc = 0;
  Int_t max_ntTpc = 999;
  Double_t close_dist_max = -1.0;
  Bool_t require_k0 = kFALSE;
  Double_t k0_mass_min = 0.45;
  Double_t k0_mass_max = 0.55;
  Long64_t scan_start = 0;
};

HelixReviewConfig gHlxCfg;

struct HelixEventData {
  UInt_t runnum = 0;
  UInt_t evnum = 0;
  Int_t ntTpc = 0;
  Int_t effective_ntTpc = 0;
  Int_t beam_flag = -1;
  std::vector<Int_t> nhtrack;
  std::vector<Int_t> pid;
  std::vector<Int_t> is_beam;
  std::vector<Int_t> is_accidental;
  std::vector<Int_t> charge;
  std::vector<Double_t> mom0;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> helix_t_min;
  std::vector<Double_t> helix_t_max;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> closeDistTpc;
  std::vector<std::vector<Double_t>> vtxTpc;
  std::vector<std::vector<Double_t>> vtyTpc;
  std::vector<std::vector<Double_t>> vtzTpc;
  std::vector<Double_t> k0_mass;
  std::vector<Double_t> lambda_mass;
};

HelixEventData gHlxEvent;

//______________________________________________________________________________
static std::string
HelixFormatSigFig(Double_t val, Int_t nsig = 3)
{
  if(!TMath::Finite(val))
    return "nan";
  char buf[48];
  snprintf(buf, sizeof(buf), "%.*g", nsig, val);
  return buf;
}

static std::string
HelixTrackSubscript(Int_t track)
{
  std::ostringstream ss;
  ss << "tr_{" << track << "}";
  return ss.str();
}

std::vector<TPolyLine3D*> gHlxHelixLines;
std::vector<TPolyMarker3D*> gHlxHitMarkers;
std::vector<TPolyLine3D*> gHlxFrame;
std::vector<TObject*> gHlxOverlay;

//______________________________________________________________________________
static std::string
DecodePidCandidates(Int_t pid_code)
{
  if(pid_code == 0)
    return "e";
  std::string names;
  if(pid_code & 0x1)
    names += "pi";
  if(pid_code & 0x2) {
    if(!names.empty())
      names += "+";
    names += "K";
  }
  if(pid_code & 0x4) {
    if(!names.empty())
      names += "+";
    names += "p";
  }
  if(names.empty())
    names = "none";
  return names;
}

//______________________________________________________________________________
static Double_t
HelixMinCloseDist(const HelixEventData& ev)
{
  Double_t best = std::numeric_limits<Double_t>::max();
  if(ev.closeDistTpc.size() >= 2U) {
    if(ev.closeDistTpc[0].size() > 1U) {
      const Double_t d = ev.closeDistTpc[0][1];
      if(TMath::Finite(d))
        best = TMath::Min(best, TMath::Abs(d));
    }
    if(ev.closeDistTpc[1].size() > 0U) {
      const Double_t d = ev.closeDistTpc[1][0];
      if(TMath::Finite(d))
        best = TMath::Min(best, TMath::Abs(d));
    }
  }
  if(best >= std::numeric_limits<Double_t>::max() * 0.5)
    return -1.0;
  return best;
}

//______________________________________________________________________________
static Bool_t
HelixEventHasAccidentalTrack(const HelixEventData& ev)
{
  for(Int_t ib : ev.is_accidental) {
    if(ib != 0)
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
static Bool_t
HelixEventHasBeamTrack(const HelixEventData& ev)
{
  for(Int_t ib : ev.is_beam) {
    if(ib != 0)
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
static Bool_t
HelixMomInRange(const HelixEventData& ev, Double_t lo, Double_t hi, Int_t track_index)
{
  if(lo < 0.0 && hi < 0.0)
    return kTRUE;
  if(ev.mom0.empty())
    return kFALSE;
  auto check_one = [&](Int_t it) -> Bool_t {
    if(it < 0 || it >= static_cast<Int_t>(ev.mom0.size()))
      return kFALSE;
    const Double_t m = ev.mom0[it];
    if(!TMath::Finite(m))
      return kFALSE;
    if(lo >= 0.0 && m < lo)
      return kFALSE;
    if(hi >= 0.0 && m > hi)
      return kFALSE;
    return kTRUE;
  };
  if(track_index >= 0)
    return check_one(track_index);
  for(Int_t it = 0; it < static_cast<Int_t>(ev.mom0.size()); ++it) {
    if(check_one(it))
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
static Bool_t
HelixK0InWindow(const HelixEventData& ev, Double_t lo, Double_t hi)
{
  for(Double_t m : ev.k0_mass) {
    if(TMath::Finite(m) && m >= lo && m <= hi)
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
Bool_t
helix_matches_filters(const HelixEventData& ev)
{
  if(ev.ntTpc < gHlxCfg.min_ntTpc || ev.ntTpc > gHlxCfg.max_ntTpc)
    return kFALSE;

  if(gHlxCfg.skip_accidental_tracks && gHlxBranches.has_is_accidental && HelixEventHasAccidentalTrack(ev))
    return kFALSE;

  if(gHlxCfg.require_beam_track) {
    if(!gHlxBranches.has_is_beam)
      return kFALSE;
    if(!HelixEventHasBeamTrack(ev))
      return kFALSE;
  }

  if((gHlxCfg.mom_min >= 0.0 || gHlxCfg.mom_max >= 0.0) && gHlxBranches.has_mom0) {
    if(!HelixMomInRange(ev, gHlxCfg.mom_min, gHlxCfg.mom_max, gHlxCfg.track_for_mom))
      return kFALSE;
  } else if((gHlxCfg.mom_min >= 0.0 || gHlxCfg.mom_max >= 0.0) && !gHlxBranches.has_mom0) {
    return kFALSE;
  }

  if(gHlxCfg.pid_mask != 0U) {
    Bool_t pid_ok = kFALSE;
    for(Int_t it = 0; it < ev.ntTpc && it < static_cast<Int_t>(ev.pid.size()); ++it) {
      const UInt_t pv = static_cast<UInt_t>(ev.pid[it]);
      if((pv & gHlxCfg.pid_mask) == gHlxCfg.pid_require) {
        pid_ok = kTRUE;
        break;
      }
    }
    if(!pid_ok)
      return kFALSE;
  }

  if(gHlxCfg.close_dist_max >= 0.0) {
    if(!gHlxBranches.has_closeDistTpc || ev.ntTpc != 2)
      return kFALSE;
    const Double_t d = HelixMinCloseDist(ev);
    if(d < 0.0 || d > gHlxCfg.close_dist_max)
      return kFALSE;
  }

  if(gHlxCfg.require_k0) {
    if(!gHlxBranches.has_k0_mass || ev.k0_mass.empty())
      return kFALSE;
    if(!HelixK0InWindow(ev, gHlxCfg.k0_mass_min, gHlxCfg.k0_mass_max))
      return kFALSE;
  }

  return kTRUE;
}

//______________________________________________________________________________
void
helix_reset_config()
{
  gHlxCfg = HelixReviewConfig{};
  std::cout << "[helix] config reset to defaults" << std::endl;
}

//______________________________________________________________________________
void
helix_status()
{
  std::cout << "\n=== Helix review config ===" << std::endl;
  std::cout << "  scan_start: " << gHlxCfg.scan_start << std::endl;
  if(gHlxCfg.mom_min >= 0.0 || gHlxCfg.mom_max >= 0.0) {
    std::cout << "  mom range: [" << gHlxCfg.mom_min << ", " << gHlxCfg.mom_max << "] GeV/c";
    if(gHlxCfg.track_for_mom >= 0)
      std::cout << "  track=" << gHlxCfg.track_for_mom;
    else
      std::cout << "  (any track)";
    std::cout << (gHlxBranches.has_mom0 ? "" : "  [branch mom0 MISSING — filter inactive/fail]")
              << std::endl;
  } else {
    std::cout << "  mom range: off" << std::endl;
  }
  if(gHlxCfg.pid_mask != 0U) {
    std::cout << "  pid filter: (pid & 0x" << std::hex << gHlxCfg.pid_mask << ") == 0x" << gHlxCfg.pid_require
              << std::dec << std::endl;
  } else {
    std::cout << "  pid filter: off" << std::endl;
  }
  std::cout << "  skip_accidental_tracks: " << (gHlxCfg.skip_accidental_tracks ? "yes" : "no")
            << (gHlxBranches.has_is_accidental ? "" : "  [is_accidental MISSING]") << std::endl;
  std::cout << "  require_beam_track: " << (gHlxCfg.require_beam_track ? "yes" : "no")
            << (gHlxBranches.has_is_beam ? "" : "  [is_beam MISSING]") << std::endl;
  std::cout << "  ntTpc range: [" << gHlxCfg.min_ntTpc << ", " << gHlxCfg.max_ntTpc << "]" << std::endl;
  if(gHlxCfg.close_dist_max >= 0.0) {
    std::cout << "  close_dist_max: " << gHlxCfg.close_dist_max << " mm"
              << (gHlxBranches.has_closeDistTpc ? "" : "  [closeDistTpc MISSING]") << std::endl;
  } else {
    std::cout << "  close_dist_max: off" << std::endl;
  }
  if(gHlxCfg.require_k0) {
    std::cout << "  require_k0: yes  mass [" << gHlxCfg.k0_mass_min << ", " << gHlxCfg.k0_mass_max << "] GeV/c^2"
              << (gHlxBranches.has_k0_mass ? "" : "  [k0_mass MISSING]") << std::endl;
  } else {
    std::cout << "  require_k0: off" << std::endl;
  }
  std::cout << "Optional branches: is_beam=" << gHlxBranches.has_is_beam
            << " is_accidental=" << gHlxBranches.has_is_accidental << " mom0=" << gHlxBranches.has_mom0
            << " dEdx=" << gHlxBranches.has_dEdx << " vtx=" << gHlxBranches.has_vtxTpc
            << " calpos=" << gHlxBranches.has_calpos << " k0_mass=" << gHlxBranches.has_k0_mass << std::endl;
  std::cout << "  display target origin (z'=z-Z_TARGET): "
            << (gHlxDisplayTargetOrigin ? "yes" : "no") << std::endl;
  std::cout << "  draw calpos overlay: " << (gHlxDrawCalpos ? "yes" : "no")
            << (gHlxBranches.has_calpos ? "" : "  [calpos MISSING]") << std::endl;
  std::cout << "  draw layer rings (bottom): " << (gHlxDrawLayerRings ? "yes" : "no") << std::endl;
  std::cout << "  prefer helix_t theta: " << (gHlxHelixPreferHitTheta ? "yes" : "no") << std::endl;
  std::cout << "===========================\n" << std::endl;
}

//______________________________________________________________________________
void
SetMomRange(Double_t min_mom, Double_t max_mom)
{
  gHlxCfg.mom_min = min_mom;
  gHlxCfg.mom_max = max_mom;
  std::cout << "[helix] SetMomRange: [" << min_mom << ", " << max_mom << "] GeV/c" << std::endl;
}

void
SetMomTrackIndex(Int_t track_index)
{
  gHlxCfg.track_for_mom = track_index;
  std::cout << "[helix] SetMomTrackIndex: " << track_index << std::endl;
}

void
SetPidFilter(UInt_t mask, UInt_t require)
{
  gHlxCfg.pid_mask = mask;
  gHlxCfg.pid_require = require;
  std::cout << "[helix] SetPidFilter: (pid & 0x" << std::hex << mask << ") == 0x" << require << std::dec
            << std::endl;
}

void
SetSkipAccidental(Bool_t on = kTRUE)
{
  gHlxCfg.skip_accidental_tracks = on;
  std::cout << "[helix] SetSkipAccidental: " << (on ? "yes" : "no") << std::endl;
}

void
SetRequireBeamTrack(Bool_t on = kTRUE)
{
  gHlxCfg.require_beam_track = on;
  std::cout << "[helix] SetRequireBeamTrack: " << (on ? "yes" : "no") << std::endl;
}

void
SetNtTpcRange(Int_t min_nt, Int_t max_nt)
{
  gHlxCfg.min_ntTpc = min_nt;
  gHlxCfg.max_ntTpc = max_nt;
  std::cout << "[helix] SetNtTpcRange: [" << min_nt << ", " << max_nt << "]" << std::endl;
}

void
SetCloseDistMax(Double_t mm)
{
  gHlxCfg.close_dist_max = mm;
  std::cout << "[helix] SetCloseDistMax: " << mm << " mm" << std::endl;
}

void
SetRequireK0(Bool_t on = kTRUE)
{
  gHlxCfg.require_k0 = on;
  std::cout << "[helix] SetRequireK0: " << (on ? "yes" : "no") << std::endl;
}

void
SetK0MassWindow(Double_t lo, Double_t hi)
{
  gHlxCfg.k0_mass_min = lo;
  gHlxCfg.k0_mass_max = hi;
  std::cout << "[helix] SetK0MassWindow: [" << lo << ", " << hi << "] GeV/c^2" << std::endl;
}

void
SetScanStart(Long64_t entry)
{
  gHlxCfg.scan_start = entry;
  std::cout << "[helix] SetScanStart: " << entry << std::endl;
}

void
SetHelixViewMargin(Double_t margin_frac = 0.22)
{
  gHlxViewMargin = TMath::Max(0.0, margin_frac);
  std::cout << "[helix] SetHelixViewMargin: " << gHlxViewMargin << std::endl;
}

void
SetHelixPreferHitTheta(Bool_t on = kTRUE)
{
  gHlxHelixPreferHitTheta = on;
  std::cout << "[helix] SetHelixPreferHitTheta: " << (on ? "yes (helix_t first)" : "no (theta_min/max first)")
            << std::endl;
}

void
SetHelixDisplayTargetOrigin(Bool_t on = kTRUE)
{
  gHlxDisplayTargetOrigin = on;
  std::cout << "[helix] SetHelixDisplayTargetOrigin: " << (on ? "yes (ROOT X = z_tpc - Z_TARGET)" : "no (ROOT X = z_tpc)")
            << std::endl;
}

void
SetHelixDrawCalpos(Bool_t on = kTRUE)
{
  gHlxDrawCalpos = on;
  std::cout << "[helix] SetHelixDrawCalpos: " << (on ? "yes" : "no") << std::endl;
}

void
SetHelixDrawLayerRings(Bool_t on = kTRUE)
{
  gHlxDrawLayerRings = on;
  std::cout << "[helix] SetHelixDrawLayerRings: " << (on ? "yes (bottom y=-310)" : "no") << std::endl;
}

void
helix_coords_help()
{
  std::cout << "\n=== TPC Helix 3D: coordinate & axis conventions ===\n"
            << "\nData in DST (what we plot after conversion):\n"
            << "  [B] TPC lab (x,y,z)  hitpos_*, calpos_*, vtxTpc  [mm]\n"
            << "      y = drift,  z = beam axis,  z=0 = TPC mechanical center\n"
            << "      target at (0, 0, Z_TARGET) with Z_TARGET = -143 mm\n"
            << "  helix_cx,cy,... are [A] internal fit coords; LocalToGlobal before draw:\n"
            << "      x = -xi,  y = zi,  z = yi + Z_TARGET\n"
            << "\nNOT in this file: experiment GEANT global (DCGeomMan etc.)\n"
            << "\nROOT SetPoint(i, X, Y, Z):\n"
            << "      X = z_tpc   (optional: z_tpc - Z_TARGET via SetHelixDisplayTargetOrigin)\n"
            << "      Y = x_tpc\n"
            << "      Z = y_tpc   (drift)\n"
            << "\nCorner arrows: red=TPC z, green=TPC x, blue=TPC y\n"
            << "Frame: octagon at (x,z)=(0,0); target cylinder at z=Z_TARGET.\n"
            << "====================================================\n" << std::endl;
}

//______________________________________________________________________________
inline Double_t
HelixPlotLongZ(Double_t z_tpc)
{
  return gHlxDisplayTargetOrigin ? (z_tpc - kTPC_Z_TARGET) : z_tpc;
}

inline void
HelixInternalToTPCGlobal(Double_t xi, Double_t yi, Double_t zi, Double_t& xg, Double_t& yg, Double_t& zg)
{
  // TPCLocalTrackHelix::LocalToGlobal / GetLocalCalPosHelix と同じ
  xg = -xi;
  yg = zi;
  zg = yi + kTPC_Z_TARGET;
}

inline void
TPCGlobalToDisplay(Double_t xg, Double_t yg, Double_t zg, Double_t& dispX, Double_t& dispY, Double_t& dispZ)
{
  dispX = HelixPlotLongZ(zg);
  dispY = xg;
  dispZ = yg;
}

//______________________________________________________________________________
void
DrawAxisLabel3D(Double_t x, Double_t y, Double_t z, const char* label, Color_t color, TView* view)
{
  if(!view || !gHlxCanvas)
    return;
  gHlxCanvas->cd();
  Double_t wc[3] = {x, y, z};
  Double_t ndc[3];
  view->WCtoNDC(wc, ndc);
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextColor(color);
  latex->SetTextSize(0.022);
  latex->SetTextAlign(22);
  latex->DrawLatex(ndc[0], ndc[1], label);
  gHlxOverlay.push_back(latex);
}

inline void
HelixPointToDisplay(Double_t cx, Double_t cy, Double_t z0, Double_t r, Double_t dz, Double_t t,
                    Double_t& dispX, Double_t& dispY, Double_t& dispZ)
{
  const Double_t xi = cx + r * TMath::Cos(t);
  const Double_t yi = cy + r * TMath::Sin(t);
  const Double_t zi = z0 + dz * r * t;
  Double_t xg, yg, zg;
  HelixInternalToTPCGlobal(xi, yi, zi, xg, yg, zg);
  TPCGlobalToDisplay(xg, yg, zg, dispX, dispY, dispZ);
}

inline void
LocalHitToDisplay(Double_t xloc, Double_t yloc, Double_t zloc, Double_t& dispX, Double_t& dispY,
                  Double_t& dispZ)
{
  TPCGlobalToDisplay(xloc, yloc, zloc, dispX, dispY, dispZ);
}

static Double_t
HelixMinDistToTargetMm(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz, Double_t tmin,
                       Double_t tmax, Int_t nSample = 64)
{
  if(nSample < 8)
    nSample = 8;
  const Double_t tgtX = 0.0;
  const Double_t tgtY = 0.0;
  const Double_t tgtZ = kTPC_Z_TARGET;
  Double_t best = std::numeric_limits<Double_t>::max();
  for(Int_t i = 0; i <= nSample; ++i) {
    const Double_t t = tmin + (tmax - tmin) * static_cast<Double_t>(i) / static_cast<Double_t>(nSample);
    const Double_t xi = cx + rr * TMath::Cos(t);
    const Double_t yi = cy + rr * TMath::Sin(t);
    const Double_t zi = z0 + dz * rr * t;
    Double_t xg, yg, zg;
    HelixInternalToTPCGlobal(xi, yi, zi, xg, yg, zg);
    const Double_t d = TMath::Sqrt(TMath::Power(xg - tgtX, 2) + TMath::Power(yg - tgtY, 2) + TMath::Power(zg - tgtZ, 2));
    if(d < best)
      best = d;
  }
  return best;
}

inline Int_t
PidLineColor(UInt_t pid)
{
  if(pid & 0x2U)
    return kBlue + 1;
  if(pid & 0x4U)
    return kRed + 1;
  if(pid & 0x1U)
    return kGreen + 2;
  return kGray + 2;
}

static Bool_t
HelixThetaDrawRange(const std::vector<Double_t>& tv, Double_t& t_lo, Double_t& t_hi)
{
  std::vector<Double_t> finite;
  finite.reserve(tv.size());
  Double_t gmin = std::numeric_limits<Double_t>::max();
  Double_t gmax = std::numeric_limits<Double_t>::lowest();
  Double_t t_first = 0.0, t_last = 0.0;
  Bool_t haveFirst = kFALSE;
  for(Double_t w : tv) {
    if(!TMath::Finite(w))
      continue;
    finite.push_back(w);
    if(!haveFirst) {
      t_first = w;
      haveFirst = kTRUE;
    }
    t_last = w;
    if(w < gmin)
      gmin = w;
    if(w > gmax)
      gmax = w;
  }
  if(!haveFirst || finite.empty())
    return kFALSE;
  const Double_t gspan = gmax - gmin;
  if(gspan < 1e-12)
    return kFALSE;

  const Double_t e_lo = TMath::Min(t_first, t_last);
  const Double_t e_hi = TMath::Max(t_first, t_last);
  const Double_t espan = e_hi - e_lo;

  const Bool_t preferEndsWide =
    (gspan > 2.0 * TMath::TwoPi() && espan + 1e-9 < gspan * 0.75 && espan > 1e-9);
  const Bool_t preferEndsOutlier =
    (espan > 1e-9 && espan + 1e-9 < gspan * 0.9 && gspan > espan * 2.5 + 0.3);

  if(preferEndsWide || preferEndsOutlier) {
    t_lo = e_lo;
    t_hi = e_hi;
  } else {
    t_lo = gmin;
    t_hi = gmax;
  }

  if(finite.size() >= 5U) {
    std::sort(finite.begin(), finite.end());
    const std::size_t n = finite.size();
    const std::size_t k0 = static_cast<std::size_t>(0.06 * static_cast<Double_t>(n - 1));
    const std::size_t k1 = static_cast<std::size_t>(0.94 * static_cast<Double_t>(n - 1));
    const Double_t p_lo = finite[k0];
    const Double_t p_hi = finite[k1];
    const Double_t pspan = p_hi - p_lo;
    const Double_t cspan = t_hi - t_lo;
    if(pspan > 1e-9 && cspan > pspan * 1.3 + 0.2) {
      t_lo = p_lo;
      t_hi = p_hi;
    }
  }

  const Double_t mrg = (t_hi - t_lo) * 0.02 + 5e-4;
  t_lo -= mrg;
  t_hi += mrg;

  if(t_hi - t_lo < 1e-12)
    return kFALSE;
  return kTRUE;
}

//______________________________________________________________________________
static Bool_t
HelixDrawThetaRange(const HelixEventData& ev, Int_t itrack, Double_t& t_lo, Double_t& t_hi)
{
  t_lo = 0.0;
  t_hi = 0.0;
  Bool_t fromHits = kFALSE;
  if(gHlxHelixPreferHitTheta && itrack < static_cast<Int_t>(ev.helix_t.size()) &&
     !ev.helix_t[itrack].empty()) {
    fromHits = HelixThetaDrawRange(ev.helix_t[itrack], t_lo, t_hi);
  }

  if(!fromHits && gHlxHelixTRangeMode != 0 && itrack < static_cast<Int_t>(ev.helix_t_min.size()) &&
     itrack < static_cast<Int_t>(ev.helix_t_max.size())) {
    const Double_t a = ev.helix_t_min[itrack];
    const Double_t b = ev.helix_t_max[itrack];
    if(TMath::Finite(a) && TMath::Finite(b) && TMath::Abs(b - a) > 1e-12) {
      const Double_t lo = TMath::Min(a, b);
      const Double_t hi = TMath::Max(a, b);
      const Double_t mrg = (hi - lo) * 0.02 + 5e-4;
      t_lo = lo - mrg;
      t_hi = hi + mrg;
      fromHits = kTRUE;
    }
  }

  if(fromHits && !gHlxHelixPreferHitTheta && gHlxHelixTRangeMode != 0 &&
     itrack < static_cast<Int_t>(ev.helix_t_min.size()) &&
     itrack < static_cast<Int_t>(ev.helix_t_max.size())) {
    const Double_t a = ev.helix_t_min[itrack];
    const Double_t b = ev.helix_t_max[itrack];
    if(TMath::Finite(a) && TMath::Finite(b) && TMath::Abs(b - a) > 1e-12) {
      const Double_t lo = TMath::Min(a, b);
      const Double_t hi = TMath::Max(a, b);
      const Double_t mrg = (hi - lo) * 0.02 + 5e-4;
      t_lo = TMath::Max(t_lo, lo - mrg);
      t_hi = TMath::Min(t_hi, hi + mrg);
      if(t_hi - t_lo < 1e-12)
        return kFALSE;
    }
  }

  if(!fromHits && itrack < static_cast<Int_t>(ev.helix_t.size()) && !ev.helix_t[itrack].empty()) {
    fromHits = HelixThetaDrawRange(ev.helix_t[itrack], t_lo, t_hi);
  }

  return fromHits;
}

//______________________________________________________________________________
static void
PadViewRange(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax,
             Double_t& oxmin, Double_t& oxmax, Double_t& oymin, Double_t& oymax, Double_t& ozmin,
             Double_t& ozmax)
{
  const Double_t m = gHlxViewMargin;
  const Double_t dx = (xmax - xmin) * m;
  const Double_t dy = (ymax - ymin) * m;
  const Double_t dz = (zmax - zmin) * m;
  oxmin = xmin - dx;
  oxmax = xmax + dx;
  oymin = ymin - dy;
  oymax = ymax + dy;
  ozmin = zmin - dz;
  ozmax = zmax + dz;
}

//______________________________________________________________________________
static void
DrawHelixPolylineSegments(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz, Double_t tmin,
                          Double_t tmax, Int_t col, Int_t nSample)
{
  if(nSample < 2)
    nSample = 2;
  const Double_t tspan = tmax - tmin;
  nSample = static_cast<Int_t>(TMath::Min(640.0, TMath::Max(96.0, static_cast<Double_t>(nSample))));
  const Double_t helixDen = static_cast<Double_t>(nSample - 1);

  std::vector<Double_t> xs(nSample);
  std::vector<Double_t> ys(nSample);
  std::vector<Double_t> zs(nSample);
  for(Int_t s = 0; s < nSample; s++) {
    const Double_t t = tmin + tspan * static_cast<Double_t>(s) / helixDen;
    HelixPointToDisplay(cx, cy, z0, rr, dz, t, xs[s], ys[s], zs[s]);
  }

  // ROOT 3D では多点 TPolyLine3D が点列に見えることがある → 2 点セグメント連結
  for(Int_t s = 0; s < nSample - 1; s++) {
    if(!TMath::Finite(xs[s]) || !TMath::Finite(xs[s + 1]))
      continue;
    TPolyLine3D* seg = new TPolyLine3D(2);
    seg->SetLineColor(col);
    seg->SetLineWidth(3);
    seg->SetLineStyle(kSolid);
    seg->SetPoint(0, xs[s], ys[s], zs[s]);
    seg->SetPoint(1, xs[s + 1], ys[s + 1], zs[s + 1]);
    seg->Draw("same");
    gHlxHelixLines.push_back(seg);
  }
}

//______________________________________________________________________________
void
clear_helix_objects()
{
  // Draw 済みオブジェクトはパッドが所有。手動 delete + Clear の二重解放を避ける。
  gHlxHelixLines.clear();
  gHlxHitMarkers.clear();
  gHlxFrame.clear();
  gHlxOverlay.clear();
  gHlxView = nullptr;
}

//______________________________________________________________________________
void
draw_tpc_frame_and_corner_axes()
{
  if(!gHlxCanvas)
    return;
  gHlxCanvas->cd();

  // check_tpc_track_3d.C と同じ固定レンジ（八角形半径・フレーム高さ）
  const Double_t xmin_orig = -311.0, xmax_orig = 311.0;
  const Double_t ymin_orig = -310.0, ymax_orig = 310.0;
  const Double_t zmin_orig = -311.0, zmax_orig = 311.0;

  const Double_t xmin = HelixPlotLongZ(zmin_orig), xmax = HelixPlotLongZ(zmax_orig);
  const Double_t ymin = xmin_orig, ymax = xmax_orig;
  const Double_t zmin = ymin_orig, zmax = ymax_orig;

  Double_t vxmin, vxmax, vymin, vymax, vzmin, vzmax;
  PadViewRange(xmin, xmax, ymin, ymax, zmin, zmax, vxmin, vxmax, vymin, vymax, vzmin, vzmax);

  gHlxView = TView::CreateView(1, 0, 0);
  gHlxView->SetRange(vxmin, vymin, vzmin, vxmax, vymax, vzmax);

  Double_t axisOriginX = xmin + (xmax - xmin) * 0.05;
  Double_t axisOriginY = ymin + (ymax - ymin) * 0.05;
  Double_t axisOriginZ = zmin + (zmax - zmin) * 0.05;
  Double_t axisLength = TMath::Min(TMath::Min(xmax - xmin, ymax - ymin), zmax - zmin) * 0.1;

  auto push_line = [](TPolyLine3D* pl) {
    pl->Draw();
    gHlxFrame.push_back(pl);
  };

  TPolyLine3D* axisX = new TPolyLine3D(2);
  axisX->SetLineColor(kRed);
  axisX->SetLineWidth(2);
  axisX->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisX->SetPoint(1, axisOriginX + axisLength, axisOriginY, axisOriginZ);
  push_line(axisX);

  Double_t arrowSize = axisLength * 0.15;
  TPolyLine3D* arrowX1 = new TPolyLine3D(3);
  arrowX1->SetLineColor(kRed);
  arrowX1->SetLineWidth(2);
  arrowX1->SetPoint(0, axisOriginX + axisLength, axisOriginY, axisOriginZ);
  arrowX1->SetPoint(1, axisOriginX + axisLength - arrowSize, axisOriginY - arrowSize * 0.3, axisOriginZ);
  arrowX1->SetPoint(2, axisOriginX + axisLength - arrowSize, axisOriginY + arrowSize * 0.3, axisOriginZ);
  push_line(arrowX1);

  TPolyLine3D* axisY = new TPolyLine3D(2);
  axisY->SetLineColor(kGreen + 2);
  axisY->SetLineWidth(2);
  axisY->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisY->SetPoint(1, axisOriginX, axisOriginY + axisLength, axisOriginZ);
  push_line(axisY);

  TPolyLine3D* arrowY1 = new TPolyLine3D(3);
  arrowY1->SetLineColor(kGreen + 2);
  arrowY1->SetLineWidth(2);
  arrowY1->SetPoint(0, axisOriginX, axisOriginY + axisLength, axisOriginZ);
  arrowY1->SetPoint(1, axisOriginX - arrowSize * 0.3, axisOriginY + axisLength - arrowSize, axisOriginZ);
  arrowY1->SetPoint(2, axisOriginX + arrowSize * 0.3, axisOriginY + axisLength - arrowSize, axisOriginZ);
  push_line(arrowY1);

  TPolyLine3D* axisZ = new TPolyLine3D(2);
  axisZ->SetLineColor(kBlue);
  axisZ->SetLineWidth(2);
  axisZ->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisZ->SetPoint(1, axisOriginX, axisOriginY, axisOriginZ + axisLength);
  push_line(axisZ);

  TPolyLine3D* arrowZ1 = new TPolyLine3D(3);
  arrowZ1->SetLineColor(kBlue);
  arrowZ1->SetLineWidth(2);
  arrowZ1->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ + axisLength);
  arrowZ1->SetPoint(1, axisOriginX - arrowSize * 0.3, axisOriginY, axisOriginZ + axisLength - arrowSize);
  arrowZ1->SetPoint(2, axisOriginX + arrowSize * 0.3, axisOriginY, axisOriginZ + axisLength - arrowSize);
  push_line(arrowZ1);

  const Double_t frameHeight = 620.0;
  const Double_t frameHalfHeight = frameHeight / 2.0;
  const Double_t frameFaceToFace = 622.0;
  const Double_t frameRadius = frameFaceToFace / 2.0;
  const Int_t nOctagonSides = 8;
  const Double_t octagonVertexRadius = frameRadius / TMath::Cos(TMath::Pi() / 8.0);

  std::vector<Double_t> octagon_x_orig(nOctagonSides);
  std::vector<Double_t> octagon_z_orig(nOctagonSides);
  for(Int_t i = 0; i < nOctagonSides; i++) {
    Double_t angle = i * TMath::Pi() / 4.0 + TMath::Pi() / 8.0;
    octagon_x_orig[i] = octagonVertexRadius * TMath::Cos(angle);
    octagon_z_orig[i] = octagonVertexRadius * TMath::Sin(angle);
  }

  for(Int_t y_sign = -1; y_sign <= 1; y_sign += 2) {
    Double_t y_orig = y_sign * frameHalfHeight;
    TPolyLine3D* octagon = new TPolyLine3D(nOctagonSides + 1);
    octagon->SetLineColor(kGray + 1);
    octagon->SetLineWidth(2);
    for(Int_t i = 0; i <= nOctagonSides; i++) {
      Int_t idx = i % nOctagonSides;
      octagon->SetPoint(i, HelixPlotLongZ(octagon_z_orig[idx]), octagon_x_orig[idx], y_orig);
    }
    push_line(octagon);
  }

  for(Int_t i = 0; i < nOctagonSides; i++) {
    TPolyLine3D* edge = new TPolyLine3D(2);
    edge->SetLineColor(kGray + 1);
    edge->SetLineWidth(2);
    edge->SetPoint(0, HelixPlotLongZ(octagon_z_orig[i]), octagon_x_orig[i], -frameHalfHeight);
    edge->SetPoint(1, HelixPlotLongZ(octagon_z_orig[i]), octagon_x_orig[i], frameHalfHeight);
    push_line(edge);
  }

  if(gHlxDrawLayerRings) {
    const Int_t nLayerSeg = 72;
    const Double_t yBot = -frameHalfHeight;
    for(Int_t layer = 0; layer < kHelixNumTpcLayers; ++layer) {
      const Double_t rPad = kHelixTpcLayerRadius[layer];
      TPolyLine3D* ring = new TPolyLine3D(nLayerSeg + 1);
      ring->SetLineColor(kTeal - 4);
      ring->SetLineStyle(3);
      ring->SetLineWidth(1);
      for(Int_t i = 0; i <= nLayerSeg; ++i) {
        const Double_t theta = 2.0 * TMath::Pi() * i / nLayerSeg;
        const Double_t x_orig = rPad * TMath::Sin(theta);
        const Double_t z_orig = rPad * TMath::Cos(theta) + kTPC_Z_TARGET;
        ring->SetPoint(i, HelixPlotLongZ(z_orig), x_orig, yBot);
      }
      push_line(ring);
    }
  }

  const Double_t targetZ_orig = kTPC_Z_TARGET;
  const Double_t targetHeight = 100.0;
  const Double_t targetHalfHeight = targetHeight / 2.0;
  const Double_t targetRadius = 40.0;
  const Int_t nCirclePoints = 64;
  const Int_t nVerticalLines = 8;

  for(Int_t y_sign = -1; y_sign <= 1; y_sign += 2) {
    Double_t y_orig = y_sign * targetHalfHeight;
    TPolyLine3D* circle = new TPolyLine3D(nCirclePoints + 1);
    circle->SetLineColor(kBlack);
    circle->SetLineWidth(2);
    for(Int_t i = 0; i <= nCirclePoints; i++) {
      Double_t theta = 2.0 * TMath::Pi() * i / nCirclePoints;
      Double_t x_orig = targetRadius * TMath::Cos(theta);
      Double_t z_orig = targetRadius * TMath::Sin(theta) + targetZ_orig;
      circle->SetPoint(i, HelixPlotLongZ(z_orig), x_orig, y_orig);
    }
    push_line(circle);
  }

  for(Int_t i = 0; i < nVerticalLines; i++) {
    Double_t theta = 2.0 * TMath::Pi() * i / nVerticalLines;
    Double_t x_orig = targetRadius * TMath::Cos(theta);
    Double_t z_orig = targetRadius * TMath::Sin(theta) + targetZ_orig;
    TPolyLine3D* line = new TPolyLine3D(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetPoint(0, HelixPlotLongZ(z_orig), x_orig, -targetHalfHeight);
    line->SetPoint(1, HelixPlotLongZ(z_orig), x_orig, targetHalfHeight);
    push_line(line);
  }

  const Double_t loff = axisLength * 0.25;
  const Double_t lxlab = axisOriginX + axisLength + loff;
  DrawAxisLabel3D(lxlab, axisOriginY, axisOriginZ, "z", kRed, gHlxView);
  DrawAxisLabel3D(axisOriginX, axisOriginY + axisLength + loff, axisOriginZ, "x", kGreen + 2, gHlxView);
  DrawAxisLabel3D(axisOriginX, axisOriginY, axisOriginZ + axisLength + loff, "y", kBlue, gHlxView);
}

//______________________________________________________________________________
static void
helix_probe_branches(TTree* tree)
{
  gHlxBranches = HelixBranchFlags{};
  if(!tree)
    return;
  gHlxBranches.has_is_beam = (tree->GetBranch("is_beam") != nullptr);
  gHlxBranches.has_is_accidental = (tree->GetBranch("is_accidental") != nullptr);
  gHlxBranches.has_charge = (tree->GetBranch("charge") != nullptr);
  gHlxBranches.has_mom0 = (tree->GetBranch("mom0") != nullptr);
  gHlxBranches.has_dEdx = (tree->GetBranch("dEdx") != nullptr);
  gHlxBranches.has_beam_flag = (tree->GetBranch("beam_flag") != nullptr);
  gHlxBranches.has_effective_ntTpc = (tree->GetBranch("effective_ntTpc") != nullptr);
  gHlxBranches.has_closeDistTpc = (tree->GetBranch("closeDistTpc") != nullptr);
  gHlxBranches.has_vtxTpc = (tree->GetBranch("vtxTpc") != nullptr);
  gHlxBranches.has_calpos =
    (tree->GetBranch("calpos_x") != nullptr) && (tree->GetBranch("calpos_y") != nullptr) &&
    (tree->GetBranch("calpos_z") != nullptr);
  gHlxBranches.has_k0_mass = (tree->GetBranch("k0_mass") != nullptr);
  gHlxBranches.has_lambda_mass = (tree->GetBranch("lambda_mass") != nullptr);
}

//______________________________________________________________________________
void
helix_set_path(const char* path)
{
  if(gHlxFile) {
    gHlxFile->Close();
    delete gHlxFile;
    gHlxFile = nullptr;
  }
  if(gHlxReader) {
    delete gHlxReader;
    gHlxReader = nullptr;
  }

  gHlxFile = TFile::Open(path);
  if(!gHlxFile || gHlxFile->IsZombie()) {
    std::cerr << "Error: Cannot open file " << path << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(gHlxFile->Get("tpc"));
  if(!tree) {
    std::cerr << "Error: Cannot find tree 'tpc' in file" << std::endl;
    return;
  }

  helix_probe_branches(tree);

  if((tree->GetBranch("helix_theta_min") != nullptr) && (tree->GetBranch("helix_theta_max") != nullptr)) {
    gHlxHelixTRangeMode = 1;
  } else if((tree->GetBranch("helix_t_min") != nullptr) && (tree->GetBranch("helix_t_max") != nullptr)) {
    gHlxHelixTRangeMode = 2;
  } else {
    gHlxHelixTRangeMode = 0;
  }

  gHlxReader = new TTreeReader(tree);
  gHlxCfg.scan_start = 0;
  gHlxCurrentEntry = -1;

  std::cout << "File opened: " << path << std::endl;
  std::cout << "Total entries: " << tree->GetEntries() << std::endl;
  if(gHlxHelixTRangeMode == 1) {
    std::cout << "  helix_theta_min / helix_theta_max: present" << std::endl;
  } else if(gHlxHelixTRangeMode == 2) {
    std::cout << "  helix_t_min / helix_t_max: present (legacy branch names)" << std::endl;
  } else {
    std::cout << "  track theta range branches: absent (use per-hit helix_t)" << std::endl;
  }
  std::cout << "Optional branches: is_beam=" << gHlxBranches.has_is_beam
            << " mom0=" << gHlxBranches.has_mom0 << " closeDist=" << gHlxBranches.has_closeDistTpc
            << " vtx=" << gHlxBranches.has_vtxTpc << " calpos=" << gHlxBranches.has_calpos
            << " k0_mass=" << gHlxBranches.has_k0_mass << std::endl;

  if(!gHlxRandom) {
    gHlxRandom = new TRandom3();
    gHlxRandom->SetSeed(0);
  }
}

//______________________________________________________________________________
Bool_t
load_helix_event(Long64_t entry = -1)
{
  if(!gHlxReader) {
    std::cerr << "Error: No file opened. Use helix_set_path() first." << std::endl;
    return kFALSE;
  }

  Long64_t nentries = gHlxReader->GetEntries(false);
  if(nentries == 0) {
    std::cerr << "Error: No entries in tree" << std::endl;
    return kFALSE;
  }

  if(entry < 0) {
    if(!gHlxRandom) {
      gHlxRandom = new TRandom3();
      gHlxRandom->SetSeed(0);
    }
    entry = gHlxRandom->Integer(nentries);
  }

  if(entry >= nentries) {
    std::cerr << "Error: Entry " << entry << " is out of range [0, " << nentries << ")" << std::endl;
    return kFALSE;
  }

  gHlxReader->Restart();

  TTreeReaderValue<UInt_t> runnum(*gHlxReader, "run_number");
  TTreeReaderValue<UInt_t> evnum(*gHlxReader, "event_number");
  TTreeReaderValue<Int_t> ntTpc(*gHlxReader, "ntTpc");
  TTreeReaderValue<std::vector<Int_t>> nhtrack(*gHlxReader, "nhtrack");
  TTreeReaderValue<std::vector<Int_t>> pid(*gHlxReader, "pid");
  TTreeReaderValue<std::vector<Double_t>> helix_cx(*gHlxReader, "helix_cx");
  TTreeReaderValue<std::vector<Double_t>> helix_cy(*gHlxReader, "helix_cy");
  TTreeReaderValue<std::vector<Double_t>> helix_z0(*gHlxReader, "helix_z0");
  TTreeReaderValue<std::vector<Double_t>> helix_r(*gHlxReader, "helix_r");
  TTreeReaderValue<std::vector<Double_t>> helix_dz(*gHlxReader, "helix_dz");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> helix_t(*gHlxReader, "helix_t");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(*gHlxReader, "hitpos_x");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_y(*gHlxReader, "hitpos_y");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(*gHlxReader, "hitpos_z");

  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pCalposX;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pCalposY;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pCalposZ;
  if(gHlxBranches.has_calpos) {
    pCalposX.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "calpos_x"));
    pCalposY.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "calpos_y"));
    pCalposZ.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "calpos_z"));
  }

  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pHelixTmin;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pHelixTmax;
  if(gHlxHelixTRangeMode == 1) {
    pHelixTmin.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_theta_min"));
    pHelixTmax.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_theta_max"));
  } else if(gHlxHelixTRangeMode == 2) {
    pHelixTmin.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_t_min"));
    pHelixTmax.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_t_max"));
  }

  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pIsBeam;
  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pIsAccidental;
  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pCharge;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pMom0;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pDEdx;
  std::unique_ptr<TTreeReaderValue<Int_t>> pBeamFlag;
  std::unique_ptr<TTreeReaderValue<Int_t>> pEffectiveNt;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pCloseDist;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pVtxX;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pVtxY;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pVtxZ;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pK0Mass;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pLambdaMass;

  if(gHlxBranches.has_is_beam)
    pIsBeam.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "is_beam"));
  if(gHlxBranches.has_is_accidental)
    pIsAccidental.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "is_accidental"));
  if(gHlxBranches.has_charge)
    pCharge.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "charge"));
  if(gHlxBranches.has_mom0)
    pMom0.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "mom0"));
  if(gHlxBranches.has_dEdx)
    pDEdx.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "dEdx"));
  if(gHlxBranches.has_beam_flag)
    pBeamFlag.reset(new TTreeReaderValue<Int_t>(*gHlxReader, "beam_flag"));
  if(gHlxBranches.has_effective_ntTpc)
    pEffectiveNt.reset(new TTreeReaderValue<Int_t>(*gHlxReader, "effective_ntTpc"));
  if(gHlxBranches.has_closeDistTpc)
    pCloseDist.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "closeDistTpc"));
  if(gHlxBranches.has_vtxTpc) {
    pVtxX.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "vtxTpc"));
    pVtxY.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "vtyTpc"));
    pVtxZ.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "vtzTpc"));
  }
  if(gHlxBranches.has_k0_mass)
    pK0Mass.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "k0_mass"));
  if(gHlxBranches.has_lambda_mass)
    pLambdaMass.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "lambda_mass"));

  gHlxReader->SetEntry(entry);

  gHlxEvent = HelixEventData{};
  gHlxEvent.runnum = *runnum;
  gHlxEvent.evnum = *evnum;
  gHlxEvent.ntTpc = *ntTpc;
  gHlxEvent.nhtrack = *nhtrack;
  gHlxEvent.pid = *pid;
  gHlxEvent.helix_cx = *helix_cx;
  gHlxEvent.helix_cy = *helix_cy;
  gHlxEvent.helix_z0 = *helix_z0;
  gHlxEvent.helix_r = *helix_r;
  gHlxEvent.helix_dz = *helix_dz;
  gHlxEvent.helix_t = *helix_t;
  gHlxEvent.hitpos_x = *hitpos_x;
  gHlxEvent.hitpos_y = *hitpos_y;
  gHlxEvent.hitpos_z = *hitpos_z;
  if(pCalposX && pCalposY && pCalposZ) {
    gHlxEvent.calpos_x = **pCalposX;
    gHlxEvent.calpos_y = **pCalposY;
    gHlxEvent.calpos_z = **pCalposZ;
  }

  if(pHelixTmin && pHelixTmax) {
    gHlxEvent.helix_t_min = **pHelixTmin;
    gHlxEvent.helix_t_max = **pHelixTmax;
  }
  if(pIsBeam)
    gHlxEvent.is_beam = **pIsBeam;
  if(pIsAccidental)
    gHlxEvent.is_accidental = **pIsAccidental;
  if(pCharge)
    gHlxEvent.charge = **pCharge;
  if(pMom0)
    gHlxEvent.mom0 = **pMom0;
  if(pDEdx)
    gHlxEvent.dEdx = **pDEdx;
  if(pBeamFlag)
    gHlxEvent.beam_flag = **pBeamFlag;
  if(pEffectiveNt)
    gHlxEvent.effective_ntTpc = **pEffectiveNt;
  if(pCloseDist)
    gHlxEvent.closeDistTpc = **pCloseDist;
  if(pVtxX && pVtxY && pVtxZ) {
    gHlxEvent.vtxTpc = **pVtxX;
    gHlxEvent.vtyTpc = **pVtxY;
    gHlxEvent.vtzTpc = **pVtxZ;
  }
  if(pK0Mass)
    gHlxEvent.k0_mass = **pK0Mass;
  if(pLambdaMass)
    gHlxEvent.lambda_mass = **pLambdaMass;

  gHlxCurrentEntry = entry;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t
helix_load_filter_snapshot(Long64_t entry, HelixEventData& out)
{
  if(!gHlxReader)
    return kFALSE;
  Long64_t nentries = gHlxReader->GetEntries(false);
  if(entry < 0 || entry >= nentries)
    return kFALSE;

  gHlxReader->Restart();
  TTreeReaderValue<Int_t> ntTpc(*gHlxReader, "ntTpc");

  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pIsBeam;
  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pIsAccidental;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pMom0;
  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> pPid;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> pCloseDist;
  std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> pK0Mass;

  if(gHlxBranches.has_is_beam)
    pIsBeam.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "is_beam"));
  if(gHlxBranches.has_is_accidental)
    pIsAccidental.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "is_accidental"));
  if(gHlxBranches.has_mom0)
    pMom0.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "mom0"));
  pPid.reset(new TTreeReaderValue<std::vector<Int_t>>(*gHlxReader, "pid"));
  if(gHlxBranches.has_closeDistTpc)
    pCloseDist.reset(new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*gHlxReader, "closeDistTpc"));
  if(gHlxBranches.has_k0_mass)
    pK0Mass.reset(new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "k0_mass"));

  gHlxReader->SetEntry(entry);

  out = HelixEventData{};
  out.ntTpc = *ntTpc;
  if(pPid)
    out.pid = **pPid;
  if(pIsBeam)
    out.is_beam = **pIsBeam;
  if(pIsAccidental)
    out.is_accidental = **pIsAccidental;
  if(pMom0)
    out.mom0 = **pMom0;
  if(pCloseDist)
    out.closeDistTpc = **pCloseDist;
  if(pK0Mass)
    out.k0_mass = **pK0Mass;
  return kTRUE;
}

//______________________________________________________________________________
void
draw_helix_vertex_marker()
{
  if(gHlxEvent.ntTpc != 2 || !gHlxBranches.has_vtxTpc)
    return;
  if(gHlxEvent.vtxTpc.size() < 1U || gHlxEvent.vtyTpc.size() < 1U || gHlxEvent.vtzTpc.size() < 1U)
    return;
  if(gHlxEvent.vtxTpc[0].size() < 2U)
    return;

  const Double_t vx = gHlxEvent.vtxTpc[0][1];
  const Double_t vy = gHlxEvent.vtyTpc[0][1];
  const Double_t vz = gHlxEvent.vtzTpc[0][1];
  if(!TMath::Finite(vx) || !TMath::Finite(vy) || !TMath::Finite(vz))
    return;

  Double_t dx, dy, dz;
  LocalHitToDisplay(vx, vy, vz, dx, dy, dz);
  TPolyMarker3D* vtx = new TPolyMarker3D(1);
  vtx->SetMarkerStyle(29);
  vtx->SetMarkerSize(1.4);
  vtx->SetMarkerColor(kMagenta + 1);
  vtx->SetPoint(0, dx, dy, dz);
  vtx->Draw("same");
  gHlxHitMarkers.push_back(vtx);
}

//______________________________________________________________________________
void
draw_helix_info_overlay()
{
  if(!gHlxCanvas)
    return;

  const Int_t nLines =
    2 + gHlxEvent.ntTpc + (gHlxBranches.has_closeDistTpc && gHlxEvent.ntTpc == 2 ? 1 : 0) +
    (gHlxBranches.has_k0_mass && !gHlxEvent.k0_mass.empty() ? 1 : 0);
  const Double_t boxH = TMath::Min(0.38, 0.06 + 0.032 * static_cast<Double_t>(nLines));
  TPaveText* pt = new TPaveText(0.02, 0.02, 0.44, 0.02 + boxH, "NDC");
  pt->SetFillColor(kWhite);
  pt->SetFillStyle(4000 + 850); // 85% 不透明
  pt->SetBorderSize(1);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.022);
  pt->SetMargin(0.04);

  std::ostringstream hdr;
  hdr << "Run " << gHlxEvent.runnum << "  Ev " << gHlxEvent.evnum << "  entry " << gHlxCurrentEntry;
  pt->AddText(hdr.str().c_str());

  std::ostringstream ntline;
  ntline << "ntTpc=" << gHlxEvent.ntTpc;
  if(gHlxBranches.has_effective_ntTpc)
    ntline << "  effective_ntTpc=" << gHlxEvent.effective_ntTpc;
  if(gHlxBranches.has_beam_flag)
    ntline << "  beam_flag=" << gHlxEvent.beam_flag;
  pt->AddText(ntline.str().c_str());

  for(Int_t it = 0; it < gHlxEvent.ntTpc; ++it) {
    std::ostringstream tl;
    tl << HelixTrackSubscript(it) << ":";
    if(gHlxBranches.has_is_beam && it < static_cast<Int_t>(gHlxEvent.is_beam.size()))
      tl << (gHlxEvent.is_beam[it] ? " beam" : " scat");
    if(gHlxBranches.has_is_accidental && it < static_cast<Int_t>(gHlxEvent.is_accidental.size()) &&
       gHlxEvent.is_accidental[it])
      tl << " acc";
    const UInt_t pv =
      (it < static_cast<Int_t>(gHlxEvent.pid.size())) ? static_cast<UInt_t>(gHlxEvent.pid[it]) : 0U;
    tl << " pid=0x" << std::hex << pv << std::dec << " (" << DecodePidCandidates(static_cast<Int_t>(pv)) << ")";
    if(gHlxBranches.has_mom0 && it < static_cast<Int_t>(gHlxEvent.mom0.size()))
      tl << " p=" << HelixFormatSigFig(gHlxEvent.mom0[it]) << " GeV/c";
    if(gHlxBranches.has_charge && it < static_cast<Int_t>(gHlxEvent.charge.size()))
      tl << " q=" << gHlxEvent.charge[it];
    if(gHlxBranches.has_dEdx && it < static_cast<Int_t>(gHlxEvent.dEdx.size()))
      tl << " dEdx=" << HelixFormatSigFig(gHlxEvent.dEdx[it]);
    pt->AddText(tl.str().c_str());
  }

  if(gHlxBranches.has_closeDistTpc && gHlxEvent.ntTpc == 2) {
    const Double_t d = HelixMinCloseDist(gHlxEvent);
    if(d >= 0.0) {
      std::ostringstream cl;
      cl << "closeDist #approx " << HelixFormatSigFig(d) << " mm";
      pt->AddText(cl.str().c_str());
    }
  }

  if(gHlxBranches.has_k0_mass && !gHlxEvent.k0_mass.empty()) {
    std::ostringstream k0s;
    k0s << "K^{0} mass:";
    for(Double_t m : gHlxEvent.k0_mass) {
      if(TMath::Finite(m))
        k0s << " " << HelixFormatSigFig(m) << " GeV/c^{2}";
    }
    pt->AddText(k0s.str().c_str());
  }

  pt->Draw("same");
  gHlxOverlay.push_back(pt);
}

//______________________________________________________________________________
void
print_helix_event_summary()
{
  std::cout << "\n=== Helix 3D summary ===" << std::endl;
  std::cout << "Run=" << gHlxEvent.runnum << " Event=" << gHlxEvent.evnum << " entry=" << gHlxCurrentEntry
            << " ntTpc=" << gHlxEvent.ntTpc;
  if(gHlxBranches.has_effective_ntTpc)
    std::cout << " effective_ntTpc=" << gHlxEvent.effective_ntTpc;
  if(gHlxBranches.has_beam_flag)
    std::cout << " beam_flag=" << gHlxEvent.beam_flag;
  std::cout << std::endl;

  for(Int_t it = 0; it < gHlxEvent.ntTpc; it++) {
    const Int_t nh = (it < static_cast<Int_t>(gHlxEvent.nhtrack.size())) ? gHlxEvent.nhtrack[it] : 0;
    const UInt_t pv = (it < static_cast<Int_t>(gHlxEvent.pid.size())) ? static_cast<UInt_t>(gHlxEvent.pid[it]) : 0U;
    std::cout << "  tr" << it << ":";
    if(gHlxBranches.has_is_beam && it < static_cast<Int_t>(gHlxEvent.is_beam.size()))
      std::cout << (gHlxEvent.is_beam[it] ? " [beam]" : " [scat]");
    if(gHlxBranches.has_is_accidental && it < static_cast<Int_t>(gHlxEvent.is_accidental.size()) &&
       gHlxEvent.is_accidental[it])
      std::cout << " [acc]";
    if(gHlxBranches.has_mom0 && it < static_cast<Int_t>(gHlxEvent.mom0.size()))
      std::cout << " mom0=" << gHlxEvent.mom0[it];
    if(gHlxBranches.has_charge && it < static_cast<Int_t>(gHlxEvent.charge.size()))
      std::cout << " charge=" << gHlxEvent.charge[it];
    if(gHlxBranches.has_dEdx && it < static_cast<Int_t>(gHlxEvent.dEdx.size()))
      std::cout << " dEdx=" << gHlxEvent.dEdx[it];
    std::cout << " nhits=" << nh << " pid=0x" << std::hex << pv << std::dec << " ("
              << DecodePidCandidates(static_cast<Int_t>(pv)) << ")" << std::endl;
  }

  if(gHlxBranches.has_closeDistTpc && gHlxEvent.ntTpc == 2) {
    const Double_t d = HelixMinCloseDist(gHlxEvent);
    if(d >= 0.0)
      std::cout << "  closeDist (min pair) = " << d << " mm" << std::endl;
  }

  if(gHlxBranches.has_k0_mass && !gHlxEvent.k0_mass.empty()) {
    std::cout << "  k0_mass (" << gHlxEvent.k0_mass.size() << "):";
    for(Double_t m : gHlxEvent.k0_mass) {
      if(TMath::Finite(m))
        std::cout << " " << m;
    }
    std::cout << " GeV/c^2" << std::endl;
  }

  std::cout << "Helix TPolyLine3D drawn: " << gHlxHelixLines.size() << std::endl;
  std::cout << "========================\n" << std::endl;
}

//______________________________________________________________________________
void
draw_helix_tracks_3d()
{
  if(!gHlxCanvas)
    return;
  gHlxCanvas->cd();

  for(Int_t itrack = 0; itrack < gHlxEvent.ntTpc; itrack++) {
    if(itrack >= static_cast<Int_t>(gHlxEvent.helix_cx.size()))
      break;

    const Double_t cx = gHlxEvent.helix_cx[itrack];
    const Double_t cy = gHlxEvent.helix_cy[itrack];
    const Double_t z0 = gHlxEvent.helix_z0[itrack];
    const Double_t rr = gHlxEvent.helix_r[itrack];
    const Double_t dz = gHlxEvent.helix_dz[itrack];
    UInt_t pidv = 0;
    if(itrack < static_cast<Int_t>(gHlxEvent.pid.size()))
      pidv = static_cast<UInt_t>(gHlxEvent.pid[itrack]);
    const Int_t col = PidLineColor(pidv);

    Double_t tmin = 0.0, tmax = 0.0;
    const Bool_t haveTrange = HelixDrawThetaRange(gHlxEvent, itrack, tmin, tmax);

    if(haveTrange && TMath::Finite(rr)) {
      const Double_t tspan = tmax - tmin;
      Int_t nHelixSample = static_cast<Int_t>(
        TMath::Min(640.0, TMath::Max(96.0, 64.0 + 80.0 * tspan / TMath::TwoPi())));
      DrawHelixPolylineSegments(cx, cy, z0, rr, dz, tmin, tmax, col, nHelixSample);
      const Double_t dTgt = HelixMinDistToTargetMm(cx, cy, z0, rr, dz, tmin, tmax);
      if(TMath::Finite(dTgt))
        std::cout << "  [helix] tr" << itrack << " min|helix-target|=" << dTgt << " mm (drawn theta range)"
                  << std::endl;
    } else {
      std::cout << "  [helix line] tr" << itrack << " skipped: haveTrange=" << haveTrange
                << " finite_r=" << TMath::Finite(rr) << std::endl;
    }

    if(itrack < static_cast<Int_t>(gHlxEvent.nhtrack.size()) && gHlxEvent.nhtrack[itrack] > 0 &&
       itrack < static_cast<Int_t>(gHlxEvent.hitpos_x.size())) {
      const Int_t nhits = gHlxEvent.nhtrack[itrack];
      const Int_t nstore = TMath::Min(nhits, static_cast<Int_t>(gHlxEvent.hitpos_x[itrack].size()));
      if(nstore > 0) {
        TPolyMarker3D* hits = new TPolyMarker3D(nstore);
        hits->SetMarkerStyle(20);
        hits->SetMarkerSize(0.42);
        hits->SetMarkerColor(col);
        for(Int_t i = 0; i < nstore; i++) {
          Double_t dx, dy, dz;
          LocalHitToDisplay(gHlxEvent.hitpos_x[itrack][i], gHlxEvent.hitpos_y[itrack][i],
                            gHlxEvent.hitpos_z[itrack][i], dx, dy, dz);
          hits->SetPoint(i, dx, dy, dz);
        }
        hits->Draw("same");
        gHlxHitMarkers.push_back(hits);
      }
    }

    if(gHlxDrawCalpos && gHlxBranches.has_calpos && itrack < static_cast<Int_t>(gHlxEvent.calpos_x.size()) &&
       itrack < static_cast<Int_t>(gHlxEvent.nhtrack.size()) && gHlxEvent.nhtrack[itrack] > 0) {
      const Int_t ncal = TMath::Min(gHlxEvent.nhtrack[itrack], static_cast<Int_t>(gHlxEvent.calpos_x[itrack].size()));
      if(ncal > 0) {
        TPolyMarker3D* cals = new TPolyMarker3D(ncal);
        cals->SetMarkerStyle(4);
        cals->SetMarkerSize(0.55);
        cals->SetMarkerColor(col);
        for(Int_t i = 0; i < ncal; i++) {
          Double_t dx, dy, dz;
          LocalHitToDisplay(gHlxEvent.calpos_x[itrack][i], gHlxEvent.calpos_y[itrack][i],
                            gHlxEvent.calpos_z[itrack][i], dx, dy, dz);
          cals->SetPoint(i, dx, dy, dz);
        }
        cals->Draw("same");
        gHlxHitMarkers.push_back(cals);
      }
    }
  }

  draw_helix_vertex_marker();
  draw_helix_info_overlay();
}

//______________________________________________________________________________
void
helix_event(Long64_t evnum = -1)
{
  if(!load_helix_event(evnum))
    return;

  if(gHlxCanvas && gHlxCanvas->IsZombie()) {
    gHlxCanvas = nullptr;
    gHlxView = nullptr;
  }

  if(!gHlxCanvas) {
    gStyle->SetCanvasPreferGL(kFALSE); // X11/SSH 転送では GL より軽いことが多い
    gHlxCanvas = new TCanvas("cHlx3d", "TPC Helix 3D View", 900, 900);
  }

  gHlxCanvas->cd();
  gHlxCanvas->Clear();
  clear_helix_objects();

  draw_tpc_frame_and_corner_axes();
  draw_helix_tracks_3d();

  gHlxCanvas->Modified();
  gHlxCanvas->Update();
  print_helix_event_summary();
}

//______________________________________________________________________________
Long64_t
helix_next(Long64_t max_scan = -1)
{
  if(!gHlxReader) {
    std::cerr << "Error: helix_set_path() first." << std::endl;
    return -1;
  }

  const Long64_t nentries = gHlxReader->GetEntries(false);
  Long64_t start = gHlxCfg.scan_start;
  if(start < 0)
    start = 0;
  if(start >= nentries) {
    std::cout << "helix_next: scan_start " << start << " >= nentries " << nentries << std::endl;
    return -1;
  }

  Long64_t end = nentries;
  if(max_scan >= 0)
    end = TMath::Min(nentries, start + max_scan);

  HelixEventData snap;
  for(Long64_t e = start; e < end; ++e) {
    if(!helix_load_filter_snapshot(e, snap))
      continue;
    if(helix_matches_filters(snap)) {
      gHlxCfg.scan_start = e + 1;
      helix_event(e);
      std::cout << "helix_next: matched entry " << e << " (next scan_start=" << gHlxCfg.scan_start << ")"
                << std::endl;
      return e;
    }
  }

  std::cout << "helix_next: no match in [" << start << ", " << end << ")" << std::endl;
  return -1;
}

//______________________________________________________________________________
void
helix_browse_filtered(Long64_t maxSteps = -1)
{
  if(!gHlxReader) {
    std::cerr << "Error: helix_set_path() してから helix_browse_filtered() を呼んでください。" << std::endl;
    return;
  }

  std::cout << "\nフィルタ付き連続表示: Enter=helix_next(), q=終了\n"
               "（entry 直指定は helix_event(n)）\n"
            << std::endl;

  char line[32];
  for(Long64_t step = 0; maxSteps < 0 || step < maxSteps; ++step) {
    if(helix_next() < 0)
      break;
    std::cout << "[browse_filtered] Enter=次, q=終了 > " << std::flush;
    if(!std::fgets(line, sizeof(line), stdin))
      break;
    if(line[0] == 'q' || line[0] == 'Q')
      break;
  }
  std::cout << "helix_browse_filtered: 終了しました。" << std::endl;
}

//______________________________________________________________________________
void
helix_browse(Long64_t maxSteps = -1)
{
  if(!gHlxReader) {
    std::cerr << "Error: helix_set_path() してから helix_browse() を呼んでください。" << std::endl;
    return;
  }

  std::cout << "\n連続表示: ターミナルで Enter を押すと次のイベント（ランダム）、"
               "q または Q を入力して Enter で終了します。\n"
               "（Ctrl+D でも終了）\n"
            << std::endl;

  char line[32];
  for(Long64_t step = 0; maxSteps < 0 || step < maxSteps; ++step) {
    helix_event(-1);
    std::cout << "[browse] Enter=次, q=終了 > " << std::flush;
    if(!std::fgets(line, sizeof(line), stdin))
      break;
    if(line[0] == 'q' || line[0] == 'Q')
      break;
  }
  std::cout << "helix_browse: 終了しました。" << std::endl;
}

//______________________________________________________________________________
// .L / ACLiC 読み込み直後に空の 3D フレームだけ表示（以前の挙動を復元）
void
helix_open_gui()
{
  if(gROOT && gROOT->IsBatch())
    return;

  if(gHlxCanvas && !gHlxCanvas->IsZombie()) {
    gHlxCanvas->cd();
    gHlxCanvas->Modified();
    gHlxCanvas->Update();
    return;
  }

  gStyle->SetCanvasPreferGL(kFALSE);
  gHlxCanvas = new TCanvas("cHlx3d", "TPC Helix 3D View", 900, 900);
  gHlxCanvas->cd();
  gHlxCanvas->Clear();
  clear_helix_objects();
  draw_tpc_frame_and_corner_axes();
  gHlxCanvas->Modified();
  gHlxCanvas->Update();
  std::cout << "[helix] GUI ready (empty frame). helix_set_path(\"file.root\"); helix_event(0);" << std::endl;
}

helix_open_gui();
