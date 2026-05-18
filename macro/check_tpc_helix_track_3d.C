// -*- C++ -*-
// DstTPCHelixTracking 出力（ツリー "tpc"）の 3D 表示。ベースは check_tpc_track_3d.C。
//
// Usage:
//   root -l
//   .L check_tpc_helix_track_3d.C+
//   helix_set_path("file.root")
//   helix_event()           // 1 イベント（ランダム）
//   helix_event(123)        // entry 固定
//   helix_browse()          // 連続: Enter で次のランダムイベント、q+Enter で終了
//   helix_browse(500)       // 最大 500 回まで（Ctrl+D でも終了）
//
// 必須枝: run_number, event_number, ntTpc, nhtrack, pid,
//         helix_cx, helix_cy, helix_z0, helix_r, helix_dz, helix_t,
//         hitpos_x, hitpos_y, hitpos_z
//
// 座標:
//   螺旋は LocalToGlobal 後の局所 (x,y,z) をツリーと同じ系で扱い、ROOT 描画には
//     (Xv,Yv,Zv) = (z, x, y)
//   を渡す。角矢印のラベルは局所軸名 z（赤）, x（緑）, y（青）— ROOT の表示軸名 X,Y,Z とは対応が逆転している点に注意。
// 螺旋の「線」は helix_event() 内で各トラック new TPolyLine3D → SetPoint → Draw("same")。
// draw_tpc_frame_and_corner_axes() の直後に同じ gHlxCanvas へ重ねる（Draw() 単独だと 3D で消える）。
// θ 区間: 枝 helix_theta_min/max（または helix_t_min/max）があれば優先、なければ helix_t から
// HelixThetaDrawRange。点数は θ の長さに応じ 72〜512 程度。

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLatex.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>

// tpc::Z_TARGET（include/TPCPadHelper.hh）と同値。相対 include に依存しないようマクロ内で定義。
static constexpr Double_t kTPC_Z_TARGET = -143.0;

TFile* gHlxFile = nullptr;
TTreeReader* gHlxReader = nullptr;
TRandom3* gHlxRandom = nullptr;
TCanvas* gHlxCanvas = nullptr;
TView* gHlxView = nullptr;
// 0: なし  1: helix_theta_min/max（現行 Dst）  2: helix_t_min/max（旧名）
static Int_t gHlxHelixTRangeMode = 0;

struct HelixEventData {
  UInt_t runnum = 0;
  UInt_t evnum = 0;
  Int_t ntTpc = 0;
  std::vector<Int_t> nhtrack;
  std::vector<Int_t> pid;
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
};

HelixEventData gHlxEvent;

std::vector<TPolyLine3D*> gHlxHelixLines;
std::vector<TPolyMarker3D*> gHlxHitMarkers;
std::vector<TPolyLine3D*> gHlxFrame;

//______________________________________________________________________________
void
DrawAxisLabel3D(Double_t x, Double_t y, Double_t z, const char* label, Color_t color, TView* view)
{
  Double_t wc[3] = {x, y, z};
  Double_t ndc[3];
  view->WCtoNDC(wc, ndc);
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextColor(color);
  latex->SetTextSize(0.022);
  latex->SetTextAlign(22);
  latex->DrawLatex(ndc[0], ndc[1], label);
}

// 螺旋: TPCLocalTrackHelix.cc の LocalPosition(par,t) で (xi,yi,zi)、
//       LocalToGlobal で局所 (x,y,z)=(-xi, zi, yi+Z_TARGET)（TPCLTrackHit::GetLocalCalPosHelix と同じ）。
// 可視化: track 3D マクロと同じ (dispX,dispY,dispZ)=(局所 z, 局所 x, 局所 y)。
inline void
HelixPointToDisplay(Double_t cx, Double_t cy, Double_t z0, Double_t r, Double_t dz, Double_t t,
                    Double_t& dispX, Double_t& dispY, Double_t& dispZ)
{
  const Double_t xi = cx + r * TMath::Cos(t);
  const Double_t yi = cy + r * TMath::Sin(t);
  const Double_t zi = z0 + dz * r * t;
  const Double_t xloc = -xi;
  const Double_t yloc = zi;
  const Double_t zloc = yi + kTPC_Z_TARGET;
  dispX = zloc;
  dispY = xloc;
  dispZ = yloc;
}

inline void
LocalHitToDisplay(Double_t xloc, Double_t yloc, Double_t zloc, Double_t& dispX, Double_t& dispY,
                  Double_t& dispZ)
{
  dispX = zloc;
  dispY = xloc;
  dispZ = yloc;
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

// helix_t から描画用 [t_lo, t_hi] を決める。単純 min/max だと外れ値・2pi 巻きで区間が巨大になり、
// サンプルが足りず螺旋が棒状に潰れる。対策: (1) 先頭〜末尾 θ 優先の条件を緩める
// (2) 6–94% タイルで外れ値を間引く (3) 端をわずかに延ばす。
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
void
clear_helix_objects()
{
  for(auto* obj : gHlxHelixLines) {
    if(obj)
      delete obj;
  }
  gHlxHelixLines.clear();

  for(auto* obj : gHlxHitMarkers) {
    if(obj)
      delete obj;
  }
  gHlxHitMarkers.clear();

  for(auto* obj : gHlxFrame) {
    if(obj)
      delete obj;
  }
  gHlxFrame.clear();
}

//______________________________________________________________________________
void
draw_tpc_frame_and_corner_axes()
{
  const Double_t xmin_orig = -311.0, xmax_orig = 311.0;
  const Double_t ymin_orig = -310.0, ymax_orig = 310.0;
  const Double_t zmin_orig = -311.0, zmax_orig = 311.0;

  const Double_t xmin = zmin_orig, xmax = zmax_orig;
  const Double_t ymin = xmin_orig, ymax = xmax_orig;
  const Double_t zmin = ymin_orig, zmax = ymax_orig;

  gHlxView = TView::CreateView(1, 0, 0);
  gHlxView->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);

  Double_t axisOriginX = xmin + (xmax - xmin) * 0.05;
  Double_t axisOriginY = ymin + (ymax - ymin) * 0.05;
  Double_t axisOriginZ = zmin + (zmax - zmin) * 0.05;
  Double_t axisLength = TMath::Min(TMath::Min(xmax - xmin, ymax - ymin), zmax - zmin) * 0.1;

  // 以下、axisX 等の変数名は「表示座標の第 1〜3 成分」を動かす矢印（局所 z,x,y 方向）
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
      octagon->SetPoint(i, octagon_z_orig[idx], octagon_x_orig[idx], y_orig);
    }
    push_line(octagon);
  }

  for(Int_t i = 0; i < nOctagonSides; i++) {
    TPolyLine3D* edge = new TPolyLine3D(2);
    edge->SetLineColor(kGray + 1);
    edge->SetLineWidth(2);
    edge->SetPoint(0, octagon_z_orig[i], octagon_x_orig[i], -frameHalfHeight);
    edge->SetPoint(1, octagon_z_orig[i], octagon_x_orig[i], frameHalfHeight);
    push_line(edge);
  }

  // ターゲット円筒: 母線=局所 y、断面=局所 x-z（中心 z≈targetZ_orig）。check_tpc_track_3d.C と同一。
  const Double_t targetZ_orig = -143.0;
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
      circle->SetPoint(i, z_orig, x_orig, y_orig);
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
    line->SetPoint(0, z_orig, x_orig, -targetHalfHeight);
    line->SetPoint(1, z_orig, x_orig, targetHalfHeight);
    push_line(line);
  }

  const Double_t loff = axisLength * 0.25;
  const Double_t lxlab = axisOriginX + axisLength + loff;
  DrawAxisLabel3D(lxlab, axisOriginY, axisOriginZ, "z", kRed, gHlxView);
  DrawAxisLabel3D(axisOriginX, axisOriginY + axisLength + loff, axisOriginZ, "x", kGreen + 2, gHlxView);
  DrawAxisLabel3D(axisOriginX, axisOriginY, axisOriginZ + axisLength + loff, "y", kBlue, gHlxView);
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

  if((tree->GetBranch("helix_theta_min") != nullptr) && (tree->GetBranch("helix_theta_max") != nullptr)) {
    gHlxHelixTRangeMode = 1;
  } else if((tree->GetBranch("helix_t_min") != nullptr) && (tree->GetBranch("helix_t_max") != nullptr)) {
    gHlxHelixTRangeMode = 2;
  } else {
    gHlxHelixTRangeMode = 0;
  }

  gHlxReader = new TTreeReader(tree);
  std::cout << "File opened: " << path << std::endl;
  std::cout << "Total entries: " << tree->GetEntries() << std::endl;
  if(gHlxHelixTRangeMode == 1) {
    std::cout << "  helix_theta_min / helix_theta_max: present" << std::endl;
  } else if(gHlxHelixTRangeMode == 2) {
    std::cout << "  helix_t_min / helix_t_max: present (legacy branch names)" << std::endl;
  } else {
    std::cout << "  track theta range branches: absent (use per-hit helix_t)" << std::endl;
  }

  if(!gHlxRandom) {
    gHlxRandom = new TRandom3();
    gHlxRandom->SetSeed(0);
  }
}

//______________________________________________________________________________
void
load_helix_event(Long64_t entry = -1)
{
  if(!gHlxReader) {
    std::cerr << "Error: No file opened. Use helix_set_path() first." << std::endl;
    return;
  }

  Long64_t nentries = gHlxReader->GetEntries(false);
  if(nentries == 0) {
    std::cerr << "Error: No entries in tree" << std::endl;
    return;
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
    return;
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

  TTreeReaderValue<std::vector<Double_t>>* pHelixTmin = nullptr;
  TTreeReaderValue<std::vector<Double_t>>* pHelixTmax = nullptr;
  if(gHlxHelixTRangeMode == 1) {
    pHelixTmin = new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_theta_min");
    pHelixTmax = new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_theta_max");
  } else if(gHlxHelixTRangeMode == 2) {
    pHelixTmin = new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_t_min");
    pHelixTmax = new TTreeReaderValue<std::vector<Double_t>>(*gHlxReader, "helix_t_max");
  }

  gHlxReader->SetEntry(entry);

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

  gHlxEvent.helix_t_min.clear();
  gHlxEvent.helix_t_max.clear();
  if(pHelixTmin && pHelixTmax) {
    gHlxEvent.helix_t_min = *(*pHelixTmin);
    gHlxEvent.helix_t_max = *(*pHelixTmax);
    delete pHelixTmin;
    delete pHelixTmax;
  }

  std::cout << "Helix event loaded: Run=" << gHlxEvent.runnum << ", Event=" << gHlxEvent.evnum
            << ", Entry=" << entry << ", ntTpc=" << gHlxEvent.ntTpc << std::endl;
}

//______________________________________________________________________________
void
helix_event(Long64_t evnum = -1)
{
  load_helix_event(evnum);
  clear_helix_objects();

  if(!gHlxCanvas) {
    gHlxCanvas = new TCanvas("cHlx3d", "TPC Helix 3D View", 800, 800);
  } else {
    gHlxCanvas->Clear();
  }

  draw_tpc_frame_and_corner_axes();

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
    Bool_t haveTrange = kFALSE;
    if(gHlxHelixTRangeMode != 0 && itrack < static_cast<Int_t>(gHlxEvent.helix_t_min.size()) &&
       itrack < static_cast<Int_t>(gHlxEvent.helix_t_max.size())) {
      const Double_t a = gHlxEvent.helix_t_min[itrack];
      const Double_t b = gHlxEvent.helix_t_max[itrack];
      if(TMath::Finite(a) && TMath::Finite(b) && TMath::Abs(b - a) > 1e-12) {
        const Double_t lo = TMath::Min(a, b);
        const Double_t hi = TMath::Max(a, b);
        const Double_t mrg = (hi - lo) * 0.02 + 5e-4;
        tmin = lo - mrg;
        tmax = hi + mrg;
        haveTrange = kTRUE;
      }
    }
    if(!haveTrange) {
      haveTrange =
        (itrack < static_cast<Int_t>(gHlxEvent.helix_t.size()) &&
         HelixThetaDrawRange(gHlxEvent.helix_t[itrack], tmin, tmax));
    }

    // 半径しきい値で線を捨てない（Dst の r が mm で小さいと従来は一切描かれなかった）
    if(haveTrange && TMath::Finite(rr)) {
      const Double_t tspan = tmax - tmin;
      Int_t nHelixSample = static_cast<Int_t>(
        TMath::Min(512.0, TMath::Max(72.0, 48.0 + 64.0 * tspan / TMath::TwoPi())));
      if(nHelixSample < 2)
        nHelixSample = 2;
      const Double_t helixDen = static_cast<Double_t>(nHelixSample - 1);
      TPolyLine3D* hel = new TPolyLine3D(nHelixSample);
      hel->SetLineColor(col);
      hel->SetLineWidth(3);
      for(Int_t s = 0; s < nHelixSample; s++) {
        const Double_t t = tmin + (tmax - tmin) * static_cast<Double_t>(s) / helixDen;
        Double_t dx, dy, dz;
        HelixPointToDisplay(cx, cy, z0, rr, dz, t, dx, dy, dz);
        hel->SetPoint(s, dx, dy, dz);
      }
      hel->Draw("same");
      gHlxHelixLines.push_back(hel);
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
        hits->SetMarkerSize(0.55);
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
  }

  gHlxCanvas->Update();

  std::cout << "\n=== Helix 3D summary ===" << std::endl;
  std::cout << "Tracks: " << gHlxEvent.ntTpc << std::endl;
  for(Int_t it = 0; it < gHlxEvent.ntTpc; it++) {
    const Int_t nh = (it < static_cast<Int_t>(gHlxEvent.nhtrack.size())) ? gHlxEvent.nhtrack[it] : 0;
    const UInt_t pv = (it < static_cast<Int_t>(gHlxEvent.pid.size())) ? static_cast<UInt_t>(gHlxEvent.pid[it]) : 0U;
    std::cout << "  tr" << it << " nhits=" << nh << " pid=0x" << std::hex << pv << std::dec << std::endl;
  }
  std::cout << "Helix TPolyLine3D drawn (this pad): " << gHlxHelixLines.size() << std::endl;
  std::cout << "========================\n" << std::endl;
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
