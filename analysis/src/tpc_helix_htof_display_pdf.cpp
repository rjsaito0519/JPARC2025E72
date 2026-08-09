// -*- C++ -*-
// Helix × HTOF event display → multi-page PDF (XZ + 3D).
//
// Usage:
//   tpc_helix_htof_display_pdf -n <N> [--start <htof_entry>] [--match-only|--no-match-only]
//     [-o <out.pdf>] <htof.root> <helix.root>
//
// Inputs:
//   htof.root  … tree "htof" from DstTPCHelixHTOF
//   helix.root … tree "tpc"  from DstTPCHelixTracking
//   Joined by (run_number, event_number). Default: --match-only.
//
// Coordinates (same as check_tpc_helix_track_3d.C):
//   Helix internal → TPC global: x=-xi, y=zi, z=yi+Z_TARGET
//   XZ: horizontal=z_tpc, vertical=x_tpc, aspect 1:1
//   3D SetPoint: X=z_tpc, Y=x_tpc, Z=y_tpc
//
// 3D view: TPC zoom (no BH2, no HTOF). Octagon + clipped layer rings + target.
// HTOF segs on XZ only (inner L=348.6 + outer L+10 mm).
// Tracks: meas = solid [tmin,tmax] on top; extrap = magenta continuous vtx→HTOF
// (helix sample with endpoints forced to DST vtx/extrap xyz).

#include "paths.h"
#include "ana_helper.h"
#include "TPCPadHelper.hh"

#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TView.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr Double_t kZTarget = -143.0;
constexpr Double_t kTargetRadius = 60.0;
constexpr Double_t kTargetHalfY = 60.0;

// check_tpc_helix_track_3d.C frame
constexpr Double_t kFrameFaceToFace = 622.0;
constexpr Double_t kFrameHeight = 620.0;
constexpr Double_t kFrameHalfHeight = kFrameHeight / 2.0;
constexpr Double_t kFrameRadius = kFrameFaceToFace / 2.0;
constexpr Double_t kOctagonVertexR = kFrameRadius / 0.92387953251; // /cos(pi/8)

constexpr Int_t kNPlanesHTOF = 8;
constexpr Double_t kHtofL = 348.6;          // inner face (= ExtrapolateToHTOF plane)
constexpr Double_t kHtofThickness = 10.0;   // 1 cm
constexpr Double_t kHtofLOuter = kHtofL + kHtofThickness;
constexpr Double_t kHtofSegWidth = 71.0;   // 70 + 1
constexpr Double_t kHtofSegHeight = 400.0;
constexpr Double_t kHtofYOffset = 4.0;
constexpr Double_t kHtofBeamWinW = 2.0 * kHtofSegWidth;
constexpr Double_t kHtofBeamWinH = 112.0;

const Color_t kExtrapColor = kMagenta + 1;

constexpr Int_t kNSegBH2 = 15;
constexpr Double_t kBH2SegWidth = 14.0;
constexpr Double_t kBH2SegGap = 0.5;
constexpr Double_t kBH2SegPitch = kBH2SegWidth + kBH2SegGap;
constexpr Double_t kBH2SegCenter = 6.0;
constexpr Double_t kBH2Thickness = 5.0;
constexpr Double_t kBH2Z = -745.7;

constexpr Double_t kXZHalf = 800.0; // XZ includes BH2
constexpr Double_t kViewMargin = 0.22;

void
usage(const char* argv0)
{
  std::cerr
    << "Usage: " << argv0
    << " -n <N> [--start <htof_entry>] [--match-only|--no-match-only]\n"
    << "  [-o <out.pdf>] <htof.root> <helix.root>\n";
}

inline void
HelixInternalToTPCGlobal(Double_t xi, Double_t yi, Double_t zi,
                         Double_t& xg, Double_t& yg, Double_t& zg)
{
  xg = -xi;
  yg = zi;
  zg = yi + kZTarget;
}

inline void
HelixPointTPC(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz, Double_t t,
              Double_t& xg, Double_t& yg, Double_t& zg)
{
  HelixInternalToTPCGlobal(cx + rr * TMath::Cos(t),
                           cy + rr * TMath::Sin(t),
                           z0 + dz * rr * t,
                           xg, yg, zg);
}

inline void
TPCGlobalToDisplay3D(Double_t xg, Double_t yg, Double_t zg,
                     Double_t& dx, Double_t& dy, Double_t& dz)
{
  dx = zg;
  dy = xg;
  dz = yg;
}

inline Double_t
BH2SegX(Double_t seg)
{
  return kBH2SegPitch * (seg - kBH2SegCenter);
}

Color_t
TrackColor(Bool_t is_beam, Bool_t is_acc)
{
  if (is_beam) return kCyan + 1;
  if (is_acc) return kGray + 1;
  return kGreen + 2;
}

std::vector<Double_t>
LayerBoundaryRadii()
{
  std::vector<Double_t> radii;
  radii.reserve(tpc::NumOfLayersTPC + 1);
  radii.push_back(tpc::padParameter[0][tpc::kRadius] - 0.5 * tpc::padParameter[0][tpc::kLength]);
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer)
    radii.push_back(tpc::padParameter[layer][tpc::kRadius]
                    + 0.5 * tpc::padParameter[layer][tpc::kLength]);
  return radii;
}

void
SampleHelixTPC(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz,
               Double_t tmin, Double_t tmax, Int_t nSample,
               std::vector<Double_t>& xg, std::vector<Double_t>& yg, std::vector<Double_t>& zg)
{
  xg.clear();
  yg.clear();
  zg.clear();
  if (nSample < 2) nSample = 2;
  if (!std::isfinite(tmin) || !std::isfinite(tmax) || tmin == tmax)
    return;
  xg.resize(nSample);
  yg.resize(nSample);
  zg.resize(nSample);
  const Double_t den = static_cast<Double_t>(nSample - 1);
  for (Int_t i = 0; i < nSample; ++i) {
    const Double_t t = tmin + (tmax - tmin) * static_cast<Double_t>(i) / den;
    HelixPointTPC(cx, cy, z0, rr, dz, t, xg[i], yg[i], zg[i]);
  }
}

Bool_t
FindNearestThetaInWindow(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz,
                         Double_t xt, Double_t yt, Double_t zt,
                         Double_t tlo, Double_t thi,
                         Double_t& tbest, Double_t& bestD2)
{
  if (!std::isfinite(xt) || !std::isfinite(yt) || !std::isfinite(zt))
    return kFALSE;
  if (!std::isfinite(tlo) || !std::isfinite(thi))
    return kFALSE;
  if (tlo > thi) std::swap(tlo, thi);
  if (thi - tlo < 1e-12) {
    Double_t xg, yg, zg;
    HelixPointTPC(cx, cy, z0, rr, dz, tlo, xg, yg, zg);
    bestD2 = (xg - xt) * (xg - xt) + (yg - yt) * (yg - yt) + (zg - zt) * (zg - zt);
    tbest = tlo;
    return kTRUE;
  }
  const Int_t nScan = TMath::Max(64, TMath::Nint(200.0 * (thi - tlo) / TMath::TwoPi()) + 1);
  bestD2 = 1e300;
  tbest = 0.5 * (tlo + thi);
  for (Int_t i = 0; i <= nScan; ++i) {
    const Double_t t = tlo + (thi - tlo) * static_cast<Double_t>(i) / nScan;
    Double_t xg, yg, zg;
    HelixPointTPC(cx, cy, z0, rr, dz, t, xg, yg, zg);
    const Double_t d2 = (xg - xt) * (xg - xt) + (yg - yt) * (yg - yt) + (zg - zt) * (zg - zt);
    if (d2 < bestD2) {
      bestD2 = d2;
      tbest = t;
    }
  }
  // Local refine around best (± one coarse step).
  const Double_t step = (thi - tlo) / nScan;
  const Double_t rlo = tbest - step;
  const Double_t rhi = tbest + step;
  for (Int_t i = 0; i <= 40; ++i) {
    const Double_t t = rlo + (rhi - rlo) * static_cast<Double_t>(i) / 40.0;
    Double_t xg, yg, zg;
    HelixPointTPC(cx, cy, z0, rr, dz, t, xg, yg, zg);
    const Double_t d2 = (xg - xt) * (xg - xt) + (yg - yt) * (yg - yt) + (zg - zt) * (zg - zt);
    if (d2 < bestD2) {
      bestD2 = d2;
      tbest = t;
    }
  }
  return bestD2 < 1e299;
}

// Map a TPC point to helix theta near the measured arc (avoid ±4π wrong-period).
// Tries [tmin,tmax], outward extensions, and ±π around mid — keeps the closest.
Bool_t
FindNearestThetaNearMeas(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz,
                         Double_t xt, Double_t yt, Double_t zt,
                         Double_t tmin, Double_t tmax,
                         Double_t& tbest, Double_t& bestD2)
{
  if (!std::isfinite(tmin) || !std::isfinite(tmax))
    return kFALSE;
  if (tmin > tmax) std::swap(tmin, tmax);
  const Double_t twopi = TMath::TwoPi();
  const Double_t mid = 0.5 * (tmin + tmax);
  const Double_t windows[][2] = {
    {tmin, tmax},
    {tmax, tmax + twopi},
    {tmin - twopi, tmin},
    {mid - TMath::Pi(), mid + TMath::Pi()},
    {tmin - 0.5 * twopi, tmax + 0.5 * twopi},
  };
  bestD2 = 1e300;
  tbest = mid;
  Bool_t any = kFALSE;
  for (const auto& w : windows) {
    Double_t tb = 0., d2 = 0.;
    if (!FindNearestThetaInWindow(cx, cy, z0, rr, dz, xt, yt, zt, w[0], w[1], tb, d2))
      continue;
    if (d2 < bestD2) {
      bestD2 = d2;
      tbest = tb;
      any = kTRUE;
    }
  }
  return any;
}

TGraph*
MakeGraphXZ(const std::vector<Double_t>& xg, const std::vector<Double_t>& zg)
{
  const Int_t n = static_cast<Int_t>(std::min(xg.size(), zg.size()));
  if (n < 2) return nullptr;
  auto* g = new TGraph(n);
  for (Int_t i = 0; i < n; ++i)
    g->SetPoint(i, zg[i], xg[i]);
  return g;
}

void
DrawHelix3D(const std::vector<Double_t>& xg,
            const std::vector<Double_t>& yg,
            const std::vector<Double_t>& zg,
            Color_t col, Style_t style = 1, Width_t width = 3)
{
  const Int_t n = static_cast<Int_t>(std::min({xg.size(), yg.size(), zg.size()}));
  for (Int_t i = 0; i < n - 1; ++i) {
    if (!std::isfinite(xg[i]) || !std::isfinite(xg[i + 1])) continue;
    Double_t d0[3], d1[3];
    TPCGlobalToDisplay3D(xg[i], yg[i], zg[i], d0[0], d0[1], d0[2]);
    TPCGlobalToDisplay3D(xg[i + 1], yg[i + 1], zg[i + 1], d1[0], d1[1], d1[2]);
    auto* seg = new TPolyLine3D(2);
    seg->SetLineColor(col);
    seg->SetLineStyle(style);
    seg->SetLineWidth(width);
    seg->SetPoint(0, d0[0], d0[1], d0[2]);
    seg->SetPoint(1, d1[0], d1[1], d1[2]);
    seg->Draw();
  }
}

void
DrawBH2XZ(Double_t time0_seg)
{
  const Double_t z0 = kBH2Z - 0.5 * kBH2Thickness;
  const Double_t z1 = kBH2Z + 0.5 * kBH2Thickness;
  for (Int_t s = 0; s < kNSegBH2; ++s) {
    const Double_t xc = BH2SegX(static_cast<Double_t>(s));
    auto* box = new TBox(z0, xc - 0.5 * kBH2SegWidth, z1, xc + 0.5 * kBH2SegWidth);
    const Bool_t is_t0 = std::isfinite(time0_seg) && (TMath::Nint(time0_seg) == s);
    box->SetFillColor(is_t0 ? kOrange + 7 : kGray + 1);
    box->SetFillStyle(is_t0 ? 1001 : 3001);
    box->SetLineColor(is_t0 ? kOrange + 9 : kGray + 2);
    box->SetLineWidth(is_t0 ? 2 : 1);
    box->Draw("l");
  }
}

void
DrawTargetXZ()
{
  auto* el = new TEllipse(kZTarget, 0., kTargetRadius, kTargetRadius);
  el->SetFillStyle(0);
  el->SetLineColor(kMagenta + 1);
  el->SetLineWidth(2);
  el->Draw();
}

Bool_t
InsideTPCOctagon(Double_t x, Double_t z)
{
  for (Int_t i = 0; i < 8; ++i) {
    const Double_t a = static_cast<Double_t>(i) * TMath::Pi() / 4.0;
    if (x * TMath::Cos(a) + z * TMath::Sin(a) > kFrameRadius + 1e-3)
      return kFALSE;
  }
  return kTRUE;
}

void
DrawLayerRingsXZ()
{
  const auto radii = LayerBoundaryRadii();
  const Int_t nSeg = 180;
  for (Double_t r : radii) {
    std::vector<Double_t> zs, xs;
    auto flush = [&]() {
      if (zs.size() < 2) {
        zs.clear();
        xs.clear();
        return;
      }
      auto* g = new TGraph(static_cast<Int_t>(zs.size()));
      for (std::size_t i = 0; i < zs.size(); ++i)
        g->SetPoint(static_cast<Int_t>(i), zs[i], xs[i]);
      g->SetLineColor(kGray + 1);
      g->SetLineStyle(3);
      g->SetLineWidth(1);
      g->Draw("L");
      zs.clear();
      xs.clear();
    };
    for (Int_t i = 0; i <= nSeg; ++i) {
      const Double_t th = 2.0 * TMath::Pi() * static_cast<Double_t>(i) / nSeg;
      const Double_t x = r * TMath::Sin(th);
      const Double_t z = r * TMath::Cos(th) + kZTarget;
      if (InsideTPCOctagon(x, z)) {
        zs.push_back(z);
        xs.push_back(x);
      } else {
        flush();
      }
    }
    flush();
  }
}

void
DrawOctagonXZ()
{
  const Int_t n = 8;
  auto* g = new TGraph(n + 1);
  for (Int_t i = 0; i <= n; ++i) {
    const Double_t a = i * TMath::Pi() / 4.0 + TMath::Pi() / 8.0;
    const Double_t x = kOctagonVertexR * TMath::Cos(a);
    const Double_t z = kOctagonVertexR * TMath::Sin(a);
    g->SetPoint(i, z, x);
  }
  g->SetLineColor(kGray + 1);
  g->SetLineWidth(2);
  g->Draw("L");
}

// Plane local (u,v) → TPC global at radial distance L along plane normal.
// u = signed tangential, v = y - y_offset. L = kHtofL (inner) or kHtofLOuter.
void
HtofLocalToGlobal(Int_t iplane, Double_t u, Double_t v, Double_t L,
                  Double_t& xg, Double_t& yg, Double_t& zg)
{
  const Double_t phi = static_cast<Double_t>(iplane) * 0.25 * TMath::Pi();
  const Double_t nx = -TMath::Sin(phi);
  const Double_t nz = -TMath::Cos(phi);
  const Double_t tx = -TMath::Cos(phi);
  const Double_t tz = TMath::Sin(phi);
  xg = L * nx + u * tx;
  zg = L * nz + u * tz;
  yg = kHtofYOffset + v;
}

struct HtofCell {
  Int_t plane;
  Int_t seg;
  Double_t u0, u1, v0, v1;
};

// Cells matching ExtrapolateToHTOF segment assignment (display geometry).
std::vector<HtofCell>
BuildHtofCells()
{
  std::vector<HtofCell> cells;
  const Double_t w = kHtofSegWidth;
  const Double_t H = kHtofSegHeight;

  for (Int_t i = 0; i < kNPlanesHTOF; ++i) {
    if (i == 0) {
      // outer full-height
      cells.push_back({0, 0, -2 * w, -w, -H, H});
      cells.push_back({0, 5, w, 2 * w, -H, H});
      // inner split by y (beam-window core left as hole — not a cell)
      cells.push_back({0, 2, -w, 0., -H, 0.});
      cells.push_back({0, 1, -w, 0., 0., H});
      cells.push_back({0, 4, 0., w, -H, 0.});
      cells.push_back({0, 3, 0., w, 0., H});
    } else {
      cells.push_back({i, 2 + 4 * i, -2 * w, -w, -H, H});
      cells.push_back({i, 3 + 4 * i, -w, 0., -H, H});
      cells.push_back({i, 4 + 4 * i, 0., w, -H, H});
      cells.push_back({i, 5 + 4 * i, w, 2 * w, -H, H});
    }
  }
  return cells;
}

void
DrawRect3D(Double_t x0, Double_t y0, Double_t z0,
           Double_t x1, Double_t y1, Double_t z1,
           Double_t x2, Double_t y2, Double_t z2,
           Double_t x3, Double_t y3, Double_t z3,
           Color_t col, Width_t width)
{
  Double_t p[4][3];
  TPCGlobalToDisplay3D(x0, y0, z0, p[0][0], p[0][1], p[0][2]);
  TPCGlobalToDisplay3D(x1, y1, z1, p[1][0], p[1][1], p[1][2]);
  TPCGlobalToDisplay3D(x2, y2, z2, p[2][0], p[2][1], p[2][2]);
  TPCGlobalToDisplay3D(x3, y3, z3, p[3][0], p[3][1], p[3][2]);
  auto* pl = new TPolyLine3D(5);
  pl->SetLineColor(col);
  pl->SetLineWidth(width);
  for (Int_t i = 0; i < 4; ++i)
    pl->SetPoint(i, p[i][0], p[i][1], p[i][2]);
  pl->SetPoint(4, p[0][0], p[0][1], p[0][2]);
  pl->Draw();
}

void
DrawLine3DGlobal(Double_t x0, Double_t y0, Double_t z0,
                 Double_t x1, Double_t y1, Double_t z1,
                 Color_t col, Width_t width)
{
  Double_t d0[3], d1[3];
  TPCGlobalToDisplay3D(x0, y0, z0, d0[0], d0[1], d0[2]);
  TPCGlobalToDisplay3D(x1, y1, z1, d1[0], d1[1], d1[2]);
  auto* pl = new TPolyLine3D(2);
  pl->SetLineColor(col);
  pl->SetLineWidth(width);
  pl->SetPoint(0, d0[0], d0[1], d0[2]);
  pl->SetPoint(1, d1[0], d1[1], d1[2]);
  pl->Draw();
}

void
DrawHtofCellFace3D(const HtofCell& c, Double_t L, Color_t col, Width_t width)
{
  Double_t x00, y00, z00, x10, y10, z10, x11, y11, z11, x01, y01, z01;
  HtofLocalToGlobal(c.plane, c.u0, c.v0, L, x00, y00, z00);
  HtofLocalToGlobal(c.plane, c.u1, c.v0, L, x10, y10, z10);
  HtofLocalToGlobal(c.plane, c.u1, c.v1, L, x11, y11, z11);
  HtofLocalToGlobal(c.plane, c.u0, c.v1, L, x01, y01, z01);
  DrawRect3D(x00, y00, z00, x10, y10, z10, x11, y11, z11, x01, y01, z01, col, width);
}

void
DrawHtofSegXZ(const std::set<Int_t>& highlight)
{
  const auto cells = BuildHtofCells();
  for (const auto& c : cells) {
    const Bool_t hi = (highlight.count(c.seg) > 0);
    const Color_t col = hi ? kRed : kAzure + 2;
    const Width_t win = hi ? 3 : 2;
    // Inner and outer mid-height edges
    for (Double_t L : {kHtofL, kHtofLOuter}) {
      Double_t x0, y0, z0, x1, y1, z1;
      HtofLocalToGlobal(c.plane, c.u0, 0., L, x0, y0, z0);
      HtofLocalToGlobal(c.plane, c.u1, 0., L, x1, y1, z1);
      auto* g = new TGraph(2);
      g->SetPoint(0, z0, x0);
      g->SetPoint(1, z1, x1);
      g->SetLineColor(col);
      g->SetLineWidth(L == kHtofL ? win : 1);
      g->SetLineStyle(L == kHtofL ? 1 : 3);
      g->Draw("L");
    }
    // Thickness connectors at u ends
    for (Double_t u : {c.u0, c.u1}) {
      Double_t xi, yi, zi, xo, yo, zo;
      HtofLocalToGlobal(c.plane, u, 0., kHtofL, xi, yi, zi);
      HtofLocalToGlobal(c.plane, u, 0., kHtofLOuter, xo, yo, zo);
      auto* g = new TGraph(2);
      g->SetPoint(0, zi, xi);
      g->SetPoint(1, zo, xo);
      g->SetLineColor(col);
      g->SetLineWidth(1);
      g->Draw("L");
    }
  }
  // Beam-window hole outline on plane 0 (XZ)
  {
    const Double_t hw = 0.5 * kHtofBeamWinW;
    Double_t x0, y0, z0, x1, y1, z1;
    HtofLocalToGlobal(0, -hw, 0., kHtofL, x0, y0, z0);
    HtofLocalToGlobal(0, hw, 0., kHtofL, x1, y1, z1);
    auto* g = new TGraph(2);
    g->SetPoint(0, z0, x0);
    g->SetPoint(1, z1, x1);
    g->SetLineColor(kGray + 2);
    g->SetLineStyle(2);
    g->SetLineWidth(1);
    g->Draw("L");
  }
}

void
DrawHtofSeg3D(const std::set<Int_t>& highlight)
{
  const auto cells = BuildHtofCells();
  for (const auto& c : cells) {
    const Bool_t hi = (highlight.count(c.seg) > 0);
    const Color_t colIn = hi ? kRed : kAzure + 2;
    const Color_t colOut = hi ? kRed - 7 : kAzure + 1;
    DrawHtofCellFace3D(c, kHtofL, colIn, hi ? 3 : 1);
    DrawHtofCellFace3D(c, kHtofLOuter, colOut, 1);
    // Thickness edges at four corners
    for (Double_t u : {c.u0, c.u1}) {
      for (Double_t v : {c.v0, c.v1}) {
        Double_t xi, yi, zi, xo, yo, zo;
        HtofLocalToGlobal(c.plane, u, v, kHtofL, xi, yi, zi);
        HtofLocalToGlobal(c.plane, u, v, kHtofLOuter, xo, yo, zo);
        DrawLine3DGlobal(xi, yi, zi, xo, yo, zo, colIn, 1);
      }
    }
  }
  // Beam-window hole outline on plane 0 inner face
  {
    const Double_t hw = 0.5 * kHtofBeamWinW;
    const Double_t hh = 0.5 * kHtofBeamWinH;
    Double_t x00, y00, z00, x10, y10, z10, x11, y11, z11, x01, y01, z01;
    HtofLocalToGlobal(0, -hw, -hh, kHtofL, x00, y00, z00);
    HtofLocalToGlobal(0, hw, -hh, kHtofL, x10, y10, z10);
    HtofLocalToGlobal(0, hw, hh, kHtofL, x11, y11, z11);
    HtofLocalToGlobal(0, -hw, hh, kHtofL, x01, y01, z01);
    DrawRect3D(x00, y00, z00, x10, y10, z10, x11, y11, z11, x01, y01, z01, kGray + 2, 1);
  }
}

void
DrawTPCCage3D()
{
  const Int_t n = 8;
  Double_t top[8][3], bot[8][3];
  for (Int_t i = 0; i < n; ++i) {
    const Double_t a = i * TMath::Pi() / 4.0 + TMath::Pi() / 8.0;
    const Double_t xg = kOctagonVertexR * TMath::Cos(a);
    const Double_t zg = kOctagonVertexR * TMath::Sin(a);
    TPCGlobalToDisplay3D(xg, kFrameHalfHeight, zg, top[i][0], top[i][1], top[i][2]);
    TPCGlobalToDisplay3D(xg, -kFrameHalfHeight, zg, bot[i][0], bot[i][1], bot[i][2]);
  }
  auto draw_loop = [&](Double_t pts[][3]) {
    auto* pl = new TPolyLine3D(n + 1);
    pl->SetLineColor(kGray + 1);
    pl->SetLineWidth(2);
    for (Int_t i = 0; i < n; ++i)
      pl->SetPoint(i, pts[i][0], pts[i][1], pts[i][2]);
    pl->SetPoint(n, pts[0][0], pts[0][1], pts[0][2]);
    pl->Draw();
  };
  draw_loop(top);
  draw_loop(bot);
  for (Int_t i = 0; i < n; ++i) {
    auto* pl = new TPolyLine3D(2);
    pl->SetLineColor(kGray + 1);
    pl->SetLineWidth(2);
    pl->SetPoint(0, top[i][0], top[i][1], top[i][2]);
    pl->SetPoint(1, bot[i][0], bot[i][1], bot[i][2]);
    pl->Draw();
  }

  // Bottom layer rings clipped to TPC octagon
  const auto radii = LayerBoundaryRadii();
  const Int_t nSeg = 180;
  const Double_t yBot = -kFrameHalfHeight;
  for (Double_t rBound : radii) {
    std::vector<Double_t> dxs, dys, dzs;
    auto flush = [&]() {
      if (dxs.size() < 2) {
        dxs.clear();
        dys.clear();
        dzs.clear();
        return;
      }
      auto* ring = new TPolyLine3D(static_cast<Int_t>(dxs.size()));
      ring->SetLineColor(kGray + 1);
      ring->SetLineStyle(3);
      ring->SetLineWidth(1);
      for (std::size_t i = 0; i < dxs.size(); ++i)
        ring->SetPoint(static_cast<Int_t>(i), dxs[i], dys[i], dzs[i]);
      ring->Draw();
      dxs.clear();
      dys.clear();
      dzs.clear();
    };
    for (Int_t i = 0; i <= nSeg; ++i) {
      const Double_t th = 2.0 * TMath::Pi() * static_cast<Double_t>(i) / nSeg;
      const Double_t xg = rBound * TMath::Sin(th);
      const Double_t zg = rBound * TMath::Cos(th) + kZTarget;
      if (!InsideTPCOctagon(xg, zg)) {
        flush();
        continue;
      }
      Double_t dx, dy, dz;
      TPCGlobalToDisplay3D(xg, yBot, zg, dx, dy, dz);
      dxs.push_back(dx);
      dys.push_back(dy);
      dzs.push_back(dz);
    }
    flush();
  }
}

void
DrawTargetCylinder3D()
{
  const Int_t n = 64;
  Double_t top[64][3], bot[64][3];
  for (Int_t i = 0; i < n; ++i) {
    const Double_t a = 2. * TMath::Pi() * i / n;
    const Double_t xg = kTargetRadius * TMath::Cos(a);
    const Double_t zg = kTargetRadius * TMath::Sin(a) + kZTarget;
    TPCGlobalToDisplay3D(xg, kTargetHalfY, zg, top[i][0], top[i][1], top[i][2]);
    TPCGlobalToDisplay3D(xg, -kTargetHalfY, zg, bot[i][0], bot[i][1], bot[i][2]);
  }
  auto* plTop = new TPolyLine3D(n + 1);
  auto* plBot = new TPolyLine3D(n + 1);
  plTop->SetLineColor(kBlack);
  plBot->SetLineColor(kBlack);
  plTop->SetLineWidth(2);
  plBot->SetLineWidth(2);
  for (Int_t i = 0; i < n; ++i) {
    plTop->SetPoint(i, top[i][0], top[i][1], top[i][2]);
    plBot->SetPoint(i, bot[i][0], bot[i][1], bot[i][2]);
  }
  plTop->SetPoint(n, top[0][0], top[0][1], top[0][2]);
  plBot->SetPoint(n, bot[0][0], bot[0][1], bot[0][2]);
  plTop->Draw();
  plBot->Draw();
  for (Int_t i = 0; i < 8; ++i) {
    const Int_t idx = i * (n / 8);
    auto* pl = new TPolyLine3D(2);
    pl->SetLineColor(kBlack);
    pl->SetLineWidth(1);
    pl->SetPoint(0, top[idx][0], top[idx][1], top[idx][2]);
    pl->SetPoint(1, bot[idx][0], bot[idx][1], bot[idx][2]);
    pl->Draw();
  }
}

void
DrawBh2ToVtxXZ(Double_t time0_seg, Double_t vx, Double_t vz)
{
  if (!std::isfinite(time0_seg) || !std::isfinite(vx) || !std::isfinite(vz))
    return;
  auto* g = new TGraph(2);
  g->SetPoint(0, kBH2Z, BH2SegX(time0_seg));
  g->SetPoint(1, vz, vx);
  g->SetLineColor(kOrange + 7);
  g->SetLineWidth(1);
  g->SetLineStyle(1);
  g->Draw("L");
}

void
DrawMeasAndExtrap(Double_t cx, Double_t cy, Double_t z0, Double_t rr, Double_t dz,
                  Double_t tmin, Double_t tmax,
                  Double_t vx, Double_t vy, Double_t vz,
                  Double_t ex, Double_t ey, Double_t ez,
                  Bool_t draw_extrap,
                  Color_t col_meas,
                  Bool_t for3d)
{
  // XZ: thinner; 3D keeps a bit thicker for visibility.
  const Width_t wLine = for3d ? 2 : 1;
  const Double_t eps = 1e-4;

  // Magenta first (under meas): continuous helix sample vtx→HTOF, endpoints = DST xyz.
  if (draw_extrap
      && std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz)
      && std::isfinite(ex) && std::isfinite(ey) && std::isfinite(ez)
      && std::isfinite(tmin) && std::isfinite(tmax)) {
    Double_t tv = 0., te = 0., d2v = 0., d2e = 0.;
    const Bool_t okV = FindNearestThetaNearMeas(cx, cy, z0, rr, dz, vx, vy, vz,
                                               tmin, tmax, tv, d2v);
    const Bool_t okE = FindNearestThetaNearMeas(cx, cy, z0, rr, dz, ex, ey, ez,
                                               tmin, tmax, te, d2e);
    if (okV && okE && TMath::Abs(tv - te) > eps) {
      Double_t lo = tv, hi = te;
      if (lo > hi) std::swap(lo, hi);
      std::vector<Double_t> exg, eyg, ezg;
      SampleHelixTPC(cx, cy, z0, rr, dz, lo, hi, 200, exg, eyg, ezg);
      if (exg.size() >= 2) {
        // Force endpoints to DST positions (avoid sample/θ mismatch at markers).
        if (TMath::Abs(tv - lo) < TMath::Abs(tv - hi)) {
          exg.front() = vx; eyg.front() = vy; ezg.front() = vz;
          exg.back() = ex; eyg.back() = ey; ezg.back() = ez;
        } else {
          exg.front() = ex; eyg.front() = ey; ezg.front() = ez;
          exg.back() = vx; eyg.back() = vy; ezg.back() = vz;
        }
        if (for3d) {
          DrawHelix3D(exg, eyg, ezg, kExtrapColor, 1, wLine);
        } else if (auto* g = MakeGraphXZ(exg, ezg)) {
          g->SetLineColor(kExtrapColor);
          g->SetLineWidth(wLine);
          g->SetLineStyle(1);
          g->Draw("L");
        }
      }
    }
  }

  // Measured helix on top: [tmin, tmax]
  std::vector<Double_t> hx, hy, hz;
  SampleHelixTPC(cx, cy, z0, rr, dz, tmin, tmax, 120, hx, hy, hz);
  if (for3d) {
    DrawHelix3D(hx, hy, hz, col_meas, 1, wLine);
  } else if (auto* g = MakeGraphXZ(hx, hz)) {
    g->SetLineColor(col_meas);
    g->SetLineWidth(wLine);
    g->SetLineStyle(1);
    g->Draw("L");
  }
}

Double_t
PlaneResidualAbsS(Double_t x, Double_t /*y*/, Double_t z, Int_t iplane)
{
  if (iplane < 0) return TMath::QuietNaN();
  const Double_t phi = static_cast<Double_t>(iplane) * 0.25 * TMath::Pi();
  const Double_t nx = -TMath::Sin(phi);
  const Double_t nz = -TMath::Cos(phi);
  const Double_t ox = kHtofL * nx;
  const Double_t oz = kHtofL * nz;
  return TMath::Abs((x - ox) * nx + (z - oz) * nz); // n_y=0
}

} // namespace

int
main(int argc, char** argv)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  Long64_t nPages = -1;
  Long64_t startEntry = 0;
  Bool_t matchOnly = kTRUE;
  std::string outPdf;
  std::string htofPath;
  std::string helixPath;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      nPages = std::atoll(argv[++i]);
    } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (std::strcmp(argv[i], "--start") == 0 && i + 1 < argc) {
      startEntry = std::atoll(argv[++i]);
    } else if (std::strcmp(argv[i], "--match-only") == 0) {
      matchOnly = kTRUE;
    } else if (std::strcmp(argv[i], "--no-match-only") == 0) {
      matchOnly = kFALSE;
    } else if (argv[i][0] == '-') {
      usage(argv[0]);
      return 1;
    } else if (htofPath.empty()) {
      htofPath = argv[i];
    } else if (helixPath.empty()) {
      helixPath = argv[i];
    } else {
      usage(argv[0]);
      return 1;
    }
  }

  if (nPages <= 0 || htofPath.empty() || helixPath.empty()) {
    usage(argv[0]);
    return 1;
  }

  TFile fHtof(htofPath.c_str(), "READ");
  TFile fHelix(helixPath.c_str(), "READ");
  if (!fHtof.IsOpen() || fHtof.IsZombie() || !fHelix.IsOpen() || fHelix.IsZombie()) {
    std::cerr << "Error: cannot open input file(s)\n";
    return 1;
  }

  TTree* tHtof = dynamic_cast<TTree*>(fHtof.Get("htof"));
  TTree* tHelix = dynamic_cast<TTree*>(fHelix.Get("tpc"));
  if (!tHtof || !tHelix) {
    std::cerr << "Error: need trees htof and tpc\n";
    return 1;
  }

  UInt_t hx_run = 0, hx_ev = 0;
  Int_t hx_nt = 0;
  std::vector<Double_t>* hx_cx = nullptr;
  std::vector<Double_t>* hx_cy = nullptr;
  std::vector<Double_t>* hx_z0 = nullptr;
  std::vector<Double_t>* hx_r = nullptr;
  std::vector<Double_t>* hx_dz = nullptr;
  std::vector<Double_t>* hx_tmin = nullptr;
  std::vector<Double_t>* hx_tmax = nullptr;
  std::vector<Int_t>* hx_is_beam = nullptr;
  std::vector<Int_t>* hx_is_acc = nullptr;
  std::vector<std::vector<Double_t>>* hx_hitx = nullptr;
  std::vector<std::vector<Double_t>>* hx_hity = nullptr;
  std::vector<std::vector<Double_t>>* hx_hitz = nullptr;

  tHelix->SetBranchAddress("run_number", &hx_run);
  tHelix->SetBranchAddress("event_number", &hx_ev);
  tHelix->SetBranchAddress("ntTpc", &hx_nt);
  tHelix->SetBranchAddress("helix_cx", &hx_cx);
  tHelix->SetBranchAddress("helix_cy", &hx_cy);
  tHelix->SetBranchAddress("helix_z0", &hx_z0);
  tHelix->SetBranchAddress("helix_r", &hx_r);
  tHelix->SetBranchAddress("helix_dz", &hx_dz);
  tHelix->SetBranchAddress("helix_theta_min", &hx_tmin);
  tHelix->SetBranchAddress("helix_theta_max", &hx_tmax);
  if (tHelix->GetBranch("is_beam"))
    tHelix->SetBranchAddress("is_beam", &hx_is_beam);
  if (tHelix->GetBranch("is_accidental"))
    tHelix->SetBranchAddress("is_accidental", &hx_is_acc);
  const Bool_t hasHits = (tHelix->GetBranch("hitpos_x") != nullptr);
  if (hasHits) {
    tHelix->SetBranchAddress("hitpos_x", &hx_hitx);
    tHelix->SetBranchAddress("hitpos_y", &hx_hity);
    tHelix->SetBranchAddress("hitpos_z", &hx_hitz);
  }

  std::map<std::pair<UInt_t, UInt_t>, Long64_t> helixIndex;
  const Long64_t nHelix = tHelix->GetEntries();
  for (Long64_t ie = 0; ie < nHelix; ++ie) {
    tHelix->GetEntry(ie);
    helixIndex[{hx_run, hx_ev}] = ie;
  }

  UInt_t hf_run = 0, hf_ev = 0;
  Int_t hf_nt = 0;
  Double_t hf_time0_seg = TMath::QuietNaN();
  std::vector<Int_t>* hf_match = nullptr;
  std::vector<Int_t>* hf_htof_seg = nullptr;
  std::vector<Int_t>* hf_cl_seg = nullptr;
  std::vector<Double_t>* hf_ex_x = nullptr;
  std::vector<Double_t>* hf_ex_y = nullptr;
  std::vector<Double_t>* hf_ex_z = nullptr;
  std::vector<Int_t>* hf_ex_plane = nullptr;
  std::vector<Double_t>* hf_vtx_x = nullptr;
  std::vector<Double_t>* hf_vtx_y = nullptr;
  std::vector<Double_t>* hf_vtx_z = nullptr;
  std::vector<Double_t>* hf_L_sec = nullptr;
  std::vector<Double_t>* hf_L_beam = nullptr;

  tHtof->SetBranchAddress("run_number", &hf_run);
  tHtof->SetBranchAddress("event_number", &hf_ev);
  tHtof->SetBranchAddress("ntTpc", &hf_nt);
  if (tHtof->GetBranch("time0_seg"))
    tHtof->SetBranchAddress("time0_seg", &hf_time0_seg);
  tHtof->SetBranchAddress("match_ok", &hf_match);
  tHtof->SetBranchAddress("htof_seg", &hf_htof_seg);
  if (tHtof->GetBranch("htof_cl_seg_matched"))
    tHtof->SetBranchAddress("htof_cl_seg_matched", &hf_cl_seg);
  tHtof->SetBranchAddress("extrap_x", &hf_ex_x);
  tHtof->SetBranchAddress("extrap_y", &hf_ex_y);
  tHtof->SetBranchAddress("extrap_z", &hf_ex_z);
  if (tHtof->GetBranch("extrap_plane"))
    tHtof->SetBranchAddress("extrap_plane", &hf_ex_plane);
  tHtof->SetBranchAddress("vtx_x", &hf_vtx_x);
  tHtof->SetBranchAddress("vtx_y", &hf_vtx_y);
  tHtof->SetBranchAddress("vtx_z", &hf_vtx_z);
  if (tHtof->GetBranch("L_sec"))
    tHtof->SetBranchAddress("L_sec", &hf_L_sec);
  if (tHtof->GetBranch("L_beam"))
    tHtof->SetBranchAddress("L_beam", &hf_L_beam);

  const Long64_t nHtof = tHtof->GetEntries();
  if (startEntry < 0 || startEntry >= nHtof) {
    std::cerr << "Error: --start out of range\n";
    return 1;
  }

  if (outPdf.empty()) {
    tHtof->GetEntry(startEntry);
    const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, static_cast<Int_t>(hf_run));
    outPdf = std::string(
      TString::Format("%s/tpc_helix_htof_display_run%05d.pdf", imgDir.Data(), hf_run).Data());
    std::cout << "Default output: " << outPdf << std::endl;
  }

  TCanvas canvas("c_hh", "Helix HTOF", 1600, 800);
  const TString outPath(outPdf.c_str());

  // 3D isotropic TPC zoom (no HTOF)
  const Double_t R = 311.0 * (1.0 + kViewMargin);
  const Double_t vxmin = -R, vxmax = R;
  const Double_t vymin = -R, vymax = R;
  const Double_t vzmin = -R, vzmax = R;

  Long64_t nDrawn = 0;
  for (Long64_t entry = startEntry; entry < nHtof && nDrawn < nPages; ++entry) {
    tHtof->GetEntry(entry);
    const auto itHx = helixIndex.find({hf_run, hf_ev});
    if (itHx == helixIndex.end())
      continue;
    tHelix->GetEntry(itHx->second);
    if (!hf_match || !hx_cx || !hx_tmin || !hx_tmax)
      continue;

    Bool_t anyMatch = kFALSE;
    for (Int_t mflag : *hf_match) {
      if (mflag == 1) {
        anyMatch = kTRUE;
        break;
      }
    }
    if (matchOnly && !anyMatch)
      continue;

    std::set<Int_t> hiSegs;
    for (std::size_t it = 0; it < hf_match->size(); ++it) {
      if ((*hf_match)[it] != 1) continue;
      if (hf_htof_seg && it < hf_htof_seg->size())
        hiSegs.insert((*hf_htof_seg)[it]);
    }

    const Int_t nUse = std::min({hf_nt, hx_nt,
                                 static_cast<Int_t>(hf_match->size()),
                                 static_cast<Int_t>(hx_cx->size()),
                                 static_cast<Int_t>(hx_tmin->size())});

    canvas.Clear();
    canvas.Divide(2, 1);

    // ===== XZ =====
    canvas.cd(1);
    gPad->SetGrid(0, 0);
    gPad->SetFixedAspectRatio(1);
    TH1F* frame = gPad->DrawFrame(-kXZHalf, -kXZHalf, kXZHalf, kXZHalf,
                                  Form("XZ Run %u Ev %u;z_{TPC} [mm];x_{TPC} [mm]",
                                       hf_run, hf_ev));
    frame->SetStats(0);

    DrawLayerRingsXZ();
    DrawOctagonXZ();
    DrawTargetXZ();
    DrawHtofSegXZ(hiSegs);
    DrawBH2XZ(hf_time0_seg);

    // BH2 (time0 seg) → vertex: straight line on XZ only (once per event).
    {
      Double_t vx0 = TMath::QuietNaN(), vz0 = TMath::QuietNaN();
      for (Int_t it = 0; it < nUse; ++it) {
        if (!hf_vtx_x || static_cast<Int_t>(hf_vtx_x->size()) <= it) continue;
        if (!std::isfinite((*hf_vtx_x)[it]) || !std::isfinite((*hf_vtx_z)[it])) continue;
        vx0 = (*hf_vtx_x)[it];
        vz0 = (*hf_vtx_z)[it];
        break;
      }
      DrawBh2ToVtxXZ(hf_time0_seg, vx0, vz0);
    }

    TString info = Form("time0_seg=%.1f ntrk=%d  green=meas  magenta=vtx->HTOF  orange=BH2->vtx",
                        hf_time0_seg, nUse);

    auto* leg = new TLegend(0.12, 0.12, 0.42, 0.12 + 0.035 * TMath::Max(3, nUse * 2 + 1));
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.022);
    {
      auto* gBh2 = new TGraph(2);
      gBh2->SetLineColor(kOrange + 7);
      gBh2->SetLineWidth(1);
      leg->AddEntry(gBh2, "BH2->vtx", "l");
    }

    for (Int_t it = 0; it < nUse; ++it) {
      const Bool_t beam = (hx_is_beam && static_cast<Int_t>(hx_is_beam->size()) > it)
                            ? ((*hx_is_beam)[it] != 0) : kFALSE;
      const Bool_t acc = (hx_is_acc && static_cast<Int_t>(hx_is_acc->size()) > it)
                           ? ((*hx_is_acc)[it] != 0) : kFALSE;
      const Color_t col = TrackColor(beam, acc);
      const Bool_t mok = ((*hf_match)[it] == 1);

      Double_t vx = TMath::QuietNaN(), vy = TMath::QuietNaN(), vz = TMath::QuietNaN();
      Double_t ex = TMath::QuietNaN(), ey = TMath::QuietNaN(), ez = TMath::QuietNaN();
      if (hf_vtx_x && static_cast<Int_t>(hf_vtx_x->size()) > it) {
        vx = (*hf_vtx_x)[it];
        vy = (*hf_vtx_y)[it];
        vz = (*hf_vtx_z)[it];
      }
      if (hf_ex_x && static_cast<Int_t>(hf_ex_x->size()) > it) {
        ex = (*hf_ex_x)[it];
        ey = (*hf_ex_y)[it];
        ez = (*hf_ex_z)[it];
      }

      const Bool_t drawEx = mok && std::isfinite(ex) && std::isfinite(vx);
      DrawMeasAndExtrap((*hx_cx)[it], (*hx_cy)[it], (*hx_z0)[it], (*hx_r)[it], (*hx_dz)[it],
                        (*hx_tmin)[it], (*hx_tmax)[it],
                        vx, vy, vz, ex, ey, ez,
                        drawEx, col, kFALSE);

      {
        auto* gMeas = new TGraph(2);
        gMeas->SetPoint(0, 0., 0.);
        gMeas->SetPoint(1, 1., 1.);
        gMeas->SetLineColor(col);
        gMeas->SetLineWidth(1);
        gMeas->SetLineStyle(1);
        TString lab = Form("tr%d", it);
        if (beam) lab += " beam";
        if (acc) lab += " acc";
        lab += " meas";
        leg->AddEntry(gMeas, lab, "l");
      }
      if (drawEx) {
        Int_t hs = -1;
        if (hf_htof_seg && static_cast<Int_t>(hf_htof_seg->size()) > it)
          hs = (*hf_htof_seg)[it];
        auto* gEx = new TGraph(2);
        gEx->SetPoint(0, 0., 0.);
        gEx->SetPoint(1, 1., 1.);
        gEx->SetLineColor(kExtrapColor);
        gEx->SetLineWidth(1);
        gEx->SetLineStyle(1);
        leg->AddEntry(gEx, Form("tr%d extrap seg%d", it, hs), "l");
      }

      if (hasHits && hx_hitx && hx_hitz && static_cast<Int_t>(hx_hitx->size()) > it) {
        const auto& px = (*hx_hitx)[it];
        const auto& pz = (*hx_hitz)[it];
        const Int_t nh = static_cast<Int_t>(std::min(px.size(), pz.size()));
        for (Int_t ih = 0; ih < nh; ++ih) {
          if (!std::isfinite(px[ih]) || !std::isfinite(pz[ih])) continue;
          auto* mk = new TMarker(pz[ih], px[ih], 20);
          mk->SetMarkerColor(col);
          mk->SetMarkerSize(0.35);
          mk->Draw();
        }
      }
      if (std::isfinite(vx) && std::isfinite(vz)) {
        auto* mk = new TMarker(vz, vx, 29);
        mk->SetMarkerColor(kBlue + 1);
        mk->SetMarkerSize(1.3);
        mk->Draw();
      }
      if (std::isfinite(ex) && std::isfinite(ez) && drawEx) {
        auto* mk = new TMarker(ez, ex, 22);
        mk->SetMarkerColor(kExtrapColor);
        mk->SetMarkerSize(1.2);
        mk->Draw();
        if (hf_htof_seg && static_cast<Int_t>(hf_htof_seg->size()) > it) {
          auto* tx = new TLatex(ez + 8., ex + 8., Form("seg%d", (*hf_htof_seg)[it]));
          tx->SetTextSize(0.025);
          tx->SetTextColor(kExtrapColor);
          tx->Draw();
        }
      }
      if (mok && drawEx) {
        Double_t Ls = TMath::QuietNaN();
        if (hf_L_sec && static_cast<Int_t>(hf_L_sec->size()) > it) Ls = (*hf_L_sec)[it];
        Int_t cl = -1, hs = -1, ipl = -1;
        if (hf_cl_seg && static_cast<Int_t>(hf_cl_seg->size()) > it) cl = TMath::Nint((*hf_cl_seg)[it]);
        if (hf_htof_seg && static_cast<Int_t>(hf_htof_seg->size()) > it) hs = TMath::Nint((*hf_htof_seg)[it]);
        if (hf_ex_plane && static_cast<Int_t>(hf_ex_plane->size()) > it) ipl = (*hf_ex_plane)[it];
        const Double_t abs_s = PlaneResidualAbsS(ex, ey, ez, ipl);
        info += Form(" | tr%d seg%d/cl%d |s|=%.3f Lsec=%.1f", it, hs, cl, abs_s, Ls);
      } else if (mok) {
        info += Form(" | tr%d match but no extrap xyz", it);
      }
    }
    auto* tex = new TLatex(0.12, 0.93, info);
    tex->SetNDC();
    tex->SetTextSize(0.022);
    tex->Draw();
    leg->Draw();

    // ===== 3D (TPC zoom, no BH2, no HTOF) =====
    canvas.cd(2);
    gPad->SetGrid(0, 0);
    auto* view = TView::CreateView(1, 0, 0);
    view->SetRange(vxmin, vymin, vzmin, vxmax, vymax, vzmax);

    DrawTPCCage3D();
    DrawTargetCylinder3D();

    for (Int_t it = 0; it < nUse; ++it) {
      const Bool_t beam = (hx_is_beam && static_cast<Int_t>(hx_is_beam->size()) > it)
                            ? ((*hx_is_beam)[it] != 0) : kFALSE;
      const Bool_t acc = (hx_is_acc && static_cast<Int_t>(hx_is_acc->size()) > it)
                           ? ((*hx_is_acc)[it] != 0) : kFALSE;
      const Color_t col = TrackColor(beam, acc);
      const Bool_t mok = ((*hf_match)[it] == 1);

      Double_t vx = TMath::QuietNaN(), vy = TMath::QuietNaN(), vz = TMath::QuietNaN();
      Double_t ex = TMath::QuietNaN(), ey = TMath::QuietNaN(), ez = TMath::QuietNaN();
      if (hf_vtx_x && static_cast<Int_t>(hf_vtx_x->size()) > it) {
        vx = (*hf_vtx_x)[it];
        vy = (*hf_vtx_y)[it];
        vz = (*hf_vtx_z)[it];
      }
      if (hf_ex_x && static_cast<Int_t>(hf_ex_x->size()) > it) {
        ex = (*hf_ex_x)[it];
        ey = (*hf_ex_y)[it];
        ez = (*hf_ex_z)[it];
      }

      const Bool_t drawEx = mok && std::isfinite(ex) && std::isfinite(vx);
      DrawMeasAndExtrap((*hx_cx)[it], (*hx_cy)[it], (*hx_z0)[it], (*hx_r)[it], (*hx_dz)[it],
                        (*hx_tmin)[it], (*hx_tmax)[it],
                        vx, vy, vz, ex, ey, ez,
                        drawEx, col, kTRUE);

      if (hasHits && hx_hitx && hx_hity && hx_hitz
          && static_cast<Int_t>(hx_hitx->size()) > it) {
        const auto& px = (*hx_hitx)[it];
        const auto& py = (*hx_hity)[it];
        const auto& pz = (*hx_hitz)[it];
        const Int_t nh = static_cast<Int_t>(std::min({px.size(), py.size(), pz.size()}));
        Int_t nGood = 0;
        for (Int_t ih = 0; ih < nh; ++ih) {
          if (std::isfinite(px[ih]) && std::isfinite(py[ih]) && std::isfinite(pz[ih]))
            ++nGood;
        }
        if (nGood > 0) {
          auto* pm = new TPolyMarker3D(nGood);
          pm->SetMarkerColor(col);
          pm->SetMarkerStyle(20);
          pm->SetMarkerSize(0.35);
          Int_t ip = 0;
          for (Int_t ih = 0; ih < nh; ++ih) {
            if (!std::isfinite(px[ih]) || !std::isfinite(py[ih]) || !std::isfinite(pz[ih]))
              continue;
            Double_t dx, dy, dz;
            TPCGlobalToDisplay3D(px[ih], py[ih], pz[ih], dx, dy, dz);
            pm->SetPoint(ip++, dx, dy, dz);
          }
          pm->Draw();
        }
      }
      if (std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz)) {
        Double_t dx, dy, dz;
        TPCGlobalToDisplay3D(vx, vy, vz, dx, dy, dz);
        auto* pm = new TPolyMarker3D(1);
        pm->SetMarkerColor(kBlue + 1);
        pm->SetMarkerStyle(29);
        pm->SetMarkerSize(1.2);
        pm->SetPoint(0, dx, dy, dz);
        pm->Draw();
      }
      if (drawEx && std::isfinite(ex) && std::isfinite(ey) && std::isfinite(ez)) {
        Double_t dx, dy, dz;
        TPCGlobalToDisplay3D(ex, ey, ez, dx, dy, dz);
        auto* pm = new TPolyMarker3D(1);
        pm->SetMarkerColor(kExtrapColor);
        pm->SetMarkerStyle(22);
        pm->SetMarkerSize(1.2);
        pm->SetPoint(0, dx, dy, dz);
        pm->Draw();
      }
    }

    auto* tex3 = new TLatex(0.10, 0.93, "3D; green=meas / magenta=vtx->HTOF (DST ends)");
    tex3->SetNDC();
    tex3->SetTextSize(0.028);
    tex3->Draw();

    auto* leg3 = new TLegend(0.55, 0.12, 0.88, 0.12 + 0.035 * TMath::Max(2, nUse * 2));
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.022);
    for (Int_t it = 0; it < nUse; ++it) {
      const Bool_t beam = (hx_is_beam && static_cast<Int_t>(hx_is_beam->size()) > it)
                            ? ((*hx_is_beam)[it] != 0) : kFALSE;
      const Bool_t acc = (hx_is_acc && static_cast<Int_t>(hx_is_acc->size()) > it)
                           ? ((*hx_is_acc)[it] != 0) : kFALSE;
      const Color_t col = TrackColor(beam, acc);
      const Bool_t mok = ((*hf_match)[it] == 1);
      auto* gMeas = new TGraph(2);
      gMeas->SetLineColor(col);
      gMeas->SetLineWidth(2);
      gMeas->SetLineStyle(1);
      TString lab = Form("tr%d", it);
      if (beam) lab += " beam";
      if (acc) lab += " acc";
      lab += " meas";
      leg3->AddEntry(gMeas, lab, "l");
      if (mok) {
        Int_t hs = -1;
        if (hf_htof_seg && static_cast<Int_t>(hf_htof_seg->size()) > it)
          hs = (*hf_htof_seg)[it];
        auto* gEx = new TGraph(2);
        gEx->SetLineColor(kExtrapColor);
        gEx->SetLineWidth(2);
        gEx->SetLineStyle(1);
        leg3->AddEntry(gEx, Form("tr%d extrap seg%d", it, hs), "l");
      }
    }
    leg3->Draw();

    canvas.cd();
    canvas.Update();
    if (nDrawn == 0)
      canvas.Print(outPath + "(");
    else
      canvas.Print(outPath);
    ++nDrawn;
  }

  if (nDrawn > 0)
    canvas.Print(outPath + ")");
  else
    std::cerr << "Warning: no pages drawn\n";
  std::cout << "Wrote " << nDrawn << " page(s) to " << outPdf << std::endl;
  return 0;
}
