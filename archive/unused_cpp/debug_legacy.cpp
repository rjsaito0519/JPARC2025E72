// TPC local-track separation debug (pad maps, tpc_event_display_pdf style)
//
//   ./bin/debug                    — gComp / SeparateClustersWithGap gallery
//   ./bin/debug --target-sep       — SeparateTracksAtTarget (beam + scatter merge)
//   ./bin/debug --target-sep -n 8  — number of merged-track pages
//   ./bin/debug --seed 42 [--random]
//   ./bin/debug --dedx
//   ./bin/debug --residual-check [--seed N] [-n N]
//       Verify ResidualVect / ResidualVectXZ are invariant under the
//       "pos - AI(absolute)" vs "AP - AI(relative)" rewrite.
//
// Default PDFs:
//   tpc_gcomp_sort_debug.pdf
//   tpc_target_sep_debug.pdf
// Build: cd myanalysis && cmake -B .build && cmake --build .build

#include "TPCPadHelper.hh"
#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TLine.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

// TPCLocalTrack.cc (SeparateClustersWithGap) file-local constants
constexpr double kMaxGapBtwClusters = 100.;
constexpr int kMaxLayerdiffBtwClusters = 32;

constexpr double kMaxResidualUnderTrack = 35.; // mm, synthetic track: 1 pad/layer on the line
constexpr double kBgPadValue = 0.;
constexpr double kPolyZmaxDefault = 0x1000;
constexpr double kGCompDisplayOffset = 10000.;

constexpr double kDispZmin = -300.;
constexpr double kDispZmax = 300.;
constexpr double kDispXmin = -300.;
constexpr double kDispXmax = 300.;
constexpr double kRankLabelOffsetX = 14.; // mm: sort rank beside pad (lab X)
constexpr int kMinTrackPads = 8;
constexpr double kTinyU0Threshold = 0.05; // |u0| below: division 1/u0 is large (risky)

// SeparateTracksAtTarget (non-vertex branch): beam box + beam/scatter transition
constexpr double kBeamBoxX = 25.;
constexpr double kBeamBoxY = 30.;
constexpr double kBeamScatterGapMm = 30.;
constexpr double kTargetSepV0Min = 0.01;
constexpr double kDimPadValue = 3000.; // COLZ for hits that would be erased

struct PadSample {
  int padid = -1;
  int layer = 0;
  int row = 0;
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double g_comp = 0.;
};

struct TrackParam {
  double x0 = 0.;
  double u0 = 0.;
  const char* tag = "";
};

struct TrackCase {
  double x0 = 0.;
  double u0 = 0.;
  double v0 = 0.;
  std::string tag;
  std::vector<PadSample> pads;
};

struct MergedTrackParam {
  double x0_fit = 0.;
  double u0_fit = 0.;
  double v0_fit = 0.;
  double x0_beam = 0.;
  double u0_beam = 0.;
  double x0_scat = 0.;
  double u0_scat = 0.;
  const char* tag = "";
};

struct TargetSepResult {
  std::vector<int> hit_order;
  std::vector<int> side; // 1 or 2 per pad index
  std::vector<bool> is_beam;
  int flip_at_sort_rank = -1;
  int n_side1 = 0;
  int n_side2 = 0;
  bool keeps_side1 = true;
  bool would_separate = false;
};

// Various slopes / x0 (steep, shallow, sign flip, near-vertical, off-center)
const TrackParam kPresetCases[] = {
    {0., 0.12, "u0 small +"},
    {0., -0.12, "u0 small -"},
    {0., 0.42, "u0 steep +"},
    {0., -0.42, "u0 steep -"},
    {0., 0.45, "u0 max +"},
    {0., -0.45, "u0 max -"},
    {35., 0.30, "x0+ u0+"},
    {-35., -0.30, "x0- u0-"},
    {25., -0.40, "x0+ u0-"},
    {-25., 0.40, "x0- u0+"},
    {0., 0.06, "u0 near-vert +"},
    {15., -0.06, "u0 near-vert -"},
    {40., 0.18, "shallow +"},
    {-40., -0.18, "shallow -"},
    {0., 0.25, "mid +"},
    {0., -0.25, "mid -"},
    // |u0| tiny: division=(1/u0,0,-1) blows up; gComp scale / sort stability check
    {0., 0.03, "u0 tiny +"},
    {0., -0.03, "u0 tiny -"},
    {0., 0.02, "u0 extreme +"},
    {0., -0.02, "u0 extreme -"},
    {0., 0.01, "u0 ~vertical +"},
    {0., -0.01, "u0 ~vertical -"},
    {0., 0.005, "u0 danger +"},
    {0., -0.005, "u0 danger -"},
    {10., 0.015, "x0 off + tiny u0"},
    {-10., -0.015, "x0 off - tiny u0"},
};

// Synthetic merged tracks: upstream beam box + downstream scatter (fit params separate)
const MergedTrackParam kMergedPresets[] = {
    {0., 0.18, 0.05, 0., 0.12, 8., 0.32, "std beam+scatter"},
    {0., 0.18, 0.05, 0., 0.12, -12., -0.28, "scatter u0-"},
    {5., 0.15, 0.08, 2., 0.10, 15., 0.40, "off-axis scatter"},
    {0., 0.20, 0.05, 0., 0.15, 0., 0.45, "steep scatter"},
    {0., 0.18, 0.002, 0., 0.12, 8., 0.32, "v0 too small (no flip)"},
    {0., 0.18, 0.05, 0., 0.12, 8., 0.32, "short scatter arm"},
    {0., 0.12, 0.06, 0., 0.08, 20., 0.22, "shallow scatter"},
    {0., 0.25, 0.05, 0., 0.20, -8., -0.35, "steep merged fit"},
};

void TPCPadTemplateXZ(TH2Poly* h)
{
  Double_t X[5];
  Double_t Y[5];
  for (Int_t l = 0; l < tpc::NumOfLayersTPC; ++l) {
    const Double_t pLength = tpc::padParameter[l][tpc::kLength];
    const Double_t st = (180. - (360. / tpc::padParameter[l][tpc::kNumOfDivision]) *
                                    tpc::padParameter[l][tpc::kNumOfPad] / 2.);
    const Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
    const Double_t dTheta = (360. / tpc::padParameter[l][tpc::kNumOfDivision]) / 180. * TMath::Pi();
    const Double_t cRad = tpc::padParameter[l][tpc::kRadius];
    const Int_t nPad = static_cast<Int_t>(tpc::padParameter[l][tpc::kNumOfPad]);
    for (Int_t j = 0; j < nPad; ++j) {
      const Double_t theta1 = j * dTheta + sTheta;
      const Double_t theta2 = (j + 1) * dTheta + sTheta;
      const Double_t r_in = cRad - (pLength / 2.);
      const Double_t r_out = cRad + (pLength / 2.);

      X[1] = r_out * TMath::Cos(theta1) + tpc::ZTarget;
      X[2] = r_out * TMath::Cos(theta2) + tpc::ZTarget;
      X[3] = r_in * TMath::Cos(theta2) + tpc::ZTarget;
      X[4] = r_in * TMath::Cos(theta1) + tpc::ZTarget;
      X[0] = X[4];

      Y[1] = r_out * TMath::Sin(theta1);
      Y[2] = r_out * TMath::Sin(theta2);
      Y[3] = r_in * TMath::Sin(theta2);
      Y[4] = r_in * TMath::Sin(theta1);
      Y[0] = Y[4];

      h->AddBin(5, X, Y);
    }
  }
  h->SetMaximum(kPolyZmaxDefault);
}

// Same as TPCLocalTrack::SeparateClustersWithGap (lines 1276-1281)
double GCompNormY(const TVector3& pos, double u0)
{
  const TVector3 tgt(0., 0., tpc::ZTarget);
  const TVector3 pos_ = pos - tgt;
  const TVector3 division(1. / u0, 0., -1.);
  const TVector3 norm = pos_.Cross(division);
  return norm.Y();
}

double ResidualXZ(double x, double z, double x0, double u0)
{
  const double z_rel = z - tpc::ZTarget;
  return TMath::Abs(u0 * z_rel - x + x0) / TMath::Hypot(u0, 1.);
}

std::vector<PadSample> CollectPadsDirectlyUnderTrack(double x0, double u0)
{
  std::vector<PadSample> best(static_cast<std::size_t>(tpc::NumOfLayersTPC));
  std::vector<double> best_dist(static_cast<std::size_t>(tpc::NumOfLayersTPC), 1e9);
  std::vector<bool> has(static_cast<std::size_t>(tpc::NumOfLayersTPC), false);

  for (Int_t padid = 0; padid < tpc::NumOfPadTPC; ++padid) {
    if (tpc::IsDead(padid)) continue;
    const TVector3 pos = tpc::getPosition(padid);
    if (!std::isfinite(pos.X()) || !std::isfinite(pos.Z())) continue;

    const Int_t layer = tpc::getLayerID(padid);
    if (layer < 0 || layer >= tpc::NumOfLayersTPC) continue;

    const double dist = ResidualXZ(pos.X(), pos.Z(), x0, u0);
    if (dist > kMaxResidualUnderTrack) continue;
    if (dist >= best_dist[layer]) continue;

    best_dist[layer] = dist;
    has[layer] = true;
    PadSample s;
    s.padid = padid;
    s.layer = layer;
    s.row = tpc::getRowID(padid);
    s.x = pos.X();
    s.y = pos.Y();
    s.z = pos.Z();
    s.g_comp = GCompNormY(pos, u0);
    best[layer] = s;
  }

  std::vector<PadSample> pads;
  pads.reserve(tpc::NumOfLayersTPC);
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    if (has[layer]) pads.push_back(best[layer]);
  }
  return pads;
}

bool IsBeamHitPos(double x, double y, double z)
{
  return TMath::Abs(x) < kBeamBoxX && TMath::Abs(y) < kBeamBoxY && z < tpc::ZTarget;
}

std::vector<PadSample> CollectPadsAlongLineZWindow(double x0, double u0, double z_lo, double z_hi,
                                                   bool beam_box_only)
{
  std::vector<PadSample> best(static_cast<std::size_t>(tpc::NumOfLayersTPC));
  std::vector<double> best_dist(static_cast<std::size_t>(tpc::NumOfLayersTPC), 1e9);
  std::vector<bool> has(static_cast<std::size_t>(tpc::NumOfLayersTPC), false);

  for (Int_t padid = 0; padid < tpc::NumOfPadTPC; ++padid) {
    if (tpc::IsDead(padid)) continue;
    const TVector3 pos = tpc::getPosition(padid);
    if (!std::isfinite(pos.X()) || !std::isfinite(pos.Z())) continue;
    if (pos.Z() < z_lo || pos.Z() > z_hi) continue;
    if (beam_box_only && !IsBeamHitPos(pos.X(), pos.Y(), pos.Z())) continue;

    const Int_t layer = tpc::getLayerID(padid);
    if (layer < 0 || layer >= tpc::NumOfLayersTPC) continue;

    const double dist = ResidualXZ(pos.X(), pos.Z(), x0, u0);
    if (dist > kMaxResidualUnderTrack) continue;
    if (dist >= best_dist[layer]) continue;

    best_dist[layer] = dist;
    has[layer] = true;
    PadSample s;
    s.padid = padid;
    s.layer = layer;
    s.row = tpc::getRowID(padid);
    s.x = pos.X();
    s.y = pos.Y();
    s.z = pos.Z();
    s.g_comp = 0.;
    best[layer] = s;
  }

  std::vector<PadSample> pads;
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    if (has[layer]) pads.push_back(best[layer]);
  }
  return pads;
}

TrackCase MakeMergedBeamScatter(const MergedTrackParam& pr, double z_scat_hi = 280.)
{
  TrackCase tc;
  tc.x0 = pr.x0_fit;
  tc.u0 = pr.u0_fit;
  tc.v0 = pr.v0_fit;
  tc.tag = pr.tag;

  // Leave a physical gap > 30 mm between last beam and first scatter pad (SeparateTracksAtTarget).
  auto beam = CollectPadsAlongLineZWindow(pr.x0_beam, pr.u0_beam, -280., tpc::ZTarget - 35., true);
  auto scat = CollectPadsAlongLineZWindow(pr.x0_scat, pr.u0_scat, tpc::ZTarget + 35., z_scat_hi,
                                          false);

  std::vector<bool> used(static_cast<std::size_t>(tpc::NumOfPadTPC), false);
  tc.pads.reserve(beam.size() + scat.size());
  for (auto p : beam) {
    if (used[p.padid]) continue;
    used[p.padid] = true;
    p.g_comp = GCompNormY(TVector3(p.x, p.y, p.z), tc.u0);
    tc.pads.push_back(p);
  }
  for (auto p : scat) {
    if (used[p.padid]) continue;
    used[p.padid] = true;
    p.g_comp = GCompNormY(TVector3(p.x, p.y, p.z), tc.u0);
    tc.pads.push_back(p);
  }
  return tc;
}

TrackCase MakeCase(double x0, double u0, const char* tag = "")
{
  TrackCase tc;
  tc.x0 = x0;
  tc.u0 = u0;
  tc.v0 = 0.;
  tc.tag = tag ? tag : "";
  tc.pads = CollectPadsDirectlyUnderTrack(x0, u0);
  for (auto& p : tc.pads) {
    p.g_comp = GCompNormY(TVector3(p.x, p.y, p.z), tc.u0);
  }
  return tc;
}

TrackCase MakeRandomCase(TRandom3& rng)
{
  for (int attempt = 0; attempt < 64; ++attempt) {
    const double x0 = rng.Uniform(-40., 40.);
    const double u0 = rng.Uniform(-0.48, 0.48);
    if (TMath::Abs(u0) < 0.05) continue;
    TrackCase tc = MakeCase(x0, u0, "random");
    if (tc.pads.size() >= kMinTrackPads) return tc;
  }
  return MakeCase(0., 0.2, "random-fallback");
}

std::vector<TrackCase> BuildCaseList(int n_cases, unsigned seed, bool random_only)
{
  std::vector<TrackCase> cases;
  cases.reserve(static_cast<std::size_t>(n_cases));

  if (!random_only) {
    for (const auto& pr : kPresetCases) {
      if (static_cast<int>(cases.size()) >= n_cases) break;
      TrackCase tc = MakeCase(pr.x0, pr.u0, pr.tag);
      if (tc.pads.size() < kMinTrackPads) {
        std::cerr << "Skip preset \"" << pr.tag << "\" x0=" << pr.x0 << " u0=" << pr.u0
                  << " (npad=" << tc.pads.size() << ")\n";
        continue;
      }
      cases.push_back(std::move(tc));
    }
  }

  TRandom3 rng(static_cast<ULong_t>(seed));
  int n_try = 0;
  while (static_cast<int>(cases.size()) < n_cases && n_try < n_cases * 30) {
    ++n_try;
    TrackCase tc = MakeRandomCase(rng);
    if (tc.pads.size() < kMinTrackPads) continue;
    cases.push_back(std::move(tc));
  }
  return cases;
}

// TPCLocalTrack::CompareDist + std::sort(hit_order, CompareDist)
void SortHitsByGComp(const std::vector<PadSample>& pads, std::vector<int>& hit_order)
{
  const int n = static_cast<int>(pads.size());
  hit_order.resize(n);
  std::iota(hit_order.begin(), hit_order.end(), 0);
  std::sort(hit_order.begin(), hit_order.end(), [&pads](int a, int b) {
    return pads[a].g_comp < pads[b].g_comp;
  });
}

// SeparateClustersWithGap flip walk (side 1/2 assignment only; no hit erasure)
int CountGapFlips(const std::vector<PadSample>& pads, const std::vector<int>& hit_order)
{
  const int n = static_cast<int>(pads.size());
  if (n == 0) return 0;

  int n_flip = 0;
  bool flip = false;
  int prev_layer = pads[hit_order[0]].layer;
  TVector3 prev_pos(0., 0., 0.);

  for (int i = 0; i < n; ++i) {
    const int id = hit_order[i];
    const TVector3 pos(pads[id].x, pads[id].y, pads[id].z);
    if (i != 0) {
      const TVector3 gap = pos - prev_pos;
      if (TMath::Abs(pads[id].layer - prev_layer) > kMaxLayerdiffBtwClusters ||
          gap.Mag() > kMaxGapBtwClusters) {
        if (!flip) ++n_flip;
        flip = true;
      }
    }
    prev_layer = pads[id].layer;
    prev_pos = pos;
  }
  return n_flip;
}

// TPCLocalTrack::SeparateTracksAtTarget — non-vertex branch (lines 1162-1199)
TargetSepResult AnalyzeSeparateTracksAtTarget(const std::vector<PadSample>& pads, double u0,
                                              double v0)
{
  TargetSepResult out;
  const int n = static_cast<int>(pads.size());
  out.side.assign(n, 1);
  out.is_beam.assign(n, false);
  if (n == 0) return out;

  SortHitsByGComp(pads, out.hit_order);

  bool flag = false;
  TVector3 prev_pos;
  bool prev_is_beam = false;

  for (int i = 0; i < n; ++i) {
    const int id = out.hit_order[i];
    const auto& p = pads[id];
    const bool is_beam = IsBeamHitPos(p.x, p.y, p.z);
    out.is_beam[id] = is_beam;

    const TVector3 pos(p.x, p.y, p.z);
    if (i != 0) {
      const TVector3 gap = pos - prev_pos;
      if (v0 > kTargetSepV0Min && gap.Mag() > kBeamScatterGapMm && prev_is_beam != is_beam) {
        if (!flag) out.flip_at_sort_rank = i;
        flag = true;
      }
    }
    out.side[id] = flag ? 2 : 1;
    prev_is_beam = is_beam;
    prev_pos = pos;
  }

  for (int id = 0; id < n; ++id) {
    if (out.side[id] == 1) ++out.n_side1;
    else ++out.n_side2;
  }
  out.would_separate = (out.n_side1 > 0 && out.n_side2 > 0);
  out.keeps_side1 = out.n_side1 >= out.n_side2;
  return out;
}

void InitPolyBackground(TH2Poly* h, double bg_value)
{
  for (int b = 1; b <= h->GetNumberOfBins(); ++b) {
    h->SetBinContent(b, bg_value);
  }
}

double GCompForDisplay(double g_comp) { return g_comp + kGCompDisplayOffset; }

void SetPolyPadMapZrange(TH2Poly* h, const std::vector<PadSample>& track_pads)
{
  h->SetMinimum(kBgPadValue);
  double zmax = kPolyZmaxDefault;
  if (!track_pads.empty()) {
    double gc_min = track_pads[0].g_comp;
    double gc_max = gc_min;
    for (const auto& p : track_pads) {
      gc_min = TMath::Min(gc_min, p.g_comp);
      gc_max = TMath::Max(gc_max, p.g_comp);
    }
    zmax = TMath::Max(kPolyZmaxDefault, GCompForDisplay(gc_max) * 1.02 + 1.);
  }
  h->SetMaximum(zmax);
}

void DrawSortRankBesidePads(const TrackCase& tc, const std::vector<int>& hit_order)
{
  const int n = static_cast<int>(tc.pads.size());
  const double side = (tc.u0 >= 0.) ? 1. : -1.;

  for (int r = 0; r < n; ++r) {
    const auto& p = tc.pads[hit_order[r]];
    double x_txt = p.x + side * kRankLabelOffsetX;
    if (x_txt > kDispXmax - 5.) x_txt = p.x - side * kRankLabelOffsetX;
    if (x_txt < kDispXmin + 5.) x_txt = p.x + side * kRankLabelOffsetX;

    TLatex lbl;
    lbl.SetTextSize(0.024);
    lbl.SetTextColor(kWhite);
    lbl.SetTextFont(82);
    lbl.DrawLatex(p.z, x_txt, Form("%d", r));
  }
}

void DrawTpcEventDisplay(const TrackCase& tc, const std::vector<int>& hit_order, int case_idx,
                         int n_flip)
{
  const char* tag_str = tc.tag.empty() ? "" : Form(" [%s]", tc.tag.c_str());
  const double inv_u0 = (TMath::Abs(tc.u0) > 1e-12) ? 1. / tc.u0 : 0.;
  const TString u0_str =
    (TMath::Abs(tc.u0) < kTinyU0Threshold) ? Form("%.5f", tc.u0) : Form("%.3f", tc.u0);
  auto* h = new TH2Poly(Form("h_tpc_gcomp_%d", case_idx),
                        Form("TPC X-Z gComp sort (case %d%s, x0=%.1f u0=%s, |1/u0|=%.0f, N=%zu);"
                             "Z [mm];X [mm]",
                             case_idx, tag_str, tc.x0, u0_str.Data(), TMath::Abs(inv_u0),
                             tc.pads.size()),
                        kDispZmin, kDispZmax, kDispXmin, kDispXmax);
  h->SetStats(0);
  TPCPadTemplateXZ(h);
  InitPolyBackground(h, kBgPadValue);

  for (const auto& p : tc.pads) {
    const Int_t bin = h->FindBin(p.z, p.x);
    if (bin > 0) h->SetBinContent(bin, GCompForDisplay(p.g_comp));
  }
  SetPolyPadMapZrange(h, tc.pads);
  h->Draw("COLZ");

  const double x1 = tc.x0 + tc.u0 * (kDispZmin - tpc::ZTarget);
  const double x2 = tc.x0 + tc.u0 * (kDispZmax - tpc::ZTarget);
  TLine track(kDispZmin, x1, kDispZmax, x2);
  track.SetLineColor(kYellow + 1);
  track.SetLineWidth(2);
  track.Draw();

  TLine tgt(tpc::ZTarget, kDispXmin, tpc::ZTarget, kDispXmax);
  tgt.SetLineColor(kRed);
  tgt.SetLineStyle(2);
  tgt.SetLineWidth(2);
  tgt.Draw();

  DrawSortRankBesidePads(tc, hit_order);

  TLatex leg;
  leg.SetNDC(1);
  leg.SetTextSize(0.028);
  leg.DrawLatex(0.05, 0.95,
                Form("COLZ = gComp+%.0f; white # = gComp sort rank (beside pad)", kGCompDisplayOffset));
  leg.DrawLatex(0.05, 0.90, Form("gap flips=%d (SeparateClustersWithGap logic)", n_flip));
  if (TMath::Abs(tc.u0) < kTinyU0Threshold) {
    leg.DrawLatex(0.05, 0.85,
                  Form("tiny |u0|: division uses 1/u0=%.1f (watch gComp scale)", inv_u0));
  }
}

void PrintCaseLog(int case_idx, const TrackCase& tc, const std::vector<int>& hit_order, int n_flip)
{
  const int n = static_cast<int>(tc.pads.size());
  double gc_min = 0., gc_max = 0.;
  if (n > 0) {
    gc_min = gc_max = tc.pads[0].g_comp;
    for (const auto& p : tc.pads) {
      gc_min = TMath::Min(gc_min, p.g_comp);
      gc_max = TMath::Max(gc_max, p.g_comp);
    }
  }
  const double inv_u0 = (TMath::Abs(tc.u0) > 1e-12) ? 1. / tc.u0 : 0.;
  int layer_inversions = 0;
  for (int r = 1; r < n; ++r) {
    if (tc.pads[hit_order[r]].layer < tc.pads[hit_order[r - 1]].layer) ++layer_inversions;
  }

  std::cout << "--- case " << case_idx;
  if (!tc.tag.empty()) std::cout << " [" << tc.tag << "]";
  std::cout << " x0=" << tc.x0 << " u0=" << tc.u0 << " |1/u0|=" << TMath::Abs(inv_u0)
            << " npad=" << n << " gap_flips=" << n_flip << " ---\n";
  std::cout << "  gComp [" << gc_min << ", " << gc_max << "] span=" << (gc_max - gc_min)
            << "  layer inversions along sort=" << layer_inversions << "\n";
  if (TMath::Abs(tc.u0) < kTinyU0Threshold) {
    std::cout << "  ** tiny |u0|: same formula as TPCLocalTrack (1/u0 in division) **\n";
  }

  const int n_print = (n <= 12) ? n : 6;
  for (int r = 0; r < n_print; ++r) {
    const int id = hit_order[r];
    std::cout << "  sort_rank=" << r << " pad=" << tc.pads[id].padid << " layer=" << tc.pads[id].layer
              << " z=" << tc.pads[id].z << " gComp=" << tc.pads[id].g_comp << "\n";
  }
  if (n > 12) {
    std::cout << "  ...\n";
    for (int r = n - 3; r < n; ++r) {
      const int id = hit_order[r];
      std::cout << "  sort_rank=" << r << " pad=" << tc.pads[id].padid
                << " layer=" << tc.pads[id].layer << " z=" << tc.pads[id].z
                << " gComp=" << tc.pads[id].g_comp << "\n";
    }
  }
}

void RunTpcGCompVisualization(const char* pdf_path, int n_cases, unsigned seed, bool random_only)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  const std::vector<TrackCase> cases = BuildCaseList(n_cases, seed, random_only);
  TCanvas canvas("c_tpc_gcomp", "TPC gComp sort debug", 900, 800);

  int n_written = 0;
  for (const auto& tc : cases) {
    std::vector<int> hit_order;
    SortHitsByGComp(tc.pads, hit_order);
    const int n_flip = CountGapFlips(tc.pads, hit_order);

    canvas.Clear();
    DrawTpcEventDisplay(tc, hit_order, n_written, n_flip);

    if (n_written == 0) {
      canvas.Print(TString::Format("%s(", pdf_path));
  } else {
      canvas.Print(pdf_path);
    }
    PrintCaseLog(n_written, tc, hit_order, n_flip);
    ++n_written;
  }

  if (n_written > 0) {
    canvas.Print(TString::Format("%s)", pdf_path));
  }
  std::cout << "Wrote " << pdf_path << " (" << n_written << " page(s), seed=" << seed
            << (random_only ? ", random only" : ", presets+random") << ")\n";
  if (n_written < n_cases) {
    std::cerr << "Warning: requested " << n_cases << " cases, wrote " << n_written << "\n";
  }
}

constexpr double kMassRelTol = 1e-4;

double PdgMassGeV(int pdg_code)
{
  TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
  return particle ? particle->Mass() : -1.;
}

void DrawTargetSepEventDisplay(const TrackCase& tc, const TargetSepResult& sep, int case_idx)
{
  const char* tag_str = tc.tag.empty() ? "" : Form(" [%s]", tc.tag.c_str());
  auto* h = new TH2Poly(Form("h_tgtsep_%d", case_idx),
                        Form("SeparateTracksAtTarget (case %d%s);Z [mm];X [mm]", case_idx,
                             tag_str),
                        kDispZmin, kDispZmax, kDispXmin, kDispXmax);
  h->SetStats(0);
  TPCPadTemplateXZ(h);
  InitPolyBackground(h, kBgPadValue);

  const bool keep1 = sep.keeps_side1;
  for (int id = 0; id < static_cast<int>(tc.pads.size()); ++id) {
    const auto& p = tc.pads[id];
    const bool kept = (sep.side[id] == 1) ? keep1 : !keep1;
    const double val = kept ? GCompForDisplay(p.g_comp) : kDimPadValue;
    const Int_t bin = h->FindBin(p.z, p.x);
    if (bin > 0) h->SetBinContent(bin, val);
  }
  SetPolyPadMapZrange(h, tc.pads);
  h->Draw("COLZ");

  const double x1 = tc.x0 + tc.u0 * (kDispZmin - tpc::ZTarget);
  const double x2 = tc.x0 + tc.u0 * (kDispZmax - tpc::ZTarget);
  TLine fit_line(kDispZmin, x1, kDispZmax, x2);
  fit_line.SetLineColor(kYellow + 1);
  fit_line.SetLineWidth(2);
  fit_line.Draw();

  TLine tgt(tpc::ZTarget, kDispXmin, tpc::ZTarget, kDispXmax);
  tgt.SetLineColor(kRed);
  tgt.SetLineStyle(2);
  tgt.SetLineWidth(2);
  tgt.Draw();

  TLine bx1(-kBeamBoxX, kDispXmin, -kBeamBoxX, kDispXmax);
  TLine bx2(kBeamBoxX, kDispXmin, kBeamBoxX, kDispXmax);
  bx1.SetLineColor(kGreen + 1);
  bx2.SetLineColor(kGreen + 1);
  bx1.SetLineStyle(3);
  bx2.SetLineStyle(3);
  bx1.Draw();
  bx2.Draw();

  const int n = static_cast<int>(tc.pads.size());
  const double side = (tc.u0 >= 0.) ? 1. : -1.;
  for (int r = 0; r < n; ++r) {
    const int id = sep.hit_order[r];
    const auto& p = tc.pads[id];
    double x_txt = p.x + side * kRankLabelOffsetX;
    if (x_txt > kDispXmax - 5.) x_txt = p.x - side * kRankLabelOffsetX;

    const bool kept = (sep.side[id] == 1) ? keep1 : !keep1;
    const bool is_flip = (r == sep.flip_at_sort_rank);
    const int s = sep.side[id];

    int text_color = kGray + 2;
    if (kept) {
      if (sep.is_beam[id]) {
        text_color = kGreen + 2;
      } else if (s == 1) {
        text_color = kCyan + 1;
      } else {
        text_color = kOrange + 7;
      }
    } else if (s == 1) {
      text_color = kGray + 1;
    }
    if (is_flip) text_color = kMagenta + 1;

    TLatex lbl;
    lbl.SetTextSize(0.022);
    lbl.SetTextFont(82);
    lbl.SetTextColor(text_color);
    lbl.DrawLatex(p.z, x_txt, Form("%d%s", r, sep.is_beam[id] ? "B" : ""));
    if (is_flip) {
      TMarker m(p.z, p.x, 29);
      m.SetMarkerColor(kMagenta + 1);
      m.SetMarkerSize(2.2);
      m.Draw();
    }
  }

  TLatex leg;
  leg.SetNDC(1);
  leg.SetTextSize(0.026);
  leg.DrawLatex(0.05, 0.95, Form("fit: x0=%.1f u0=%.3f v0=%.3f (v0>%.2f for beam/scat flip)",
                                 tc.x0, tc.u0, tc.v0, kTargetSepV0Min));
  leg.DrawLatex(0.05, 0.90,
                Form("side1=%d side2=%d  keep %s  separate=%s  flip@sort_rank=%d",
                     sep.n_side1, sep.n_side2, keep1 ? "side1" : "side2",
                     sep.would_separate ? "yes" : "no", sep.flip_at_sort_rank));
  leg.DrawLatex(0.05, 0.85,
                "#color: cyan=side1 kept  orange=side2 kept  green+B=beam  gray=erased  "
                "magenta=flip");
}

void PrintTargetSepLog(int case_idx, const TrackCase& tc, const TargetSepResult& sep)
{
  std::cout << "--- target-sep case " << case_idx;
  if (!tc.tag.empty()) std::cout << " [" << tc.tag << "]";
  std::cout << " fit x0=" << tc.x0 << " u0=" << tc.u0 << " v0=" << tc.v0 << " npad="
            << tc.pads.size() << " n_beam=";
  int n_beam = 0;
  for (bool b : sep.is_beam) if (b) ++n_beam;
  std::cout << n_beam << " ---\n";
  std::cout << "  side1=" << sep.n_side1 << " side2=" << sep.n_side2
            << " keep=" << (sep.keeps_side1 ? "side1" : "side2")
            << " flip_rank=" << sep.flip_at_sort_rank << "\n";
  const int n = static_cast<int>(tc.pads.size());
  for (int r = 0; r < n && r < 8; ++r) {
    const int id = sep.hit_order[r];
    std::cout << "  rank=" << r << " side=" << sep.side[id] << (sep.is_beam[id] ? " beam" : " scat")
              << " pad=" << tc.pads[id].padid << " z=" << tc.pads[id].z << "\n";
  }
  if (n > 8) std::cout << "  ...\n";
}

std::vector<TrackCase> BuildMergedCaseList(int n_cases)
{
  std::vector<TrackCase> cases;
  for (const auto& pr : kMergedPresets) {
    if (static_cast<int>(cases.size()) >= n_cases) break;
    double z_scat_hi = 280.;
    if (std::strstr(pr.tag, "short")) z_scat_hi = 80.;
    TrackCase tc = MakeMergedBeamScatter(pr, z_scat_hi);
    if (tc.pads.size() < kMinTrackPads) {
      std::cerr << "Skip merged \"" << pr.tag << "\" npad=" << tc.pads.size() << "\n";
      continue;
    }
    cases.push_back(std::move(tc));
  }
  return cases;
}

void RunTargetSepVisualization(const char* pdf_path, int n_cases)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  const std::vector<TrackCase> cases = BuildMergedCaseList(n_cases);
  TCanvas canvas("c_tgt_sep", "SeparateTracksAtTarget debug", 900, 800);

  int n_written = 0;
  for (const auto& tc : cases) {
    const TargetSepResult sep = AnalyzeSeparateTracksAtTarget(tc.pads, tc.u0, tc.v0);

    canvas.Clear();
    DrawTargetSepEventDisplay(tc, sep, n_written);

    if (n_written == 0) {
      canvas.Print(TString::Format("%s(", pdf_path));
    } else {
      canvas.Print(pdf_path);
    }
    PrintTargetSepLog(n_written, tc, sep);
    ++n_written;
  }

  if (n_written > 0) {
    canvas.Print(TString::Format("%s)", pdf_path));
  }
  std::cout << "Wrote " << pdf_path << " (" << n_written << " page(s))\n";
}

int RunDedxValidation()
{
  const double mpi = 1000. * PdgMassGeV(kPiMinus);
  const bool ok = TMath::Abs(mpi - 139.57039) / mpi < kMassRelTol;
  std::cout << "dedx smoke: pi mass " << (ok ? "OK" : "NG") << std::endl;
  return ok ? 0 : 1;
}

// --- ResidualVect equivalence check --------------------------------------
// V1: current TPCLocalTrack.cc (pos - AI_absolute)
// V2: rewrite (AP - AI_relative, AI = u * (u . AP))
// Mathematically identical; this verifies numerical equivalence.

TVector3 ResidualVectV1(const double par[4], const TVector3& pos)
{
  const TVector3 x0(par[0], par[1], tpc::ZTarget);
  const TVector3 x1(par[0] + par[2], par[1] + par[3], tpc::ZTarget + 1.);
  const TVector3 u = (x1 - x0).Unit();
  const TVector3 AP = pos - x0;
  const double dist_AX = u.Dot(AP);
  const TVector3 AI_abs(x0.x() + u.x() * dist_AX,
                        x0.y() + u.y() * dist_AX,
                        x0.z() + u.z() * dist_AX);
  return pos - AI_abs;
}

TVector3 ResidualVectV2(const double par[4], const TVector3& pos)
{
  const TVector3 x0(par[0], par[1], tpc::ZTarget);
  const TVector3 x1(par[0] + par[2], par[1] + par[3], tpc::ZTarget + 1.);
  const TVector3 u = (x1 - x0).Unit();
  const TVector3 AP = pos - x0;
  const TVector3 AI = u * u.Dot(AP);
  return AP - AI;
}

TVector3 ResidualVectXZV1(const double par[4], const TVector3& pos)
{
  const TVector3 x0(par[0], pos.y(), tpc::ZTarget);
  const TVector3 x1(par[0] + par[2], pos.y(), tpc::ZTarget + 1.);
  const TVector3 u = (x1 - x0).Unit();
  const TVector3 AP = pos - x0;
  const double dist_AX = u.Dot(AP);
  const TVector3 AI_abs(x0.x() + u.x() * dist_AX,
                        x0.y() + u.y() * dist_AX,
                        x0.z() + u.z() * dist_AX);
  return pos - AI_abs;
}

TVector3 ResidualVectXZV2(const double par[4], const TVector3& pos)
{
  const TVector3 x0(par[0], pos.y(), tpc::ZTarget);
  const TVector3 x1(par[0] + par[2], pos.y(), tpc::ZTarget + 1.);
  const TVector3 u = (x1 - x0).Unit();
  const TVector3 AP = pos - x0;
  const TVector3 AI = u * u.Dot(AP);
  return AP - AI;
}

int RunResidualCheck(unsigned seed, int n_samples)
{
  TRandom3 rng(static_cast<ULong_t>(seed));
  const int n = std::max(n_samples, 1);

  double max_diff_3d = 0., max_diff_xz = 0.;
  double sum_diff_3d = 0., sum_diff_xz = 0.;
  TVector3 worst_pos_3d, worst_pos_xz;
  double worst_par_3d[4] = {0.}, worst_par_xz[4] = {0.};
  double worst_v1_mag_3d = 0., worst_v1_mag_xz = 0.;

  for (int i = 0; i < n; ++i) {
    double par[4] = {rng.Uniform(-50., 50.),
                     rng.Uniform(-50., 50.),
                     rng.Uniform(-0.5, 0.5),
                     rng.Uniform(-0.5, 0.5)};
    const TVector3 pos(rng.Uniform(-300., 300.),
                       rng.Uniform(-300., 300.),
                       rng.Uniform(-300., 300.));

    const TVector3 v1 = ResidualVectV1(par, pos);
    const TVector3 v2 = ResidualVectV2(par, pos);
    const double d3 = (v1 - v2).Mag();
    sum_diff_3d += d3;
    if (d3 > max_diff_3d) {
      max_diff_3d = d3;
      worst_pos_3d = pos;
      std::copy(std::begin(par), std::end(par), std::begin(worst_par_3d));
      worst_v1_mag_3d = v1.Mag();
    }

    const TVector3 v1xz = ResidualVectXZV1(par, pos);
    const TVector3 v2xz = ResidualVectXZV2(par, pos);
    const double dxz = (v1xz - v2xz).Mag();
    sum_diff_xz += dxz;
    if (dxz > max_diff_xz) {
      max_diff_xz = dxz;
      worst_pos_xz = pos;
      std::copy(std::begin(par), std::end(par), std::begin(worst_par_xz));
      worst_v1_mag_xz = v1xz.Mag();
    }
  }

  const double tol = 1e-9;
  const bool ok = max_diff_3d < tol && max_diff_xz < tol;

  std::cout << "--- ResidualVect / ResidualVectXZ rewrite check ---\n"
            << "  samples=" << n << " seed=" << seed << " tol=" << tol << "\n"
            << "  ResidualVect   max|v1-v2|=" << max_diff_3d
            << "  mean=" << (sum_diff_3d / n) << "  (typical |v|=" << worst_v1_mag_3d << ")\n"
            << "    worst par=[" << worst_par_3d[0] << "," << worst_par_3d[1] << ","
            << worst_par_3d[2] << "," << worst_par_3d[3] << "] pos=(" << worst_pos_3d.x()
            << "," << worst_pos_3d.y() << "," << worst_pos_3d.z() << ")\n"
            << "  ResidualVectXZ max|v1-v2|=" << max_diff_xz
            << "  mean=" << (sum_diff_xz / n) << "  (typical |v|=" << worst_v1_mag_xz << ")\n"
            << "    worst par=[" << worst_par_xz[0] << "," << worst_par_xz[1] << ","
            << worst_par_xz[2] << "," << worst_par_xz[3] << "] pos=(" << worst_pos_xz.x()
            << "," << worst_pos_xz.y() << "," << worst_pos_xz.z() << ")\n"
            << "result: " << (ok ? "OK" : "NG") << std::endl;
  return ok ? 0 : 1;
}

} // namespace

int main(int argc, char** argv)
{
  bool run_dedx = false;
  bool run_target_sep = false;
  bool run_residual_check = false;
  bool random_only = false;
  int n_cases = -1; // default set per mode below
  unsigned seed = static_cast<unsigned>(time(nullptr));
  Int_t run_number = 0;
  std::string pdf_path;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--dedx") == 0) {
      run_dedx = true;
    } else if (std::strcmp(argv[i], "--target-sep") == 0) {
      run_target_sep = true;
    } else if (std::strcmp(argv[i], "--residual-check") == 0) {
      run_residual_check = true;
    } else if (std::strcmp(argv[i], "--random") == 0) {
      random_only = true;
    } else if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      n_cases = std::atoi(argv[++i]);
    } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
      seed = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
    } else if (std::strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      run_number = std::atoi(argv[++i]);
    } else if ((std::strcmp(argv[i], "-o") == 0 || std::strcmp(argv[i], "--pdf") == 0) && i + 1 < argc) {
      pdf_path = argv[++i];
    } else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
      std::cout << "Usage:\n"
                   "  debug [-n N] [--seed N] [-r 0] [--random]     gComp / SeparateClustersWithGap\n"
                   "  debug --target-sep [-n N] [-r 0]            SeparateTracksAtTarget (merged)\n"
                   "  debug --residual-check [-n N] [--seed N]    ResidualVect rewrite check\n"
                   "  debug --dedx\n";
      return 0;
    } else {
      std::cerr << "Unknown: " << argv[i] << "\n";
      return 2;
    }
  }

  if (run_dedx) return RunDedxValidation();

  if (run_residual_check) {
    const int n = (n_cases > 0) ? n_cases : 100000;
    return RunResidualCheck(seed, n);
  }

  if (run_target_sep) {
    if (n_cases < 1) {
      n_cases = static_cast<int>(sizeof(kMergedPresets) / sizeof(kMergedPresets[0]));
    }
    if (pdf_path.empty()) {
      const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
      pdf_path = TString::Format("%s/tpc_target_sep_debug.pdf", img_dir.Data()).Data();
    }
    std::cout << "TPCPadHelper ZTarget=" << tpc::ZTarget << " mm  (beam box |x|<" << kBeamBoxX
              << " |y|<" << kBeamBoxY << " z<Ztgt)\n";
    std::cout << "Output: " << pdf_path << "  cases=" << n_cases << "\n";
    RunTargetSepVisualization(pdf_path.c_str(), n_cases);
    return 0;
  }

  if (pdf_path.empty()) {
    const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
    pdf_path = TString::Format("%s/tpc_gcomp_sort_debug.pdf", img_dir.Data()).Data();
  }

  if (n_cases < 1) {
    n_cases = static_cast<int>(sizeof(kPresetCases) / sizeof(kPresetCases[0]));
  }
  std::cout << "TPCPadHelper ZTarget=" << tpc::ZTarget << " mm\n";
  std::cout << "Output: " << pdf_path << "  cases=" << n_cases << "  seed=" << seed << "\n";
  RunTpcGCompVisualization(pdf_path.c_str(), n_cases, seed, random_only);
  return 0;
}
