#ifndef TPC_GAIN_TREE_H
#define TPC_GAIN_TREE_H

// Pad gain factors from the tpc tree of a Helix DST of pion-beam runs (tpc_gain_calib --tree).
//
// Sample
//   track: >= kMinTrackHits hits, kMinMom < p < kMaxMom, and a pion dE/dx window: the track dE/dx
//          (truncated mean) is within +-kPionWindow of the pion Bethe-Bloch expectation and below the
//          midpoint between the pion and proton expectations. The expectation is c * HypTPCdEdx(P10),
//          with the scale c taken from the data (median of dE/dx / Bethe over MIP-like negative tracks).
//   hit:   cluster size 1, |alpha| below kAlphaInner / kAlphaOuter with alpha = theta_diff folded into
//          [-pi/2, pi/2) (theta_diff and theta_diff +- pi give the same path), |y| < kMaxAbsY,
//          dip factor sqrt(1+dz^2) < kMaxDipFactor.
//   value: dE per kOuterPadLength of path (track_cluster_de / pathhit_cor), divided by the pion
//          Bethe-Bloch ratio between the track momentum and kRefMom.
// Factors
//   pad:   the MPV of the value is equalized within each region (inner layers <= kLastInnerLayer, outer
//          layers): factor = (median pad MPV of the region) / (pad MPV). Pads need >= min_hits hits and
//          an MPV in (kMpvMin, kMpvMax). Outer pads keep the current overall level.
//   inner: one more factor for all inner pads so that, on the same track, the truncated mean dE/dx of the
//          inner hits equals that of the outer hits (all cluster sizes, as the track dE/dx uses them).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "TPCPadHelper.hh"

namespace tpc_gain_tree {

constexpr Int_t kMinTrackHits = 10;
constexpr Double_t kMinMom = 0.2;  // GeV/c
constexpr Double_t kMaxMom = 1.6;
constexpr Double_t kPionWindow = 0.25;
constexpr Double_t kAlphaInner = 0.1;  // rad, layers <= kLastInnerLayer (9 mm pads)
constexpr Double_t kAlphaOuter = 0.2;  // rad, outer layers (12.5 mm pads)
constexpr Int_t kLastInnerLayer = 9;
constexpr Double_t kMaxAbsY = 50.0;  // mm
constexpr Double_t kMaxDipFactor = 1.15;
constexpr Double_t kOuterPadLength = 12.5;  // mm
constexpr Double_t kRefMom = 1.0;  // GeV/c
constexpr Double_t kMpvMin = 50.0;
constexpr Double_t kMpvMax = 350.0;
// Tracks for the dE/dx scale c
constexpr Double_t kScaleMinMom = 0.3;
constexpr Double_t kScaleMaxMom = 1.2;
constexpr Double_t kScaleMinDedx = 10.0;
constexpr Double_t kScaleMaxDedx = 30.0;
// Tracks for the inner / outer ratio
constexpr Double_t kRatioMinMom = 0.3;
constexpr Double_t kRatioMaxMom = 0.7;
constexpr Int_t kRatioMinInner = 4;
constexpr Int_t kRatioMinOuter = 8;
constexpr Double_t kTruncatedMeanRatio = 0.8;  // as TPCLocalTrackHelix::GetdEdx
// Value histogram of one pad
constexpr Int_t kHistBins = 100;
constexpr Double_t kHistMax = 800.0;

constexpr Double_t kPionMass = 139.57039;     // MeV/c^2
constexpr Double_t kProtonMass = 938.2720813;  // MeV/c^2

// Kinematics.cc HypTPCdEdx(0 = P10, mass, beta) [MeV/cm] as a function of momentum [GeV/c]
inline Double_t BetheP10(Double_t mass, Double_t p_gev) {
  const Double_t p = 1000.0 * std::fabs(p_gev);
  const Double_t e = std::hypot(mass, p);
  const Double_t beta2 = (p / e) * (p / e);
  const Double_t gamma2 = 1.0 / (1.0 - beta2);
  const Double_t rho = 1e-3 * (0.9 * 1.662 + 0.1 * 0.6672);
  const Double_t zoa = 17.2 / 37.6;
  const Double_t I = 0.9 * 188.0 + 0.1 * 41.7;
  const Double_t a = 0.9 * 0.19714 + 0.1 * 0.09253, k = 0.9 * 2.9618 + 0.1 * 3.6257;
  const Double_t x0 = 0.9 * 1.7635 + 0.1 * 1.6263, x1 = 0.9 * 4.4855 + 0.1 * 3.9716;
  const Double_t C = 0.9 * 11.9480 + 0.1 * 9.5243, d0 = 0.0;
  const Double_t me = 0.5109989461;
  const Double_t X = std::log10(std::sqrt(beta2 * gamma2));
  Double_t delta = 0.0;
  if (X <= x0) delta = d0 * std::pow(10.0, 2.0 * (X - x0));
  else if (X < x1) delta = 2.0 * std::log(10.0) * X - C + a * std::pow(x1 - X, k);
  else delta = 2.0 * std::log(10.0) * X - C;
  const Double_t wmax = 2.0 * me * beta2 * gamma2 /
                        ((me / mass + 1.0) * (me / mass + 1.0) + 2.0 * (me / mass) * (std::sqrt(gamma2) - 1.0));
  return rho * 0.307075 * zoa / beta2 *
         (0.5 * std::log(2.0 * me * beta2 * gamma2 * wmax * 1e12 / (I * I)) - beta2 - 0.5 * delta);
}

inline Double_t Median(std::vector<Double_t> v) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Peak of a Landau-like spectrum: parabola fit to the log of the counts within +-20% of the (3-bin smoothed)
// maximum. Returns <= 0 if it cannot be determined.
inline Double_t PeakMpv(const TH1D* h) {
  const Int_t nb = h->GetNbinsX();
  std::vector<Double_t> sm(nb + 2, 0.0);
  for (Int_t i = 1; i <= nb; ++i) {
    Double_t s = 0.0;
    Int_t n = 0;
    for (Int_t j = i - 1; j <= i + 1; ++j) {
      if (j < 1 || j > nb) continue;
      s += h->GetBinContent(j);
      ++n;
    }
    sm[i] = s / 3.0;
  }
  const Int_t kmax = static_cast<Int_t>(std::max_element(sm.begin() + 1, sm.begin() + nb + 1) - sm.begin());
  const Double_t c0 = h->GetBinCenter(kmax);
  Double_t S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, T0 = 0, T1 = 0, T2 = 0;
  Int_t n = 0;
  for (Int_t i = 1; i <= nb; ++i) {
    const Double_t c = h->GetBinCenter(i);
    const Double_t y = h->GetBinContent(i);
    if (c <= 0.8 * c0 || c >= 1.2 * c0 || y <= 0) continue;
    const Double_t x = c - c0, ly = std::log(y);
    S0 += 1; S1 += x; S2 += x * x; S3 += x * x * x; S4 += x * x * x * x;
    T0 += ly; T1 += x * ly; T2 += x * x * ly;
    ++n;
  }
  if (n < 4) return c0;
  // Solve the 3x3 normal equations for ly = A + B x + C x^2.
  const Double_t M[3][3] = {{S0, S1, S2}, {S1, S2, S3}, {S2, S3, S4}};
  const Double_t det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                       M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  if (std::fabs(det) < 1e-300) return c0;
  const Double_t detB = M[0][0] * (T1 * M[2][2] - M[1][2] * T2) - T0 * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                        M[0][2] * (M[1][0] * T2 - T1 * M[2][0]);
  const Double_t detC = M[0][0] * (M[1][1] * T2 - T1 * M[2][1]) - M[0][1] * (M[1][0] * T2 - T1 * M[2][0]) +
                        T0 * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  const Double_t B = detB / det, C = detC / det;
  if (!(C < 0)) return c0;
  return c0 - B / (2.0 * C);
}

inline Double_t TruncatedMean(std::vector<Double_t> v) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const Int_t nt = std::max(1, static_cast<Int_t>(v.size() * kTruncatedMeanRatio));
  Double_t s = 0.0;
  for (Int_t i = 0; i < nt; ++i) s += v[i];
  return s / nt;
}

inline Bool_t IsInner(Int_t layer) { return layer <= kLastInnerLayer; }

using PadKey = std::pair<Int_t, Int_t>;

struct Result {
  Double_t dedx_scale = 0.0;          // c: track dE/dx / Bethe-Bloch(pion)
  Long64_t n_tracks_read = 0;
  Long64_t n_tracks_pass = 0;
  Long64_t n_hits_used = 0;
  std::map<PadKey, TH1D*> pad_hist;   // value histograms of the selected hits (owned)
  std::map<PadKey, Double_t> pad_mpv; // fitted pads only
  Double_t median_mpv_inner = 0.0;
  Double_t median_mpv_outer = 0.0;
  Double_t inner_scale = 1.0;         // extra factor on all inner pads
  Double_t ratio_before = 0.0;        // median inner/outer truncated mean, current gain
  Double_t ratio_pads = 0.0;          // ... with the pad factors only
  Double_t ratio_after = 0.0;         // ... with pad factors and inner_scale
  Long64_t n_ratio_tracks = 0;
  std::map<PadKey, Double_t> factor;  // final factor of every pad seen in the sample
  std::vector<Double_t> ratio_before_v, ratio_after_v;

  ~Result() {
    for (auto& kv : pad_hist) delete kv.second;
  }
};

// Returns kFALSE if the tree or a branch is missing.
inline Bool_t Run(TFile* f, Int_t min_hits, Result& res) {
  TTree* tree = dynamic_cast<TTree*>(f->Get("tpc"));
  if (!tree) {
    std::cerr << "Error: tree 'tpc' not found in " << f->GetName() << std::endl;
    return kFALSE;
  }
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<Double_t>> rv_mom(reader, "mom0");
  TTreeReaderValue<std::vector<Int_t>> rv_charge(reader, "charge");
  TTreeReaderValue<std::vector<Double_t>> rv_dedx(reader, "dEdx");
  TTreeReaderValue<std::vector<Double_t>> rv_dzf(reader, "dz_factor");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_de(reader, "track_cluster_de");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_cs(reader, "track_cluster_size");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_td(reader, "theta_diff");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_ly(reader, "hitlayer");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_row(reader, "track_cluster_row_center");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_y(reader, "track_cluster_y_center");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_path(reader, "pathhit_cor");

  struct Trk {
    Float_t p, dedx;
    Int_t q;
  };
  std::vector<Trk> trks;
  // selected hits for the pad MPVs
  std::vector<Int_t> hit_trk, hit_key;
  std::vector<Float_t> hit_val;
  // all hits of the tracks used for the inner / outer ratio
  std::vector<Int_t> rt_trk, rt_off;
  std::vector<Int_t> rt_key;
  std::vector<Float_t> rt_x;  // dE per mm with the current gain

  const Double_t bethe_ref = BetheP10(kPionMass, kRefMom);
  const Long64_t n_entries = tree->GetEntries();
  Long64_t iev = 0;
  while (reader.Next()) {
    if (++iev % 100000 == 0) std::cout << "  tree: " << iev << " / " << n_entries << " events" << std::endl;
    const Int_t nt = static_cast<Int_t>(rv_mom->size());
    for (Int_t it = 0; it < nt; ++it) {
      ++res.n_tracks_read;
      if (it >= static_cast<Int_t>(rv_de->size())) continue;
      const auto& de = (*rv_de)[it];
      const Int_t nh = static_cast<Int_t>(de.size());
      const Double_t p = std::fabs((*rv_mom)[it]);
      if (nh < kMinTrackHits || !(p > kMinMom && p < kMaxMom)) continue;
      const Int_t itrk = static_cast<Int_t>(trks.size());
      trks.push_back({static_cast<Float_t>(p), static_cast<Float_t>((*rv_dedx)[it]), (*rv_charge)[it]});
      const Double_t dzf = (*rv_dzf)[it];
      const Double_t norm = kOuterPadLength * bethe_ref / BetheP10(kPionMass, p);
      const auto& cs = (*rv_cs)[it];
      const auto& td = (*rv_td)[it];
      const auto& ly = (*rv_ly)[it];
      const auto& row = (*rv_row)[it];
      const auto& y = (*rv_y)[it];
      const auto& path = (*rv_path)[it];
      Int_t n_in = 0, n_out = 0;
      for (Int_t ih = 0; ih < nh; ++ih) {
        if (!(path[ih] > 0) || !(de[ih] > 0)) continue;
        const Int_t L = static_cast<Int_t>(ly[ih]);
        (IsInner(L) ? n_in : n_out)++;
        if (static_cast<Int_t>(cs[ih]) != 1) continue;
        const Double_t alpha = std::fmod(std::fmod(td[ih] + 0.5 * TMath::Pi(), TMath::Pi()) + TMath::Pi(), TMath::Pi()) -
                               0.5 * TMath::Pi();
        if (std::fabs(alpha) >= (IsInner(L) ? kAlphaInner : kAlphaOuter)) continue;
        if (std::fabs(y[ih]) >= kMaxAbsY || dzf >= kMaxDipFactor) continue;
        hit_trk.push_back(itrk);
        hit_key.push_back(L * 1000 + static_cast<Int_t>(row[ih]));
        hit_val.push_back(static_cast<Float_t>(de[ih] / path[ih] * norm));
      }
      if (p > kRatioMinMom && p < kRatioMaxMom && n_in >= kRatioMinInner && n_out >= kRatioMinOuter) {
        rt_trk.push_back(itrk);
        rt_off.push_back(static_cast<Int_t>(rt_x.size()));
        for (Int_t ih = 0; ih < nh; ++ih) {
          if (!(path[ih] > 0) || !(de[ih] > 0)) continue;
          rt_key.push_back(static_cast<Int_t>(ly[ih]) * 1000 + static_cast<Int_t>(row[ih]));
          rt_x.push_back(static_cast<Float_t>(de[ih] / path[ih]));
        }
      }
    }
  }
  rt_off.push_back(static_cast<Int_t>(rt_x.size()));

  // dE/dx scale relative to the pion Bethe-Bloch
  std::vector<Double_t> sc;
  for (const auto& t : trks) {
    if (t.q < 0 && t.p > kScaleMinMom && t.p < kScaleMaxMom && t.dedx > kScaleMinDedx && t.dedx < kScaleMaxDedx)
      sc.push_back(t.dedx / BetheP10(kPionMass, t.p));
  }
  res.dedx_scale = Median(sc);
  std::vector<Bool_t> pass(trks.size(), kFALSE);
  for (size_t i = 0; i < trks.size(); ++i) {
    const Double_t e_pi = res.dedx_scale * BetheP10(kPionMass, trks[i].p);
    const Double_t e_p = res.dedx_scale * BetheP10(kProtonMass, trks[i].p);
    pass[i] = std::fabs(trks[i].dedx / e_pi - 1.0) < kPionWindow && trks[i].dedx < 0.5 * (e_pi + e_p);
    if (pass[i]) ++res.n_tracks_pass;
  }

  // pad MPVs
  for (size_t i = 0; i < hit_val.size(); ++i) {
    if (!pass[hit_trk[i]]) continue;
    const PadKey key{hit_key[i] / 1000, hit_key[i] % 1000};
    TH1D*& h = res.pad_hist[key];
    if (!h) {
      h = new TH1D(Form("hGainTree_L%02d_R%03d", key.first, key.second),
                   Form("L%02d R%03d;dE per %.1f mm at %.1f GeV/c (pion);hits", key.first, key.second, kOuterPadLength,
                        kRefMom),
                   kHistBins, 0.0, kHistMax);
      h->SetDirectory(nullptr);
    }
    h->Fill(hit_val[i]);
    ++res.n_hits_used;
  }
  std::vector<Double_t> mpv_in, mpv_out;
  for (const auto& kv : res.pad_hist) {
    if (kv.second->GetEntries() < min_hits) continue;
    const Double_t m = PeakMpv(kv.second);
    if (!(m > kMpvMin && m < kMpvMax)) continue;
    res.pad_mpv[kv.first] = m;
    (IsInner(kv.first.first) ? mpv_in : mpv_out).push_back(m);
  }
  res.median_mpv_inner = Median(mpv_in);
  res.median_mpv_outer = Median(mpv_out);
  std::map<PadKey, Double_t> pad_factor;
  for (const auto& kv : res.pad_mpv) {
    const Double_t ref = IsInner(kv.first.first) ? res.median_mpv_inner : res.median_mpv_outer;
    pad_factor[kv.first] = ref / kv.second;
  }

  // inner / outer truncated mean on the same track
  auto ratios = [&](Double_t inner_scale, Bool_t use_pad) {
    std::vector<Double_t> r;
    for (size_t k = 0; k + 1 < rt_off.size(); ++k) {
      if (!pass[rt_trk[k]]) continue;
      std::vector<Double_t> xin, xout;
      for (Int_t j = rt_off[k]; j < rt_off[k + 1]; ++j) {
        const PadKey key{rt_key[j] / 1000, rt_key[j] % 1000};
        Double_t f = 1.0;
        if (use_pad) {
          const auto it = pad_factor.find(key);
          if (it != pad_factor.end()) f = it->second;
        }
        if (IsInner(key.first)) xin.push_back(rt_x[j] * f * inner_scale);
        else xout.push_back(rt_x[j] * f);
      }
      if (static_cast<Int_t>(xin.size()) < kRatioMinInner || static_cast<Int_t>(xout.size()) < kRatioMinOuter) continue;
      r.push_back(TruncatedMean(xin) / TruncatedMean(xout));
    }
    return r;
  };
  res.ratio_before_v = ratios(1.0, kFALSE);
  res.ratio_before = Median(res.ratio_before_v);
  res.ratio_pads = Median(ratios(1.0, kTRUE));
  res.inner_scale = (res.ratio_pads > 0) ? 1.0 / res.ratio_pads : 1.0;
  res.ratio_after_v = ratios(res.inner_scale, kTRUE);
  res.ratio_after = Median(res.ratio_after_v);
  res.n_ratio_tracks = static_cast<Long64_t>(res.ratio_after_v.size());

  for (const auto& kv : res.pad_hist) {
    const auto it = pad_factor.find(kv.first);
    Double_t f = (it != pad_factor.end()) ? it->second : 1.0;
    if (IsInner(kv.first.first)) f *= res.inner_scale;
    res.factor[kv.first] = f;
  }
  return kTRUE;
}

}  // namespace tpc_gain_tree

#endif  // TPC_GAIN_TREE_H
