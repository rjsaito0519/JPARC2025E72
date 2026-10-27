#ifndef TPC_PAD_FILL_H
#define TPC_PAD_FILL_H

// Fill of pads that a calibration fit did not update (TPCPRM per-pad parameters).
//
// Sources are fitted pads that agree with their (layer, AsAd) group: the group median and
// MAD are computed from all fitted pads of the group (at least kMinGroup pads), and the
// caller's predicate decides whether a value is consistent with them. Center-frame pads
// are never sources or targets.
//
// A target pad (not fitted, not on the center frame) gets
//   1. the mean of rows row-1 and row+1 if both are sources in the same AsAd as the pad;
//   2. otherwise the median of the sources of its (layer, AsAd) group (>= kMinGroup);
//   3. otherwise nothing (the caller keeps the previous value).
// Neighbours across an AsAd boundary are not used: offsets jump there.

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <Rtypes.h>

#include "TPCPadHelper.hh"

namespace tpc_fill {

using PadKey = std::pair<Int_t, Int_t>;

enum class FillKind { kNeighbor, kAsadMedian };

struct GroupStat {
  Double_t median = 0.0;
  Double_t mad = 0.0;  // median absolute deviation (not scaled)
  Int_t n = 0;
};

struct FillResult {
  std::map<PadKey, Double_t> value;
  std::map<PadKey, FillKind> kind;
  Int_t n_neighbor = 0;
  Int_t n_asad_median = 0;
  Int_t n_left = 0;  // targets that could not be filled
};

constexpr Int_t kMinGroup = 3;

inline Double_t Median(std::vector<Double_t> v) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const size_t mid = v.size() / 2;
  return (v.size() % 2 == 0) ? 0.5 * (v[mid - 1] + v[mid]) : v[mid];
}

// values  : current parameter of every pad that may be a source or a target
// fitted  : pads updated by the fit in this pass
// excluded: fitted pads that must not be used as sources (e.g. forced outliers)
// consistent(value, group): source vetting against the (layer, AsAd) group of fitted pads
inline FillResult FillUnfittedPads(const std::map<PadKey, Double_t>& values, const std::set<PadKey>& fitted,
                                   const std::set<PadKey>& excluded,
                                   const std::function<Bool_t(Double_t, const GroupStat&)>& consistent) {
  auto group_of = [](const PadKey& k) { return std::make_pair(k.first, tpc::GetASADId(k.first, k.second)); };

  std::map<PadKey, std::vector<Double_t>> group_vals;
  for (const auto& kv : values) {
    const PadKey& k = kv.first;
    if (!fitted.count(k) || excluded.count(k) || tpc::IsPadOnCenterFrame(k.first, k.second)) continue;
    if (!std::isfinite(kv.second)) continue;
    group_vals[group_of(k)].push_back(kv.second);
  }
  std::map<PadKey, GroupStat> group_stat;
  for (const auto& kv : group_vals) {
    if (static_cast<Int_t>(kv.second.size()) < kMinGroup) continue;
    GroupStat st;
    st.n = static_cast<Int_t>(kv.second.size());
    st.median = Median(kv.second);
    std::vector<Double_t> dev;
    dev.reserve(kv.second.size());
    for (Double_t v : kv.second) dev.push_back(std::fabs(v - st.median));
    st.mad = Median(dev);
    group_stat[kv.first] = st;
  }

  std::map<PadKey, Double_t> source;
  std::map<PadKey, std::vector<Double_t>> group_sources;
  for (const auto& kv : values) {
    const PadKey& k = kv.first;
    if (!fitted.count(k) || excluded.count(k) || tpc::IsPadOnCenterFrame(k.first, k.second)) continue;
    if (!std::isfinite(kv.second)) continue;
    const auto it = group_stat.find(group_of(k));
    if (it == group_stat.end() || !consistent(kv.second, it->second)) continue;
    source[k] = kv.second;
    group_sources[group_of(k)].push_back(kv.second);
  }

  FillResult res;
  for (const auto& kv : values) {
    const PadKey& k = kv.first;
    if (fitted.count(k) || tpc::IsPadOnCenterFrame(k.first, k.second)) continue;
    const Int_t asad = tpc::GetASADId(k.first, k.second);
    const PadKey lo{k.first, k.second - 1};
    const PadKey hi{k.first, k.second + 1};
    const auto it_lo = source.find(lo);
    const auto it_hi = source.find(hi);
    if (it_lo != source.end() && it_hi != source.end() && tpc::GetASADId(lo.first, lo.second) == asad &&
        tpc::GetASADId(hi.first, hi.second) == asad) {
      res.value[k] = 0.5 * (it_lo->second + it_hi->second);
      res.kind[k] = FillKind::kNeighbor;
      ++res.n_neighbor;
      continue;
    }
    const auto it_g = group_sources.find(group_of(k));
    if (it_g != group_sources.end() && static_cast<Int_t>(it_g->second.size()) >= kMinGroup) {
      res.value[k] = Median(it_g->second);
      res.kind[k] = FillKind::kAsadMedian;
      ++res.n_asad_median;
      continue;
    }
    ++res.n_left;
  }
  return res;
}

}  // namespace tpc_fill

#endif  // TPC_PAD_FILL_H
