// m_hit_order / EraseHit / SortHitOrder — mirrors production code
//
// Faithful to:
//   SortHitOrder     src/TPCLocalTrackHelix.cc (CompareTheta on m_hit_t)
//   EraseHit ORIG    ref/src/TPCLocalTrackHelix.cc  (rank-style update)
//   EraseHit NEW     src/TPCLocalTrackHelix.cc      (permutation update)
//
// Model: m_hit_t[i] <-> theta[i], m_hit_order <-> hit_order (same vectors)
//
//   ./bin/debug
//
// Backup: analysis/src/unused/debug_residual_bk.cpp

#include <TRandom3.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr int kDefaultCases = 10;

// --- production: SortHitOrder (src ~2887) ---
// CompareTheta(a,b) => m_hit_t[a] < m_hit_t[b]
void SortHitOrderProd(std::vector<int>& hit_order, const std::vector<double>& hit_t)
{
  hit_order.clear();
  for (std::size_t i = 0; i < hit_t.size(); ++i) {
    hit_order.push_back(static_cast<int>(i));
  }
  std::sort(hit_order.begin(), hit_order.end(),
            [&](int a, int b) { return hit_t[a] < hit_t[b]; });
}

// --- production ORIG: ref EraseHit (~1907) ---
// order = m_hit_order[delete_hit]; erase parallel at delete_hit; decrement > order
struct EraseStepOrig {
  int order_read = -1;
  std::vector<int> after_erase_pos;
  std::vector<int> after_decrement;
};

EraseStepOrig EraseHitOrig(std::vector<double>& hit_t, std::vector<int>& hit_order,
                           int delete_hit)
{
  EraseStepOrig step;
  step.order_read = hit_order[delete_hit];
  hit_t.erase(hit_t.begin() + delete_hit);
  hit_order.erase(hit_order.begin() + delete_hit);
  step.after_erase_pos = hit_order;
  for (int& v : hit_order) {
    if (v > step.order_read) {
      --v;
    }
  }
  step.after_decrement = hit_order;
  return step;
}

// --- production NEW: src EraseHit (~1981) ---
struct EraseStepNew {
  std::vector<int> removed_values;
  std::vector<int> after_reindex;
};

EraseStepNew EraseHitNew(std::vector<double>& hit_t, std::vector<int>& hit_order,
                         int delete_hit)
{
  EraseStepNew step;
  std::vector<int> new_order;
  new_order.reserve(hit_order.size());
  for (int idx : hit_order) {
    if (idx == delete_hit) {
      step.removed_values.push_back(idx);
      continue;
    }
    new_order.push_back(idx > delete_hit ? idx - 1 : idx);
  }
  hit_t.erase(hit_t.begin() + delete_hit);
  hit_order.swap(new_order);
  step.after_reindex = hit_order;
  return step;
}

bool PermMatchesSortHitOrder(const std::vector<int>& hit_order,
                             const std::vector<double>& hit_t)
{
  std::vector<int> expect;
  SortHitOrderProd(expect, hit_t);
  return hit_order == expect;
}

void PrintIntVec(const char* label, const std::vector<int>& v)
{
  std::cout << label << " = [";
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i) std::cout << ", ";
    std::cout << v[i];
  }
  std::cout << "]\n";
}

// storage #N = m_hit_array[N] の添字（配列は θ で並べ替えない）
void PrintHits(const std::vector<double>& hit_t, int mark_delete = -1)
{
  std::cout << "  [ヒット一覧]  #N = storage 番号, θ = そのヒットの角度\n";
  for (std::size_t n = 0; n < hit_t.size(); ++n) {
    std::cout << "    #" << n << "  θ=" << std::fixed << std::setprecision(2)
              << hit_t[n];
    if (static_cast<int>(n) == mark_delete) {
      std::cout << "  ← 削除";
    }
    std::cout << "\n";
  }
}

// P = SortHitOrder の結果。「θ 順で何番目か → どの #N か」
void PrintThetaTrack(const char* title, const std::vector<int>& perm,
                     const std::vector<double>& hit_t)
{
  std::cout << "  [" << title << "]  左から θ 小→大\n    ";
  for (std::size_t rank = 0; rank < perm.size(); ++rank) {
    if (rank) std::cout << " → ";
    const int n = perm[rank];
    std::cout << "#" << n << "(θ=" << std::fixed << std::setprecision(2)
              << hit_t[n] << ")";
  }
  std::cout << "\n";
  std::cout << "    ";
  PrintIntVec("P", perm);
}

struct PathCheck {
  bool matches_rebuild = false;
};

PathCheck RunProductionPath(const char* path_name, bool use_orig_erase,
                            const std::vector<double>& hit_t_in,
                            const std::vector<int>& hit_order_in,
                            int delete_hit)
{
  std::cout << "\n  --- " << path_name << " ---\n";

  std::vector<double> hit_t = hit_t_in;
  std::vector<int> hit_order = hit_order_in;

  PathCheck chk;
  if (use_orig_erase) {
    const EraseStepOrig step = EraseHitOrig(hit_t, hit_order, delete_hit);
    std::cout << "  処理 (ref オリジナル):\n"
              << "    1. order ← P の位置 " << delete_hit << " の値 = "
              << step.order_read << "\n"
              << "       (※ hit #" << delete_hit
              << " の θ 順位ではない)\n"
              << "    2. P から位置 " << delete_hit << " を削除 → ";
    PrintIntVec("P", step.after_erase_pos);
    std::cout << "    3. 残り要素で 値 > " << step.order_read << " なら -1 → ";
    PrintIntVec("P", step.after_decrement);
  } else {
    const EraseStepNew step = EraseHitNew(hit_t, hit_order, delete_hit);
    std::cout << "  処理 (src 現行):\n"
              << "    1. P から値が #" << delete_hit << " の要素を除去\n"
              << "    2. 残りの #番号を再付番 (>#" << delete_hit << " は -1) → ";
    PrintIntVec("P", step.after_reindex);
  }

  PrintHits(hit_t);
  PrintThetaTrack("EraseHit 後の θ 順", hit_order, hit_t);

  std::vector<int> rebuild;
  SortHitOrderProd(rebuild, hit_t);
  PrintThetaTrack("正解 (SortHitOrder し直し)", rebuild, hit_t);

  chk.matches_rebuild = (hit_order == rebuild);
  std::cout << "  判定: " << (chk.matches_rebuild ? "OK" : "NG");
  if (!chk.matches_rebuild) {
    std::cout << "  (上の2行が一致すべき)";
  }
  std::cout << "\n";
  return chk;
}

void RunCase(int case_id, const std::vector<double>& hit_t_in, int delete_hit)
{
  std::cout << "\n========================================\n"
            << " case " << case_id << "\n"
            << "========================================\n";

  std::vector<double> hit_t = hit_t_in;
  std::vector<int> hit_order;
  SortHitOrderProd(hit_order, hit_t);

  std::cout << "  SortHitOrder 済み。P は「θ 順 → #N」の対応 (rank 配列 R ではない)\n";
  PrintHits(hit_t, delete_hit);
  PrintThetaTrack("EraseHit 前の θ 順", hit_order, hit_t);
  std::cout << "  削除: storage #" << delete_hit << "  (θ="
            << std::fixed << std::setprecision(2) << hit_t[delete_hit] << ")\n";

  const PathCheck orig =
      RunProductionPath("ref EraseHit (ORIG)", true, hit_t, hit_order, delete_hit);
  const PathCheck newer =
      RunProductionPath("src EraseHit (NEW)", false, hit_t, hit_order, delete_hit);

  std::cout << "\n  case " << case_id << " result:  ORIG="
            << (orig.matches_rebuild ? "OK" : "NG")
            << "  NEW=" << (newer.matches_rebuild ? "OK" : "NG") << "\n";
}

}  // namespace

int main(int argc, char** argv)
{
  int n_random = kDefaultCases - 1;
  unsigned seed = static_cast<unsigned>(time(nullptr));

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-n" && i + 1 < argc) {
      n_random = std::max(0, std::atoi(argv[++i]) - 1);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
    }
  }

  TRandom3 rng(seed);

  std::cout << "m_hit_order / EraseHit 検証\n"
            << "  #N   = m_hit_array の添字 (storage 番号)\n"
            << "  P    = SortHitOrder 後の並び [θ最小の#N, 次に小さい#N, ...]\n"
            << "  判定 = EraseHit 後の P が SortHitOrder し直しと一致するか\n"
            << "  ORIG = ref オリジナル / NEW = src 現行\n"
            << "  case 0=固定例 + " << n_random << " 乱数  seed=" << seed << "\n";

  // Fixed example from discussion: theta = [-3.1, 3.9, 0.1, 1.2], erase hit 2
  RunCase(0, {-3.1, 3.9, 0.1, 1.2}, 2);

  int n_orig_ok = 0;
  int n_new_ok = 0;
  const int total = 1 + n_random;

  // Re-run fixed case for counting (case 0)
  {
    std::vector<double> t = {-3.1, 3.9, 0.1, 1.2};
    std::vector<int> o;
    SortHitOrderProd(o, t);
    std::vector<double> t1 = t;
    std::vector<int> o1 = o;
    EraseHitOrig(t1, o1, 2);
    if (PermMatchesSortHitOrder(o1, t1)) ++n_orig_ok;
    std::vector<double> t2 = t;
    std::vector<int> o2 = o;
    EraseHitNew(t2, o2, 2);
    if (PermMatchesSortHitOrder(o2, t2)) ++n_new_ok;
  }

  for (int c = 1; c < total; ++c) {
    const int n = rng.Integer(3) + 4;
    std::vector<double> hit_t(n);
    for (int i = 0; i < n; ++i) {
      hit_t[i] = rng.Uniform(-3.5, 3.5);
    }
    const int delete_hit = rng.Integer(n);
    RunCase(c, hit_t, delete_hit);

    std::vector<int> o;
    SortHitOrderProd(o, hit_t);
    std::vector<double> t1 = hit_t;
    std::vector<int> o1 = o;
    EraseHitOrig(t1, o1, delete_hit);
    if (PermMatchesSortHitOrder(o1, t1)) ++n_orig_ok;
    std::vector<double> t2 = hit_t;
    std::vector<int> o2 = o;
    EraseHitNew(t2, o2, delete_hit);
    if (PermMatchesSortHitOrder(o2, t2)) ++n_new_ok;
  }

  std::cout << "\n========================================\n"
            << " summary (" << total << " cases)\n"
            << "   ORIG (ref) matches rebuild SortHitOrder: " << n_orig_ok << "/"
            << total << "\n"
            << "   NEW  (src) matches rebuild SortHitOrder: " << n_new_ok << "/"
            << total << "\n"
            << "========================================\n";

  return (n_new_ok == total) ? 0 : 1;
}
