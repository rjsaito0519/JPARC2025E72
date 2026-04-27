#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

double TruncatedMeanLoop(const std::vector<double>& sorted_values, double fraction)
{
  int n = static_cast<int>(sorted_values.size() * fraction);
  if (n <= 0 || n > static_cast<int>(sorted_values.size())) return 0.;

  double sum = 0.;
  for (int i = 0; i < n; ++i) sum += sorted_values[i];
  return sum / static_cast<double>(n);
}

double TruncatedMeanCumulative(const std::vector<double>& sorted_values, double fraction)
{
  int n = static_cast<int>(sorted_values.size() * fraction);
  if (n <= 0 || n > static_cast<int>(sorted_values.size())) return 0.;

  std::vector<double> cumulative_sum(sorted_values.size() + 1, 0.);
  for (std::size_t i = 0; i < sorted_values.size(); ++i) {
    cumulative_sum[i + 1] = cumulative_sum[i] + sorted_values[i];
  }
  return cumulative_sum[n] / static_cast<double>(n);
}

bool CompareOne(const std::vector<double>& sorted_values, double fraction, const std::string& tag, double eps)
{
  const double loop = TruncatedMeanLoop(sorted_values, fraction);
  const double cumulative = TruncatedMeanCumulative(sorted_values, fraction);
  const double diff = std::abs(loop - cumulative);
  const bool ok = (diff <= eps);

  std::cout << tag
            << "  frac=" << std::fixed << std::setprecision(1) << fraction
            << "  loop=" << std::setprecision(8) << loop
            << "  cumulative=" << cumulative
            << "  |diff|=" << diff
            << "  -> " << (ok ? "OK" : "NG")
            << std::endl;
  return ok;
}

} // namespace

int main()
{
  const double eps = 1e-12;
  bool all_ok = true;

  // 例: dEdx 相当のサンプル（下位から使う前提なのでソートする）
  std::vector<double> dedx = {
    1.23, 0.95, 1.10, 0.78, 1.42, 0.88, 1.05, 1.31, 0.67, 1.56,
    0.73, 0.99, 1.12, 0.84, 1.45, 1.27, 0.91, 1.18, 0.81, 1.36
  };
  std::sort(dedx.begin(), dedx.end());

  std::cout << "=== Compare: legacy loop vs cumulative-sum ===" << std::endl;
  all_ok &= CompareOne(dedx, 1.0, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.9, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.8, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.7, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.6, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.5, "dedx", eps);
  all_ok &= CompareOne(dedx, 0.4, "dedx", eps);

  // 境界ケース確認
  std::vector<double> tiny = {2.0, 1.0};
  std::sort(tiny.begin(), tiny.end());
  std::cout << "\n=== Boundary checks ===" << std::endl;
  all_ok &= CompareOne(tiny, 0.4, "tiny", eps); // n=0 になり得るケース
  all_ok &= CompareOne(tiny, 0.5, "tiny", eps);
  all_ok &= CompareOne(tiny, 1.0, "tiny", eps);

  std::cout << "\nResult: " << (all_ok ? "ALL OK" : "MISMATCH FOUND") << std::endl;

  return all_ok ? 0 : 1;
}

