#include "jet_tools/include/math_helpers.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <random>

namespace jet_tools {

double delta_phi(double a, double b) {
  constexpr double pi = 3.14159265358979323846;
  double d = std::fmod(a - b, 2 * pi);
  if (d > pi)
    d -= 2 * pi;
  if (d <= -pi)
    d += 2 * pi;
  return d;
}

double delta_r(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = delta_phi(phi1, phi2);
  const double deta = eta1 - eta2;
  return std::sqrt(dphi * dphi + deta * deta);
}

double minkowski_dot(const P4 &a, const P4 &b) {
  return a.E() * b.E() - (a.Px() * b.Px() + a.Py() * b.Py() + a.Pz() * b.Pz());
}

std::vector<double> finite_values(const std::vector<double> &values) {
  std::vector<double> out;
  out.reserve(values.size());
  for (double value : values) {
    if (std::isfinite(value))
      out.push_back(value);
  }
  return out;
}

double mean_value(const std::vector<double> &values) {
  if (values.empty())
    return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0;
  for (double value : values)
    sum += value;
  return sum / static_cast<double>(values.size());
}

double median_value(std::vector<double> values) {
  if (values.empty())
    return std::numeric_limits<double>::quiet_NaN();
  const std::size_t n = values.size();
  const std::size_t mid = n / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  double med = values[mid];
  if (n % 2 == 0) {
    std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
    med = 0.5 * (med + values[mid - 1]);
  }
  return med;
}

double quantile_value(std::vector<double> values, double q) {
  if (values.empty())
    return std::numeric_limits<double>::quiet_NaN();
  if (q <= 0.0)
    return *std::min_element(values.begin(), values.end());
  if (q >= 1.0)
    return *std::max_element(values.begin(), values.end());
  std::sort(values.begin(), values.end());
  const double index = q * static_cast<double>(values.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(index));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(index));
  if (lo == hi)
    return values[lo];
  const double t = index - static_cast<double>(lo);
  return (1.0 - t) * values[lo] + t * values[hi];
}

double sample_stddev_value(const std::vector<double> &values) {
  if (values.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();
  const double mean = mean_value(values);
  double var = 0.0;
  for (double value : values) {
    const double diff = value - mean;
    var += diff * diff;
  }
  var /= static_cast<double>(values.size() - 1);
  return std::sqrt(var);
}

double bootstrap_err(const std::vector<double> &values, std::uint32_t seed,
                     int n_samples, BootstrapStatFn stat_fn, void *ctx) {
  if (!stat_fn || n_samples <= 1)
    return std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> finite = finite_values(values);
  if (finite.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();

  std::mt19937 rng(seed);
  std::uniform_int_distribution<std::size_t> dist(0, finite.size() - 1);
  std::vector<double> sample(finite.size());
  std::vector<double> stats;
  stats.reserve(static_cast<std::size_t>(n_samples));

  for (int sample_id = 0; sample_id < n_samples; ++sample_id) {
    for (std::size_t index = 0; index < finite.size(); ++index)
      sample[index] = finite[dist(rng)];
    const double stat = stat_fn(sample, ctx);
    if (std::isfinite(stat))
      stats.push_back(stat);
  }
  return sample_stddev_value(stats);
}

}  // namespace jet_tools
