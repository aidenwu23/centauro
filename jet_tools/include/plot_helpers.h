#pragma once

// Shared plot and statistics helpers for QA scripts.

#include <cstdint>
#include <vector>

class TH2;

namespace jet_tools {

struct AxisRange {
  double min = 0.0;
  double max = 0.0;
};

double pad_max(double max_value, double fallback, double scale = 1.1);

AxisRange pad_range_signed(double min_value, double max_value,
                           double fallback_min, double fallback_max,
                           double scale = 1.1);

void clamp_filled_ymax(TH2 &hist);
void clamp_filled_yrange(TH2 &hist);

double bootstrap_mean_err(const std::vector<double>& v, std::uint32_t seed,
                          int n_samples = 300);

double bootstrap_err_stat(const std::vector<double>& v, std::uint32_t seed,
                          double (*stat)(std::vector<double>),
                          int n_samples = 300);

double bootstrap_quantile_err(const std::vector<double>& v, std::uint32_t seed,
                              double q, int n_samples = 300);

double bootstrap_rate_err(const std::vector<double>& v, std::uint32_t seed,
                          double threshold, int n_samples = 300);

}  // namespace jet_tools
