#include "jet_tools/include/plot_helpers.h"
#include "jet_tools/include/math_helpers.h"

#include <cmath>
#include <limits>

#include <TAxis.h>
#include <TH2.h>

namespace jet_tools {

namespace {

double mean_stat(const std::vector<double>& sample, void *) {
  return mean_value(sample);
}

struct StatCtx {
  double (*stat)(std::vector<double>);
};

double stat_adapter(const std::vector<double>& sample, void *ctx_ptr) {
  const auto *ctx = static_cast<const StatCtx *>(ctx_ptr);
  if (!ctx || !ctx->stat)
    return std::numeric_limits<double>::quiet_NaN();
  return ctx->stat(sample);
}

struct QuantileCtx {
  double q = 0.5;
};

double quantile_stat(const std::vector<double>& sample, void *ctx_ptr) {
  const auto *ctx = static_cast<const QuantileCtx *>(ctx_ptr);
  const double q = ctx ? ctx->q : 0.5;
  return quantile_value(sample, q);
}

struct RateCtx {
  double threshold = 0.0;
};

double rate_stat(const std::vector<double>& sample, void *ctx_ptr) {
  const auto *ctx = static_cast<const RateCtx *>(ctx_ptr);
  const double threshold = ctx ? ctx->threshold : 0.0;
  std::size_t n_pass = 0;
  for (double value : sample) {
    if (value > threshold)
      ++n_pass;
  }
  return static_cast<double>(n_pass) / static_cast<double>(sample.size());
}

}  // namespace

double pad_max(double max_value, double fallback, double scale) {
  if (std::isfinite(max_value) && max_value > 0.0)
    return max_value * scale;
  return fallback;
}

AxisRange pad_range_signed(double min_value, double max_value,
                           double fallback_min, double fallback_max,
                           double scale) {
  if (std::isfinite(min_value) && std::isfinite(max_value) &&
      min_value < max_value) {
    return AxisRange{min_value * scale, max_value * scale};
  }
  return AxisRange{fallback_min, fallback_max};
}

void clamp_filled_ymax(TH2 &hist) {
  // Pad y max to avoid clipping filled bins.
  const int nx = hist.GetNbinsX();
  const int ny = hist.GetNbinsY();
  double max_y = -std::numeric_limits<double>::infinity();
  for (int ybin = 1; ybin <= ny; ++ybin) {
    bool has_content = false;
    for (int xbin = 1; xbin <= nx; ++xbin) {
      if (hist.GetBinContent(xbin, ybin) != 0.0) {
        has_content = true;
        break;
      }
    }
    if (has_content) {
      max_y = hist.GetYaxis()->GetBinUpEdge(ybin);
    }
  }
  if (!std::isfinite(max_y))
    return;
  const double y_max_current = hist.GetYaxis()->GetXmax();
  const double y_max_padded = pad_max(max_y, y_max_current);
  hist.GetYaxis()->SetRangeUser(hist.GetYaxis()->GetXmin(), y_max_padded);
}

void clamp_filled_yrange(TH2 &hist) {
  // Pad y min/max to avoid clipping filled bins.
  const int nx = hist.GetNbinsX();
  const int ny = hist.GetNbinsY();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  for (int ybin = 1; ybin <= ny; ++ybin) {
    bool has_content = false;
    for (int xbin = 1; xbin <= nx; ++xbin) {
      if (hist.GetBinContent(xbin, ybin) != 0.0) {
        has_content = true;
        break;
      }
    }
    if (has_content) {
      if (min_y == std::numeric_limits<double>::infinity()) {
        min_y = hist.GetYaxis()->GetBinLowEdge(ybin);
      }
      max_y = hist.GetYaxis()->GetBinUpEdge(ybin);
    }
  }
  if (!std::isfinite(min_y) || !std::isfinite(max_y))
    return;
  const auto range = pad_range_signed(
      min_y, max_y, hist.GetYaxis()->GetXmin(), hist.GetYaxis()->GetXmax());
  hist.GetYaxis()->SetRangeUser(range.min, range.max);
}

double bootstrap_mean_err(const std::vector<double>& v, std::uint32_t seed,
                          int n_samples) {
  return bootstrap_err(v, seed, n_samples, &mean_stat);
}

double bootstrap_err_stat(const std::vector<double>& v, std::uint32_t seed,
                          double (*stat)(std::vector<double>),
                          int n_samples) {
  StatCtx ctx;
  ctx.stat = stat;
  return bootstrap_err(v, seed, n_samples, &stat_adapter, &ctx);
}

double bootstrap_quantile_err(const std::vector<double>& v, std::uint32_t seed,
                              double q, int n_samples) {
  QuantileCtx ctx;
  ctx.q = q;
  return bootstrap_err(v, seed, n_samples, &quantile_stat, &ctx);
}

double bootstrap_rate_err(const std::vector<double>& v, std::uint32_t seed,
                          double threshold, int n_samples) {
  RateCtx ctx;
  ctx.threshold = threshold;
  return bootstrap_err(v, seed, n_samples, &rate_stat, &ctx);
}

}  // namespace jet_tools
