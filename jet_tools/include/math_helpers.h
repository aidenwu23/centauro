#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

#include "jet_tools/include/types.h"

namespace jet_tools {

// Wrap delta phi into (-pi, pi].
double delta_phi(double a, double b);

// Delta-R in eta-phi.
double delta_r(double eta1, double phi1, double eta2, double phi2);

// Minkowski dot product for 4-vectors (metric +,-,-,-).
double minkowski_dot(const P4 &a, const P4 &b);

// Filter finite values from an input vector.
std::vector<double> finite_values(const std::vector<double> &values);

// Mean of a finite-valued vector.
double mean_value(const std::vector<double> &values);

// Median via nth-element on a local copy.
double median_value(std::vector<double> values);

// Quantile with linear interpolation on a sorted local copy.
double quantile_value(std::vector<double> values, double q);

// Sample standard deviation (n-1 denominator).
double sample_stddev_value(const std::vector<double> &values);

// Callback type for bootstrap statistic evaluation.
using BootstrapStatFn = double (*)(const std::vector<double> &sample, void *ctx);

// Generic bootstrap standard error for a statistic callback.
double bootstrap_err(const std::vector<double> &values, std::uint32_t seed,
                     int n_samples, BootstrapStatFn stat_fn, void *ctx = nullptr);

}  // namespace jet_tools
