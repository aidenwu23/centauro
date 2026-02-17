#pragma once

#include <vector>

#include "jet_tools/include/config.h"
#include "jet_tools/include/types.h"

namespace jet_tools {

struct CalibMaps {
  std::vector<double> crel[2];
  std::vector<double> cabs[2];
  std::vector<double> cdis_pt_y[2];
  std::vector<double> cabs_final[2];
  std::vector<double> eta_edges;
  std::vector<double> pT_edges;
  std::vector<double> pT_dis_edges;
  std::vector<double> y_edges;
};

// Step per bin, frozen outside.
double interpolate_bins(const std::vector<double> &edges,
                        const std::vector<double> &vals, double x);

CalibMaps derive_calibration_maps(const CalibConfig &cfg,
                                  const std::vector<MatchRecord> &records);

}  // namespace jet_tools
