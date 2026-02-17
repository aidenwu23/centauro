#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace jet_tools {

// CLI-facing options (parse_args fills these).
struct Args {
  std::string input{};
  std::string output{};

  double train_frac = 0.5;
  bool renorm = false;

  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  std::optional<std::size_t> N_min{};
  std::optional<double> max_abs_eta{};

  std::vector<double> eta_bins_opt{};
  std::vector<double> pT_bins_opt{};

  // Reco hygiene / matching.
  std::optional<int> min_constituents_reco{};
  std::optional<double> min_match_pT{};
  std::optional<double> max_match_dR{};
  std::optional<double> min_match_pt_ratio{};
  std::optional<double> max_match_pt_ratio{};
  std::optional<double> min_cst_pT{};
  std::optional<double> max_cst_pT{};
  std::optional<double> max_constituent_eta{};

  // Frame and scattered-e handling.
  bool lab_frame = false;
  bool drop_scat_electron = true;
  bool exclusive_match = false;
  bool antikt_in_breit = false;
};

// Pipeline-facing config (main copies Args -> CalibConfig).
struct CalibConfig {
  // Split + calibration behavior.
  double train_frac = 0.3;
  bool renorm = false;

  // Dataset / selection.
  std::size_t N_min = 50;
  double max_abs_eta = 3;
  bool lab_frame = false;
  bool drop_scat_electron = true;

  // Binning.
  std::vector<double> eta_bins{-2.4, -1.6, -0.8, 0.0, 0.8, 1.6, 2.4};
  std::vector<double> pT_bins{5, 7, 9, 12, 20, 30, 60}; // GeV

  // Jet-pair matching gates (used after clustering).
  double min_match_pT = 4.0;
  double max_match_eta = 3.0;
  double max_match_dR = 0.13;
  double min_match_pt_ratio = 0.3;
  double max_match_pt_ratio = 3.0;

  // Reco input and post-cluster hygiene.
  // Set to >=2 for stricter matching studies.
  std::size_t min_constituents_reco = 1;
  double min_constituent_pT = 0.0;
  double max_constituent_pT = std::numeric_limits<double>::infinity();
  double max_constituent_eta = std::numeric_limits<double>::infinity();
};

// Stable 64-bit hash (FNV-1a) over a serialized config string.
std::uint64_t fnv1a64(std::string_view s);

// Serialize CalibConfig to a deterministic string for hashing/metadata.
std::string serialize_config(const CalibConfig &cfg);

} // namespace jet_tools
