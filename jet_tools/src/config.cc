#include "jet_tools/include/config.h"

#include <sstream>

namespace jet_tools {

// Simple hash for config provenance tracking (ensures we can check if settings changed).
std::uint64_t fnv1a64(std::string_view s) {
  std::uint64_t h = 14695981039346656037ULL;
  for (unsigned char c : s) {
    h ^= static_cast<std::uint64_t>(c);
    h *= 1099511628211ULL;
  }
  return h;
}

std::string serialize_config(const CalibConfig &cfg) {
  std::ostringstream oss;
  oss << "train_frac=" << cfg.train_frac;
  oss << ";lab_frame=" << (cfg.lab_frame ? 1 : 0);
  oss << ";renorm=" << (cfg.renorm ? 1 : 0);
  oss << ";N_min=" << cfg.N_min;
  oss << ";max_abs_eta=" << cfg.max_abs_eta;
  oss << ";min_constituents_reco=" << cfg.min_constituents_reco;
  oss << ";min_match_pT=" << cfg.min_match_pT;
  oss << ";max_match_eta=" << cfg.max_match_eta;
  oss << ";max_match_dR=" << cfg.max_match_dR;
  oss << ";min_match_pt_ratio=" << cfg.min_match_pt_ratio;
  oss << ";max_match_pt_ratio=" << cfg.max_match_pt_ratio;
  oss << ";min_constituent_pT=" << cfg.min_constituent_pT;
  oss << ";max_constituent_pT=" << cfg.max_constituent_pT;
  oss << ";max_constituent_eta=" << cfg.max_constituent_eta;
  oss << ";drop_scat_electron=" << (cfg.drop_scat_electron ? 1 : 0);
  oss << ";eta_bins=";
  for (std::size_t i = 0; i < cfg.eta_bins.size(); ++i) {
    if (i)
      oss << ",";
    oss << cfg.eta_bins[i];
  }
  oss << ";pT_bins=";
  for (std::size_t i = 0; i < cfg.pT_bins.size(); ++i) {
    if (i)
      oss << ",";
    oss << cfg.pT_bins[i];
  }
  return oss.str();
}

} // namespace jet_tools
