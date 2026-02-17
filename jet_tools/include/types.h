#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <Math/Vector4D.h>

namespace jet_tools {

using P4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

struct Provenance {
  std::string input_path;
  std::vector<std::string> cli_args;
  std::uint64_t rng_seed = 0;
  std::uint64_t cfg_hash = 0;
};

struct MatchRecord {
  std::size_t event_index = 0;
  std::size_t truth_index = 0;
  std::size_t reco_index = 0;
  bool is_train = false;
  int alg = 0;  // 0 = anti-kt, 1 = Centauro
  double Q2 = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
  double pT_truth = 0.0;
  double eta_truth = 0.0;
  double phi_truth = 0.0;
  double px_truth = 0.0;
  double py_truth = 0.0;
  double pT_reco_raw = 0.0;
  double pT_calib = 0.0;
  double resp_calib = 0.0;
  double eta_reco = 0.0;
  double phi_reco = 0.0;
  double dR = 0.0;
  double z_inv_truth = std::numeric_limits<double>::quiet_NaN();
  double z_inv_reco = std::numeric_limits<double>::quiet_NaN();
};

struct TruthRecord {
  std::size_t event_index = 0;
  bool is_train = false;
  int alg = 0;
  double Q2 = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
  double pT_truth = 0.0;
  double eta_truth = 0.0;
  double phi_truth = 0.0;
  double px_truth = 0.0;
  double py_truth = 0.0;
  double z_inv_truth = std::numeric_limits<double>::quiet_NaN();
  bool is_matched = false;
};

struct RecoRecord {
  std::size_t event_index = 0;
  bool is_train = false;
  int alg = 0;
  double Q2 = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
  double pT_reco_raw = 0.0;
  double pT_calib = 0.0;
  double eta_reco = 0.0;
  double phi_reco = 0.0;
  double z_inv_reco = std::numeric_limits<double>::quiet_NaN();
  bool is_matched = false;
};

struct SimpleJet {
  double pT = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double E = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  std::vector<int> constituents;
  std::vector<double> constituent_px;
  std::vector<double> constituent_py;
  std::vector<double> constituent_pz;
  std::vector<double> constituent_E;
  std::size_t n_constituents = 0;
  std::size_t n_reco_constituents = 0;
  double z_inv = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace jet_tools
