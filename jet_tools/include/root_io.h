#pragma once

#include <RtypesCore.h>

#include <cstddef>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "jet_tools/include/types.h"

class TFile;
class TTree;

namespace jet_tools {

// Per-event jet storage keyed by event_id.
using EventJets = std::unordered_map<ULong64_t, std::vector<SimpleJet>>;

// Truth-to-reco pair used by same-frame match outputs.
struct TruthRecoMatchRow {
  ULong64_t event = 0;
  int truth_index = -1;
  int reco_index = -1;
  double dR = std::numeric_limits<double>::quiet_NaN();
};

// Nearest anti-kt jet for one Centauro jet in cross-frame matching.
struct CrossFrameBestMatchRow {
  ULong64_t event = 0;
  int centauro_index = -1;
  int antikt_index = -1;
  double dR = std::numeric_limits<double>::quiet_NaN();
};

// All anti-kt candidates inside the dR cut for one Centauro jet.
struct CrossFrameAllMatchRow {
  ULong64_t event = 0;
  int centauro_index = -1;
  std::vector<int> antikt_indices;
};

void ensure_parent(const std::string& path);

TTree* get_required_tree(TFile& file, const char* tree_name, const char* tag);

struct JetTreeReadOptions {
  // Require px/py/pz/E and derive pT, eta, phi from that 4-momentum.
  bool require_four_momentum = false;
  // Require n_constituents.
  bool require_n_constituents = false;
  // Require both constituent_indices and constituent_E vectors.
  bool require_constituent_vectors = false;
};

// Read one jet tree schema into EventJets with strict requirements from options.
void read_jet_tree(TTree& tree, std::size_t max_events, EventJets& events,
                   const std::string& tag,
                   const JetTreeReadOptions& options = JetTreeReadOptions{});

// Read same-frame match rows (event, truth_index, reco_index, dR).
std::vector<TruthRecoMatchRow> read_match_tree(TTree& tree, std::size_t max_events, const std::string& tag);

// Define same-frame match tree branches.
void setup_truth_reco_match_tree(TTree& tree, TruthRecoMatchRow& row);
// Define cross-frame best-match branches.
void setup_best_match_tree(TTree& tree, CrossFrameBestMatchRow& row);
// Define cross-frame all-match branches.
void setup_all_match_tree(TTree& tree, CrossFrameAllMatchRow& row);

}  // namespace jet_tools
