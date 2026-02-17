#include "jet_tools/include/root_io.h"

#include "jet_tools/include/event_progress.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace jet_tools {

namespace {

// Branch checks.
bool has_branch(TTree& tree, const char* name) {
  return tree.GetBranch(name) != nullptr;
}

void fail_missing_branch(TTree& tree, const char* name, const std::string& tag) {
  std::cerr << "[" << tag << "] tree '" << tree.GetName()
            << "' missing required branch '" << name << "'\n";
  std::exit(1);
}

void require_branch(TTree& tree, const char* name, const std::string& tag) {
  if (!has_branch(tree, name)) {
    fail_missing_branch(tree, name, tag);
  }
}

std::string make_read_message(const std::string& tag, const char* tree_name,
                              std::size_t index, std::size_t total) {
  return std::string("[") + tag + "] read " + tree_name + " " +
         std::to_string(index) + "/" + std::to_string(total);
}

}  // namespace

void ensure_parent(const std::string& path) {
  const std::filesystem::path p(path);
  if (p.has_parent_path()) {
    std::filesystem::create_directories(p.parent_path());
  }
}

TTree* get_required_tree(TFile& file, const char* tree_name, const char* tag) {
  auto* tree = dynamic_cast<TTree*>(file.Get(tree_name));
  if (tree != nullptr) {
    return tree;
  }
  std::cerr << "[" << tag << "] missing required tree '" << tree_name << "'\n";
  std::exit(1);
}

// Unified jet tree reader, different purposes use different flags.
void read_jet_tree(TTree& tree, std::size_t max_events, EventJets& events,
                   const std::string& tag, const JetTreeReadOptions& options) {
  require_branch(tree, "event_id", tag);

  // Accept either full 4-momentum or direct pT/eta/phi branches.
  const bool have_px = has_branch(tree, "px");
  const bool have_py = has_branch(tree, "py");
  const bool have_pz = has_branch(tree, "pz");
  const bool have_E = has_branch(tree, "E");
  const bool have_four_momentum = have_px && have_py && have_pz && have_E;

  const bool have_pT = has_branch(tree, "pT");
  const bool have_eta = has_branch(tree, "eta");
  const bool have_phi = has_branch(tree, "phi");
  const bool have_kinematics = have_pT && have_eta && have_phi;

  // Enforce 4-momentum when a caller needs exact momentum components.
  if (options.require_four_momentum && !have_four_momentum) {
    if (!have_px) {
      fail_missing_branch(tree, "px", tag);
    }
    if (!have_py) {
      fail_missing_branch(tree, "py", tag);
    }
    if (!have_pz) {
      fail_missing_branch(tree, "pz", tag);
    }
    fail_missing_branch(tree, "E", tag);
  }

  // Reject trees that expose neither supported schema.
  if (!have_four_momentum && !have_kinematics) {
    std::cerr << "[" << tag << "] tree '" << tree.GetName()
              << "' must provide either px/py/pz/E or pT/eta/phi\n";
    std::exit(1);
  }

  const bool have_n_constituents = has_branch(tree, "n_constituents");
  const bool have_n_constituent_count = have_n_constituents;

  if (options.require_n_constituents && !have_n_constituent_count) {
    fail_missing_branch(tree, "n_constituents", tag);
    std::exit(1);
  }

  // Only read constituent vectors when both paired branches are present.
  const bool have_constituents = has_branch(tree, "constituents");
  const bool have_constituent_E = has_branch(tree, "constituent_E");
  const bool have_constituent_vectors = have_constituents && have_constituent_E;

  if (options.require_constituent_vectors && !have_constituent_vectors) {
    if (!have_constituents) {
      fail_missing_branch(tree, "constituents", tag);
    }
    fail_missing_branch(tree, "constituent_E", tag);
  }

  ULong64_t event_id = 0;
  double px = std::numeric_limits<double>::quiet_NaN();
  double py = std::numeric_limits<double>::quiet_NaN();
  double pz = std::numeric_limits<double>::quiet_NaN();
  double E = std::numeric_limits<double>::quiet_NaN();
  double pT = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  int n_constituents = -1;
  std::vector<int>* constituents = nullptr;
  std::vector<double>* constituent_E = nullptr;

  tree.SetBranchAddress("event_id", &event_id);

  // Bind both schemas when present so we can prefer stored kinematics.
  if (have_four_momentum) {
    tree.SetBranchAddress("px", &px);
    tree.SetBranchAddress("py", &py);
    tree.SetBranchAddress("pz", &pz);
    tree.SetBranchAddress("E", &E);
  }
  if (have_kinematics) {
    tree.SetBranchAddress("pT", &pT);
    tree.SetBranchAddress("eta", &eta);
    tree.SetBranchAddress("phi", &phi);
  }

  if (have_n_constituents) {
    tree.SetBranchAddress("n_constituents", &n_constituents);
  }

  if (have_constituent_vectors) {
    tree.SetBranchAddress("constituents", &constituents);
    tree.SetBranchAddress("constituent_E", &constituent_E);
  }

  const std::size_t entries = static_cast<std::size_t>(tree.GetEntries());
  EventLimitGate event_gate(max_events);
  ProgressTicker progress;

  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(make_read_message(tag, tree.GetName(), i, entries));
    }
    tree.GetEntry(static_cast<Long64_t>(i));
    if (!event_gate.allow(event_id)) {
      continue;
    }

    SimpleJet jet;
    jet.px = std::numeric_limits<double>::quiet_NaN();
    jet.py = std::numeric_limits<double>::quiet_NaN();
    jet.pz = std::numeric_limits<double>::quiet_NaN();
    jet.E = std::numeric_limits<double>::quiet_NaN();
    jet.pT = std::numeric_limits<double>::quiet_NaN();
    jet.eta = std::numeric_limits<double>::quiet_NaN();
    jet.phi = std::numeric_limits<double>::quiet_NaN();
    jet.n_constituents = 0;

    // Prefer stored pT/eta/phi when available.
    if (have_four_momentum) {
      jet.px = px;
      jet.py = py;
      jet.pz = pz;
      jet.E = E;
    }
    if (have_kinematics) {
      jet.pT = pT;
      jet.eta = eta;
      jet.phi = phi;
    } else {
      jet.pT = std::hypot(px, py);
      jet.eta = (jet.pT > 0.0) ? std::asinh(pz / jet.pT) : std::numeric_limits<double>::quiet_NaN();
      jet.phi = std::atan2(py, px);
    }

    // Copy constituent count when the input branch exists.
    if (have_n_constituent_count) {
      jet.n_constituents = static_cast<std::size_t>(std::max(0, n_constituents));
    }

    if (have_constituent_vectors) {
      if (constituents != nullptr) {
        jet.constituents = *constituents;
      }
      if (constituent_E != nullptr) {
        jet.constituent_E = *constituent_E;
      }
    }

    events[event_id].push_back(std::move(jet));
  }

  progress.finish(make_read_message(tag, tree.GetName(), entries, entries));
}

// Shared reader for same-frame match trees.
std::vector<TruthRecoMatchRow> read_match_tree(TTree& tree, std::size_t max_events,
                                               const std::string& tag) {
  require_branch(tree, "event", tag);
  require_branch(tree, "truth_index", tag);
  require_branch(tree, "reco_index", tag);
  require_branch(tree, "dR", tag);

  TruthRecoMatchRow row;
  tree.SetBranchAddress("event", &row.event);
  tree.SetBranchAddress("truth_index", &row.truth_index);
  tree.SetBranchAddress("reco_index", &row.reco_index);
  tree.SetBranchAddress("dR", &row.dR);

  std::vector<TruthRecoMatchRow> rows;
  const std::size_t entries = static_cast<std::size_t>(tree.GetEntries());
  rows.reserve(entries);

  EventLimitGate event_gate(max_events);
  ProgressTicker progress;
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(make_read_message(tag, tree.GetName(), i, entries));
    }
    tree.GetEntry(static_cast<Long64_t>(i));
    if (!event_gate.allow(row.event)) {
      continue;
    }
    rows.push_back(row);
  }

  progress.finish(make_read_message(tag, tree.GetName(), entries, entries));
  return rows;
}

void setup_truth_reco_match_tree(TTree& tree, TruthRecoMatchRow& row) {
  tree.Branch("event", &row.event);
  tree.Branch("truth_index", &row.truth_index);
  tree.Branch("reco_index", &row.reco_index);
  tree.Branch("dR", &row.dR);
}

void setup_best_match_tree(TTree& tree, CrossFrameBestMatchRow& row) {
  tree.Branch("event", &row.event);
  tree.Branch("centauro_index", &row.centauro_index);
  tree.Branch("antikt_index", &row.antikt_index);
  tree.Branch("dR", &row.dR);
}

void setup_all_match_tree(TTree& tree, CrossFrameAllMatchRow& row) {
  tree.Branch("event", &row.event);
  tree.Branch("centauro_index", &row.centauro_index);
  tree.Branch("antikt_indices", &row.antikt_indices);
}

}  // namespace jet_tools
