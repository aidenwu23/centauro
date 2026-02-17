/* Run:
./run.sh programs/same_frame/match_jets.cc \
  -i data/jets/raw/jets_0-999_native_cuts.root -o data/jets/matched/same_frame_matches.root

*/
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>

#include "jet_tools/include/math_helpers.h"
#include "jet_tools/include/root_io.h"
#include "jet_tools/include/types.h"

namespace {

struct Args {
  std::string input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double match_dR = 0.2;
  double min_match_pt_ratio = 0.3;
  double max_match_pt_ratio = 3.0;
  bool use_lab = false;
};

struct EventJetPair {
  std::vector<jet_tools::SimpleJet> truth;
  std::vector<jet_tools::SimpleJet> reco;
};

using EventMap = std::unordered_map<ULong64_t, EventJetPair>;

void usage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " -i CLUSTERED.root -o OUTPUT.root [--max-events N]\n"
            << "       [--match-dR VAL] [--min_match_pt_ratio VAL] [--max_match_pt_ratio VAL]\n"
            << "       [--use-lab]\n";
}

Args parse_args(int argc, char *argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[match_jets_sf] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--match-dR" && i + 1 < argc) {
      args.match_dR = std::atof(argv[++i]);
    } else if (v == "--min_match_pt_ratio" && i + 1 < argc) {
      args.min_match_pt_ratio = std::atof(argv[++i]);
    } else if (v == "--max_match_pt_ratio" && i + 1 < argc) {
      args.max_match_pt_ratio = std::atof(argv[++i]);
    } else if (v == "--use-lab") {
      args.use_lab = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!std::isfinite(args.match_dR) || args.match_dR <= 0.0) {
    std::cerr << "[match_jets_sf] match-dR must be finite and > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.min_match_pt_ratio) || args.min_match_pt_ratio <= 0.0) {
    std::cerr << "[match_jets_sf] min_match_pt_ratio must be finite and > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.max_match_pt_ratio) ||
      args.max_match_pt_ratio < args.min_match_pt_ratio) {
    std::cerr << "[match_jets_sf] max_match_pt_ratio must be >= min_match_pt_ratio\n";
    std::exit(1);
  }
  return args;
}

std::vector<std::pair<std::size_t, std::size_t>> match_truth_reco_pairs(
    const std::vector<jet_tools::SimpleJet> &truth_jets,
    const std::vector<jet_tools::SimpleJet> &reco_jets, double max_match_dR,
    double min_match_pt_ratio, double max_match_pt_ratio) {
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  if (truth_jets.empty() || reco_jets.empty()) {
    return pairs;
  }

  std::vector<std::size_t> reco_used(reco_jets.size(), 0);
  std::vector<std::size_t> truth_order(truth_jets.size());
  std::iota(truth_order.begin(), truth_order.end(), 0);
  std::sort(truth_order.begin(), truth_order.end(), [&](std::size_t a, std::size_t b) {
    return truth_jets[a].pT > truth_jets[b].pT;
  });

  const double max_dR2 = max_match_dR * max_match_dR;
  pairs.reserve(std::min(truth_jets.size(), reco_jets.size()));

  // Match truth jets in descending pT, each reco jet used at most once.
  for (const std::size_t truth_index : truth_order) {
    const auto &truth_jet = truth_jets[truth_index];
    double best_dR2 = max_dR2;
    std::size_t best_reco_index = reco_jets.size();

    for (std::size_t reco_index = 0; reco_index < reco_jets.size(); ++reco_index) {
      if (reco_used[reco_index]) {
        continue;
      }

      const auto &reco_jet = reco_jets[reco_index];
      const double dphi = jet_tools::delta_phi(truth_jet.phi, reco_jet.phi);
      const double deta = truth_jet.eta - reco_jet.eta;
      const double dR2 = dphi * dphi + deta * deta;
      if (!(std::isfinite(dR2) && dR2 < best_dR2)) {
        continue;
      }

      const double pT_ratio = (truth_jet.pT > 0.0)
                                  ? (reco_jet.pT / truth_jet.pT)
                                  : std::numeric_limits<double>::quiet_NaN();
      if (!std::isfinite(pT_ratio) || pT_ratio < min_match_pt_ratio ||
          pT_ratio > max_match_pt_ratio) {
        continue;
      }

      best_dR2 = dR2;
      best_reco_index = reco_index;
    }

    if (best_reco_index != reco_jets.size()) {
      reco_used[best_reco_index] = 1;
      pairs.push_back({truth_index, best_reco_index});
    }
  }

  return pairs;
}

void fill_matches_for_alg(const EventMap &events, double max_match_dR,
                          double min_match_pt_ratio, double max_match_pt_ratio,
                          TTree &tree, jet_tools::TruthRecoMatchRow &row,
                          const std::string &progress_tag) {
  std::vector<ULong64_t> event_ids;
  event_ids.reserve(events.size());
  for (const auto &kv : events) {
    event_ids.push_back(kv.first);
  }
  std::sort(event_ids.begin(), event_ids.end());

  for (std::size_t i = 0; i < event_ids.size(); ++i) {
    if (i % 100 == 0) {
      std::cout << "[match_jets_sf] " << progress_tag << " event " << i << "/"
                << event_ids.size() << "\r" << std::flush;
    }

    const ULong64_t event_id = event_ids[i];
    const auto it = events.find(event_id);
    if (it == events.end()) {
      continue;
    }

    const auto &truth_jets = it->second.truth;
    const auto &reco_jets = it->second.reco;
    if (truth_jets.empty() || reco_jets.empty()) {
      continue;
    }

    const auto pairs = match_truth_reco_pairs(truth_jets, reco_jets, max_match_dR,
                                              min_match_pt_ratio, max_match_pt_ratio);

    for (const auto &pair : pairs) {
      const std::size_t truth_index = pair.first;
      const std::size_t reco_index = pair.second;
      if (truth_index >= truth_jets.size() || reco_index >= reco_jets.size()) {
        continue;
      }
      const auto &truth_jet = truth_jets[truth_index];
      const auto &reco_jet = reco_jets[reco_index];
      const double dphi = jet_tools::delta_phi(truth_jet.phi, reco_jet.phi);
      const double deta = truth_jet.eta - reco_jet.eta;
      const double dR2 = dphi * dphi + deta * deta;

      row.event = event_id;
      row.truth_index = static_cast<int>(truth_index);
      row.reco_index = static_cast<int>(reco_index);
      row.dR = std::sqrt(dR2);
      tree.Fill();
    }
  }

  std::cout << "[match_jets_sf] " << progress_tag << " event " << event_ids.size()
            << "/" << event_ids.size() << "\n";
}

} // namespace

int main(int argc, char *argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  TFile fin(args.input.c_str(), "READ");
  if (!fin.IsOpen()) {
    std::cerr << "[match_jets_sf] failed to open input " << args.input << "\n";
    return 1;
  }

  const char *antikt_truth_name = args.use_lab ? "LabFrameAntiktTruth" : "BreitFrameAntiktTruth";
  const char *antikt_reco_name = args.use_lab ? "LabFrameAntiktReco" : "BreitFrameAntiktReco";
  const char *centauro_truth_name = args.use_lab ? "LabFrameCentauroTruth" : "BreitFrameCentauroTruth";
  const char *centauro_reco_name = args.use_lab ? "LabFrameCentauroReco" : "BreitFrameCentauroReco";

  auto *t_antikt_truth =
      jet_tools::get_required_tree(fin, antikt_truth_name, "match_jets_sf");
  auto *t_antikt_reco =
      jet_tools::get_required_tree(fin, antikt_reco_name, "match_jets_sf");
  auto *t_centauro_truth =
      jet_tools::get_required_tree(fin, centauro_truth_name, "match_jets_sf");
  auto *t_centauro_reco =
      jet_tools::get_required_tree(fin, centauro_reco_name, "match_jets_sf");

  jet_tools::EventJets antikt_truth_events;
  jet_tools::EventJets antikt_reco_events;
  jet_tools::EventJets centauro_truth_events;
  jet_tools::EventJets centauro_reco_events;
  jet_tools::JetTreeReadOptions jet_read_options;
  jet_read_options.require_four_momentum = true;
  jet_tools::read_jet_tree(*t_antikt_truth, args.max_events, antikt_truth_events,
                           "match_jets_sf", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_reco, args.max_events, antikt_reco_events,
                           "match_jets_sf", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_truth, args.max_events, centauro_truth_events,
                           "match_jets_sf", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_reco, args.max_events, centauro_reco_events,
                           "match_jets_sf", jet_read_options);

  EventMap antikt_events;
  EventMap centauro_events;
  for (const auto& event_entry : antikt_truth_events) {
    antikt_events[event_entry.first].truth = event_entry.second;
  }
  for (const auto& event_entry : antikt_reco_events) {
    antikt_events[event_entry.first].reco = event_entry.second;
  }
  for (const auto& event_entry : centauro_truth_events) {
    centauro_events[event_entry.first].truth = event_entry.second;
  }
  for (const auto& event_entry : centauro_reco_events) {
    centauro_events[event_entry.first].reco = event_entry.second;
  }

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[match_jets_sf] failed to open output " << args.output << "\n";
    return 1;
  }

  TTree antikt_matches("antikt_matches", "truth-reco anti-kt matches");
  TTree centauro_matches("centauro_matches", "truth-reco centauro matches");

  jet_tools::TruthRecoMatchRow antikt_row;
  jet_tools::TruthRecoMatchRow centauro_row;
  jet_tools::setup_truth_reco_match_tree(antikt_matches, antikt_row);
  jet_tools::setup_truth_reco_match_tree(centauro_matches, centauro_row);

  fill_matches_for_alg(antikt_events, args.match_dR, args.min_match_pt_ratio,
                       args.max_match_pt_ratio, antikt_matches, antikt_row,
                       "anti-kt");
  fill_matches_for_alg(centauro_events, args.match_dR, args.min_match_pt_ratio,
                       args.max_match_pt_ratio, centauro_matches, centauro_row,
                       "centauro");

  antikt_matches.Write();
  centauro_matches.Write();

  fout.Close();
  fin.Close();
  return 0;
}
