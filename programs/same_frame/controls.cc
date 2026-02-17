/* Run:
./run.sh programs/same_frame/controls.cc \
  -m data/jets/matched/same_frame_matches.root -j data/jets/raw/jets_0-999_native_cuts.root \
  -o data/graphs/controls.root

./run.sh programs/same_frame/controls.cc \
  -m data/jets/matched/same_frame_matches.root -j data/jets/raw/jets_0-999_native_cuts.root \
  -o data/graphs/controls.root \
  --abs-eta-max 3 --pT-min 5 --match-dR 0.2
*/
#include <Rtypes.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>

#include "jet_tools/include/event_progress.h"
#include "jet_tools/include/root_io.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

//=====================================================================================================
// Types And Data Models.
//=====================================================================================================
struct Args {
  std::string matches_input;
  std::string jets_input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double max_eta = std::numeric_limits<double>::infinity();
  double min_pT = 0.0;
  double max_match_dR = 0.2;
  bool use_lab = false;
};

struct JetKey {
  ULong64_t event = 0;
  int index = -1;
};

struct JetKeyHash {
  std::size_t operator()(const JetKey& key) const noexcept {
    const std::size_t h1 = std::hash<ULong64_t>{}(key.event);
    const std::size_t h2 = std::hash<int>{}(key.index);
    return h1 ^ (h2 + 0x9e3779b9U + (h1 << 6U) + (h1 >> 2U));
  }
};

struct JetKeyEq {
  bool operator()(const JetKey& a, const JetKey& b) const noexcept {
    return a.event == b.event && a.index == b.index;
  }
};

using EventJets = jet_tools::EventJets;

struct SummaryRow {
  int alg = -1;
  // By the current definition of matching, truth_matched = reco_matched.
  // The current split is kept just in case a different definition of matching is needed.
  ULong64_t truth_total = 0;    // Total truth jets
  ULong64_t truth_matched = 0;  // Truth jets with a reco pair
  ULong64_t reco_total = 0;     // Total reco jets
  ULong64_t reco_matched = 0;   // Total reco jets with a truth pair
  double truth_efficiency = 0.0;  // truth_matched / truth_total
  double reco_efficiency = 0.0;   // reco_matched / reco_total
};

struct DomainEventCounts {
  std::unordered_map<ULong64_t, int> truth_per_event;
  std::unordered_map<ULong64_t, int> reco_per_event;
  std::unordered_map<ULong64_t, int> truth_matched_per_event;
  std::unordered_map<ULong64_t, int> reco_matched_per_event;
};

struct AlgData {
  std::string label;
  SummaryRow summary;
  DomainEventCounts counts;

  // Truth and reco jet kinematics.
  std::vector<double> truth_pT;
  std::vector<double> truth_eta;
  std::vector<double> truth_eta_lab;  // Occasionally used for debugging.
  std::vector<double> reco_pT;
  std::vector<double> reco_eta;

  // Kinematics only for matched truth jets, usually plotted as one of the axis.
  std::vector<double> truth_pT_matched;
  std::vector<double> truth_eta_matched;
  std::vector<double> match_dR;

  std::vector<double> truth_jets_per_event_pT;
  std::vector<double> truth_jets_per_event_eta;

  std::vector<double> truth_nconst;
  std::vector<double> reco_nconst;

  std::vector<std::size_t> truth_total_counts_pT;   // All jets (matched and unmatched) that pass the pT cut.
  std::vector<std::size_t> truth_matched_counts_pT; // Matched jets that pass the pT cut.
  std::vector<std::size_t> reco_total_counts_pT;
  std::vector<std::size_t> reco_matched_counts_pT;

  std::vector<std::size_t> truth_total_counts_eta;  // All jets (matched and unmatched) that pass the eta cut.
  std::vector<std::size_t> truth_matched_counts_eta;// Matched jets that pass the eta cut.
  std::vector<std::size_t> reco_total_counts_eta;
  std::vector<std::size_t> reco_matched_counts_eta;
};

//=====================================================================================================
// CLI Helpers.
//=====================================================================================================
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -m MATCHES.root -j JETS.root -o OUTPUT.root [--max-events N]\n"
            << "       [--abs-eta-max VAL] [--pT-min VAL] [--match-dR VAL] [--use-lab]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-m" || v == "--matches-input" || v == "--matches") && i + 1 < argc) {
      args.matches_input = argv[++i];
    } else if ((v == "-j" || v == "--jets-input") && i + 1 < argc) {
      args.jets_input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[controls] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--abs-eta-max" && i + 1 < argc) {
      args.max_eta = std::atof(argv[++i]);
    } else if (v == "--pT-min" && i + 1 < argc) {
      args.min_pT = std::atof(argv[++i]);
    } else if (v == "--match-dR" && i + 1 < argc) {
      args.max_match_dR = std::atof(argv[++i]);
    } else if (v == "--use-lab") {
      args.use_lab = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.matches_input.empty() || args.jets_input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (std::isfinite(args.max_eta) && args.max_eta <= 0.0) {
    std::cerr << "[controls] abs-eta-max must be > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.min_pT) || args.min_pT < 0.0) {
    std::cerr << "[controls] pT-min must be finite and >= 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.max_match_dR) || args.max_match_dR <= 0.0) {
    std::cerr << "[controls] match-dR must be finite and > 0\n";
    std::exit(1);
  }
  return args;
}

//=====================================================================================================
// Selection And Binning Helpers.
//=====================================================================================================
bool pass_jet_cuts(const jet_tools::SimpleJet& jet, const Args& args) {
  if (!std::isfinite(jet.pT) || jet.pT < args.min_pT) {
    return false;
  }
  if (!std::isfinite(jet.eta)) {
    return false;
  }
  if (std::isfinite(args.max_eta) && std::abs(jet.eta) > args.max_eta) {
    return false;
  }
  return true;
}

std::vector<double> make_uniform_edges(double min_v, double max_v, int n_bins) {
  std::vector<double> edges;
  if (!(max_v > min_v) || n_bins < 1) {
    return edges;
  }
  edges.reserve(static_cast<std::size_t>(n_bins + 1));
  const double step = (max_v - min_v) / static_cast<double>(n_bins);
  for (int i = 0; i <= n_bins; ++i) {
    edges.push_back(min_v + step * static_cast<double>(i));
  }
  return edges;
}

std::vector<double> build_pT_edges(double min_pT) {
  const std::vector<double> base = {5.0, 7.0, 8.0, 12.0, 20.0, 60.0};
  if (!(min_pT > base.front())) {
    return base;
  }

  std::vector<double> edges;
  edges.push_back(min_pT);
  for (double edge : base) {
    if (edge > min_pT) {
      edges.push_back(edge);
    }
  }
  if (edges.size() < 2) {
    edges.push_back(min_pT + 1.0);
  }
  return edges;
}

int find_bin(const std::vector<double>& edges, double x) {
  if (!std::isfinite(x) || edges.size() < 2) {
    return -1;
  }
  for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
    const bool is_last = (i + 2 == edges.size());
    if (x >= edges[i] && (x < edges[i + 1] || (is_last && x <= edges[i + 1]))) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

//=====================================================================================================
// Histogram Fill Helpers.
//=====================================================================================================
void fill_hist_values(TH1D& hist, const std::vector<double>& values) {
  for (double value : values) {
    if (std::isfinite(value)) {
      hist.Fill(value);
    }
  }
}

void fill_truth_eta_lab(const EventJets& truth_jets, const Args& args,
                        std::vector<double>& truth_eta_lab) {
  truth_eta_lab.clear();
  for (const auto& event_entry : truth_jets) {
    for (const auto& truth : event_entry.second) {
      if (!pass_jet_cuts(truth, args)) {
        continue;
      }
      truth_eta_lab.push_back(truth.eta);
    }
  }
}

void fill_event_hists(const DomainEventCounts& counts, TH1D& h_truth, TH1D& h_reco,
                      TH1D& h_reco_unmatched) {
  std::unordered_set<ULong64_t> event_ids;
  event_ids.reserve(counts.truth_per_event.size() + counts.reco_per_event.size());
  for (const auto& kv : counts.truth_per_event) {
    event_ids.insert(kv.first);
  }
  for (const auto& kv : counts.reco_per_event) {
    event_ids.insert(kv.first);
  }

  auto get_count = [](const std::unordered_map<ULong64_t, int>& values, ULong64_t event_id) {
    const auto it = values.find(event_id);
    return (it == values.end()) ? 0 : it->second;
  };

  for (ULong64_t event_id : event_ids) {
    const int n_truth = get_count(counts.truth_per_event, event_id);
    const int n_reco = get_count(counts.reco_per_event, event_id);
    const int n_reco_matched = get_count(counts.reco_matched_per_event, event_id);
    const int n_reco_unmatched = std::max(0, n_reco - n_reco_matched);
    h_truth.Fill(static_cast<double>(n_truth));
    h_reco.Fill(static_cast<double>(n_reco));
    h_reco_unmatched.Fill(static_cast<double>(n_reco_unmatched));
  }
}

//=====================================================================================================
// Core Data Collection.
//=====================================================================================================
void collect_alg_data(const std::string& label, int alg, const EventJets& truth_jets,
                      const EventJets& reco_jets,
                      const std::vector<jet_tools::TruthRecoMatchRow>& matches,
                      const Args& args, const std::vector<double>& pT_edges,
                      const std::vector<double>& eta_edges, AlgData& out) {
  out = AlgData{};
  out.label = label;
  out.summary.alg = alg;

  out.truth_total_counts_pT.assign(pT_edges.size() > 1 ? pT_edges.size() - 1 : 0, 0);
  out.truth_matched_counts_pT.assign(pT_edges.size() > 1 ? pT_edges.size() - 1 : 0, 0);
  out.reco_total_counts_pT.assign(pT_edges.size() > 1 ? pT_edges.size() - 1 : 0, 0);
  out.reco_matched_counts_pT.assign(pT_edges.size() > 1 ? pT_edges.size() - 1 : 0, 0);

  out.truth_total_counts_eta.assign(eta_edges.size() > 1 ? eta_edges.size() - 1 : 0, 0);
  out.truth_matched_counts_eta.assign(eta_edges.size() > 1 ? eta_edges.size() - 1 : 0, 0);
  out.reco_total_counts_eta.assign(eta_edges.size() > 1 ? eta_edges.size() - 1 : 0, 0);
  out.reco_matched_counts_eta.assign(eta_edges.size() > 1 ? eta_edges.size() - 1 : 0, 0);

  for (const auto& event_entry : truth_jets) {
    int truth_count = 0;
    for (const auto& truth : event_entry.second) {
      if (truth.n_constituents > 0) {
        out.truth_nconst.push_back(static_cast<double>(truth.n_constituents));
      }
      if (!pass_jet_cuts(truth, args)) {
        continue;
      }
      ++truth_count;
      ++out.summary.truth_total;
      out.truth_pT.push_back(truth.pT);
      out.truth_eta.push_back(truth.eta);
      out.truth_jets_per_event_pT.push_back(truth.pT);
      out.truth_jets_per_event_eta.push_back(truth.eta);

      const int pT_bin = find_bin(pT_edges, truth.pT);
      if (pT_bin >= 0) {
        out.truth_total_counts_pT[static_cast<std::size_t>(pT_bin)] += 1;
      }
      const int eta_bin = find_bin(eta_edges, truth.eta);
      if (eta_bin >= 0) {
        out.truth_total_counts_eta[static_cast<std::size_t>(eta_bin)] += 1;
      }
    }
    out.counts.truth_per_event[event_entry.first] = truth_count;
  }

  for (const auto& event_entry : reco_jets) {
    int reco_count = 0;
    for (const auto& reco : event_entry.second) {
      if (reco.n_constituents > 0) {
        out.reco_nconst.push_back(static_cast<double>(reco.n_constituents));
      }
      if (!pass_jet_cuts(reco, args)) {
        continue;
      }
      ++reco_count;
      ++out.summary.reco_total;
      out.reco_pT.push_back(reco.pT);
      out.reco_eta.push_back(reco.eta);

      const int pT_bin = find_bin(pT_edges, reco.pT);
      if (pT_bin >= 0) {
        out.reco_total_counts_pT[static_cast<std::size_t>(pT_bin)] += 1;
      }
      const int eta_bin = find_bin(eta_edges, reco.eta);
      if (eta_bin >= 0) {
        out.reco_total_counts_eta[static_cast<std::size_t>(eta_bin)] += 1;
      }
    }
    out.counts.reco_per_event[event_entry.first] = reco_count;
  }

  std::unordered_set<JetKey, JetKeyHash, JetKeyEq> matched_truth_keys;
  std::unordered_set<JetKey, JetKeyHash, JetKeyEq> matched_reco_keys;
  matched_truth_keys.reserve(matches.size());
  matched_reco_keys.reserve(matches.size());

  jet_tools::ProgressTicker progress;
  for (std::size_t i = 0; i < matches.size(); ++i) {
    if (progress.should_report(i)) {
      const std::string message = std::string("[controls] collect ") + label + " " +
                                  std::to_string(i) + "/" + std::to_string(matches.size());
      progress.report(message);
    }

    const auto& match = matches[i];
    if (match.truth_index < 0 || match.reco_index < 0) {
      continue;
    }
    if (!std::isfinite(match.dR) || match.dR > args.max_match_dR) {
      continue;
    }

    const auto truth_it = truth_jets.find(match.event);
    const auto reco_it = reco_jets.find(match.event);
    if (truth_it == truth_jets.end() || reco_it == reco_jets.end()) {
      continue;
    }

    const std::size_t truth_index = static_cast<std::size_t>(match.truth_index);
    const std::size_t reco_index = static_cast<std::size_t>(match.reco_index);
    if (truth_index >= truth_it->second.size() || reco_index >= reco_it->second.size()) {
      continue;
    }

    const auto& truth = truth_it->second[truth_index];
    const auto& reco = reco_it->second[reco_index];
    if (!pass_jet_cuts(truth, args) || !pass_jet_cuts(reco, args)) {
      continue;
    }

    const JetKey truth_key{match.event, match.truth_index};
    const JetKey reco_key{match.event, match.reco_index};
    const bool inserted_truth = matched_truth_keys.insert(truth_key).second;
    const bool inserted_reco = matched_reco_keys.insert(reco_key).second;

    if (inserted_truth) {
      out.truth_pT_matched.push_back(truth.pT);
      out.truth_eta_matched.push_back(truth.eta);
      out.match_dR.push_back(match.dR);
      ++out.counts.truth_matched_per_event[match.event];

      const int pT_bin = find_bin(pT_edges, truth.pT);
      if (pT_bin >= 0) {
        out.truth_matched_counts_pT[static_cast<std::size_t>(pT_bin)] += 1;
      }
      const int eta_bin = find_bin(eta_edges, truth.eta);
      if (eta_bin >= 0) {
        out.truth_matched_counts_eta[static_cast<std::size_t>(eta_bin)] += 1;
      }
    }

    if (inserted_reco) {
      ++out.counts.reco_matched_per_event[match.event];

      const int pT_bin = find_bin(pT_edges, reco.pT);
      if (pT_bin >= 0) {
        out.reco_matched_counts_pT[static_cast<std::size_t>(pT_bin)] += 1;
      }
      const int eta_bin = find_bin(eta_edges, reco.eta);
      if (eta_bin >= 0) {
        out.reco_matched_counts_eta[static_cast<std::size_t>(eta_bin)] += 1;
      }
    }

  }

  const std::string done_message = std::string("[controls] collect ") + label + " " +
                                   std::to_string(matches.size()) + "/" +
                                   std::to_string(matches.size());
  progress.finish(done_message);

  out.summary.truth_matched = static_cast<ULong64_t>(matched_truth_keys.size());
  out.summary.reco_matched = static_cast<ULong64_t>(matched_reco_keys.size());
  out.summary.truth_efficiency =
      (out.summary.truth_total > 0)
          ? static_cast<double>(out.summary.truth_matched) / static_cast<double>(out.summary.truth_total)
          : 0.0;
  out.summary.reco_efficiency =
      (out.summary.reco_total > 0)
          ? static_cast<double>(out.summary.reco_matched) / static_cast<double>(out.summary.reco_total)
          : 0.0;
}

//=====================================================================================================
// Summary And Output Helpers.
//=====================================================================================================
void write_canvas(TFile& fout, TDirectory* dir, TCanvas& canvas) {
  if (dir) {
    dir->cd();
  } else {
    fout.cd();
  }
  canvas.Write();
  fout.cd();
}

void set_style_pair(TH1D& a, TH1D& b) {
  a.SetStats(false);
  b.SetStats(false);
  a.SetLineColor(kBlue + 1);
  b.SetLineColor(kRed + 1);
  a.SetLineWidth(2);
  b.SetLineWidth(2);
}

std::size_t count_norm_events(const AlgData (&data)[2]) {
  std::unordered_set<ULong64_t> event_ids;
  for (const auto& d : data) {
    for (const auto& kv : d.counts.truth_per_event) {
      event_ids.insert(kv.first);
    }
    for (const auto& kv : d.counts.reco_per_event) {
      event_ids.insert(kv.first);
    }
  }
  return event_ids.size();
}

void fill_binomial_bin(TH1D& hist, int bin, std::size_t num, std::size_t den) {
  const double value = (den > 0) ? static_cast<double>(num) / static_cast<double>(den) : 0.0;
  const double err =
      (den > 0)
          ? std::sqrt(value * (1.0 - value) / static_cast<double>(den))
          : 0.0;
  hist.SetBinContent(bin, value);
  hist.SetBinError(bin, err);
}

void fill_control_hists(const AlgData (&data)[2], TH1D& h_eff_pT_antikt,
                        TH1D& h_eff_pT_cent, TH1D& h_eff_eta_antikt,
                        TH1D& h_eff_eta_cent,
                        TH1D& h_reco_efficiency_pT_antikt,
                        TH1D& h_reco_efficiency_pT_cent,
                        TH1D& h_reco_efficiency_eta_antikt,
                        TH1D& h_reco_efficiency_eta_cent,
                        TH1D& h_truth_unmatched_frac_pT_antikt,
                        TH1D& h_truth_unmatched_frac_pT_cent,
                        TH1D& h_truth_unmatched_frac_eta_antikt,
                        TH1D& h_truth_unmatched_frac_eta_cent,
                        TH1D& h_reco_unmatched_frac_pT_antikt,
                        TH1D& h_reco_unmatched_frac_pT_cent,
                        TH1D& h_reco_unmatched_frac_eta_antikt,
                        TH1D& h_reco_unmatched_frac_eta_cent) {
  for (int alg = 0; alg < 2; ++alg) {
    const AlgData& d = data[alg];
    for (std::size_t b = 0; b < d.truth_total_counts_pT.size(); ++b) {
      const std::size_t truth_total = d.truth_total_counts_pT[b];
      const std::size_t truth_matched = d.truth_matched_counts_pT[b];
      const std::size_t reco_total = d.reco_total_counts_pT[b];
      const std::size_t reco_matched = d.reco_matched_counts_pT[b];
      const std::size_t truth_unmatched =
          (truth_total >= truth_matched) ? (truth_total - truth_matched) : 0;
      const std::size_t reco_unmatched =
          (reco_total >= reco_matched) ? (reco_total - reco_matched) : 0;

      const int root_bin = static_cast<int>(b + 1);
      if (alg == 0) {
        fill_binomial_bin(h_eff_pT_antikt, root_bin, truth_matched, truth_total);
        fill_binomial_bin(h_reco_efficiency_pT_antikt, root_bin, reco_matched, reco_total);
        fill_binomial_bin(h_truth_unmatched_frac_pT_antikt, root_bin, truth_unmatched, truth_total);
        fill_binomial_bin(h_reco_unmatched_frac_pT_antikt, root_bin, reco_unmatched, reco_total);
      } else {
        fill_binomial_bin(h_eff_pT_cent, root_bin, truth_matched, truth_total);
        fill_binomial_bin(h_reco_efficiency_pT_cent, root_bin, reco_matched, reco_total);
        fill_binomial_bin(h_truth_unmatched_frac_pT_cent, root_bin, truth_unmatched,
                          truth_total);
        fill_binomial_bin(h_reco_unmatched_frac_pT_cent, root_bin, reco_unmatched, reco_total);
      }
    }

    for (std::size_t b = 0; b < d.truth_total_counts_eta.size(); ++b) {
      const std::size_t truth_total = d.truth_total_counts_eta[b];
      const std::size_t truth_matched = d.truth_matched_counts_eta[b];
      const std::size_t reco_total = d.reco_total_counts_eta[b];
      const std::size_t reco_matched = d.reco_matched_counts_eta[b];
      const std::size_t truth_unmatched =
          (truth_total >= truth_matched) ? (truth_total - truth_matched) : 0;
      const std::size_t reco_unmatched =
          (reco_total >= reco_matched) ? (reco_total - reco_matched) : 0;

      const int root_bin = static_cast<int>(b + 1);
      if (alg == 0) {
        fill_binomial_bin(h_eff_eta_antikt, root_bin, truth_matched, truth_total);
        fill_binomial_bin(h_reco_efficiency_eta_antikt, root_bin, reco_matched, reco_total);
        fill_binomial_bin(h_truth_unmatched_frac_eta_antikt, root_bin, truth_unmatched, truth_total);
        fill_binomial_bin(h_reco_unmatched_frac_eta_antikt, root_bin, reco_unmatched, reco_total);
      } else {
        fill_binomial_bin(h_eff_eta_cent, root_bin, truth_matched, truth_total);
        fill_binomial_bin(h_reco_efficiency_eta_cent, root_bin, reco_matched, reco_total);
        fill_binomial_bin(h_truth_unmatched_frac_eta_cent, root_bin, truth_unmatched, truth_total);
        fill_binomial_bin(h_reco_unmatched_frac_eta_cent, root_bin, reco_unmatched, reco_total);
      }
    }
  }
}

void fill_eff_dR_hists(const AlgData (&data)[2], TH1D& h_eff_dR_thresh_antikt,
                       TH1D& h_eff_dR_thresh_cent) {
  for (int alg = 0; alg < 2; ++alg) {
    const AlgData& d = data[alg];
    std::size_t truth_total = 0;
    for (std::size_t value : d.truth_total_counts_pT) {
      truth_total += value;
    }

    std::vector<double> dR = d.match_dR;
    std::sort(dR.begin(), dR.end());
    TH1D* hist = (alg == 0) ? &h_eff_dR_thresh_antikt : &h_eff_dR_thresh_cent;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
      const double threshold = hist->GetXaxis()->GetBinUpEdge(bin);
      const auto it = std::upper_bound(dR.begin(), dR.end(), threshold);
      const std::size_t matched = static_cast<std::size_t>(it - dR.begin());
      fill_binomial_bin(*hist, bin, matched, truth_total);
    }
  }
}

void write_summary_tree(const AlgData (&data)[2], TDirectory* hist_dir, TFile& fout) {
  if (hist_dir) {
    hist_dir->cd();
  } else {
    fout.cd();
  }

  TTree summary("summary", "paper summary");
  SummaryRow row;
  summary.Branch("alg", &row.alg);
  summary.Branch("truth_total", &row.truth_total);
  summary.Branch("truth_matched", &row.truth_matched);
  summary.Branch("reco_total", &row.reco_total);
  summary.Branch("reco_matched", &row.reco_matched);
  summary.Branch("truth_efficiency", &row.truth_efficiency);
  summary.Branch("reco_efficiency", &row.reco_efficiency);

  row = data[0].summary;
  summary.Fill();
  row = data[1].summary;
  summary.Fill();
  summary.Write();

  fout.cd();
}

void write_outputs(const Args& args, const AlgData (&data)[2], TFile& fout) {
  const std::vector<double> pT_edges = build_pT_edges(args.min_pT);
  const double max_eta_eff = std::isfinite(args.max_eta) ? args.max_eta : 3.0;
  const std::vector<double> eta_edges = make_uniform_edges(-max_eta_eff, max_eta_eff, 12);

  TDirectory* hist_dir = fout.mkdir("hist");

  TDirectory* antikt_dir = fout.mkdir("antikt");
  TDirectory* antikt_raw_dir = antikt_dir ? antikt_dir->mkdir("raw") : nullptr;
  TDirectory* antikt_raw_truth_dir = antikt_raw_dir ? antikt_raw_dir->mkdir("truth") : nullptr;
  TDirectory* antikt_raw_reco_dir = antikt_raw_dir ? antikt_raw_dir->mkdir("reco") : nullptr;
  TDirectory* antikt_raw_match_dir = antikt_raw_dir ? antikt_raw_dir->mkdir("matching") : nullptr;
  TDirectory* antikt_raw_unmatched_dir = antikt_raw_dir ? antikt_raw_dir->mkdir("unmatched") : nullptr;

  TDirectory* centauro_dir = fout.mkdir("centauro");
  TDirectory* centauro_raw_dir = centauro_dir ? centauro_dir->mkdir("raw") : nullptr;
  TDirectory* centauro_raw_truth_dir = centauro_raw_dir ? centauro_raw_dir->mkdir("truth") : nullptr;
  TDirectory* centauro_raw_reco_dir = centauro_raw_dir ? centauro_raw_dir->mkdir("reco") : nullptr;
  TDirectory* centauro_raw_match_dir = centauro_raw_dir ? centauro_raw_dir->mkdir("matching") : nullptr;
  TDirectory* centauro_raw_unmatched_dir = centauro_raw_dir ? centauro_raw_dir->mkdir("unmatched") : nullptr;

  if (hist_dir) {
    hist_dir->cd();
  }

  //===================================================================================================
  // Histogram Definitions.
  //===================================================================================================
  TH1D h_truth_pT_antikt("truth_jet_pT_antikt",
                         "Truth jet p_{T,Breit} (anti-k_{t});p_{T,Breit} [GeV];jets", 200, 0.0, 200.0);
  TH1D h_truth_pT_cent("truth_jet_pT_centauro",
                       "Truth jet p_{T,Breit} (Centauro);p_{T,Breit} [GeV];jets", 200, 0.0, 200.0);
  TH1D h_truth_eta_antikt("truth_jet_eta_breit_antikt",
                          "Truth jet #eta_{Breit} (anti-k_{t});#eta_{Breit};jets", 120, -6.0, 6.0);
  TH1D h_truth_eta_cent("truth_jet_eta_breit_centauro",
                        "Truth jet #eta_{Breit} (Centauro);#eta_{Breit};jets", 120, -6.0, 6.0);
  TH1D h_truth_eta_lab_antikt("truth_jet_eta_lab_antikt",
                              "Truth jet #eta_{lab} (anti-k_{t});#eta_{lab};jets", 120, -6.0, 6.0);
  TH1D h_truth_eta_lab_cent("truth_jet_eta_lab_centauro",
                            "Truth jet #eta_{lab} (Centauro);#eta_{lab};jets", 120, -6.0, 6.0);
  TH1D h_truth_pT_matched_antikt(
      "truth_jet_pT_matched_antikt",
      "Matched truth jet p_{T,Breit} per event (anti-k_{t});p_{T,Breit} [GeV];jets/event", 200, 0.0, 200.0);
  TH1D h_truth_pT_matched_cent(
      "truth_jet_pT_matched_centauro",
      "Matched truth jet p_{T,Breit} per event (Centauro);p_{T,Breit} [GeV];jets/event", 200, 0.0, 200.0);
  TH1D h_truth_eta_matched_antikt(
      "truth_jet_eta_matched_antikt",
      "Matched truth jet #eta_{Breit} per event (anti-k_{t});#eta_{Breit};jets/event", 120, -6.0, 6.0);
  TH1D h_truth_eta_matched_cent(
      "truth_jet_eta_matched_centauro",
      "Matched truth jet #eta_{Breit} per event (Centauro);#eta_{Breit};jets/event", 120, -6.0, 6.0);

  TH1D h_reco_pT_antikt("reco_jet_pT_antikt",
                        "Reco jet p_{T,Breit} (anti-k_{t});p_{T,Breit} [GeV];jets", 200, 0.0, 200.0);
  TH1D h_reco_pT_cent("reco_jet_pT_centauro",
                      "Reco jet p_{T,Breit} (Centauro);p_{T,Breit} [GeV];jets", 200, 0.0, 200.0);
  TH1D h_reco_eta_antikt("reco_jet_eta_antikt",
                         "Reco jet #eta_{Breit} (anti-k_{t});#eta_{Breit};jets", 120, -6.0, 6.0);
  TH1D h_reco_eta_cent("reco_jet_eta_centauro",
                       "Reco jet #eta_{Breit} (Centauro);#eta_{Breit};jets", 120, -6.0, 6.0);

  TH1D h_match_dR_antikt("dr_antikt", "#DeltaR_{Breit} matched (anti-k_{t});#DeltaR_{Breit};jets",
                         120, 0.0, args.max_match_dR);
  TH1D h_match_dR_cent("dr_centauro", "#DeltaR_{Breit} matched (Centauro);#DeltaR_{Breit};jets",
                       120, 0.0, args.max_match_dR);

  TH1D h_njets_event_truth_antikt(
      "njets_event_truth_antikt", "Truth anti-k_{t} jets per event (Breit);N_{jets};events", 30, 0.0, 30.0);
  TH1D h_njets_event_truth_cent(
      "njets_event_truth_centauro", "Truth Centauro jets per event (Breit);N_{jets};events", 30, 0.0, 30.0);
  TH1D h_njets_event_reco_antikt(
      "njets_event_reco_antikt",
      "Reco anti-k_{t} jets per event (Breit, baseline selection);N_{jets};events", 30, 0.0, 30.0);
  TH1D h_njets_event_reco_cent(
      "njets_event_reco_centauro",
      "Reco Centauro jets per event (Breit, baseline selection);N_{jets};events", 30, 0.0, 30.0);
  TH1D h_njets_event_reco_unmatched_antikt(
      "njets_event_reco_unmatched_antikt",
      "Unmatched reco anti-k_{t} jets per event (Breit);N_{jets}^{unmatched};events", 30, 0.0, 30.0);
  TH1D h_njets_event_reco_unmatched_cent(
      "njets_event_reco_unmatched_centauro",
      "Unmatched reco Centauro jets per event (Breit);N_{jets}^{unmatched};events", 30, 0.0, 30.0);

  TH1D h_jets_per_event_truth_pT_antikt(
      "jets_event_truth_pT_antikt",
      "Truth jets per event vs p_{T,Breit} (anti-k_{t});p_{T,Breit} [GeV];jets/event",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_jets_per_event_truth_pT_cent(
      "jets_event_truth_pT_centauro",
      "Truth jets per event vs p_{T,Breit} (Centauro);p_{T,Breit} [GeV];jets/event",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_jets_per_event_truth_eta_antikt(
      "jets_event_truth_eta_antikt",
      "Truth jets per event vs #eta_{Breit} (anti-k_{t});#eta_{Breit};jets/event",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_jets_per_event_truth_eta_cent(
      "jets_event_truth_eta_centauro",
      "Truth jets per event vs #eta_{Breit} (Centauro);#eta_{Breit};jets/event",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());

  TH1D h_nconst_jet_truth_antikt(
      "nconst_jet_truth_antikt", "Truth anti-k_{t} particles per jet (Breit);N_{const};jets", 60, 0.0, 60.0);
  TH1D h_nconst_jet_truth_cent(
      "nconst_jet_truth_centauro", "Truth Centauro particles per jet (Breit);N_{const};jets", 60, 0.0, 60.0);
  TH1D h_nconst_jet_reco_antikt(
      "nconst_jet_reco_antikt", "Reco anti-k_{t} particles per jet (Breit);N_{const};jets", 60, 0.0, 60.0);
  TH1D h_nconst_jet_reco_cent(
      "nconst_jet_reco_centauro", "Reco Centauro particles per jet (Breit);N_{const};jets", 60, 0.0, 60.0);

  TH1D h_eff_pT_antikt(
      "eff_match_pT_antikt",
      "Matching efficiency vs truth p_{T,Breit} (anti-k_{t});p_{T,Breit}^{truth} [GeV];efficiency",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_eff_pT_cent(
      "eff_match_pT_centauro",
      "Matching efficiency vs truth p_{T,Breit} (Centauro);p_{T,Breit}^{truth} [GeV];efficiency",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_eff_eta_antikt(
      "eff_match_eta_antikt",
      "Matching efficiency vs truth #eta_{Breit} (anti-k_{t});#eta_{Breit}^{truth};efficiency",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_eff_eta_cent(
      "eff_match_eta_centauro",
      "Matching efficiency vs truth #eta_{Breit} (Centauro);#eta_{Breit}^{truth};efficiency",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());

  TH1D h_reco_efficiency_pT_antikt(
      "reco_efficiency_match_pT_antikt",
      "Reco-side matching efficiency vs reco p_{T,Breit} (anti-k_{t});p_{T,Breit}^{reco} [GeV];reco efficiency",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_reco_efficiency_pT_cent(
      "reco_efficiency_match_pT_centauro",
      "Reco-side matching efficiency vs reco p_{T,Breit} (Centauro);p_{T,Breit}^{reco} [GeV];reco efficiency",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_reco_efficiency_eta_antikt(
      "reco_efficiency_match_eta_antikt",
      "Reco-side matching efficiency vs reco #eta_{Breit} (anti-k_{t});#eta_{Breit}^{reco};reco efficiency",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_reco_efficiency_eta_cent(
      "reco_efficiency_match_eta_centauro",
      "Reco-side matching efficiency vs reco #eta_{Breit} (Centauro);#eta_{Breit}^{reco};reco efficiency",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());

  TH1D h_eff_dR_thresh_antikt(
      "eff_match_dR_thresh_antikt",
      "Matching efficiency vs #DeltaR_{Breit} threshold (anti-k_{t});max #DeltaR_{Breit};efficiency",
      40, 0.0, args.max_match_dR);
  TH1D h_eff_dR_thresh_cent(
      "eff_match_dR_thresh_centauro",
      "Matching efficiency vs #DeltaR_{Breit} threshold (Centauro);max #DeltaR_{Breit};efficiency",
      40, 0.0, args.max_match_dR);

  TH1D h_truth_unmatched_frac_pT_antikt(
      "truth_unmatched_frac_pT_antikt",
      "Truth unmatched fraction vs truth p_{T,Breit} (anti-k_{t});p_{T,Breit}^{truth} [GeV];fraction",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_truth_unmatched_frac_pT_cent(
      "truth_unmatched_frac_pT_centauro",
      "Truth unmatched fraction vs truth p_{T,Breit} (Centauro);p_{T,Breit}^{truth} [GeV];fraction",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_truth_unmatched_frac_eta_antikt(
      "truth_unmatched_frac_eta_antikt",
      "Truth unmatched fraction vs truth #eta_{Breit} (anti-k_{t});#eta_{Breit}^{truth};fraction",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_truth_unmatched_frac_eta_cent(
      "truth_unmatched_frac_eta_centauro",
      "Truth unmatched fraction vs truth #eta_{Breit} (Centauro);#eta_{Breit}^{truth};fraction",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_reco_unmatched_frac_pT_antikt(
      "reco_unmatched_frac_pT_antikt",
      "Reco unmatched fraction vs reco p_{T,Breit} (anti-k_{t});p_{T,Breit}^{reco} [GeV];fraction",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_reco_unmatched_frac_pT_cent(
      "reco_unmatched_frac_pT_centauro",
      "Reco unmatched fraction vs reco p_{T,Breit} (Centauro);p_{T,Breit}^{reco} [GeV];fraction",
      static_cast<int>(pT_edges.size()) - 1, pT_edges.data());
  TH1D h_reco_unmatched_frac_eta_antikt(
      "reco_unmatched_frac_eta_antikt",
      "Reco unmatched fraction vs reco #eta_{Breit} (anti-k_{t});#eta_{Breit}^{reco};fraction",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());
  TH1D h_reco_unmatched_frac_eta_cent(
      "reco_unmatched_frac_eta_centauro",
      "Reco unmatched fraction vs reco #eta_{Breit} (Centauro);#eta_{Breit}^{reco};fraction",
      static_cast<int>(eta_edges.size()) - 1, eta_edges.data());

  //===================================================================================================
  // Histogram Styling.
  //===================================================================================================
  set_style_pair(h_truth_pT_antikt, h_truth_pT_cent);
  set_style_pair(h_truth_eta_antikt, h_truth_eta_cent);
  set_style_pair(h_truth_eta_lab_antikt, h_truth_eta_lab_cent);
  set_style_pair(h_truth_pT_matched_antikt, h_truth_pT_matched_cent);
  set_style_pair(h_truth_eta_matched_antikt, h_truth_eta_matched_cent);
  set_style_pair(h_reco_pT_antikt, h_reco_pT_cent);
  set_style_pair(h_reco_eta_antikt, h_reco_eta_cent);
  set_style_pair(h_match_dR_antikt, h_match_dR_cent);
  set_style_pair(h_njets_event_truth_antikt, h_njets_event_truth_cent);
  set_style_pair(h_njets_event_reco_antikt, h_njets_event_reco_cent);
  set_style_pair(h_njets_event_reco_unmatched_antikt, h_njets_event_reco_unmatched_cent);
  set_style_pair(h_jets_per_event_truth_pT_antikt, h_jets_per_event_truth_pT_cent);
  set_style_pair(h_jets_per_event_truth_eta_antikt, h_jets_per_event_truth_eta_cent);
  set_style_pair(h_nconst_jet_truth_antikt, h_nconst_jet_truth_cent);
  set_style_pair(h_nconst_jet_reco_antikt, h_nconst_jet_reco_cent);
  set_style_pair(h_eff_pT_antikt, h_eff_pT_cent);
  set_style_pair(h_eff_eta_antikt, h_eff_eta_cent);
  set_style_pair(h_reco_efficiency_pT_antikt, h_reco_efficiency_pT_cent);
  set_style_pair(h_reco_efficiency_eta_antikt, h_reco_efficiency_eta_cent);
  set_style_pair(h_eff_dR_thresh_antikt, h_eff_dR_thresh_cent);
  set_style_pair(h_truth_unmatched_frac_pT_antikt, h_truth_unmatched_frac_pT_cent);
  set_style_pair(h_truth_unmatched_frac_eta_antikt, h_truth_unmatched_frac_eta_cent);
  set_style_pair(h_reco_unmatched_frac_pT_antikt, h_reco_unmatched_frac_pT_cent);
  set_style_pair(h_reco_unmatched_frac_eta_antikt, h_reco_unmatched_frac_eta_cent);

  //===================================================================================================
  // Histogram Filling.
  //===================================================================================================
  fill_hist_values(h_truth_pT_antikt, data[0].truth_pT);
  fill_hist_values(h_truth_pT_cent, data[1].truth_pT);
  fill_hist_values(h_truth_eta_antikt, data[0].truth_eta);
  fill_hist_values(h_truth_eta_cent, data[1].truth_eta);
  fill_hist_values(h_truth_eta_lab_antikt, data[0].truth_eta_lab);
  fill_hist_values(h_truth_eta_lab_cent, data[1].truth_eta_lab);
  fill_hist_values(h_truth_pT_matched_antikt, data[0].truth_pT_matched);
  fill_hist_values(h_truth_pT_matched_cent, data[1].truth_pT_matched);
  fill_hist_values(h_truth_eta_matched_antikt, data[0].truth_eta_matched);
  fill_hist_values(h_truth_eta_matched_cent, data[1].truth_eta_matched);
  fill_hist_values(h_reco_pT_antikt, data[0].reco_pT);
  fill_hist_values(h_reco_pT_cent, data[1].reco_pT);
  fill_hist_values(h_reco_eta_antikt, data[0].reco_eta);
  fill_hist_values(h_reco_eta_cent, data[1].reco_eta);
  fill_hist_values(h_match_dR_antikt, data[0].match_dR);
  fill_hist_values(h_match_dR_cent, data[1].match_dR);
  fill_hist_values(h_jets_per_event_truth_pT_antikt, data[0].truth_jets_per_event_pT);
  fill_hist_values(h_jets_per_event_truth_pT_cent, data[1].truth_jets_per_event_pT);
  fill_hist_values(h_jets_per_event_truth_eta_antikt, data[0].truth_jets_per_event_eta);
  fill_hist_values(h_jets_per_event_truth_eta_cent, data[1].truth_jets_per_event_eta);
  fill_hist_values(h_nconst_jet_truth_antikt, data[0].truth_nconst);
  fill_hist_values(h_nconst_jet_truth_cent, data[1].truth_nconst);
  fill_hist_values(h_nconst_jet_reco_antikt, data[0].reco_nconst);
  fill_hist_values(h_nconst_jet_reco_cent, data[1].reco_nconst);

  fill_event_hists(data[0].counts, h_njets_event_truth_antikt, h_njets_event_reco_antikt,
                   h_njets_event_reco_unmatched_antikt);
  fill_event_hists(data[1].counts, h_njets_event_truth_cent, h_njets_event_reco_cent,
                   h_njets_event_reco_unmatched_cent);

  fill_control_hists(
      data, h_eff_pT_antikt, h_eff_pT_cent, h_eff_eta_antikt, h_eff_eta_cent,
      h_reco_efficiency_pT_antikt, h_reco_efficiency_pT_cent,
      h_reco_efficiency_eta_antikt, h_reco_efficiency_eta_cent,
      h_truth_unmatched_frac_pT_antikt, h_truth_unmatched_frac_pT_cent,
      h_truth_unmatched_frac_eta_antikt, h_truth_unmatched_frac_eta_cent,
      h_reco_unmatched_frac_pT_antikt, h_reco_unmatched_frac_pT_cent,
      h_reco_unmatched_frac_eta_antikt, h_reco_unmatched_frac_eta_cent);
  fill_eff_dR_hists(data, h_eff_dR_thresh_antikt, h_eff_dR_thresh_cent);

  //===================================================================================================
  // Normalization And Summary.
  //===================================================================================================
  const std::size_t norm_events = count_norm_events(data);
  if (norm_events > 0) {
    const double scale = 1.0 / static_cast<double>(norm_events);
    h_truth_pT_matched_antikt.Scale(scale);
    h_truth_pT_matched_cent.Scale(scale);
    h_truth_eta_matched_antikt.Scale(scale);
    h_truth_eta_matched_cent.Scale(scale);
    h_jets_per_event_truth_pT_antikt.Scale(scale);
    h_jets_per_event_truth_pT_cent.Scale(scale);
    h_jets_per_event_truth_eta_antikt.Scale(scale);
    h_jets_per_event_truth_eta_cent.Scale(scale);
  }

  write_summary_tree(data, hist_dir, fout);

  //===================================================================================================
  // Histogram Writes.
  //===================================================================================================
  if (hist_dir) {
    hist_dir->cd();
  } else {
    fout.cd();
  }
  h_truth_pT_antikt.Write();
  h_truth_pT_cent.Write();
  h_truth_eta_antikt.Write();
  h_truth_eta_cent.Write();
  h_truth_eta_lab_antikt.Write();
  h_truth_eta_lab_cent.Write();
  h_truth_pT_matched_antikt.Write();
  h_truth_pT_matched_cent.Write();
  h_truth_eta_matched_antikt.Write();
  h_truth_eta_matched_cent.Write();
  h_reco_pT_antikt.Write();
  h_reco_pT_cent.Write();
  h_reco_eta_antikt.Write();
  h_reco_eta_cent.Write();
  h_match_dR_antikt.Write();
  h_match_dR_cent.Write();
  h_njets_event_truth_antikt.Write();
  h_njets_event_truth_cent.Write();
  h_njets_event_reco_antikt.Write();
  h_njets_event_reco_cent.Write();
  h_njets_event_reco_unmatched_antikt.Write();
  h_njets_event_reco_unmatched_cent.Write();
  h_jets_per_event_truth_pT_antikt.Write();
  h_jets_per_event_truth_pT_cent.Write();
  h_jets_per_event_truth_eta_antikt.Write();
  h_jets_per_event_truth_eta_cent.Write();
  h_nconst_jet_truth_antikt.Write();
  h_nconst_jet_truth_cent.Write();
  h_nconst_jet_reco_antikt.Write();
  h_nconst_jet_reco_cent.Write();
  h_eff_pT_antikt.Write();
  h_eff_pT_cent.Write();
  h_eff_eta_antikt.Write();
  h_eff_eta_cent.Write();
  h_reco_efficiency_pT_antikt.Write();
  h_reco_efficiency_pT_cent.Write();
  h_reco_efficiency_eta_antikt.Write();
  h_reco_efficiency_eta_cent.Write();
  h_eff_dR_thresh_antikt.Write();
  h_eff_dR_thresh_cent.Write();
  h_truth_unmatched_frac_pT_antikt.Write();
  h_truth_unmatched_frac_pT_cent.Write();
  h_truth_unmatched_frac_eta_antikt.Write();
  h_truth_unmatched_frac_eta_cent.Write();
  h_reco_unmatched_frac_pT_antikt.Write();
  h_reco_unmatched_frac_pT_cent.Write();
  h_reco_unmatched_frac_eta_antikt.Write();
  h_reco_unmatched_frac_eta_cent.Write();
  fout.cd();

  //===================================================================================================
  // Canvas Writes.
  //===================================================================================================
  gStyle->SetOptStat(0);

  auto write_single_canvas = [&](TDirectory* dir, const char* name, const char* title,
                                 TH1D& hist, const char* draw_opt) {
    TCanvas canvas(name, title, 900, 600);
    hist.Draw(draw_opt);
    write_canvas(fout, dir, canvas);
  };

  write_single_canvas(antikt_raw_truth_dir, "c_truth_jet_pT", "Truth jet p_{T,Breit}",
                      h_truth_pT_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_truth_jet_eta", "Truth jet #eta_{Breit}",
                      h_truth_eta_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_truth_jet_eta_lab", "Truth jet #eta_{lab}",
                      h_truth_eta_lab_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_truth_jet_pT_matched",
                      "Matched truth jet p_{T,Breit} per event", h_truth_pT_matched_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_truth_jet_eta_matched",
                      "Matched truth jet #eta_{Breit} per event", h_truth_eta_matched_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_njets_event_truth",
                      "Truth anti-k_{t} jets per event (Breit)", h_njets_event_truth_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_jets_event_truth_pT",
                      "Truth jets per event vs p_{T,Breit}", h_jets_per_event_truth_pT_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_jets_event_truth_eta",
                      "Truth jets per event vs #eta_{Breit}",
                      h_jets_per_event_truth_eta_antikt, "HIST");
  write_single_canvas(antikt_raw_truth_dir, "c_nconst_jet_truth",
                      "Truth anti-k_{t} particles per jet (Breit)", h_nconst_jet_truth_antikt, "HIST");

  write_single_canvas(antikt_raw_reco_dir, "c_reco_jet_pT", "Reco jet p_{T,Breit}",
                      h_reco_pT_antikt, "HIST");
  write_single_canvas(antikt_raw_reco_dir, "c_reco_jet_eta", "Reco jet #eta_{Breit}",
                      h_reco_eta_antikt, "HIST");
  write_single_canvas(antikt_raw_reco_dir, "c_njets_event_reco",
                      "Reco anti-k_{t} jets per event (Breit)", h_njets_event_reco_antikt, "HIST");
  write_single_canvas(antikt_raw_reco_dir, "c_nconst_jet_reco",
                      "Reco anti-k_{t} particles per jet (Breit)", h_nconst_jet_reco_antikt, "HIST");

  write_single_canvas(antikt_raw_match_dir, "c_match_dR", "Matched #DeltaR_{Breit}",
                      h_match_dR_antikt, "HIST");
  write_single_canvas(antikt_raw_match_dir, "c_eff_match_pT",
                      "Matching efficiency vs truth p_{T,Breit}", h_eff_pT_antikt, "E1");
  write_single_canvas(antikt_raw_match_dir, "c_eff_match_eta",
                      "Matching efficiency vs truth #eta_{Breit}", h_eff_eta_antikt, "E1");
  write_single_canvas(antikt_raw_match_dir, "c_reco_efficiency_match_pT",
                      "Reco-side matching efficiency vs reco p_{T,Breit}",
                      h_reco_efficiency_pT_antikt, "E1");
  write_single_canvas(antikt_raw_match_dir, "c_reco_efficiency_match_eta",
                      "Reco-side matching efficiency vs reco #eta_{Breit}",
                      h_reco_efficiency_eta_antikt, "E1");
  write_single_canvas(antikt_raw_match_dir, "c_eff_match_dR_thresh",
                      "Matching efficiency vs #DeltaR_{Breit} threshold",
                      h_eff_dR_thresh_antikt, "E1");

  write_single_canvas(antikt_raw_unmatched_dir, "c_njets_event_reco_unmatched",
                      "Unmatched reco jets per event (Breit)",
                      h_njets_event_reco_unmatched_antikt, "HIST");
  write_single_canvas(antikt_raw_unmatched_dir, "c_truth_unmatched_frac_pT",
                      "Truth unmatched fraction vs p_{T,Breit}",
                      h_truth_unmatched_frac_pT_antikt, "E1");
  write_single_canvas(antikt_raw_unmatched_dir, "c_truth_unmatched_frac_eta",
                      "Truth unmatched fraction vs #eta_{Breit}",
                      h_truth_unmatched_frac_eta_antikt, "E1");
  write_single_canvas(antikt_raw_unmatched_dir, "c_reco_unmatched_frac_pT",
                      "Reco unmatched fraction vs p_{T,Breit}",
                      h_reco_unmatched_frac_pT_antikt, "E1");
  write_single_canvas(antikt_raw_unmatched_dir, "c_reco_unmatched_frac_eta",
                      "Reco unmatched fraction vs #eta_{Breit}",
                      h_reco_unmatched_frac_eta_antikt, "E1");

  write_single_canvas(centauro_raw_truth_dir, "c_truth_jet_pT", "Truth jet p_{T,Breit}",
                      h_truth_pT_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_truth_jet_eta", "Truth jet #eta_{Breit}",
                      h_truth_eta_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_truth_jet_eta_lab", "Truth jet #eta_{lab}",
                      h_truth_eta_lab_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_truth_jet_pT_matched",
                      "Matched truth jet p_{T,Breit} per event", h_truth_pT_matched_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_truth_jet_eta_matched",
                      "Matched truth jet #eta_{Breit} per event", h_truth_eta_matched_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_njets_event_truth",
                      "Truth Centauro jets per event (Breit)", h_njets_event_truth_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_jets_event_truth_pT",
                      "Truth jets per event vs p_{T,Breit}", h_jets_per_event_truth_pT_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_jets_event_truth_eta",
                      "Truth jets per event vs #eta_{Breit}",
                      h_jets_per_event_truth_eta_cent, "HIST");
  write_single_canvas(centauro_raw_truth_dir, "c_nconst_jet_truth",
                      "Truth Centauro particles per jet (Breit)", h_nconst_jet_truth_cent, "HIST");

  write_single_canvas(centauro_raw_reco_dir, "c_reco_jet_pT", "Reco jet p_{T,Breit}",
                      h_reco_pT_cent, "HIST");
  write_single_canvas(centauro_raw_reco_dir, "c_reco_jet_eta", "Reco jet #eta_{Breit}",
                      h_reco_eta_cent, "HIST");
  write_single_canvas(centauro_raw_reco_dir, "c_njets_event_reco",
                      "Reco Centauro jets per event (Breit)", h_njets_event_reco_cent, "HIST");
  write_single_canvas(centauro_raw_reco_dir, "c_nconst_jet_reco",
                      "Reco Centauro particles per jet (Breit)", h_nconst_jet_reco_cent, "HIST");

  write_single_canvas(centauro_raw_match_dir, "c_match_dR", "Matched #DeltaR_{Breit}",
                      h_match_dR_cent, "HIST");
  write_single_canvas(centauro_raw_match_dir, "c_eff_match_pT",
                      "Matching efficiency vs truth p_{T,Breit}", h_eff_pT_cent, "E1");
  write_single_canvas(centauro_raw_match_dir, "c_eff_match_eta",
                      "Matching efficiency vs truth #eta_{Breit}", h_eff_eta_cent, "E1");
  write_single_canvas(centauro_raw_match_dir, "c_reco_efficiency_match_pT",
                      "Reco-side matching efficiency vs reco p_{T,Breit}",
                      h_reco_efficiency_pT_cent, "E1");
  write_single_canvas(centauro_raw_match_dir, "c_reco_efficiency_match_eta",
                      "Reco-side matching efficiency vs reco #eta_{Breit}",
                      h_reco_efficiency_eta_cent, "E1");
  write_single_canvas(centauro_raw_match_dir, "c_eff_match_dR_thresh",
                      "Matching efficiency vs #DeltaR_{Breit} threshold",
                      h_eff_dR_thresh_cent, "E1");

  write_single_canvas(centauro_raw_unmatched_dir, "c_njets_event_reco_unmatched",
                      "Unmatched reco jets per event (Breit)",
                      h_njets_event_reco_unmatched_cent, "HIST");
  write_single_canvas(centauro_raw_unmatched_dir, "c_truth_unmatched_frac_pT",
                      "Truth unmatched fraction vs p_{T,Breit}",
                      h_truth_unmatched_frac_pT_cent, "E1");
  write_single_canvas(centauro_raw_unmatched_dir, "c_truth_unmatched_frac_eta",
                      "Truth unmatched fraction vs #eta_{Breit}",
                      h_truth_unmatched_frac_eta_cent, "E1");
  write_single_canvas(centauro_raw_unmatched_dir, "c_reco_unmatched_frac_pT",
                      "Reco unmatched fraction vs p_{T,Breit}",
                      h_reco_unmatched_frac_pT_cent, "E1");
  write_single_canvas(centauro_raw_unmatched_dir, "c_reco_unmatched_frac_eta",
                      "Reco unmatched fraction vs #eta_{Breit}",
                      h_reco_unmatched_frac_eta_cent, "E1");
}

}  // namespace

//=====================================================================================================
// Program Entry Point.
//=====================================================================================================
int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  TFile f_matches(args.matches_input.c_str(), "READ");
  if (!f_matches.IsOpen()) {
    std::cerr << "[controls] failed to open matches input " << args.matches_input << "\n";
    return 1;
  }

  TFile f_jets(args.jets_input.c_str(), "READ");
  if (!f_jets.IsOpen()) {
    std::cerr << "[controls] failed to open jets input " << args.jets_input << "\n";
    return 1;
  }

  const char* antikt_truth_name = args.use_lab ? "LabFrameAntiktTruth" : "BreitFrameAntiktTruth";
  const char* antikt_reco_name = args.use_lab ? "LabFrameAntiktReco" : "BreitFrameAntiktReco";
  const char* centauro_truth_name = args.use_lab ? "LabFrameCentauroTruth" : "BreitFrameCentauroTruth";
  const char* centauro_reco_name = args.use_lab ? "LabFrameCentauroReco" : "BreitFrameCentauroReco";

  auto* t_antikt_truth = jet_tools::get_required_tree(f_jets, antikt_truth_name, "controls");
  auto* t_antikt_reco = jet_tools::get_required_tree(f_jets, antikt_reco_name, "controls");
  auto* t_centauro_truth = jet_tools::get_required_tree(f_jets, centauro_truth_name, "controls");
  auto* t_centauro_reco = jet_tools::get_required_tree(f_jets, centauro_reco_name, "controls");

  auto* t_antikt_matches = jet_tools::get_required_tree(f_matches, "antikt_matches", "controls");
  auto* t_centauro_matches = jet_tools::get_required_tree(f_matches, "centauro_matches", "controls");

  jet_tools::JetTreeReadOptions jet_read_options;
  jet_read_options.require_n_constituents = true;

  EventJets antikt_truth;
  EventJets antikt_reco;
  EventJets centauro_truth;
  EventJets centauro_reco;
  jet_tools::read_jet_tree(*t_antikt_truth, args.max_events, antikt_truth, "controls", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_reco, args.max_events, antikt_reco, "controls", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_truth, args.max_events, centauro_truth, "controls", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_reco, args.max_events, centauro_reco, "controls", jet_read_options);

  const std::vector<jet_tools::TruthRecoMatchRow> antikt_matches =
      jet_tools::read_match_tree(*t_antikt_matches, args.max_events, "controls");
  const std::vector<jet_tools::TruthRecoMatchRow> centauro_matches =
      jet_tools::read_match_tree(*t_centauro_matches, args.max_events, "controls");

  const std::vector<double> pT_edges = build_pT_edges(args.min_pT);
  const double max_eta_eff = std::isfinite(args.max_eta) ? args.max_eta : 3.0;
  const std::vector<double> eta_edges = make_uniform_edges(-max_eta_eff, max_eta_eff, 12);

  AlgData data[2];
  collect_alg_data("antikt", 0, antikt_truth, antikt_reco, antikt_matches, args, pT_edges, eta_edges, data[0]);
  collect_alg_data("centauro", 1, centauro_truth, centauro_reco, centauro_matches, args, pT_edges, eta_edges, data[1]);

  if (args.use_lab) {
    data[0].truth_eta_lab = data[0].truth_eta;
    data[1].truth_eta_lab = data[1].truth_eta;
  } else {
    auto* t_antikt_truth_lab = dynamic_cast<TTree*>(f_jets.Get("LabFrameAntiktTruth"));
    auto* t_centauro_truth_lab = dynamic_cast<TTree*>(f_jets.Get("LabFrameCentauroTruth"));
    if (t_antikt_truth_lab != nullptr && t_centauro_truth_lab != nullptr) {
      EventJets antikt_truth_lab;
      EventJets centauro_truth_lab;
      jet_tools::read_jet_tree(*t_antikt_truth_lab, args.max_events, antikt_truth_lab,
                               "controls", jet_read_options);
      jet_tools::read_jet_tree(*t_centauro_truth_lab, args.max_events, centauro_truth_lab,
                               "controls", jet_read_options);
      fill_truth_eta_lab(antikt_truth_lab, args, data[0].truth_eta_lab);
      fill_truth_eta_lab(centauro_truth_lab, args, data[1].truth_eta_lab);
    } else {
      data[0].truth_eta_lab = data[0].truth_eta;
      data[1].truth_eta_lab = data[1].truth_eta;
    }
  }

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[controls] failed to open output " << args.output << "\n";
    return 1;
  }

  write_outputs(args, data, fout);

  fout.Close();
  f_matches.Close();
  f_jets.Close();
  std::cout << "[controls] wrote " << args.output << "\n";
  return 0;
}
