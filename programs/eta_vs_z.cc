/*
Run:
./run.sh programs/eta_vs_z.cc \
  -i data/jets/raw/eta_vs_z_jets.root \
  -o data/graphs/eta_vs_z.root

./run.sh programs/eta_vs_z.cc \
  -i data/jets/raw/eta_vs_z_jets.root \
  -o data/graphs/eta_vs_z.root --max-events 50000 --min-jet-pT 2
*/

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <Rtypes.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

namespace {

struct Args {
  std::string input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double min_jet_pT = 0.0;
  std::size_t min_constituents = 1;
};

struct FillStats {
  std::size_t scanned = 0;
  std::size_t kept = 0;
  std::size_t q2_filled = 0;
  std::size_t events_seen = 0;
};

void usage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " -i CLUSTERED.root -o OUTPUT.root [--max-events N]\n"
            << "       [--min-jet-pT VAL] [--min-constituents N]\n";
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
        std::cerr << "[eta_vs_z] --max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--min-jet-pT" && i + 1 < argc) {
      args.min_jet_pT = std::atof(argv[++i]);
    } else if (v == "--min-constituents" && i + 1 < argc) {
      const long long raw = std::stoll(argv[++i]);
      if (raw <= 0) {
        std::cerr << "[eta_vs_z] --min-constituents must be >= 1\n";
        std::exit(1);
      }
      args.min_constituents = static_cast<std::size_t>(raw);
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!std::isfinite(args.min_jet_pT) || args.min_jet_pT < 0.0) {
    std::cerr << "[eta_vs_z] --min-jet-pT must be finite and >= 0\n";
    std::exit(1);
  }
  return args;
}

void ensure_parent(const std::string &output_path) {
  const std::filesystem::path out(output_path);
  const std::filesystem::path parent = out.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
}

TH2D make_eta_z_hist(const std::string &name, const std::string &title,
                     int z_bins, int eta_bins) {
  return TH2D(name.c_str(), title.c_str(), z_bins, 0.0, 1.0, eta_bins, -6.0,
              6.0);
}

std::string q2_tag(double value) {
  return std::to_string(static_cast<int>(value));
}

struct HistGroup {
  TH2D all_truth;
  TH2D all_reco;
  std::vector<TH2D> q2_truth;
  std::vector<TH2D> q2_reco;
  TH2D q2_over_truth;
  TH2D q2_over_reco;

  HistGroup(const std::string &prefix, const std::string &label, int z_bins,
            int eta_bins, const std::vector<double> &q2_edges,
            double q2_over_min)
      : all_truth(make_eta_z_hist(
            prefix + "_eta_z_truth",
            label +
                " truth jets: #eta_{Breit} vs z_{jet} (inclusive);z_{jet};#eta_{Breit}",
            z_bins, eta_bins)),
        all_reco(make_eta_z_hist(
            prefix + "_eta_z_reco",
            label +
                " reco jets: #eta_{Breit} vs z_{jet} (inclusive);z_{jet};#eta_{Breit}",
            z_bins, eta_bins)),
        q2_over_truth(make_eta_z_hist(
            prefix + "_eta_z_truth_Q2_" + q2_tag(q2_over_min) + "_inf",
            label + " truth jets: #eta_{Breit} vs z_{jet} (Q^{2} #geq " +
                q2_tag(q2_over_min) + " GeV^{2});z_{jet};#eta_{Breit}",
            z_bins, eta_bins)),
        q2_over_reco(make_eta_z_hist(
            prefix + "_eta_z_reco_Q2_" + q2_tag(q2_over_min) + "_inf",
            label + " reco jets: #eta_{Breit} vs z_{jet} (Q^{2} #geq " +
                q2_tag(q2_over_min) + " GeV^{2});z_{jet};#eta_{Breit}",
            z_bins, eta_bins)) {
    q2_truth.reserve(q2_edges.size() - 1U);
    q2_reco.reserve(q2_edges.size() - 1U);
    for (std::size_t i = 0; i + 1 < q2_edges.size(); ++i) {
      const double lo = q2_edges[i];
      const double hi = q2_edges[i + 1U];
      const std::string range =
          std::to_string(static_cast<int>(lo)) + "_" +
          std::to_string(static_cast<int>(hi));
      const std::string title_range =
          std::to_string(static_cast<int>(lo)) + "-" +
          std::to_string(static_cast<int>(hi));

      q2_truth.emplace_back(make_eta_z_hist(
          prefix + "_eta_z_truth_Q2_" + range,
          label + " truth jets: #eta_{Breit} vs z_{jet} (Q^{2} " + title_range +
              " GeV^{2});z_{jet};#eta_{Breit}",
          z_bins, eta_bins));

      q2_reco.emplace_back(make_eta_z_hist(
          prefix + "_eta_z_reco_Q2_" + range,
          label + " reco jets: #eta_{Breit} vs z_{jet} (Q^{2} " + title_range +
              " GeV^{2});z_{jet};#eta_{Breit}",
          z_bins, eta_bins));
    }
  }
};

struct SampleHists {
  TH2D *all = nullptr;
  std::vector<TH2D> *q2 = nullptr;
  TH2D *q2_over = nullptr;
};

SampleHists select_hists(HistGroup &group, bool is_reco) {
  if (is_reco) {
    return {&group.all_reco, &group.q2_reco, &group.q2_over_reco};
  }
  return {&group.all_truth, &group.q2_truth, &group.q2_over_truth};
}

bool find_q2_bin(double q2, const std::vector<double> &edges, std::size_t &out) {
  for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
    if (q2 >= edges[i] && q2 < edges[i + 1U]) {
      out = i;
      return true;
    }
  }
  return false;
}

void write_hist_group(HistGroup &group, TDirectory *dir) {
  if (dir) {
    dir->cd();
  }
  group.all_truth.Write();
  group.all_reco.Write();
  for (auto &hist : group.q2_truth) {
    hist.Write();
  }
  for (auto &hist : group.q2_reco) {
    hist.Write();
  }
  group.q2_over_truth.Write();
  group.q2_over_reco.Write();
}

FillStats fill_from_tree(TTree &tree, const Args &args, const std::string &tag,
                         const std::vector<double> &q2_edges,
                         double q2_over_min, SampleHists hists) {
  FillStats stats;

  ULong64_t event_id = 0;
  double pT = 0.0;
  double eta = std::numeric_limits<double>::quiet_NaN();
  double z_inv = std::numeric_limits<double>::quiet_NaN();
  double q2 = std::numeric_limits<double>::quiet_NaN();
  int n_constituents = 0;

  if (!tree.GetBranch("event_id") || !tree.GetBranch("pT") ||
      !tree.GetBranch("eta") || !tree.GetBranch("n_constituents") ||
      !tree.GetBranch("z_inv") || !tree.GetBranch("Q2")) {
    std::cerr << "[eta_vs_z] tree '" << tree.GetName()
              << "' missing required branches\n";
    std::exit(1);
  }

  tree.SetBranchAddress("event_id", &event_id);
  tree.SetBranchAddress("pT", &pT);
  tree.SetBranchAddress("eta", &eta);
  tree.SetBranchAddress("n_constituents", &n_constituents);
  tree.SetBranchAddress("z_inv", &z_inv);
  tree.SetBranchAddress("Q2", &q2);

  const std::size_t entries = static_cast<std::size_t>(tree.GetEntries());
  ULong64_t last_event_id = std::numeric_limits<ULong64_t>::max();

  for (std::size_t i = 0; i < entries; ++i) {
    tree.GetEntry(static_cast<Long64_t>(i));
    ++stats.scanned;

    // Count events by event_id transitions while streaming jet rows once.
    if (event_id != last_event_id) {
      last_event_id = event_id;
      ++stats.events_seen;
      if (stats.events_seen % 100 == 0) {
        std::cout << "[eta_vs_z] " << tag << " event " << stats.events_seen
                  << "\r" << std::flush;
      }
      if (args.max_events != std::numeric_limits<std::size_t>::max() &&
          stats.events_seen > args.max_events) {
        break;
      }
    }

    if (!std::isfinite(pT) || pT < args.min_jet_pT) {
      continue;
    }
    if (!std::isfinite(eta)) {
      continue;
    }
    if (n_constituents < static_cast<int>(args.min_constituents)) {
      continue;
    }
    if (!std::isfinite(z_inv) || z_inv < 0.0 || z_inv > 1.0) {
      continue;
    }

    hists.all->Fill(z_inv, eta);
    ++stats.kept;

    if (!std::isfinite(q2) || q2 <= 0.0) {
      continue;
    }
    std::size_t q2_index = 0;
    if (find_q2_bin(q2, q2_edges, q2_index)) {
      hists.q2->at(q2_index).Fill(z_inv, eta);
      ++stats.q2_filled;
    }
    if (q2 >= q2_over_min) {
      hists.q2_over->Fill(z_inv, eta);
    }
  }

  std::cout << "[eta_vs_z] " << tag << " event " << stats.events_seen << "\n";
  return stats;
}

} // namespace

int main(int argc, char *argv[]) {
  const Args args = parse_args(argc, argv);
  ensure_parent(args.output);

  TFile fin(args.input.c_str(), "READ");
  if (!fin.IsOpen()) {
    std::cerr << "[eta_vs_z] failed to open input " << args.input << "\n";
    return 1;
  }

  auto *t_cent_truth = dynamic_cast<TTree *>(fin.Get("BreitFrameCentauroTruth"));
  auto *t_cent_reco = dynamic_cast<TTree *>(fin.Get("BreitFrameCentauroReco"));
  auto *t_antikt_truth = dynamic_cast<TTree *>(fin.Get("BreitFrameAntiktTruth"));
  auto *t_antikt_reco = dynamic_cast<TTree *>(fin.Get("BreitFrameAntiktReco"));

  if (!t_cent_truth || !t_cent_reco || !t_antikt_truth || !t_antikt_reco) {
    std::cerr << "[eta_vs_z] missing required Breit-frame jet trees\n";
    return 1;
  }

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[eta_vs_z] failed to open output " << args.output << "\n";
    return 1;
  }

  gStyle->SetOptStat(0);

  constexpr int kZBins = 120;
  constexpr int kEtaBins = 240;
  constexpr double kQ2OverflowMin = 800.0;
  const std::vector<double> q2_edges{
      100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0, 6400.0, 12800.0};

  HistGroup h_centauro("centauro", "Centauro", kZBins, kEtaBins, q2_edges,
                       kQ2OverflowMin);
  HistGroup h_antikt("antikt", "anti-k_{t}", kZBins, kEtaBins, q2_edges,
                     kQ2OverflowMin);

  const FillStats s_cent_truth =
      fill_from_tree(*t_cent_truth, args, "centauro truth", q2_edges,
                     kQ2OverflowMin, select_hists(h_centauro, false));
  const FillStats s_cent_reco =
      fill_from_tree(*t_cent_reco, args, "centauro reco", q2_edges,
                     kQ2OverflowMin, select_hists(h_centauro, true));
  const FillStats s_antikt_truth =
      fill_from_tree(*t_antikt_truth, args, "antikt truth", q2_edges,
                     kQ2OverflowMin, select_hists(h_antikt, false));
  const FillStats s_antikt_reco =
      fill_from_tree(*t_antikt_reco, args, "antikt reco", q2_edges,
                     kQ2OverflowMin, select_hists(h_antikt, true));

  TDirectory *cent_dir = fout.mkdir("centauro");
  TDirectory *antikt_dir = fout.mkdir("antikt");
  write_hist_group(h_centauro, cent_dir);
  write_hist_group(h_antikt, antikt_dir ? antikt_dir : &fout);

  fout.Close();
  fin.Close();

  const std::size_t total_kept = s_cent_truth.kept + s_cent_reco.kept +
                                 s_antikt_truth.kept + s_antikt_reco.kept;
  const std::size_t total_q2 = s_cent_truth.q2_filled + s_cent_reco.q2_filled +
                               s_antikt_truth.q2_filled +
                               s_antikt_reco.q2_filled;

  std::cout << "[eta_vs_z] kept jets: " << total_kept << "\n";
  std::cout << "[eta_vs_z] Q2-sliced fills: " << total_q2 << "\n";
  std::cout << "[eta_vs_z] wrote " << args.output << "\n";
  return 0;
}
