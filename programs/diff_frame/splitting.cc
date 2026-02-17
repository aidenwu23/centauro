/* Run:
./run.sh programs/diff_frame/splitting.cc \
  -m data/jets/matched/diff_frame_matches.root -j data/jets/raw/jets_0-999_lab_cuts.root \
  -o data/graphs/diff_frame_splitting.root

./run.sh programs/diff_frame/splitting.cc \
  -m data/jets/matched/diff_frame_matches.root -j data/jets/raw/jets_0-999_lab_cuts.root \
  -o data/graphs/diff_frame_splitting.root \
  --max-events 50000 --max-dR 0.3

--max-dR = maximum dR allowed between a pair of jets.
*/
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
#include <vector>

#include <Rtypes.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include "jet_tools/include/event_progress.h"
#include "jet_tools/include/root_io.h"

namespace {

struct Args {
  std::string matches_input;
  std::string jets_input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double max_dR = std::numeric_limits<double>::infinity();
};

struct BestMatchRow {
  ULong64_t event = 0;
  int centauro_index = -1;
  int antikt_index = -1;
  double dR = std::numeric_limits<double>::quiet_NaN();
};

struct AllMatchRow {
  ULong64_t event = 0;
  int centauro_index = -1;
  std::vector<int>* antikt_indices = nullptr;
};

using EventJets = jet_tools::EventJets;

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -m MATCHES.root -j JETS.root -o OUTPUT.root [--max-events N]\n"
            << "       [--max-dR VAL]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-m" || v == "--matches-input") && i + 1 < argc) {
      args.matches_input = argv[++i];
    } else if ((v == "-j" || v == "--jets-input") && i + 1 < argc) {
      args.jets_input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[splitting_cf] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--max-dR" && i + 1 < argc) {
      args.max_dR = std::atof(argv[++i]);
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.matches_input.empty() || args.jets_input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!std::isfinite(args.max_dR) && args.max_dR != std::numeric_limits<double>::infinity()) {
    std::cerr << "[splitting_cf] max-dR must be finite\n";
    std::exit(1);
  }
  if (std::isfinite(args.max_dR) && args.max_dR <= 0.0) {
    std::cerr << "[splitting_cf] max-dR must be > 0\n";
    std::exit(1);
  }
  return args;
}

struct DomainHists {
  TH1D* h_ratio = nullptr;
  TH2D* h_ratio_eta_lab = nullptr;
  TH2D* h_ratio_eta_breit = nullptr;
  TH1D* h_nconst_antikt = nullptr;
  TH1D* h_nconst_centauro = nullptr;
  TH1D* h_soc_count = nullptr;
  TH1D* h_soc_energy = nullptr;
  TH1D* h_centauro_match_multiplicity = nullptr;
};

struct SoCResult {
  double soc_count = std::numeric_limits<double>::quiet_NaN();
  double soc_energy = std::numeric_limits<double>::quiet_NaN();
};

SoCResult compute_soc(const jet_tools::SimpleJet& centauro,
                      const jet_tools::SimpleJet& antikt) {
  SoCResult out;
  const std::size_t centauro_size = std::min(centauro.constituents.size(), centauro.constituent_E.size());
  const std::size_t antikt_size = std::min(antikt.constituents.size(), antikt.constituent_E.size());

  if (centauro_size == 0 || antikt_size == 0) {
    return out;
  }

  std::size_t centauro_total_count = 0;
  std::size_t shared_count = 0;
  double centauro_total_E = 0.0;
  double shared_E = 0.0;

  // Compute count and energy overlap with direct constituent ID matching.
  for (std::size_t ci = 0; ci < centauro_size; ++ci) {
    const int centauro_id = centauro.constituents[ci];
    const double centauro_E = centauro.constituent_E[ci];

    if (!std::isfinite(centauro_E) || centauro_E <= 0.0) {
      continue;
    }

    // Increment denominator (Centauro).
    ++centauro_total_count;
    centauro_total_E += centauro_E;

    // If constituent IDs match, increment numerator (shared).
    for (std::size_t ai = 0; ai < antikt_size; ++ai) {
      const int antikt_id = antikt.constituents[ai];
      if (centauro_id != antikt_id) {
        continue;
      }
      const double antikt_E = antikt.constituent_E[ai];
      if (!std::isfinite(antikt_E) || antikt_E <= 0.0) {
        continue;
      }
      
      ++shared_count;
      shared_E += centauro_E;
      break;
    }
  }

  if (centauro_total_count > 0) {
    out.soc_count =
        static_cast<double>(shared_count) / static_cast<double>(centauro_total_count);
  }
  if (std::isfinite(centauro_total_E) && centauro_total_E > 0.0) {
    out.soc_energy = shared_E / centauro_total_E;
  }
  return out;
}

void process_best_matches(TTree& tree, const EventJets& centauro_jets,
                          const EventJets& antikt_jets, std::size_t max_events,
                          double max_dR, DomainHists& h) {
  BestMatchRow row;
  if (!tree.GetBranch("event") || !tree.GetBranch("centauro_index") ||
      !tree.GetBranch("antikt_index") || !tree.GetBranch("dR")) {
    std::cerr << "[splitting_cf] tree '" << tree.GetName() << "' missing required branches\n";
    std::exit(1);
  }

  tree.SetBranchAddress("event", &row.event);
  tree.SetBranchAddress("centauro_index", &row.centauro_index);
  tree.SetBranchAddress("antikt_index", &row.antikt_index);
  tree.SetBranchAddress("dR", &row.dR);

  const std::size_t entries = static_cast<std::size_t>(tree.GetEntries());
  jet_tools::EventLimitGate event_gate(max_events);
  jet_tools::ProgressTicker progress;
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      const std::string message = std::string("[splitting_cf] best ") + tree.GetName() + " " + 
                                  std::to_string(i) + "/" + std::to_string(entries);
      progress.report(message);
    }
    tree.GetEntry(static_cast<Long64_t>(i));

    if (!event_gate.allow(row.event)) {
      continue;
    }
    if (!std::isfinite(row.dR)) {
      continue;
    }
    if (std::isfinite(max_dR) && row.dR > max_dR) {
      continue;
    }

    // Grab event ids.
    const auto centauro_it = centauro_jets.find(row.event);
    const auto antikt_it = antikt_jets.find(row.event);
    if (centauro_it == centauro_jets.end() || antikt_it == antikt_jets.end()) {
      continue;
    }

    if (row.centauro_index < 0 || row.antikt_index < 0) {
      continue;
    }
    const std::size_t centauro_index = static_cast<std::size_t>(row.centauro_index);
    const std::size_t antikt_index = static_cast<std::size_t>(row.antikt_index);
    if (centauro_index >= centauro_it->second.size() || antikt_index >= antikt_it->second.size()) {
      continue;
    }

    // Grab jet info.
    const auto& centauro = centauro_it->second[centauro_index];
    const auto& antikt = antikt_it->second[antikt_index];
    if (!std::isfinite(centauro.pT) || centauro.pT <= 0.0 || !std::isfinite(antikt.pT)) {
      continue;
    }

    // Transverse momenta ratio.
    const double pT_ratio = antikt.pT / centauro.pT;
    if (!std::isfinite(pT_ratio)) {
      continue;
    }

    // FIXME: switch the antikt.eta to a boosted centauro.eta, antikt.eta probably works as an estimate 
    // for now, would need to store the boosted centauro lab eta before/during matching?
    h.h_ratio->Fill(pT_ratio);
    if (std::isfinite(antikt.eta)) {
      h.h_ratio_eta_lab->Fill(antikt.eta, pT_ratio);
    }
    if (std::isfinite(centauro.eta)) {
      h.h_ratio_eta_breit->Fill(centauro.eta, pT_ratio);
    }

    h.h_nconst_antikt->Fill(static_cast<double>(antikt.n_constituents));
    h.h_nconst_centauro->Fill(static_cast<double>(centauro.n_constituents));

    const SoCResult soc = compute_soc(centauro, antikt);
    // Shared constituent count / centauro total constituent count
    if (std::isfinite(soc.soc_count)) {
      h.h_soc_count->Fill(soc.soc_count);
    }
    // Shared energy amount / centauro total energy
    if (std::isfinite(soc.soc_energy)) {
      h.h_soc_energy->Fill(soc.soc_energy);
    }
  }
  progress.finish(std::string("[splitting_cf] best ") + tree.GetName() + " " +
                  std::to_string(entries) + "/" + std::to_string(entries));
}

// TODO: expand process_all_matches to not only deal with multiplicity.
void process_all_matches(TTree& tree, std::size_t max_events, DomainHists& h) {
  AllMatchRow row;
  if (!tree.GetBranch("event") || !tree.GetBranch("centauro_index") ||
      !tree.GetBranch("antikt_indices")) {
    std::cerr << "[splitting_cf] tree '" << tree.GetName() << "' missing required branches\n";
    std::exit(1);
  }

  tree.SetBranchAddress("event", &row.event);
  tree.SetBranchAddress("centauro_index", &row.centauro_index);
  tree.SetBranchAddress("antikt_indices", &row.antikt_indices);

  const std::size_t entries = static_cast<std::size_t>(tree.GetEntries());
  jet_tools::EventLimitGate event_gate(max_events);
  jet_tools::ProgressTicker progress;
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      const std::string message = std::string("[splitting_cf] all ") +
                                  tree.GetName() + " " + std::to_string(i) +
                                  "/" + std::to_string(entries);
      progress.report(message);
    }
    tree.GetEntry(static_cast<Long64_t>(i));

    if (!event_gate.allow(row.event)) {
      continue;
    }
    if (row.centauro_index < 0 || row.antikt_indices == nullptr) {
      continue;
    }

    const std::size_t multiplicity = row.antikt_indices->size();
    h.h_centauro_match_multiplicity->Fill(static_cast<double>(multiplicity));
  }
  progress.finish(std::string("[splitting_cf] all ") + tree.GetName() + " " +
                  std::to_string(entries) + "/" + std::to_string(entries));
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  TFile f_matches(args.matches_input.c_str(), "READ");
  if (!f_matches.IsOpen()) {
    std::cerr << "[splitting_cf] failed to open matches input " << args.matches_input << "\n";
    return 1;
  }

  TFile f_jets(args.jets_input.c_str(), "READ");
  if (!f_jets.IsOpen()) {
    std::cerr << "[splitting_cf] failed to open jets input " << args.jets_input << "\n";
    return 1;
  }

  //=======================================================================================
  // Grab matches and copy the jet values.
  //=======================================================================================
  auto* t_best_truth = jet_tools::get_required_tree(f_matches, "best_matches_truth", "splitting_cf");
  auto* t_all_truth = jet_tools::get_required_tree(f_matches, "all_matches_truth", "splitting_cf");
  auto* t_best_reco = jet_tools::get_required_tree(f_matches, "best_matches_reco", "splitting_cf");
  auto* t_all_reco = jet_tools::get_required_tree(f_matches, "all_matches_reco", "splitting_cf");

  auto* t_centauro_truth = jet_tools::get_required_tree(f_jets, "BreitFrameCentauroTruth", "splitting_cf");
  auto* t_antikt_truth = jet_tools::get_required_tree(f_jets, "LabFrameAntiktTruth", "splitting_cf");
  auto* t_centauro_reco = jet_tools::get_required_tree(f_jets, "BreitFrameCentauroReco", "splitting_cf");
  auto* t_antikt_reco = jet_tools::get_required_tree(f_jets, "LabFrameAntiktReco", "splitting_cf");

  EventJets centauro_truth;
  EventJets antikt_truth;
  EventJets centauro_reco;
  EventJets antikt_reco;
  jet_tools::JetTreeReadOptions jet_read_options;
  jet_read_options.require_n_constituents = true;
  jet_read_options.require_constituent_vectors = true;
  jet_tools::read_jet_tree(*t_centauro_truth, args.max_events, centauro_truth, "splitting_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_truth, args.max_events, antikt_truth, "splitting_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_reco, args.max_events, centauro_reco, "splitting_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_reco, args.max_events, antikt_reco, "splitting_cf", jet_read_options);

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[splitting_cf] failed to open output " << args.output << "\n";
    return 1;
  }

  //=======================================================================================
  // Make plots.
  //=======================================================================================
  TDirectory* dir_ratio = fout.mkdir("ratio");
  TDirectory* dir_eta_lab = fout.mkdir("eta_lab");
  TDirectory* dir_eta_breit = fout.mkdir("eta_breit");
  TDirectory* dir_nconst = fout.mkdir("nconst");
  TDirectory* dir_match = fout.mkdir("match");
  TDirectory* dir_shared_over_centauro = fout.mkdir("shared_over_centauro");
  if (!dir_ratio || !dir_eta_lab || !dir_eta_breit || !dir_nconst || !dir_match || !dir_shared_over_centauro) {
    std::cerr << "[splitting_cf] failed to create output directories\n";
    return 1;
  }

  TH1D h_ratio_lab_reco("pT_ratio_antikt_centauro_lab_reco",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} (reco);ratio;jets", 120, 0.0, 3.0);
  TH1D h_ratio_lab_truth("pT_ratio_antikt_centauro_lab_truth",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} (truth);ratio;jets", 120, 0.0, 3.0);
  TH2D h_ratio_eta_lab_reco("pT_ratio_antikt_centauro_vs_eta_lab_reco",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} vs #eta_{lab} (reco);#eta_{lab};ratio",
      120, -6.0, 6.0, 120, 0.0, 10.0);
  TH2D h_ratio_eta_lab_truth("pT_ratio_antikt_centauro_vs_eta_lab_truth",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} vs #eta_{lab} (truth);#eta_{lab};ratio",
      120, -6.0, 6.0, 120, 0.0, 10.0);
  TH2D h_ratio_eta_breit_reco("pT_ratio_antikt_centauro_vs_eta_breit_reco",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} vs #eta_{Breit} (reco);#eta_{Breit};ratio",
      120, -6.0, 6.0, 120, 0.0, 10.0);
  TH2D h_ratio_eta_breit_truth("pT_ratio_antikt_centauro_vs_eta_breit_truth",
      "p_{T,lab}^{anti-k_{t}}/p_{T,lab}^{Centauro} vs #eta_{Breit} (truth);#eta_{Breit};ratio",
      120, -6.0, 6.0, 120, 0.0, 10.0);
  TH1D h_nconst_antikt_reco("nconst_antikt_reco", "Nconst (anti-k_{t}, reco);Nconst;jets", 120, 0.0,
      120.0);
  TH1D h_nconst_centauro_truth("nconst_centauro_truth", "Nconst (Centauro, truth);Nconst;jets", 120, 0.0,
      120.0);
  TH1D h_nconst_antikt_truth("nconst_antikt_truth", "Nconst (anti-k_{t}, truth);Nconst;jets", 120, 0.0,
      120.0);
  TH1D h_nconst_centauro_reco("nconst_centauro_reco", "Nconst (Centauro, reco);Nconst;jets", 120, 0.0,
      120.0);
  TH1D h_soc_count_reco("soc_count_reco", "SoC count (reco);SoC count;jets", 60, 0.0, 1.0);
  TH1D h_soc_count_truth("soc_count_truth", "SoC count (truth);SoC count;jets", 60, 0.0, 1.0);
  TH1D h_soc_energy_reco("soc_energy_reco", "SoC energy (reco);SoC energy;jets", 60, 0.0, 1.0);
  TH1D h_soc_energy_truth("soc_energy_truth", "SoC energy (truth);SoC energy;jets", 60, 0.0, 1.0);
  TH1D h_centauro_mult_reco("centauro_match_multiplicity_reco",
      "Centauro match multiplicity (reco);N(anti-k_{t}) per Centauro;Centauro", 10, 0.5, 10.5);
  TH1D h_centauro_mult_truth("centauro_match_multiplicity_truth",
      "Centauro match multiplicity (truth);N(anti-k_{t}) per Centauro;Centauro", 10, 0.5, 10.5);

  DomainHists truth_hists{&h_ratio_lab_truth,
                          &h_ratio_eta_lab_truth,
                          &h_ratio_eta_breit_truth,
                          &h_nconst_antikt_truth,
                          &h_nconst_centauro_truth,
                          &h_soc_count_truth,
                          &h_soc_energy_truth,
                          &h_centauro_mult_truth};
  DomainHists reco_hists{&h_ratio_lab_reco,
                         &h_ratio_eta_lab_reco,
                         &h_ratio_eta_breit_reco,
                         &h_nconst_antikt_reco,
                         &h_nconst_centauro_reco,
                         &h_soc_count_reco,
                         &h_soc_energy_reco,
                         &h_centauro_mult_reco};

  process_best_matches(*t_best_truth, centauro_truth, antikt_truth, args.max_events,
                       args.max_dR, truth_hists);
  process_best_matches(*t_best_reco, centauro_reco, antikt_reco, args.max_events,
                       args.max_dR, reco_hists);
  process_all_matches(*t_all_truth, args.max_events, truth_hists);
  process_all_matches(*t_all_reco, args.max_events, reco_hists);

  auto write_in_dir = [&](TDirectory* dir) {
    if (dir) {
      dir->cd();
    } else {
      fout.cd();
    }
  };

  write_in_dir(dir_ratio);
  h_ratio_lab_reco.Write();
  h_ratio_lab_truth.Write();

  write_in_dir(dir_eta_lab);
  h_ratio_eta_lab_reco.Write();
  h_ratio_eta_lab_truth.Write();

  write_in_dir(dir_eta_breit);
  h_ratio_eta_breit_reco.Write();
  h_ratio_eta_breit_truth.Write();

  write_in_dir(dir_nconst);
  h_nconst_antikt_reco.Write();
  h_nconst_centauro_reco.Write();
  h_nconst_antikt_truth.Write();
  h_nconst_centauro_truth.Write();

  write_in_dir(dir_match);
  h_centauro_mult_reco.Write();
  h_centauro_mult_truth.Write();

  write_in_dir(dir_shared_over_centauro);
  h_soc_count_reco.Write();
  h_soc_count_truth.Write();
  h_soc_energy_reco.Write();
  h_soc_energy_truth.Write();

  fout.Close();
  f_matches.Close();
  f_jets.Close();
  std::cout << "[splitting_cf] wrote " << args.output << "\n";
  return 0;
}
