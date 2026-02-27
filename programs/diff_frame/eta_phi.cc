/*
Run:
./run.sh programs/diff_frame/eta_phi.cc \
  -m data/jets/matched/diff_frame_matches.root \
  -j data/jets/raw/jets_0-999_lab_cuts.root \
  -e data/events/filtered/merged_0-999.root \
  -o data/graphs/diff_frame_eta_phi.root

Select one truth event and one reco event using the highest average matched-pair constituent count 
(so distributions are more visible?), then plot constituent eta-phi distributions for anti-kt (lab)
and Centauro (Breit) jets.

Also select an event pair with the highest number of matched jets.
*/

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Rtypes.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/ReconstructedParticleCollection.h>

#include "jet_tools/include/electron_veto.h"
#include "jet_tools/include/event_progress.h"
#include "jet_tools/include/root_io.h"

namespace {

using EventJets = jet_tools::EventJets;
using BestRowsByEvent = std::unordered_map<ULong64_t, std::vector<jet_tools::CrossFrameBestMatchRow>>;
using AllRowsByEvent = std::unordered_map<ULong64_t, std::vector<jet_tools::CrossFrameAllMatchRow>>;

struct Args {
  std::string matches_input;
  std::string jets_input;
  std::string events_input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
};

struct SelectedEvent {
  bool valid = false;
  ULong64_t event = 0;
  std::size_t pair_count = 0;
  std::size_t centauro_jet_count = 0;
  std::size_t antikt_jet_count = 0;
  double cent_avg = std::numeric_limits<double>::quiet_NaN();
  double antikt_avg = std::numeric_limits<double>::quiet_NaN();
  double score = -std::numeric_limits<double>::infinity();
};

struct DomainSelection {
  SelectedEvent selected;
  BestRowsByEvent rows_by_event;
  std::vector<ULong64_t> event_order;
};

struct DomainHists {
  TH2D* h_antikt = nullptr;
  TH2D* h_centauro = nullptr;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -m MATCHES.root -j JETS.root -e EVENTS.root -o OUTPUT.root\n"
            << "       [--max-events N]\n";
}

//-----------------------------------------------------------------------------------------
// I/O helpers.
//-----------------------------------------------------------------------------------------
Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-m" || v == "--matches-input") && i + 1 < argc) {
      args.matches_input = argv[++i];
    } else if ((v == "-j" || v == "--jets-input") && i + 1 < argc) {
      args.jets_input = argv[++i];
    } else if ((v == "-e" || v == "--events-input") && i + 1 < argc) {
      args.events_input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[eta_phi_cf] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.matches_input.empty() || args.jets_input.empty() ||
      args.events_input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

bool has_collection(const podio::Frame& frame, const std::string& name) {
  for (const auto& n : frame.getAvailableCollections()) {
    if (n == name) {
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------------------
// Compute eta and phi by searching up the particle using an index.
//-----------------------------------------------------------------------------------------
double compute_eta(const edm4eic::ReconstructedParticleCollection& particles, int index) {
  if (index < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const std::size_t ui = static_cast<std::size_t>(index);
  if (ui >= particles.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto mom = particles[ui].getMomentum();
  const double pT = std::hypot(mom.x, mom.y);
  if (!std::isfinite(pT) || pT <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double eta = std::asinh(mom.z / pT);
  return std::isfinite(eta) ? eta : std::numeric_limits<double>::quiet_NaN();
}

double compute_phi(const edm4eic::ReconstructedParticleCollection& particles, int index) {
  if (index < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const std::size_t ui = static_cast<std::size_t>(index);
  if (ui >= particles.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto mom = particles[ui].getMomentum();
  if (!std::isfinite(mom.x) || !std::isfinite(mom.y)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double phi = std::atan2(mom.y, mom.x);
  return std::isfinite(phi) ? phi : std::numeric_limits<double>::quiet_NaN();
}

//-----------------------------------------------------------------------------------------
// Histogram fillers.
//-----------------------------------------------------------------------------------------
// Let the color density represent the raw number of particles at (eta,phi)
void fill_constituent_eta_phi_count(const edm4eic::ReconstructedParticleCollection& particles,
                                    const jet_tools::SimpleJet& jet, TH2D& hist) {
  for (int index : jet.constituent_indices) {
    const double eta = compute_eta(particles, index);
    const double phi = compute_phi(particles, index);
    if (!std::isfinite(eta) || !std::isfinite(phi)) {
      continue;
    }
    hist.Fill(eta, phi);
  }
}

// Let the color density represent the energy of the particles at (eta,phi).
void fill_constituent_eta_phi_energy(const edm4eic::ReconstructedParticleCollection& particles,
                                     const jet_tools::SimpleJet& jet, TH2D& hist) {
  const std::size_t entries = std::min(jet.constituent_indices.size(), jet.constituent_E.size());

  for (std::size_t i = 0; i < entries; ++i) {
    const double energy = jet.constituent_E[i];
    if (!std::isfinite(energy) || energy <= 0.0) {
      continue;
    }
    const double eta = compute_eta(particles, jet.constituent_indices[i]);
    const double phi = compute_phi(particles, jet.constituent_indices[i]);
    if (!std::isfinite(eta) || !std::isfinite(phi)) {
      continue;
    }
    hist.Fill(eta, phi, energy);
  }
}

// Wrap both the count-based and energy-based histograms.
void fill_constituent_eta_phi_both(const edm4eic::ReconstructedParticleCollection& particles,
                                   const jet_tools::SimpleJet& jet, DomainHists& count_hists,
                                   DomainHists& energy_hists, bool is_antikt) {
  TH2D* const count_hist = is_antikt ? count_hists.h_antikt : count_hists.h_centauro;
  TH2D* const energy_hist = is_antikt ? energy_hists.h_antikt : energy_hists.h_centauro;
  if (count_hist) {
    fill_constituent_eta_phi_count(particles, jet, *count_hist);
  }
  if (energy_hist) {
    fill_constituent_eta_phi_energy(particles, jet, *energy_hist);
  }
}

// Sort match entries by event number.
DomainSelection group_best_rows_by_event(const std::vector<jet_tools::CrossFrameBestMatchRow>& rows) {
  DomainSelection out;
  out.rows_by_event.reserve(rows.size());

  for (const auto& row : rows) {
    auto it = out.rows_by_event.find(row.event);
    if (it == out.rows_by_event.end()) {
      out.event_order.push_back(row.event);
      it = out.rows_by_event.emplace(row.event, std::vector<jet_tools::CrossFrameBestMatchRow>{}).first;
    }
    it->second.push_back(row);
  }

  return out;
}

// Find the event where jets have the greatest number of constituents.
DomainSelection build_domain_selection_by_constituent_count(
    const std::vector<jet_tools::CrossFrameBestMatchRow>& rows,
    const EventJets& centauro_jets, const EventJets& antikt_jets) {
  DomainSelection out = group_best_rows_by_event(rows);
  jet_tools::ProgressTicker progress;
  const std::size_t entries = out.event_order.size();

  // Loop through candidate events.
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(std::string("[eta_phi_cf] score events ") + std::to_string(i) + "/" 
          + std::to_string(entries));
    }

    const ULong64_t event_id = out.event_order[i];
    const auto rows_it = out.rows_by_event.find(event_id);
    if (rows_it == out.rows_by_event.end()) {
      continue;
    }

    // Find Centauro and antikt jets per event.
    const auto cent_it = centauro_jets.find(event_id);
    const auto antikt_it = antikt_jets.find(event_id);
    if (cent_it == centauro_jets.end() || antikt_it == antikt_jets.end()) {
      continue;
    }

    double cent_total = 0.0;
    double antikt_total = 0.0;
    std::size_t denom = 0;

    // Loop through all matched pairs and read the jet indices.
    for (const auto& row : rows_it->second) {
      if (!std::isfinite(row.dR) || row.centauro_index < 0 || row.antikt_index < 0) {
        continue;
      }
      const std::size_t cent_idx = static_cast<std::size_t>(row.centauro_index);
      const std::size_t antikt_idx = static_cast<std::size_t>(row.antikt_index);
      if (cent_idx >= cent_it->second.size() || antikt_idx >= antikt_it->second.size()) {
        continue;
      }

      // Accumulate constituent counts for matched jets.
      cent_total += static_cast<double>(cent_it->second[cent_idx].n_constituents);
      antikt_total += static_cast<double>(antikt_it->second[antikt_idx].n_constituents);

      // Increment denominator for matched-pair averages.
      ++denom;
    }
    if (denom == 0) {
      continue;
    }

    const double cent_avg = cent_total / static_cast<double>(denom);
    const double antikt_avg = antikt_total / static_cast<double>(denom);

    // Use both the centauro average and anti-kt average.
    const double score = 0.5 * (cent_avg + antikt_avg);
    if (!std::isfinite(score)) {
      continue;
    }

    // Only keep if greater value than previous best score.
    if (!out.selected.valid || score > out.selected.score) {
      out.selected.valid = true;
      out.selected.event = event_id;
      out.selected.pair_count = denom;
      out.selected.cent_avg = cent_avg;
      out.selected.antikt_avg = antikt_avg;
      out.selected.centauro_jet_count = cent_it->second.size();
      out.selected.antikt_jet_count = antikt_it->second.size();
      out.selected.score = score;
    }
  }
  progress.finish(std::string("[eta_phi_cf] score events ") + std::to_string(entries) + "/" 
      + std::to_string(entries));
  return out;
}

// Find the event with the most number of jets.
DomainSelection build_domain_selection_by_jet_count(
    const std::vector<jet_tools::CrossFrameBestMatchRow>& rows,
    const EventJets& centauro_jets, const EventJets& antikt_jets) {
  // Reuse the per-event grouping so event iteration stays consistent with input order.
  DomainSelection out = group_best_rows_by_event(rows);
  jet_tools::ProgressTicker progress;
  const std::size_t entries = out.event_order.size();

  // Score each candidate event by total jet multiplicity (Centauro + anti-kt).
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(std::string("[eta_phi_cf] score jet-count events ") +
                      std::to_string(i) + "/" + std::to_string(entries));
    }

    const ULong64_t event_id = out.event_order[i];
    const auto cent_it = centauro_jets.find(event_id);
    const auto antikt_it = antikt_jets.find(event_id);

    if (cent_it == centauro_jets.end() || antikt_it == antikt_jets.end()) {
      continue;
    }

    // Count only valid matched pairs so selected events have usable match content.
    std::size_t valid_pairs = 0;
    const auto rows_it = out.rows_by_event.find(event_id);
    if (rows_it != out.rows_by_event.end()) {
      for (const auto& row : rows_it->second) {
        if (!std::isfinite(row.dR) || row.centauro_index < 0 || row.antikt_index < 0) {
          continue;
        }
        const std::size_t cent_idx = static_cast<std::size_t>(row.centauro_index);
        const std::size_t antikt_idx = static_cast<std::size_t>(row.antikt_index);
        if (cent_idx >= cent_it->second.size() || antikt_idx >= antikt_it->second.size()) {
          continue;
        }

        ++valid_pairs;
      }
    }
    // Skip events that have no usable best-match rows.
    if (valid_pairs == 0) {
      continue;
    }

    // Jet-count score for this event: N_centauro + N_antikt.
    const std::size_t centauro_jet_count = cent_it->second.size();
    const std::size_t antikt_jet_count = antikt_it->second.size();
    const double score = static_cast<double>(centauro_jet_count + antikt_jet_count);
    if (!std::isfinite(score)) {
      continue;
    }

    if (!out.selected.valid || score > out.selected.score) {
      // This selector is jet-count based, so constituent-average fields are not defined.
      out.selected.valid = true;
      out.selected.event = event_id;
      out.selected.pair_count = valid_pairs;
      out.selected.centauro_jet_count = centauro_jet_count;
      out.selected.antikt_jet_count = antikt_jet_count;
      out.selected.cent_avg = std::numeric_limits<double>::quiet_NaN();
      out.selected.antikt_avg = std::numeric_limits<double>::quiet_NaN();
      out.selected.score = score;
    }
  }
  progress.finish(std::string("[eta_phi_cf] score jet-count events ") +
                  std::to_string(entries) + "/" + std::to_string(entries));
  return out;
}

std::optional<podio::Frame> load_event_frame(podio::ROOTReader& reader,
                                             ULong64_t event_id) {
  const std::size_t entries_total = reader.getEntries("events");
  const std::size_t index = static_cast<std::size_t>(event_id);
  if (index >= entries_total) {
    return std::nullopt;
  }
  auto data = reader.readEntry("events", index);
  if (!data) {
    return std::nullopt;
  }
  return podio::Frame(std::move(data));
}

AllRowsByEvent group_all_rows_by_event(
    const std::vector<jet_tools::CrossFrameAllMatchRow>& rows) {
  AllRowsByEvent out;
  out.reserve(rows.size());
  for (const auto& row : rows) {
    out[row.event].push_back(row);
  }
  return out;
}

// Fill histograms for the best-match rows (e.g. binary pairs).
void fill_best_domain_hists_for_event(const char* tag, const podio::Frame& frame,
                                      const SelectedEvent& selected,
                                      const BestRowsByEvent& rows_by_event,
                                      const EventJets& centauro_jets,
                                      const EventJets& antikt_jets,
                                      const std::string& antikt_collection_name,
                                      const std::string& centauro_collection_name,
                                      DomainHists& count_hists,
                                      DomainHists& energy_hists) {
  if (!selected.valid) {
    return;
  }
  if ((!count_hists.h_antikt && !energy_hists.h_antikt) ||
      (!count_hists.h_centauro && !energy_hists.h_centauro)) {
    return;
  }
  if (!has_collection(frame, antikt_collection_name) ||
      !has_collection(frame, centauro_collection_name)) {
    std::cerr << "[eta_phi_cf] missing required collections for " << tag
              << " event " << selected.event << "\n";
    return;
  }

  const auto& antikt_particles = frame.get<edm4eic::ReconstructedParticleCollection>(antikt_collection_name);
  const auto& centauro_particles = frame.get<edm4eic::ReconstructedParticleCollection>(centauro_collection_name);

  const auto rows_it = rows_by_event.find(selected.event);
  const auto cent_it = centauro_jets.find(selected.event);
  const auto antikt_it = antikt_jets.find(selected.event);
  if (rows_it == rows_by_event.end() || cent_it == centauro_jets.end() || antikt_it == antikt_jets.end()) {
    return;
  }

  const std::size_t entries = rows_it->second.size();
  jet_tools::ProgressTicker progress;

  // Loop through matched pairs in the selected event.
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(std::string("[eta_phi_cf] fill ") + tag + " pairs " +
                      std::to_string(i) + "/" + std::to_string(entries));
    }
    const auto& row = rows_it->second[i];
    if (row.centauro_index < 0 || row.antikt_index < 0) {
      continue;
    }

    const std::size_t cent_idx = static_cast<std::size_t>(row.centauro_index);
    const std::size_t antikt_idx = static_cast<std::size_t>(row.antikt_index);
    if (cent_idx >= cent_it->second.size() || antikt_idx >= antikt_it->second.size()) {
      continue;
    }

    fill_constituent_eta_phi_both(
        antikt_particles, antikt_it->second[antikt_idx], count_hists, energy_hists, true);
    fill_constituent_eta_phi_both(
        centauro_particles, cent_it->second[cent_idx], count_hists, energy_hists, false);
  }
  progress.finish(std::string("[eta_phi_cf] fill ") + tag + " pairs " +
                  std::to_string(entries) + "/" + std::to_string(entries));
}

// Fill histograms for (centauro, all antikt within dR).
void fill_all_domain_hists_for_event(const char* tag, const podio::Frame& frame,
                                     const SelectedEvent& selected,
                                     const AllRowsByEvent& rows_by_event,
                                     const EventJets& centauro_jets,
                                     const EventJets& antikt_jets,
                                     const std::string& antikt_collection_name,
                                     const std::string& centauro_collection_name,
                                     DomainHists& count_hists,
                                     DomainHists& energy_hists) {
  if (!selected.valid) {
    return;
  }
  if ((!count_hists.h_antikt && !energy_hists.h_antikt) ||
      (!count_hists.h_centauro && !energy_hists.h_centauro)) {
    return;
  }
  if (!has_collection(frame, antikt_collection_name) ||
      !has_collection(frame, centauro_collection_name)) {
    std::cerr << "[eta_phi_cf] missing required collections for " << tag
              << " event " << selected.event << "\n";
    return;
  }

  // Grab particles.
  const auto& antikt_particles = frame.get<edm4eic::ReconstructedParticleCollection>(antikt_collection_name);
  const auto& centauro_particles = frame.get<edm4eic::ReconstructedParticleCollection>(centauro_collection_name);

  // Find per-event data.
  const auto rows_it = rows_by_event.find(selected.event);
  const auto cent_it = centauro_jets.find(selected.event);
  const auto antikt_it = antikt_jets.find(selected.event);
  if (rows_it == rows_by_event.end() || cent_it == centauro_jets.end() ||
      antikt_it == antikt_jets.end()) {
    return;
  }

  const std::size_t entries = rows_it->second.size();
  jet_tools::ProgressTicker progress;

  // Loop through all match rows.
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      progress.report(std::string("[eta_phi_cf] fill ") + tag + " rows " +
                      std::to_string(i) + "/" + std::to_string(entries));
    }
    const auto& row = rows_it->second[i];
    if (row.centauro_index < 0) {
      continue;
    }

    const std::size_t cent_idx = static_cast<std::size_t>(row.centauro_index);
    if (cent_idx >= cent_it->second.size()) {
      continue;
    }

    fill_constituent_eta_phi_both(
        centauro_particles, cent_it->second[cent_idx], count_hists, energy_hists, false);

    for (int antikt_index : row.antikt_indices) {
      if (antikt_index < 0) {
        continue;
      }
      const std::size_t antikt_idx = static_cast<std::size_t>(antikt_index);
      if (antikt_idx >= antikt_it->second.size()) {
        continue;
      }
      fill_constituent_eta_phi_both(
          antikt_particles, antikt_it->second[antikt_idx], count_hists, energy_hists, true);
    }
  }
  progress.finish(std::string("[eta_phi_cf] fill ") + tag + " rows " +
                  std::to_string(entries) + "/" + std::to_string(entries));
}

void write_leaf_plots(TDirectory* dir, TH2D& h_antikt, TH2D& h_centauro,
                      const char* canvas_name, const char* canvas_title) {
  if (!dir) {
    return;
  }
  dir->cd();
  h_antikt.Write();
  h_centauro.Write();
  TCanvas canvas(canvas_name, canvas_title, 1200, 500);
  canvas.Divide(2, 1);
  canvas.cd(1);
  h_antikt.Draw("COLZ");
  canvas.cd(2);
  h_centauro.Draw("COLZ");
  canvas.Write();
}

void write_selection_outputs(
    TFile& fout, const char* selection_dir_name,
    const SelectedEvent& truth_selected, const SelectedEvent& reco_selected,
    const BestRowsByEvent& truth_best_rows_by_event,
    const BestRowsByEvent& reco_best_rows_by_event,
    const AllRowsByEvent& truth_all_rows_by_event,
    const AllRowsByEvent& reco_all_rows_by_event,
    const EventJets& centauro_truth, const EventJets& antikt_truth,
    const EventJets& centauro_reco, const EventJets& antikt_reco,
    const std::optional<podio::Frame>& truth_frame,
    const std::optional<podio::Frame>& reco_frame,
    const jet_tools::JetToolsCollections& collections) {
  constexpr double min_eta = -6.0;
  constexpr double max_eta = 6.0;
  constexpr double min_phi = -3.2;
  constexpr double max_phi = 3.2;
  constexpr int n_eta_bins = 240;
  constexpr int n_phi_bins = 128;

  const std::string prefix = std::string(selection_dir_name);
  auto name_of = [&](const char* base) {
    return prefix + "_" + base;
  };

  //---------------------------------------------------------------------------------------
  // Histogram declarations
  //---------------------------------------------------------------------------------------
  TH2D h_count_truth_best_antikt(name_of("antikt_eta_phi_count_truth_best").c_str(),
      "Count-weighted truth anti-k_{t} constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_truth_best_centauro(name_of("centauro_eta_phi_count_truth_best").c_str(),
      "Count-weighted truth Centauro constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_reco_best_antikt(name_of("antikt_eta_phi_count_reco_best").c_str(),
      "Count-weighted reco anti-k_{t} constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_reco_best_centauro(name_of("centauro_eta_phi_count_reco_best").c_str(),
      "Count-weighted reco Centauro constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_truth_all_antikt(name_of("antikt_eta_phi_count_truth_all").c_str(),
      "Count-weighted truth anti-k_{t} constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_truth_all_centauro(name_of("centauro_eta_phi_count_truth_all").c_str(),
      "Count-weighted truth Centauro constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_reco_all_antikt(name_of("antikt_eta_phi_count_reco_all").c_str(),
      "Count-weighted reco anti-k_{t} constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_count_reco_all_centauro(name_of("centauro_eta_phi_count_reco_all").c_str(),
      "Count-weighted reco Centauro constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);

  TH2D h_energy_truth_best_antikt(name_of("antikt_eta_phi_energy_truth_best").c_str(),
      "Energy-weighted truth anti-k_{t} constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_truth_best_centauro(name_of("centauro_eta_phi_energy_truth_best").c_str(),
      "Energy-weighted truth Centauro constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_reco_best_antikt(name_of("antikt_eta_phi_energy_reco_best").c_str(),
      "Energy-weighted reco anti-k_{t} constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_reco_best_centauro(name_of("centauro_eta_phi_energy_reco_best").c_str(),
      "Energy-weighted reco Centauro constituent #eta-#phi (best rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_truth_all_antikt(name_of("antikt_eta_phi_energy_truth_all").c_str(),
      "Energy-weighted truth anti-k_{t} constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_truth_all_centauro(name_of("centauro_eta_phi_energy_truth_all").c_str(),
      "Energy-weighted truth Centauro constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_reco_all_antikt(name_of("antikt_eta_phi_energy_reco_all").c_str(),
      "Energy-weighted reco anti-k_{t} constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);
  TH2D h_energy_reco_all_centauro(name_of("centauro_eta_phi_energy_reco_all").c_str(),
      "Energy-weighted reco Centauro constituent #eta-#phi (all rows);#eta;#phi",
      n_eta_bins, min_eta, max_eta, n_phi_bins, min_phi, max_phi);

  TH2D* all_hists[] = {
      &h_count_truth_best_antikt, &h_count_truth_best_centauro,
      &h_count_reco_best_antikt, &h_count_reco_best_centauro,
      &h_count_truth_all_antikt, &h_count_truth_all_centauro,
      &h_count_reco_all_antikt, &h_count_reco_all_centauro,
      &h_energy_truth_best_antikt, &h_energy_truth_best_centauro,
      &h_energy_reco_best_antikt, &h_energy_reco_best_centauro,
      &h_energy_truth_all_antikt, &h_energy_truth_all_centauro,
      &h_energy_reco_all_antikt, &h_energy_reco_all_centauro};
  for (TH2D* hist : all_hists) {
    hist->SetStats(0);
  }

  DomainHists count_truth_best{&h_count_truth_best_antikt, &h_count_truth_best_centauro};
  DomainHists count_reco_best{&h_count_reco_best_antikt, &h_count_reco_best_centauro};
  DomainHists count_truth_all{&h_count_truth_all_antikt, &h_count_truth_all_centauro};
  DomainHists count_reco_all{&h_count_reco_all_antikt, &h_count_reco_all_centauro};
  DomainHists energy_truth_best{&h_energy_truth_best_antikt, &h_energy_truth_best_centauro};
  DomainHists energy_reco_best{&h_energy_reco_best_antikt, &h_energy_reco_best_centauro};
  DomainHists energy_truth_all{&h_energy_truth_all_antikt, &h_energy_truth_all_centauro};
  DomainHists energy_reco_all{&h_energy_reco_all_antikt, &h_energy_reco_all_centauro};

  //---------------------------------------------------------------------------------------
  // Fill Histograms
  //---------------------------------------------------------------------------------------
  if (truth_frame.has_value()) {
    fill_best_domain_hists_for_event(
        "truth_best", *truth_frame, truth_selected, truth_best_rows_by_event,
        centauro_truth, antikt_truth, "GeneratedParticles",
        "GeneratedBreitFrameParticles", count_truth_best, energy_truth_best);
    fill_all_domain_hists_for_event(
        "truth_all", *truth_frame, truth_selected, truth_all_rows_by_event,
        centauro_truth, antikt_truth, "GeneratedParticles",
        "GeneratedBreitFrameParticles", count_truth_all, energy_truth_all);
  }
  if (reco_frame.has_value()) {
    fill_best_domain_hists_for_event(
        "reco_best", *reco_frame, reco_selected, reco_best_rows_by_event,
        centauro_reco, antikt_reco, collections.reco_lab_collection,
        "ReconstructedBreitFrameParticles", count_reco_best, energy_reco_best);
    fill_all_domain_hists_for_event(
        "reco_all", *reco_frame, reco_selected, reco_all_rows_by_event,
        centauro_reco, antikt_reco, collections.reco_lab_collection,
        "ReconstructedBreitFrameParticles", count_reco_all, energy_reco_all);
  }

  //---------------------------------------------------------------------------------------
  // Write Outputs
  //---------------------------------------------------------------------------------------
  TDirectory* dir_selection = fout.mkdir(selection_dir_name);
  TDirectory* dir_count = dir_selection ? dir_selection->mkdir("count") : nullptr;
  TDirectory* dir_energy = dir_selection ? dir_selection->mkdir("energy") : nullptr;
  TDirectory* dir_count_truth_best = dir_count ? dir_count->mkdir("truth_best") : nullptr;
  TDirectory* dir_count_reco_best = dir_count ? dir_count->mkdir("reco_best") : nullptr;
  TDirectory* dir_count_truth_all = dir_count ? dir_count->mkdir("truth_all") : nullptr;
  TDirectory* dir_count_reco_all = dir_count ? dir_count->mkdir("reco_all") : nullptr;
  TDirectory* dir_energy_truth_best = dir_energy ? dir_energy->mkdir("truth_best") : nullptr;
  TDirectory* dir_energy_reco_best = dir_energy ? dir_energy->mkdir("reco_best") : nullptr;
  TDirectory* dir_energy_truth_all = dir_energy ? dir_energy->mkdir("truth_all") : nullptr;
  TDirectory* dir_energy_reco_all = dir_energy ? dir_energy->mkdir("reco_all") : nullptr;

  write_leaf_plots(dir_count_truth_best, h_count_truth_best_antikt, h_count_truth_best_centauro,
                   name_of("c_count_truth_best_eta_phi").c_str(),
                   "count truth best eta-phi");
  write_leaf_plots(dir_count_reco_best, h_count_reco_best_antikt, h_count_reco_best_centauro,
                   name_of("c_count_reco_best_eta_phi").c_str(),
                   "count reco best eta-phi");
  write_leaf_plots(dir_count_truth_all, h_count_truth_all_antikt, h_count_truth_all_centauro,
                   name_of("c_count_truth_all_eta_phi").c_str(),
                   "count truth all eta-phi");
  write_leaf_plots(dir_count_reco_all, h_count_reco_all_antikt, h_count_reco_all_centauro,
                   name_of("c_count_reco_all_eta_phi").c_str(),
                   "count reco all eta-phi");
  write_leaf_plots(dir_energy_truth_best, h_energy_truth_best_antikt, h_energy_truth_best_centauro,
                   name_of("c_energy_truth_best_eta_phi").c_str(),
                   "energy truth best eta-phi");
  write_leaf_plots(dir_energy_reco_best, h_energy_reco_best_antikt, h_energy_reco_best_centauro,
                   name_of("c_energy_reco_best_eta_phi").c_str(),
                   "energy reco best eta-phi");
  write_leaf_plots(dir_energy_truth_all, h_energy_truth_all_antikt, h_energy_truth_all_centauro,
                   name_of("c_energy_truth_all_eta_phi").c_str(),
                   "energy truth all eta-phi");
  write_leaf_plots(dir_energy_reco_all, h_energy_reco_all_antikt, h_energy_reco_all_centauro,
                   name_of("c_energy_reco_all_eta_phi").c_str(),
                   "energy reco all eta-phi");
}

void print_constituent_selection_summary(const char* tag, const SelectedEvent& sel) {
  if (!sel.valid) {
    std::cout << "[eta_phi_cf] " << tag << " no valid matched event found\n";
    return;
  }
  std::cout << "[eta_phi_cf] " << tag
            << " best_event=" << sel.event
            << " pairs=" << sel.pair_count
            << " cent_avg=" << sel.cent_avg
            << " antikt_avg=" << sel.antikt_avg
            << " score=" << sel.score << "\n";
}

void print_jet_count_selection_summary(const char* tag, const SelectedEvent& sel) {
  if (!sel.valid) {
    std::cout << "[eta_phi_cf] " << tag << " no valid matched event found\n";
    return;
  }
  std::cout << "[eta_phi_cf] " << tag
            << " best_event=" << sel.event
            << " pairs=" << sel.pair_count
            << " centauro_jet_count=" << sel.centauro_jet_count
            << " antikt_jet_count=" << sel.antikt_jet_count
            << " total_jets=" << (sel.centauro_jet_count + sel.antikt_jet_count)
            << " score=" << sel.score << "\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  TFile f_matches(args.matches_input.c_str(), "READ");
  if (!f_matches.IsOpen()) {
    std::cerr << "[eta_phi_cf] failed to open matches input " << args.matches_input << "\n";
    return 1;
  }

  TFile f_jets(args.jets_input.c_str(), "READ");
  if (!f_jets.IsOpen()) {
    std::cerr << "[eta_phi_cf] failed to open jets input " << args.jets_input << "\n";
    return 1;
  }

  //---------------------------------------------------------------------------------------
  // Load Trees And Rows
  //---------------------------------------------------------------------------------------
  auto* t_best_truth = jet_tools::get_required_tree(f_matches, "best_matches_truth", "eta_phi_cf");
  auto* t_all_truth = jet_tools::get_required_tree(f_matches, "all_matches_truth", "eta_phi_cf");
  auto* t_best_reco = jet_tools::get_required_tree(f_matches, "best_matches_reco", "eta_phi_cf");
  auto* t_all_reco = jet_tools::get_required_tree(f_matches, "all_matches_reco", "eta_phi_cf");
  auto* t_centauro_truth = jet_tools::get_required_tree(f_jets, "BreitFrameCentauroTruth", "eta_phi_cf");
  auto* t_antikt_truth = jet_tools::get_required_tree(f_jets, "LabFrameAntiktTruth", "eta_phi_cf");
  auto* t_centauro_reco = jet_tools::get_required_tree(f_jets, "BreitFrameCentauroReco", "eta_phi_cf");
  auto* t_antikt_reco = jet_tools::get_required_tree(f_jets, "LabFrameAntiktReco", "eta_phi_cf");

  const auto best_truth_rows = jet_tools::read_cf_best_match_tree(*t_best_truth, args.max_events, "eta_phi_cf");
  const auto all_truth_rows = jet_tools::read_cf_all_match_tree(*t_all_truth, args.max_events, "eta_phi_cf");
  const auto best_reco_rows = jet_tools::read_cf_best_match_tree(*t_best_reco, args.max_events, "eta_phi_cf");
  const auto all_reco_rows = jet_tools::read_cf_all_match_tree(*t_all_reco, args.max_events, "eta_phi_cf");

  EventJets centauro_truth;
  EventJets antikt_truth;
  EventJets centauro_reco;
  EventJets antikt_reco;
  jet_tools::JetTreeReadOptions jet_read_options;
  jet_read_options.require_n_constituents = true;
  jet_read_options.require_constituent_vectors = true;
  jet_tools::read_jet_tree(*t_centauro_truth, args.max_events, centauro_truth, "eta_phi_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_truth, args.max_events, antikt_truth, "eta_phi_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_centauro_reco, args.max_events, centauro_reco, "eta_phi_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_antikt_reco, args.max_events, antikt_reco, "eta_phi_cf", jet_read_options);

  const DomainSelection truth_sel_constituent =
      build_domain_selection_by_constituent_count(best_truth_rows, centauro_truth, antikt_truth);
  const DomainSelection reco_sel_constituent =
      build_domain_selection_by_constituent_count(best_reco_rows, centauro_reco, antikt_reco);
  const DomainSelection truth_sel_jet_count =
      build_domain_selection_by_jet_count(best_truth_rows, centauro_truth, antikt_truth);
  const DomainSelection reco_sel_jet_count =
      build_domain_selection_by_jet_count(best_reco_rows, centauro_reco, antikt_reco);
  const AllRowsByEvent all_truth_rows_by_event = group_all_rows_by_event(all_truth_rows);
  const AllRowsByEvent all_reco_rows_by_event = group_all_rows_by_event(all_reco_rows);

  print_constituent_selection_summary("truth highest_constituent_count", truth_sel_constituent.selected);
  print_constituent_selection_summary("reco highest_constituent_count", reco_sel_constituent.selected);
  print_jet_count_selection_summary("truth highest_jet_count", truth_sel_jet_count.selected);
  print_jet_count_selection_summary("reco highest_jet_count", reco_sel_jet_count.selected);

  podio::ROOTReader reader;
  reader.openFile(args.events_input);

  const std::optional<podio::Frame> truth_constituent_frame =
      truth_sel_constituent.selected.valid
          ? load_event_frame(reader, truth_sel_constituent.selected.event)
          : std::nullopt;
  const std::optional<podio::Frame> reco_constituent_frame =
      reco_sel_constituent.selected.valid
          ? load_event_frame(reader, reco_sel_constituent.selected.event)
          : std::nullopt;
  const std::optional<podio::Frame> truth_jet_count_frame =
      truth_sel_jet_count.selected.valid
          ? load_event_frame(reader, truth_sel_jet_count.selected.event)
          : std::nullopt;
  const std::optional<podio::Frame> reco_jet_count_frame =
      reco_sel_jet_count.selected.valid
          ? load_event_frame(reader, reco_sel_jet_count.selected.event)
          : std::nullopt;

  if (truth_sel_constituent.selected.valid && !truth_constituent_frame.has_value()) {
    std::cerr << "[eta_phi_cf] failed to read truth selected event "
              << truth_sel_constituent.selected.event
              << " (highest_constituent_count)\n";
  }
  if (reco_sel_constituent.selected.valid && !reco_constituent_frame.has_value()) {
    std::cerr << "[eta_phi_cf] failed to read reco selected event "
              << reco_sel_constituent.selected.event
              << " (highest_constituent_count)\n";
  }
  if (truth_sel_jet_count.selected.valid && !truth_jet_count_frame.has_value()) {
    std::cerr << "[eta_phi_cf] failed to read truth selected event "
              << truth_sel_jet_count.selected.event
              << " (highest_jet_count)\n";
  }
  if (reco_sel_jet_count.selected.valid && !reco_jet_count_frame.has_value()) {
    std::cerr << "[eta_phi_cf] failed to read reco selected event "
              << reco_sel_jet_count.selected.event
              << " (highest_jet_count)\n";
  }

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[eta_phi_cf] failed to open output " << args.output << "\n";
    return 1;
  }

  const jet_tools::JetToolsCollections collections;
  write_selection_outputs(
      fout, "highest_constituent_count", truth_sel_constituent.selected,
      reco_sel_constituent.selected, truth_sel_constituent.rows_by_event,
      reco_sel_constituent.rows_by_event, all_truth_rows_by_event, all_reco_rows_by_event,
      centauro_truth, antikt_truth, centauro_reco, antikt_reco,
      truth_constituent_frame, reco_constituent_frame, collections);
  write_selection_outputs(
      fout, "highest_jet_count", truth_sel_jet_count.selected,
      reco_sel_jet_count.selected, truth_sel_jet_count.rows_by_event,
      reco_sel_jet_count.rows_by_event, all_truth_rows_by_event, all_reco_rows_by_event,
      centauro_truth, antikt_truth, centauro_reco, antikt_reco,
      truth_jet_count_frame, reco_jet_count_frame, collections);

  fout.Close();
  f_matches.Close();
  f_jets.Close();
  std::cout << "[eta_phi_cf] wrote " << args.output << "\n";
  return 0;
}
