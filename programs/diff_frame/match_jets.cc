/* 
Run:
./run.sh programs/diff_frame/match_jets.cc \
  -j data/jets/raw/jets_0-999_lab_cuts.root -e data/events/filtered/merged_0-999.root \
  -o data/jets/matched/diff_frame_matches.root --match-dR 1

./run.sh programs/diff_frame/match_jets.cc \
  -j data/jets/raw/jets_0-999_lab_cuts.root -e data/events/filtered/merged_0-999.root \
  -o data/jets/matched/diff_frame_matches.root --max-events 1000 --match-dR 0.3
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Math/Boost.h>
#include <Math/Rotation3D.h>
#include <Math/Vector3D.h>
#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/InclusiveKinematicsCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include "jet_tools/include/beam_helpers.h"
#include "jet_tools/include/electron_veto.h"
#include "jet_tools/include/math_helpers.h"
#include "jet_tools/include/root_io.h"
#include "jet_tools/include/transform_breit.h"
#include "jet_tools/include/types.h"

namespace {

using P4 = jet_tools::P4;
using Vec3 = ROOT::Math::XYZVector;

constexpr double kElectronMass = 0.000511;
constexpr double kProtonMass = 0.938272;
constexpr double kCrossingAngle = -0.025;
const std::array<double, 3> kElectronNominalPz{-5.0, -10.0, -18.0};
const std::array<double, 3> kProtonNominalPz{41.0, 100.0, 275.0};

struct Args {
  std::string jets_input;
  std::string events_input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double match_dR = 0.4;
};

struct DomainEventJets {
  std::vector<jet_tools::SimpleJet> centauro_breit;
  std::vector<jet_tools::SimpleJet> antikt_lab;
};

// <event_id, jet_data>
using DomainMap = std::unordered_map<ULong64_t, DomainEventJets>;

void usage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " -j CLUSTERED.root -e EVENTS.root -o OUTPUT.root [--max-events N]\n"
            << "       [--match-dR VAL]\n";
}

Args parse_args(int argc, char *argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-j" || v == "--jets-input") && i + 1 < argc) {
      args.jets_input = argv[++i];
    } else if ((v == "-e" || v == "--events-input") && i + 1 < argc) {
      args.events_input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[match_jets_cf] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--match-dR" && i + 1 < argc) {
      args.match_dR = std::atof(argv[++i]);
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.jets_input.empty() || args.events_input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!std::isfinite(args.match_dR) || args.match_dR <= 0.0) {
    std::cerr << "[match_jets_cf] match-dR must be finite and > 0\n";
    std::exit(1);
  }
  return args;
}

bool has_collection(const podio::Frame &frame, const std::string &name) {
  for (const auto &n : frame.getAvailableCollections()) {
    if (n == name) {
      return true;
    }
  }
  return false;
}

bool build_transform(const podio::Frame &frame, bool use_reco,
                     const jet_tools::JetToolsCollections &collections,
                     const jet_tools::ElectronVetoCuts &cuts,
                     ROOT::Math::Boost &boost, ROOT::Math::Rotation3D &rot) {
  if (!has_collection(frame, "MCParticles")) {
    return false;
  }
  const auto &mc = frame.get<edm4hep::MCParticleCollection>("MCParticles");
  if (mc.empty()) {
    return false;
  }

  P4 p_in_raw;
  P4 k_in_raw;
  if (!jet_tools::find_beam(mc, 2212, p_in_raw) ||
      !jet_tools::find_beam(mc, 11, k_in_raw)) {
    return false;
  }

  const Vec3 p_vec(p_in_raw.Px(), p_in_raw.Py(), p_in_raw.Pz());
  const Vec3 k_vec(k_in_raw.Px(), k_in_raw.Py(), k_in_raw.Pz());
  const P4 p_in = jet_tools::round_beam(p_vec, kProtonMass, kProtonNominalPz, kCrossingAngle);
  const P4 k_in = jet_tools::round_beam(k_vec, kElectronMass, kElectronNominalPz, 0.0);

  const std::string particles_name = use_reco ? collections.reco_lab_collection : "GeneratedParticles";
  if (!has_collection(frame, particles_name)) {
    return false;
  }
  const auto &particles = frame.get<edm4eic::ReconstructedParticleCollection>(particles_name);
  if (particles.empty()) {
    return false;
  }

  const int scat_idx = jet_tools::find_scattered_electron_index(frame, particles, use_reco, collections, cuts);
  if (scat_idx < 0 || static_cast<std::size_t>(scat_idx) >= particles.size()) {
    return false;
  }

  const auto &scat = particles[static_cast<std::size_t>(scat_idx)];
  const auto m = scat.getMomentum();
  double e_f = scat.getEnergy();
  if (!std::isfinite(e_f) || e_f <= 0.0) {
    e_f = std::sqrt(m.x * m.x + m.y * m.y + m.z * m.z + kElectronMass * kElectronMass);
  }
  const P4 k_out(m.x, m.y, m.z, e_f);

  jet_tools::BreitXSource x_source = jet_tools::BreitXSource::Derived;
  double x_input = std::numeric_limits<double>::quiet_NaN();
  const std::string kin_name =
      use_reco ? collections.electron_kinematics : "InclusiveKinematicsTruth";
  if (has_collection(frame, kin_name)) {
    const auto &kin = frame.get<edm4eic::InclusiveKinematicsCollection>(kin_name);
    if (!kin.empty()) {
      const double x = kin[0].getX();
      if (std::isfinite(x) && x > 0.0) {
        x_source = jet_tools::BreitXSource::Input;
        x_input = x;
      }
    }
  }

  return jet_tools::transform_breit(p_in, k_in, k_out, boost, rot, x_source, x_input);
}

void match_domain_event(ULong64_t event_id, const DomainEventJets &jets,
                        const ROOT::Math::Boost &boost,
                        const ROOT::Math::Rotation3D &rot, double max_dR2,
                        jet_tools::CrossFrameBestMatchRow &best_row,
                        jet_tools::CrossFrameAllMatchRow &all_row,
                        TTree &best_tree, TTree &all_tree) {
  // Loop through all Centauro jets (per event).
  for (std::size_t ci = 0; ci < jets.centauro_breit.size(); ++ci) {
    const jet_tools::SimpleJet cent_lab = jet_tools::breit_to_lab(jets.centauro_breit[ci], boost, rot);

    std::vector<std::pair<double, int>> accepted;
    accepted.reserve(jets.antikt_lab.size());

    for (std::size_t ai = 0; ai < jets.antikt_lab.size(); ++ai) {
      const auto &aj = jets.antikt_lab[ai];
      const double deta = aj.eta - cent_lab.eta;
      const double dphi = jet_tools::delta_phi(aj.phi, cent_lab.phi);
      const double dR2 = deta * deta + dphi * dphi;

      // Skip if not a valid value or beyond dR threshold.
      if (!std::isfinite(dR2) || dR2 > max_dR2) {
        continue;
      }
      accepted.push_back({dR2, static_cast<int>(ai)});
    }

    // Sort candidates by increasing dR to this Centauro jet.
    std::sort(accepted.begin(), accepted.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    
    // First fill all matches within dR.
    all_row.event = event_id;
    all_row.centauro_index = static_cast<int>(ci);
    all_row.antikt_indices.clear();
    all_row.antikt_indices.reserve(accepted.size());
    for (const auto &pair : accepted) {
      all_row.antikt_indices.push_back(pair.second);
    }
    all_tree.Fill();

    // If at least one anti-kt candidate passed the dR cut, the first sorted entry is the nearest one,
    // so use it as the best match for this Centauro jet.
    if (!accepted.empty()) {
      best_row.event = event_id;
      best_row.centauro_index = static_cast<int>(ci);
      best_row.antikt_index = accepted.front().second;
      best_row.dR = std::sqrt(accepted.front().first);
      best_tree.Fill();
    }
  }
}

} // namespace

int main(int argc, char *argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  TFile f_jets(args.jets_input.c_str(), "READ");
  if (!f_jets.IsOpen()) {
    std::cerr << "[match_jets_cf] failed to open jets input " << args.jets_input
              << "\n";
    return 1;
  }

  //========================================================================================
  // Grab and read jet information.
  //========================================================================================
  auto *t_truth_cent = jet_tools::get_required_tree(
      f_jets, "BreitFrameCentauroTruth", "match_jets_cf");
  auto *t_truth_antikt = jet_tools::get_required_tree(
      f_jets, "LabFrameAntiktTruth", "match_jets_cf");
  auto *t_reco_cent = jet_tools::get_required_tree(
      f_jets, "BreitFrameCentauroReco", "match_jets_cf");
  auto *t_reco_antikt = jet_tools::get_required_tree(
      f_jets, "LabFrameAntiktReco", "match_jets_cf");

  // Group jets by event.
  jet_tools::EventJets truth_cent_jets;
  jet_tools::EventJets truth_antikt_jets;
  jet_tools::EventJets reco_cent_jets;
  jet_tools::EventJets reco_antikt_jets;
  jet_tools::JetTreeReadOptions jet_read_options;
  jet_read_options.require_four_momentum = true;
  jet_tools::read_jet_tree(*t_truth_cent, args.max_events, truth_cent_jets,
                           "match_jets_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_truth_antikt, args.max_events, truth_antikt_jets,
                           "match_jets_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_reco_cent, args.max_events, reco_cent_jets,
                           "match_jets_cf", jet_read_options);
  jet_tools::read_jet_tree(*t_reco_antikt, args.max_events, reco_antikt_jets,
                           "match_jets_cf", jet_read_options);

  DomainMap truth_events;
  DomainMap reco_events;
  for (const auto& event_entry : truth_cent_jets) {
    truth_events[event_entry.first].centauro_breit = event_entry.second;
  }
  for (const auto& event_entry : truth_antikt_jets) {
    truth_events[event_entry.first].antikt_lab = event_entry.second;
  }
  for (const auto& event_entry : reco_cent_jets) {
    reco_events[event_entry.first].centauro_breit = event_entry.second;
  }
  for (const auto& event_entry : reco_antikt_jets) {
    reco_events[event_entry.first].antikt_lab = event_entry.second;
  }

  podio::ROOTReader reader;
  reader.openFile(args.events_input);
  const std::size_t entries_total = reader.getEntries("events");
  const std::size_t entries = std::min(entries_total, args.max_events);

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[match_jets_cf] failed to open output " << args.output << "\n";
    return 1;
  }

  // Best matches and all matches within dR <= 0.3.
  TTree best_truth("best_matches_truth", "best truth cross-frame matches");
  TTree all_truth("all_matches_truth", "all truth cross-frame matches");
  TTree best_reco("best_matches_reco", "best reco cross-frame matches");
  TTree all_reco("all_matches_reco", "all reco cross-frame matches");

  //========================================================================================
  // Matched entries.
  //========================================================================================
  jet_tools::CrossFrameBestMatchRow best_row_truth;
  jet_tools::CrossFrameBestMatchRow best_row_reco;
  jet_tools::CrossFrameAllMatchRow all_row_truth;
  jet_tools::CrossFrameAllMatchRow all_row_reco;

  jet_tools::setup_best_match_tree(best_truth, best_row_truth);
  jet_tools::setup_all_match_tree(all_truth, all_row_truth);
  jet_tools::setup_best_match_tree(best_reco, best_row_reco);
  jet_tools::setup_all_match_tree(all_reco, all_row_reco);

  const jet_tools::JetToolsCollections collections;
  const jet_tools::ElectronVetoCuts cuts;
  const double max_dR2 = args.match_dR * args.match_dR;

  for (std::size_t i = 0; i < entries; ++i) {
    if (i % 100 == 0) {
      std::cout << "[match_jets_cf] event " << i << "/" << entries << "\r"
                << std::flush;
    }

    const ULong64_t event_id = static_cast<ULong64_t>(i);

    // Find the jet entries for event_id.
    const auto truth_iterator = truth_events.find(event_id);
    const auto reco_iterator = reco_events.find(event_id);

    const bool usable_truth = (truth_iterator != truth_events.end() &&
                             !truth_iterator->second.centauro_breit.empty() &&
                             !truth_iterator->second.antikt_lab.empty());
    const bool usable_reco = (reco_iterator != reco_events.end() &&
                            !reco_iterator->second.centauro_breit.empty() &&
                            !reco_iterator->second.antikt_lab.empty());
    if (!usable_truth && !usable_reco) {
      continue;
    }

    auto data = reader.readEntry("events", i);
    if (!data) {
      continue;
    }
    podio::Frame frame(std::move(data));

    // If transform(s) are valid, boost a Centauro jet back to lab frame, and proceed with matching.
    if (usable_truth) {
      ROOT::Math::Boost boost;
      ROOT::Math::Rotation3D rot;
      if (build_transform(frame, false, collections, cuts, boost, rot)) {
        // Use truth_iterator->second to grab the jet_data from DomainMap <event_id, jet_data>.
        match_domain_event(event_id, truth_iterator->second, boost, rot, max_dR2,
                           best_row_truth, all_row_truth, best_truth, all_truth);
      }
    }

    if (usable_reco) {
      ROOT::Math::Boost boost;
      ROOT::Math::Rotation3D rot;
      if (build_transform(frame, true, collections, cuts, boost, rot)) {
        match_domain_event(event_id, reco_iterator->second, boost, rot, max_dR2,
                           best_row_reco, all_row_reco, best_reco, all_reco);
      }
    }
  }
  std::cout << "[match_jets_cf] event " << entries << "/" << entries << "\n";

  best_truth.Write();
  all_truth.Write();
  best_reco.Write();
  all_reco.Write();

  fout.Close();
  f_jets.Close();
  return 0;
}
